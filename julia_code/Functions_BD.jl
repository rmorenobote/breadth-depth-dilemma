module Functions_BD

export divisors
export uni_utility_flat
export uni_utility_Bb
export utility_estimate
export utility_sim_Bb
export dist_max
export dist_max_Gauss
export utility_uniform_Gauss
export utility_sim_Gauss
export policy_iteration
export policy_iteration_Bb

using Distributions, SpecialFunctions

function divisors(n)
    out = Int64[]
    for i in 1:n
        if rem(n,i) == 0
            push!(out,i)
        end
    end
    out
end

function sum_of_powers(C,M)
    L = Int64(C/M)
    out = 0
    for s in 1:L
        out += s^big(M)
    end
    out
end

#Utility for flat prior, analytic result
function uni_utility_flat(C,M)
    L = Int64(C/M)
    out = 0
    p_sum = sum_of_powers(C,M)
    out = (C+M)/(C+2*M)-p_sum/((L+2)*(L+1)^big(M))
    out
end

prob_Bb(x,n,a,b) = exp(loggamma(n+1) + loggamma(x+a) + loggamma(n-x+b) + loggamma(a+b) - (loggamma(x+1) + loggamma(n-x+1) + loggamma(a) + loggamma(b) + loggamma(n+a+b)))

function cumulative_Bb(n,L,α,β)
    out = 0
    for i in 0:n
        out += prob_Bb(i,L,α,β)
    end
    out
end

function uni_utility_Bb(C,M,α,β,alice = false)
    L = Int64(C/M)
    out = 0
    cumulative1 = 0
    cumulative2 = 0
    for n in 0:L
        cumulative1 += prob_Bb(n,L,α,β)#pdf(BetaBinomial(L,α,β),n)
        if n >= 1
            cumulative2 += prob_Bb(n-1,L,α,β)#pdf(BetaBinomial(L,α,β),n-1)
        end
        if alice == true
            out+=(cumulative1^M - cumulative2^M)*max((n+α)/(L+α+β),α/(α+β))            
        else
            out+=(cumulative1^M - cumulative2^M)*((n+α)/(L+α+β))
        end
    end
    out
end

#Function that estimates the utility of distribution L by simulating outcomes
function utility_sim_Bb(L,α,β,N_sim = 1E5,alice=false)
    mean_new = 0
    var_new = 0
    indices = Any[]
    for j in 1:length(L)
        if L[j]>=1
            push!(indices,j)
        end
    end
    M = length(indices)
    for k in 1:N_sim
        mean_old = mean_new
        var_old = var_new
        post = zeros(M)
        p = zeros(M)
        for i in 1:M
            p[i] = rand(Beta(α,β))
            n = 0
            for j in 1:L[i]
                if rand() < p[i]
                    n += 1
                end
            end
            post[i] = (α + n)/(L[i]+α+β)
        end
        idx = findmax(post)[2]
        a = p[idx]
        if alice == true
            x = max(a,α/(α+β))
        else
            x = a
        end
        mean_new = mean_old + (x-mean_old)/k
        if k > 1
            var_new = var_old + (x-mean_old)*(x-mean_new)/(k-1)
        end
    end
    mean_new, var_new
end

## Markov-Chain Monte Carlo utility estimate
function utility_estimate(L,α,β,n=10000,alice=false)
    out = 0
    indices = Any[]
    ns = Int64[]
    empty = 0
    for j in 1:length(L)
        if L[j]>=1
            push!(indices,j)
            push!(ns,rand(0:L[j]))
        else
            empty += 1
        end
    end
    for i in 1:n
        prop_j = rand(1:length(indices))
        ratio = 0
        which = 0
        if ns[prop_j] == 0
            if rand() <= 0.5
                ratio = α*L[indices][prop_j]/(L[indices][prop_j]+β-1)
                which = 1
            end
        elseif ns[prop_j] == L[indices][prop_j]
            if rand()<= 0.5
                ratio = L[indices][prop_j]*β/(L[indices][prop_j]+α-1)
                which = -1
            end
        else
            if rand() <= 0.5
                ratio = ((ns[prop_j]+α)/(L[indices][prop_j]-ns[prop_j]+β-1))*((L[indices][prop_j]-ns[prop_j])/(ns[prop_j]+1))
                which = 1
            else
                ratio = ((L[indices][prop_j]-ns[prop_j]+β)/(ns[prop_j]+α-1))*((ns[prop_j])/(L[indices][prop_j]-ns[prop_j]+1))
                which = -1
            end
        end
        alpha = min(1,ratio)
        if rand()<= alpha
            ns[prop_j]+=which
        end
        if alice == true
            out += max(maximum((ns.+α)./(L[indices].+(α+β))),α/(α+β))/n
        else
            out += maximum((ns.+α)./(L[indices].+(α+β)))/n
        end
    end
    out
end

#Alternative way to estimate utility
function ind_utility_estimate(L,α,β,n=100)
    out2 = 0
    indices = Any[]
    reps = 1
    for j in 1:length(L)
        if L[j]>=1
            push!(indices,j)
            reps *= (L[j]+1)
        end
    end

    for n_i in 1:n
        out = 0
        for i in 1:reps
            ns = zeros(Int64,length(indices))
            prod = 1
            for j in 1:length(indices)
                ns[j] = rand(0:L[j])
                prod *= factorial(L[j])*factorial(α+β-1)*factorial(ns[j]+α-1)*factorial(L[j]-ns[j]+β-1)/(factorial(L[j]-ns[j])*factorial(ns[j])*factorial(α-1)*factorial(β-1)*factorial(L[j]+α+β-1))
            end

            out += prod*maximum((ns.+α)./(L[indices].+(α+β)))
        end
        out2+=out/n
    end
    out2
end

#Function that produces distribution perturbations
function perturbation(L)
    if L[end] >= 1
        indices_to_take = [length(L)]
    else
        indices_to_take = Any[]
    end
    #Determine which are the indices from which we can take 
    #and which are the ones we can put samples into.
    #Distributions have to be always ordered.
    indices_to_put = [1]
    for j in 1:length(L)-1
        if L[j]-L[j+1]>=1
            push!(indices_to_take,j)
            push!(indices_to_put,j+1)
        end
    end
    prop_take = 0
    prop_put = 0
    while abs(prop_take-prop_put)<1
        #Propose two random indices, one from which we take a sample
        #and the other two put it into
        take,put = rand(indices_to_take),rand(indices_to_put)
        if put-take == 1 && L[take]-L[put]==1
            prop_take,prop_put = 0,0
        else
            prop_take,prop_put = take,put
        end
    end
    L_new = copy(L)
    L_new[prop_take] -= 1
    L_new[prop_put] += 1
    L_new
end

#Gradient-descent method to find optimal distribution
function dist_max(α,β,L_0,n=10000,util_precision=1E6,alice = false)
    #We sort distributions to consider only distinct distributions,
    #and to propose perturbations that keep this sorting.
    L_old = sort(L_0,rev=true)
    u_old = 0
    u_old = utility_estimate(L_old,α,β,Int64(util_precision),alice)
    for i in 1:n
        L_new = perturbation(L_old)
        u_new_1 = utility_estimate(L_new,α,β,Int64(1E4),alice)
        #if rough estimate of new distribution is better,
        if u_new_1-u_old > 0
            #then calculate it better and compare it. 
            #Important to recalculate old utility in order to not get stuck in spurious results
            u_old = utility_estimate(L_old,α,β,Int64(util_precision),alice)
            u_new = utility_estimate(L_new,α,β,Int64(util_precision),alice)
            if u_new > 1.0005*u_old 
                L_old = L_new
                u_old = u_new
            end
        end
    end
    L_old,u_old
end 

function greedy_policy_Bb(policy, ϵ, C)
    terminate = false
    L = Int64[]
    t = 1
    current_a = C
    #Explore with probability ϵ the C first possibilities
    if rand() < ϵ
        current_a = policy[1]
    else
        current_a = rand(1:C)
    end
    #push!(L,current_a)
    while terminate == false
        b = C - sum(L)
        #print("b = ", b, "\n")
        next_a = 1
        if b > 0
            #Explore allocations for all times
            if rand() < ϵ || t > length(policy)
                next_a = rand(1:current_a)
            else
                #Can't choose more than previous
                next_a = min(policy[t],current_a)               
            end
            #Action is either determined by policy or the amount of capacity we have left
            current_a = Int64(min(next_a,b))
            push!(L,current_a)
            t += 1
        else
            #push!(L,0)
            terminate = true
        end
    end
    L
end

function policy_iteration_Bb(seed, ϵ, C, α, β, n_episodes, alice = false)
    Q = Float64[]
    policies = Any[]
    Qn = Float64[]
    counter = Int64[]
    policy = seed
    for ep in 1:n_episodes
        p = zeros(C)
        #Initialize options
        for i in 1:C
            p[i] = rand(Beta(α,β))
        end
        R = 0
        #ϵ-greedy policy
        L = greedy_policy_Bb(policy,ϵ,C)
        M = length(L)
        post = zeros(M)
        #n = zeros(M)
        #Simulate with policy
        for j in 1:M
            #current_a = L[j]
            n = 0
            for k in 1:L[j]
                if rand() < p[j]
                    n += 1
                end
            end
            post[j] = (α + n)/(L[j]+α+β)
        end 
        id = findmax(post)[2]
        #If we are allowed to pick from unsampled options
        if alice == true
            #Choose the option with the best expected probability of success, and that is reward
            R = max(p[id],α/(α+β))
            #If we can only pick from the options we sampled
        else
            R = p[id]
        end
        #Only update q-value functions for the visited action, which is given by 
        idx = findall(x -> x == L, policies)
        id = 0
        #If it hasn't been evaluated before, add it
        if length(idx) == 0
            push!(policies, L)
            id = length(policies)
            push!(Q, 0)
            push!(Qn, 0)
            push!(counter,1)
        else
            id = idx[1]
            counter[id] += 1
        end   
        Q[id] = Qn[id] + (R - Qn[id])/counter[id]
        #Improve policy after each episode
        best_idx = findmax(Q)[2]
        policy = policies[best_idx]
        #Reinitilize Qn
        Qn = copy(Q)
    end 
    policy, policies, Q, counter, length(policies)
end

####################################################################################
########################## Gaussian distributed samples ############################
####################################################################################

#Cumulative probability function
cum_flat(x,t,σ) = (1+(x/t)*(erf(x/(sqrt(2*σ^2*t)))-erf((x-t)/(sqrt(2*σ^2*t))))+erf((x-t)/(sqrt(2*σ^2*t)))+sqrt(2*σ^2/(t*pi))*(exp(-x^2/(2*σ^2*t))-exp(-(x-t)^2/(2*σ^2*t))))/2

#Probability density function
prob_flat(x,t,σ)=(erf(x/(σ*sqrt(2*t)))-erf(sqrt(t)*(x/t - 1)/(σ*sqrt(2))))/(2*t)

#Analytical, numerical integration of utility
function utility_uniform_Gauss(C,N,σ,N_points = 1E5)
    L = C/N
    x_max = 1000
    Δx = 2*x_max/N_points
    x = collect(-x_max:Δx:x_max)
    total_x = length(x)
    out = 0
    cumulative = 0
    for i in 1:total_x
        cumulative = cum_flat(x[i]*sqrt(2*σ^2*L),L,σ)
        out += cumulative^(N-1)*(2*cumulative - (1 + erf(x[i]-sqrt(L/(2*σ^2)))))*Δx
    end
    out *= sqrt(σ^2*N^3/(2*C))
    out
end

#Monte-Carlo estimation of utility of L distribution
function utility_sim_Gauss(L,σ,N_sim = 1E5,alice=false)
    indices = Any[]
    for j in 1:length(L)
        if L[j]>=1
            push!(indices,j)
        end
    end
    M = length(indices)
    out = 0
    mean_new = 0
    var_new = 0
    for k in 1:N_sim
        mu_ut = 0
        post_old = -99999
        for i in 1:M
            mu = rand()
            x = 0
            for j in 1:L[i]
                x += rand(Normal(mu,σ))
            end
            post = x/L[i] + σ*(exp(-x^2/(2*σ^2*L[i])) - exp(-(sqrt(L[i]/(2*σ^2)) - x/sqrt(2*σ^2*L[i]))^2))/(sqrt(2*pi*L[i]^3)*prob_flat(x,L[i],σ))
            if post > post_old
                post_old = post
                mu_ut = mu
            end
        end
        if alice == true
            max_mu = max(mu_ut,0.5)
        else
            max_mu = mu_ut
        end
        mean_new += max_mu
    end
    mean_new/N_sim
end

#Function that finds the best distribution for Gaussian distributed samples
function dist_max_Gauss(L_0,σ,n=100,util_precision=1E5,alice = false)
    #We sort distributions to consider only distinct distributions,
    #and to propose perturbations that keep this sorting.
    L_old = sort(L_0,rev=true)
    u_old = 0
    u_old = utility_sim_Gauss(L_old,σ,Int64(util_precision),alice)
    for i in 1:n
        L_new = perturbation(L_old)
        #print("A perturbation \n")
        u_new_1 = utility_sim_Gauss(L_new,σ,Int64(1E4),alice)
        #if rough estimate of new distribution is better,
        if u_new_1-u_old > 0
            #then calculate it better and compare it. 
            #Important to recalculate old utility in order to not get stuck in spurious results
            u_old = utility_sim_Gauss(L_old,σ,Int64(util_precision),alice)
            u_new = utility_sim_Gauss(L_new,σ,Int64(util_precision),alice)
            if u_new > 1.0005*u_old 
                L_old = L_new
                u_old = u_new
            end
        end
    end
    L_old,u_old
end 

####################################################################################
################################# Dynamic allocation ###############################
####################################################################################

function greedy_policy(policy, ϵ, C)
    terminate = false
    t = 1
    t_max = C
    actions = Int64[]
    R = 0
    current_a = C
    #Explore with probability ϵ the C first possibilities
    if rand() < ϵ
        current_a = policy[1]
    else
        current_a = rand(1:C)
    end
    while terminate == false
        b = C - sum(actions)
        #print("b = ", b, "\n")
        next_a = 1
        if b > 0
            #Explore allocations for all times
            if rand() < ϵ || policy[t] == 0.0 || t > length(policy)
                next_a = rand(1:current_a)
            else
                #Can't choose more than previous
                next_a = min(policy[t],current_a)               
            end
            #Action is either determined by policy or the amount of capacity we have left
            current_a = Int64(min(next_a,b))
            push!(actions,current_a)
        else
            push!(actions,0)
            t_max = t
            terminate = true
        end
    end
    actions
end

function policy_iteration(seed, ϵ, C, α, β, n_episodes, alice = false)
    Q = Float64[]
    policies = Any[]
    Qn = Float64[]
    counter = Int64[]
    policy = seed
    for ep in 1:n_episodes
        n = zeros(C)
        L = zeros(C)
        p = zeros(C)
        #Initialize options
        for i in 1:C
            p[i] = rand(Beta(α,β))
        end
        R = 0
        #ϵ-greedy policy
        actions = greedy_policy(policy,ϵ,C)
        t_max = length(actions)
        #Simulate with policy
        for t in 1:t_max
            current_a = actions[t]
            if t < t_max
                for i in 1:current_a
                    L[i] += 1
                    if rand() < p[i]
                        n[i] += 1
                    end
                end
                #Sort them at every time step
                indices = sortperm(n,rev = true)
                n = n[indices]
                L = L[indices]
                p = p[indices] 
                t += 1
            else
                #We are allowed to pick from unsampled options
                if alice == true
                    #Choose the option with the best expected probability of success, and that is reward
                    id = findmax( (n .+α) ./(L .+ (α + β)))[2]
                    R = p[id]
                #We can only pick from the options we sampled in the first wave
                else
                    ids = findall(x -> L[x] > 0, 1:C)
                    id = findmax( (n[ids] .+α) ./(L[ids] .+ (α + β)))[2]
                    R = p[id]
                end
            end
        end
        #Only update q-value functions for the visited action, which is given by 
        idx = findall(x -> x == actions, policies)
        id = 0
        #If it hasn't been evaluated before, add it
        if length(idx) == 0
            push!(policies, actions)
            id = length(policies)
            push!(Q, 0)
            push!(Qn, 0)
            push!(counter,1)
        else
            id = idx[1]
            counter[id] += 1
        end   
        Q[id] = Qn[id] + (R - Qn[id])/counter[id]
        #Improve policy after each episode
        best_idx = findmax(Q)[2]
        policy = policies[best_idx]
        #Reinitilize Qn
        Qn = copy(Q)
    end 
    policy, policies, Q, counter, length(policies)
end

end