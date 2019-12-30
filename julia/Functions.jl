module Functions

export divisors
export uni_utility_flat
export uni_utility
export utility_estimate
export dist_max

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

function Bb(n,L,α,β)
    prod1=1
    prod2=1
    den=1
    if α>1
        for i in 1:α-1
            prod1 *= n+α-i
        end
    end
    if β>1
        for i in 1:β-1
            prod2 *= L-n+β-i
        end
    end
    for j in 1:α+β-1
        den*= L+α+β-j
    end
    prod1*prod2*factorial(α+β-1)/(den*factorial(α-1)*factorial(β-1))
end

function cumulative_Bb(n,L,α,β)
    out = 0
    for i in 0:n
        out += Bb(i,L,α,β)
    end
    out
end

function uni_utility(C,M,α,β)
    L = Int64(C/M)
    out = 0
    for n in 0:L
        out+=(cumulative_Bb(n,L,α,β)^M - cumulative_Bb(n-1,L,α,β)^M)*((n+α)/(L+α+β))
    end
    out
end

## Markov-Chain Monte Carlo utility estimate
function utility_estimate(L,α,β,n=10000)
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
    out += maximum((ns.+α)./(L[indices].+(α+β)))/n
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
function dist_max(α,β,L_0,n=10000,util_precision=1E6)
    #We sort distributions to consider only distinct distributions,
    #and to propose perturbations that keep this sorting.
    L_old = sort(L_0,rev=true)
    u_old = 0
    u_old = utility_estimate(L_old,α,β,Int64(util_precision))
    for i in 1:n
        L_new = perturbation(L_old)
        u_new_1 = utility_estimate(L_new,α,β,Int64(1E4))
        #if rough estimate of new distribution is better,
        if u_new_1-u_old > 0
            #then calculate it better and compare it. 
            #Important to recalculate old utility in order to not get stuck in spurious results
            u_old = utility_estimate(L_old,α,β,Int64(util_precision))
            u_new = utility_estimate(L_new,α,β,Int64(util_precision))
            if u_new > 1.0005*u_old 
                L_old = L_new
                u_old = u_new
            end
        end
    end
    L_old,u_old
end 

end