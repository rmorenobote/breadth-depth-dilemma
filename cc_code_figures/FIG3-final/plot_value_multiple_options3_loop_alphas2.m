%goes with the code ...._loop_alphas

clear;
Y_th = load('value_actions_th2th_loop_alphas_opt.m'); 
mean_post_vec = sort( unique(Y_th(:,1)) ); %%alpha/(alpha+beta), mean posterior
num_mean_post = length( mean_post_vec );

for k=1:num_mean_post
    
%clf;
N_vec_th = sort( unique(Y_th(:,2)) );
num_N_th = length( N_vec_th );
    for i=1:num_N_th 
       index_th = find( Y_th(:,1) == mean_post_vec(k) & Y_th(:,2) == N_vec_th(i) ); 

       %exact
       M_index_th3(i) = Y_th(index_th,3);
       frac_index_th3(i) = Y_th(index_th,3) / N_vec_th(i); 
    end
    
    %finding the critical C (C before a drop):
    i_max = min( find( frac_index_th3(:) < 1 ) ) - 1;
    C_critical_vec(k) = N_vec_th(i_max);

    %finding the exponent (we use values above 200 to avoid transients):
    %power law fit
    N_vec_range = N_vec_th(1000:2000); %selecting large capacity only
    M_index_range = M_index_th3(1000:2000);
    X_reg = [log(N_vec_range) ones(length(log(N_vec_range)),1) ]; %adding column of ones
    [w,CI] = regress(log(M_index_range'),X_reg);
    exponent_vec(k) =  w(1);
    exponent_CI_vec(k) = CI(1,2)-w(1); %error estimate

end


figure(3)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 2.5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

%%% subplot 1 %%% 
pos3 = [0.15 0.2 0.3 0.6];
subplot('Position',pos3)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
plot(mean_post_vec,C_critical_vec,'k.');
plot(mean_post_vec,C_critical_vec,'k-');
xticks([0 0.25 0.5 0.75 1]);
ylim([0 50]);
xlabel('\alpha/(\alpha + \beta)');
ylabel('critical capacity')


%%% subplot 2 %%%
pos4 = [0.6 0.2 0.3 0.6];
subplot('Position',pos4)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
%errorbar(mean_post_vec,exponent_vec,exponent_CI_vec,'k-');
%shaded error bars
x = mean_post_vec;
y = exponent_vec';
dy = exponent_CI_vec';
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.7 .7 .7],'linestyle','none');
%
plot(mean_post_vec,exponent_vec,'k.');
plot(mean_post_vec,exponent_vec,'k-');
xticks([0 0.25 0.5 0.75 1]);
ylim([0.4 0.75]);
xlabel('\alpha/(\alpha + \beta)');
ylabel('exponent');


%print pdf
print('fig3_criticalC_power_vs_beta','-dpdf')

