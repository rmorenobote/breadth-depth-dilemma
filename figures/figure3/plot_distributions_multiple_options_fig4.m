

%plots stats of the trajectories

clear;
%clf;

% flat
 X = load('distributions7_mixed.m');
 Z = load('value_actions_gradient7_mixed.m');

% - skewed
% X = load('distributions72_mixed.m');
% Z = load('value_actions_gradient72_mixed.m');

% + skewed
%X = load('distributions73_mixed.m');
%Z = load('value_actions_gradient73_mixed.m');



N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );



figure(3)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 2.5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

pos3 = [0.15 0.2 0.3 0.6];
subplot('Position',pos3)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   %number of 0s
   index2 = find(X(index,3) == 0);
   
   num_zeros_vec(i) = length(index2);
   frac_actions_vec(i) = 1 - length(index2)/N_vec(i);
   
  
   
   frac_max_asymptotic(i) = (sqrt(N_vec(i)^3-N_vec(i)) + N_vec(i))/(N_vec(i) - 2) / N_vec(i);
end
%plot(N_vec,num_zeros_vec);
plot(N_vec,frac_actions_vec,'Color',[0. 0. 1]); %fractions of 0s
%plot(N_vec,frac_actions_vec2); %fractions os 1s
%plot(N_vec(10:num_N),frac_max_asymptotic(10:num_N),'k');
xlabel('Capacity');
ylabel('fraction of sampled options')

%legends
text(0.55,0.70,'+ skewed / poor','Units','normalized','Color',[0 0 1],'FontSize',8,...
    'FontName','Times New Roman');
text(0.55,0.80,'- skewed / rich','Units','normalized','Color',[0.5 0.2 0],'FontSize',8,...
    'FontName','Times New Roman');
text(0.55,0.90,'uniform / flat','Units','normalized','Color',[0 0 0],'FontSize',8,...
    'FontName','Times New Roman','fontweight', 'bold');
text(0.55,1.00,'beta prior:','Units','normalized','Color',[0 0 0],'FontSize',8,...
    'FontName','Times New Roman');



%%% subplot 2 %%%
pos4 = [0.6 0.2 0.3 0.6];
subplot('Position',pos4)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
plot(Z(:,1),100*Z(:,6)./Z(:,9)-100,'Color',[0. 0. 1]); %numerical optimal square root sim, with correction becuase non-intenger sqrt(C) 
    %in the new version, the comparison rule 9 is uniform for N<=7, and
    %then square root with correction for inexact square root
%plot(Z(:,1),100*Z(:,6)./Z(:,8)-100,'b--'); %max uniform
%plot(Z(:,1),100*Z(:,6)./Z(:,7)-100,'r--'); %optimal square root, th, no correction becuase of non-integer sqrt(C) 
xlabel('Capacity');
ylabel('% improvement');








%print pdf
print('fig3_reward_improv_fraction_zeros','-dpdf')


