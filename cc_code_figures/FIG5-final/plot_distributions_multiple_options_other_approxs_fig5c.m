

%plots stats of the trajectories

clear;
%clf;


%final plots 
%gradients (put also inputs for the associated approxs below):

%%% flat environment

X = load('distributions7_other_approxs4.m');
Z = load('value_actions_gradient7_other_approxs4.m');
Z_breadth = load('value_actions_gradient7_other_approxs_breadth4.m');
Z_depth = load('value_actions_gradient7_other_approxs_depth4.m');
Z_random = load('value_actions_gradient7_other_approxs_random4.m');
Z_triangular = load('value_actions_gradient7_other_approxs_triangular4.m');


N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );


figure(3)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

pos1 = [0.15 0.65 0.3 0.3];
subplot('Position',pos1)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
plot(Z(:,1),Z(:,6),'r'); %found by gradient
plot(Z(:,1),Z(:,7),'b--'); %optimal square root reward
plot(Z(:,1),Z(:,9),'k--'); %uniform, simulation, with correction becuase non-intenger sqrt(C)
ylim([0. 1])
xlabel('Capacity');
ylabel('Reward')

pos2 = [0.6 0.65 0.3 0.3];
subplot('Position',pos2)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
plot(Z(:,1),100*Z(:,6)./Z(:,9)-100,'k-'); %numerical optimal square root sim, with correction becuase non-intenger sqrt(C) 
    %in the new version ...5 (not in ...4), the comparison rule 9 is uniform for N<=7, and
    %then square root with correction for inexact square root
xlabel('Capacity');
ylabel('% reward gain');


pos3 = [0.15 0.2 0.3 0.3];
subplot('Position',pos3)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
%comparisons to optimal (found by gradient descent):
plot(Z(:,1),100*Z(:,9)./Z(:,6)-100,'k-'); %compared square root (with corrections) approx low to optimal
plot(Z(:,1),100*Z_breadth(:,2)./Z(:,6)-100,'Color',[1 0 0]); %breadth comparison 
plot(Z(:,1),100*Z_depth(:,2)./Z(:,6)-100,'Color',[1 0 1]); %depth comparison
plot(Z(:,1),100*Z_random(:,2)./Z(:,6)-100,'Color',[1 0.5 0]); %random comparison
plot(Z(:,1),100*Z_triangular(:,2)./Z(:,6)-100,'Color',[0.5 0.5 0.5]); %triangular comparison
xlabel('Capacity');
ylabel('% reward loss');

%legends
text(0.55,1.00,'triangular','Units','normalized','Color',[0.5 0.5 0.5],'FontSize',6,...
    'FontName','Times New Roman');
text(0.55,0.95,'square root','Units','normalized','Color',[0 0 0],'FontSize',6,...
    'FontName','Times New Roman');
text(0.55,0.90,'random','Units','normalized','Color',[1 0.5 0],'FontSize',6,...
    'FontName','Times New Roman');
text(0.55,0.85,'pure breadth','Units','normalized','Color',[1 0 0],'FontSize',6,...
    'FontName','Times New Roman');
text(0.55,0.80,'pure depth','Units','normalized','Color',[1 0 1],'FontSize',6,...
    'FontName','Times New Roman');




pos4 = [0.6 0.2 0.3 0.3];
subplot('Position',pos4)
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
plot(N_vec,frac_actions_vec); %fractions of 0s
plot(N_vec(10:num_N),frac_max_asymptotic(10:num_N),'k');
xlabel('Capacity');
ylabel('fraction of zeros')


%print pdf
print('fig5_approx_panel_fig5c','-dpdf')


