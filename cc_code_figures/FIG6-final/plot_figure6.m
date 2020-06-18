

%plots stats of the trajectories

clear;
%clf;


%final plots 
%gradients (put also inputs for the associated approxs below):

%%% flat environment

%%% All capacities included here from files below
X = load('distributions7_dyn_alloc5_other_approxs_mixed.m');
Z = load('value_actions_gradient7_dyn_alloc5_other_approxs_mixed.m');
Z_breadth = load('value_actions_gradient7_dyn_alloc5_other_approxs_breadth_mixed.m');
Z_depth = load('value_actions_gradient7_dyn_alloc5_other_approxs_depth_mixed.m');
Z_random = load('value_actions_gradient7_dyn_alloc5_other_approxs_random_mixed.m');
Z_triangular = load('value_actions_gradient7_dyn_alloc5_other_approxs_triangular_mixed.m');
Z_sq_root = load('value_actions_gradient7_dyn_alloc5_other_approxs_sq_root_mixed.m');

%%% Capacity from 1 to 13
% X = load('distributions7_dyn_alloc5_other_approxs23.m');
% Z = load('value_actions_gradient7_dyn_alloc5_other_approxs23.m');
% Z_breadth = load('value_actions_gradient7_dyn_alloc5_other_approxs_breadth23.m');
% Z_depth = load('value_actions_gradient7_dyn_alloc5_other_approxs_depth23.m');
% Z_random = load('value_actions_gradient7_dyn_alloc5_other_approxs_random23.m');
% Z_triangular = load('value_actions_gradient7_dyn_alloc5_other_approxs_triangular23.m');
% Z_sq_root = load('value_actions_gradient7_dyn_alloc5_other_approxs_sq_root23.m');

%%% Capacity from 14 to 20
% X = load('distributions7_dyn_alloc5_other_approxs24.m');
% Z = load('value_actions_gradient7_dyn_alloc5_other_approxs24.m');
% Z_breadth = load('value_actions_gradient7_dyn_alloc5_other_approxs_breadth24.m');
% Z_depth = load('value_actions_gradient7_dyn_alloc5_other_approxs_depth24.m');
% Z_random = load('value_actions_gradient7_dyn_alloc5_other_approxs_random24.m');
% Z_triangular = load('value_actions_gradient7_dyn_alloc5_other_approxs_triangular24.m');
% Z_sq_root = load('value_actions_gradient7_dyn_alloc5_other_approxs_sq_root24.m');

%%% Capacity from 14 to 20
% X = load('distributions7_dyn_alloc5_other_approxs24.m');
% Z = load('value_actions_gradient7_dyn_alloc5_other_approxs24.m');
% Z_breadth = load('value_actions_gradient7_dyn_alloc5_other_approxs_breadth24.m');
% Z_depth = load('value_actions_gradient7_dyn_alloc5_other_approxs_depth24.m');
% Z_random = load('value_actions_gradient7_dyn_alloc5_other_approxs_random24.m');
% Z_triangular = load('value_actions_gradient7_dyn_alloc5_other_approxs_triangular24.m');
% Z_sq_root = load('value_actions_gradient7_dyn_alloc5_other_approxs_sq_root24.m');

%%% Capacity up to 2500
% X = load('distributions7_dyn_alloc5_other_approxs2.m');
% Z = load('value_actions_gradient7_dyn_alloc5_other_approxs2.m');
% Z_breadth = load('value_actions_gradient7_dyn_alloc5_other_approxs_breadth2.m');
% Z_depth = load('value_actions_gradient7_dyn_alloc5_other_approxs_depth2.m');
% Z_random = load('value_actions_gradient7_dyn_alloc5_other_approxs_random2.m');
% Z_triangular = load('value_actions_gradient7_dyn_alloc5_other_approxs_triangular2.m');
% Z_sq_root = load('value_actions_gradient7_dyn_alloc5_other_approxs_sq_root2.m');


N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );



%%% Optimal Allocations
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

hold on;
for i=1:20 %%num_N 
   
   subplot(5,5,i); 
   
   index = find( X(:,1) == N_vec(i) ); 
   
   bar(-sort(-X(index,3)),'k');
   hold on;
   plot(1+X(index,2),-sort(-X(index,3)),'k.');
   
   axis tight
   if i <= 8
        xticks([1 2 3 4 5 6 7 8 9 10]);
   end
   if i > 8
        xticks([1 5 10 15 20]);
   end
   yticks([0 1 2 3 4 5 6 8 10]);
   
   if i>3
       text(0.5,0.70,['C=',num2str(i)],'Units','normalized','Color',[0., 0., 0.],'FontSize',8,...
        'FontName','Times New Roman');
   end
   
   set(gca,'fontsize',9)
   set(gca, 'FontName', 'Times New Roman')
end
%xlabel('option #');
%ylabel('# samples')

%print pdf
print('fig6_optimal_allocations_fig6c','-dpdf')




%%% Optimal Waves
figure(3)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

hold on;
L = floor(sqrt(num_N+1))+1;
%L = min(L,18);
for i=1:20 %%num_N 
   
   %subplot(L,L,i);
   subplot(5,5,i);
   
   index = find( X(:,1) == N_vec(i) ); 
   
   p1 = bar(-sort(-X(index,4)),'k');
   set(p1,'FaceColor','[0.6 0.2 0.2]');  
   
   hold on;
   plot(1+X(index,2),-sort(-X(index,4)),'k.');
   
   axis tight
   if i <= 8
        xticks([1 2 3 4 5 6 7 8 9 10]);
   end
   if i > 8
        xticks([1 5 10 15 20]);
   end
   if i <= 8
        yticks([1 2 3 4 5 6 7 8 9 10]);
   end
   if i > 8
        yticks([0 2 4 6 8 10]);
   end
   
   
   if i>3
       text(0.5,0.70,['C=',num2str(i)],'Units','normalized','Color',[0., 0., 0.],'FontSize',8,...
        'FontName','Times New Roman');
   end
   
   set(gca,'fontsize',9)
   set(gca, 'FontName', 'Times New Roman')
end
%xlabel('wave #');
%ylabel('# samples')

%print pdf
print('fig6_optimal_waves_fig6b','-dpdf')



figure(4)

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
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   %number of 0s
   index2 = find(X(index,3) == 0);
   
   num_zeros_vec(i) = length(index2);
   frac_actions_vec(i) = 1 - length(index2)/N_vec(i);
end
plot(N_vec,frac_actions_vec,'k'); %fractions of 0s
xlabel('Capacity');
ylabel('fraction of sampled options')


pos4 = [0.6 0.2 0.3 0.3];
subplot('Position',pos4)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
%comparisons to optimal (found by gradient descent):
plot(Z(:,1),100*Z_breadth(:,2)./Z(:,6)-100,'Color',[1 0 0]); %breadth comparison 
plot(Z(:,1),100*Z_depth(:,2)./Z(:,6)-100,'Color',[1 0 1]); %depth comparison
plot(Z(:,1),100*Z_random(:,2)./Z(:,6)-100,'Color',[1 0.5 0]); %random comparison
plot(Z(:,1),100*Z_triangular(:,2)./Z(:,6)-100,'Color',[0.5 0.5 0.5]); %triangular comparison
plot(Z(:,1),100*Z_sq_root(:,2)./Z(:,6)-100,'k'); %compared to square root
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



%print pdf
print('fig6_fractions&reward_fig6_panels_a&b','-dpdf')







