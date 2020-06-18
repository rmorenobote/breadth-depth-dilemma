%this codes takes outputs from code ...capacity_fixed.cc

%plots stats of the trajectories

clear;
%clf;

%%%unit variance
X = load('value_actions_gauss_capacity_fixed2.m');
 
%%%low variance
X2 = load('value_actions_gauss_capacity_fixed22.m');

%%%high variance
X3 = load('value_actions_gauss_capacity_fixed23.m');

N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );


figure(5)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

pos1 = [0.15 0.65 0.3 0.3];
subplot('Position',pos1)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
color(4,:) = [0.8, 0.8, 0.8];
color(10,:) = [0.6, 0.6, 0.6];
color(100,:) = [0., 0., 0.];
for i=[4,10,100]

   index = find( X(:,1) == i ); 
   
   plot(X(index,2),X(index,4),'k.'); %aerage reward, but reward is not generated randomly and averaged
   plot(X(index,2),X(index,4),'-','Color',color(i,:)); %aerage reward, but reward is not generated randomly and averaged
end
line([0 N_vec(13)],[2./3. 2./3.],'Color','k','LineStyle','--')
ylim([0.5 1])
text(0.8,0.70,'C=4','Units','normalized','Color',[0.8, 0.8, 0.8],'FontSize',8,...
    'FontName','Times New Roman');
text(0.8,0.80,'C=10','Units','normalized','Color',[0.6, 0.6, 0.6],'FontSize',8,...
    'FontName','Times New Roman');
text(0.8,0.90,'C=100','Units','normalized','Color',[0., 0., 0.],'FontSize',8,...
    'FontName','Times New Roman');
%xticks([0 20 40 60 80 100])
xticks([1 2 4 10 20 50 100])
yticks([0.5 0.6 0.7 0.8 0.9 1]) 
%xticklabels({'3\pi'})
xlabel('M, sampled alternatives');
ylabel('average reward')


pos2 = [0.6 0.65 0.3 0.3];
subplot('Position',pos2)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   %data, simulations
   max_vec(i) = max( X(index,4) ); 
   M_index_aux(i) = find( X(index,4) == max( X(index,4) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
end
n_max = 23;
plot(N_vec(1:n_max),N_vec(1:n_max),'--','Color',[0 0.5 0.1]); %just a slope one straight line
plot(N_vec,M_index,'k.');
plot(N_vec,M_index,'k-');
%power law fit
N_vec_range = N_vec(12:23); %selecting large capacity only
M_index_range = M_index(12:23);
X_reg = [log(N_vec_range) ones(length(log(N_vec_range)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_range'),X_reg);
w(1)
plot(N_vec_range,exp(w(2))*N_vec_range.^w(1),'r--');
%xticks([1 10 100 10^3 10^4 10^5]) 
xticks([10 10^3 10^5]) 
yticks([10 10^3 10^5]) 
xlabel('Capacity');
ylabel('optimal M')


pos3 = [0.15 0.2 0.3 0.3];
subplot('Position',pos3)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,4) ); 
   M_index_aux(i) = find( X(index,4) == max( X(index,4) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
end
plot(N_vec,frac_index,'k.');
plot(N_vec,frac_index,'k-');
x = [9 9];
y = [0 1];
pl = line(x,y); %vertical linea at 10
pl.Color = 'r';
pl.LineStyle = '--';
xticks([9 10^2 10^3]) 
xticklabels({'9' '10^2', '10^2'}) 
xlabel('Capacity');
ylabel('optimal M / C')


%%%% plot here all variances
pos4 = [0.6 0.2 0.3 0.3];
subplot('Position',pos4)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
%unit variance
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,4) ); 
   M_index_aux(i) = find( X(index,4) == max( X(index,4) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
end
plot(N_vec,frac_index,'k.');
plot(N_vec,frac_index,'k-');

%low variance
X = X2;
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,4) ); 
   M_index_aux(i) = find( X(index,4) == max( X(index,4) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
end
plot(N_vec,frac_index,'b.');
plot(N_vec,frac_index,'b-');

%high variance
X = X3;
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,4) ); 
   M_index_aux(i) = find( X(index,4) == max( X(index,4) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
end
plot(N_vec,frac_index,'r.');
plot(N_vec,frac_index,'r-');

pl.Color = 'r';
pl.LineStyle = '--';
xticks([10^1 10^2 10^3]) 
xticklabels({'10^1', '10^2', '10^3'}) 
xlabel('Capacity');
ylabel('optimal M / C')




%print pdf
print('figSI1_rewards_scaling','-dpdf')



