

%plots stats of the trajectories

clear;
%clf;

X = load('value_actions2.m');
X_th = load('value_actions_th2.m'); 

Y_th = load('value_actions_th2th.m'); 

N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );

N_vec_th = sort( unique(Y_th(:,1)) );
num_N_th = length( N_vec_th );


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
   
   plot(X(index,2),X(index,5),'k.'); %average reward, but reward is not generated randomly and averaged
   plot(X(index,2),X_th(index,5),'-','Color',color(i,:)); %exact flat prior
   %plot(X(index,2),X_th(index,6),'--'); %exact
   %plot(X(index,2),X(index,7),'k.'); %average reward, computed from stats
end
line([0 N_vec(13)],[2./3. 2./3.],'Color','k','LineStyle','--')
ylim([0.5 1])
text(0.8,0.70,'C=4','Units','normalized','Color',[0.8, 0.8, 0.8],'FontSize',8,...
    'FontName','Times New Roman');
text(0.8,0.80,'C=10','Units','normalized','Color',[0.6, 0.6, 0.6],'FontSize',8,...
    'FontName','Times New Roman');
text(0.8,0.90,'C=100','Units','normalized','Color',[0., 0., 0.],'FontSize',8,...
    'FontName','Times New Roman');
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
   max_vec(i) = max( X(index,5) ); 
   M_index_aux(i) = find( X(index,5) == max( X(index,5) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   
   index_th = find( X_th(:,1) == N_vec(i) ); 
   
   %exact flat prior
   max_vec_th(i) = max( X_th(index_th,5) ); 
   M_index_th_aux(i) = find( X_th(index_th,5) == max( X_th(index_th,5) ) ); 
   M_index_th(i) = X(index_th(M_index_th_aux(i)),2);
   
   %exact
   max_vec_th2(i) = max( X_th(index_th,6) ); 
   M_index_th_aux2(i) = find( X_th(index_th,6) == max( X_th(index_th,6) ) ); 
   M_index_th2(i) = X(index_th(M_index_th_aux2(i)),2);
   
   M_max_asymptotic(i) = (sqrt(N_vec(i)^3 - N_vec(i)) + N_vec(i)) / (N_vec(i) - 2);
   %M_max_asymptotic(i) = (sqrt(2/3*N_vec(i)));

end
for i=1:num_N_th 
   index_th = find( Y_th(:,1) == N_vec_th(i) ); 
    
   %exact
   max_vec_th3(i) = max( Y_th(index_th,5) ); 
   M_index_th_aux3(i) = find( Y_th(index_th,5) == max( Y_th(index_th,5) ) ); 
   M_index_th3(i) = Y_th(index_th(M_index_th_aux3(i)),2);
    
end
n_max = 20;
plot(N_vec(1:n_max),N_vec(1:n_max),'--','Color',[0 0.5 0.1]); %just a slope one straight line
%plot(N_vec,M_index);
plot(N_vec,M_index,'k.');
plot(N_vec,M_index_th,'k-');
%plot(N_vec,M_index_th2,'-');
%plot(N_vec_th,M_index_th3,'-.');
%plot(N_vec(10:num_N),M_max_asymptotic(10:num_N),'b--');
%power law fit
N_vec_range = N_vec(11:20); %selecting large capacity only
M_index_range = M_index(11:20);
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
   
   max_vec(i) = max( X(index,5) ); 
   M_index_aux(i) = find( X(index,5) == max( X(index,5) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
   
   index_th = find( X_th(:,1) == N_vec(i) ); 
   
   %exact flat prior
   max_vec_th(i) = max( X_th(index_th,5) ); 
   M_index_aux_th(i) = find( X_th(index_th,5) == max( X_th(index_th,5) ) ); 
   M_index_th(i) = X(index_th(M_index_aux_th(i)),2);
   frac_index_th(i) = M_index_th(i) / N_vec(i); 
   
   %exact
   max_vec_th2(i) = max( X_th(index_th,6) ); 
   M_index_aux_th2(i) = find( X_th(index_th,6) == max( X_th(index_th,6) ) ); 
   M_index_th2(i) = X(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
   
   frac_max_asymptotic(i) = (sqrt(N_vec(i)^3-N_vec(i)) + N_vec(i))/(N_vec(i) - 2) / N_vec(i);
end
for i=1:num_N_th 
   index_th = find( Y_th(:,1) == N_vec_th(i) ); 
    
   %exact
   max_vec_th3(i) = max( Y_th(index_th,5) ); 
   M_index_th_aux3(i) = find( Y_th(index_th,5) == max( Y_th(index_th,5) ) ); 
   M_index_th3(i) = Y_th(index_th(M_index_th_aux3(i)),2);
   frac_index_th3(i) = M_index_th3(i) / N_vec_th(i); 
    
end
%plot(N_vec,frac_index);
plot(N_vec,frac_index,'k.');
plot(N_vec,frac_index_th,'k-');
%plot(N_vec,frac_index_th2,'-');
%plot(N_vec_th,frac_index_th3,'-.');
%plot(N_vec(10:num_N),frac_max_asymptotic(10:num_N),'k');
x = [7 7];
y = [0 1];
pl = line(x,y); %vertical linea at 7
pl.Color = 'r';
pl.LineStyle = '--';
xticks([7 10^3 10^5]) 
xticklabels({'7', '10^3', '10^5'}) 
xlabel('Capacity');
ylabel('optimal M / C')

%%%%%%%%%%%%%
%subplot 4 %%
%%%%%%%%%%%%%

%flat/uniform beta
pos4 = [0.6 0.2 0.3 0.3];
subplot('Position',pos4)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,5) ); 
   M_index_aux(i) = find( X(index,5) == max( X(index,5) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
   
   index_th = find( X_th(:,1) == N_vec(i) ); 
   
   %exact flat prior
   max_vec_th(i) = max( X_th(index_th,5) ); 
   M_index_aux_th(i) = find( X_th(index_th,5) == max( X_th(index_th,5) ) ); 
   M_index_th(i) = X(index_th(M_index_aux_th(i)),2);
   frac_index_th(i) = M_index_th(i) / N_vec(i); 
   
   %exact
   max_vec_th2(i) = max( X_th(index_th,6) ); 
   M_index_aux_th2(i) = find( X_th(index_th,6) == max( X_th(index_th,6) ) ); 
   M_index_th2(i) = X(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
   
   frac_max_asymptotic(i) = (sqrt(N_vec(i)^3-N_vec(i)) + N_vec(i))/(N_vec(i) - 2) / N_vec(i);
end
plot(N_vec,frac_index,'k.');
plot(N_vec,frac_index_th,'k-');

N_vec_range = N_vec(11:18); %selecting large capacity only
M_index_range = M_index(11:18);
X_reg = [log(N_vec_range) ones(length(log(N_vec_range)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_range'),X_reg);
w(1)
CI

%Gaussian like beta 
clear;
X = load('value_actions2b.m'); % gaussian like beta
X_th = load('value_actions_th2b.m'); 

N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );

for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,5) ); 
   M_index_aux(i) = find( X(index,5) == max( X(index,5) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
   
   index_th = find( X_th(:,1) == N_vec(i) ); 
      
   %exact
   max_vec_th2(i) = max( X_th(index_th,6) ); 
   M_index_aux_th2(i) = find( X_th(index_th,6) == max( X_th(index_th,6) ) ); 
   M_index_th2(i) = X(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
   
end
plot(N_vec,frac_index,'.','Color',[0 0.5 0]);
plot(N_vec,frac_index_th2,'Color',[0 0.5 0]);

N_vec_range = N_vec(11:18); %selecting large capacity only
M_index_range = M_index(11:18);
X_reg = [log(N_vec_range) ones(length(log(N_vec_range)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_range'),X_reg);
w(1)
CI

%Negatively skeweed beta 
clear;
X = load('value_actions2c.m'); % gaussian like beta
X_th = load('value_actions_th2c.m'); 

N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );

for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,5) ); 
   M_index_aux(i) = find( X(index,5) == max( X(index,5) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
   
   index_th = find( X_th(:,1) == N_vec(i) ); 
      
   %exact
   max_vec_th2(i) = max( X_th(index_th,6) ); 
   M_index_aux_th2(i) = find( X_th(index_th,6) == max( X_th(index_th,6) ) ); 
   M_index_th2(i) = X(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
   
end
plot(N_vec,frac_index,'.','Color',[0.5 0.2 0]);
plot(N_vec,frac_index_th2,'Color',[0.5 0.2 0]);
xticks([10 10^3 10^5]) 
xlabel('Capacity');
ylabel('optimal M / C')

N_vec_range = N_vec(11:18); %selecting large capacity only
M_index_range = M_index(11:18);
X_reg = [log(N_vec_range) ones(length(log(N_vec_range)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_range'),X_reg);
w(1)
CI

%Positively skeweed beta 
clear;
X = load('value_actions2d.m'); % gaussian like beta
X_th = load('value_actions_th2d.m'); 

N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );

for i=1:num_N
   index = find( X(:,1) == N_vec(i) ); 
   
   max_vec(i) = max( X(index,5) ); 
   M_index_aux(i) = find( X(index,5) == max( X(index,5) ) ); 
   M_index(i) = X(index(M_index_aux(i)),2);
   frac_index(i) = M_index(i) / N_vec(i); 
   
   index_th = find( X_th(:,1) == N_vec(i) ); 
      
   %exact
   max_vec_th2(i) = max( X_th(index_th,6) ); 
   M_index_aux_th2(i) = find( X_th(index_th,6) == max( X_th(index_th,6) ) ); 
   M_index_th2(i) = X(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
   
end
plot(N_vec,frac_index,'b.');
plot(N_vec,frac_index_th2,'b-');

x = [7 7];
y = [0 1];
pl = line(x,y); %vertical linea at 7
pl.Color = 'r';
pl.LineStyle = '--';

N_vec_range = N_vec(11:18); %selecting large capacity only
M_index_range = M_index(11:18);
X_reg = [log(N_vec_range) ones(length(log(N_vec_range)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_range'),X_reg);
w(1)
CI

%legends
text(0.55,0.70,'+ skewed / poor','Units','normalized','Color',[0 0 1],'FontSize',8,...
    'FontName','Times New Roman');
text(0.55,0.80,'- skewed / rich','Units','normalized','Color',[0.5 0.2 0],'FontSize',8,...
    'FontName','Times New Roman');
text(0.55,0.90,'Gaussian-like','Units','normalized','Color',[0 0.5 0],'FontSize',8,...
    'FontName','Times New Roman');
text(0.55,1.00,'beta prior:','Units','normalized','Color',[0 0 0],'FontSize',8,...
    'FontName','Times New Roman');

xticks([10 10^3 10^5]) 
xlabel('Capacity');
ylabel('optimal M / C')


%%%%%%%%%%%%%%%%%
%NEW way of computing powers of fits with the ...th2th code and many points: 
%%%%%%%%%%%%%%%%%%

figure(11)

%flat
clear;

X_th = load('value_actions_th2th_flat.m'); 

N_vec = sort( unique(X_th(:,1)) );
num_N = length( N_vec );

hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
for i=1:num_N
   index_th = find( X_th(:,1) == N_vec(i) ); 
   
   %exact
   max_vec_th2(i) = max( X_th(index_th,5) ); 
   M_index_aux_th2(i) = find( X_th(index_th,5) == max( X_th(index_th,5) ) ); 
   M_index_th2(i) = X_th(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
end
plot(N_vec,frac_index_th2,'k');
xticks([10 10^3 10^5]) 
xlabel('Capacity');
ylabel('optimal M / C')

X_reg = [log(N_vec) ones(length(log(N_vec)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_th2'),X_reg);
w(1)
CI

%rich
clear;

X_th = load('value_actions_th2th_rich.m'); 

N_vec = sort( unique(X_th(:,1)) );
num_N = length( N_vec );

hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
for i=1:num_N
   index_th = find( X_th(:,1) == N_vec(i) ); 
   
   %exact
   max_vec_th2(i) = max( X_th(index_th,5) ); 
   M_index_aux_th2(i) = find( X_th(index_th,5) == max( X_th(index_th,5) ) ); 
   M_index_th2(i) = X_th(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
end
plot(N_vec,frac_index_th2,'Color',[0.5 0.2 0]);
xticks([10 10^3 10^5]) 
xlabel('Capacity');
ylabel('optimal M / C')

X_reg = [log(N_vec) ones(length(log(N_vec)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_th2'),X_reg);
w(1)
CI

%poor
clear;

X_th = load('value_actions_th2th_poor.m'); 

N_vec = sort( unique(X_th(:,1)) );
num_N = length( N_vec );

hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XScale', 'log')
for i=1:num_N
   index_th = find( X_th(:,1) == N_vec(i) ); 
   
   %exact
   max_vec_th2(i) = max( X_th(index_th,5) ); 
   M_index_aux_th2(i) = find( X_th(index_th,5) == max( X_th(index_th,5) ) ); 
   M_index_th2(i) = X_th(index_th(M_index_aux_th2(i)),2);
   frac_index_th2(i) = M_index_th2(i) / N_vec(i); 
end
plot(N_vec,frac_index_th2,'Color',[0 0 1]);
xticks([10 10^3 10^5]) 
xlabel('Capacity');
ylabel('optimal M / C')

X_reg = [log(N_vec) ones(length(log(N_vec)),1) ]; %adding column of ones
[w,CI] = regress(log(M_index_th2'),X_reg);
w(1)
CI



%print pdf
print('fig1_scaling','-dpdf')
