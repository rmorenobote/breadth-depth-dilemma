

%plots optimal distributions

clear;
%clf;

%final plots
X = load('distributions6.m'); %flat beta prior

%X = load('distributions62.m'); %negatively skewed, rich environment
 
%X = load('distributions63.m'); %positively skewed, needle in a haystack,
%poor environemnt


N_vec = sort( unique(X(:,1)) );
num_N = length( N_vec );

figure(2)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

hold on;
L = floor(sqrt(num_N+1))+1;
%L = min(L,18);
for i=1:num_N 
   
   subplot(L,L,i);
   
   index = find( X(:,1) == N_vec(i) ); 
   
   p1 = bar(-sort(-X(index,3)),'k');
   set(p1,'FaceColor','b');  
   %%%%% colors: flat = k;
   %%%%% - skeweed = [0.5 0.2 0];
   %%%%% + skeweed = b
   hold on;
   plot(1+X(index,2),-sort(-X(index,3)),'k.');
   
   axis tight
   if i <= 8
        xticks([1 2 3 4 5 6 7 8 9 10]);
   end
   if i > 8
        xticks([1 5 10 15 20]);
   end
   yticks([0 1 2 3 4 5]);
   
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
print('fig2_opt_distributions','-dpdf')


