%-------------------------------------------------------------------------
% Script for illustrating the recovered volatility process 
% and the confidence interval
% Implementation: Maria Kulikova
%-------------------------------------------------------------------------
function Illustrate_Volatility(std_vol,std_vol_UB,std_vol_LB,Datas,Data_Name)

if nargin<3, Data_Name = 'Data'; end;
N_total  = size(std_vol,2);

points   = N_total/24;                        % plot each two years;
xti      = round(linspace(2,N_total,points));
xti_lab  = Datas(xti,:);

figure;
  x_ax = 1:length(std_vol_LB);   
  X_plot  = [x_ax, fliplr(x_ax)];
  Y_plot  = [std_vol_LB, fliplr(std_vol_UB)];
  fill(X_plot,Y_plot,1,'facecolor',[0.9 0.9 0.9],'edgecolor','none'); 
   hold on;

  plot(1:N_total,std_vol,'k-','LineWidth',2);
  title(Data_Name,'FontSize',14);
  xlabel('Time'); ylabel('Data')

  set(gca,'XTick',xti,'XTickLabel',xti_lab,'FontSize',10,'XGrid','on','XLim',[1 N_total]); 
  rotateXLabels(gca,45);

end
