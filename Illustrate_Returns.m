function Illustrate_Returns(cc_returns,Datas,Data_Name)

if nargin<3, Data_Name = 'Data'; end;
N_total  = size(cc_returns,2);

points   = N_total/24;                        % plot each two years;
xti      = round(linspace(2,N_total,points));
xti_lab  = Datas(xti,:);

figure;
  plot(1:N_total,cc_returns,'k-','LineWidth',0.5);
  title(Data_Name,'FontSize',14);
  xlabel('Time'); ylabel('Data')

  set(gca,'XTick',xti,'XTickLabel',xti_lab,'FontSize',10,'XGrid','on','XLim',[1 N_total]); 
  rotateXLabels(gca,45);

end
