%-------------------------------------------------------------------------
% Script for illustrating the log returns and recovered volatility process
% Conditional variance: \sigma_t = exp(1/2 \hat h_t)
% Implementation: Maria Kulikova
%-------------------------------------------------------------------------
function Illustrate_ReturnsVol(cc_returns,std_vol,Datas,Data_Name)

if nargin<3, Data_Name = 'Data'; end;
N_total  = size(cc_returns,2);

points   = N_total/24;                        % plot each 24 points;
xti      = round(linspace(2,N_total,points)); 
xti_lab  = Datas(xti,:);

figure;
  plot(1:N_total,abs(cc_returns),'k-','LineWidth',0.5);
  title(Data_Name,'FontSize',14);
  xlabel('Time'); ylabel('Data')
 
  hold on;
  plot(1:N_total,std_vol,'b-','LineWidth',2);

  set(gca,'XTick',xti,'XTickLabel',xti_lab,'FontSize',10,'XGrid','on','XLim',[1 N_total]); 
  rotateXLabels(gca,45);
end
