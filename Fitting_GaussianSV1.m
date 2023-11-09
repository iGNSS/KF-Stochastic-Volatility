% ------------------------------------------------------------------- 
% Script for fitting the univariate SV model by the state-space approach.
% Implementation: Maria Kulikova
% ------------------------------------------------------------------- 
%  Model: This is the SV model in [1] with the Gaussian assumption.
%   log (y^2_t) = h_t + \varepsilon_t, \varepsilon_t \sim N(0,\sigma^2_{varepsilon})
%           h_t = omega + phi h_{t-1} + eta_t, eta_t \sim N(0,\sigma^2_{eta}) 
%  Conditional std (volatility): \sigma_t = exp(1/2 h_t)
%  Parameters: |phi|<1, \sigma^2_{eta}>0, \sigma^2_varepsilon>0, omega is any
% ------------------------------------------------------------------- 
% Casted into the State-space form: 
%        x_k = F x_{k-1} + B ControlInput_k + G noise1_k, noise1_k \sim N(0,Q)
%        z_k = const_data + H x_k + noise2_k, noise2_k \sim N(0,R) 
% ------------------------------------------------------------------- 
% References: 
%  [1]  Harvey A., Ruiz E., Shephard N. (1994)
%       Multivariate stochastic variance models. 
%       The Review of Economic Studies. 1:61(2), 247-264.
%       DOI: https://doi.org/10.2307/2297980
% ------------------------------------------------------------------- 

clc; clear all; close all;  warning off;
MC_runs = 2; % Number of Monte Carlo runs in optimization

% ------- Load the data ----------------------------------------------
  Data_Name = 'SomeDataExample';      
  fprintf(1,'1. Load data and clean it from NaN and Inf. Sample:  %s \n',Data_Name);
  [prices, text,~] =  xlsread([Data_Name '.xlsx']);
  date_returns = char(text{:}); Close_price = prices(:,1);
  [Close_price,date_returns] = Clean_Data(Close_price,date_returns);

% --- Compute the cc_returns ------------------------------------------- 
  cc_returns    = log(Close_price(2:end)')-log(Close_price(1:end-1)');       % cc retunt series should be a row for our codes
  date_returns(1,:)  = [];                                                   % for cc_returns the data timeline is one less

% ---- Find the mean-adjusted returns ------------------------------------
 cc_returns_mean = mean(cc_returns,2);                % Subtracting the mean ensures that 
 Measurements = log((cc_returns-cc_returns_mean).^2); % there are no y_t identically equal to zero 

% --- Load  Model to be examined -----------------------------------------
p = pwd; cd('Models/'); 
[system_matrices,model_parameters,P0,x0,parameters_options,const_data] = Gaussian_SV1; cd(p);
 Measurements = Measurements - const_data; Input_signal = ones(1,length(Measurements));

% --- Choose the KF implementation to use for log LF calculation----------
p = pwd; cd('Methods-KF'); 
   handle_funs{1} = @Riccati_KF_Joseph;      % Conventional Joseph stabilized implementation
   % handle_funs{2} = you can add any other KF method in the same way 
cd(p); 
Number_Methods = size(handle_funs,2);       % number of methods to be tested
filters        = cell(Number_Methods,1);    % prelocate for efficiency

% --- Choose the optimization method to use ------------------------------
opt = optimset(@fmincon); 
opt = optimset(opt,'TolFun',1e-7,'TolX',1e-5,'Display','off','Algorithm','interior-point'); 
% 'MaxFunEvals',10000,'MaxIter',10000 and 'GradObj','off', Alg.'trust-region-reflective', 'interior-point' and Hessian: https://www.mathworks.com/help/optim/ug/hessian.html#bsapedt 

%% --- Monte Carlo runs --------------------------------------------------------------------
 for exp_number = 1:MC_runs  
    fprintf(1,'Fitting SV model #%d: \n',exp_number); 
     % ----- Optimization ------------------
     [sumA,sumB,Aeq,Beq,LB,UB,nonlcon] = deal(parameters_options{:}); 
     theta_initial = 0.01 * randi([60  99],1,3);
     fprintf(1,'      Model parameters: \t      '); fprintf(1,'%s \t    ',model_parameters{:}); fprintf(1,' \n');
     fprintf(1,'      Initial parametes:\t '); fprintf(1,'%8.4f\t ',theta_initial); fprintf(1,' \n');
     fprintf(1,'      Fitting model... \n');

    for i=1:Number_Methods;
     [optimum,min_LLF,flagexit,output,~,~,Hess] = fmincon(handle_funs{i},theta_initial,sumA,sumB,Aeq,Beq,LB,UB,nonlcon,opt,model_parameters,system_matrices,{P0,x0},Measurements,Input_signal);

     % --- Perform the filtering process at the optimum point ----
       filters{i}.legend = func2str(handle_funs{i}); 
       [LLF,predX,predDP,hatX,hatDP] = feval(handle_funs{i},optimum,model_parameters,system_matrices,{P0,x0},Measurements,Input_signal);
       filters{i}.predX(:,:,exp_number)  = predX;
       filters{i}.predDP(:,:,exp_number) = predDP;
       filters{i}.hatX(:,:,exp_number)  = predX;
       filters{i}.hatDP(:,:,exp_number) = predDP;
       filters{i}.yk(:,:,exp_number)    = Measurements;
       filters{i}.hatTheta(:,exp_number)  = optimum;
       filters{i}.initialTheta(:,exp_number)  = theta_initial;

      fprintf(1,'%d.%22s\t ',i,filters{i}.legend); 
      fprintf(1,'%8.4f\t ',optimum); fprintf(1,' \n');
    end;
end;
% ---- Find the a posteriori mean of the parameters ---------------------
for i=1:Number_Methods
 filters{i}.SE = std(filters{i}.hatTheta,0,2);
 mean_optimum = mean(filters{i}.hatTheta,2);
 [LLF,predX,predDP,hatX,hatDP] = feval(handle_funs{i},mean_optimum,model_parameters,system_matrices,{P0,x0},Measurements,Input_signal);
 filters{i}.neg_LLF = LLF;
end;
% --- Smooth at the a posteriori mean of the parameters -----------------
 [~,~,~,~,~,smX,smDP] = KFsmoother(mean_optimum,model_parameters,system_matrices,{P0,x0},Measurements,Input_signal);
% --- Print the results -------------------------------------------------  
fprintf(1,'--------------------- \n'); fprintf(1,'  Filter Implementations:    \t');       
fprintf(1,'%s \t   ',model_parameters{:}); fprintf(1,' | '); fprintf(1,'SE(%s) \t',model_parameters{:}); fprintf(1,' | ');  fprintf(1,'neg LLF \n');
for i=1:Number_Methods
 fprintf(1,'%d.%22s\t ',i,filters{i}.legend); 
 fprintf(1,'%8.4f\t',mean(filters{i}.hatTheta,2)); 
 fprintf(1,' | ');  fprintf(1,'%8.6f\t',filters{i}.SE);
 fprintf(1,' | '); fprintf(1,'%8.2f \n',mean(filters{i}.neg_LLF));
end;
% --- Plot the results -------------------------------------------------  
% Conditional std (volatility): \sigma_t = exp(1/2 h_t)
 std_vol = exp(1/2 * smX); % the smoothed volatility at the a posteriori mean of the parameter vector 
 Illustrate_ReturnsVol(abs(cc_returns),std_vol,date_returns,'Absolute values of first difference of logged index and smoothed estimate of standard deviation')
% --- smoothed volatility and its confidence interval ----------------
 alpha_val = 0.05;                % confidence at 5%
 cVal = norminv(1-alpha_val/2);   % confidence interval for the mean when sigma is known
 SE = sqrt(smDP);         
 UB = smX(:,2:end) + cVal*SE;  % upper bound
 LB = smX(:,2:end) - cVal*SE;  % lower bound
 std_vol_UB  = exp(1/2 * UB);
 std_vol_LB  = exp(1/2 * LB);
 Illustrate_Volatility(std_vol(:,2:end),std_vol_UB,std_vol_LB,date_returns,'Smoothed estimate of standard deviation (with confidence interval)')
