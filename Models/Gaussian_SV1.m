% ------------------------------------------------------------------- 
%  Model: SV model in [1] with the Gaussian assumption before the transformation.
%   Transformation yields Equations 6a and 6b in the cited paper, i.e. 
%   log (y^2_t) = -1.27 + h_t + \varepsilon_t, \varepsilon_t \sim N(0,pi^2/2 )
%           h_t = omega + phi h_{t-1} + eta_t, eta_t \sim N(0,\sigma^2_{eta}) 
%  Conditional std (volatility): \sigma_t = exp(1/2 h_t)
%  Parameters: |phi|<1, \sigma^2_{eta}>0, omega is any
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
function [system_matrices,model_parameters,P0,x0,parameters_options,const_data] = Gaussian_SV1

% ---- model parameters  -------------------------------------
phi   = sym('phi','real');      % matrix F
covQ  = sym('covQ','real');     % \sigma^2_{eta} is the noise in process eq.
omega   = sym('omega','real');  % scaled volatility parameter

model_parameters = {'phi','covQ','omega'};

% ---- parameters' restrictions: sumA* parameters < sumB--------------
sumA =  [1,  0, 0;
        -1,  0, 0;
         0, -1, 0;];  % the number of columns is the number of param 
                      % number of rows is the number of constrains
sumB =  [1;1;0];      % it is a column with the size equals to the number of constratins       

Aeq = []; Beq = []; LB = []; UB = []; nonlcon = [];

parameters_options  = {sumA,sumB,Aeq,Beq,LB,UB,nonlcon};

% ---- process equation --------------------------------------
Fsys = [phi];  Gsys = [1]; Qsys = [covQ]; Bsys = [omega];

% ---- measurement equation -----------------------------------
Hsys = [1]; Rsys = pi^2/2;  const_data = -1.27;

system_matrices = {Fsys,Bsys,Gsys,Qsys,Hsys,Rsys};

% ---- filter initials  -------------------------------------
P0 = covQ/(1 - phi^2);    % if model is stationary,
x0 = omega/(1 - phi);     % then we know the mean and var

end
