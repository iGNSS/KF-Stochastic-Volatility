% ------------------------------------------------------------------- 
% Standard Kalman Filter Implementation 
%         Method: Conventional implementation
%           Type: Covariance filtering
%      Recursion: Riccati recursion
%           Form: Two stages, a priori form
%        Authors: R.E. Kalman
% Implementation: Maria Kulikova 
% ------------------------------------------------------------------- 
% Model: x_k = F x_{k-1} + B ControlInput_k + G noise1_k, noise1_k \sim N(0,Q)
%        z_k = H x_k + noise2_k, noise2_k \sim N(0,R) 
% ------------------------------------------------------------------- 
% References:
%   Kalman R.E. 
% 1. A new approach to linear filtering and prediction problems. 
%    Journal of basic Engineering. 1960 Mar, 82(1):35-45.
%
%   See also implementation in
% 2. Grewal, M.S., Andrews, A.P. (2015) 
%    Kalman filtering: theory and practice using MATLAB 
%    Prentice-Hall, New Jersey, 4th edn. 
% ------------------------------------------------------------------- 
% Input:
%     matrices        - system matrices F,H,Q etc
%     initials_filter - initials x0,P0
%     measurements    - measurements (where y(t_k) is the k-th column)
% Output:
%     neg_LLF     - negative log LF
%     predX       - a priori estimates (history) 
%     predDP      - diag of the predicted error covariance (history)
%     hatX        - filtered estimate (history) 
%     hatDP       - diag of the filtered error covariance (history)
% ------------------------------------------------------------------- 

function [neg_LLF,predX,predDP,hatX,hatDP] = Riccati_KF_standard(theta,parameters,matrices,initials_filter,Smeasurements,control_input)
   [symF,symB,symG,symQ,symH,symR] = deal(matrices{:});         % get system matrices
                  [symP,symX] = deal(initials_filter{:});       % get initials for the filter 

   % values at the current \theta
   [F,B,G,Q,H,R,P,X,measurements] = Substitute(parameters,theta,symF,symB,symG,symQ,symH,symR,symP,symX,Smeasurements); 

        [m,n]   = size(H);                    % dimensions
       N_total  = size(measurements,2);       % number of measurements
         hatX   = zeros(n,N_total);           % prelocate for efficiency
         hatDP  = zeros(n,N_total);           % prelocate for efficiency
         predX  = zeros(n,N_total+1);         % prelocate for efficiency
         predDP = zeros(n,N_total+1);         % prelocate for efficiency
       neg_LLF  = 1/2*m*log(2*pi)*N_total;    % set initial value for the neg Log LF
 predX(:,1)  = X;  predDP(:,1) = diag(P);     % save initials at the first entry
 for k = 1:N_total               
      [X,P,ek,Rek] = kf_update(X,P,measurements(:,k),H,R);     
      
      neg_LLF = neg_LLF+1/2*log((det(Rek)))+1/2*ek'/Rek*ek;
      hatX(:,k)  = X; hatDP(:,k) = diag(P); 

      [X,P]      = kf_predict(X,P,F,B,control_input(:,k),G,Q); 
    predX(:,k+1) = X; predDP(:,k+1) = diag(P); 
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Time update: a priori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P] = kf_predict(X,P,F,B,u,G,Q)
    X = F*X + B*u;            % Predicted state estimate  
    P = F*P*F' + G*Q*G';      % Predicted error covariance 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement update: a posteriori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P,residual,cov_residual] = kf_update(X,P,z,H,R)
  residual     = z - H*X;                 % residual
  cov_residual = R + H*P*H';              % residual covariance 
  Kalman_gain  = P*H'*inv(cov_residual);  % Filter gain

  X = X + Kalman_gain*residual;           % Filtered state estimate
  P = P - Kalman_gain*H*P;                % Filtered error covariance

end