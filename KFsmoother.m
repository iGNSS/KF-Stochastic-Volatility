% ------------------------------------------------------------------- 
% Kalman Filter Implementation + smoothing
% ------------------------------------------------------------------- 

function [neg_LLF,predX,predDP,hatX,hatDP,smoothX,smoothDP] = KFsmoother(theta,parameters,matrices,initials_filter,Smeasurements,control_input)
   [symF,symB,symG,symQ,symH,symR] = deal(matrices{:});    % get system matrices
                 [symP,symX] = deal(initials_filter{:});  % get initials for the filter 

   % values at the current \theta
   [F,B,G,Q,H,R,P,X,measurements] = Substitute(parameters,theta,symF,symB,symG,symQ,symH,symR,symP,symX,Smeasurements); 
   
        [m,n]  = size(H);                    % dimensions
       N_total = size(measurements,2);       % number of measurements
         hatX  = zeros(n,N_total);           % prelocate for efficiency
         hatDP = zeros(n,N_total);           % prelocate for efficiency
         predX  = zeros(n,N_total+1);        % prelocate for efficiency
         predDP = zeros(n,N_total+1);        % prelocate for efficiency
       neg_LLF = 1/2*m*log(2*pi)*N_total;    % set initial value for the neg Log LF
 predX(:,1)  = X; predDP(:,1) = diag(P);     % save initials at the first entry
 for k = 1:N_total               
      [X,P,ek,Rek] = kf_update_Joseph(X,P,measurements(:,k),H,R);     
      
      neg_LLF = neg_LLF+1/2*log((det(Rek)))+1/2*ek'/Rek*ek;
      hatX(:,k) = X; hatDP(:,k) = diag(P); 
      PP(:,:,k) = P;

      [X,P]      = kf_predict(X,P,F,B,control_input(:,k),G,Q); 
    predX(:,k+1) = X; predDP(:,k+1) = diag(P); 
  end;
   [smoothX,smoothDP] = rts_smooth(hatX,PP,F,B,control_input,G*Q*G');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P] = kf_predict(X,P,F,B,u,G,Q)
    X = F*X + B*u;            % Predicted state estimate  
    P = F*P*F' + G*Q*G';      % Predicted error covariance 
end

function [X,P,residual,cov_residual] = kf_update_Joseph(X,P,z,H,R)

  residual     = z - H*X;                   % residual
  cov_residual = R + H*P*H';                % residual covariance 
  Kalman_gain  = P*H'/cov_residual;         % Filter Gain 

  X = X + Kalman_gain*residual;             % Filtered state estimate
  A = (eye(size(P,1)) - Kalman_gain*H);
  P = A*P*A' + Kalman_gain*R*Kalman_gain';  % Filtered error covariance
end


function [M,smoothDP,D] = rts_smooth(M,P,A,B,u,Q)

  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
  end;

  %
  % Extend A and Q if they are NxN matrices
  %
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);
  end;
  if size(Q,3)==1
    Q = repmat(Q,[1 1 size(M,2)]);
  end;

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  for k=(size(M,2)-1):-1:1
    P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
    D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - A(:,:,k) * M(:,k) - B*u(:,k));
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
    smoothDP(:,k) = diag(P(:,:,k));
  end;
end

