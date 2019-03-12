% Construct the diffusion maps kalman filter
% ***************************************************************@

function [psi_hat, y_est] = dmk(psi, lambda, yM, deltaT)
%DMK constructs and applies the diffusion maps kalman filter to the noisy
% measurements 'yM', based on the eigenvertors 'psi' and eignevalues
% 'lambda'. The function returned the estimated system state 'psi_hat' and
% estimated clean measurement 'y_est'.

DMdim = size(psi,2);

%% Kalman parameters:
F = eye(DMdim) + diag(lambda)*deltaT;       % DMK dynamics matrix
H = yM * psi;                               % DMK measurement mapping
Q = diag(var(psi*diag(lambda)));            % estimation of process covariance
R = diag(var(yM.'));                        % estimation of measurement covariance

%% Initializations:
psi_hat      = nan(DMdim, size(yM,2));
P_est        = nan(DMdim, DMdim, size(yM,2));
psi_hat(:,1) = psi(1,:).';
P_est(:,:,1) = eye(DMdim);

%% Construct DMK:
for ii = 2:(length(yM))
    % Gain matrix update:
    K = (F*P_est(:,:,ii-1)*F.' + Q)*H.' / (H*F*P_est(:,:,ii-1)*F.'*H.' + H*Q*H.' + R);
    
    % State estimation and state covariance update:
    psi_hat(:,ii) = F*psi_hat(:,ii-1) + K*(yM(:,ii) - H*F*psi_hat(:,ii-1));
    P_est(:,:,ii) = (eye(DMdim) - K*H) * (F*P_est(:,:,ii-1)*F.' + Q);
    
    % Enforcing symmetry and PSD in P_est for stability:
    P_est(:,:,ii) = (P_est(:,:,ii)+P_est(:,:,ii).')/2;
    [Vtmp,Dtmp]   = eig(P_est(:,:,ii));
    P_est(:,:,ii) = Vtmp*abs(Dtmp)/Vtmp;
end

% Lift the new representation back to the measurement space:
y_est = (H * psi_hat);

end