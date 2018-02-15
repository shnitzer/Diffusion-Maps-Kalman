% MATLAB code implementation of the non-Linear object tracking example 
% from: T. Shnitzer and R. Talmon, J.J. Slotine, "Diffusion maps Kalman 
% filter", submitted to IEEE Transactions on Signal Processing.
% ***************************************************************@
% This implementation generates the underlying diffusion processes
% and the corresponding measurements of the non-linear tracking problem and 
% recovers the underlying processes using the proposed DMK method.
% Author: Tal Shnitzer.
% Created:  2/14/18.
% ***************************************************************@

%% Configuration
% Signal generation parameters:
len      = 1000;                % signal length
deltaT   = 0.01;                % time step
noiseStd = [0.5,1,1.5,2,2.5];   % relative standard deviation of noise
DriftMat = [-0.1,0; 0, -0.1];   % drift dynamics matrix
DriftVec = [0.001; 3];          % drift dynamics offset vector (c_1, c_2)
DMdim    = 2;                   % dimensions of the diffusion maps coordinates to use in the Kalman filter
iteNum   = 20;                  % number of process generation iterations for RMSE calulations
jumpErr  = 60;                  % processes with changes larger than jumpErr will be ignored

AlgNo    = 4;                   % number of compared algorithms (dmk, noisy measurements, optimal Kalman, EKF)

%% Initialization

rErr      = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = dmk, (:,2,:) = measurement error, (:,3,:) = EKF, (:,4,:) = particle filter
phiErr    = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = dmk, (:,2,:) = measurement error, (:,3,:) = EKF, (:,4,:) = particle filter
theta1Err = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = dmk, (:,2,:) = measurement error, (:,3,:) = EKF, (:,4,:) = particle filter
theta2Err = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = dmk, (:,2,:) = measurement error, (:,3,:) = EKF, (:,4,:) = particle filter

avgSNR = zeros(length(noiseStd),iteNum);

%% Generate / Load Process

h = waitbar(0, 'Please wait while the errors for different noise-STD are computed');
for nstd = 1:length(noiseStd)
    waitbar(nstd/length(noiseStd), h);
    noiseStd_curr = noiseStd(nstd);

    for iter = 1:iteNum
        
        goodProcess = 0;
        % Generating a process without discontinuities caused by arctan or
        % discontinuities larger than 'jumpErr':
        while ~goodProcess
            % Generate new underlying process:
            InitLoc    = 5*randn(2,1)+0.1;                                   % initial process location
            DriftRate  = drift(DriftVec,DriftMat);                           % set drift parameters
            DiffRate   = @(t,X) [sqrt(2), 0; 0 sqrt(2)];                     % set diffusion parameters
            SDE        = sde(DriftRate, DiffRate, 'StartState', InitLoc);    % define SDE
            [theta, t] = SDE.simulate(len-1, 'DeltaTime', deltaT);           % simulate process
            
            % Generate measurements:
            phiT = atand(theta(:,1)./theta(:,2)).';       % clean angle values
            rT   = sqrt(theta(:,1).^2 + theta(:,2).^2).'; % clean radius values
            
            % Ignoring processes which contain significant discontinuities due
            % to arctan and drastic chanes
            if ~any(abs(diff(phiT))>60) && ~any(theta(:,2)<0)
                goodProcess = 1;
            end
        end
        phiM = phiT + noiseStd_curr*std(phiT) * randn(size(phiT)); % noisy angle
        rM   = rT   + noiseStd_curr*std(rT)   * randn(size(rT));   % noisy radius
        
        avgSNR(nstd,iter) = ( var(rT) + var(phiT) ) / ( noiseStd_curr.^2*var(phiT) + noiseStd_curr.^2*var(rT) );
        
        yT = [phiT; rT];
        yM = [phiM; rM];
        
        tt = 100:size(yM,2); % samples to consider - ignoring the first samples in the error calculations due to initialization effect errors
        
        %% Constrcut diffusion maps with the modified mahalanobis distance
        
        % Compute the modified mahalanobis distance for the (noisy) measurements:
        mahDist       = modified_mahalanobis(yM);
        
        % Compute diffusion maps coordinates and eigenvalues:
        [psi, lambda] = diffusion_maps(mahDist, DMdim);
        
        % Apply DMK:
        [psi_hat, yDMK_est] = dmk(psi, lambda, yM, deltaT, tt);
        
        % Computing the errors for DMK:
        dmk_data.Tr_y  = yT; dmk_data.Tr_theta = theta;
        dmk_data.Est_y = yDMK_est;
        aI = 1;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter), ...
            theta1Err(nstd,aI,iter), theta2Err(nstd,aI,iter)] = error_calc(dmk_data, tt);
        
        %% Compute estimations of compared algorithms and noisy measurements
        
        % Measurement error:
        meas_data.Tr_y  = yT; meas_data.Tr_theta = theta;
        meas_data.Est_y = yM;
        aI = 2;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter), ...
            theta1Err(nstd,aI,iter), theta2Err(nstd,aI,iter)] = error_calc(meas_data, tt);
        
        % Extended Kalman Filter:
        y_est_ekf = ekf( yM, DriftMat, DriftVec, deltaT, noiseStd_curr*std(yT,[],2), InitLoc );
        % Extended Kalman Filter error:
        ekf_data.Tr_y  = yT; ekf_data.Tr_theta = theta;
        ekf_data.Est_y = y_est_ekf;
        aI = 3;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter), ...
            theta1Err(nstd,aI,iter), theta2Err(nstd,aI,iter)] = error_calc(ekf_data, tt);
        
        % Particle Filter:
        theta_est_PF = particle_filter( yM, DriftMat, DriftVec, deltaT, noiseStd_curr*std(yT,[],2), InitLoc );
        % Particle Filter error:
        pf_data.Tr_y      = yT; pf_data.Tr_theta = theta;
        pf_data.Est_theta = theta_est_PF.';
        aI = 4;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter), ...
            theta1Err(nstd,aI,iter), theta2Err(nstd,aI,iter)] = error_calc(pf_data, tt);
        
    end
end
close(h);

%% Displaying the results and comparing to other algorithms:

plot_rmse_err( theta1Err, theta2Err, avgSNR, {'DMK','Meas.','EKF','PF'} )

figure
ax(1) = subplot(2,1,1); plot(tt,yT(2,tt),'--b',tt,yDMK_est(2,tt),'g','LineWidth',1);
xlabel('samples'); ylabel('r')
hold on
subplot(2,1,1); scatter(tt,yM(2,tt),10,[0.7,0.7,0.7],'x');
legend('Clean measurements','DMK estimation','Noisy measurements');
hold off;

ax(2) = subplot(2,1,2); plot(tt,yT(1,tt),'--b',tt,yDMK_est(1,tt),'g','LineWidth',1);
xlabel('samples'); ylabel('\phi')
hold on
subplot(2,1,2); scatter(tt,yM(1,tt),10,[0.7,0.7,0.7],'x');
legend('Clean measurements','DMK estimation','Noisy measurements');
hold off;

suptitle('DMK filtering trace example')
linkaxes(ax,'x')



