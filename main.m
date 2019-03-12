% MATLAB code implementation of the non-linear object tracking example 1
% from: T. Shnitzer and R. Talmon, J.J. Slotine, "Diffusion maps Kalman 
% Filter for a Class of Systems with Gradient Flows", submitted to IEEE 
% Transactions on Signal Processing.
% ***************************************************************@
% This implementation generates the underlying diffusion processes
% and the corresponding measurements of the non-linear tracking problem 1
% and recovers the underlying processes using the proposed DMK method.
% Author: Tal Shnitzer.
% Created:  3/10/19.
% ***************************************************************@

function main
%MAIN creates Fig. 1 and Fig. 3 in the paper


gendata = 0; % 0 - use saved data to produce the figure in the paper, 1 - generate a new realization

DMdim   = 2;              % dimensions of the diffusion maps coordinates to use in the Kalman filter
deltaT  = 0.01;           % time step

if gendata
    len        = 1000;    % process length
    procStd    = sqrt(2); % standard deviation of the process noise
    noiseStd   = 1;       % standard deviation of the measurement noise
    InitLoc    = 1*randn(2,1)+[1; 5];        % initial process location
    DriftRate1 = @(t,X) -0.5*(X-1).^3+(X-1); % set drift parameters
    DriftRate2 = @(t,X) -0.5*(X-6).^3+(X-6); % set drift parameters
    DiffRate1  = @(t,X) procStd;             % set diffusion parameters
    DiffRate2  = @(t,X) procStd;             % set diffusion parameters
    SDE1       = sde(DriftRate1, DiffRate1, 'StartState', InitLoc(1)); % define SDE
    [thet1, ~] = SDE1.simulate(len-1, 'DeltaTime', deltaT);            % simulate process
    SDE2       = sde(DriftRate2, DiffRate2, 'StartState', InitLoc(2)); % define SDE
    [thet2, ~] = SDE2.simulate(len-1, 'DeltaTime', deltaT);            % simulate process
    
    theta = [thet1, thet2];
    % Generate measurements:
    phiT = atan(theta(:,1)./theta(:,2)).';       % clean angle values
    rT   = sqrt(theta(:,1).^2 + theta(:,2).^2).'; % clean radius values
    phiM = phiT + noiseStd*std(phiT) * randn(size(phiT)); % noisy angle
    rM   = rT   + noiseStd*std(rT)   * randn(size(rT));   % noisy radius
    yT   = [phiT; rT];
    yM   = [phiM; rM];
else
    load('data.mat');
    InitLoc = theta(1,:).';   % true initial location for the particle filter
end

% Plot state trajectory:
figure
plot(theta(:,1),theta(:,2),'k')
grid on; axis equal
xlabel('$$x_1$$','Interpreter','latex','FontSize',16)
ylabel('$$x_2$$','Interpreter','latex','FontSize',16)

%% Constrcut diffusion maps with the modified mahalanobis distance

% Compute the modified mahalanobis distance for the (noisy) measurements:
mahDist       = modified_mahalanobis(yM);

% Compute diffusion maps coordinates and eigenvalues:
[psi, lambda] = diffusion_maps(mahDist, DMdim);

%% DMK framework:

% Apply DMK:
[~, yDMK_est] = dmk(psi, lambda, yM, deltaT);

%% Compute estimations of the particle filter:

DriftRate1 = @(t,X) -0.5*(X-1).^3+(X-1); % set true drift parameters
DriftRate2 = @(t,X) -0.5*(X-6).^3+(X-6); % set true drift parameters

% Particle Filter:
y_est_pf = particle_filter( yM, DriftRate1, DriftRate2, deltaT, noiseStd(1)*std(yT,[],2), procStd, InitLoc );

%% Plot trajectory examples - for each SNR:
    
tt = 100:size(yM,2); % samples to consider - ignoring the first samples in the error calculations due to initialization effect errors

figure
scatter(tt*deltaT,yM(2,tt),20,[0.6,0.6,0.6],'x');
hold on
xlabel('t [sec]','FontSize',14); ylabel('$$r$$','Interpreter','latex','FontSize',16)
plot(tt*deltaT,yT(2,tt),':k','LineWidth',2);
plot(tt*deltaT,y_est_pf(2,tt),'Color',[0.5,0.5,0.5],'LineWidth',2);
plot(tt*deltaT,yDMK_est(2,tt),'b','LineWidth',1);
lgd = legend('Noisy meas.','Clean meas.','PF estimation','DMK estimation');
lgd.FontSize = 12;
if ~gendata; xlim([2,5]); ylim([2,10]); end
hold off;
set(gcf,'Position',[45,90,560,285]);

figure
scatter(tt*deltaT,yM(1,tt),20,[0.6,0.6,0.6],'x');
xlabel('t [sec]','FontSize',14); ylabel('$$\phi$$','Interpreter','latex','FontSize',16)
hold on
plot(tt*deltaT,yT(1,tt),':k','LineWidth',2);
plot(tt*deltaT,y_est_pf(1,tt),'Color',[0.5,0.5,0.5],'LineWidth',2);
plot(tt*deltaT,yDMK_est(1,tt),'b','LineWidth',1);
if ~gendata; xlim([2,5]); ylim([-0.6,0.8]); end
hold off;
set(gcf,'Position',[45,90,560,285]);

end

