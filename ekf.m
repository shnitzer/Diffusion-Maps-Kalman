% Compute the estimation of the Extended Kalman filter (with given system matrices)
% ***************************************************************@

function [ y_est ] = ekf( yM, DriftMat, DriftVec, deltaT, noiseSTD, InitLoc )
%EKF computes the Extended Kalman Filter estimation given the system
% matrices

% System dynamics (considering the time step)
DriftMat = DriftMat*deltaT + eye(size(DriftMat));
DriftVec = DriftVec*deltaT;

% System dynamics function:
f = @(y) [ atand( (DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1))) / (DriftVec(2)+DriftMat(2,2)*y(2)*cosd(y(1))) );...
    sqrt( (DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1)))^2 + (DriftVec(2)+DriftMat(2,2)*y(2)*cosd(y(1)))^2 )];
% System Measurement function:
h     = @(y) [ y(1); y(2) ];

% Derivatives of system dynamics and measurement functions:
diff_f = @(y) [ y(2)*( DriftMat(1,1)*DriftVec(2)*cosd(y(1))+DriftVec(1)*DriftMat(2,2)*sind(y(1))+DriftMat(1,1)*DriftMat(2,2)*y(2) ) / ( DriftVec(1)^2+y(2)*(DriftMat(1,1)*sind(y(1))*(2*DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1)))+2*DriftVec(2)*DriftMat(2,2)*cosd(y(1))+DriftMat(2,2)^2*y(2)*(cosd(y(1)))^2)+DriftVec(2)^2 ),...
    ( DriftMat(1,1)*DriftVec(2)*sind(y(1))-DriftVec(1)*DriftMat(2,2)*cosd(y(1)) ) / ( DriftVec(1)^2+y(2)*(DriftMat(1,1)*sind(y(1))*(2*DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1)))+2*DriftVec(2)*DriftMat(2,2)*cosd(y(1))+DriftMat(2,2)^2*y(2)*(cosd(y(1)))^2)+DriftVec(2)^2 );...
    ( DriftMat(1,1)*cosd(y(1))*(DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1))) - DriftMat(2,2)*sind(y(1))*(DriftVec(2)+DriftMat(2,2)*y(2)*cosd(y(1))) ) / sqrt((DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1)))^2 + (DriftVec(2)+DriftMat(2,2)*y(2)*cosd(y(1)))^2),...
    ( DriftMat(1,1)*sind(y(1))*(DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1))) + DriftMat(2,2)*cosd(y(1))*(DriftVec(2)+DriftMat(2,2)*y(2)*cosd(y(1))) ) / sqrt((DriftVec(1)+DriftMat(1,1)*y(2)*sind(y(1)))^2 + (DriftVec(2)+DriftMat(2,2)*y(2)*cosd(y(1)))^2) ];
diff_h = @(y) [ 1, 0; 0 1 ];

Q = [2, 0; 0, 2]*deltaT;                 % process noise covariance
R = [noiseSTD(1)^2, 0; 0 noiseSTD(2)^2]; % true measurement noise covariance

y_est = nan(size(yM));
P_est = nan(size(yM,1),size(yM,1),size(yM,2));

y_est(:,1)   = InitLoc;
P_est(:,:,1) = eye(size(yM,1));

% Extended Kalman Filter:
obj                             = extendedKalmanFilter(f,h,InitLoc,'ProcessNoise',Q,'MeasurementNoise',R);
obj.StateTransitionJacobianFcn  = diff_f;
obj.MeasurementJacobianFcn      = diff_h;

for ii = 2:size(yM,2)
    predict(obj);
    [y_est(:,ii),P_est(:,:,ii)] = correct(obj,yM(:,ii));
end

