function [ y_est ] = particle_filter( yM, DriftRate1, DriftRate2, deltaT, noiseSTD, procStd, InitLoc )
%PARTICLE_FILTER implementation of the particle filter for the non-linear
% object tracking example 1
% Input:
% yM                  -   the measurements
% DriftMat, DriftVec  -   true drift parameter 
% deltaT              -   true time step
% noiseSTD            -   true noise std
% InitLoc             -   true initial location
% Output:   y_est     -   estimation of the cleaned measurments

% True covariance matrices:
Q = [procStd^2, 0; 0, procStd^2]*deltaT;
R = [noiseSTD(1)^2, 0; 0, noiseSTD(2)^2];

% True state and measurement equations:
procFcn = @(x) [DriftRate1(1,x(:,1))*deltaT+x(:,1), DriftRate2(1,x(:,2))*deltaT+x(:,2)];
measFcn = @(x) [(atan(x(:,1)./x(:,2))), sqrt((x(:,1)).^2+(x(:,2)).^2)];


% Construct the particle filter using the true system equations:
pf                          = robotics.ParticleFilter;
pf.StateEstimationMethod    = 'mean';
pf.ResamplingMethod         = 'systematic';
pf.StateTransitionFcn       = @(pf,x) mvnrnd(procFcn(x), Q);
pf.MeasurementLikelihoodFcn = @(pf,x,y) mvnpdf([y(1), y(2)]-measFcn(x), [0, 0], R);

initialize(pf,1000,InitLoc,eye(size(yM,1)));

x_est        = nan(size(yM));
y_est        = nan(size(yM));
x_est(:,1)   = InitLoc;
y_est(:,1)   = measFcn(x_est(:,1).');

% Calculate the state and measurement estimations using the particle
% filter:
for ii = 2:size(yM,2)
    predict(pf);
    correct(pf,yM(:,ii));
    x_est(:,ii) = getStateEstimate(pf);
    y_est(:,ii) = measFcn(x_est(:,ii).');
end

