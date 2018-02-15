function [ theta_est ] = particle_filter( yM, DriftMat, DriftVec, deltaT, noiseSTD, InitLoc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

DriftMat = DriftMat*deltaT + eye(size(yM,1));
DriftVec = DriftVec*deltaT;

Q = [sqrt(2), 0; 0, sqrt(2)]*deltaT;
R = [noiseSTD(1)^2, 0; 0, noiseSTD(2)^2];

procFcn = @(theta) [theta(:,1), theta(:,2)]*DriftMat + DriftVec.';
measFcn = @(theta) [(atand(theta(:,1)./theta(:,2))), sqrt((theta(:,1)).^2+(theta(:,2)).^2)];

pf                          = robotics.ParticleFilter;
pf.StateEstimationMethod    = 'mean';
pf.ResamplingMethod         = 'systematic';
pf.StateTransitionFcn       = @(pf,theta) mvnrnd(procFcn(theta), Q);
pf.MeasurementLikelihoodFcn = @(pf,theta,y) mvnpdf([y(1), y(2)]-measFcn(theta), [0, 0], R);

initialize(pf,1000,InitLoc,eye(size(yM,1)));

theta_est        = nan(size(yM));
theta_est(:,1)   = InitLoc;

for ii = 2:size(yM,2)
    predict(pf);
    correct(pf,yM(:,ii));
    theta_est(:,ii) = getStateEstimate(pf);
end

