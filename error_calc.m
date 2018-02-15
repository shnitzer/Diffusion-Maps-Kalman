% Compute measurement error
% ***************************************************************@

function [rErr, phiErr, theta1Err, theta2Err] = error_calc(data, tt)
%MEAS_ERROR computes the measurement error and the error in the underlying
% coordinates (using the true measurement functions). 'tt' containes the
% samples which will be included in these error calculations.

Tr_y     = data.Tr_y;
Tr_theta = data.Tr_theta;
if isfield(data,'Est_y')
    Est_y  = data.Est_y;
    rErr   = sqrt(mean((Est_y(2,tt)-Tr_y(2,tt)).^2))/std(Tr_y(2,tt));
    phiErr = sqrt(mean((Est_y(1,tt)-Tr_y(1,tt)).^2))/std(Tr_y(1,tt));
    
    theta1_est = (Est_y(2,tt).*sind(Est_y(1,tt))).';
    theta2_est = (Est_y(2,tt).*cosd(Est_y(1,tt))).';
    
    theta1Err = sqrt(mean((theta1_est-Tr_theta(tt,1)).^2))/std(Tr_theta(tt,1));
    theta2Err = sqrt(mean((theta2_est-Tr_theta(tt,2)).^2))/std(Tr_theta(tt,2));
    
elseif isfield(data,'Est_theta')
    Est_theta = data.Est_theta;
    
    theta1Err = sqrt(mean((Est_theta(tt,1)-Tr_theta(tt,1)).^2))/std(Tr_theta(tt,1));
    theta2Err = sqrt(mean((Est_theta(tt,2)-Tr_theta(tt,2)).^2))/std(Tr_theta(tt,2));
    
    r_est   = sqrt(Est_theta(tt,1).^2 + Est_theta(tt,2).^2).';
    phi_est = atand(Est_theta(tt,1)./Est_theta(tt,2)).';
    
    rErr   = sqrt(mean((r_est-Tr_y(2,tt)).^2))/std(Tr_y(2,tt));
    phiErr = sqrt(mean((phi_est-Tr_y(1,tt)).^2))/std(Tr_y(1,tt));
end

