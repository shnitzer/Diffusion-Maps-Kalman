% Estimate local covariance matrices and the modified mahalanobis distance
% ***************************************************************@

function [mahDist] = modified_mahalanobis(yM)
%MODIFIED_MAHALANOBIS calculates estimated covariances and the modified
% mahalanobis distance for the input data yM (size: variables x samples)

%% Configuration
ncov     = 15; % size of neighborhood for covariance 
finalDim = 2;  % final data dimension

%% Covariance estimation

inv_c = zeros(size(yM,1), size(yM,1), size(yM,2)); 
for i = 1+ncov:length(yM)-ncov
    
    % Estimate covariance in short time windows
    win = yM(:, i-ncov:i+ncov-1);
    c   = cov(win');
    
    % Denoise via projection on "known" # of dimensions
    [U, S, V]    = svd(c);
    inv_c(:,:,i) = V(:,1:finalDim) / (S(1:finalDim,1:finalDim)) * U(:,1:finalDim)';
    
end

% Complete missing covariance matrices (beginning and end) by duplication
for i = 1:ncov
    inv_c(:,:,i) = inv_c(:,:,1+ncov);
end
for i = (length(yM)-ncov+1):length(yM)
    inv_c(:,:,i) = inv_c(:,:,length(yM)-ncov);
end

%% Mahalanobis distance calculation

data    = yM.';
mahDist = zeros(size(yM,2));

for i = 1:size(yM,2)
    mahDist(:,i) = sum((bsxfun(@minus,data,data(i,:))*inv_c(:,:,i)).*bsxfun(@minus,data,data(i,:)),2);
end
mahDist = (mahDist + mahDist.');

end