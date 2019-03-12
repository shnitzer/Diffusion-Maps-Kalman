% Construct the diffusion maps kernel and eigenvectors
% ***************************************************************@

function [psi, lambda] = diffusion_maps(dis, DMdim)
%DIFFUSION_MAPS constructs the diffusion maps kernel based on the distance
% matrix 'dis' and returns the (non-trivial) eigenvectors 'psi' and 
% eigenvalues 'lambda' of dimension DMdim.

%% Construct diffusion maps

ep = 2*median(dis(:));  % kernel scale - may require some tuning for different data
W  = exp(-dis/ep);      % distance kernel
D  = diag(1./sum(W,2)); % kernel normalization term
A  = D*W;               % normalized diffusion maps kernel

[psi, mu]        = eigs(A,DMdim+1);
[muSrtd, srtInd] = sort(diag(mu),'descend');

lambda = (2/ep)*log(muSrtd(2:end));
psi    = psi(:,srtInd(2:end));

