function [pc, score, eigenvalues, tsquare, loadings, Xmean] = xASL_im_PCA(dataIn)
%xASL_im_PCA Perform a Principal Component Analysis.
%
% FORMAT:       [pc, score, eigenvalues, tsquare, loadings, Xmean] = xASL_im_PCA(dataIn)
% 
% INPUT:        dataIn      - ...
%
% OUTPUT:       pc          - ...
%               score       - ...
%               eigenvalues - ...
%               tsquare     - ...
%               loadings    - ...
%               Xmean       - ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Perform a Principal Component Analysis.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL


% Data matrix X, with M rows and N columns. Each row m is a set of measurements with N features.
[nM,nN] = size(dataIn);
nR = min(nM-1,nN);     % max possible rank of dataIn

% Average over all measurements
Xmean = mean(dataIn,1);

% Zero-mean for each feature
dataIn = (dataIn - Xmean(ones(nM,1),:));

% U E W^T = X, svd(X) = [U,E,W]
[~,eigenvalues,pc] = svd(dataIn./sqrt(nM-1),0);

% Scores T = XW
score = dataIn*pc;

% W*sqrt of eigenvalues (eigenvectors scales by the variances) = loadings
if nargout >= 5
	loadings = pc*eigenvalues;
end

if nargout < 3, return; end
eigenvalues = diag(eigenvalues).^2;
if (nR<nN)
   eigenvalues = [eigenvalues(1:nR); zeros(nN-nR,1)];
   score(:,nR+1:end) = 0;
end

if nargout < 4, return; end
tmp = sqrt(diag(1./eigenvalues(1:nR)))*score(:,1:nR)';
tsquare = sum(tmp.*tmp)';

return;
