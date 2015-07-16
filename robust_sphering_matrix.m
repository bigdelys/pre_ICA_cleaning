function [robustSphering, mixing, covarianceMatrix] = robust_sphering_matrix(X)
% [robustSphering mixing] = robust_sphering_matrix(X);
% X is channel x times data, e.g. EEG.data

[C,S] = size(X);
X = X';
blocksize = 10;
blocksize = max(blocksize,ceil((C*C*S*8*3*2)/hlp_memfree));

% calculate the sample covariance matrices U (averaged in blocks of blocksize successive samples)
U = zeros(length(1:blocksize:S),C*C);
for k=1:blocksize
    range = min(S,k:blocksize:(S+k-1));
    U = U + reshape(bsxfun(@times,reshape(X(range,:),[],1,C),reshape(X(range,:),[],C,1)),size(U));
end

% get the mixing matrix M
covarianceMatrix = real(reshape(block_geometric_median(U/blocksize),C,C));
mixing = sqrtm(covarianceMatrix);
robustSphering = inv(mixing);