function [covarianceDirectlyFromGeometricMedian geometricMedian convergenceHistory] = robust_covariance_from_direct_geometric_median(data, useAllData)
% [covarianceDirectlyFromGeometricMedian geometricMedian]= robust_covariance_from_diect_geometric_median(data)
% data should be channels x time points.

if nargin < 2
    useAllData = false;
end;

% calculate the minimum number of points required to calculate the covariance
minNumberOfPoints = min((size(data,1) ^ 2) * 30, size(data,2));
covarianceDirectlyFromGeometricMedianInitial = [];

if  size(data,2) > minNumberOfPoints
    % first run it on a much shorter data    
    [covarianceDirectlyFromGeometricMedian geometricMedian convergenceHistory] = robust_covariance_from_direct_geometric_median(data(:,round(linspace(1, size(data,2), minNumberOfPoints))), true);
    if ~useAllData
        return;
    end;
else
    geometricMedian = mean(data');      
end;

try   
    geometricMedian = geometric_median(data', 'initialGuess', geometricMedian);
    
    % remove the medial from all
    data = bsxfun(@minus, data, geometricMedian');
    pointCovariance = zeros(size(data,2), size(data,1)^2);
    for i=1:size(data,2)
        pointCovariance(i,:) = vec(data(:,i) * data(:,i)');
    end;
    
    if isempty(covarianceDirectlyFromGeometricMedianInitial)
        covarianceDirectlyFromGeometricMedianInitial = mean(pointCovariance);
    end;
    
    [vectorizedCovariance convergenceHistory]= geometric_median(pointCovariance, 'initialGuess', covarianceDirectlyFromGeometricMedianInitial(:)');
    covarianceDirectlyFromGeometricMedian = reshape(vectorizedCovariance, [size(data,1) size(data,1)]);
catch
    fprintf('There was a issue (likely not enough memory) when trying to use all the data, falling back to use a sufficiently large number of samples.\n');
end;