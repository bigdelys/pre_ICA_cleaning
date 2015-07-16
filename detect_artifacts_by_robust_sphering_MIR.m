function [isFrameAnArtifact rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEG, showPlot, threshold)
% [isFrameAnArtifact rejectionWindows mutualInformationReductionRobustSphering]= detect_artifacts_by_robust_sphering_MIR(EEG)
% isFrameAnArtifact is a boolean vector with the size of number of frames in EEG.data.

if nargin < 2
    showPlot = false;
end;

if nargin < 3
    threshold = 1.9; % 2.1 or 2.5 too
end;

if isempty(EEG.icachansind)
    EEG.icachansind = 1:size(EEG.data,1);
end;

% try
%     [covarianceDirectlyFromGeometricMedian geometricMedian] = robust_covariance_from_direct_geometric_median(EEG.data(EEG.icachansind,1:15:end));
% catch
%     [covarianceDirectlyFromGeometricMedian geometricMedian] = robust_covariance_from_direct_geometric_median(EEG.data(EEG.icachansind,1:50:end));
% end;
% robustSphering = real(inv(sqrtm(covarianceDirectlyFromGeometricMedian)));

robustSphering = robust_sphering_matrix(EEG.data(EEG.icachansind,:));

windowDuration = 1;
mutualInformationReductionRobustSphering = mututal_info_reduction_time_course(EEG, robustSphering, 'windowDuration',1);

%%

medianMIR = median(mutualInformationReductionRobustSphering);
madMIR = median(abs(mutualInformationReductionRobustSphering - medianMIR));

isFrameAnArtifact = mutualInformationReductionRobustSphering <= medianMIR - threshold * madMIR;

% expnd the rejection part to enclose any frame contributed to the MIR window asociated with
% rejected section.
isFrameAnArtifact = logical(moving_average(windowDuration * EEG.srate, double(isFrameAnArtifact)));

%%

% create a format compatible with eegplot artifact plotting

raisingEdge = find(diff([0 isFrameAnArtifact']) > 0);
fallingEdge = find(diff([isFrameAnArtifact' 0]) < 0);

rejectionWindows = [];
for i=1:length(fallingEdge)
    rejectionWinowStart = raisingEdge(i);
    rejectionWinowEnd = fallingEdge(i);
    rejectionWindows =cat(1, rejectionWindows, [rejectionWinowStart rejectionWinowEnd]);
end;

%%
if showPlot
    eegplot(robustSphering * EEG.data(EEG.icachansind,:), 'srate', EEG.srate, 'winrej', rejectionWindows,  'events' , EEG.event);
end;