function [isFrameAnArtifact rejectionWindows]= eeg_clean_data_by_probability_robust(EEG, showPlot, minRemoveRatio, ignoreBadChannels)
%    [isFrameAnArtifact rejectionWindows]= eeg_clean_data_by_probability_robust(EEG, showPlot, minRemoveRatio)
%    minRemoveRatio: at least this much of data will be marked as artifact.
%    Usage:
%		
%		[isFrameAnArtifact rejectionWindows]= eeg_clean_data_by_probability_robust(EEG, false);
% 		if isempty(EEG.icachansind) % making sure that EEG.icachansind which contains channels to be used for ICA is not empty
% 			EEG.icachansind  =1:size(EEG.data, 1);
% 		end;
% 
% 		% assuming your data is 2D, and not epoched
% 		cleanData = EEG.data(EEG.icachansind, ~isFrameAnArtifact);
% 
% 		% do your ICA here and get wts and sph matrices, for example
% 		% [wts sph] = binica(cleanData, ...);
% 
% 		iWts = pinv(wts*sph);
% 		scaling = repmat(sqrt(mean(iWts.^2))', [1 size(wts,2)]);
% 		wts = wts.*scaling;
% 
% 		EEG.icawinv = pinv(wts*sph);
% 		EEG.icasphere = sph;
% 		EEG.icaweights = wts;
% 		EEG.icaact = [];
% 		EEG = eeg_checkset(EEG);

clear cleanData

if nargin < 2
    showPlot = false;
end;

if nargin < 3
    minRemoveRatio = 0.05; % at least mark 5% of the least likely parts as artifgacts.
end;

if nargin < 4
    ignoreBadChannels = true;
end;

% ignore bad channels

goodChanneId = 1:EEG.nbchan;

if ignoreBadChannels    
    if isempty(EEG.icaweights)
        badChannel = eeg_detect_bad_channels(EEG);
        goodChanneId(badChannel) = [];
    end;
    
    % also ignore channels not present in EEG.icachansind
    if ~isempty(EEG.icachansind)
        goodChanneId = intersect(goodChanneId, EEG.icachansind);
    end;
end;

data = double(eeg_getdatact(EEG));
data = data(goodChanneId,:);

if isempty(EEG.icaweights)
   % covarianceDirectlyFromGeometricMedian = robust_covariance_from_direct_geometric_median(data(:,1:round(EEG.srate / 4):end));
   % sphere =inv(real(sqrtm(covarianceDirectlyFromGeometricMedian)));
     sphere = robust_sphering_matrix(data);
    
    data = real(sphere*data);      % actually decorrelate the electrode signals
else
    data = real((EEG.icaweights * EEG.icasphere) * data);
end;


%%
logLikelihood = zeros(size(data));
for i=1:size(data,1)
    if mod(round(100*i/size(data,1)),10) == 0
        fprintf('%d%%..', round(100*i/size(data,1)));
    end;
    triedRank = tiedrank(data(i,:)) / size(data, 2);
    twoSidedPvalue = min(triedRank, 1 - triedRank);
    logLikelihood(i,:) = -log(twoSidedPvalue);
end;

fprintf('\n');
meanLogLikelihood= mean(logLikelihood, 1);

windowTimeLenght = 200;%in ms
windowFrameLenght = round((EEG.srate * windowTimeLenght/1000));
windowFrame = round((-windowFrameLenght/2):(windowFrameLenght/2));
smoothMeanLogLikelihood =  moving_average(windowFrameLenght, meanLogLikelihood)';


isArtifactWindowCenter = find(smoothMeanLogLikelihood > 2.1 | smoothMeanLogLikelihood > quantile(smoothMeanLogLikelihood, 1-minRemoveRatio));

% add two sides on the window
artifactFrames = repmat(windowFrame, length(isArtifactWindowCenter), 1) + repmat(isArtifactWindowCenter, 1, length(windowFrame));
artifactFrames = max(artifactFrames, 1);
artifactFrames = min(artifactFrames, length(smoothMeanLogLikelihood));
artifactFrames = unique(artifactFrames(:));

isFrameAnArtifact = false(1, length(smoothMeanLogLikelihood));
isFrameAnArtifact(artifactFrames) = true;

raisingEdge = find(diff([0 isFrameAnArtifact]) > 0);
fallingEdge = find(diff([isFrameAnArtifact 0]) < 0);

rejectionWindows = [];
for i=1:length(fallingEdge)
    rejectionWinowStart = raisingEdge(i);
    rejectionWinowEnd = fallingEdge(i);
    rejectionWindows =cat(1, rejectionWindows, [rejectionWinowStart rejectionWinowEnd]);
end;

%%

if showPlot
    eegplot(data, 'srate', EEG.srate, 'winrej', rejectionWindows,  'events', EEG.event);
end;