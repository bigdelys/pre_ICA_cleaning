function [badChannelsBothMethods channelBadWindowRatio badChannelsByCorrelation badChannelsBasedOnAmplitude normalizeAmplitudeDistanceToMean ] = eeg_detect_bad_channels(EEG, ratioTimeBadTolerated, minAbsCorrelationThresholdAccepted)
% [badChannelsBothMethods channelBadWindowRatio badChannelsByCorrelation badChannelsBasedOnAmplitude normalizeAmplitudeDistanceToMean ] = eeg_detect_bad_channels(EEG, ratioTimeBadTolerated, minAbsCorrelationThresholdAccepted)
% Detects bad channels and places them in badChannels output variable. the last three variables are optional.
% Usage:
% badChannels = eeg_detect_bad_channels(EEG);
%
% This function uses methods for detecting bad channels:
% method 1: low correlation with other channels in 1 s windows
% method 2: too low or high amplitude.

%% detect bad channel by their correlation in 1 second windows

fprintf('Detecting channels with low correlations to others in 1 s windows...\n');

timeWindowForCalculatingCorrelation = 1; % in seconds

if nargin<2
    ratioTimeBadTolerated = 0.01;
end;

if nargin<3
    minAbsCorrelationThresholdAccepted = 0.4;
end;

data = eeg_getdatact(EEG);

% make sure channel data is 2D, conatenates epochs if present
data = data(:,:);

numberOfChannels = size(data,1);
numberOfFramesInTimeWindow = round(timeWindowForCalculatingCorrelation * EEG.srate);
numberOfTimeWindows = floor(size(data, 2) / numberOfFramesInTimeWindow);
halfNumberOfFrames = round(numberOfFramesInTimeWindow / 2);


rejectedChannels = zeros(numberOfTimeWindows, numberOfChannels);

for i=2:(numberOfTimeWindows-2) % ignore last two time windows, so we dont go out of range by chance
    eegPortion = data(:, (i*(numberOfFramesInTimeWindow-1) - halfNumberOfFrames): (i*(numberOfFramesInTimeWindow-1) + halfNumberOfFrames));
    correlationBetweenChannelsInWindow = corrcoef(eegPortion');
    correlationBetweenChannelsInWindow = correlationBetweenChannelsInWindow - diag(diag(correlationBetweenChannelsInWindow));
    absCorrelationBetweenChannelsInWindow = abs(correlationBetweenChannelsInWindow);
    % maxChannelAbsCorrelationInWindow = max(absCorrelationBetweenChannelsInWindow);
    maxChannelAbsCorrelationInWindow = quantile(absCorrelationBetweenChannelsInWindow,0.98);
    rejectedChannels(i,:)  = maxChannelAbsCorrelationInWindow < minAbsCorrelationThresholdAccepted;
end;

channelBadWindowRatio = mean(rejectedChannels,1);
badChannelsByCorrelation = find(channelBadWindowRatio > ratioTimeBadTolerated);

%% detect bad channels by having an unusually highor low amplitude (quantified by Standard Deviation)

fprintf('Detecting channels with unusually high or low (> 3 std) amplitude...\n');
channelStandardDeviation =  0.7413 *iqr(data');

% Estimate the std of the distribution of channel amplitudes (estimated by their standrad deviations)
% using the inter-quantile distance which is more robust to outliers
estimateStandardDeviationOfChannelAmplitude =  0.7413 * iqr(channelStandardDeviation);

% median used instead of mean to be more robust to bad channel STDs
normalizeAmplitudeDistanceToMean = (channelStandardDeviation-median(channelStandardDeviation)) / estimateStandardDeviationOfChannelAmplitude;

% channels with amplitude far from 3 std. from mean channel powers are unusuall (bad).
badChannelsBasedOnAmplitude = find(abs(normalizeAmplitudeDistanceToMean) > 5);


% combine bad channels detected from both methods
badChannelsBothMethods = union(badChannelsByCorrelation, badChannelsBasedOnAmplitude);