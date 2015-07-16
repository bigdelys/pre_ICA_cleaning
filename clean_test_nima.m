function [signal,sample_mask] = clean_test_nima(signal,maxbad,stddev)
% Substitute windows with high-power activation by repaired content.
% [Signal,Mask] = clean_test(Signal,MaxBadChannels)
%
% In:
%   Signal         : Continuous data set, assumed to be appropriately high-passed (e.g. >1Hz or
%                    0.5Hz - 2.0Hz transition band)
%
%   MaxBadChannels : The maximum number or fraction of bad channels that a retained window may still
%                    contain (more than this and it is removed). Reasonable range is 0.05 (very clean
%                    output) to 0.3 (very lax cleaning of only coarse artifacts). Default: 0.15.
%
%   StandardDevCutoff: StdDev cutoff for repairs. Data segments whose variance is beyond this cutoff 
%                      from the distribution of variance across the recording are considered missing data.
%                      The reasonable range is 3 (fairly aggressive), 4 (fairly safe), or 5 (very safe,
%                      likely not harming any EEG). Default: 3.
%
% Out:
%   Signal : data set with bad time periods repaired.
%
%   Mask   : mask that is 1 for untouched samples and 0 for repaired samples
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-11-23

if ~exist('maxbad','var') || isempty(maxbad) maxbad = 0.15; end
if ~exist('stddev','var') || isempty(stddev) stddev = 3; end

% first determine the breakage mask
%[dummy,sample_mask] = clean_windows(signal,maxbad);
% instead we use Nima's amplitude-based
isFrameAnArtifact = eeg_clean_data_by_probability_robust(signal, false);
sample_mask = ~isFrameAnArtifact;

% generate a repaired version of the data set
repaired = clean_asr(signal,stddev);

% substitute repaired content into the original signal
signal.data(:,~sample_mask) = repaired.data(:,~sample_mask);

