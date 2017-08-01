function [PIC, outputStructure] = pre_ica_cleaning_index(nonCleanedEEG, icaWeight, icaSphere)
% function [PICI outputStructure] = = pre_ica_cleaning_index(nonCleanedEEG, icaWeight, icaSphere)
%
% Calculates an index (Pre ICA Cleaning Index, or PICI) that can be used to compare the quality of ICA cleaning (e.g. artifact rejection)
% using median windowed mutual infromation method. In this method, Mutual information reduction
% index, or MIR is calculated in short ~1 s windows of the original (full, non rejected or otherwise cleaned data) and
% the median of these MIR values is calculated for both an ICA decomposition (often calculated on a
% separate, clean data) and an sphering matrix calculated based on a geometric-median calculated
% covariance matrix.
%
% To use this function you should run the same ICA algorthms (e.g Infomax) on the 'clean' portion of
% the data, and then use provide this ICA matrix along with FULL NON-REJECTED DATA in EEG input
% variable to this function. Please DO NOT send cleaned data in nonCleanedEEG variable to this function as
% this would certainly produce misleading results. The cleaning method that produces the highest
% PICI value is more suitable to be used in later analysis. Bad channels should be removed from
% nonCleanedEEG and all other EEG datasets before using this function.
%
% Also, this function may NOT be used to compare ICA solution from different number of channels (e.g. different bad channels removed) or
% ICA components. You need to always compare the same channel subset (unfortunetaly using PCA to
% reduce dimension will prevent you from this function either because Mutual Information reduction index cannot be used in this case).
%
% Again, if you use this function on cleaned (e.g. artifact rejected) data, the resulting PICI would completely
% meaningless and the 'quality jedgement' produced by the function would be certainly wrong.
%
% Usage Example:
%
% Let originalEEG be a contineous EEG dataset, and we want to compare the cleaning performance of
% methods 1 and 2:
%
% >> cleanEEG1 = clean_by_method_1(originalEEG);
% >> cleanEEG2 = clean_by_method_2(originalEEG);
%
% then we do ICA for data cleaned by these two methods
% >> cleanEEG1WithICA = do_ica(cleanEEG1);
% >> cleanEEG2WithICA = do_ica(cleanEEG2);
%
% then
% >> method1PIC = pre_ica_cleaning_index(nonCleanedEEG, cleanEEG1WithICA.icaweights, cleanEEG1WithICA.icasphere);
% >> method2PICI= pre_ica_cleaning_index(nonCleanedEEG, cleanEEG2WithICA.icaweights, cleanEEG2WithICA.icasphere);
%
% and then look at  method1ICQ - method2ICQ. Even a small difference in PICI (e.g. 1%) could makes a
% significant difference in resulting ICA quality.
%
%   Written by Nima Bigdely-Shamlo, Swartz Center. Copyright 2013, UCSD.

if nargin < 3
	icaSphere = 1;
end;

if nargin < 2
	icaWeight = nonCleanedEEG.icaweights;
	icaSphere = nonCleanedEEG.icasphere;
end;

if nargin < 4
	windowDuration = 1; % 1 second
end;

if isempty(nonCleanedEEG.icachansind)
    fprintf('EEG.icachansind is empty, assuming all channels were used in the ICA.\n');
    nonCleanedEEG.icachansind = 1:size(nonCleanedEEG.data,1);
end;

mutualInformationReductionTimeCourseOfICA = mututal_info_reduction_time_course(nonCleanedEEG, icaWeight * icaSphere , 'windowDuration', windowDuration);
mutualInfoMedianFromICA = median(mutualInformationReductionTimeCourseOfICA);

% now calculate mutual information using robust (geometric-median) sphering
fprintf('Calculating robust sphering matrix using geometric median (takes ~1 minute)...\n');
robustSphere = robust_sphering_matrix(nonCleanedEEG.data);

% get MIR of robust sphering
mutualReductionTimeCourseOfRobustSphering = mututal_info_reduction_time_course(nonCleanedEEG, robustSphere , 'windowDuration', windowDuration);
mutualInfoMedianFromRobustSphering= median(mutualReductionTimeCourseOfRobustSphering);

PIC = 100 * (mutualInfoMedianFromICA - mutualInfoMedianFromRobustSphering) / mutualInfoMedianFromRobustSphering;

if PIC < 0
	qualityJudgement = 'Definately Bad';
elseif PIC > 1.5
	qualityJudgement = 'Could be Good';
elseif PIC < 1.5
	qualityJudgement = 'Might Not be Very Good';
end;

fprintf('PIC = %3.2f%%, which is ''%s''.\n', PIC, qualityJudgement);
fprintf('PIC''s Rule of Thumb: PIC < 0 is bad since robust sphering does better, PIC > 1.5%% could be good, 0 < PIC < 1.5%% might not be very good.\n', PIC, qualityJudgement);

if nargout > 1
	outputStructure.mutualInformationReductionTimeCourseOfICA = mutualInformationReductionTimeCourseOfICA;
	outputStructure.mutualInfoMedianFromRobustSphering = mutualInfoMedianFromRobustSphering;
	outputStructure.robustSphere = robustSphere;
	outputStructure.pic = PIC;
end;