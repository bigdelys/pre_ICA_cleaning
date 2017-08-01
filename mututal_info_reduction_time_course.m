function mutualInformationReduction = mututal_info_reduction_time_course(EEG, unmixingMatrix, varargin)
% mutualInformationReduction = mututal_info_reduction_time_course(EEG, unmixingMatrix, {name, value pairs})

inputOption = finputcheck(varargin, ...
    { 'windowDuration'     'real'   [0 Inf]       2; ... % in seconds
    'interestFrame'         'real'   []    [1 size(EEG.data,2)] ; ...
    'progressbar'         'boolean'  [false true]    false ; ...
    'windowOverlapRatio'   'real'  [0 1]   0.5 });

windowLength = inputOption.windowDuration * EEG.srate ; % in frames
windowStep = round(windowLength * (1-inputOption.windowOverlapRatio)); % in frames

currentWindowStart = inputOption.interestFrame(1);
i = 1;
mutualInformationReductionSum = zeros(inputOption.interestFrame(2) - inputOption.interestFrame(1) + 1, 1);
numberOfWindowsForFrame = zeros(inputOption.interestFrame(2) - inputOption.interestFrame(1) + 1, 1);

numberOfWindows = ceil(length(numberOfWindowsForFrame) / windowStep);
if inputOption.progressbar
progress('init','');
end;

while (currentWindowStart + windowLength -1) <= max(inputOption.interestFrame(2))
    
	if inputOption.progressbar && mod(i,10) == 0
	progress(i/numberOfWindows, sprintf('progress: %d %%', round(100*i/numberOfWindows)));
	end;
	
    windowFrame = currentWindowStart:(currentWindowStart + windowLength -1);
    windowMutualInfoReduction = getMIR(unmixingMatrix, EEG.data(EEG.icachansind,windowFrame));
    
    windowFrameForMir = windowFrame -  inputOption.interestFrame(1) +1;
    mutualInformationReductionSum(windowFrameForMir) = mutualInformationReductionSum(windowFrameForMir) + windowMutualInfoReduction;
    numberOfWindowsForFrame(windowFrameForMir) = numberOfWindowsForFrame(windowFrameForMir) + 1;
    
    mutualInformationReductionForWindow(i) = windowMutualInfoReduction;
    
    currentWindowStart = currentWindowStart + windowStep;
    i = i +1;
end;

if inputOption.progressbar
pause(0.1)
progress('close');
end;

mutualInformationReduction =real(mutualInformationReductionSum ./ numberOfWindowsForFrame);
mutualInformationReduction(isnan(mutualInformationReduction)) = 0;

end