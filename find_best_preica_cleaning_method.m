function result = find_best_preica_cleaning_method(EEG, EEGHandCleaned, fileToSave, methodNumbersToRun)
% find the best method for rejecting parts of the data before performing
% ICA

handCleanedExists = nargin > 1 && ~isempty(EEGHandCleaned);
if nargin<3
    fileToSave = [tempname '.mat'];
    fprintf(['Result will be saved to file ' fileToSave]);
end;

if nargin<4
    methodNumbersToRun = [];
end;

if isempty(EEG.icachansind)
    EEG.icachansind = 1:size(EEG.data,1);
end;

if handCleanedExists && isempty(EEGHandCleaned.icachansind)
    EEGHandCleaned.icachansind = 1:size(EEG.data,1);
end;

% EEG.icachansind should be similar between the two
if handCleanedExists &&  ~isequal(EEG.icachansind, EEGHandCleaned.icachansind)
    error('icachansind differs!');
end;

EEG.data = double(EEG.data);
if handCleanedExists
EEGHandCleaned.data = double(EEGHandCleaned.data);
end;

for methodCounter = 1:(26+uint8(handCleanedExists)) % the last method uses hand-cleaned data
    
    if ~isempty(methodNumbersToRun) && ~ismember(methodCounter, methodNumbersToRun)
        continue;
    end;
    
    tic;
    result.method(methodCounter).ratioRemoved = 0;
    w = nan;
    s = nan;
    switch methodCounter
        case 1
            result.method(methodCounter).label = 'Arno''s spectrum thresholding';
            EEGRejected = pop_rejcont(EEG, 'elecrange', EEG.icachansind ,'freqlimit',[20 40] ,'threshold',10,'epochlength',0.5,'contiguous',4,'addlength',0.25, 'taper', 'hamming');
            
            try
                [w, s] = cudaica(double(EEGRejected.data(EEG.icachansind,:)), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = (size(EEG.data,2) - size(EEGRejected.data,2)) / size(EEG.data,2);
        case 2
            result.method(methodCounter).label = 'Nima''s amplitude-based frame rejection';
            [cleanData bad_frames bad_frames_marked_zero]= eeg_badframes(EEG);
            try
                [w, s] = cudaica(cleanData(EEG.icachansind,:), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = sum(~bad_frames_marked_zero) / length(bad_frames_marked_zero);
        case 3
            result.method(methodCounter).label = 'Nima''s frame and MIR rejection';
            [cleanData bad_frames bad_frames_marked_zero]= eeg_badframes(EEG);
            [isFrameAnArtifact rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEG, false, 2.1);
            
            isFrameAnArtifact = isFrameAnArtifact | ~bad_frames_marked_zero';
            try
                [w, s] = cudaica(EEG.data(EEG.icachansind, ~isFrameAnArtifact), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = sum(isFrameAnArtifact) / length(isFrameAnArtifact);
        case 4
            result.method(methodCounter).label = 'Nima''s amplitude-based window rejection';
            
            isFrameAnArtifact = eeg_clean_data_by_probability_robust(EEG, false);
            try
                [w, s] = cudaica(EEG.data(EEG.icachansind, ~isFrameAnArtifact), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = sum(isFrameAnArtifact) / length(isFrameAnArtifact);
            
        case 5
            result.method(methodCounter).label = 'Nima''s amplitude and MIR-based window rejection';
            
            isFrameAnArtifact = eeg_clean_data_by_probability_robust(EEG, false);
            [isFrameAnArtifactByMIR rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEG, false, 2.1);
            isFrameAnArtifact = isFrameAnArtifact(:) | isFrameAnArtifactByMIR(:);
            
            try
                [w, s] = cudaica(EEG.data(EEG.icachansind, ~isFrameAnArtifact), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = sum(isFrameAnArtifact) / length(isFrameAnArtifact);
        case 6
            result.method(methodCounter).label = 'Christian''s amplitude-based window rejection';
            [EEGcleaned, Mask]= clean_windows(EEG);
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = sum(~Mask) / length(Mask);
        case 7
            result.method(methodCounter).label = 'Christian''s amplitude and MIR-based window rejection';
            [EEGcleaned, Mask]= clean_windows(EEG);
            [isFrameAnArtifactByMIR rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEG, false, 2.1);
            
            Mask = Mask(:) & ~isFrameAnArtifactByMIR(:);
            try
                [w, s] = cudaica(EEG.data(EEG.icachansind, Mask), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = sum(~Mask) / length(Mask);
        case 8
            result.method(methodCounter).label = 'Nima-Christian Combo then MIR-based rejection'; % MIR rejectiuon is still on the original here.
            [EEGcleaned, Mask]= clean_test_nima(EEG);
            [isFrameAnArtifactByMIR rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEG, false, 2.1);
            
            Mask = ~isFrameAnArtifactByMIR(:);
            try
                [w, s] = cudaica(EEGcleaned.data(EEG.icachansind, Mask), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = sum(~Mask) / length(Mask);
        case 9
            result.method(methodCounter).label = 'Regular Sphering';
            s = real(inv(sqrtm(cov(EEG.data(EEG.icachansind,:)'))));
            w = eye(size(s));
        case 10
            result.method(methodCounter).label = 'Robust Sphering';
            
            %[covarianceDirectlyFromGeometricMedian geometricMedian] = robust_covariance_from_diect_geometric_median(EEG.data(EEG.icachansind,1:10:end));
            %s = real(inv(sqrtm(covarianceDirectlyFromGeometricMedian)));
            s = robust_sphering_matrix(EEG.data(EEG.icachansind,:));
            w = eye(size(s));
            
            %             case 11
            %                 result.method(methodCounter).label = 'Hand-Cleaned with Robust Sphering';
            %
            %                 [covarianceDirectlyFromGeometricMedian geometricMedian] = robust_covariance_from_diect_geometric_median(EEGHandCleaned.data(EEG.icachansind,1:10:end));
            %                 robustSphering = real(inv(sqrtm(covarianceDirectlyFromGeometricMedian)));
            %                 spheredData = robustSphering * double(EEGHandCleaned.data(EEG.icachansind,:));
            %                 spheredData = bsxfun(@minus, spheredData, mean(spheredData,2));
            %                 try
            %                     [w, sIca] = cudaica(spheredData, 'extended', 3, 'sphering', 'off');
            %                     clear spheredData;
            %                     s = sIca * robustSphering;
            %                 catch, end;
            %                 result.method(methodCounter).ratioRemoved = size(EEGHandCleaned.data,2) / size(EEG.data,2);
        case 11
            result.method(methodCounter).label = 'Arno''s spectrum thresholding with MIR';
            EEGRejected = pop_rejcont(EEG, 'elecrange', EEG.icachansind ,'freqlimit',[20 40] ,'threshold',10,'epochlength',0.5,'contiguous',4,'addlength',0.25, 'taper', 'hamming');
            
            isFrameAnArtifactByMIR = detect_artifacts_by_robust_sphering_MIR(EEGRejected, false, 2.1);
            Mask = ~isFrameAnArtifactByMIR(:);
            
            try
                [w, s] = cudaica(double(EEGRejected.data(EEG.icachansind, Mask)), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = (size(EEG.data,2) - size(EEGRejected.data,2)) / size(EEG.data,2);
        case 12
            result.method(methodCounter).label = 'Original';
            try
                [w, s] = cudaica(EEG.data(EEG.icachansind,:), 'extended', 3);
            catch, end;
        case 13
            result.method(methodCounter).label = 'MIR rejection';
            [isFrameAnArtifact ]= detect_artifacts_by_robust_sphering_MIR(EEG, false, 2.1);
            try
                [w, s] = cudaica(EEG.data(EEG.icachansind, ~isFrameAnArtifact), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = sum(isFrameAnArtifact) / length(isFrameAnArtifact);
        case 14
            result.method(methodCounter).label = 'Christian''s Combo';
            [EEGcleaned, Mask]= clean_test(EEG);
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = 0;
        case 15
            result.method(methodCounter).label = 'Christian''s Combo then MIR-based rejection'; % MIR rejectiuon is still on the original here.
            [EEGcleaned, Mask]= clean_test(EEG);
            [isFrameAnArtifactByMIR rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEG, false, 2.1);
            
            Mask = ~isFrameAnArtifactByMIR(:);
            try
                [w, s] = cudaica(EEGcleaned.data(EEG.icachansind, Mask), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = sum(~Mask) / length(Mask);
        case 16
            result.method(methodCounter).label = 'ASR with 3 std';
            EEGcleaned = clean_asr(EEG, 3);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
        case 17
            result.method(methodCounter).label = 'ASR then MIR-based rejection on changed';
            EEGcleaned = clean_asr(EEG, 3);
            
            [isFrameAnArtifactByMIR rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEGcleaned, false, 2.1);
            
            Mask = ~isFrameAnArtifactByMIR(:);
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, Mask), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
        case 18
            result.method(methodCounter).label = 'Christian''s Combo then MIR-based rejection on changed';
            [EEGcleaned, Mask]= clean_test(EEG);
            [isFrameAnArtifactByMIR rejectionWindows mutualInformationReductionRobustSphering, robustSphering]= detect_artifacts_by_robust_sphering_MIR(EEGcleaned, false, 2.1);
            
            Mask = ~isFrameAnArtifactByMIR(:);
            try
                [w, s] = cudaica(EEGcleaned.data(EEG.icachansind, Mask), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = sum(~Mask) / length(Mask);
            
            % we need to add a hristian Combo but using Nima's amplitude instead here, also a
            % version with MIR
        case 19
            result.method(methodCounter).label = 'Nima-Christian Combo';
            [EEGcleaned, Mask]= clean_test_nima(EEG);
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = 0;
        case 20
            result.method(methodCounter).label = 'ASR-15';
            EEGcleaned = clean_asr(EEG, 15);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
                  result.method(methodCounter).ratioRemoved = 0;
        case 21
            result.method(methodCounter).label = 'ASR with 50 std';
            EEGcleaned = clean_asr(EEG, 50);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
        case 22
            result.method(methodCounter).label = 'ASR with 70 std';
            EEGcleaned = clean_asr(EEG, 70);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
        case 23
            result.method(methodCounter).label = 'Squeeze alpha 1';
            EEGcleaned = squeeze_EEG_amplitude(EEG, 1);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
        case 24
            result.method(methodCounter).label = 'Squeeze alpha 1.5';
            EEGcleaned = squeeze_EEG_amplitude(EEG, 1.5);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
            
        case 25
            result.method(methodCounter).label = 'Squeeze alpha 0.7';
            EEGcleaned = squeeze_EEG_amplitude(EEG, 0.7);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
            
        case 26
            result.method(methodCounter).label = 'Squeeze alpha 2';
            EEGcleaned = squeeze_EEG_amplitude(EEG, 2);
            
            try
                [w, s] = cudaica(EEGcleaned.data(EEGcleaned.icachansind, :), 'extended', 3);
            catch, end;
            
            result.method(methodCounter).ratioRemoved = mean(any(EEGcleaned.data(EEGcleaned.icachansind, :) ~= EEG.data(EEGcleaned.icachansind, :)));
            
            
        case 27
            result.method(methodCounter).label = 'Hand-Cleaned';
            try
                [w, s] = cudaica(EEGHandCleaned.data(EEGHandCleaned.icachansind,:), 'extended', 3);
            catch, end;
            result.method(methodCounter).ratioRemoved = size(EEGHandCleaned.data,2) / size(EEG.data,2);
    end;
    
    result.method(methodCounter).w = w;
    result.method(methodCounter).s = s;
    result.method(methodCounter).chanlocs = EEG.chanlocs(EEG.icachansind);
    
    if any(isnan(vec(w*s)))
        mutualInformationReductionTimeCourse = 0;
    else
        mutualInformationReductionTimeCourse = mututal_info_reduction_time_course(EEG, result.method(methodCounter).w * result.method(methodCounter).s , 'windowDuration',1);
    end;
    
    result.method(methodCounter).mutualInformationReductionMedian = median(mutualInformationReductionTimeCourse);
    result.method(methodCounter).timeIncludingICA = toc;
    
    save(fileToSave, 'result');
    fprintf(['Result saved to file ' fileToSave]);
end;

