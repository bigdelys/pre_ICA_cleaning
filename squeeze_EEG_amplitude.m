function squeezedEEG = squeeze_EEG_amplitude(EEG, alpha)
% squeezedData = squeeze_EEG_amplitude(EEG, alpha)

if nargin < 2
    alpha = 1;
end;

frameAmplitude = vec(sum(EEG.data .^2, 1).^0.5);
robustStd = std_from_mad(frameAmplitude);

squeezedFrameAmplitude= 1./(1+exp(-alpha * frameAmplitude / (4 * robustStd)));

factor  = squeezedFrameAmplitude ./ frameAmplitude;

squeezedEEG = EEG;
squeezedEEG.dat = bsxfun(@times, EEG.data, factor'); 