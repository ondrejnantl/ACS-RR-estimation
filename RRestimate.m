function breathVAL = RRestimate(fs, PPG, ECG)
% Estimate of respiratory rate [breaths per minute] from a one minute 
% section of signal
% -------------------------------------------------------------------------
% This function estimates the respiratory rate as dominant frequency within
% respiratory frequency range from spectrum created by multiplication  
% the spectra of approximations of respiratory signals from ECG and PPG
% 
% Respiratory signal is gained by approximating baseline wander of biosignal 
% by interpolating consequent peaks of biosignal by cubic spline
% -------------------------------------------------------------------------
% Input: 
% fs - sampling rate of biosignals as integer, must be the same for both 
%      input signals
% ECG - input ECG signal as vector with length Nx1, where N is the number 
%       of samples
% PPG - input PPG signal as vector with length Nx1, where N is the number
%       of samples
% Output:
% breathVAL - estimated average respiratory rate as double
% -------------------------------------------------------------------------
% Author: Ondřej Nantl
% =========================================================================
% calculating length of signal and Hamming window for spectra
len = length(PPG);
wind = hamming(len,'periodic');

% subtracting baseline offset
ECG = ECG - mean(ECG);
PPG = PPG - mean(PPG);

% estimating R peaks distance for interpolation - PPG
PPGSpec = abs(fft(wind.*PPG))./len; % calculating spectrum of PPG 
w = linspace(0,fs-fs/len,len);
idL = find(w >= 0.5,1,'first');
idR = find(w <= 5,1,'last');
[~,HRID] = findpeaks(PPGSpec(idL:idR),'SortStr','descend','NPeaks',1);
HR = w(HRID+idL-1); % finding the peak - approximation of average heart rate
peakDist = round(fs/HR)-15; % subtracting 15 samples from average R peak distance

% calculating breathing frequency axis
wR = 60.*linspace(0,fs-fs/(len),len);

% estimating baseline wander by calculating peak envelope of ECG signal
bwECG = envelope(ECG,peakDist,'peak');
% estimating its spectrum
bwSpecECG = (abs(fft(wind.*bwECG)))./len;

% estimating baseline wander by calculating peak envelope of PPG signal
bwPPG = envelope(PPG,peakDist,'peak');
% estimating its spectrum
bwSpecPPG = (abs(fft(wind.*bwPPG)))./len;

% multipling of the obtained spectra
bwSpec = bwSpecPPG.*bwSpecECG./len;

% finding peak in breathing frequency range and obtaining breathing
% frequency value
idL = find(wR >= 0.067*60,1,'first');
idR = find(wR <= 1.08*60,1,'last');
[~,breathVALID] = findpeaks(bwSpec(idL:idR),'SortStr','descend','NPeaks',1);
breathVAL = wR(breathVALID+idL-1);
end