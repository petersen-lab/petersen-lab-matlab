function [isOsc, oscFreq, maxPower] = detectLFPOsc(lfp, freqRange, options)
% [isOsc, oscFreq, maxPower] = detectLFPOsc(lfp, freqRange, <samplingFreq>)
%
% Detect if the continuous signal is oscillatory within a given frequency
% range.
%
% Args:
%   lfp (numeric, required, positional): a shape-(M, N) numeric array of
%     LFP recordings with rows corresponding to individual probe channels
%     and columns corresponding to individual samples.
%   freqRange (numeric, required, positional): a shape-(1, 2) numeric array
%     of frequency range of interest.
%   samplingFreq (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar with sampling frequency (default=500).
%
% Returns:
%   isOsc (logical): a shape-(M, 1) logical array indicating if recorded
%     LFPs corresponding to individual probe channels are oscillatory.
%   oscFreq (numeric): a shape-(M, 1) numeric array of oscillation
%     frequencies of individual LFP channels.
%   maxPower (numeric): a shape-(M, 1) numeric array of power amplitudes of
%     the oscillatory power peak for individual LFP channels.
%
% Authors:
%   Martynas Dervinis

arguments
  lfp (:,:) {mustBeNumeric}
  freqRange (1,2) {mustBePositive}
  options.samplingFreq (1,1) {mustBePositive} = 500
end

nCh = size(lfp,1);
isOsc = false(nCh,1);
oscFreq = zeros(nCh,1);
maxPower = zeros(nCh,1);
for iCh = 1:nCh

  % Compute the power spectrum of the signal
  inds = 1:(floor(numel(lfp(iCh,:))/2)*2);
  N = length(lfp(iCh,inds));
  LFP = fft(lfp(iCh,inds) - mean(lfp(iCh,inds)));
  P2 = abs(LFP/N);
  P1 = P2(1:N/2+1);
  P1(2:end-1) = 2*P1(2:end-1);
  f = options.samplingFreq*(0:(N/2))/N;

  % Smooth the power
  df = f(2) - f(1);
  smoothWindow = 1/df;
  P1_smooth = movmean(P1,smoothWindow);

  % Plot the power spectrum
  %figure; loglog(f,P1); hold on
  %loglog(f,P1_smooth); hold off
  %figure; plot(f,P1); hold on
  %plot(f,P1_smooth); hold off
  %xlim([0 12])

  % Identify peak frequency
  P1_smooth_focus = P1_smooth(f>=freqRange(1) & f<=freqRange(2));
  f_focus = f(f>=freqRange(1) & f<=freqRange(2));
  [maxPower(iCh),indexMax] = max(P1_smooth_focus);
  oscFreq(iCh) = f_focus(indexMax);

  % Determine if the signal is oscillatory
  if oscFreq(iCh) > f_focus(1) + 0.05*(diff(freqRange))
    isOsc(iCh) = true;
  end
end