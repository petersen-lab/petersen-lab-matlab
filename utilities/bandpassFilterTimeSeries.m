function timeSeries = bandpassFilterTimeSeries(timeSeries, options)
% filtTimeSeries = bandpassFilterTimeSeries(timeSeries, <options>)
%
% Band-pass filter a row matrix containing time series vectors.
%
% Args:
%   timeSeries (numeric, required, positional): a shape-(N, M) numeric
%     array of continuous time series vectors (rows).
%   sampleRate (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%    corresponding to sampling frequency in Hz (default=500).
%   frequencyRange (numeric, optional, keyword): a shape-(1, 2) numeric
%    array with the band-pass filtering frequency range.
%
% Returns:
%   filtTimeSeries (numeric): a shape-(N, M) numeric array of band-pass
%     filtered time series vectors (rows). Rows containing only NaNs are
%     skipped while rows containing some NaNs elicit an error.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  timeSeries (:,:) {mustBeNumeric}
  options.sampleRate (1,1) {mustBeNumeric,mustBePositive} = 500
  options.frequencyRange (1,2) {mustBeNumeric,mustBeVector,mustBePositive} = [4 12]
end

% Filter time series data
Wn_theta = [options.frequencyRange(1)/(options.sampleRate/2) options.frequencyRange(2)/(options.sampleRate/2)]; % Band-pass frequencies normalised by the Nyquist frequency
[btheta,atheta] = butter(3,Wn_theta); % Apply butterworth filter
for row = 1:size(timeSeries,1)
  if sum(isnan(timeSeries(row,:))) && sum(isnan(timeSeries(row,:))) < numel(timeSeries(row,:))
    error('Input data contains NaNs');
  elseif ~sum(isnan(timeSeries(row,:)))
    timeSeries(row,:) = filtfilt(btheta,atheta,double(timeSeries(row,:)));
  end
end