function chSpikeTimes = getChannelSpikeTimes(spikeTimes, maxWaveformCh, nCh)
% chSpikeTimes = getChannelSpikeTimes(spikeTimes, maxWaveformCh, nCh)
%
% Function takes in spikes times of units or MUAs and outputs spike times
% organised based on probe recording channels instead.
%
% Args:
%   spikeTimes (cell, required, positional): a shape-(N, 1) cell array of
%     spike time numeric arrays corresponding to individual units (or
%     MUAs).
%   maxWaveformCh (numeric, required, positional): a shape-(N, 1) numeric
%     array of recording channels with units' maximal waveform amplitude.
%   nCh (numeric, required, positional): a shape-(1, 1) numeric scalar
%     indicating the number of channels that spike times should be
%     re-organised to.
%
% Returns:
%   chSpikeTimes (cell): a shape-(nCh, 1) cell array of spike times
%     corresponding to individual recording channels. Certain channels can
%     be empty, if no spiking is associated with those channels.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (:,1) {mustBeA(spikeTimes,'cell'),mustBeNonempty}
  maxWaveformCh (:,1) {mustBeNumeric,mustBePositive}
  nCh (1,1) {mustBeNumeric,mustBePositive}
end

chSpikeTimes = cell(nCh,1);
for ch = 1:nCh
  chInds = maxWaveformCh == ch;
  if any(chInds)
    chSpikeTimes{ch} = sort(concatenateCells(spikeTimes(chInds)));
  end
end