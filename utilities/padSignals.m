function [paddedSignal1, paddedSignal2] = padSignals(signal1, signal2)
% [paddedSignal1, paddedSignal2] = padSignals(signal1, signal2)
%
% Function takes in two vector signals and pads the shorter one with zeros
% to make both signals equal in length.
%
% Args:
%   signal1 (numeric, required, positional): a shape-(1, N) numeric array;
%     could be a spike count vector.
%   signal2 (numeric, required, positional): a shape-(1, M) numeric array;
%     could be a spike count vector.
%
% Returns:
%   paddedSignal1 (numeric): a shape-(1, N) or shape-(1, M) numeric array
%     that is a zero-end-padded signal1 array.
%   paddedSignal2 (numeric): a shape-(1, N) or shape-(1, M) numeric array
%     that is a zero-end-padded signal2 array.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  signal1 {mustBeVector}
  signal2 {mustBeVector}
end

% Parse inputs
signal1 = signal1(:)';
signal2 = signal2(:)';

% Pad signals
if numel(signal1) > numel(signal2)
  signal2 = [signal2 zeros(1, numel(signal1)-numel(signal2))];
elseif numel(signal2) > numel(signal1)
  signal1 = [signal1 zeros(1, numel(signal2)-numel(signal1))];
end

% Assign output
paddedSignal1 = signal1;
paddedSignal2 = signal2;