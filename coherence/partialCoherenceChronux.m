% Inputs: data1, data2, dataR, type, params, win
%   type - 'segpb' (all three signals treated as bin counts) or 'segc' (all three signals treated as continuous)
% Outputs: freq 
%          parC - the partial coherence
%             C - coherence (i.e. of data1 and data2, disregarding dataR)
% Author: Michael Okun.
function [freq, parC, C] = partialCoherenceChronux(data1, data2, dataR, type, params, win)

assert(max(data1) - min(data1)>0 && max(data2) - min(data2)>0, 'data1 or/and data2 is constant?!')

switch type
  case 'segpb'
    assert(nargin >= 6, 'partialCoherenceChronux: win parameter not provided')
    [c, ~, S12, S1, S2, freq] = coherencysegpb(data1, data2, win, params);     % when data1 = data2, S12 = S1 = S2
    [~, ~, S1r, ~, Sr] = coherencysegpb(data1, dataR, win, params);
    [~, ~, Sr2] = coherencysegpb(dataR, data2, win, params); 
  case 'segc'
    assert(nargin >= 6, 'partialCoherenceChronux: win parameter not provided')
    [c, ~, S12, S1, S2, freq] = coherencysegc(torow(data1)', torow(data2)', win, params, true);
    [~, ~, S1r, ~, Sr] = coherencysegc(torow(data1)', torow(dataR)', win, params);
    [~, ~, Sr2] = coherencysegc(torow(dataR)', torow(data2)', win, params); 
%   case 'segcpb'
%     assert(nargin >= 6, 'partialCoherenceChronux: win parameter not provided')
%     [~, ~, S12, S1, S2, freq] = coherencysegcpb(torow(data1)', torow(data2)', win, params, true);
%     [~, ~, S1r, ~, Sr] = coherencysegc(torow(data1)', torow(dataR)', win, params);
%     [~, ~, Sr2] = coherencysegpb(torow(dataR)', torow(data2)', win, params);     
  otherwise
    error('unknown type');
end
R12 = S12 ./ sqrt(S1.*S2);
R1r = S1r ./ sqrt(S1.*Sr);
Rr2 = Sr2 ./ sqrt(Sr.*S2);

parC = abs(R12 - R1r.*Rr2).^2 ./ (1-abs(R1r).^2) ./ (1-abs(Rr2).^2);
parC = sqrt(parC); % The above equation actually gives the partial coherence squared
%(e.g. consider the case when dataR is completely unrelated)

% % Alternative formulation, based on 
% % Multivariate Partial Coherence Analysis for Identification of Neuronal Connectivity from Multiple Electrode Array Recordings 
% % (2014 IEEE Conference on Biomedical Engineering and Sciences, 8 - 10 December 2014, Miri, Sarawak, Malaysia)
% % https://doi.org/10.1109/IECBES.2014.7047613
% S12cond = S12 - S1r.*Sr2 ./ Sr;
% S1cond  = S1  - S1r.*conj(S1r) ./ Sr; % conj(S1r) is the same as Sr1
% S2cond  = S2  - conj(Sr2).*Sr2 ./ Sr;
% parC2 = abs(S12cond).^2 ./  S1cond ./ S2cond;
% parC2 = sqrt(parC2);

C = sqrt(abs(S12).^2 ./ S1 ./ S2); 

if max(dataR) - min(dataR) == 0 % Dividing by Sr produces NaN above
  parC = C;
end

assert( all(C >=0) && all(C <= 1) && all(parC >=0) && all(parC <= 1+1e-9) && max(abs(C-c)) < 1e-6) % because numerical inaccuracies parC can get slightly above 1, see example below
return


%% some code to test the above:
x = randn(1, 1e5);
y = randn(size(x));
z1 = randn(size(x));
z2 = randn(size(x));

params.Fs = 1;
params.tapers = [3/20 20 3]; % 3 tapers, 20s windows

% [freq, parC12, C12] = partialCoherenceChronux(x - z1 - z2, y + z1 + z2, z1+z2, 'segc', params, 20);
% [freq, parC1, C1] = partialCoherenceChronux(x - z1 - z2, y + z1 + z2, z1, 'segc', params, 20);

[freq, parC, C] = partialCoherenceChronux(z1 + z2, z1, z2, 'segc', params, 20); % parC > C
% Cases when partial correlation > correlation involve so called 'suppressor' variables