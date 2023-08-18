function [recentredPhaseData, phaseRange] = recentrePhase(phaseData, phaseCentre)
% phaseData = recentrePhase(phaseData, phaseCentre)
%
% Function converts phase values to values within a chosen phase range
% centred around input variable phaseCentre: [phaseCentre-pi
% phaseCentre+pi].
%
% Args:
%   phaseData (numeric, required, positional): a shape-(M, N) numeric array
%     of phase values expressed in radians.
%   phaseCentre (numeric, required, positional): a shape-(1, 1) numeric
%     scalar representing the centre of the new phase range in radians.
%
% Returns:
%   recentredPhaseData (numeric): a shape-(M, N) numeric array of phase
%     values expressed in radians recentred around the phaseCentre.
%   phaseRange (numeric): a shape-(1,2) numeric array describing the range
%     of phase values ([lowerLimit upperLimit]) in recentredPhaseData.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  phaseData {mustBeNumeric}
  phaseCentre (1,1) {mustBeNumeric,mustBeNonNan,mustBeReal}
end

% Define the new range of phase values
upperLimit = phaseCentre + pi;
lowerLimit = phaseCentre - pi;

% Convert phase data into values limited by the desired range
for i = 1:numel(phaseData)
  if phaseData(i) > upperLimit
    while phaseData(i) > upperLimit
      phaseData(i) = phaseData(i) - 2*pi;
    end
  elseif phaseData(i) < lowerLimit
    while phaseData(i) < lowerLimit
      phaseData(i) = phaseData(i) + 2*pi;
    end
  end
end

% Assign output
recentredPhaseData = phaseData;
phaseRange = [lowerLimit, upperLimit];