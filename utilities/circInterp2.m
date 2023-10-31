function interpolatedMatrix = circInterp2(matrixContainingNaNs)
% interpolatedMatrix = circInterp2(matrixContainingNaNs)
%
% Interpolates missing circular values (NaNs) of a matrix by averaging the
% nearest neighbours. If the neighbours are also missing, the value is not
% interpolated and remains NaN.
%
% Args:
%   matrixContainingNaNs (numeric, required, positional): a shape-(M, N)
%     numeric array of circular values like a phase in radians.
%
% Returns:
%   interpolatedMatrix (numeric): a shape-(M, N) numeric array of
%     interpolated circular values.
%
% Dependencies:
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%   dervinism/circStatNP (https://github.com/dervinism/circStatNP).
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  matrixContainingNaNs (:,:) {mustBeNumeric}
end

% Parse input
nanLocs = isnan(matrixContainingNaNs);
interpolatedMatrix = matrixContainingNaNs;
if ~any(nanLocs)
  return
end

% Fill in all NaNs with the mean of their neighbours
[nanRows,nanCols] = find(nanLocs);
for iNan = 1:sum(sum(nanLocs))
  neighbourRows = max([1 nanRows(iNan)-1]):min([nanRows(iNan)+1 size(matrixContainingNaNs,1)]);
  neighbourCols = max([1 nanCols(iNan)-1]):min([nanCols(iNan)+1 size(matrixContainingNaNs,2)]);
  suroundSize = numel(neighbourRows)*numel(neighbourCols)-1;
  if sum(sum(~isnan(matrixContainingNaNs(neighbourRows,neighbourCols))))/suroundSize >= 0.5
    interpolatedMatrix(nanRows(iNan),nanCols(iNan)) = ...
      datamean(reshape(matrixContainingNaNs(neighbourRows,neighbourCols),[],1), 'circular');
  end
end