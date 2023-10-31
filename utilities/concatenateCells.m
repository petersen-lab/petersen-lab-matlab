function concatenatedMatrix = concatenateCells(cellArray, dimension)
% concatenatedMatrix = concatenateCells(cellArray, dimension)
%
% Concatenates elements of a one-dimensional cell array into a matrix and
% outputs it.
%
% Args:
%   cellArray (cell): a shape-(M, N) cell array where M = 1 or N = 1. The
%     elements are expected to be numeric.
%   dimension (numeric, optional): a shape-(1, 1) numeric type variable
%     indicating the dimension along which to concatenate. 1 is default.
%
% Returns:
%   concatenatedMatrix (numeric): a shape-(M, N) numeric matrix containing
%     concatentated elements of the input cell array.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  cellArray cell
  dimension (1,1) {mustBeNumeric,mustBePositive} = 1
end

% Concatenate
nElements = numel(cellArray);
concatenatedMatrix = cellArray{1};
assert(isnumeric(concatenatedMatrix), 'Elements of the supplied cell array must be numeric.');
for iElement = 2:nElements
  if isempty(concatenatedMatrix)
    concatenatedMatrix = cellArray{iElement};
  elseif ~isempty(cellArray{iElement})
    assert(isnumeric(cellArray{iElement}), 'Elements of the supplied cell array must be numeric.');
    if dimension == 1
      assert(size(concatenatedMatrix,2) == size(cellArray{iElement},2), ...
        'Elements in the supplied cell array have unequal number of columns');
      concatenatedMatrix = [concatenatedMatrix; cellArray{iElement}]; %#ok<*AGROW>
    elseif dimension == 2
      assert(size(concatenatedMatrix,1) == size(cellArray{iElement},1), ...
        'Elements in the supplied cell array have unequal number of rows');
      concatenatedMatrix = [concatenatedMatrix, cellArray{iElement}];
    end
  end
end

assert(isnumeric(concatenatedMatrix));