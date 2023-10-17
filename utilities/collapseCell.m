function [collapsedCellArray, nonEmptyInds] = collapseCell(cellArray, options)
% collapsedCellArray = collapseCell(cellArray, <options>)
%
% Collapses a 2-D cell array along the specified dimension without
% preserving any dimensionality of individual numeric arrays contained
% within each cell. Cells of the resulting 1-D cell array contain numeric
% row vectors with original numeric arrays flattened along the specified
% dimension (the same as the cell array collapse dimension).
%
% Args:
%   cellArray (cell, required, positional): a shape-(M, N) cell array of
%     numeric arrays.
%   dim (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     indicating the direction along which the cell array should be
%     collapsed (and individual cells to be flattened). By default the cell
%     array is collapsed along the first dimension. Individual numeric
%     arrays within each cell are always flattened along the same dimension
%     as the collapse dimension.
%   sortElements (logical, optional, keyword): a shape-(1, K) character
%     array to indicate whether numeric values within the elements of the
%     concatenated cell array should also be sorted in a particular order.
%     Available choices include: 'ascend', 'descend', 'none' (default).
%
% Returns:
%   collapsedCellArray (cell): a shape-(M, 1) or shape-(1, N) cell array
%     with row vectors of concatenated and flattened original cells.
%   nonEmptyInds (numeric): a shape-(1, L) numeric array of indices
%     refering to non-empty elements of collapsedCellArray.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  cellArray (:,:) {mustBeA(cellArray,'cell')}
  options.dim (1,1) {mustBeNumeric,mustBePositive,mustBeLessThan(options.dim,3)} = 1
  options.sortElements (1,:) {mustBeMember(options.sortElements,{'ascend','descend','none'})} = 'none'
end

% Create output containers
if options.dim == 1
  nOutputCells = size(cellArray,1);
  cells2collapse = size(cellArray,2);
  collapsedCellArray = cell(nOutputCells,1);
elseif options.dim == 2
  nOutputCells = size(cellArray,2);
  cells2collapse = size(cellArray,1);
  collapsedCellArray = cell(nOutputCells,1)';
end
nonEmptyInds = zeros(1,nOutputCells);

% Collapse the cell array
for iCellOut = 1:nOutputCells
  for iCellIn = 1:cells2collapse
    if options.dim == 1
      flattenedElement = reshape(cellArray{iCellOut,iCellIn}, ...
        [1,numel(cellArray{iCellOut,iCellIn})]);
    elseif options.dim == 2
      flattenedElement = reshape(cellArray{iCellIn,iCellOut}', ...
        [1,numel(cellArray{iCellIn,iCellOut})]);
    end
    if strcmpi(options.sortElements,'none')
      collapsedCellArray{iCellOut} = [collapsedCellArray{iCellOut} ...
        flattenedElement];
    else
      collapsedCellArray{iCellOut} = sort([collapsedCellArray{iCellOut} ...
        flattenedElement], options.sortElements);
    end
  end
  if ~isempty(collapsedCellArray{iCellOut})
    nonEmptyInds(iCellOut) = 1;
  end
end

% Find non-empty cell indices
nonEmptyInds = find(nonEmptyInds);