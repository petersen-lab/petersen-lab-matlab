function [reducedCellArray, reducedCellArrayIdx] = selectCellValues(cellArray, cutoffs)
% [reducedCellArray, reducedCellArrayIdx] = selectCellValues(cellArray, cutoffs)
%
% function selects values of individual cells in a cell array within given
% cutoffs.
%
% Args:
%   cellArray (cell, required, positional): a shape-(K, 1) cell
%     array of shape-(1, N) numeric arrays of spike times where N
%     corresponds to individual values (e.g., spike times). N might be
%     different in different cells.
%   cutoffs (numeric, required, positional): a shape-(L, 2) numeric array
%     of cutoffs for selecting cell values. L represents the number of
%     cutoff intervals, while the second dimension corresponds to the lower
%     and the upper cutoff values. Only elements falling within these
%     cutoff values will remain within individual cells. The rest of the
%     values will be removed.
%
% Returns:
%   reducedCellArray (cell): a shape-(K, 1) cell array of shape-(1, M)
%     numeric arrays of spike times where M corresponds to individual
%     values (e.g., spike times). M might be different in different cells.
%     This cell array contains only selected values from cellArray after
%     applying the cutoffs.
%   reducedCellArrayIdx (numeric): a shape-(K, 1) cell array of linear
%     indices to elements within individual cells of the cellArray
%     producing corresponding elements within individual cells of
%     reducedCellArray
%     (reducedCellArray{k} = cellArray{k}(reducedCellArrayIdx)).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  cellArray (:,1) {mustBeA(cellArray,'cell')}
  cutoffs (:,2) {mustBeNumeric}
end

% Select cell values within given cutoffs
nCells = numel(cellArray);
nCutoffs = size(cutoffs,1);
reducedCellArray = cell(nCells,1);
reducedCellArrayIdx = cell(nCells,1);
for iCell = 1:nCells
  reducedCell = [];
  reducedIdx = [];
  for cutoff = 1:nCutoffs
    idx = find(cellArray{iCell} >= cutoffs(cutoff,1) ...
      & cellArray{iCell} <= cutoffs(cutoff,2));
    reducedIdx = [reducedIdx, idx(:)'];
    cutoffValues = cellArray{iCell}(cellArray{iCell} >= cutoffs(cutoff,1) ...
      & cellArray{iCell} <= cutoffs(cutoff,2));
    reducedCell = [reducedCell, cutoffValues(:)']; %#ok<*AGROW> 
  end
  reducedCellArray{iCell} = reducedCell;
  reducedCellArrayIdx{iCell} = reducedIdx;
end