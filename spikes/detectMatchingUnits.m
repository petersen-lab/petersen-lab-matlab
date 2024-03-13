function [matchingPairs, totalScoreThr, uniqueUnitIDs, unitIDArray, ...
  recSessions, unitIDTable] = detectMatchingUnits(matchTableFilename, options)
% [matchingPairs, totalScoreThr, uniqueUnitIDs, unitIDArray, recSessions, ...
%   unitIDTable] = detectMatchingUnits(matchTableFilename, <save>)
%
% Function detects matching units in the MatchTable output of the UnitMatch
% pipeline (https://github.com/EnnyvanBeest/UnitMatch).
%
% Args:
%   matchTableFilename (char, required, positional): a shape-(1, N)
%     character array containing the MAT filename with MatchTable (i.e.,
%     UnitMatch.mat; give the full path though).
%   save (char, optional, keyword): a shape-(1, M) character array with the
%     name of the file to save the output variables of the function. If
%     left empty, the output variables will not be saves (default).
%
% Returns:
%   matchingPairs (logical): a shape-(L, 1) logical array marking matching
%     unit pairs. The number of elements in the vector is equal to the
%     number of rows in MatchTable. True or matching pairs exceed
%     totalScoreThr.
%   totalScoreThr (numeric): a shape-(1, 1) numeric scalar corresponding to
%     the total score threshold as defined in van Beest et al. (2023)
%     Tracking neurons across days with high-density probes
%     (https://doi.org/10.1101/2023.10.12.562040).
%   uniqueUnitIDs (numeric): a shape-(K, 1) numeric array of unique unit
%     IDs. These IDs are not the same as the ones in MatchTable.
%   unitIDArray (numeric): a shape-(K, J) numeric array of original unit
%     IDs (unique for individual sections but not globally). Rows of this
%     matrix corresponds to unique IDs and columns correspond to
%     consequetive recording sessions.
%   recSessions (cell): a shape-(J, I) cell array of character arrays
%     containing names (IDs) of recording sessions.
%   unitIDTable (table): a shape-(K, J+1) table corresponding to
%     unitIDArray but with an extra column with unique unit IDs.
%
% Dependencies:
%   UnitMatch pipeline (https://github.com/EnnyvanBeest/UnitMatch).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  matchTableFilename (1,:) {mustBeA(matchTableFilename, 'char')}
  options.save (1,:) {mustBeA(options.save, 'char')} = ''
end

% Parameters
nBins = 100;
binSize = 1/nBins;
bins = binSize/2:binSize:1;

% Calculate same unit and neighbouring pair total score distributions
load(matchTableFilename); %#ok<*LOAD> 
sameUnits = MatchTable.RecSes1 == MatchTable.RecSes2 & ...
  MatchTable.ID1 == MatchTable.ID2; %#ok<*USENS> 
sameUnitTotalScores = MatchTable.TotalScore(sameUnits);
sameUnitDistro = hist(sameUnitTotalScores, bins); %#ok<*HIST>
normalisedSameUnitDistro = sameUnitDistro./sum(sameUnitDistro);

neighbourUnits = MatchTable.RecSes1 == MatchTable.RecSes2 & ...
  MatchTable.ID1 ~= MatchTable.ID2;
neighbourUnitTotalScores = MatchTable.TotalScore(neighbourUnits);
neighbourUnitDistro = hist(neighbourUnitTotalScores, bins);
normalisedNeighbourUnitDistro = neighbourUnitDistro./sum(neighbourUnitDistro);

% Find the total score threshold
meanSame = mean(sameUnitTotalScores,'omitnan');
meanNeighbour = mean(neighbourUnitTotalScores,'omitnan');
thrApproximateArea = bins >= meanNeighbour & bins <= meanSame & ...
  normalisedSameUnitDistro > normalisedNeighbourUnitDistro;
potentialThrInd = find(thrApproximateArea,1);
totalScoreThr = max([0.5 bins(potentialThrInd)]);

% figure; plot(bins, normalisedSameUnitDistro); hold on
% plot(bins, normalisedNeighbourUnitDistro);
% yLim = ylim;
% plot([totalScoreThr totalScoreThr], [0 yLim(2)], 'r:')

% Reorder the match table
reducedMatchTable = [double(table2array(MatchTable(:,1:4))) ...
  double(table2array(MatchTable(:,8)))];
matchingPairs = reducedMatchTable(:,5) >= totalScoreThr;
matchingPairsInds = find(matchingPairs);
reducedMatchTable = reducedMatchTable(matchingPairs,:);
nPairs = size(reducedMatchTable,1);
nCols = size(reducedMatchTable,2);
for pair = 1:nPairs
  if reducedMatchTable(pair,4) < reducedMatchTable(pair,3)
    reducedMatchTable(pair,:) = reducedMatchTable(pair, [2 1 4 3 5]);
  end
end

% Eliminate one-to-many matches across sections and section gaps
for pair = 1:nPairs
  if reducedMatchTable(pair,3) == reducedMatchTable(pair,4)
    reducedMatchTable(pair,:) = NaN(1,nCols);
  elseif abs(reducedMatchTable(pair,3) - reducedMatchTable(pair,4)) > 1
    reducedMatchTable(pair,:) = NaN(1,nCols);
  else
    relevantEntries = (reducedMatchTable(:,1) == reducedMatchTable(pair,1) & ...
      reducedMatchTable(:,3) == reducedMatchTable(pair,3) & ...
      reducedMatchTable(:,4) == reducedMatchTable(pair,4)) | ...
      (reducedMatchTable(:,2) == reducedMatchTable(pair,1) & ...
      reducedMatchTable(:,3) == reducedMatchTable(pair,3) & ...
      reducedMatchTable(:,3) == reducedMatchTable(pair,4));
    relevantEntries(pair) = false;
    if any(reducedMatchTable(relevantEntries,5) >= reducedMatchTable(pair,5))
      reducedMatchTable(pair,:) = NaN(1,nCols);
    end
  end
end

% Eliminate many-to-one matches across sections
reducedMatchTable = reducedMatchTable(end:-1:1, [2 1 4 3 5]);
for pair = 1:nPairs
  if isnan(reducedMatchTable(pair,1))
    continue
  else
    relevantEntries = (reducedMatchTable(:,1) == reducedMatchTable(pair,1) & ...
      reducedMatchTable(:,3) == reducedMatchTable(pair,3) & ...
      reducedMatchTable(:,4) == reducedMatchTable(pair,4)) | ...
      (reducedMatchTable(:,2) == reducedMatchTable(pair,1) & ...
      reducedMatchTable(:,3) == reducedMatchTable(pair,3) & ...
      reducedMatchTable(:,3) == reducedMatchTable(pair,4));
    relevantEntries(pair) = false;
    if any(reducedMatchTable(relevantEntries,5) >= reducedMatchTable(pair,5))
      reducedMatchTable(pair,:) = NaN(1,nCols);
    end
  end
end
reducedMatchTable = reducedMatchTable(end:-1:1, [2 1 4 3 5]);

% Final non-unique match elimination
reducedMatchTableExt = [reducedMatchTable; reducedMatchTable(:, [2 1 4 3 5])];
for pair = 1:nPairs
  if isnan(reducedMatchTableExt(pair,1))
    continue
  else
    relevantEntries = (reducedMatchTableExt(:,1) == reducedMatchTableExt(pair,1) & ...
      reducedMatchTableExt(:,3) == reducedMatchTableExt(pair,3) & ...
      reducedMatchTableExt(:,4) == reducedMatchTableExt(pair,4)) | ...
      (reducedMatchTableExt(:,2) == reducedMatchTableExt(pair,1) & ...
      reducedMatchTableExt(:,3) == reducedMatchTableExt(pair,3) & ...
      reducedMatchTableExt(:,3) == reducedMatchTableExt(pair,4));
    relevantEntries(pair) = false;
    if any(reducedMatchTableExt(relevantEntries,5) >= reducedMatchTableExt(pair,5))
      reducedMatchTableExt(pair,:) = NaN(1,nCols);
    end
  end
end
nanInds = isnan(reducedMatchTableExt(:,5));
nanInds = unique([find(nanInds(1:nPairs)); find(nanInds(nPairs+1:end))]);
reducedMatchTable(nanInds,:) = NaN(numel(nanInds),nCols);

% Find matching pairs across recording sessions
matchingPairs(matchingPairsInds(isnan(reducedMatchTable(:,5)))) = false;
reducedMatchTable = reducedMatchTable(~isnan(reducedMatchTable(:,5)),:);

% Assign unique IDs
% Initial IDs
entries = NaN(size(reducedMatchTable,1),6);
unitCount = 0;
for pair = 1:size(reducedMatchTable,1)
  entry = reducedMatchTable(pair,1:4);
  if any(entries(:,1) == entry(1) & entries(:,3) == entry(3))
    entry = [entry unique(entries(find(entries(:,1) == entry(1) & entries(:,3) == entry(3)),5)) ...
      unique(entries(find(entries(:,1) == entry(1) & entries(:,3) == entry(3)),5))]; %#ok<*AGROW,*FNDSB> 
  else
    unitCount = unitCount + 1;
    entry = [entry unitCount unitCount];
  end
  entries(pair,:) = entry;
end
% Eliminate crossovers or overlaps
for pair = 1:size(reducedMatchTable,1)
  entry = reducedMatchTable(pair,1:4);
  inds = find((entries(:,2) == entry(2) & entries(:,4) == entry(4)) | ...
    (entries(:,1) == entry(2) & entries(:,3) == entry(4)));
  IDs = unique([entries(inds,5); entries(inds,6)]);
  smallestID = unique(min(IDs));
  entries(inds,5) = smallestID;
  entries(inds,6) = smallestID;
end
uniqueUnitIDs = unique(entries(:,5:6));

% Identify sessions and original unit IDs for all matching units
sessions = unique(MatchTable.RecSes1);
nSessions = numel(sessions);
unitIDArray = NaN(max(uniqueUnitIDs),nSessions);
for pair = 1:size(entries,1)
  unitIDArray(entries(pair,5),entries(pair,3)) = entries(pair,1);
  unitIDArray(entries(pair,5),entries(pair,4)) = entries(pair,2);
end
unitIDArray = unitIDArray(uniqueUnitIDs,:);
uniqueUnitIDs = (1:numel(uniqueUnitIDs))';
recSessions = cell(nSessions,1);
recSessionsStrArray = "uniqueUnitIDs";
for session = 1:nSessions
  recSessions{session} = UMparam.RawDataPaths{session}.name(1:end-4);
  recSessionsStrArray = [recSessionsStrArray, ...
    string(UMparam.RawDataPaths{session}.name(1:end-4))];
end

% Create the unit ID table
unitIDTable = array2table([uniqueUnitIDs unitIDArray], ...
  'VariableNames',recSessionsStrArray);

% Save the function output
if ~isempty(options.save)
  if endsWith(options.save, '.mat')
    folder = fileparts(options.save);
    if ~exist(folder,'dir')
      mkdir(folder);
    end
    save(options.save, 'matchingPairs','totalScoreThr','uniqueUnitIDs', ...
      'unitIDArray','recSessions', '-v7.3');
  end
end