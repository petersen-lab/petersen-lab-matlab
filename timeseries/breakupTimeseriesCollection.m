function timeseriesFileList = breakupTimeseriesCollection(file, options)
% timeseriesFilList = breakupTimeseriesCollection(file, <basename>)
%
% Function takes in a timeseries collection MAT file (CellExplorer format:
% https://cellexplorer.org/data-structure/) or a timeseries collection
% Matlab structure that has already been loaded or setup, breaks it up into
% its constituent individual timeseries and seves them in seprate MAT files.
%
% Args:
%   file (char | struct, required, positional): a shape-(1, N)
%     character array containing a MAT filename with a timeseries
%     collection. Alternatively, it is a shape-(1, 1) scalar structure
%     loaded from that type of file or set up accordingly.
%   basename (char, optional, keyword): a shape-(1, M) character array
%     with the new basename for saving individual timeseries data.
%     Individual filenames will be appended by corresponding channel names
%     from the timesieres collection. By default, the basename is kept the
%     same as the timeseries collection file. In case an actual structure
%     variable is supplied as an input, timeseries files are saved in the
%     current working directory by default.
%
% Returns:
%   timeseriesFileList (cell): a shape-(L, 1) cell array of MAT files
%     containing individual timeseries vectors from a breakup collection.
%     Individual files follow CellExplorer timeseries format
%     (https://cellexplorer.org/datastructure/data-structure-and-format/#time-series).
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  file (:,:) {mustBeCharOrListedType(file,'struct')}
  options.basename (1,:) {mustBeA(options.basename,'char')} = '';
end

% Parse input
if isempty(options.basename)
  if ischar(file)
    options.basename = strrep(file, '.timeseriesCollection.mat', '');
    [basefolder, options.basename] = fileparts(options.basename);
    options.basename = fullfile(basefolder, options.basename);
  else
    options.basename = fullfile(pwd, 'basename');
  end
end

% Load the file
if ischar(file)
  fileContents = load(file);
else
  fileContents = file;
end

% Get collection container name
if ischar(file)
  fn = fieldnames(fileContents);
  fileContents = fileContents.(fn{1});
else
  fn = {inputname(1)};
end

% Check whether file contents follow the format requirements
if any(~isfield(fileContents, {'data', 'timestamps', 'precision', 'units', ...
  'nChannels', 'channelNames', 'sr', 'nSamples', 'description', 'processingInfo'}))
  error('The file format does not match CellExplorer timeseries collection format');
end

% Break up collection and save individual timeseries
nSeries = fileContents.nChannels;
timeseriesFileList = cell(nSeries,1);
for iSeries = 1:nSeries
  timeseries = fileContents;
  timeseries.data = fileContents.data(:,iSeries);
  timeseries.units = fileContents.units{iSeries};
  timeseries.nChannels = 1;
  timeseries.channelNames = [fileContents.channelNames{iSeries} '_' fn{1}];
  timeseries.description = [ ...
    'Column ' num2str(iSeries) ' of ' fn{1} ...
    '.data timeseries collection. The full collection description: "' ...
    fileContents.description '"'];
  timeseriesFileList{iSeries} = [options.basename '.' timeseries.channelNames '.timeseries.mat'];
  eval([timeseries.channelNames ' = timeseries;']);
  save(timeseriesFileList{iSeries}, timeseries.channelNames, '-v7.3');
end