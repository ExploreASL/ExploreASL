function imPar = ExploreASL_ImportConfig(StudyRoot)
%ExploreASL_ImportConfig Configures the import parameters used by ExploreASL_Import
%
% FORMAT:      imPar = ExploreASL_ImportConfig(StudyRoot)
% 
% INPUT:       root of study folder containing DICOMs, e.g. '//MyDisk/MyStudy'
%
% OUTPUT:      imPar which is input to ExploreASL_Import
%
% DESCRIPTION: Please read the help of ExploreASL_Import for more
%              information.
%              Insert your `JSON` study file into the
%              `ExploreASL/Development/ConfigFiles` directory.
%              Call this script with the name of the `JSON` file, e.g.:
%              imPar = ExploreASL_ImportConfig('incoming');
%
% EXAMPLE:     imPar = ExploreASL_ImportConfig('incoming');
% __________________________________
% Copyright 2015-2021 ExploreASL

% Initialize ExploreASL
[x] = ExploreASL_Initialize;

if strcmp(StudyRoot(end),'\') || strcmp(StudyRoot(end),'/')
    StudyRoot = StudyRoot(1:end-1); % bugfix
end

% Get file path, name and extension
[fpath, fname, fext] = fileparts(StudyRoot);

% Set imPar fields accordingly
imPar.studyID = [fname fext];
imPar.AnalysisRoot = fpath;
imPar.RawRoot = fpath;

%% -----------------------------------------------------------------------------
%% Initialize the study-specific parameters
%% -----------------------------------------------------------------------------

imPar.folderHierarchy = {}; % must define this per study; use a cell array of regular expressions. One cell per directory level.
imPar.tokenOrdering = []; % must match imPar.folderHierarchy: 1==subject, 2=visit, 3==session, 4==scan (if visit or session are omitted, they will be skipped)
imPar.tokenScanAliases = [];
imPar.tokenVisitAliases = [];
imPar.tokenSessionAliases = [];
imPar.bMatchDirectories  = false;

% Define configFiles folder
configFiles = fullfile(x.MyPath, 'Development', 'ConfigFiles');

% Load imPar file
imParUpdate = spm_jsonread(fullfile(configFiles, [StudyRoot '.json']));

% Update imPar template with fields from imParUpdate
imParUpdateFieldNames = fieldnames(imParUpdate);
for iField = 1:size(imParUpdateFieldNames,1)
    imPar.(imParUpdateFieldNames{iField}) = imParUpdate.(imParUpdateFieldNames{iField});
end


end
