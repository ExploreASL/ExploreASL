function imPar = ExploreASL_ImportConfig(StudyRoot)
%ExploreASL_ImportConfig Configures the import parameters used by ExploreASL_Import
%
% FORMAT:      imPar = ExploreASL_ImportConfig(StudyRoot)
% 
% INPUT:       root of study folder containing DICOMs, e.g. '//MyDisk/MyStudy'
% OUTPUT:      imPar which is input to ExploreASL_Import
%
% DESCRIPTION: Please read the help of ExploreASL_Import for more information
% EXAMPLE:     n/a
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

% Select CustomScripts directory
if usejava('desktop')
    customScripts = uigetdir(x.MyPath,'Please select the CustomScripts directory...');
else
    customScripts = input('Please insert the CustomScripts directory: ');
end

% Load imPar file
imPar = spm_jsonread(fullfile(customScripts, 'ConfigFiles', [StudyRoot '.json']));


end
