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

% Print warning
fprintf('\n========================================== WARNING ===========================================\n');
fprintf(['\nThe support for this ExploreASL functionality will expire with a future release.\n'...
         'If you want to be able to use this function within a future release, please pull\n'...
         'the public <a href="https://github.com/ExploreASL/CustomScripts" rel="nofollow">CustomScripts</a> repository.\n']);
fprintf('\n==============================================================================================\n\n');

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
imParUpdate = spm_jsonread(fullfile(configFiles, [imPar.studyID '.json']));

% Update imPar template with fields from imParUpdate
imParUpdateFieldNames = fieldnames(imParUpdate);
for iField = 1:size(imParUpdateFieldNames,1)
    if ~isempty(imParUpdate.(imParUpdateFieldNames{iField}))
        if ~isempty(strfind(imParUpdateFieldNames{iField},'Aliases'))
            % Retrieve pairs
            it = 1;
            for iPair = 1:round(length(imParUpdate.(imParUpdateFieldNames{iField}))/2)
                imPar.(imParUpdateFieldNames{iField}){iPair,1} = imParUpdate.(imParUpdateFieldNames{iField}){it};
                it = it+1;
                imPar.(imParUpdateFieldNames{iField}){iPair,2} = imParUpdate.(imParUpdateFieldNames{iField}){it};
                it = it+1;
            end
        else
            imPar.(imParUpdateFieldNames{iField}) = imParUpdate.(imParUpdateFieldNames{iField});
        end
    end
end


end
