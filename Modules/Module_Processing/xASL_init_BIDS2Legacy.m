function [x] = xASL_init_BIDS2Legacy(x)
%xASL_init_BIDS2Legacy Prepare relevant fields in xASL x struct for BIDS to Legacy conversion
%
% FORMAT: [x] = xASL_init_BIDS2Legacy(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x          - ExploreASL x structure
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Checks the necessary input directories and loads the BIDS directory structure that should be done only once for all subjects and not
%                 for each subject separately
%
% 1. Check basic directories
% 2. Start with checking dataset_description.json & rawdata
% - 1. The input is dataset_description.json in the rawdata folder
% - 2. The input is dataPar.json or sourceStructure.json - have to look for a rawdata folder
% 3. Check if a dataPar is provided, otherwise use the defaults
% 4. Load the BIDS structure of all subjects
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_BIDS2Legacy(x);
% __________________________________
% Copyright (c) 2015-2023 ExploreASL

% Admin
if nargin < 1 || isempty(x)
	error('x-struct is a required input');
end


%% 1. Check basic directories
if ~isfield(x,'dir')
	error('Missing directories field...');
end

if isempty(x.dir.DatasetRoot)
	error('x.dir.DatasetRoot is a required parameter.');
end

% Verify that the rawdata subfolder exists
if ~exist(fullfile(x.dir.DatasetRoot,'rawdata'), 'dir')
    warning(['Rawdata folder missing: ' fullfile(x.dir.DatasetRoot,'rawdata')]);
    return;
end

% Check derivatives of ExploreASL
if ~isfield(x.dir,'xASLDerivatives')
    error('Missing xASL derivatives field...');
end





%% 4. Load the BIDS structure of all subjects & translate to ExploreASL legacy
x.modules.bids2legacy.BIDS = bids.layout(fullfile(x.dir.DatasetRoot,'rawdata'));

% Remove any pre-existing fields
x.SUBJECTS = {};
x.dataset.TotalSubjects = {};

% Now we translate the BIDS naming convention to ExploreASL's legacy name convention
for iSubj=1:length(x.modules.bids2legacy.BIDS.subjects)
    subjectName = x.modules.bids2legacy.BIDS.subjects(iSubj).name;
    visitName = x.modules.bids2legacy.BIDS.subjects(iSubj).session;
    if ~strcmp(visitName(1:4), 'ses-')
        warning(['Illegal BIDS session definition for ' subjectName '_' visitName]);
    else
        visitName = visitName(5:end);
    end
    
    x.dataset.TotalSubjects{iSubj} = [subjectName '_' visitName];
end

x.SUBJECTS = x.dataset.TotalSubjects;



%% 5. Create the derivatives directory
if exist(x.dir.xASLDerivatives, 'dir')
    fprintf('%s\n', [x.dir.xASLDerivatives ' already exists']);
    fprintf('%s\n', 'Note that all pre-existing derivative subject folders will be overwritten,');
    fprintf('%s\n', 'unless BIDS2Legacy lock files already exist for a subject');
else
    xASL_adm_CreateDir(x.dir.xASLDerivatives);
end


%% 6. Load BIDS configuration for file renaming
x.modules.bids2legacy.bidsPar = xASL_bids_Config;


    
% This is the line needed by xASL_init_Iteration for BIDS2Legacy
x.D.ROOT = x.dir.DatasetRoot;

% SESSIONS DUMMY
x.SESSIONS = {''};


end