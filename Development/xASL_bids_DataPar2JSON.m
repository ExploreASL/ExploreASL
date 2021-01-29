function xASL_bids_DataPar2JSON(DataParPath)
%xASL_bids_DataPar2JSON Take DataParPath ASL-legacy parameters and move them to JSON sidecars per BIDS
%
% FORMAT: [x] = xASL_bids_FromDataPar2JSON(DataParPath)
% 
% INPUT:
%   DataParPath - path to data parameter JSON file with xASL legacy format (REQUIRED)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes all parameters from the DataPar & moves them into all lower-level JSONs, per BIDS inheritance
% Note that this function assumes that the DataPar file is in the ROOT folder of the study, that contains all the JSON sidecars.
% Also note that this function will recursively create JSON files (if non-existing) for all NIfTI files, so is supposed to run on raw data only.
%
% This function runs the following steps:
% 1) Load parent DataParPath JSON file in xASL legacy format
% 2) Get list of NIfTIs (i.e. "children" that will get the parameters)
% 3) Load & add JSON child (if exist) to memory
% 4) Load & add parms.mat child (if exist (legacy)) to memory
% 5) Add parent fields to memory
% 6) Save (& overwrite if existed) new JSON from memory
% 7) Delete parms.mat if existed
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_bids_FromDataPar2JSON('/MyStudy/DataParameterFile.json');
% __________________________________
% Copyright 2015-2020 ExploreASL

%% 1) Load DataPar file
if nargin<1 || isempty(DataParPath) || ~exist(DataParPath, 'file')
    error('Invalid input argument or non-existing DataPar file');
end

DataPar = xASL_import_json(DataParPath);

%% 2) Get list of NIfTIs
AnalysisDir = fileparts(DataParPath);

fprintf('Converting *_parms.mat to *.json & implementing BIDS inheritance:   ');
FileList = xASL_adm_GetFileList(AnalysisDir, '^(T1|FLAIR|ASL4D|M0)\.nii$','FPListRec',[0 Inf]);

Fields2Skip = {'Quality' 'DELETETEMP' 'subject_regexp' 'name' 'exclusion' 'exclusionReason' 'SESSIONS' 'ROOT'};
% These fields are environment parameters, not ASL-specific parameters

for iList=1:length(FileList)
    xASL_TrackProgress(iList, length(FileList));
    [Fpath, Ffile] = xASL_fileparts(FileList{iList});
    PathJSON = fullfile(Fpath, [Ffile '.json']);
    PathMAT = fullfile(Fpath, [Ffile '_parms.mat']);
    
    %% 3) Load & add JSON child if exist
    if exist(PathJSON, 'file')
        JSON = xASL_import_json(PathJSON);
    else
        JSON = struct;
    end
    
    %% 4) Load & add parms.mat child if exist (legacy)
    if exist(PathMAT, 'file')
        mat = load(PathMAT,'-mat');

        % BIDS timing corrections to SI units
        if isfield(mat.parms,'RepetitionTime')
            mat.parms.RepetitionTime = mat.parms.RepetitionTime/1000;
        end
        if isfield(mat.parms,'EchoTime')
            mat.parms.EchoTime = mat.parms.EchoTime/1000;
        end        
        
        JSON = xASL_bids_InsertJSONFields(mat.parms, JSON, Fields2Skip);
    end
    
    %% 5) Add parent fields
    JSON = xASL_bids_InsertJSONFields(DataPar, JSON, Fields2Skip);
    
    %% 6) Save (& overwrite if existed) new JSON
    spm_jsonwrite(PathJSON, JSON);
    %% 7) Delete parms.mat if existed
    xASL_delete(PathMAT);
end

fprintf('\n');

end
