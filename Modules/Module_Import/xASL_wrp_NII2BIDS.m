function x = xASL_wrp_NII2BIDS(x)
%xASL_wrp_NII2BIDS Run the NII2BIDS conversion.
%
% FORMAT: x = xASL_wrp_NII2BIDS(x)
% 
% INPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   x.modules.import.imPar - JSON file with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   x               - ExploreASL x structure (STRUCT)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run the NII2BIDS conversion.
%
% 1. Load the study parameters + dataset description
% 2. Create the study description output and verify that all is there
% 3. Go through all subjects and check all the M0 and ASLs and modify the JSONs
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_wrp_NII2BIDS(x);
%
% __________________________________
% Copyright 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    %% Run the NII2BIDS conversion
    
    % Make sure that logging is still active
    diary(x.dir.diaryFile);
    
    % Print feedback
    xASL_adm_BreakString('NIFTI to BIDS CONVERSION');
    
    % We do not iterate over subjects anymore, since this is done in xASL_Iteration now
    iSubject = strcmp(x.SUBJECT,x.SUBJECTS);
    subjectName = x.SUBJECTS{iSubject};
    
    % Check if the temp folder exists
    existTempRoot = xASL_exist(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL','temp'),'dir');
    if ~existTempRoot
        error('The temp directory does not exist. Please run DICOM to NIfTI on your sourcedata first...');
    end
    % Loads the general configuration necessary for the conversion and BIDS saving
	bidsPar = xASL_bids_Config();
	
	%% 1. Load the study parameters + dataset description
	if ~exist(x.dir.studyPar,'file')
		warning('Could not find the studyPar.json file...');
		studyParAll = struct;
	else
		studyParAll = xASL_io_ReadDataPar(x.dir.studyPar, true);
	end
	
	% The name always has to be assigned as it is used in the DatasetDescription
	
	if isfield(studyParAll,'StudyPars')
		% For multi-studyParameter studyPar - add this to the first studyPar
		if ~isfield(studyParAll.StudyPars{1},'Name')
			studyParAll.StudyPars{1}.Name = x.modules.import.imPar.studyID;
		end
		studyParFirstParameters = studyParAll.StudyPars{1};
	else
		if ~isfield(studyParAll,'Name')
			studyParAll.Name = x.modules.import.imPar.studyID;
		end
		studyParFirstParameters = studyParAll;
	end
	
	%% 2. Create the study description output and verify that all is there
	datasetDescription = xASL_bids_CreateDatasetDescriptionTemplate(studyParFirstParameters, x.Version);
	
	% Make the output directory and save the description
    if ~xASL_exist(x.modules.import.imPar.BidsRoot,'dir')
        xASL_adm_CreateDir(x.modules.import.imPar.BidsRoot);
    end
	
	xASL_io_WriteJson(fullfile(x.modules.import.imPar.BidsRoot,[bidsPar.datasetDescription.filename '.json']),datasetDescription);
	
	%% 3. Go through all subjects and check all the M0 and ASLs and modify the JSONs
	% This step should be completely automatic, just taking the info filled above and using it to convert to full BIDS.
    
	% Go through all subjects in this directory (but proceed only with the current subject)
	listSubjectsSessions = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],false,[],true);
    for iSubjectSession = 1:length(listSubjectsSessions)
        % Only run it for the current subject (maybe we can do this more elegantly in the future)
        if ~isempty(regexpi(listSubjectsSessions{iSubjectSession}, subjectName, 'once'))
            x = xASL_wrp_NII2BIDS_Subject(x, bidsPar, studyParAll, listSubjectsSessions{iSubjectSession});
        end
    end
    
    % Move log files of current subject from datasetroot & temp to derivatives/ExploreASL 
    logFilesDatasetRoot = xASL_adm_GetFileList(x.dir.DatasetRoot,'^import.+$');
    logFilesTempRoot = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,'^import.+$');
    allLogFiles = vertcat(logFilesDatasetRoot,logFilesTempRoot);
    if ~isempty(allLogFiles)
        for importFile=1:size(allLogFiles,1)
            % Check for subject name - only move them when they have not been moved already
            if ~isempty(regexpi(allLogFiles{importFile}, subjectName, 'once'))
                [~,thisFileMeta,thisExtensionMeta] = xASL_fileparts(allLogFiles{importFile,1});
                xASL_Move(allLogFiles{importFile,1},fullfile(x.modules.import.imPar.DerivativesRoot,'ExploreASL',[thisFileMeta thisExtensionMeta]),1);
            end
        end
    end
    
    % Delete temp folder of all the subjects that might be there
    tempDirs = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],[],[],true);
    for iTempDir=1:numel(tempDirs)
        if ~isempty(regexpi(tempDirs{iTempDir}, subjectName, 'once'))
            xASL_delete(fullfile(tempDirs{iTempDir}), true);
        end
    end
    
    % Delete temp directory if it is empty
    dirsInTemp = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],[],[],true);
    filesInTemp = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],[],[],false);
    if isempty(dirsInTemp) && isempty(filesInTemp)
        xASL_delete(x.modules.import.imPar.TempRoot, true);
    end
end