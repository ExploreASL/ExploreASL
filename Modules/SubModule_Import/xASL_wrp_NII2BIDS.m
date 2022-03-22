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
% Copyright 2015-2022 ExploreASL


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
		studyPar = struct;
	else
		studyPar = xASL_io_ReadDataPar(x.dir.studyPar, true);
	end
	
	% The name always has to be assigned
	if ~isfield(studyPar,'Name')
		studyPar.Name = x.modules.import.imPar.studyID;
	end
	
	%% 2. Create the study description output and verify that all is there
	datasetDescription = xASL_bids_CreateDatasetDescriptionTemplate(studyPar, x.Version);
	
	% Make the output directory and save the description
    if ~xASL_exist(x.modules.import.imPar.BidsRoot,'dir')
        xASL_adm_CreateDir(x.modules.import.imPar.BidsRoot);
    end
	
	spm_jsonwrite(fullfile(x.modules.import.imPar.BidsRoot,[bidsPar.datasetDescription.filename '.json']),datasetDescription);
	
	%% 3. Go through all subjects and check all the M0 and ASLs and modify the JSONs
	% This step should be completely automatic, just taking the info filled above and using it to convert to full BIDS.
    
	% Go through all subjects
	listSubjectsSessions = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],false,[],true);
    for iSubjectSession = 1:length(listSubjectsSessions)
        % Only run it for the current subject (maybe we can do this more elegantly in the future)
        if ~isempty(regexpi(listSubjectsSessions{iSubjectSession},subjectName))
            x = xASL_wrp_NII2BIDS_Subject(x, bidsPar, studyPar, listSubjectsSessions{iSubjectSession});
        end
    end
    
    % Move log files of current subject from datasetroot & temp to derivatives/ExploreASL if they aren't there already
    logFilesDatasetRoot = xASL_adm_GetFileList(x.dir.DatasetRoot,'^import.+$');
    logFilesTempRoot = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,'^import.+$');
    allLogFiles = vertcat(logFilesDatasetRoot,logFilesTempRoot);
    if ~isempty(allLogFiles)
        for importFile=1:size(allLogFiles,1)
            % Check for subject name
            if ~isempty(regexpi(allLogFiles{importFile},subjectName))
                [~,thisFileMeta,thisExtensionMeta] = xASL_fileparts(allLogFiles{importFile,1});
                xASL_Move(allLogFiles{importFile,1},fullfile(x.modules.import.imPar.DerivativesRoot,'ExploreASL',[thisFileMeta thisExtensionMeta]),1);
            end
        end
    end
    
    % Delete temp folder of current subject
    tempDirs = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],[],[],true);
    for iTempDir=1:numel(tempDirs)
        if ~isempty(regexpi(tempDirs{iTempDir},subjectName))
            xASL_delete(fullfile(tempDirs{iTempDir}), true);
        end
    end
    
    % Delete temp directory if it is empty
    dirsInTemp = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],[],[],true);
    filesInTemp = xASL_adm_GetFileList(x.modules.import.imPar.TempRoot,[],[],[],false);
    if isempty(dirsInTemp) && isempty(filesInTemp)
        xASL_delete(x.modules.import.imPar.TempRoot, true);
    end
    
    % Update x.opts.DatasetRoot
    x = xASL_imp_UpdateDatasetRoot(x);

end



