function x = xASL_imp_NII2BIDS(x, imPar, studyPath, studyParPath)
%xASL_imp_NII2BIDS Run the NII2BIDS conversion.
%
% FORMAT: x = xASL_imp_NII2BIDS(x, imPar, studyPath, studyParPath)
% 
% INPUT:
%   x               - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar           - JSON file with structure with import parameters (REQUIRED, STRUCT)
%   studyPath       - Path to the study directory containing the 'sourcedata' directory with the DICOM files (REQUIRED, CHAR ARRAY)
%   studyParPath    - Path to the JSON file with the BIDS parameters relevant for the whole study (REQUIRED, CHAR ARRAY)
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
% EXAMPLE:     xASL_imp_NII2BIDS(x, imPar, studyPath, studyParPath);
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Run the NII2BIDS conversion
    diary(fullfile(studyPath,'xASL_module_Import.log'));
    fprintf('================================== NIFTI to BIDS CONVERSION ==================================\n');
    
    % Create lock directory for import if it does not exist already
    xASL_adm_CreateDir(x.modules.import.dir.lockImport);
    
    % Create NII2BIDS sub-directory
    xASL_adm_CreateDir(x.modules.import.dir.lockNII2BIDS);
    
    % Check if the temp folder exists
    existTempRoot = xASL_exist(fullfile(studyPath,'derivatives','ExploreASL','temp'),'dir');
    if ~existTempRoot
        error('The temp directory does not exist. Please run DICOM to NIfTI on your sourcedata first...');
    end

    % Loads the general configuration necessary for the conversion and BIDS saving
	bidsPar = xASL_bids_Config();
	
	%% 1. Load the study parameters + dataset description
	if ~exist(studyParPath,'file')
		warning('Study-par file is not provided.');
		studyPar = struct;
	else
		studyPar = xASL_io_ReadDataPar(studyParPath, true);
	end
	
	% The Name has to be always assigned
	if ~isfield(studyPar,'Name')
		studyPar.Name = imPar.studyID;
	end
	
	%% 2. Create the study description output and verify that all is there
	datasetDescription = xASL_bids_CreateDatasetDescriptionTemplate(studyPar, x.Version);
	
	% Make the output directory and save the description
    if xASL_exist(imPar.BidsRoot,'dir')
        fprintf('The rawdata BIDS directory exists already...\n');
    else
        xASL_adm_CreateDir(imPar.BidsRoot);
    end
	
	spm_jsonwrite(fullfile(imPar.BidsRoot,[bidsPar.datasetDescription.filename '.json']),datasetDescription);
	
	%% 3. Go through all subjects and check all the M0 and ASLs and modify the JSONs
	% This step should be completely automatic, just taking the info filled above and using it to convert to full BIDS.
    
	% Go through all subjects
	listSubjectsSessions = xASL_adm_GetFileList(imPar.TempRoot,[],false,[],true);
    for iSubjectSession = 1:length(listSubjectsSessions)
        x = xASL_imp_NII2BIDS_Subject(x,imPar,bidsPar,studyPar,listSubjectsSessions{iSubjectSession});
    end
    
    % Copy log files
    importMetaFiles = xASL_adm_GetFileList(imPar.TempRoot,'^import.+$');
    for importFile=1:size(importMetaFiles,1)
        [~,thisFileMeta,thisExtensionMeta] = xASL_fileparts(importMetaFiles{importFile,1});
        xASL_Copy(importMetaFiles{importFile,1},fullfile(studyPath,[thisFileMeta thisExtensionMeta]));
    end
    
    % Delete temp folder
    xASL_delete(imPar.TempRoot, true);

end



