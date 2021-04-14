function xASL_bids_NII2BIDS(studyParPath, imPar)
%xASL_bids_NII2BIDS Run the NII2BIDS conversion.
%
% FORMAT: xASL_bids_NII2BIDS(studyParPath, imPar)
% 
% INPUT:
%   studyParPath    - path to the JSON file with the BIDS parameters relevant for the whole study
%   imPar           - JSON file with structure with import parameters
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run the NII2BIDS conversion.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_bids_NII2BIDS(studyParPath, imPar);
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Loads the general configuration necessary for the conversion and BIDS saving
	bidsPar = xASL_bids_Config();
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load the study parameters + dataset description
	if ~exist(studyParPath,'file')
		warning('Study-par file is not provided.');
		studyPar = struct;
	else
		studyPar = xASL_io_ReadDataPar(studyParPath);
	end
	
	% The Name has to be always assigned
	if ~isfield(studyPar,'Name')
		studyPar.Name = imPar.studyID;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Create the study description output and verify that all is there
	datasetDescription = xASL_bids_CreateDatasetDescriptionTemplate(studyPar);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Make the output directory and save the description
	
	if ~exist(fullfile(imPar.BidsRoot),'dir')
		mkdir(fullfile(imPar.BidsRoot));
	end
	
	spm_jsonwrite(fullfile(imPar.BidsRoot,[bidsPar.datasetDescription.filename '.json']),datasetDescription);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Go through all subjects and check all the M0 and ASLs and modify the JSONs
	% This step should be completely automatic, just taking the info filled above and using it to convert to full BIDS.
	
	% Go through all subjects
	listSubjects = xASL_adm_GetFileList(imPar.AnalysisRoot,[],false,[],true);
    for iSubject = 1:length(listSubjects)
        xASL_bids_NII2BIDS_Subject(imPar,bidsPar,studyPar,listSubjects,iSubject);
    end
    
    % Copy log files
    importMetaFiles = xASL_adm_GetFileList(imPar.AnalysisRoot,'^import.+$');
    for importFile=1:size(importMetaFiles,1)
        [~,thisFileMeta,thisExtensionMeta] = xASL_fileparts(importMetaFiles{importFile,1});
        xASL_Copy(importMetaFiles{importFile,1},fullfile(studyPath,[thisFileMeta thisExtensionMeta]));
    end
    
    % Delete analysis folder
    xASL_delete(imPar.AnalysisRoot, true);

end



