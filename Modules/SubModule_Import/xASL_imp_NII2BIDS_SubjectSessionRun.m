function xASL_imp_NII2BIDS_SubjectSessionRun(imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun)
%xASL_imp_NII2BIDS_SubjectSessionRun NII2BIDS conversion for a single sessions, single run.
%
% FORMAT: xASL_imp_NII2BIDS_SubjectSessionRun(bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun)
% 
% INPUT:
% imPar          - JSON file with structure with import parameter (STRUCT, REQUIRED)
% bidsPar             - Output of xASL_imp_Config (STRUCT, REQUIRED)
% studyPar            - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
% subjectSessionLabel - subject-session label (CHAR ARRAY, REQUIRED)
% inSessionPath       - input session path (CHAR ARRAY, PATH, REQUIRED)
% outSessionPath      - output session path (CHAR ARRAY, PATH, REQUIRED)
% listRuns            - list of runs (INTEGER, REQUIRED)
% iRun                - run number (INTEGER, REQUIRED)
%
% OUTPUT:
% n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single sessions, single run.
% 
% 1. Define the pathnames
% 2. Load the JSONs and NIfTI information
% 3. BIDSify ASL
% 4. Prepare the link to M0 in ASL.json
% 5. BIDSify M0
% 6. Save all ASL files (JSON, NIFTI, CONTEXT) to the BIDS directory
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

%% 1. Define the pathnames
if length(listRuns)
	aslLabel = ['ASL4D_' num2str(iRun)];
	aslOutLabel = fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectSessionLabel '_run-' num2str(iRun)]);
	aslOutLabelRelative = fullfile(bidsPar.strPerfusion,['sub-' subjectSessionLabel '_run-' num2str(iRun)]);
else
	aslLabel = 'ASL4D';
	aslOutLabel = fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectSessionLabel]);
	aslOutLabelRelative = fullfile(bidsPar.strPerfusion,['sub-' subjectSessionLabel]);
end

%% 2. Load the JSONs and NIfTI information
if exist(fullfile(inSessionPath, [aslLabel '.json']),'file')
	jsonDicom = spm_jsonread(fullfile(inSessionPath, [aslLabel '.json']));
else
	error('Missing file: %s\n',fullfile(inSessionPath, [aslLabel '.json']));
end
if xASL_exist(fullfile(inSessionPath, [aslLabel '.nii']),'file')
	headerASL = xASL_io_ReadNifti(fullfile(inSessionPath, [aslLabel '.nii']));
else
	error('Missing file: %s\n\',fullfile(inSessionPath, [aslLabel '.nii']));
end
	
%% 3. BIDSify ASL
% Merge the information from DICOM, manually entered parameters and BIDSify
jsonLocal = xASL_bids_BIDSifyASLJSON(jsonDicom, studyPar, headerASL);

%% 4. Prepare the link to M0 in ASL.json	
% Define the M0 type	
[jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, fullfile(inSessionPath,'M0.nii'), fullfile(bidsPar.strPerfusion,['sub-' subjectSessionLabel]));	
	
%% 5. BIDSify M0	
% Check the M0 files only for the first ASL-image run within a session
if iRun == 1 % For first run process also M0, for the second run, only reference it
	for iReversedPE = 1:2 % Check if AP/PA coding is given for M0
		% Load the M0 JSON if existing
		if iReversedPE == 1
			pathM0In = fullfile(inSessionPath,'M0');
		else
			pathM0In = fullfile(inSessionPath,'M0PERev');
		end
			
		% Load the M0 JSON
		if xASL_exist([pathM0In '.json'])
			jsonM0 = spm_jsonread([pathM0In '.json']);
		else
			jsonM0 = struct;
		end
			
		if iReversedPE == 1
			if xASL_exist(fullfile(inSessionPath,'M0PERev.nii'))
				strPEDirection = '_dir-ap';

				jsonM0.PhaseEncodingDirection = 'j-';
				jsonLocal.PhaseEncodingDirection = 'j-';
			else
				strPEDirection = '';
			end
			if bJsonLocalM0isFile
				jsonLocal.M0 = [jsonLocal.M0 strPEDirection '_' bidsPar.strM0scan '.nii.gz'];
			end
			% Define the path to the respective ASL
			jsonM0.IntendedFor = [aslOutLabelRelative '_asl.nii.gz'];
			pathM0Out = fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectSessionLabel strPEDirection '_' bidsPar.strM0scan]);
		else
			jsonM0.PhaseEncodingDirection = 'j';
			jsonM0.IntendedFor = fullfile(bidsPar.strPerfusion,['sub-' subjectSessionLabel '_dir-ap' '_' bidsPar.strM0scan '.nii.gz']);
			pathM0Out = fullfile(outSessionPath,bidsPar.strFmap,['sub-' subjectSessionLabel '_dir-pa' '_' bidsPar.strM0scan]);
		end
		
		% Create the directory for the reversed PE if needed
		if iReversedPE == 2 && xASL_exist([pathM0In '.nii']) && ~exist(fullfile(outSessionPath,bidsPar.strFmap),'dir')
			mkdir(fullfile(outSessionPath,bidsPar.strFmap));
		end
		
		% If M0, then copy M0 and add ASL path to the IntendedFor
		if xASL_exist([pathM0In '.nii'])
			jsonM0 = xASL_bids_BIDSifyM0(jsonM0, jsonLocal, studyPar, pathM0In, pathM0Out, headerASL);
			
			% Save JSON to new dir
			jsonM0 = xASL_bids_VendorFieldCheck(jsonM0);
			jsonM0 = xASL_bids_JsonCheck(jsonM0,'M0');
			spm_jsonwrite([pathM0Out '.json'],jsonM0);
		end
	end
else
	if bJsonLocalM0isFile
		jsonLocal.M0 = [jsonLocal.M0 '.nii.gz'];
	end
end
	
	
	
    
	
%% 6. Save all ASL files (JSON, NIFTI, CONTEXT) to the BIDS directory
jsonLocal = xASL_bids_BIDSifyASLNII(jsonLocal, bidsPar, fullfile(inSessionPath,[aslLabel '.nii']), aslOutLabel);
jsonLocal = xASL_bids_VendorFieldCheck(jsonLocal);
[jsonLocal,bidsReport] = xASL_bids_JsonCheck(jsonLocal,'ASL');
spm_jsonwrite([aslOutLabel '_asl.json'],jsonLocal);

% Export report file for ASL dependencies
if exist('bidsReport','var')
    if ~isempty(fieldnames(bidsReport))
        spm_jsonwrite(fullfile(fileparts(imPar.BidsRoot),'bidsReportASL.json'),bidsReport);
    end
end

end
