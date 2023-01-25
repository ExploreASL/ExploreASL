function xASL_imp_NII2BIDS_RunPerf(imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun)
%xASL_imp_NII2BIDS_RunPerf NII2BIDS conversion for a single sessions, single run.
%
% FORMAT: xASL_imp_NII2BIDS_RunPerf(imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun)
% 
% INPUT:
% imPar               - JSON file with structure with import parameter (STRUCT, REQUIRED)
% bidsPar             - Output of xASL_imp_Config (STRUCT, REQUIRED)
% studyPar            - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
% subjectSessionLabel - subject-session label (CHAR ARRAY, REQUIRED)
% inSessionPath       - input session path (CHAR ARRAY, PATH, REQUIRED)
% outSessionPath      - output session path (CHAR ARRAY, PATH, REQUIRED)
% listRuns            - list of runs (INTEGER, REQUIRED)
% iRun                - run number (INTEGER, REQUIRED)
%
% OUTPUT:
% Files               - Saves bids_report and conditionally ASL & M0.json files
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
% EXAMPLE:     xASL_imp_NII2BIDS_RunPerf(imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun);
%
% __________________________________
% Copyright 2015-2022 ExploreASL

    %% 1. Define the pathnames
	if (length(listRuns)>1) || (str2num(listRuns{1}(5))>1)
		aslLegacyLabel = 'ASL4D';
		subjectSessionLabel = ['sub-' subjectSessionLabel];
        runLabel = ['_run-' num2str(iRun)];
	else
		aslLegacyLabel = 'ASL4D';
		subjectSessionLabel = ['sub-' subjectSessionLabel];
        runLabel = '';
	end
	aslOutLabel = fullfile(outSessionPath,bidsPar.stringPerfusion, [subjectSessionLabel runLabel]);
	aslOutLabelRelative = fullfile(bidsPar.stringPerfusion, [subjectSessionLabel runLabel]);
    
    % Current scan name
    [~, scanName] = xASL_fileparts([aslOutLabel '_' bidsPar.stringASL '.json']);
    
    % Print perfusion name
    fprintf('scan %s ...\n', scanName);

    %% 2. Load the JSONs and NIfTI information
    if exist(fullfile(inSessionPath, [aslLegacyLabel '.json']),'file')
        jsonDicom = xASL_io_ReadJson(fullfile(inSessionPath, [aslLegacyLabel '.json']));
    else
        warning('Missing JSON file: %s\n',fullfile(inSessionPath, [aslLegacyLabel '.json']));
		
		% Try to check if a file with a similar name is there
		jsonDicom = xASL_adm_GetFileList(inSessionPath, ['^' aslLegacyLabel '.*\.json$']);
		if ~isempty(jsonDicom)
			fprintf('File with an incorrect but similar name was found: %s \n',jsonDicom{1});
		end
		
        return;
    end
    if xASL_exist(fullfile(inSessionPath, [aslLegacyLabel '.nii']),'file')
        headerASL = xASL_io_ReadNifti(fullfile(inSessionPath, [aslLegacyLabel '.nii']));
		
		% For GE scanner and 3D ASL file, we try to reformat to 4D correctly
		if length(headerASL.dat.dim) == 3 && isfield(jsonDicom,'Manufacturer') && strcmp(jsonDicom.Manufacturer,'GE') &&...
			isfield(jsonDicom,'NumberOfTemporalPositions') && jsonDicom.NumberOfTemporalPositions > 2 && ...
			mod(headerASL.dat.dim(3),jsonDicom.NumberOfTemporalPositions) == 0
		
			fprintf('GE ASL data incorrectly read as 3D. Reshaping...\n');
			% Load the image data and reshape accordingly
		    imASL = xASL_io_Nifti2Im(fullfile(inSessionPath, [aslLegacyLabel '.nii']));
			imASL = reshape(imASL, size(imASL,1), size(imASL,2), size(imASL,3)/jsonDicom.NumberOfTemporalPositions, jsonDicom.NumberOfTemporalPositions);
			
			% Center the new image correctly
			newMat = headerASL.mat0;
			newMat(1:3,4) = [0;-22;0] - headerASL.mat0(1:3,1:3)*([size(imASL,1); size(imASL,2); size(imASL,3)]/2);
			
			% However, the matrix is given with a wrong shift
			xASL_io_SaveNifti(fullfile(inSessionPath, [aslLegacyLabel '.nii']), fullfile(inSessionPath, [aslLegacyLabel '.nii']), imASL, [], [], newMat);
			headerASL = xASL_io_ReadNifti(fullfile(inSessionPath, [aslLegacyLabel '.nii']));
		end
    else
        warning('Missing NIfTI file: %s\n\', fullfile(inSessionPath, [aslLegacyLabel '.nii']));
		
		% Try to check if a file with a similar name is there
		headerASL = xASL_adm_GetFileList(inSessionPath, ['^' aslLegacyLabel '.*\.nii$']);
		if ~isempty(headerASL)
			fprintf('File with an incorrect but similar name was found: %s \n', headerASL{1});
		end
		
        return;
	end
	
    %% 3. BIDSify ASL
    % Merge the information from DICOM, manually entered parameters and BIDSify
    jsonLocal = xASL_bids_BIDSifyASLJSON(jsonDicom, studyPar, headerASL);

    %% 4. Prepare the link to M0 in ASL.json	
    % Define the M0 type	
    [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(...
        studyPar, bidsPar, jsonLocal, fullfile(inSessionPath,'M0.nii'), fullfile(bidsPar.stringPerfusion,[subjectSessionLabel runLabel]));	

    %% 5. BIDSify M0	
    % Check the M0 files 
    
	for iReversedPE = 1:2 % Check if AP/PA coding is given for M0
		% Load the M0 JSON if existing
		if iReversedPE == 1
			pathM0In = fullfile(inSessionPath,'M0');
		else
			pathM0In = fullfile(inSessionPath,'M0PERev');
		end
		
		% Load the M0 JSON
		if xASL_exist([pathM0In '.json'])
			jsonM0 = xASL_io_ReadJson([pathM0In '.json']);
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
				jsonLocal.M0 = [jsonLocal.M0 strPEDirection '_' bidsPar.stringM0scan '.nii.gz'];
            end
            
			% Define the path to the respective ASL
			jsonM0.IntendedFor = [aslOutLabelRelative '_' bidsPar.stringASL '.nii.gz'];
            
            % Determine output name
            aslLegacyLabel = 'ASL4D';
            bidsm0scanLabel = [subjectSessionLabel strPEDirection runLabel '_' bidsPar.stringM0scan];
			pathM0Out = fullfile(outSessionPath,bidsPar.stringPerfusion,bidsm0scanLabel);
		else
			jsonM0.PhaseEncodingDirection = 'j';
            strPEDirectionPA = '_dir-pa';
            strPEDirectionAP = '_dir-ap';
            % Determine output name
            aslLegacyLabel = 'ASL4D';
            bidsm0scanLabelPA = [subjectSessionLabel strPEDirectionPA runLabel '_' bidsPar.stringM0scan];
            bidsm0scanLabelAP = [subjectSessionLabel strPEDirectionAP runLabel '_' bidsPar.stringM0scan];
            jsonM0.IntendedFor = fullfile(bidsPar.stringPerfusion,[bidsm0scanLabelAP '.nii.gz']);
            pathM0Out = fullfile(outSessionPath,bidsPar.stringFmap,bidsm0scanLabelPA);
		end
		
		% Create the directory for the reversed PE if needed
		if iReversedPE == 2 && xASL_exist([pathM0In '.nii'])
			xASL_adm_CreateDir(fullfile(outSessionPath,bidsPar.stringFmap));
		end
		
		% If M0, then copy M0 and add ASL path to the IntendedFor
		if xASL_exist([pathM0In '.nii'])
			jsonM0 = xASL_bids_BIDSifyM0(jsonM0, jsonLocal, studyPar, pathM0In, pathM0Out, headerASL);
			
			% Save JSON to new dir
			jsonM0 = xASL_bids_VendorFieldCheck(jsonM0);
			jsonM0 = xASL_bids_JsonCheck(jsonM0,'M0');
			xASL_io_WriteJson([pathM0Out '.json'],jsonM0);
		end
	end



    %% 6. Save all ASL files (JSON, NIFTI, CONTEXT) to the BIDS directory
    jsonLocal = xASL_bids_BIDSifyASLNII(jsonLocal, bidsPar, fullfile(inSessionPath,[aslLegacyLabel '.nii']), aslOutLabel);
    jsonLocal = xASL_bids_VendorFieldCheck(jsonLocal);
    [jsonLocal,bidsReport] = xASL_bids_JsonCheck(jsonLocal,'ASL');
    xASL_io_WriteJson([aslOutLabel '_' bidsPar.stringASL '.json'],jsonLocal);

    % Export report file for ASL dependencies
    if exist('bidsReport','var')
        if ~isempty(fieldnames(bidsReport))
            xASL_io_WriteJson(fullfile(fileparts(imPar.BidsRoot),'derivatives','ExploreASL','log', ...
                ['bids_report_' subjectSessionLabel runLabel '.json']), bidsReport);
        end
    end


end


