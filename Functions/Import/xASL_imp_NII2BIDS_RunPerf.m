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
% Copyright 2015-2023 ExploreASL

    %% 1. Define the pathnames
	if length(listRuns)>1 || xASL_str2num(listRuns{1}(5:end))>1
		aslLegacyLabel = 'ASL4D';
		subjectSessionLabel = ['sub-' subjectSessionLabel];
        runLabel = ['_run-' listRuns{iRun}(5:end)];
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
    jsonASL = xASL_bids_BIDSifyASLJSON(jsonDicom, studyPar, headerASL);

    %% 4. Prepare the link to M0 in ASL.json	
    % Define the M0 type	
    [jsonASL, bjsonASLM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(...
        studyPar, bidsPar, jsonASL, fullfile(inSessionPath,'M0.nii'), fullfile(bidsPar.stringPerfusion,[subjectSessionLabel runLabel]));	

    %% 5. BIDSify M0	
    % Check the M0 files 
    
	for iReversedPE = 1:2 % Check if AP/PA coding is given for M0
		% Obtain the path to M0
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
			% Check first the normal M0 without reversed PE direction
			if xASL_exist(fullfile(inSessionPath, 'M0PERev.nii'))
				% Check if the encoding direction was specified in DICOM
				if ~isfield(jsonM0, 'PhaseEncodingDirection') && ~isfield(jsonASL, 'PhaseEncodingDirection')
					% If none of the fields exist in DICOM, then assign the default
					jsonM0.PhaseEncodingDirection = 'j-';
					jsonASL.PhaseEncodingDirection = 'j-';
					fprintf('Phase-encoding direction for ASL and M0 is not specified, using the default AP direction.\n');
				elseif isfield(jsonM0, 'PhaseEncodingDirection') && isfield(jsonASL, 'PhaseEncodingDirection')
					% Both ASL and M0 have the phase-encoding direction field in DICOM, check if this is consistent
					if ~strcmp(jsonM0.PhaseEncodingDirection, jsonASL.PhaseEncodingDirection)
						warning('Phase-encoding direction differs between ASL and M0');
					end
				elseif isfield(jsonM0, 'PhaseEncodingDirection')
					% Only M0 has the field
					jsonASL.PhaseEncodingDirection = jsonM0.PhaseEncodingDirection;
					fprintf('Phase-encoding direction was not specified for ASL, using the same value as for M0.\n');
				else
					% Only ASL has the field
					jsonM0.PhaseEncodingDirection = jsonASL.PhaseEncodingDirection;
					fprintf('Phase-encoding direction was not specified for M0, using the same value as for ASL.\n');
				end
				
				% Define the direction tag for the BIDS field name based on the DICOM value
				switch (jsonM0.PhaseEncodingDirection)
					case 'i'
						strPEDirection = '_dir-rl';
					case 'i-'
						strPEDirection = '_dir-lr';
					case 'j'
						strPEDirection = '_dir-pa';
					case 'j-'
						strPEDirection = '_dir-ap';
					case 'k'
						strPEDirection = '_dir-si';
					case 'k-'
						strPEDirection = '_dir-is';
				end

			else
				strPEDirection = '';
			end
			if bjsonASLM0isFile
				jsonASL.M0 = [jsonASL.M0 strPEDirection '_' bidsPar.stringM0scan '.nii.gz'];
            end
            
			% Define the path to the respective ASL
			jsonM0.IntendedFor = [aslOutLabelRelative '_' bidsPar.stringASL '.nii.gz'];
            
            % Determine output name
            aslLegacyLabel = 'ASL4D';
            bidsm0scanLabel = [subjectSessionLabel strPEDirection runLabel '_' bidsPar.stringM0scan];
			pathM0Out = fullfile(outSessionPath,bidsPar.stringPerfusion,bidsm0scanLabel);
		else
			% The value is missing for the M0 with reversed PE, we assign the same default as for the ASL and M0 scans
			if ~isfield(jsonM0, 'PhaseEncodingDirection')
				jsonM0.PhaseEncodingDirection = 'j';
				fprintf('Phase-encoding direction for reversed-PE M0 is not specified, using the default PA direction.\n');
			end

			% Define the direction tag and inversed-direction tag for the BIDS field name based on the DICOM value
			switch(jsonM0.PhaseEncodingDirection)
				case 'j'
					strPEDirectionNorm    = '_dir-pa';
					strPEDirectionReverse = '_dir-ap';
				case 'j-'
					strPEDirectionNorm    = '_dir-ap';
					strPEDirectionReverse = '_dir-pa';
				case 'i'
					strPEDirectionNorm    = '_dir-rl';
					strPEDirectionReverse = '_dir-lr';
				case 'i-'
					strPEDirectionNorm    = '_dir-lr';
					strPEDirectionReverse = '_dir-rl';
				case 'k'
					strPEDirectionNorm    = '_dir-si';
					strPEDirectionReverse = '_dir-is';
				case 'k-'
					strPEDirectionNorm    = '_dir-is';
					strPEDirectionReverse = '_dir-si';
			end

            % Determine output name
            aslLegacyLabel = 'ASL4D';
            bidsm0scanLabelNorm    = [subjectSessionLabel strPEDirectionNorm    runLabel '_' bidsPar.stringM0scan];
            bidsm0scanLabelReverse = [subjectSessionLabel strPEDirectionReverse runLabel '_' bidsPar.stringM0scan];
            jsonM0.IntendedFor = fullfile(bidsPar.stringPerfusion, [bidsm0scanLabelReverse '.nii.gz']);
            pathM0Out = fullfile(outSessionPath, bidsPar.stringFmap, bidsm0scanLabelNorm);
		end
		
		% Create the directory for the reversed PE if needed
		if iReversedPE == 2 && xASL_exist([pathM0In '.nii'])
			xASL_adm_CreateDir(fullfile(outSessionPath,bidsPar.stringFmap));
		end
		
		% If M0, then copy M0 and add ASL path to the IntendedFor
		if xASL_exist([pathM0In '.nii'])
			jsonM0 = xASL_bids_BIDSifyM0(jsonM0, jsonASL, studyPar, pathM0In, pathM0Out, headerASL);
			
			% Save JSON to new dir
			jsonM0 = xASL_bids_VendorFieldCheck(jsonM0);
			jsonM0 = xASL_bids_JsonCheck(jsonM0,'M0');
			xASL_io_WriteJson([pathM0Out '.json'],jsonM0);
		end
	end



    %% 6. Save all ASL files (JSON, NIFTI, CONTEXT) to the BIDS directory
    jsonASL = xASL_bids_BIDSifyASLNII(jsonASL, bidsPar, fullfile(inSessionPath,[aslLegacyLabel '.nii']), aslOutLabel);
    jsonASL = xASL_bids_VendorFieldCheck(jsonASL);
    [jsonASL,bidsReport] = xASL_bids_JsonCheck(jsonASL,'ASL');
    xASL_io_WriteJson([aslOutLabel '_' bidsPar.stringASL '.json'],jsonASL);

    % Export report file for ASL dependencies
    if exist('bidsReport','var')
        if ~isempty(fieldnames(bidsReport))
            xASL_io_WriteJson(fullfile(fileparts(imPar.BidsRoot),'derivatives','ExploreASL','log', ...
                ['bids_report_' subjectSessionLabel runLabel '.json']), bidsReport);
        end
    end


end
