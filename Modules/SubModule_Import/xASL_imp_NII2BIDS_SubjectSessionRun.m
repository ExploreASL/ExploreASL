function xASL_imp_NII2BIDS_SubjectSessionRun(imPar, bidsPar, studyPar, subjectSessionLabel, nameSubject, nameSession, inSessionPath, outSessionPath, listRuns, iRun)
%xASL_imp_NII2BIDS_SubjectSessionRun NII2BIDS conversion for a single sessions, single run.
%
% FORMAT: xASL_imp_NII2BIDS_SubjectSessionRun(imPar, bidsPar, studyPar, subjectSessionLabel, nameSubject, nameSession, inSessionPath, outSessionPath, listRuns, iRun)
% 
% INPUT:
% imPar               - JSON file with structure with import parameter (STRUCT, REQUIRED)
% bidsPar             - Output of xASL_imp_Config (STRUCT, REQUIRED)
% studyPar            - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
% subjectSessionLabel - subject-session label (CHAR ARRAY, REQUIRED)
% nameSubject         - name of the subject (CELL ARRAY, REQUIRED)
% nameSession         - name of the session (CELL ARRAY, REQUIRED)
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
jsonDicom = spm_jsonread(fullfile(inSessionPath, [aslLabel '.json']));
headerASL = xASL_io_ReadNifti(fullfile(inSessionPath, [aslLabel '.nii']));
	
%% 3. BIDSify ASL
% Merge the information from DICOM, manually entered parameters and BIDSify
jsonLocal = xASL_bids_BIDSifyASLJSON(jsonDicom, studyPar, headerASL);

	
	
	
	
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Define the M0 type
    [jsonLocal, bJsonLocalM0isFile] = xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, inSessionPath, subjectSessionLabel);
		
	if iRun == 1
        for nn = 1:2
            if nn == 1
                nnStrIn = '';
                if xASL_exist(fullfile(imPar.TempRoot,nameSubject,nameSession,'M0PERev.nii'))
                    nnStrOut = '_dir-ap';

                    tagPhaseEncodingDirection = 'j-';
                    jsonLocal.PhaseEncodingDirection = 'j-';
                    tagIntendedFor = [];
                    tagTotalReadoutTime = studyPar.TotalReadoutTime;

                    if bJsonLocalM0isFile
                        jsonLocal.M0 = [jsonLocal.M0 nnStrOut '_' bidsPar.strM0scan '.nii.gz'];
                    end
                else
                    if bJsonLocalM0isFile
                        jsonLocal.M0 = [jsonLocal.M0 '_' bidsPar.strM0scan '.nii.gz'];
                    end
                    nnStrOut = '';
                    tagPhaseEncodingDirection = [];
                    tagIntendedFor = [];
                    tagTotalReadoutTime = [];
                end
            else
                nnStrIn = 'PERev';
                nnStrOut = '_dir-pa';
                tagPhaseEncodingDirection = 'j';
                tagIntendedFor = fullfile(bidsPar.strPerfusion,['sub-' subjectSessionLabel '_dir-ap' '_' bidsPar.strM0scan '.nii.gz']);

                if isfield(studyPar,'TotalReadoutTime')
                    tagTotalReadoutTime = studyPar.TotalReadoutTime;
                else
                    tagTotalReadoutTime = [];
                end
            end

            % If M0, then copy M0 and add ASL path to the IntendedFor
            if xASL_exist(fullfile(imPar.TempRoot,nameSubject,nameSession,['M0' nnStrIn '.nii']))
                jsonM0 = spm_jsonread(fullfile(inSessionPath,['M0' nnStrIn '.json']));
                imM0   = xASL_io_Nifti2Im(fullfile(inSessionPath,['M0' nnStrIn '.json']));

                if ~isempty(regexpi(jsonDicom.Manufacturer,'Philips'))
                    jsonM0.scaleFactor = xASL_adm_GetPhilipsScaling(jsonM0,xASL_io_ReadNifti(fullfile(inSessionPath,['M0' nnStrIn '.nii'])));
                else
                    jsonM0.scaleFactor = 0;
                end

                if jsonM0.scaleFactor
                    imM0 = imM0 .* jsonM0.scaleFactor;
                end

                % Check echo time, for vectors
                if isfield(jsonM0,'EchoTime') && length(jsonM0.EchoTime)>1
                    % Remove zero entries
                    jsonM0.EchoTime = jsonM0.EchoTime(jsonM0.EchoTime ~= 0);
                end

                jsonM0Write = jsonM0;

                if isfield(jsonLocal,'SliceTiming')
                    % Issue a warning if the SliceTiming was already existing for M0, but still overwrite with ASL one
                    if isfield(jsonM0Write,'SliceTiming')
                        warning('SliceTiming already existed for M0, overwriting with ASL');
                    end

                    if headerASL.dat.dim(3) == size(imM0,3)
                        % Either copy if the save number of slices in M0 as in ASL
                        jsonM0Write.SliceTiming = jsonLocal.SliceTiming;
                    else
                        % Or recalculate for M0 if the number of slices differ
                        jsonM0Write.SliceTiming = ((0:(size(imM0,3)-1))')*(jsonLocal.SliceTiming(2)-jsonLocal.SliceTiming(1));
                    end
                else
                    if isfield(jsonM0Write,'SliceTiming')
                        jsonM0Write = rmfield(jsonM0Write,'SliceTiming');
                        warning('Removing pre-existing SliceTiming from M0, as there was no SliceTiming for ASL');
                    end
                end

                %if isfield(studyPar,'RepetitionTime')
                %	jsonM0Write.RepetitionTimePreparation = studyPar.RepetitionTime;
                %else
                jsonM0Write.RepetitionTimePreparation = jsonM0.RepetitionTime;
                %end

                jsonM0Write.IntendedFor = [aslOutLabelRelative '_asl.nii.gz'];

                if ~isempty(tagPhaseEncodingDirection)
                    jsonM0Write.PhaseEncodingDirection = tagPhaseEncodingDirection;
                end

                if ~isempty(tagIntendedFor)
                    jsonM0Write.IntendedFor = tagIntendedFor;
                end

                if ~isempty(tagTotalReadoutTime)
                    jsonM0Write.TotalReadoutTime = tagTotalReadoutTime;
                end

                if nn == 2 && ~exist(fullfile(outSessionPath,'fmap'),'dir')
                    mkdir(fullfile(outSessionPath,'fmap'));
                end

                % if scaling modified then save instead of move
                if jsonM0.scaleFactor || size(imM0,4) == 1
                    if nn == 1
                        xASL_io_SaveNifti(fullfile(inSessionPath,['M0' nnStrIn '.nii']),fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectSessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
                        % Delete original Nifti
                        xASL_delete(fullfile(inSessionPath,['M0' nnStrIn '.nii']));
                    else
                        xASL_io_SaveNifti(fullfile(inSessionPath,['M0' nnStrIn '.nii']),fullfile(outSessionPath,'fmap',['sub-' subjectSessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
                        % Delete original Nifti
                        xASL_delete(fullfile(inSessionPath,['M0' nnStrIn '.nii']));
                    end
                else
                    % Move the M0
                    if nn == 1
                        xASL_Move(fullfile(inSessionPath,['M0' nnStrIn '.nii']),...
                            fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectSessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),1);
                    else
                        xASL_Move(fullfile(inSessionPath,['M0' nnStrIn '.nii']),...
                            fullfile(outSessionPath,'fmap',['sub-' subjectSessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),1);
                    end
                end
                % Save JSON to new dir
                jsonM0Write = xASL_bids_VendorFieldCheck(jsonM0Write);
                jsonM0Write = xASL_bids_JsonCheck(jsonM0Write,'M0');
                if nn == 1
                    spm_jsonwrite(fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectSessionLabel nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
                else
                    spm_jsonwrite(fullfile(outSessionPath,'fmap',['sub-' subjectSessionLabel nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
                end
            end
        end
    else
        if bJsonLocalM0isFile
            jsonLocal.M0 = [jsonLocal.M0 '.nii.gz'];
        end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TODO
    % Remove the AslContext field and save it as a separate file
	filenameTSV = [aslOutLabel '_' bidsPar.strAslContext '.tsv'];
	[pathTSV,~,~] = fileparts(filenameTSV);
	if ~exist(pathTSV,'dir')
		mkdir(pathTSV);
	end
    fContext = fopen(filenameTSV,'w+');
    fwrite(fContext,sprintf('volume_type\n'));
    fwrite(fContext,jsonLocal.ASLContext);
    fclose(fContext);

    jsonLocal = rmfield(jsonLocal,'ASLContext');
	
	%%%%%%TODO%%%%% Save the JSON and NII to final location for ASL
	if jsonLocal.scaleFactor || length(headerASL.dat.dim) < 4 || headerASL.dat.dim(4) == 1
		imNii = xASL_io_Nifti2Im(fullfile(inSessionPath,[aslLabel '.nii']));
		
		if jsonLocal.scaleFactor
			imNii = imNii .* jsonLocal.scaleFactor;
		end
		
		% Scaling changed, so we have to save again OR
		% The fourth dimension is 1, so we have to write the file again, to make sure the
		xASL_io_SaveNifti(fullfile(inSessionPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],imNii,[],1,[]);
		% Delete original Nifti
		xASL_delete(fullfile(inSessionPath,[aslLabel '.nii']));
	else
		% Move the ASL
		xASL_Move(fullfile(inSessionPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],1);
	end
	
    %% Save all ASL files (JSON, NIFTI, CONTEXT) to the BIDS directory
    jsonLocal = xASL_bids_VendorFieldCheck(jsonLocal);
    jsonLocal = xASL_bids_JsonCheck(jsonLocal,'ASL');
    spm_jsonwrite([aslOutLabel '_asl.json'],jsonLocal);

end
