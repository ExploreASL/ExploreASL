function [imPar, bidsPar, studyPar, subjectLabel, sessionLabel, listSubjects, fSes, inSessionPath, outSessionPath, nSes, iSubject] = xASL_imp_NII2BIDS_SubjectSessionRun(imPar, bidsPar, studyPar, subjectLabel, sessionLabel, listSubjects, fSes, inSessionPath, outSessionPath, nSes, iSubject, kk, mm)
%xASL_imp_NII2BIDS_SubjectSessionRun NII2BIDS conversion for a single sessions, single run.
%
% FORMAT: [imPar, bidsPar, studyPar, subjectLabel, sessionLabel, listSubjects, fSes, inSessionPath, outSessionPath, nSes, iSubject] = xASL_imp_NII2BIDS_SubjectSessionRun(imPar, bidsPar, studyPar, subjectLabel, sessionLabel, listSubjects, fSes, inSessionPath, outSessionPath, nSes, iSubject, kk, mm)
% 
% INPUT:
% imPar            - JSON file with structure with import parameter (STRUCT, REQUIRED)
% bidsPar          - Output of xASL_imp_Config (STRUCT, REQUIRED)
% studyPar         - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
% subjectLabel     - subject label (CHAR ARRAY, REQUIRED)
% sessionLabel     - session label (CHAR ARRAY, REQUIRED)
% listSubjects     - list of subjects (CELL ARRAY, REQUIRED)
% fSes             - f session (CELL ARRAY, REQUIRED)
% inSessionPath    - input session path (CHAR ARRAY, PATH, REQUIRED)
% outSessionPath   - output session path (CHAR ARRAY, PATH, REQUIRED)
% nSes             - Number of sessions (INTEGER, REQUIRED)
% iSubject         - Subject ID (INTEGER, REQUIRED)
% kk               - session number (INTEGER, REQUIRED)
% mm               - run number (INTEGER, REQUIRED)
%
% OUTPUT:
% imPar            - imPar struct
% bidsPar          - bidsPar struct
% studyPar         - studyPar struct
% subjectLabel     - subject label
% sessionLabel     - session label
% listSubjects     - list of subjects
% fSes             - f session
% inSessionPath    - input session path
% outSessionPath   - output session path
% nSes             - Number of sessions
% iSubject         - Subject ID
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single sessions, single run.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Run
    if nSes
        aslLabel = ['ASL4D_' num2str(mm)];
        aslOutLabel = fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel '_run-' num2str(mm)]);
        aslOutLabelRelative = fullfile(bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel '_run-' num2str(mm)]);
    else
        aslLabel = 'ASL4D';
        aslOutLabel = fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel]);
        aslOutLabelRelative = fullfile(bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel]);
    end

    % Load the JSON
    jsonDicom = spm_jsonread(fullfile(inSessionPath,[aslLabel '.json']));
    imNii = xASL_io_Nifti2Im(fullfile(inSessionPath,[aslLabel '.nii']));

    if ~isempty(regexpi(jsonDicom.Manufacturer,'Philips'))
        scaleFactor = xASL_adm_GetPhilipsScaling(jsonDicom,xASL_io_ReadNifti(fullfile(inSessionPath,[aslLabel '.nii'])));
    else
        scaleFactor = 0;
    end

    % Check if the Phoenix protocol is present and parse it
    if isfield(jsonDicom,'PhoenixProtocol') && ~isempty(regexpi(jsonDicom.Manufacturer,'Siemens'))
        PhoenixParsed = xASL_bids_PhoenixProtocolReader(jsonDicom.PhoenixProtocol);
        jsonDicom.PhoenixAnalyzed = xASL_bids_PhoenixProtocolAnalyzer(PhoenixParsed);
    end

    % TotalAcquiredPairs - use only the maximum values from DICOM
    if isfield(jsonDicom,'TotalAcquiredPairs') && ~isempty(jsonDicom.TotalAcquiredPairs) && length(jsonDicom.TotalAcquiredPairs)>1
        jsonDicom.TotalAcquiredPairs = max(jsonDicom.TotalAcquiredPairs);
    end

    if scaleFactor
        imNii = imNii .* scaleFactor;
    end

    if scaleFactor || (size(imNii,4) == 1)
        % Scaling changed, so we have to save again OR
        % The fourth dimension is 1, so we have to write the file again, to make sure the
        xASL_io_SaveNifti(fullfile(inSessionPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],imNii,[],1,[]);
        % Delete original Nifti
        xASL_delete(fullfile(inSessionPath,[aslLabel '.nii']));
    else
        % Move the ASL
        xASL_Move(fullfile(inSessionPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],1);
    end

    % Take all the manually predefined fields from studyPar
    jsonLocal = studyPar;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In case of LookLocker and manually defined multiple flip angle, use this as a priority
    if isfield(studyPar,'LookLocker') && ~isempty(studyPar.LookLocker) && studyPar.LookLocker
        if isfield(studyPar,'FlipAngle') && ~isempty(studyPar.FlipAngle) && length(studyPar.FlipAngle)>1
            if isfield(jsonDicom,'FlipAngle') && ~isequal(jsonDicom.FlipAngle,jsonLocal.FlipAngle)
                jsonDicom.FlipAngle = jsonLocal.FlipAngle;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Overwrite differing fields with those from Dicom, but report all differences
    strDifferentFields = '';
    for fn = fieldnames(jsonDicom)'
        if isfield(jsonLocal,fn{1})
            % If the field is there, then report different fields
            if ~isequal(jsonLocal.(fn{1}),jsonDicom.(fn{1}))
                % Just if this is not only a different vector orientation
                if ~isnumeric(jsonLocal.(fn{1})) ||...
                        (size(jsonLocal.(fn{1}),1)>1 && size(jsonLocal.(fn{1}),1) >1) ||...
                        ~isequal((jsonLocal.(fn{1}))',jsonDicom.(fn{1}))
                    strDifferentFields = [strDifferentFields ' ' fn{1}];
                end
                if strcmp(fn{1},'TotalAcquiredPairs')
                    jsonDicom.(fn{1}) = jsonLocal.(fn{1});
                end
            end
        end
        % Prioritize the DICOM values in general case
        jsonLocal.(fn{1}) = jsonDicom.(fn{1});
    end
    % Report if certain fields were different as a warning
    if ~isempty(strDifferentFields)
        warning('The following user defined and DICOM fields differ:\n %s \n',strDifferentFields);
    end

    % Check repetition time
    if isfield(studyPar,'RepetitionTimePreparation')
        % RTP from studyPar has the highest priority
        jsonLocal.RepetitionTimePreparation = studyPar.RepetitionTimePreparation;
    elseif isfield(studyPar,'RepetitionTime')
        % RT from studyPar comes next
        jsonLocal.RepetitionTimePreparation = studyPar.RepetitionTime;
    elseif isfield(jsonDicom,'RepetitionTime')
        % RT from the DICOM has the lowest priority
        jsonLocal.RepetitionTimePreparation = jsonDicom.RepetitionTime;
        % Check if TR is a vector - replace by the maximum then
        if length(jsonLocal.RepetitionTimePreparation) > 1
            jsonLocal.RepetitionTimePreparation = max(jsonLocal.RepetitionTimePreparation);
            warning('TR was a vector. Taking the maximum only.');
        end
    end

    % Check echo time, for vectors
    if isfield(jsonLocal,'EchoTime') && length(jsonLocal.EchoTime)>1
        % Remove zero entries
        jsonLocal.EchoTime = jsonLocal.EchoTime(jsonLocal.EchoTime ~= 0);
    end

    % Convert the field name of the old LabelingType
    if isfield(jsonLocal,'LabelingType')
        if isfield(jsonLocal,'ArterialSpinLabelingType') && ~isequal(jsonLocal.LabelingType,jsonLocal.ArterialSpinLabelingType)
            warning('Both LabelingType and ArterialSpinLabelingType are defined with a different value');
        else
            jsonLocal.ArterialSpinLabelingType = jsonLocal.LabelingType;
        end
    end

    % The Labeling defined in a private GE field has a priority
    if isfield(jsonLocal,'GELabelingDuration') && ~isempty(jsonLocal.GELabelingDuration)
        if isfield(jsonLocal,'LabelingDuration') && ~isequal(jsonLocal.GELabelingDuration,jsonLocal.LabelingDuration)
            warning('Labeling duration mismatch with GE private field.');
        end
        jsonLocal.LabelingDuration = jsonLocal.GELabelingDuration;

        % GELabelingDuration comes together with the PostLabelingDelay defined in the standard DICOM field called InversionTime
        if isfield(jsonLocal,'InversionTime') && ~isempty(jsonLocal.InversionTime)
            % Verify if this doesn't differ from the predefined file, but the DICOM field has priority
            if isfield(jsonLocal,'PostLabelingDelay') && ~isequal(jsonLocal.PostLabelingDelay,jsonLocal.InversionTime)
                warning('PostLabelingDelay mismatch with the GE-DICOM value in Inversion time.');
            end
            jsonLocal.PostLabelingDelay = jsonLocal.InversionTime;
        end
    end

    % Free info about the sequence, now just the scanner type+software
    if isfield(jsonDicom,'ManufacturersModelName')
        jsonLocal.PulseSequenceDetails = jsonDicom.ManufacturersModelName;
    end
    if isfield(jsonDicom,'SoftwareVersions')
        if ~isfield(jsonLocal,'PulseSequenceDetails')
            jsonLocal.PulseSequenceDetails = '';
        end
        if ~isempty(jsonLocal.PulseSequenceDetails)
            jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails '-'];
        end
        jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails jsonDicom.SoftwareVersions];
    end

    % Process all the data and automatically fill in the missing parameters
    if ~isfield(jsonLocal,'MRAcquisitionType')
        error('MRAcquisitionType has to be defined either in the data or studyPar');
    end

    if strcmpi(jsonLocal.MRAcquisitionType,'2D')
        jsonLocal.PulseSequenceType = '2D_EPI';
    else
        if strcmpi(jsonLocal.Manufacturer,'GE') || strcmpi(jsonLocal.Manufacturer,'GE_WIP') || strcmpi(jsonLocal.Manufacturer,'GE_product')
            jsonLocal.PulseSequenceType = '3D_spiral';
        else
            jsonLocal.PulseSequenceType = '3D_GRASE';
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FME sequence check
    if isfield(jsonDicom,'SeriesDescription')
        isHadamardFME = ~isempty(regexp(char(jsonDicom.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
    else
        isHadamardFME = false;
    end
    if isHadamardFME
        if isfield(jsonLocal,'EchoTime') && isfield(jsonLocal,'PostLabelingDelay')
			% From the import, the length of EchoTime should correspond to the number of volumes
            if length(jsonLocal.EchoTime)~=length(jsonLocal.PostLabelingDelay)
                % So here, we first make sure that each PLD is repeated for the whole block of echo-times
				numberTEs = length(uniquetol(jsonLocal.EchoTime,0.001)); % Obtain the number of echo times
				repeatedPLDs = repmat(jsonLocal.PostLabelingDelay(:),1,numberTEs)';
				repeatedPLDs = repeatedPLDs(:);
                
				if length(repeatedPLDs) > length(jsonLocal.EchoTime) || mod(length(jsonLocal.EchoTime),length(repeatedPLDs)) ~= 0
					warning('Did not succeed in repeating PLDs for each TE for Hadamard sequence import');
				end
				% Make sure that number of volumes can be divided by the repeated PLDs
                jsonLocal.PostLabelingDelay = repeatedPLDs;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Overwrite differing fields with those from the Phoenix protocol, but report all differences
    if isfield(jsonDicom,'PhoenixAnalyzed') && ~isempty(jsonDicom.PhoenixAnalyzed)
        strDifferentFields = '';
        for fn = fieldnames(jsonDicom.PhoenixAnalyzed)'
            if isfield(jsonLocal,fn{1})
                % If the field is there, then report different fields
                if ~isequal(jsonLocal.(fn{1}),jsonDicom.PhoenixAnalyzed.(fn{1}))
                    % Just if this is not only a different vector orientation
                    if ~isnumeric(jsonLocal.(fn{1})) ||...
                            (size(jsonLocal.(fn{1}),1)>1 && size(jsonLocal.(fn{1}),1) >1) ||...
                            ~isequal((jsonLocal.(fn{1}))',jsonDicom.PhoenixAnalyzed.(fn{1}))
                        strDifferentFields = [strDifferentFields ' ' fn{1}];
                    end
                end
            end
            % Prioritize the Phoenix field values in the general case
            if ~(strcmp(fn{1},'EchoTime') && isHadamardFME)
                jsonLocal.(fn{1}) = jsonDicom.PhoenixAnalyzed.(fn{1});
            end
        end
        % Report if certain fields were different as a warning
        if ~isempty(strDifferentFields)
            warning('The following user-defined/DICOM fields and DICOM-Phoenix fields differ:\n %s \n',strDifferentFields);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BSup sanity check
	if ~isfield(jsonLocal,'BackgroundSuppression')
		warning('BackgroundSuppression field should be define in BIDS, setting the default to false');
		jsonLocal.BackgroundSuppression = false;
	end
	
	if jsonLocal.BackgroundSuppression == false
		% remove pulsenumbers and timings if BSup is OFF
		if isfield(jsonLocal,'BackgroundSuppressionNumberPulses')
			jsonLocal = rmfield(jsonLocal,'BackgroundSuppressionNumberPulses');
		end
		if isfield(jsonLocal,'BackgroundSuppressionPulseTime')
			jsonLocal = rmfield(jsonLocal,'BackgroundSuppressionPulseTime');
		end
	else
		% If times are given, but not the number of pulses, then assign the length
		if isfield(jsonLocal,'BackgroundSuppressionPulseTime')
			if isfield(jsonLocal,'BackgroundSuppressionNumberPulses')
				if jsonLocal.BackgroundSuppressionNumberPulses ~= length(jsonLocal.BackgroundSuppressionPulseTime)
					fprintf('Warning: Number of pulses and their timings do not match.');
				end
			else
				jsonLocal.BackgroundSuppressionNumberPulses = length(jsonLocal.BackgroundSuppressionPulseTime);
			end
		end
	end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SliceReadoutTime from the manual entry is prioritized
    if isfield (studyPar,'SliceReadoutTime')
        if isfield(studyPar,'SliceTiming') && ~isequal(studyPar.SliceTiming,studyPar.SliceReadoutTime)
            sprintf('Warning, difference in SliceTiming and SliceReadoutTime');
        end
        studyPar.SliceTiming = studyPar.SliceReadoutTime;
    end

    % Fill in extra parameters based on the JSON from the data
    if jsonLocal.MRAcquisitionType(1) == '2'
        % Take the studyPar as a prior source of SliceTiming since this is mostly wrong in DICOM otherwise
        if isfield(studyPar,'SliceTiming')
            jsonLocal.SliceTiming = studyPar.SliceTiming;
        end

        % The Siemens field is rather reliable though
        if isfield(jsonLocal,'SiemensSliceTime') && ~isempty(jsonLocal.SiemensSliceTime)
            jsonLocal.SliceTiming = jsonLocal.SiemensSliceTime;
        end

        % If the length of SliceTiming fits to the number of slices, do nothing
        if length(jsonLocal.SliceTiming) ~= size(imNii,3)
            % if the length of studyPar.sliceTiming is higher than 1 and the difference non-zero then use this
            if length(jsonLocal.SliceTiming) > 1 && abs(jsonLocal.SliceTiming(2)-jsonLocal.SliceTiming(1)) > 0
                jsonLocal.SliceTiming = jsonLocal.SliceTiming(2)-jsonLocal.SliceTiming(1);
            end

            if abs(jsonLocal.SliceTiming) > 0
                jsonLocal.SliceTiming = ((0:(size(imNii,3)-1))')*jsonLocal.SliceTiming;
            end
        end
    else
        % 3D sequences should not have a SliceTiming or have it defined as zero
        if isfield(jsonLocal,'SliceTiming')
            jsonLocal = rmfield(jsonLocal,'SliceTiming');
        end
    end

    % Import the number of averages
    if isfield(jsonLocal,'NumberOfAverages') && (max(jsonLocal.NumberOfAverages) > 1)
        if isfield(studyPar,'TotalAcquiredPairs')
            if max(jsonLocal.NumberOfAverages) ~= studyPar.TotalAcquiredPairs
                warning('Discrepancy in the number of averages');
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for FME Hadamard sequences
    if isfield(jsonLocal,'SeriesDescription')
        isHadamardFME = ~isempty(regexp(char(jsonLocal.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
        if isHadamardFME
            startDetails = regexp(char(jsonLocal.SeriesDescription),'\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once');
			jsonLocal.HadamardType = xASL_str2num(jsonLocal.SeriesDescription(startDetails:startDetails+1),'auto');
			jsonLocal.HadamardNumberPLD = xASL_str2num(jsonLocal.SeriesDescription(startDetails+3:startDetails+4),'auto');
            jsonLocal.HadamardNumberTE = xASL_str2num(jsonLocal.SeriesDescription(startDetails+10:startDetails+11),'auto');
            fprintf('FME sequence, Hadamard-%d encoded images, %d PLDs, %d TEs\n', jsonLocal.HadamardType, jsonLocal.HadamardNumberPLD, jsonLocal.HadamardNumberTE);
        end
    else
        isHadamardFME = false;
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Define the M0 type
    [studyPar, bidsPar, jsonLocal, inSessionPath, subjectLabel, sessionLabel, bJsonLocalM0isFile] = ...
        xASL_imp_NII2BIDS_Subject_DefineM0Type(studyPar, bidsPar, jsonLocal, inSessionPath, subjectLabel, sessionLabel);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If Post-labeling delay or labeling duration is longer than 1, but shorter then number of volumes then repeat it
    listFieldsRepeat = {'PostLabelingDelay', 'LabelingDuration','VascularCrushingVENC','FlipAngle','RepetitionTimePreparation'};
    for iRepeat = 1:length(listFieldsRepeat)
        if isfield(jsonLocal,(listFieldsRepeat{iRepeat})) && (length(jsonLocal.(listFieldsRepeat{iRepeat})) > 1) && (size(imNii,4) ~= length(jsonLocal.(listFieldsRepeat{iRepeat})))
			if mod(size(imNii,4),length(jsonLocal.(listFieldsRepeat{iRepeat})))
				error('Cannot find a match between the %s and the 4th dimension of the NIFTI.\n',listFieldsRepeat{iRepeat});
			else
				jsonLocal.(listFieldsRepeat{iRepeat}) = repmat(jsonLocal.(listFieldsRepeat{iRepeat})(:),[size(imNii,4)/length(jsonLocal.(listFieldsRepeat{iRepeat})) 1]);
			end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reformat ASLcontext field
    % Remove ',' and ';' at the end
    if (jsonLocal.ASLContext(end) == ';') || (jsonLocal.ASLContext(end) == ',')
        jsonLocal.ASLContext = jsonLocal.ASLContext(1:(end-1));
    end

    % Replace all ',' and ';' by \n
    jsonLocal.ASLContext = strrep(jsonLocal.ASLContext, ' ','');
    jsonLocal.ASLContext = strrep(jsonLocal.ASLContext, ';',',');
    lengthASLContext = sum(jsonLocal.ASLContext == ',')+1;
    jsonLocal.ASLContext = strrep(jsonLocal.ASLContext, ',',sprintf('\n'));

    % Check if the length is the same
    if size(imNii,4) ~= lengthASLContext
        % Check if we can simply repeat it
        if mod(size(imNii,4),lengthASLContext)
            error('Cannot find a match between the ASLContext and the 4th dimension of the NIFTI');
        else
            numRepeat = size(imNii,4)/lengthASLContext;
            tmpStr = jsonLocal.ASLContext;
            for iRepeat = 2:numRepeat
                jsonLocal.ASLContext = sprintf('%s\n%s',jsonLocal.ASLContext,tmpStr);
            end
        end
    end
    jsonLocal.ASLContext = sprintf('%s\n',jsonLocal.ASLContext);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If Post-labeling delay or labeling duration is longer than 1 then we need to verify this against
    % the ASLContext - these fields should be 0 for m0scans
    listFieldsZero = {'PostLabelingDelay', 'LabelingDuration'};

    % Find m0scans indices in ASLContext
    ASLContextCell = strsplit(jsonLocal.ASLContext,'\n'); % Split to cells by line-end
    ASLContextM0Index = regexp(ASLContextCell,'^m0scan'); % Find m0scans
    ASLContextM0Index = cellfun(@(x)~isempty(x),ASLContextM0Index); % Create a vector out of it
    ASLContextM0Index = ASLContextM0Index(1:size(imNii,4)); % Remove the last empty field

    % Go through all variables, check those that have length bigger than 1
    for iRepeat = 1:length(listFieldsZero)
        if isfield(jsonLocal,listFieldsZero{iRepeat}) && length(jsonLocal.(listFieldsZero{iRepeat})) > 1
            % Make sure the vector is a row vector
            jsonLocal.(listFieldsZero{iRepeat}) = jsonLocal.(listFieldsZero{iRepeat})(:)';
            % Check the all m0scans have zeros
            if ~isequal(ASLContextM0Index,jsonLocal.(listFieldsZero{iRepeat})==0)
                % If not, then set to zeros and report a warning
                jsonLocal.(listFieldsZero{iRepeat})(ASLContextM0Index) = 0;
                warning(['Had to set non-zero values for m0scan to zero in ' listFieldsZero{iRepeat}]);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Verify that TotalAcquiredPairs is a reasonable number

    % First count the number of controls, labels, and deltaMs
    ASLContextControlIndex = cellfun(@(x)~isempty(x),regexpi(ASLContextCell,'^control')); % Create a vector out of it
    ASLContextControlIndex = ASLContextControlIndex(1:size(imNii,4)); % Remove the last empty field
    ASLContextLabelIndex = cellfun(@(x)~isempty(x),regexpi(ASLContextCell,'^label')); % Create a vector out of it
    ASLContextLabelIndex = ASLContextLabelIndex(1:size(imNii,4)); % Remove the last empty field
    ASLContextDeltaMIndex = cellfun(@(x)~isempty(x),regexpi(ASLContextCell,'^deltam')); % Create a vector out of it
    ASLContextDeltaMIndex = ASLContextDeltaMIndex(1:size(imNii,4)); % Remove the last empty field

    % If TotalAcquiredPairs is 1, but more control/label pairs od deltaMs are present, then set this to the correct
    % number
    if ~isfield(jsonLocal,'TotalAcquiredPairs') || jsonLocal.TotalAcquiredPairs == 1
        if sum(ASLContextDeltaMIndex) > 0
            if sum(ASLContextControlIndex) == 0 && sum(ASLContextControlIndex) == 0
                jsonLocal.TotalAcquiredPairs = sum(ASLContextDeltaMIndex);
            else
                warning('Cannot calculate TotalAcquiredPairs when both controls and deltaMs are present');
            end
        elseif sum(ASLContextControlIndex) > 0
            if sum(ASLContextControlIndex) == sum(ASLContextLabelIndex)
                jsonLocal.TotalAcquiredPairs = sum(ASLContextControlIndex);
            else
                warning('Cannot calculte TotalAcquiredPairs when control and label numbers differ');
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove the AslContext field and save it as a separate file
    fContext = fopen([aslOutLabel '_' bidsPar.strAslContext '.tsv'],'w+');
    fwrite(fContext,sprintf('volume_type\n'));
    fwrite(fContext,jsonLocal.ASLContext);
    fclose(fContext);

    jsonLocal = rmfield(jsonLocal,'ASLContext');

    if isfield(jsonLocal,'BolusCutOffFlag') && jsonLocal.BolusCutOffFlag
        if isfield(jsonLocal,'BolusCutOffTechnique') && strcmpi(jsonLocal.BolusCutOffTechnique,'Q2TIPS')
            if ~isfield(jsonLocal,'BolusCutOffDelayTime') || length(jsonLocal.BolusCutOffDelayTime)~=2
                warning('Q2TIPS BolusCutOff has to have 2 values defined');
            end
        end
    end

    if mm == 1
        for nn = 1:2
            if nn == 1
                nnStrIn = '';
                if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk},'M0PERev.nii'))
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
                tagIntendedFor = fullfile(bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel '_dir-ap' '_' bidsPar.strM0scan '.nii.gz']);

                if isfield(studyPar,'TotalReadoutTime')
                    tagTotalReadoutTime = studyPar.TotalReadoutTime;
                else
                    tagTotalReadoutTime = [];
                end
            end

            % If M0, then copy M0 and add ASL path to the IntendedFor
            if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk},['M0' nnStrIn '.nii']))
                jsonM0 = spm_jsonread(fullfile(inSessionPath,['M0' nnStrIn '.json']));
                imM0   = xASL_io_Nifti2Im(fullfile(inSessionPath,['M0' nnStrIn '.json']));

                if ~isempty(regexpi(jsonDicom.Manufacturer,'Philips'))
                    scaleFactor = xASL_adm_GetPhilipsScaling(jsonM0,xASL_io_ReadNifti(fullfile(inSessionPath,['M0' nnStrIn '.nii'])));
                else
                    scaleFactor = 0;
                end

                if scaleFactor
                    imM0 = imM0 .* scaleFactor;
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

                    if size(imNii,3) == size(imM0,3)
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
                if scaleFactor || size(imM0,4) == 1
                    if nn == 1
                        xASL_io_SaveNifti(fullfile(inSessionPath,['M0' nnStrIn '.nii']),fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
                        % Delete original Nifti
                        xASL_delete(fullfile(inSessionPath,['M0' nnStrIn '.nii']));
                    else
                        xASL_io_SaveNifti(fullfile(inSessionPath,['M0' nnStrIn '.nii']),fullfile(outSessionPath,'fmap',['sub-' subjectLabel sessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
                        % Delete original Nifti
                        xASL_delete(fullfile(inSessionPath,['M0' nnStrIn '.nii']));
                    end
                else
                    % Move the M0
                    if nn == 1
                        xASL_Move(fullfile(inSessionPath,['M0' nnStrIn '.nii']),...
                            fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),1);
                    else
                        xASL_Move(fullfile(inSessionPath,['M0' nnStrIn '.nii']),...
                            fullfile(outSessionPath,'fmap',['sub-' subjectLabel sessionLabel nnStrOut '_' bidsPar.strM0scan '.nii.gz']),1);
                    end
                end
                % Save JSON to new dir
                jsonM0Write = xASL_bids_VendorFieldCheck(jsonM0Write);
                jsonM0Write = xASL_bids_JsonCheck(jsonM0Write,'M0');
                if nn == 1
                    spm_jsonwrite(fullfile(outSessionPath,bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
                else
                    spm_jsonwrite(fullfile(outSessionPath,'fmap',['sub-' subjectLabel sessionLabel nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
                end
            end
        end
    else
        if bJsonLocalM0isFile
            jsonLocal.M0 = [jsonLocal.M0 '.nii.gz'];
        end
    end
    % Save JSON to new dir
    jsonLocal = xASL_bids_VendorFieldCheck(jsonLocal);
    jsonLocal = xASL_bids_JsonCheck(jsonLocal,'ASL');
    spm_jsonwrite([aslOutLabel '_asl.json'],jsonLocal);


end



