function xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, listSubjects, iSubject)
%xASL_imp_NII2BIDS_Subject Run NII to ASL-BIDS for one individual subject.
%
% FORMAT: xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, listSubjects, iSubject)
% 
% INPUT:
%   imPar           - JSON file with structure with import parameter (REQUIRED, STRUCT)
%   bidsPar         - Output of xASL_imp_Config (REQUIRED, STRUCT)
%   studyPar        - JSON file with the BIDS parameters relevant for the whole study (REQUIRED, STRUCT)
%   listSubjects    - List of subjects (REQUIRED, LIST)
%   iSubject        - Current subject number (REQUIRED, INTEGER)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run NII to ASL-BIDS for one individual subject.
%
% 1. Initialize
% 2. Process all the anatomical files
% 3. Process the perfusion files
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, listSubjects, iSubject);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1. Initialize
    
    subjectLabel = xASL_adm_CorrectName(listSubjects{iSubject},2);

	% Make a subject directory
	if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel]),'dir')
		mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel]));
	end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2. Process all the anatomical files
    % Go throught the list of anat files
    for iAnatType = bidsPar.listAnatTypes

        % Check if it exists
        anatPath = '';
        if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},[iAnatType{1},'.nii']),'file')
            anatPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},iAnatType{1});
        end

        if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},[iAnatType{1} '_1'],[iAnatType{1},'.nii']),'file')
            anatPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},[iAnatType{1} '_1'],iAnatType{1});
        end

        % If anatomical file of this type exist, then BIDSify its structure
        if ~isempty(anatPath)

            % Create the anatomical directory
            if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'),'dir')
                mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'));
            end

            % Move the NiFTI file
            xASL_Move([anatPath '.nii'],fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat',...
                ['sub-' subjectLabel '_' iAnatType{1} '.nii.gz']),1);

            % Load the JSON
            jsonAnat = spm_jsonread([anatPath,'.json']);

			% If RepetitionTimePreparation is equal to RepetitionTime, then remove RepetitionTimePreparation
			if isfield(jsonAnat,'RepetitionTime') && isfield(jsonAnat,'RepetitionTimePreparation') &&...
					isequal(jsonAnat.RepetitionTime,jsonAnat.RepetitionTimePreparation)
				jsonAnat = rmfield(jsonAnat,'RepetitionTimePreparation');
			end
			
            % Save the JSON
            jsonAnat = xASL_bids_VendorFieldCheck(jsonAnat);
            jsonAnat = xASL_bids_JsonCheck(jsonAnat,'');
            spm_jsonwrite(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat',['sub-' subjectLabel '_' iAnatType{1} '.json']),jsonAnat);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3. Process the perfusion files
    fSes = xASL_adm_GetFileList(fullfile(imPar.AnalysisRoot,listSubjects{iSubject}),'^ASL.+$',false,[],true);

    % Go through all sessions
    for kk = 1:length(fSes)

        % Make a subject directory
        if length(fSes)>1
            sessionLabel = ['ses-' fSes{kk}(5:end)];

            if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel),'dir')
                mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel));
                mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel,'asl'));
            end
            inSessionPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk});
            outSessionPath = fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel);

            % Need to add the underscore so that it doesn't need to be added automatically and can be skipped for empty session
            sessionLabel = ['_' sessionLabel];
        else
            % Session label is skipped
            sessionLabel = '';

            % Only one session - no session labeling
            if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel]),'dir')
                mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel]));
                mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],bidsPar.strPerfusion));
            end
            inSessionPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk});
            outSessionPath = fullfile(imPar.BidsRoot,['sub-' subjectLabel]);
        end

        % Check if there are multiple runs per session
        fRuns = xASL_adm_GetFileList(inSessionPath,'^ASL4D_\d.nii+$',false,[],false);
        nSes = length(fRuns);

        for mm = 1:(max(nSes,1))
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
                    jsonLocal.(fn{1}) = jsonDicom.PhoenixAnalyzed.(fn{1});
                end
                % Report if certain fields were different as a warning
                if ~isempty(strDifferentFields)
                    warning('The following user-defined/DICOM fields and DICOM-Phoenix fields differ:\n %s \n',strDifferentFields);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BSup sanity check
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Define the M0 type
            % Type of an M0 image
            bJsonLocalM0isFile = 0;
            if ~isfield(studyPar,'M0') || isempty(studyPar.M0) || strcmpi(studyPar.M0,'separate_scan')
                if isfield(studyPar,'M0PositionInASL4D') && (max(studyPar.M0PositionInASL4D(:))>0)
                    jsonLocal.M0 = true;
                    jsonLocal.M0Type = bidsPar.strM0Included;
                elseif xASL_exist(fullfile(inSessionPath,'M0.nii'))
                    jsonLocal.M0 = fullfile(bidsPar.strPerfusion,['sub-' subjectLabel sessionLabel]);
                    jsonLocal.M0Type = bidsPar.strM0Separate;
                    bJsonLocalM0isFile = 1;
                else
                    if ~isempty(strfind(jsonLocal.ASLContext,bidsPar.strM0scan))
                        jsonLocal.M0 = true;
                        jsonLocal.M0Type = bidsPar.strM0Included;
                    else
                        jsonLocal.M0 = false;
                        jsonLocal.M0Type = bidsPar.strM0Absent;
                    end
                end
            else
                if strcmpi(studyPar.M0,'UseControlAsM0')
                    jsonLocal.M0 = bidsPar.strM0Absent;
                else
                    if strcmpi(studyPar.M0,'no_background_suppression')
                        jsonLocal.M0 = bidsPar.strM0Absent;
                    else
                        jsonLocal.M0 = studyPar.M0;
                        if isnumeric(studyPar.M0)
                            jsonLocal.M0Type = bidsPar.strM0Estimate;
                            jsonLocal.M0Estimate = studyPar.M0;
                        elseif xASL_exist(fullfile(inSessionPath,'M0.nii'))
                            jsonLocal.M0Type = bidsPar.strM0Separate;
                        elseif ~isempty(strfind(jsonLocal.ASLContext,bidsPar.strM0scan))
                            jsonLocal.M0Type = bidsPar.strM0Included;
                        else
                            jsonLocal.M0Type = bidsPar.strM0Absent;
                        end
                    end
                end
            end

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
    end
end


