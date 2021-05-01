function [parms, pathDcmDictOut] = xASL_bids_Dicom2JSON(imPar, pathIn, pathJSON, dcmExtFilter, bUseDCMTK, pathDcmDictIn)
%xASL_bids_Dicom2JSON Go through the DICOM or PAR/REC files, parses the header and saves it in JSON
%
% FORMAT: [parms pathDcmDictOut] = xASL_bids_Dicom2Parms(imPar, pathIn[, pathJSON, dcmExtFilter, bUseDCMTK, pathDcmDictIn])
%
% INPUT:
%        imPar               - struct with import parameters (REQUIRED)
%        pathIn (PATH)       - path to the RAW files (REQUIRED)
%        pathJSON (PATH)     - cell with paths to the JSON file for saving parameters (OPTIONAL, DEFAULT = don't save)
%        dcmExtFilter (STR)  - wildcards specifying the allowed extensions for the RAW files
%        bUseDCMTK (BOOL)    - if yes, then use DCMTK instead of dicominfo
%        pathDcmDictIn (STR) - path to the dicom dictionary in case DCMTK fails and DICOMINFO is used
%
% OUTPUT:
%        parms               - structure containing the parsed parameters
%        pathDcmDictOut      - if dicom dict for dicominfo is initialized then clear this path, otherwise return unchanged pathDcmDictIn
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      The function goes through the pathIn files, reads the DICOM or PAR/REC files 
%                   and parses their headers. It extracts the DICOM parameters important for ASL, 
%                   makes sure they are in the correct format, if missing then replaces with default 
%                   value, it also checks if the parameters are consistent across DICOM files for a 
%                   single sequence.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          ...
%
% REFERENCES:
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% ----------------------------------------------------------------------------------
	% Admin
	% ----------------------------------------------------------------------------------
	
	if nargin<2 || isempty(pathJSON)
		pathJSON = cell(1,1);
	end
	if nargin<3 || isempty(dcmExtFilter)
		dcmExtFilter='^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$'; % the last one is because some convertors save files without extension, but there would be a dot/period before a bunch of numbers
	end
	if nargin<4 || isempty(bUseDCMTK)
		bUseDCMTK = true; % use this by default
    elseif ~bUseDCMTK && isempty(which('dicomdict'))
        error('Dicomdict missing, image processing probably not installed, try DCMTK instead');
	end
	if nargin<5
		pathDcmDictIn = [];
	end
	
	pathDcmDictOut = pathDcmDictIn;
    
	%% ----------------------------------------------------------------------------------
	% Set up the default values
	% ----------------------------------------------------------------------------------
	
	DcmParDefaults.RepetitionTime               = NaN;
	DcmParDefaults.EchoTime                     = NaN;
	DcmParDefaults.NumberOfAverages             = 1;   % no temporal positions in 3D, as default for non-Philips scan. CAVE!!!!
	DcmParDefaults.NumberOfTemporalPositions    = 1;   % no temporal positions in 3D, as default for non-Philips scan. CAVE!!!!
	DcmParDefaults.RescaleSlope                 = 1;   % RescaleSlope; added by Paul to get rid of misleading RescaleSlopeOriginal
	DcmParDefaults.RescaleSlopeOriginal         = NaN; % RescaleSlopeOriginal; will be set to RescaleSlope if missing
	DcmParDefaults.MRScaleSlope                 = 1;   % MRScaleSlope
	DcmParDefaults.RescaleIntercept             = 0;   % RescaleIntercept (although this one is standard dicom)
	DcmParDefaults.AcquisitionTime              = 0;   % AcquisitionTime
    
	% TopUp parameters
	DcmParDefaults.AcquisitionMatrix            = NaN;
	DcmParDefaults.EffectiveEchoSpacing         = NaN;
	DcmParDefaults.AssetRFactor                 = NaN;
	DcmParDefaults.MRSeriesWaterFatShift        = NaN;
	DcmParDefaults.MRSeriesEPIFactor            = NaN;
	DcmParDefaults.BandwidthPerPixelPhaseEncode = NaN;
	DcmParDefaults.RWVIntercept                 = NaN;
	DcmParDefaults.RWVSlope                     = NaN;
	
    % Image parameters
	DcmParDefaults.Rows                         = NaN;
	DcmParDefaults.Columns                      = NaN;
	DcmParDefaults.TemporalPositionIdentifier   = NaN;
	DcmParDefaults.PhilipsNumberTemporalScans   = NaN;
	DcmParDefaults.GELabelingDuration           = NaN;
	DcmParDefaults.InversionTime                = NaN;
	
	DcmSkipNan = {'Rows' 'Columns' 'TemporalPositionIdentifier' 'PhilipsNumberTemporalScans' ...
		          'GELabelingDuration' 'InversionTime' 'RWVIntercept' 'RWVSlope'};
	
	DcmComplexFieldFirst = {'PulseSequenceName' 'GELabelingType'  'SiemensSliceTime' 'PhoenixProtocol' 'InPlanePhaseEncodingDirection'};
	DcmComplexFieldAll = {'ComplexImageComponent' 'AcquisitionContrast' 'ImageType' 'PhilipsLabelControl'};
	
	DcmFieldList = {'RepetitionTime', 'EchoTime', 'NumberOfAverages', 'RescaleSlope', ...
					'RescaleSlopeOriginal', 'MRScaleSlope', 'RescaleIntercept', 'AcquisitionTime', ...
					'AcquisitionMatrix', 'Rows', 'Columns', 'NumberOfAverages', 'NumberOfTemporalPositions', ...
                    'RWVIntercept', 'RWVSlope'};
	
	bVendor = 'Unknown';
	
	%% ----------------------------------------------------------------------------------
	% Recreate the parameter file from raw data
	% ----------------------------------------------------------------------------------
	for iJSON = 1:length(pathJSON)
		if ~isempty(pathJSON{iJSON})
			if imPar.bVerbose; fprintf('Recreating parameter file: %s\n',pathJSON{iJSON}); end
			
			% Make a list of the instanceNumbers
			tmpJSON = spm_jsonread(pathJSON{iJSON});
			if isfield(tmpJSON,'InstanceNumber')
				instanceNumberList(iJSON) = xASL_str2num(tmpJSON.InstanceNumber);
			else
				instanceNumberList(iJSON) = 0;
			end
		end
		% Create a cell of parms
		parms{iJSON} = struct();
	end
	
	% Reoder JSONs by increasing instanceNumber to allow easy categorizing to the correct range.
	[instanceNumberList,sortInstance] = sort(instanceNumberList,'ascend');
	pathJSON = pathJSON(sortInstance);
	
	if exist(pathIn, 'dir')
		FileList            = xASL_adm_GetFileList(pathIn, dcmExtFilter, 'FPList', [0 Inf]); % we assume all the dicoms are in the same folder
		for iF=1:length(FileList)
			[~, fname, ext] = fileparts(FileList{iF});
			FileList{iF}        = [fname ext];
		end
	else
		[pathIn, fname, ext] = fileparts(pathIn);
		FileList = {[fname ext]};
	end
	
    TryDicominfo = true; % this is only set to false below upon succesful DcmtkRead
    
    if ~isempty(FileList)
        iMrFileAll = zeros(length(parms));
        
        % Check all dicom files throughout the sequence to validate that they are the same - because we are fast now!
        nFiles  = length(FileList);
        
        for iFile = 1:nFiles
            if TryDicominfo && iFile>1
                continue;
                % with dicominfo, reading is very slow, so we only read 1 dicom
            elseif ismac() && iFile>1
                % spm_jsonread is very slow, which is used in
                % xASL_io_DmtkRead for macOS, so we also only read one file
                % for macOS
                continue;
			else
				if nFiles > 300
					if mod(iFile,50) == 0
						xASL_TrackProgress(iFile, nFiles);
					end
				else
					xASL_TrackProgress(iFile, nFiles);
				end
            end
            
			ifname = FileList{iFile};
			filepath = fullfile(pathIn, ifname); % this is a file by definition, according to the xASL_adm_GetFileList command above
			
			%% ----------------------------------------------------------------------------------
			% Use DCMTK library to read the DICOM header to temp
			% ----------------------------------------------------------------------------------
            
            if bUseDCMTK
                try
                    temp = xASL_io_DcmtkRead(filepath, 0);
                    TryDicominfo = false;
                catch ME
                    warning(['xASL_adm_Dicom2JSON: xASL_io_DcmtkRead failed for ' filepath ', trying dicominfo']);
                    fprintf('%s\n', ['Message: ' ME.message]);
                end
            end
            
			if TryDicominfo
				% The dicom-dict has not yet been initialized and set
				if ~isempty(pathDcmDictIn)
					dicomdict('set', pathDcmDictIn);
					pathDcmDictOut = '';
				end
				
				DictionaryDCM = dicomdict('get_current'); % save current dictionary
				
				try % 1) use current dictionary
					temp = dicominfo(filepath);
                catch ME
					try % 2) use default dictionary
                        fprintf('%s\n', ['Warning: ' ME.message]);
						dicomdict('factory'); % temporarily set to default dictionary
						temp = dicominfo(filepath); % retry reading dicom
                    catch ME
						try % 3) change UseDictionaryVR setting with current dictionary
							fprintf('%s\n', ['Warning: ' ME.message]);
                            dicomdict('set',DictionaryDCM);
							temp = dicominfo(filepath,'UseDictionaryVR', true);
                        catch ME
							try % 4) change UseDictionaryVR setting with default dictionary
                                fprintf('%s\n', ['Warning: ' ME.message]);
								dicomdict('factory');
								temp = dicominfo(filepath,'UseDictionaryVR', true);
                            catch ME
								warning('dicominfo also did not work, check this!');
                                fprintf('%s\n', ['Message: ' ME.message]);
								dicomdict('set',DictionaryDCM); % reset dictionary
								continue;
							end
						end
					end
				end
				
				dicomdict('set',DictionaryDCM); % reset dictionary
			end
			
			%% -----------------------------------------------------------------------------
			% Obtain the instance number and JSON index
			%% -----------------------------------------------------------------------------
			if isfield(temp,'InstanceNumber') && ~isempty(temp.InstanceNumber)
				currentInstanceNumber = temp.InstanceNumber;
			else
				currentInstanceNumber = 0;
			end
			
			parmsIndex = 1;
			for iInst = 1:length(instanceNumberList)
				if currentInstanceNumber > instanceNumberList(iInst)
					parmsIndex = iInst;
				end
			end
			
			%% -----------------------------------------------------------------------------
            % Take information from enhanced DICOM, if exists
            %% -----------------------------------------------------------------------------
			if isfield(temp, 'MediaStorageSOPClassUID')
				if strcmp(temp.MediaStorageSOPClassUID,'1.2.840.10008.5.1.4.1.1.4.1')==1 % Enhanced MR Image Storage
					bEnhancedMR = true;
					iMrFileAll(parmsIndex) = iMrFileAll(parmsIndex)+1;
					iMrFile = iMrFileAll(parmsIndex);
				elseif strcmp(temp.MediaStorageSOPClassUID,'1.2.840.10008.5.1.4.1.1.4')==1 % MR Image Storage
					bEnhancedMR = false; % default
					iMrFileAll(parmsIndex) = iMrFileAll(parmsIndex)+1;
					iMrFile = iMrFileAll(parmsIndex);
				else
% 					continue; % THIS SEEMS STRANGE >>>>>>>>>> RESULTS IN EMPTY PARMS, BY TRYING TO READ IMRFILE =0
                    iMrFileAll(parmsIndex) = iMrFileAll(parmsIndex)+1;
					iMrFile = iMrFileAll(parmsIndex);
                    bEnhancedMR = false; % default
                end
			else
				bEnhancedMR = false; % default
				iMrFileAll(parmsIndex) = iMrFileAll(parmsIndex)+1;
				iMrFile = iMrFileAll(parmsIndex);
			end
			
			% Deal with enhanced DICOM format imported through DICOMINFO and not DCMTK
            if bEnhancedMR && ~TryDicominfo
                % fprintf('Enhanced DICOM detected, but no dicominfo selected, skipping obtaining parameters\n');
            elseif bEnhancedMR && TryDicominfo
                % for simplicity, take the first value from the enhanced
                % sequences and store them in the temp-struct as if it is a
                % classic structure.

                if ~isfield(temp, 'PerFrameFunctionalGroupsSequence')
                    warning('xASL_adm_Dicom2JSON: Enhanced DICOM but PerFrameFunctionalGroupsSequence not found');
                else
                    itemNames = fieldnames(temp.PerFrameFunctionalGroupsSequence);
                    SequenceFields = temp.PerFrameFunctionalGroupsSequence.(itemNames{1}); % here we assume the same ScaleSlope for all images in the sequence

                    % Fill temp dicominfo with these fields
                    Fields1Are = fields(SequenceFields);
                    for iField1=1:length(Fields1Are)
                        if ischar(SequenceFields.(Fields1Are{iField1}))
                            temp.(Fields1Are{iField1}) = SequenceFields.(Fields1Are{iField1});
                        elseif isstruct(SequenceFields.(Fields1Are{iField1}))
                            Fields2Are = fields(SequenceFields.(Fields1Are{iField1}));
                            SequenceFields2 = SequenceFields.(Fields1Are{iField1}).(Fields2Are{1});

                            Fields3Are = fields(SequenceFields2);
                            for iField3=1:length(Fields3Are)
                                temp.(Fields3Are{iField3}) = SequenceFields2.(Fields3Are{iField3});
                            end
                        end
                    end
                end
            end
			
			
			% Do this only once, and do not reset the manufacturer for other files (assume the same is for all files within the directory)
			if ~exist('dcmfields','var')
				%% -----------------------------------------------------------------------------
				% Identify vendor and select the important header parameters to read
				% -----------------------------------------------------------------------------
				if  isfield(temp, 'Manufacturer')
					manufacturer    = lower(temp.Manufacturer);
					if ~isempty(strfind(manufacturer,'ge'))
						bVendor = 'GE';
					elseif ~isempty(strfind(manufacturer,'philips'))
						bVendor = 'Philips';
					elseif ~isempty(strfind(manufacturer,'siemens'))
						bVendor = 'Siemens';
					else
						warning('xASL_adm_Dicom2JSON: Manufacturer unknown for %s', filepath);
					end
				else
					warning('xASL_adm_Dicom2JSON: Manufacturer unknown for %s', filepath);
				end
				
				dcmfields = DcmFieldList;
				
				switch bVendor
					case 'GE'
						dcmfields(end+1:end+4) = {'AssetRFactor', 'EffectiveEchoSpacing'...
							'GELabelingDuration' 'InversionTime' }; % (0043,1083) (0043,102c)
					case 'Philips'
						dcmfields(end+1:end+4) = {'MRSeriesWaterFatShift', 'MRSeriesEPIFactor'...
							'TemporalPositionIdentifier'  'PhilipsNumberTemporalScans'}; % (2001,1022) (2001,1013)
					case 'Siemens'
						dcmfields(end+1) = {'BandwidthPerPixelPhaseEncode'}; % (0019,1028)
					otherwise
						% do nothing
				end
			end
			
			%% -----------------------------------------------------------------------------
			% Obtain the selected DICOM parameters from the header
			% Write the new parameter to the list (or put the default value)
			% -----------------------------------------------------------------------------
			for iField=1:length(dcmfields)
				fieldname = dcmfields{iField};
				
                % Extract value of current field
                if  isfield(temp, fieldname)
                    thevalue = temp.(fieldname);
                    thevalue = thevalue(~isnan(thevalue));
                else
                    thevalue = [];
                end
				
                % Convert string-numbers to numbers if necessary
                if ~isempty(thevalue)
                    if length(xASL_str2num(thevalue))>1
                        tmpTheValue = nonzeros(thevalue);
                        tmpTheValue = tmpTheValue(1);
                    else
                        tmpTheValue = thevalue;
                    end
                    t_parms{parmsIndex}.(fieldname)(iMrFile) = xASL_str2num(tmpTheValue);
                else
					if imPar.bVerbose
						if iMrFile==1, fprintf('%s\n',['Parameter ' fieldname ' not found, default used']); end
					end
					t_parms{parmsIndex}.(fieldname)(iMrFile) = DcmParDefaults.(fieldname);
                end
            end
			
			c_all_parms{parmsIndex} = struct;
			% The more complex fields - strings and arrays are saved in cell
			for iField=1:length(DcmComplexFieldAll)
				if isfield(temp,DcmComplexFieldAll{iField}) && ~isempty(temp.(DcmComplexFieldAll{iField}))
					c_all_parms{parmsIndex}.(DcmComplexFieldAll{iField}){iMrFile} = temp.(DcmComplexFieldAll{iField});
				end
			end
			
			c_first_parms{parmsIndex} = struct;
			for iField=1:length(DcmComplexFieldFirst)
				if isfield(temp,DcmComplexFieldFirst{iField}) && ~isempty(temp.(DcmComplexFieldFirst{iField}))
					c_first_parms{parmsIndex}.(DcmComplexFieldFirst{iField}) = temp.(DcmComplexFieldFirst{iField});
				end
			end
		end
		for indexInstance = 1:length(parms)
			parmsIndex = 1;
			for iInst = 1:length(instanceNumberList)
				if instanceNumberList(indexInstance) > 0
					parmsIndex = iInst;
				end
			end
			%% If no files were found previously (just directories etc.) then the manufacturer won't be identified and
			% dcmfields won't be assigned
			if exist('dcmfields','var')
				% -----------------------------------------------------------------------------
				% Dealing with empty or inconsistent field values
				% -----------------------------------------------------------------------------
				
				% Limit AcquisitionTime to one value
				if  isfield(t_parms{parmsIndex},'AcquisitionTime')
					for iL=1:length(t_parms{parmsIndex})
						t_parms{parmsIndex}(iL).AcquisitionTime = xASL_str2num(t_parms{parmsIndex}(iL).AcquisitionTime);
					end
					
					if length(t_parms{parmsIndex})>1
						if imPar.bVerbose; fprintf('%s\n','Parameter AcquisitionTime has multiple values:'); end
						%                     t_parms{parmsIndex}.AcquisitionTime
						if imPar.bVerbose; fprintf('%s\n','using minumum value:'); end
						%                     t_parms{parmsIndex}.AcquisitionTime = min(t_parms{parmsIndex}.AcquisitionTime)
						tempAcquisitionTime = t_parms{parmsIndex}(1).AcquisitionTime;
						t_parms{parmsIndex} = rmfield(t_parms{parmsIndex},'AcquisitionTime');
						t_parms{parmsIndex}(1).AcquisitionTime = tempAcquisitionTime;
					end
                end
				
                % Convert number to time format (check that current struct exists / is not empty)
                if ~isempty(t_parms{parmsIndex})
                    t_parms{parmsIndex}(1).AcquisitionTime = xASL_adm_ConvertNr2Time(t_parms{parmsIndex}(1).AcquisitionTime);
                else
                    fprintf('\nWarning: t_parms{parmsIndex} is empty...\n');
                end
				
				% Checks if the field values were the same for all dicoms and keep only one from the same value
				for iField=1:length(dcmfields)
					fieldname = dcmfields{iField};
					if  isfield(t_parms{parmsIndex},fieldname)
						parms{parmsIndex}.(fieldname) = t_parms{parmsIndex}.(fieldname);
						% Remove (set to NaN) also those that differ only minimally
						if length(parms{parmsIndex}.(fieldname))>1 && isnumeric(parms{parmsIndex}.(fieldname))
							iDiff = abs(parms{parmsIndex}.(fieldname) - parms{parmsIndex}.(fieldname)(1))./parms{parmsIndex}.(fieldname)(1);
							iDiff(1) = 1;
							parms{parmsIndex}.(fieldname)(iDiff<0.001) = NaN;
						end
						% There's one or more NaNs
						nNaN = sum(isnan(parms{parmsIndex}.(fieldname)));
						if nNaN > 0
							% Only NaNs
							if nNaN == length(parms{parmsIndex}.(fieldname))
								parms{parmsIndex}.(fieldname) = NaN;
							else
								parms{parmsIndex}.(fieldname) = parms{parmsIndex}.(fieldname)(~isnan(parms{parmsIndex}.(fieldname)));
							end
						end
					end
				end
			end
			
			% Remove fields that are NaN
			for iField=1:length(DcmSkipNan)
				if isfield(parms,DcmSkipNan{iField}) && sum(isnan(parms{parmsIndex}.(DcmSkipNan{iField})))
					parms{parmsIndex} = rmfield(parms,DcmSkipNan{iField});
				end
			end
			
			% The more complex fields - strings and arrays are saved in cell
			for iField=1:length(DcmComplexFieldAll)
				if isfield(c_all_parms{parmsIndex},DcmComplexFieldAll{iField})
					listEmptyFields = find(cellfun(@isempty,c_all_parms{parmsIndex}.(DcmComplexFieldAll{iField})));
					if ~isempty(listEmptyFields)
						fprintf('Field %s contains empty fields, skipping\n',DcmComplexFieldAll{iField});
					else
						c_all_unique = unique(c_all_parms{parmsIndex}.(DcmComplexFieldAll{iField}));
						if length(c_all_unique) == 1
							parms{parmsIndex}.(DcmComplexFieldAll{iField}) = c_all_unique;
						else
							parms{parmsIndex}.(DcmComplexFieldAll{iField}) = c_all_parms{parmsIndex}.(DcmComplexFieldAll{iField});
						end
					end
				end
			end
			
			for iField=1:length(DcmComplexFieldFirst)
				if isfield(c_first_parms{parmsIndex},DcmComplexFieldFirst{iField}) && ~isempty(c_first_parms{parmsIndex}.(DcmComplexFieldFirst{iField}))
					parms{parmsIndex}.(DcmComplexFieldFirst{iField}) = c_first_parms{parmsIndex}.(DcmComplexFieldFirst{iField});
				end
			end
			
			% Check for valid RescaleSlope value
			if isfield(parms{parmsIndex},'RescaleSlopeOriginal') && max(isnan(parms{parmsIndex}.RescaleSlopeOriginal))
				parms{parmsIndex}.RescaleSlopeOriginal = parms{parmsIndex}.RescaleSlope;
			end
			
			
			%% -----------------------------------------------------------------------------
			% Check and purge the parameters
			% -----------------------------------------------------------------------------
			% Check whether multiple scale slopes exist, this happens in 1
			% Philips software version, and should give an error (later we can
			% make this a warning, or try to deal with this)
			
			
			
			%% Calculate the TopUp parameters
			if isfield(parms{parmsIndex},'AcquisitionMatrix')
				parms{parmsIndex}.AcquisitionMatrix = double(parms{parmsIndex}.AcquisitionMatrix(1));
			end
			
			switch bVendor
				case 'GE'
					if isfield(parms{parmsIndex},'AssetRFactor') && isfield(parms{parmsIndex},'EffectiveEchoSpacing') && isfield(parms{parmsIndex},'AcquisitionMatrix')
						parms{parmsIndex}.EffectiveEchoSpacing = double(parms{parmsIndex}.EffectiveEchoSpacing);
						parms{parmsIndex}.AssetRFactor = double(parms{parmsIndex}.AssetRFactor);
						parms{parmsIndex}.AcquisitionMatrix = double(parms{parmsIndex}.AcquisitionMatrix);
						
						parms{parmsIndex}.EffectiveEchoSpacing = (parms{parmsIndex}.EffectiveEchoSpacing .* parms{parmsIndex}.AssetRFactor) / 10^6;
						parms{parmsIndex}.TotalReadoutTime = (parms{parmsIndex}.AcquisitionMatrix-1) .* parms{parmsIndex}.EffectiveEchoSpacing;
						parms{parmsIndex} = rmfield(parms{parmsIndex}, 'AssetRFactor');
					end
				case 'Philips'
					if isfield(parms{parmsIndex},'MRSeriesWaterFatShift') && isfield(parms{parmsIndex},'MRSeriesEPIFactor') && isfield(parms{parmsIndex},'AcquisitionMatrix')
						parms{parmsIndex}.MRSeriesWaterFatShift = double(parms{parmsIndex}.MRSeriesWaterFatShift);
						parms{parmsIndex}.MRSeriesEPIFactor = double(parms{parmsIndex}.MRSeriesEPIFactor);
						parms{parmsIndex}.AcquisitionMatrix = double(parms{parmsIndex}.AcquisitionMatrix);
						
						EffectiveEchoSpacingPhilips = parms{parmsIndex}.MRSeriesWaterFatShift/(434.215 * (parms{parmsIndex}.MRSeriesEPIFactor+1));
						parms{parmsIndex}.TotalReadoutTime = EffectiveEchoSpacingPhilips*(parms{parmsIndex}.MRSeriesEPIFactor-1);
						parms{parmsIndex}.EffectiveEchoSpacing = parms{parmsIndex}.TotalReadoutTime/(parms{parmsIndex}.AcquisitionMatrix(1)-1);
						parms{parmsIndex} = rmfield(parms{parmsIndex}, {'MRSeriesWaterFatShift' 'MRSeriesEPIFactor'});
					end
				case 'Siemens'
					if isfield(parms{parmsIndex},'BandwidthPerPixelPhaseEncode') && (~isnan(parms{parmsIndex}.BandwidthPerPixelPhaseEncode)) && isfield(parms{parmsIndex},'InPlanePhaseEncodingDirection')
						parms{parmsIndex}.BandwidthPerPixelPhaseEncode = double(parms{parmsIndex}.BandwidthPerPixelPhaseEncode);
						
						if isfield(parms{parmsIndex},'AcquisitionMatrix') && ~isempty(parms{parmsIndex}.AcquisitionMatrix) && ~sum(isnan(parms{parmsIndex}.AcquisitionMatrix))
							if length(parms{parmsIndex}.AcquisitionMatrix) == 1
								parms{parmsIndex}.ReconMatrixPE = parms{parmsIndex}.AcquisitionMatrix;
							elseif strcmp(parms{parmsIndex}.InPlanePhaseEncodingDirection,'COL')
								parms{parmsIndex}.ReconMatrixPE = parms{parmsIndex}.AcquisitionMatrix(2);
							else
								parms{parmsIndex}.ReconMatrixPE = parms{parmsIndex}.AcquisitionMatrix(1);
							end
						elseif strcmp(parms{parmsIndex}.InPlanePhaseEncodingDirection,'COL')
							parms{parmsIndex}.ReconMatrixPE = double(parms{parmsIndex}.Rows);
						elseif strcmp(parms{parmsIndex}.InPlanePhaseEncodingDirection,'ROW')
							parms{parmsIndex}.ReconMatrixPE = double(parms{parmsIndex}.Columns);
						else
							error('Unknown InPlanePhaseEncodingDirection');
						end
						
						parms{parmsIndex}.EffectiveEchoSpacing = 1/(parms{parmsIndex}.BandwidthPerPixelPhaseEncode*parms{parmsIndex}.ReconMatrixPE);
						parms{parmsIndex}.TotalReadoutTime = (parms{parmsIndex}.ReconMatrixPE-1) * parms{parmsIndex}.EffectiveEchoSpacing;
					end
				otherwise
					% skip
			end
			
			%% First remove non-finite values
            if  isfield(parms{parmsIndex},'EchoTime')
                parms{parmsIndex}.EchoTime              = parms{parmsIndex}.EchoTime(isfinite(parms{parmsIndex}.EchoTime));
            end
            if  isfield(parms{parmsIndex},'RepetitionTime')
                parms{parmsIndex}.RepetitionTime        = parms{parmsIndex}.RepetitionTime(isfinite(parms{parmsIndex}.RepetitionTime));
            end
			if  isfield(parms{parmsIndex},'MRScaleSlope')
				parms{parmsIndex}.MRScaleSlope          = parms{parmsIndex}.MRScaleSlope(isfinite(parms{parmsIndex}.MRScaleSlope));
            end
			if  isfield(parms{parmsIndex},'RescaleSlopeOriginal')
				parms{parmsIndex}.RescaleSlopeOriginal  = parms{parmsIndex}.RescaleSlopeOriginal(isfinite(parms{parmsIndex}.RescaleSlopeOriginal));
            end
			if  isfield(parms{parmsIndex},'RescaleIntercept')
				parms{parmsIndex}.RescaleIntercept      = parms{parmsIndex}.RescaleIntercept(isfinite(parms{parmsIndex}.RescaleIntercept));
			end
			
			% In case more than one value is given, then keep only the value that is not equal to 1. Or set to 1 if all are 1
			parmNameToCheck = {'MRScaleSlope','RescaleSlopeOriginal','RescaleSlope','RWVSlope'};
			for parmNameInd = 1:length(parmNameToCheck)
				parmName = parmNameToCheck{parmNameInd};
				if  isfield(parms,parmName) && (length(parms{parmsIndex}.(parmName))>1)
					
					indNonOne = find(parms{parmsIndex}.(parmName)~=1);
					if isempty(indNonOne)
						parms{parmsIndex}.(parmName) = 1;
					else
						parms{parmsIndex}.(parmName) = parms{parmsIndex}.(parmName)(indNonOne);
					end
					
				end
			end
			
			% In case multiple different scale slopes are given, report a warning
            if isfield(parms{parmsIndex},'MRScaleSlope') && isfield(parms{parmsIndex},'RescaleSlopeOriginal') && isfield(parms{parmsIndex},'RescaleSlope')
                if length(parms{parmsIndex}.MRScaleSlope)>1  || length(parms{parmsIndex}.RescaleSlopeOriginal)>1 || length(parms{parmsIndex}.RescaleSlope)>1 ||...
                        (isfield(parms{parmsIndex},'RWVSlope') && length(parms{parmsIndex}.RWVSlope)>1)
                    warning('xASL_adm_Dicom2JSON: Multiple scale slopes exist for a single scan!');
                    %parms{parmsIndex} = rmfield(parms{parmsIndex},'MRScaleSlope');
                    %parms{parmsIndex} = rmfield(parms{parmsIndex},'RescaleSlope');
                    %parms{parmsIndex} = rmfield(parms{parmsIndex},'RescaleSlopeOriginal');
                    %parms{parmsIndex} = rmfield(parms{parmsIndex},'RescaleIntercept');
                    %parms{parmsIndex} = rmfield(parms{parmsIndex},'RWVSlope');
                end
            else
                fprintf('Warning: MRScaleSlope, RescaleSlopeOriginal or RescaleSlope not found...\n');
            end
			
			%% Save the info in JSON file
			
			% Loads the JSON parms
			if exist(pathJSON{parmsIndex},'file')
				JSONParms = spm_jsonread(pathJSON{parmsIndex});
			else
				JSONParms = [];
			end
						
			% Merges them with parms
			parms{parmsIndex} = xASL_bids_parms2BIDS(parms{parmsIndex}, JSONParms, 1, 0);
			
			% Saves the JSON file
			spm_jsonwrite(pathJSON{parmsIndex}, parms{parmsIndex});
		end
    end


end

%% Obtain a value from a structure 
function val = GetDicomValue(I, fieldname, default)
if nargin>2
	val = default;
else
	val = [];
end

if isfield(I, fieldname)
	val = I.(fieldname);
end
end
