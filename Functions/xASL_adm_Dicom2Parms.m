function [parms, pathDcmDictOut] = xASL_adm_Dicom2Parms(imPar, inp, parmsfile, dcmExtFilter, bUseDCMTK, pathDcmDictIn)
% Goes through the DICOM or PAR/REC files, parses the header and saves it in MAT.
%
% FORMAT: [parms pathDcmDictOut] = xASL_adm_Dicom2Parms(inp[, parmsfile, dcmExtFilter, bUseDCMTK, pathDcmDictIn])
%
% INPUT:
%        inp (PATH)         - path to the RAW files
%        parmsfile (PATH)   - path to the MAT file for saving parsed parameters, if empty or not given, then don't save
%        dcmExtFilter (STR) - wildcards specifying the allowed extensions for the RAW files
%        bUseDCMTK (BOOL)   - if yes, then use DCMTK instead of dicominfo
%        pathDcmDictIn (STR)- path to the dicom dictionary in case DCMTK fails and DICOMINFO is used
% OUTPUT:
%        parms              - structure containing the parsed parameters
%        pathDcmDictOut     - if dicom dict for dicominfo is initialized then clear this path, otherwise return unchanged pathDcmDictIn
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
%        The function goes through the INP files, reads the DICOM or PAR/REC files and parses their headers.
%        It extracts the DICOM parameters important for ASL, makes sure they are in the correct format, if missing then 
%        replaces with default value, it also checks if the parameters are consistent across DICOM files for a single sequence.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
% __________________________________
% Copyright @ 2015-2019 ExploreASL

    %% ----------------------------------------------------------------------------------
	% Admin
	% ----------------------------------------------------------------------------------
	
	if nargin<2 || isempty(parmsfile)
		parmsfile = [];
	end
	if nargin<3 || isempty(dcmExtFilter)
		dcmExtFilter='^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$'; % the last one is because some convertors save files without extension, but there would be a dot/period before a bunch of numbers
	end
	if nargin<4 || isempty(bUseDCMTK)
		bUseDCMTK = true; % use this by default
    elseif bUseDCMTK && isempty(which('dicomdict'))
        error('Dicomdict missing, image processing probably not installed, try DCMTK instead');
	end
	if nargin<5
		pathDcmDictIn = [];
	end
	
	pathDcmDictOut = pathDcmDictIn;
	
	%% ----------------------------------------------------------------------------------
	% Set up the default values
	% ----------------------------------------------------------------------------------
	
	DcmParDefaults.RepetitionTime            = NaN;
	DcmParDefaults.EchoTime                  = NaN;
	DcmParDefaults.NumberOfAverages          = 1;   % no temporal positions in 3D, as default for non-Philips scan. CAVE!!!!
	DcmParDefaults.NumberOfTemporalPositions = 1;   % no temporal positions in 3D, as default for non-Philips scan. CAVE!!!!
	DcmParDefaults.RescaleSlope              = 1;   % RescaleSlope; added by Paul to get rid of misleading RescaleSlopeOriginal
	DcmParDefaults.RescaleSlopeOriginal      = NaN; % RescaleSlopeOriginal; will be set to RescaleSlope if missing
	DcmParDefaults.MRScaleSlope              = 1;   % MRScaleSlope
	DcmParDefaults.RescaleIntercept          = 0;   % RescaleIntercept (although this one is standard dicom)
	DcmParDefaults.AcquisitionTime           = 0;   % AcquisitionTime
	% TopUp parameters
	DcmParDefaults.AcquisitionMatrix         = NaN;
	DcmParDefaults.EffectiveEchoSpacing      = NaN;
	DcmParDefaults.AssetRFactor  = NaN;
	DcmParDefaults.MRSeriesWaterFatShift = NaN;
	DcmParDefaults.MRSeriesEPIFactor = NaN;
	DcmParDefaults.BandwidthPerPixelPhaseEncode = NaN;
	
	bVendor = 'Unknown';
	
	%% ----------------------------------------------------------------------------------
	% Recreate the parameter file from raw data
	% ----------------------------------------------------------------------------------
	if ~isempty(parmsfile)
		if imPar.bVerbose; fprintf('Recreating parameter file: %s\n',parmsfile); end
	end
	parms = struct();
	
	if exist(inp, 'dir')
		FileList            = xASL_adm_GetFileList(inp, dcmExtFilter, 'FPList', [0 Inf]); % we assume all the dicoms are in the same folder
		for iF=1:length(FileList)
			[~, fname, ext] = fileparts(FileList{iF});
			FileList{iF}        = [fname ext];
		end
	else
		[inp, fname, ext] = fileparts(inp);
		FileList = {[fname ext]};
	end
	
	%% ----------------------------------------------------------------------------------
	% Quick & Dirty fix for multiple FLAIR/T1w NIfTIs/scale slopes,
	% this should not give an error, but multiple ASL scale slopes should give an error
	% Check if FLAIR occurs
	isASL = true; % by default, assume these are ASL DICOMs
	if ~isempty(findstr(FileList{1},'FLAIR')) || ~isempty(findstr(FileList{1},'T1'))
		isASL = false;
	end
	
    TryDicominfo = true; % this is only set to false below upon succesful DcmtkRead
    
	if ~isempty(FileList)
		iMrFile = 0;
		
		% Check All dicom files throughout the sequence to validate that they are the same - because we are fast now!
		nFiles  = length(FileList);
		
		for iFile = 1:nFiles
            if TryDicominfo && iFile>1
                continue;
                % with dicominfo, reading is very slow, so we only read 1 dicom
            end
            
			ifname = FileList{iFile};
			filepath = fullfile(inp, ifname); % this is a file by definition, according to the xASL_adm_GetFileList command above
			
			%% ----------------------------------------------------------------------------------
			% Use DCMTK library to read the DICOM header to temp
			% ----------------------------------------------------------------------------------
            
            if bUseDCMTK
                try
                    temp = xASL_io_DcmtkRead(filepath, 0);
                    TryDicominfo = false;
				catch
                    warning(['xASL_adm_Dicom2Parms: xASL_io_DcmtkRead failed for ' filepath ', trying dicominfo']);
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
				catch
					try % 2) use default dictionary
						dicomdict('factory'); % temporarily set to default dictionary
						temp            = dicominfo(filepath); % retry reading dicom
					catch
						try % 3) change UseDictionaryVR setting with current dictionary
							dicomdict('set',DictionaryDCM);
							temp = dicominfo(filepath,'UseDictionaryVR', true);
						catch
							try % 4) change UseDictionaryVR setting with default dictionary
								dicomdict('factory');
								temp = dicominfo(filepath,'UseDictionaryVR', true);
							catch
								warning('xASL_adm_Dicom2Parms: dicominfo also did not work, check this!');
								dicomdict('set',DictionaryDCM); % reset dictionary
								continue;
							end
						end
					end
				end
				
				dicomdict('set',DictionaryDCM); % reset dictionary
			end
			
			%% -----------------------------------------------------------------------------
            % Take information from enhanced DICOM, if exists
            %% -----------------------------------------------------------------------------
			if isfield(temp, 'MediaStorageSOPClassUID')
				if strcmp(temp.MediaStorageSOPClassUID,'1.2.840.10008.5.1.4.1.1.4.1')==1 % Enhanced MR Image Storage
					bEnhancedMR = true;
					iMrFile = iMrFile+1;
				elseif strcmp(temp.MediaStorageSOPClassUID,'1.2.840.10008.5.1.4.1.1.4')==1 % MR Image Storage
					bEnhancedMR = false; % default
					iMrFile = iMrFile+1;
				else
% 					continue; % THIS SEEMS STRANGE >>>>>>>>>> RESULTS IN EMPTY PARMS, BY TRYING TO READ IMRFILE =0
                    iMrFile = iMrFile+1;
                    bEnhancedMR = false; % default
                end
			else
				bEnhancedMR = false; % default
				iMrFile = iMrFile+1;
			end
			
			% Deal with enhanced DICOM format imported through DICOMINFO and not DCMTK
            if bEnhancedMR && ~TryDicominfo
                %warning('Enhanced DICOM detected, but no dicominfo selected, skipping obtaining parameters');
            elseif bEnhancedMR && TryDicominfo
                % for simplicity, take the first value from the enhanced
                % sequences and store them in the temp-struct as if it is a
                % classic structure.

                if ~isfield(temp, 'PerFrameFunctionalGroupsSequence')
                    warning('Enhanced DICOM but PerFrameFunctionalGroupsSequence not found');
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
						warning('Manufacturer unknown for %s', filepath);
					end
				else
					warning('Manufacturer unknown for %s', filepath);
				end
				
				dcmfields = {'RepetitionTime', 'EchoTime', 'NumberOfAverages', 'RescaleSlope', ...
					'RescaleSlopeOriginal', 'MRScaleSlope', 'RescaleIntercept', 'AcquisitionTime', ...
					'AcquisitionMatrix'};
				
				switch bVendor
					case 'GE'
						dcmfields(end+1:end+2) = {'AssetRFactor', 'EffectiveEchoSpacing'}; % (0043,1083) (0043,102c)
					case 'Philips'
						dcmfields(end+1:end+2) = {'MRSeriesWaterFatShift', 'MRSeriesEPIFactor'}; % (2001,1022) (2001,1013)
					case 'Siemens'
						dcmfields(end+1) = {'BandwidthPerPixelPhaseEncode'}; % (0019,1028)
					otherwise
						% do nothing
				end
			end
			
			%% -----------------------------------------------------------------------------
			% Obtain the selected DICOM parameters from the header
			% Write the new parameter to the list (or put the default value
			% -----------------------------------------------------------------------------
			for iField=1:length(dcmfields)
				fieldname = dcmfields{iField};
				
				if  isfield(temp, fieldname)
					thevalue = temp.(fieldname);
					thevalue = thevalue(~isnan(thevalue));
				else
					thevalue = [];
				end
				
				if ~isempty(thevalue)
					if length(xASL_str2num(thevalue))>1
						tmpTheValue = nonzeros(thevalue);
						tmpTheValue = tmpTheValue(1);
					else
						tmpTheValue = thevalue;
					end
					
					t_parms.(fieldname)(iMrFile) = xASL_str2num(tmpTheValue);
				else
					if imPar.bVerbose; if iMrFile==1, fprintf('%s\n',['Parameter ' fieldname ' not found, default used']); end; end
					t_parms.(fieldname)(iMrFile) = DcmParDefaults.(fieldname);
				end
			end
		end
		
		%% If no files were found previously (just directories etc.) then the manufacturer won't be identified and 
		% dcmfields won't be assigned
		if exist('dcmfields','var')
			% -----------------------------------------------------------------------------
			% Dealing with empty or inconsistent field values
			% -----------------------------------------------------------------------------
			
            % Limit AcquisitionTime to one value
            if  isfield(t_parms,'AcquisitionTime')
                for iL=1:length(t_parms)
                    t_parms(iL).AcquisitionTime = xASL_str2num(t_parms(iL).AcquisitionTime);
                end

                if length(t_parms)>1
                    if imPar.bVerbose; fprintf('%s\n','Parameter AcquisitionTime has multiple values:'); end
%                     t_parms.AcquisitionTime
                    if imPar.bVerbose; fprintf('%s\n','using minumum value:'); end
%                     t_parms.AcquisitionTime = min(t_parms.AcquisitionTime)
                    tempAcquisitionTime = t_parms(1).AcquisitionTime;
                    t_parms = rmfield(t_parms,'AcquisitionTime');
                    t_parms(1).AcquisitionTime = tempAcquisitionTime;
                end
            end            
            
            t_parms(1).AcquisitionTime = xASL_adm_ConvertNr2Time(t_parms(1).AcquisitionTime);
            
			% Checks if the field values were the same for all dicoms and keep only one from the same value
			for iField=1:length(dcmfields)
				fieldname = dcmfields{iField};
				if  isfield(t_parms,fieldname)
					parms.(fieldname) = unique(t_parms.(fieldname));
					% There's one or more NaNs
					nNaN = sum(isnan(parms.(fieldname)));
					if nNaN > 0
						% Only NaNs
						if nNaN == length(parms.(fieldname))
							parms.(fieldname) = NaN;
						else
							parms.(fieldname) = parms.(fieldname)(~isnan(parms.(fieldname)));
						end
					end
				end
			end
		end
		
		% Check for valid RescaleSlope value
		if isfield(parms,'RescaleSlopeOriginal') && max(isnan(parms.RescaleSlopeOriginal))
			parms.RescaleSlopeOriginal = parms.RescaleSlope;
		end
		
		
		%% -----------------------------------------------------------------------------
		% Check and purge the parameters
		% -----------------------------------------------------------------------------
		% Check whether multiple scale slopes exist, this happens in 1
		% Philips software version, and should give an error (later we can
		% make this a warning, or try to deal with this)
		

		
        %% Calculate the TopUp parameters
        if isfield(parms,'AcquisitionMatrix')
            parms.AcquisitionMatrix = double(parms.AcquisitionMatrix(1));
        end

        switch bVendor
            case 'GE'
                if isfield(parms,'AssetRFactor') && isfield(parms,'EffectiveEchoSpacing') && isfield(parms,'AcquisitionMatrix')
                    parms.EffectiveEchoSpacing = double(parms.EffectiveEchoSpacing);
                    parms.AssetRFactor = double(parms.AssetRFactor);
                    parms.AcquisitionMatrix = double(parms.AcquisitionMatrix);
                    
                    parms.EffectiveEchoSpacing = (parms.EffectiveEchoSpacing * parms.AssetRFactor) / 10^6;
                    parms.TotalReadoutTime = (parms.AcquisitionMatrix-1) * parms.EffectiveEchoSpacing;
                    parms = rmfield(parms, 'AssetRFactor');
                end
            case 'Philips'
                if isfield(parms,'MRSeriesWaterFatShift') && isfield(parms,'MRSeriesEPIFactor') && isfield(parms,'AcquisitionMatrix')
                    parms.MRSeriesWaterFatShift = double(parms.MRSeriesWaterFatShift);
                    parms.MRSeriesEPIFactor = double(parms.MRSeriesEPIFactor);
                    parms.AcquisitionMatrix = double(parms.AcquisitionMatrix);
                    
                    EffectiveEchoSpacingPhilips = parms.MRSeriesWaterFatShift/(434.215 * (parms.MRSeriesEPIFactor+1));
                    parms.TotalReadoutTime = EffectiveEchoSpacingPhilips*(parms.MRSeriesEPIFactor-1);
                    parms.EffectiveEchoSpacing = parms.TotalReadoutTime/(parms.AcquisitionMatrix(1)-1);
                    parms = rmfield(parms, {'MRSeriesWaterFatShift' 'MRSeriesEPIFactor'});
                end
            case 'Siemens'
                if isfield(parms,'BandwidthPerPixelPhaseEncode') && isfield(parms,'InPlanePhaseEncodingDirection')
                    parms.BandwidthPerPixelPhaseEncode = double(parms.BandwidthPerPixelPhaseEncode);
                    
                    if strcmp(parms.InPlanePhaseEncodingDirection,'COL')
                       parms.ReconMatrixPE = double(parms.Rows);
                    elseif strcmp(parms.InPlanePhaseEncodingDirection,'ROW')
                       parms.ReconMatrixPE = double(parms.Columns);
                    else
                       error('Unknown InPlanePhaseEncodingDirection');
                    end
                    
                    parms.EffectiveEchoSpacing = 1/(parms.BandwidthPerPixelPhaseEncode*parms.ReconMatrixPE);
                    parms.TotalReadoutTime = (parms.ReconMatrixPE-1) * parms.EffectiveEchoSpacing;
                end
            otherwise
                % skip
        end
                    
		%% First remove non-finite values
        parms.EchoTime              = parms.EchoTime(isfinite(parms.EchoTime));
		parms.RepetitionTime        = parms.RepetitionTime(isfinite(parms.RepetitionTime));
		
		if  isfield(parms,'MRScaleSlope')
			parms.MRScaleSlope          = parms.MRScaleSlope(isfinite(parms.MRScaleSlope));
		end
		
		if  isfield(parms,'RescaleSlopeOriginal')
			parms.RescaleSlopeOriginal  = parms.RescaleSlopeOriginal(isfinite(parms.RescaleSlopeOriginal));
		end
		
		if  isfield(parms,'RescaleIntercept')
			parms.RescaleIntercept      = parms.RescaleIntercept(isfinite(parms.RescaleIntercept));
		end
		
		if  (length(parms.MRScaleSlope)>1 && parms.MRScaleSlope(2)~=1) || length(parms.RescaleSlopeOriginal)>1 || length(parms.RescaleIntercept)>1
			if  isASL % quickfix, see above
				warning('xASL_adm_Dicom2Parms: Multiple scale slopes exist for a single scan!');
                warning(['Could not perform dicom2nii conversion for ' parmsfile]);
                return;
			end
		end
		
		if ~isempty(parmsfile)
			% save as mat-file
			save(parmsfile, 'parms');
		end
% 	end
	
end

if ~isempty(parmsfile)
	X = load(parmsfile,'-mat');
% 	parms = X.parms % show me the values >>>>>>>>>> TOO VERBOSE
else
% 	parms % show me the values >>>>>>>>>> TOO VERBOSE
end

return

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
return
