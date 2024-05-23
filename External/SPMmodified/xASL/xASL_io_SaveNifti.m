function xASL_io_SaveNifti(pathOrigNifti, pathNewNifti, imNew, nBits, bGZip, changeMat, bCopyOrigJson, JsonFields, bLegacy2BIDS, bOverwrite, changeMat0)
% Save a file to a Nifti format, while taking the parameters from another file
%
% FORMAT: xASL_io_SaveNifti(pathOrigNifti, pathNewNifti, imNew[, nBits, bGZip, changeMat, bCopyOrigJson, JsonFields, bLegacy2BIDS, bOverwrite, changeMat0])
%
% INPUT:
%   pathOrigNifti  Path to the original Nifti file to take parameters from (REQUIRED)
%   pathNewNifti   Name of the file to save the results to (REQUIRED)
%   imNew          Image matrix to save (REQUIRED)
%                  The dimension must correspond to the dimension of pathOrigNifti
%   nBits           Number of bits to save the result in - 8,16,32
%                  (OPTIONAL, by DEFAULT it checks which bits
%                  representation is enough.
%                  For bit conversion, 32 is best precision for sensitive data, 16 is still
%                  OK and saves some space, 8 is most economic but should only be used for
%                  masks, since it cannot contain a large data range
%   bGZip          Gzip the result (OPTIONAL, DEFAULT 1)
%   changeMat      New orientation matrix 4x4 (OPTIONAL, DEFAULT same as previous)
%   bCopyOrigJson  Copy a JSON sidecar file (if existed for the original file) as sidecar of the new NIfTI (OPTIONAL, DEFAULT = false)
%   JsonFields     Save additional JsonFields to a new JSON sidecar (or merged with the original Json content)
%                  These fields overwrite fields with the same fieldnames in the original Json
%                  (OPTIONAL, DEFAULT = struct() empty)
%   bLegacy2BIDS   Convert saved JSON from Legacy to BIDS. This only applies to the JsonFields passed on the input. The JSON read from 
%                  the origNifti is in BIDS already and is never converted. This means that the saved JSON is always in BIDS, this only 
%                  solves the problem if the input is in Legacy or already in BIDS (OPTIONAL, DEFAULT = true)
%   bOverwrite     Vector with 3 Booleans for overwriting the following pre-existing destination files:
%                  1) NIfTI file, 2) JSON file, 3) mat-orientation motion file (for 4D NIfTIs) (OPTIONAL, DEFAULT = [1 1 1]; 
%                  1 -> [1 1 1] and 0 -> [0 0 0] and creates a warning; [1 0] - error
%   changeMat0     New orientation matrix mat0 4x4. Note that mat0 should always stay as it was originally from the scanner, 
%                  but we can exceptionally change it, e.g. if mat0 was incorrectly converted from the DICOM header (OPTIONAL, DEFAULT same as previous)
%                  
%
% JSON saving options are:
% 1. bCopyOrigJson = false & JsonFields = empty  -> we don't save a new JSON sidecar
% 2. bCopyOrigJson = true  & JsonFields = empty  -> we copy the original JSON sidecar only
% 3. bCopyOrigJson = false & JsonFields = filled -> we create a new JSON sidecar with new fields only
% 4. bCopyOrigJson = true  & JsonFields = filled -> we merge the original JSON sidecar with the new fields
%                                                   the new fields have priority
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It loads the pathOrigNifti, takes all the parameters from it, and creates a new Nifti file with
%              these parameters, but new image matrix from imNew. It saves the result in pathNewNifti.
%              It runs the following steps:
% 
%              1. Unzipping and manage name input file
%              2. Determine the bit precision
%              3. Create new NIfTI
%              4. Remove redundant .mat orientation files
%              5. Manage scale slopes
%              6. Save new NIfTI
%
% EXAMPLE: xASL_io_SaveNifti('c:\User\path\old.nii', 'c:\User\path\new.nii', im)
%          xASL_io_SaveNifti('c:\User\path\old.nii', 'c:\User\path\new.nii', im, [], 0)
%          xASL_io_SaveNifti('c:\User\path\old.nii', 'c:\User\path\new.nii', im, 32, 0)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2024 ExploreASL

%% ----------------------------------------------------------------------------
%% Admin
if nargin < 3
	error('Needs at least three input parameters.');
end

if isempty(imNew) || min(size(imNew))==0 % QA
    error(['Empty image, could not create: ' pathNewNifti]);
end

if nargin < 5 || isempty(bGZip)
	bGZip = 1; % zip by default
end

if nargin < 6 || isempty(changeMat)
	changeMat = [];
else
	if ~isequal(size(changeMat),[4 4])
		error('changeMat has to be 4x4');
	end
end

if nargin < 7 || isempty(bCopyOrigJson)
    bCopyOrigJson = false;
end

if nargin < 8 || isempty(JsonFields)
    JsonFields = struct;
elseif ~isstruct(JsonFields)
    error('JsonFields input variable should be a struct');
end

if nargin < 9 || isempty(bLegacy2BIDS)
	bLegacy2BIDS = true;
end

if nargin < 10 || isempty(bOverwrite)
	bOverwrite = [true true true];
elseif sum(isinf(bOverwrite)) > 0 || sum(isnan(bOverwrite)) > 0
	error('bOverwrite should not contain INF or NaN');
elseif length(bOverwrite) == 1
	bOverwrite(1,[2,3]) = bOverwrite(1);
	warning('bOverwrite should have a length of 3');
elseif length(bOverwrite) == 2
	error('bOverwrite should have a length of 3');
end

if nargin < 11 || isempty(changeMat0)
	changeMat0 = [];
else
	if ~isequal(size(changeMat0),[4 4])
		error('changeMat0 has to be 4x4');
	end
end

% JSON saving, see detailed explanation above.
% If bCopyOrigJson is true OR JsonFields is not empty, we set bSaveJson to true
% So we only save a JSON in either or both of these cases
if bCopyOrigJson || ~isempty(fields(JsonFields))
    bSaveJson = true;
else
    bSaveJson = false;
end


% If the absolute path is missing and filename is given only, then add the current path to the absolute path
% Do this both for the new filename and the original filename
[PathNew, NameNew, ExtNew] = xASL_fileparts(pathNewNifti);
if isempty(PathNew)
	% If a file-name only, then add the full current path to avoid ambiguity
	pathNewNifti = fullfile(pwd(), [NameNew ExtNew]);
end

[PathOri, NameOri, ExtOri] = xASL_fileparts(pathOrigNifti);
if isempty(PathOri)
	% If a file-name only, then add the full current path to avoid ambiguity
	pathOrigNifti = fullfile(pwd(), [NameOri ExtOri]);
end

% Define the JSON sidecar paths
pathNewJson = fullfile(PathNew, [NameNew '.json']);
pathOrigJson = fullfile(PathOri, [NameOri '.json']);

% Create temporary name for new NIFTI, since if pathOrigNifti & pathNewNifti
% are the same, this will work better
tempName = fullfile(PathNew, [NameNew '_temp.nii']);
tempMat = fullfile(PathNew, [NameNew '_temp.mat']);
newMat = fullfile(PathNew, [NameNew '.mat']);

% Remove temp files in case they exist from a previous crash
if xASL_exist(tempName)
    warning(['Temporary file already existed, perhaps from a previous crash? Removing: ' tempName]);
    xASL_delete(tempName);
    xASL_delete(tempMat);
end

%% ====================================================================================
%% 1. Unzipping and manage name input file
% First unzip original Nifti if needed
xASL_io_ReadNifti(pathOrigNifti);

% Then change name if needed
pathOrigNifti = xASL_adm_ZipFileNameHandling(pathOrigNifti);

% This will make sure that the newly created nifti file has the correct name (otherwise it may complain that it wants to read .nii.gz which doesn't exist anymore
newNifti = xASL_io_ReadNifti(pathOrigNifti);
newNifti.dat.fname = tempName;


%% ====================================================================================
%% 2. Determine the bit precision
bInteger8 = min(min(min(min(min(min(min(uint8(imNew)==imNew)))))));
bInteger16 = min(min(min(min(min(min(min(int16(imNew)==imNew)))))));

if nargin < 4 || isempty(nBits)
    % Automatically adapt according to data input bitdepth
    if bInteger8
        nBits = 8;
    elseif bInteger16
        nBits = 16;
    else
        % this image had more precision, save as single precision
        nBits = 32;
    end
end

switch nBits
    case 32
        % Convert to single floating point, larger data & slower
        % processing but virtually no rounding errors
        newNifti.dat.dtype = 'FLOAT32-LE';
        imNew = single(imNew);

    case 16
        % Convert to int16, to save space & speed up processing.
        % Leads to minimal/negligible rounding errors, if original
        % image wasn't 16 bit
        newNifti.dat.dtype = 'INT16-LE';

        if bInteger16
            % this image is integer16 already, doesn't need scale slope
            imNew = int16(imNew);
        else
            bImFinite = imNew(isfinite(imNew));
            ScaleSlope16 = max([-min(bImFinite(:)) max(bImFinite(:))]) / 2^15; % define new range (this was signed, stays signed, hence 2^16-1)
            imNew = imNew ./ ScaleSlope16;
            imNew = int16(imNew);
        end

    case 8
        % Convert to UINT8-LE, e.g. for masks.
        % Leads to minimal/negligible rounding errors, if original
        % image wasn't 16 bit
        newNifti.dat.dtype = 'UINT8-LE';

        if bInteger8
            % this image is UINT8 already, doesn't need scale slope
            imNew = uint8(imNew);
        else
            InterceptN   = floor(min(imNew(:)));
            bImFinite    = imNew(isfinite(imNew));
            ScaleSlope16 = (max(bImFinite(:)) - min(bImFinite(:)))/255; % define new range (this is converted to unsigned)
            if  ScaleSlope16==0 % theoretical case
                ScaleSlope16 = 1;
            end
            imNew = imNew-InterceptN;
            imNew = round(imNew ./ ScaleSlope16);
            % Because a nifti viewer will do
            % (RawValue+Intercept)*ScaleSlope
        end

    otherwise
        error('Unknown bit-choice.');
end


%% ====================================================================================
%% 3. Create new NIfTI
if ~isempty(changeMat)
	newNifti.mat  = changeMat;
end

if ~isempty(changeMat0)
	warning('Note that mat0 (the original NIfTI orientation) is changed upon request. This can jeopardize re-running ExploreASL when done in the processing module!');
	newNifti.mat0 = changeMat0;
end

newNifti.dat.scl_slope = 1;
newNifti.dat.scl_inter = 0;
newNifti.dat.dim       = [size(imNew,1) size(imNew,2) size(imNew,3) size(imNew,4) size(imNew,5) size(imNew,6) size(imNew,7)];

% Create new NIFTI
xASL_adm_CreateDir(fileparts(newNifti.dat.fname));
create(newNifti);

newNifti.dat(:,:,:,:,:) = imNew;


%% ====================================================================================
%% 4. Remove redundant .mat orientation files
% Start with removing any pre-existing destination .mat files
% Even if we don't create a new one, we don't want the wrong .mat motion orientation file with a new NIfTI

if exist(newMat, 'file')
	if bOverwrite(3)
        xASL_delete(newMat);
	end % Otherwise don't do anything
end

if size(imNew,4)==1
    xASL_delete(tempMat);
end

if exist(tempMat, 'file')
	motionMat = load(tempMat);
	% Remove .mat if dimensions do not fit
	if size(newNifti.dat,4) ~= size(motionMat.mat,3)
		xASL_delete(tempMat);
	end
end


%% ====================================================================================
%% 5. Manage scale slopes
if exist('ScaleSlope16', 'var')
    newNifti = xASL_io_ReadNifti(tempName);
    newNifti.dat.scl_slope = ScaleSlope16;
    create(newNifti);
end
if exist('InterceptN', 'var')
    newNifti = xASL_io_ReadNifti(tempName);
    newNifti.dat.scl_inter = InterceptN;
    create(newNifti);
end


%% ====================================================================================
%% 6. Save NIfTI (& -mat sidecar)
% Remove existing NIfTI
if xASL_exist(pathNewNifti)
	if bOverwrite(1)
        xASL_delete(pathNewNifti);
	end % Otherwise don't do anything
end

if bOverwrite(1) || ~xASL_exist(pathNewNifti, 'file')
    % if we want to overwrite an existing NIfTI-file (or none existed)
    xASL_Move(tempName, pathNewNifti, 1, 0);
end

if exist(tempMat, 'file')
    % if a mat-file exists that we need
    if bOverwrite(3) || ~exist(newMat, 'file')
        % if we want to overwrite an existing mat-file (or none existed)
        xASL_Move(tempMat, newMat, 1, 0);
    end
end


%% ====================================================================================
%% 7. Create JSON-sidecar
% 0. Remove JSON sidecar if it already exists
% Even if we don't save a new one, then we still don't want the wrong sidecar to a new NIfTI
% But don't delete the NewJson in case the new and reference file are equivalent
if ~strcmp(pathNewNifti, pathOrigNifti)
	if bOverwrite(2)
        xASL_delete(pathNewJson);
	end % Otherwise don't do anything
end

if bSaveJson
    % a. First, we load the reference sidecar JSON file, if it exists
	if bCopyOrigJson && xASL_exist(pathOrigJson, 'file')
        json = xASL_io_ReadJson(pathOrigJson);
    elseif bCopyOrigJson
        warning([pathOrigJson ' did not exist']);
        json = struct;
    else
        json = struct;
	end

	% b. Second, we convert provided JSON fields to BIDS if requested
	if bLegacy2BIDS && ~isempty(JsonFields) && ~isempty(fieldnames(JsonFields))
		JsonFields = xASL_bids_parms2BIDS(JsonFields, [], 1);
	end

    % c. Third, we add any provided JSON fields
    fieldNames = fields(JsonFields);
    for iField = 1:length(fieldNames)
        json.(fieldNames{iField}) = JsonFields.(fieldNames{iField});
    end

    % d. Fourth, we save the JSON
    if bOverwrite(2) || ~xASL_exist(pathNewJson, 'file')
        xASL_io_WriteJson(pathNewJson, json, 1);
    end
end


end