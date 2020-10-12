function xASL_io_SaveNifti(pathOrigNifti, pathNewNifti, imNew, nBits, bGZip, changeMat)
% Save a file to a Nifti format, while taking the parameters from another file
%
% FORMAT: xASL_io_SaveNifti(pathOrigNifti, pathNewNifti, imNew[, nBits, bGZip, newMat])
%
% INPUT:
%   pathOrigNifti  Path to the original Nifti file to take parameters from (REQUIRED)
%   pathNewNifti   Name of the file to save the results to (REQUIRED)
%   imNew          Image matrix to save (REQUIRED)
%                  The dimension must correspond to the dimension of pathOrigNifti
%   nBits           Number of bits to save the result in - 8,16,32
%                  (OPTIONAL, by DEFAULT it checks if 16 bits representation is enough or 32 are needed not to
%                  loose the precition of imNew
%                  For bit conversion, 32 is best precision for sensitive data, 16 is still
%                  OK and saves some space, 8 is most economic but should only be used for
%                  masks, since it cannot contain a large data range
%   bGZip          Gzip the result (OPTIONAL, DEFAULT 1)
%   changeMat      New orientation matrix 4x4 (OPTIONAL, DEFAULT same as previous)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It loads the pathOrigNifti, takes all the parameters from it, and creates a new Nifti file with
%              these parameters, but new image matrix from imNew. It saves the result in pathNewNifti.
%
% EXAMPLE: xASL_io_SaveNifti('c:\User\path\old.nii', 'c:\User\path\new.nii', im)
%          xASL_io_SaveNifti('c:\User\path\old.nii', 'c:\User\path\new.nii', im, [], 0)
%          xASL_io_SaveNifti('c:\User\path\old.nii', 'c:\User\path\new.nii', im, 32, 0)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

% Admin
if nargin < 3
	error('Needs at least three input parameters.');
end

if min(size(imNew))==0 % QA
    error(['Empty image: ' pathOrigNifti]);
end

if nargin < 5 || isempty(bGZip)
	bGZip   = 1; % zip by default
end

if nargin < 6 || isempty(changeMat)
	changeMat = [];
else
	if ~isequal(size(changeMat),[4 4])
		error('changeMat has to be 4x4');
	end
end

% If the absolute path is missing and filename is given only, then add the current path to the absolute path
% Do this both for the new filename and the original filename
[pathstr, name0, ext0] = fileparts(pathNewNifti);
if isempty(pathstr)
	% If a file-name only, then add the full current path to avoid ambiguity
	pathNewNifti = fullfile(xASL_adm_UnixPath(pwd()),[name0 ext0]);
end

[pathstr, name0, ext0] = fileparts(pathOrigNifti);
if isempty(pathstr)
	% If a file-name only, then add the full current path to avoid ambiguity
	pathOrigNifti = fullfile(xASL_adm_UnixPath(pwd()),[name0 ext0]);
end

% Create temporary name for new NIFTI, since if pathOrigNifti & pathNewNifti
% are the same, this will work better
tempName = [pathNewNifti(1:end-4) '_temp.nii'];
tempMat = [pathNewNifti(1:end-4) '_temp.mat'];
newMat = [pathNewNifti(1:end-4) '.mat'];

% Remove temp files in case they exist from a previous crash
if xASL_exist(tempName)
    warning(['Temporary file already existed, removing: ' tempName]);
    xASL_delete(tempName);
end
if xASL_exist(tempMat)
    warning(['Temporary file already existed, removing: ' tempMat]);
    xASL_delete(tempMat);
end

% First unzip original Nifti if needed
xASL_io_ReadNifti(pathOrigNifti);

% Then change name if needed
pathOrigNifti = xASL_adm_ZipFileNameHandling(pathOrigNifti);

% this will make sure that the newly created nifti file has the correct name (otherwise it may complain that it wants to read .nii.gz which doesn't exist anymore
newNifti = xASL_io_ReadNifti(pathOrigNifti);
newNifti.dat.fname = tempName;

% Determines the bit precision
if nargin < 4 || isempty(nBits)
    bImInt16 = int16(imNew)==imNew;

    if  min(bImInt16(:))==1 && newNifti.dat.scl_slope==1
        % this image is integer16 already, doesn't need FLOAT32
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
        bImInt16 = int16(imNew)==imNew;

        if  min(bImInt16(:))==1
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
        bImInt8 = uint8(imNew)==imNew;

        if  min(bImInt8(:))==1
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
        error('xASL_io_SaveNifti: Unknown bit-choice.');
end

if ~isempty(changeMat)
	newNifti.mat = changeMat;
end

newNifti.dat.scl_slope = 1;
newNifti.dat.scl_inter = 0;
newNifti.dat.dim       = [size(imNew,1) size(imNew,2) size(imNew,3) size(imNew,4) size(imNew,5) size(imNew,6) size(imNew,7)];

% Create new NIFTI
xASL_adm_CreateDir(fileparts(newNifti.dat.fname));
create(newNifti);

newNifti.dat(:,:,:,:,:) = imNew;

if  exist('ScaleSlope16','var')
    newNifti = xASL_io_ReadNifti(tempName);
    newNifti.dat.scl_slope = ScaleSlope16;
    create(newNifti);
end
if  exist('InterceptN','var')
    newNifti = xASL_io_ReadNifti(tempName);
    newNifti.dat.scl_inter = InterceptN;
    create(newNifti);
end

xASL_Move(tempName,pathNewNifti,1,0);

if  exist(tempMat,'file')
    xASL_Move(tempMat,newMat,1,0);
end


%% Remove redundant .mat orientation files
if size(imNew,4)==1
    xASL_delete(newMat);
end

if exist(newMat,'file')
	tmpMat = load(newMat);
	% Remove .mat if dimensions do not fit
	if size(newNifti.dat,4) ~= size(tmpMat.mat,3)
		xASL_delete(newMat);
	end
end


%% Always avoid having two of the same files, of which one copy is zipped
% E.g. in a rerun
if  exist([pathNewNifti '.gz'],'file')
    delete([pathNewNifti '.gz']);
end

if  bGZip
    xASL_adm_GzipNifti(pathNewNifti);
end

if  strcmp(pathNewNifti(end-3:end),'.nii') && exist(pathNewNifti,'file') && exist([pathNewNifti '.gz'],'file')
    delete([pathNewNifti '.gz']);
end

end
