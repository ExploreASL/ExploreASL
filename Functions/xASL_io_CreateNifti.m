function xASL_io_CreateNifti(pathNewNifti, imNew, resMat, nBits, bGZip)
% Creates a new Nifti image.
%
% FORMAT: xASL_io_CreateNifti(pathNewNifti, imNew, resMat, nBits, bGZip)
%
% INPUT:
%   pathNewNifti   Name of the file to save the results to (REQUIRED)
%   imNew          Image matrix to save (REQUIRED)
%                  The dimension must correspond to the dimension of pathOrigNifti 
%   resMat         Either a 1x3 vector of the voxel size or a 4x4 matrix describing the transformation matrix (OPTIONAL, DEFAULT [1 1 1])
%   nBits           Number of bits to save the result in - 8,16,32 
%                  (OPTIONAL, by DEFAULT it checks if 16 bits representation is enough or 32 are needed not to 
%                  loose the precition of imNew
%                  For bit conversion, 32 is best precision for sensitive data, 16 is still
%                  OK and saves some space, 8 is most economic but should only be used for
%                  masks, since it cannot contain a large data range
%   bGZip          Gzip the result (OPTIONAL, DEFAULT 1)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It loads the pathOrigNifti, takes all the parameters from it, and creates a new Nifti file with
%              these parameters, but new image matrix from imNew. It saves the result in pathNewNifti.
%
% EXAMPLE: xASL_io_CreateNifti('c:\User\path\new.nii', im)
%          xASL_io_CreateNifti('c:\User\path\new.nii', im, [3 3 7])
%          xASL_io_CreateNifti('c:\User\path\new.nii', im, [], 32, 0)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2016-00-00 HM

if nargin < 2
	error('xASL_io_CreateNifti: Minimum of two arguments is needed.');
end

if nargin < 3 || isempty(resMat)
	resMat = [1 1 1];
end

if nargin < 4 || isempty(nBits)
	nBits = 32;
end
	
if nargin < 5 || isempty(bGZip) 
    bGZip   = 1; % zip by default
end

newNifti = nifti();
newNifti.dat = file_array(pathNewNifti);%,dim,dtype,offset,scl_slope,scl_inter,permission)

if min(size(imNew))==0 
    error('xASL_io_CreateNifti: Empty image');
end
   
switch nBits
    case 32
        % Convert to single floating point, larger data & slower
        % processing but virtually no rounding errors 
        newNifti.dat.dtype = 'FLOAT32-LE';
        imNew              = single(imNew);

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
            bImFinite    = imNew(isfinite(imNew));
            ScaleSlope16 = max([-min(bImFinite(:)) max(bImFinite(:))]) / 2^15; % define new range (this was signed, stays signed)
            imNew        = imNew ./ ScaleSlope16;
            imNew        = int16(imNew);
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
        error('xASL_io_CreateNifti: Unknown bit-choice.');
end

newNifti.dat.scl_slope = 1;
newNifti.dat.scl_inter = 0;
if length(resMat) == 3
	newNifti.mat = [resMat(1) 0 0 0;0 resMat(2) 0 0; 0 0 resMat(3) 0; 0 0 0 1];
else
	if isequal(size(resMat),[4 4])
		newNifti.mat = resMat;
	else
		error('xASL_io_CreateNifti: Incorrect size of resMat.Either voxel size or transformation matrix has to be given');
	end
end
newNifti.mat0    = newNifti.mat;
newNifti.dat.dim = [size(imNew,1) size(imNew,2) size(imNew,3) size(imNew,4) size(imNew,5) size(imNew,6) size(imNew,7)];


% Create new NIFTI
xASL_adm_CreateDir(fileparts(newNifti.dat.fname));
create(newNifti);

newNifti.dat(:,:,:,:,:) = imNew;

if  exist('ScaleSlope16','var')
    newNifti.dat.scl_slope = ScaleSlope16;
    create(newNifti);
end
if  exist('InterceptN','var')
    newNifti.dat.scl_inter   = InterceptN;
    create(newNifti);
end

% Always avoid having two of the same files, of which one copy is zipped
% E.g. in a rerun
if  exist([pathNewNifti '.gz'],'file')
    delete([pathNewNifti '.gz']);
end

if  bGZip
    pathNewNifti = xASL_adm_ZipFileNameHandling(pathNewNifti);
    gzip(pathNewNifti);
    delete(pathNewNifti);
end

if  strcmp(pathNewNifti(end-3:end),'.nii') && exist(pathNewNifti,'file') && exist([pathNewNifti '.gz'],'file')
    delete([pathNewNifti '.gz']);
end

end