function imOut = xASL_io_Nifti2Im(niftiIn, ImageSize)
% Load image matrix from NIfTI given as a file path or preloaded
%
% FORMAT: imOut = xASL_io_Nifti2Im(niftiIn [, ImageSize])
%
% INPUT:
%   niftiIn   - can be one of the following (REQUIRED):
%                 1) path to NIfTIfile to load image matrix from
%                 2) NIfTI object to load image matrix from
%                 3) image matrix (to simply pass through)
%   ImageSize - if NIfTI doesnt exist or is corrupt,
%               will create a dummy image matrix with this size
%               (OPTIONAL, DEFAULT=none)
%
% OUTPUT:
%   imOut     - image matrix with single precision
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads a NIfTI image matrix with flexible input
%              (as explained under INPUT: niftiIn). It does the following.
%              1) Try to load a NIfTI
%              2) If NIfTI successfully loaded, try to load the NIfTI image
%              3) If the above didnt work, try to create a dummy image
%              4) Convert to single precision data format
%              5) Also able to load NIfTI as .nii.mat format
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: imOut = xASL_io_Nifti2Im('/analysis/Sub-001/ASL_1/CBF.nii');
% __________________________________
% Copyright 2015-2019 ExploreASL

% Admin
if nargin < 1
	error('xASL_io_Nifti2Im: Needs at least one input.');
end
if nargin<2
	ImageSize = [];
end
if iscell(niftiIn)
    niftiIn = niftiIn{1};
end

niiMat = false; % default

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 1) Load the NIfTI
% -----------------------------------------------------------------------------------------------------------------------------------------------------
if numel(niftiIn) == 1 % Assume this is a NIfTI
    nii = niftiIn;
elseif isnumeric(niftiIn) || islogical(niftiIn) && numel(niftiIn)>1000 % assume this is an image already
    imOut = niftiIn;
else % assume this is a path, we try to open
    try
        [Fpath, Ffile] = xASL_fileparts(niftiIn);
        niiPathMat = fullfile(Fpath, [Ffile '.nii.mat']);
        if exist(niiPathMat,'file')
            niiMat = true; % load mat below
        else
            nii = xASL_io_ReadNifti(niftiIn);
        end
    catch ME
        warning(['Could not load ' niftiIn]);
        fprintf('%s\n',['Message: ' ME.message]);
    end
end

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 2) Load the image
% -----------------------------------------------------------------------------------------------------------------------------------------------------
if ~exist('imOut','var') && exist('nii','var') && ~niiMat
    % if succesfully loaded the NIfTI
    try
        imOut = nii.dat(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
    catch ME
        warning(['Corrupt ' niftiIn]);
        fprintf('%s\n',['Message: ' ME.message]);
    end
end

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 3) Try creating dummy output if above crashed
% -----------------------------------------------------------------------------------------------------------------------------------------------------
if ~exist('imOut','var') && ~niiMat
	try % Try if we can find the supposed NIfTI size
		ImageSize = size(nii.dat);
	catch
		warning(['Cannot find image size ' niftiIn]);
	end

	% Create an empty dummy image
	if ~isempty(ImageSize)
		imOut = zeros(ImageSize);
		imOut(:) = NaN;
	end
end

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 4) Convert to single precision data format
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPM opens to double precision data format by default, which takes up
% unnecessary large memory
if exist('imOut','var') && ~isa(imOut,'single') && ~niiMat
    imOut = single(imOut);
elseif ~exist('imOut','var') && ~niiMat
    imOut = [];
    fprintf('Something went wrong opening NifTI\n');
end

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 5) Load mat file if it exists
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% In this case we don't convert to any other format, but just load the .nii.mat

if niiMat
    try
        niiMat = load(niiPathMat,'-mat');
        FieldsMat = fields(niiMat);
        if isempty(FieldsMat)
            error(['No image/variable found in ' niiPathMat]);
        elseif length(FieldsMat)>1
            error(['Too many images/variables found in ' niiPathMat]);
        else
            imOut = niiMat.(FieldsMat{1});
        end
    catch ME
        warning(['Corrupt ' niiPathMat]);
        fprintf('%s\n',['Message: ' ME.message]);
    end
end        


% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 6) Check scaling
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% In the import stage, the scaling can go wrong. In this case, we throw a
% warning & try to correct this

if exist('imOut', 'var') && exist('nii', 'var')
	if nii.dat.scl_slope>16
		% this if-clause speeds up in all other cases that don't have this potential Philips import issue
		MaxIm = max(imOut(isfinite(imOut)));
		if MaxIm>1e9
			if exist('Fpath', 'var') && exist('Ffile', 'var')
				if ~isempty(regexp(Ffile,'^.*(T1|FLAIR).*$'))
					warning('%s\n%s','Found a structural image with unusually large NIfTI scale slope and image data, resetting the maximum to 4096 upon reading:',niftiIn);
					imOut = imOut.*4096./MaxIm;
				else
					warning('%s\n%s','Found unusually large NIfTI scale slope and image data, check if all processing and quantification went correctly:',niftiIn);
				end
			end
		end
	end
end


end
