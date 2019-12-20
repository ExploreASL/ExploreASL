function xASL_adm_ManageMoCoMat(PathIn)
% Manage the .mat orientation sidecar that SPM produces from realign
%
% FORMAT: xASL_adm_ManageMoCoMat(PathIn)
%
% INPUT:
%   PathIn  - path to NIfTI file (REQUIRED)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function manages the orientation matrices that SPM puts
%              in an external .mat sidecar file if there are more than 1
%              volumes. The first volume should be equal to the orientation
%              header of the NIfTI, if not, we assume that the NIfTI
%              header is correct.
%              This function performs several checks & corrects if
%              necessary, combined with throwing a warning:
%              A) the nVolumes in .mat & .nii image should be equal, if not, delete sidecar
%              B) .mat should have more than one volume, if not delete sidecar
%              C) If there are illegal numbers in the diagonal of the .mat
%                 orientation matrices (here only checked for zeros or non
%                 finite values) then the .mat is removed
%              D) If this is true for the first volume only, the .mat is
%                 retained but the first volume orientation is overwritten
%                 with a zero matrix
%
% EXAMPLE: xASL_adm_ManageMoCoMat('/analysis/Sub-001/ASL_1/ASL4D.nii');
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

% ------------------------------------------------------------------------------------------------------------
% Admin, input arguments & file loading
% ------------------------------------------------------------------------------------------------------------
if nargin < 1 || isempty(PathIn)
	error('Needs at least one input.');
end

[Fpath, Ffile, Fext] = xASL_fileparts(PathIn);

if ~(strcmp(Fext,'.nii') || strcmp(Fext,'.nii.gz'))
    warning('Input was no NIfTI file, skipping...');
    return;
end

PathMat = fullfile(Fpath, [Ffile '.mat']);
if ~exist(PathMat,'file')
    return;
    % only go through this function when there are both nii & mat files
end

NiftiObject = nifti(PathIn); % here xASL_io_ReadNifti would create an infinite loop
mat = load(PathMat,'-mat');

% ------------------------------------------------------------------------------------------------------------
% Manage the .mat orientation sidecar
% ------------------------------------------------------------------------------------------------------------
if size(mat.mat,3)~=size(NiftiObject.dat,4)
	% ------------------------------------------------------------------------------------------------------------
    % A) The dimensions are not equal
	% ------------------------------------------------------------------------------------------------------------
    warning('.MAT sidecar had different size than the image. Deleting the sidecar');
    fprintf('%s\n', ['Deleting ' PathMat]);
    clear mat
    xASL_delete(PathMat);
	return;
end
	
if size(mat.mat,3)==1 % We test only the sidecar dimensions because the Nifti has the same dimensions as tested previously
	% ------------------------------------------------------------------------------------------------------------
    % B) There's only a single volume in the sidecar
	% ------------------------------------------------------------------------------------------------------------
    warning('.MAT sidecar contains only 1 volume. Using the NifTI header instead and deleting the sidecar');
    clear mat
    xASL_delete(PathMat);
	return;
end

% ------------------------------------------------------------------------------------------------------------
% C) check any illegal numbers in the diagonal
% ------------------------------------------------------------------------------------------------------------
MaskEye = logical(zeros(4,4));
MaskEye(1:3,1:3) = eye(3,3); % create mask for diagonal
MaskEyeFull = repmat(MaskEye,[1 1 size(mat.mat,3)]);
MaskEyeFull(:,:,1) = 0; % ignore the first volume

DiagonalMat = squeeze(mat.mat(MaskEyeFull)); % get diagonal data
if sum(~isfinite(DiagonalMat) | DiagonalMat==0)>0
	warning('.MAT sidecar volumes contained illegal diagnonal. Deleting the sidecar');
	clear mat
	xASL_delete(PathMat);
	return;
end

% ------------------------------------------------------------------------------------------------------------
% D) Check the same as in C for the first volume
% ------------------------------------------------------------------------------------------------------------
MaskEyeFull = repmat(MaskEye,[1 1 size(mat.mat,3)]);
MaskEyeFull(:,:,2:end) = 0;
DiagonalMat = squeeze(mat.mat(MaskEyeFull)); % get diagonal data
	
if sum(~isfinite(DiagonalMat) | DiagonalMat==0)>1
	% If the whole matrix is zero, then it's ok, when only the diagonal, then report an error
	if sum(sum(abs(mat.mat(:,:,1)))) > 0
		warning('Corrupt first volume of .MAT sidecar, replacing it with a zero matrix');
		mat = mat.mat;
		mat(:,:,1) = 0;
		save(PathMat,'mat');
	end
end

end
