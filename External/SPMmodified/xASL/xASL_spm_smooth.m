function xASL_spm_smooth(pathIn, fwhmSmooth, pathNew)
% ExploreASL SPM wrapper for SPM smooth function to smooth a volume with a given Gaussian kernel
%
% FORMAT: xASL_spm_smooth(pathIn, fwhmSmooth[, pathNew])
%
% INPUT:
%   pathIn     - path to NifTI image you want to smooth. Must be a single file, but NIfTI file can contain multiple volumes (REQUIRED)
%                can be a .NII or NII.GZ
%   fwhmSmooth - Full-Width-Half-Maximum (FWHM, mm) of 1x3 Gaussian kernel you want to smooth by  (REQUIRED)
%   pathNew    - the output file name, by default SPM smooth prefixes a 's' to the input path (OPTIONAL)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This SPM wrapper runs SPM's smooth function, which spatially smooths the input image with a Gaussian kernel.
% In the case of multiple volumes (i.e. a 4D NIfTI), each 3D volume is spatially smoothed separately.
% Note that smoothnesses combine with Pythagoras' rule (i.e. sum quadratically)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_smooth('/MyStudy/Subject1/T1.nii.gz', [8 8 8]);
%          xASL_spm_smooth('/MyStudy/Subject1/T1.nii.gz', [8 8 8],'newName.nii');
% __________________________________
% Copyright 2015-2019 ExploreASL

% ------------------------------------------------------------------------------------------
% Admin

if nargin < 2
	error('xASL_spm_smooth: Needs at least two input arguments.');
end

if length(fwhmSmooth) ~= 3
	error('xASL_spm_smooth: The smoothing kernel must have size 3.');
end

if nargin<3 || isempty(pathNew)
    pathNew = [];
end

% Load the file
[Fpath, Ffile, Fext] = xASL_fileparts(pathIn);
nii = xASL_io_ReadNifti(pathIn);
Dims4 = size(nii.dat,4);
Dims5 = size(nii.dat,5);

% ------------------------------------------------------------------------------------------
% Save multiple volumes as individual NIfTIs
% ------------------------------------------------------------------------------------------
for i4=1:Dims4
	for i5=1:Dims5
		% Save to 3D volumes and keep all the names
		if i4*i5 == 1
			Path2Smooth{1} = pathIn;
			Dim4(1) = 1;
			Dim5(1) = 1;
			NN = 2;
		else
			Path2Smooth{NN,1} = fullfile(Fpath,[Ffile '_' num2str(NN) '.nii']);
			xASL_io_SaveNifti(pathIn,Path2Smooth{NN},nii.dat(:,:,:,i4,i5),[], false);

			Dim4(NN) = i4;
			Dim5(NN) = i5;
			NN  = NN+1;
		end
	end
end

% ------------------------------------------------------------------------------------------
% Run the Matlab job for all NIfTI volumes
% ------------------------------------------------------------------------------------------
matlabbatch{1}.spm.spatial.smooth.fwhm   = double(fwhmSmooth);
matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
matlabbatch{1}.spm.spatial.smooth.im     = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
matlabbatch{1}.spm.spatial.smooth.data   = Path2Smooth;

spm_jobman('run',matlabbatch);

fprintf('%s\n',['Smoothed ' Ffile Fext ' with [' num2str(fwhmSmooth(1),3) ' ' num2str(fwhmSmooth(2),3) ' ' num2str(fwhmSmooth(3),3) '] mm FWHM Gaussian filter']);

% ------------------------------------------------------------------------------------------
% Reconstruct image from multiple volumes
% ------------------------------------------------------------------------------------------
if length(Path2Smooth)>1
	% Allocated the smoothed image
	imSmoothed = zeros(size(nii.dat));
    for iVolume=1:length(Path2Smooth)
		% Load the smoothed image
        [Fpath2, Ffile2, Fext2] = xASL_fileparts(Path2Smooth{iVolume});
        SmoothedPath = fullfile(Fpath2,['s' Ffile2 Fext2]);
        imSmoothed(:,:,:,Dim4(iVolume),Dim5(iVolume)) = xASL_io_Nifti2Im(SmoothedPath);

		% Delete the (two) temporary images
        xASL_delete(SmoothedPath);
        if iVolume>1
            xASL_delete(Path2Smooth{iVolume});
        end
	end

	% Save the completed image
    TempName = fullfile(Fpath, ['s' Ffile Fext]);
    xASL_io_SaveNifti(pathIn,TempName,imSmoothed,[], false);
end

% ------------------------------------------------------------------------------------------
% Rename temporary SPM file into output file
% ------------------------------------------------------------------------------------------
if ~isempty(pathNew) % otherwise don't rename, keep prefix "s" for smoothed

    TempName = fullfile(Fpath, ['s' Ffile Fext]);

    if ~strcmp(pathNew,TempName)
        xASL_Move(TempName,pathNew,1);
    end
end

end
