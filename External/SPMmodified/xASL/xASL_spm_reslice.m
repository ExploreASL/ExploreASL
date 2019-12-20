function xASL_spm_reslice(refPath, srcPath, srcAffinePath, bInvAffine, bQuality, NewName, InterpolationPar)
% ExploreASL wrapper for SPM reslice function (coreg-write).
%
% FORMAT: xASL_spm_reslice(refPath, srcPath[, srcAffinePath, bInvAffine, bQuality, NewName, InterpolationPar])
%
% INPUT:
%   refPath             - path to reference space you want to resample the source image into (REQUIRED)
%                         entered either as 'path.nii' or as 'path.nii,1'
%                         can be a cell as well, but only the first cell is used
%   srcPath             - path to source image you want to resample (REQUIRED)
%                         entered either as 'path.nii' or as 'path.nii,1'
%                         can be a cell as well, but only the first cell is used
%   srcAffinePath       - path to the '*_sn.mat' that describes the Affine transformation matrix related
%                         to the source image and describing the src->ref transformation.
%   bInvAffine          - Do the inverse of the affine transform (OPTIONAL, default = FALSE)
%                         if the affine goes from src->ref, then FALSE. For ref-src, we need to use the inversion = TRUE
%   bQuality            - Quality setting affecting the speed of calculation (OPTIONAL, default = TRUE)
%   NewName             - Output name of the resampled NIfTI (OPTIONAL, default = uses srcPath but prefixes the input filename with 'r')
%   InterpolationPar    - interpolation method used for resampling (integer 0-7), higher is usually better quality but also slower
%                         options are 0) Nearest Neighbor 1) Trilinear 2-7) 2nd to 7th degree B-spline (OPTIONAL, default = 2)
%                         For bQuality==false, it resets it to maximum of 1.
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This wrapper runs SPM's reslice function (a.k.a. coregister: reslice) which resamples a source image into the space of a reference image,
% taking into account any orientation differences between the two images that are defined in the orientation matrix in the NIfTI header.
% When the source image contains multiple volumes, they are all resampled.
% The source image will get the same orientation matrix as the reference image,
% as it is now in the same space as the reference image. This can be useful when two images are to be compared voxel-by-voxel, e.g. when
% overlaying a CBF image over a structural image, or when wanting to use a mask from a different space. When after running this function, the
% reference and source images are not in alignment, they need to be registered first (i.e. xASL_spm_register).
% Resampling/reslicing needs an interpolation method to know from which voxels of the source image, a voxel in the new image will be computed.
% A simplistic explanation is that this determines the number of surrounding neighborhood voxels it uses to determine the new voxel value.
% The example syntax below would reslice/resample the CBF image to the T1w space (assuming the Affine was done to register CBF to T1w)
% It also works with the external .mat file of the source file that has the same name as the source file. It also can optionally take a _sn.mat containing the
% affine transformation information.
%
% EXAMPLE: xASL_spm_reslice('/TestDataSet/Sub-001/c1T1.nii', {'/TestDataSet/Sub-001/CBF.nii'}, {'xxx_sn.mat'},0,[], 0, '/TestDataSet/Sub-001/CBF_resampled.nii')
%          xASL_spm_reslice('/TestDataSet/Sub-001/c1T1.nii,3', '/TestDataSet/Sub-001/CBF.nii,2', '/TestDataSet/Sub-001/CBF_sn.mat')
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-08-01 HJM

    % ------------------------------------------------------------------------------------------
    % Admin, manage input
	% ------------------------------------------------------------------------------------------
	if (nargin < 2) || isempty(refPath) || isempty(srcPath)
		error('Needs two input arguments with input and output paths.');
	end

	% Input can be char or cell, in case of cell, take only the first input
	if iscell(srcPath)
		if length(srcPath) > 1
			warning('srcPath contains multiple cells, using the first only.');
		end
		srcPath = srcPath{1};
	end

	if iscell(refPath)
		if length(refPath) > 1
			warning('refPath contains multiple cells, using the first only.');
		end
		refPath = refPath{1};
	end

    [~, refFile, refExt] = xASL_fileparts(refPath);
    [~, srcFile, srcExt] = xASL_fileparts(srcPath);
    refFile = [refFile refExt];
    srcFile = [srcFile srcExt];
    refPath = xASL_spm_admin(refPath);
    srcPath = xASL_spm_admin(srcPath);

	if nargin < 3 || isempty(srcAffinePath)
		srcAffinePath = [];
    elseif iscell(srcAffinePath)
		srcAffinePath = srcAffinePath{1};
	end
	if nargin < 4 || isempty(bInvAffine)
		bInvAffine = false;
	end
	if nargin < 5 || isempty(bQuality)
		bQuality= true; % default quality is high
	end
	if nargin < 6
		NewName = [];
	end
	if iscell(NewName)
		NewName = NewName{1};
	end
	if nargin < 7 || isempty(InterpolationPar)
		InterpolationPar = 2;
	end

	% For low-quality setting, set it to maximum of 1
	if bQuality==false && InterpolationPar>1
		InterpolationPar = 1;
	end

	% ------------------------------------------------------------------------------------------
	% Prepare the files
	% ------------------------------------------------------------------------------------------
    % Always delete the temporary file with a prefix 'r'
	% Prepare this path and delete the result file if exists
	[Fpath, Ffile, Fext] = xASL_fileparts(srcPath{1}(1:end-2));
    rPath = fullfile(Fpath, ['r' Ffile Fext]);
    if xASL_exist(rPath)
        fprintf('%s\n', ['Warning, we delete ' rPath]);
        xASL_delete(rPath);
    end

	if isempty(NewName)
		NewName = rPath;
	end

	% And checks if the sn-mat exists and then loads it
	if ~isempty(srcAffinePath) && exist(srcAffinePath,'file')
		SnMat = load(srcAffinePath);
	else
		SnMat = [];
	end

	% ------------------------------------------------------------------------------------------
    % Manage the dimensionality of source NIfTIs and add the affine matrix
	% ------------------------------------------------------------------------------------------
    srcNii = xASL_io_ReadNifti(srcPath{1}(1:end-2));
    Sz = size(srcNii.dat);

	% Load the extra transformation matrices saved in a mat.file.
	if exist([srcPath{1}(1:end-5) 'mat'],'file')
		srcMat = load([srcPath{1}(1:end-5) 'mat']);
		srcMat = srcMat.mat;
		if ~isequal(srcNii.hdr.dim(5:8),[size(srcMat,3),size(srcMat,4),size(srcMat,5),size(srcMat,6)])
			warning('Image dimension and nD realignment MAT orientation dimensions dont match, not applying them');
			srcMat = [];
		end
	else
		srcMat = [];
    end

    if length(Sz)>3; fprintf('%s\n','Preparing reslice:...'); end % skip this otherwise, or create a bVerbose

    if ~min(srcNii.hdr.dim(2:length(Sz)+1)==Sz)
		warning(['Invalid NIfTI size: ' srcPath{1}]);
	elseif length(Sz)>3 % reslice multiple volumes separately
		DimO = srcNii.hdr.dim(5:8);
		iO = 1;

		% Pre-allocate memory
		FileName = cell(DimO(1)*DimO(2)*DimO(3)*DimO(4),1);
		rFileName = cell(DimO(1)*DimO(2)*DimO(3)*DimO(4),1);
		VolO = zeros(DimO(1)*DimO(2)*DimO(3)*DimO(4),4);

		% Create list of separate 3D volumes
		for iN1=1:DimO(1)
            if length(Sz)>3; xASL_TrackProgress(iN1,DimO(1)); end
			for iN2=1:DimO(2)
				for iN3=1:DimO(3)
					for iN4=1:DimO(4)
						VolO(iO,1:4) = [iN1 iN2 iN3 iN4];

						ImVolume = srcNii.dat(:,:,:,VolO(iO,1),VolO(iO,2),VolO(iO,3),VolO(iO,4));
						FileName{iO,1} = fullfile(Fpath, [Ffile '_' num2str(VolO(iO,1)) '_' num2str(VolO(iO,2)) '_' num2str(VolO(iO,3)) '_' num2str(VolO(iO,4)) Fext]);
						rFileName{iO,1} = fullfile(Fpath, ['r' Ffile '_' num2str(VolO(iO,1)) '_' num2str(VolO(iO,2)) '_' num2str(VolO(iO,3)) '_' num2str(VolO(iO,4)) Fext]);

						if isempty(SnMat)
							xASL_io_SaveNifti(srcPath{1}, FileName{iO,1}, ImVolume, [], false);
						else
							if bInvAffine
								mat = SnMat.VF.mat*SnMat.Affine;
							else
								mat = SnMat.VG.mat/SnMat.Affine;%SnMat.VG.mat*inv(SnMat.Affine);
							end
							if isempty(srcMat)
								locMat = srcNii.mat;
							else
								locMat = srcMat(:,:,iO);
							end

							% There is an external Mat file provided or potentially different internal MAT (different from the Affine trans definition) - take this into account
							if bInvAffine
								% Doing the inverse
								mat = (mat/SnMat.VG.mat)*locMat;%mat*inv(SnMat.VG.mat)*locMat;
							else
								mat = (mat/SnMat.VF.mat)*locMat;%mat*inv(SnMat.VF.mat)*locMat;
							end
							xASL_io_SaveNifti(srcPath{1}, FileName{iO,1}, ImVolume, [], false, mat);
						end
						iO              = iO+1;
					end
				end
			end
		end
	else
		% Simple 3D file - still check for Affine trans presence
		if isempty(SnMat)
			% Simple rigid - do nothing
			FileName = srcPath;
			rFileName = {fullfile(Fpath, ['r' Ffile Fext])};
		else
			if bInvAffine
				% Doing the inverse transformation
				mat = SnMat.VF.mat*SnMat.Affine;
			else
				mat = SnMat.VG.mat/SnMat.Affine;%SnMat.VG.mat*inv(SnMat.Affine);
			end
			% If there's not external mat, then take the internal MAT - because that can be different from the original MAT of the affine transformation
			if isempty(srcMat)
				srcMat = srcNii.mat;
			end
			% There is an external Mat file provided or potentially different internal MAT - take this into account
			if bInvAffine
				% Doing the inverse
				mat = (mat/SnMat.VG.mat)*srcMat;%mat*inv(SnMat.VG.mat)*srcMat;
			else
				mat = (mat/SnMat.VF.mat)*srcMat;%mat = mat*inv(SnMat.VF.mat)*srcMat;
			end
			xASL_io_SaveNifti(fullfile(Fpath,[Ffile Fext]),fullfile(Fpath,['r' Ffile Fext]),srcNii.dat(:,:,:,1),[],0,mat);
			FileName = {fullfile(Fpath, ['r' Ffile Fext ',1'])};
			rFileName = {fullfile(Fpath, ['rr' Ffile Fext])};
		end
    end

    if length(Sz)>3; xASL_TrackProgress(DimO(1),DimO(1)); end

	% ------------------------------------------------------------------------------------------
	% Run the reslicing
	% ------------------------------------------------------------------------------------------
	%  Default SPM settings
    matlabbatch{1}.spm.spatial.coreg.write.ref                  = refPath;
    matlabbatch{1}.spm.spatial.coreg.write.source               = FileName;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp      = InterpolationPar;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap        = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask        = 0;
	matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix      = 'r';

	%  Screenprint & run
	fprintf('\n%s\n','------------------------------------------------------------------------------------------');
	InterpMeth  = {'Nearest Neighbor' 'Trilinear' '2nd B-spline' '3rd degree Bspline' '4th degree Bspline' '5th degree Bspline' '6th degree Bspline' '7th degree Bspline'};
    fprintf('%s\n',['Reslicing ' srcFile ' to ' refFile ' using ' InterpMeth{InterpolationPar+1} ' interpolation']);

%     if length(Sz)>3; warning('off','all'); end
%     warning('off','QFORM0 representation has been rounded.'); -> has no
%     identifier
    spm_jobman('run',matlabbatch);
%     if length(Sz)>3; warning('on','all'); end

	% ------------------------------------------------------------------------------------------
    % Put volumes back together
	% ------------------------------------------------------------------------------------------
	if length(Sz)>3; fprintf('%s\n','Putting volumes back together:...'); end

    if length(Sz)>3
		for iO=1:length(rFileName)
            if length(Sz)>3; xASL_TrackProgress(iO, length(rFileName)); end
			NewIm(:,:,:,VolO(iO,1),VolO(iO,2),VolO(iO,3),VolO(iO,4)) = xASL_io_Nifti2Im(rFileName{iO,1});
        end
        if length(Sz)>3; xASL_TrackProgress(1, 1); end
		TempName = fullfile(Fpath, ['r' Ffile Fext]);
		xASL_io_SaveNifti(rFileName{1,1}, TempName, NewIm, [], 0);

		for iO=1:length(rFileName)
			xASL_delete(FileName{iO,1});
            [Fpath, Ffile] = xASL_fileparts(FileName{iO,1});
            xASL_delete(fullfile(Fpath, [Ffile '.mat']));

            xASL_delete(rFileName{iO,1});
            [Fpath, Ffile] = xASL_fileparts(rFileName{iO,1});
            xASL_delete(fullfile(Fpath, [Ffile '.mat']));
		end
		rFileName{1,1} = TempName;
	end

	% ------------------------------------------------------------------------------------------
	% Rename temporary SPM file(s) into output file(s)
	% ------------------------------------------------------------------------------------------
	% We used the affine transform - then delete the renamed file
	if (length(Sz)<=3) && (~isempty(SnMat))
		xASL_delete(FileName{1,1});
	end

	% Check if it needs renaming to the new name
	if ~strcmp(NewName,rFileName{1,1})
		xASL_Move(rFileName{1,1},NewName,1);
	end

end
