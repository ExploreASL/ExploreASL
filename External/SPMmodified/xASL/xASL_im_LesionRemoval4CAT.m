function [ Ycls, LesionImOut] = xASL_im_LesionRemoval4CAT( Ycls, PathIn )
%xASL_im_LesionRemoval4CAT For all lesion masks in the anatomical directory, remove
%them from the current segmentations

%%%   If there are no lesions found, the images are untouched

[x.SUBJECTDIR, x.P.STRUCT , ~]     = xASL_fileparts(PathIn);
x.P.STRUCT    = 'T1';
x.P.FLAIR     = 'FLAIR';

Lesion_T1_list      = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.STRUCT '.*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);
Lesion_FLAIR_list   = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.FLAIR  '.*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);

% Initialize the rLesion list to have the same size as the original (or empty)
rLesion_FLAIR_list  = Lesion_FLAIR_list;

% First resample the lesion masks in FLAIR space to the T1 space
for iL=1:length(Lesion_FLAIR_list)
	[Fpath, Ffile, Fext]        = xASL_fileparts(Lesion_FLAIR_list{iL});
	rLesion_FLAIR_list{iL}      = fullfile(Fpath,['r' Ffile Fext]);
	xASL_spm_reslice( PathIn, Lesion_FLAIR_list{iL}, [], [], x.Quality,'Lesion_FLAIR 2 T1 space', rLesion_FLAIR_list{iL}, 2 );
end

%%% ---------------------------------------------------------------------
%%% Check if there are lesions
LesionImOut                             = false(size(Ycls{1})); % default

if  isempty(Lesion_T1_list) && isempty(Lesion_FLAIR_list)
	fprintf('%s\n','No Lesion*.nii specified, cost function masking skipped');
else

	%%% ---------------------------------------------------------------------
	%%% Save probability maps before DARTEL & before lesion masking
	cPathOri{1}   = fullfile( x.SUBJECTDIR, ['c1' x.P.STRUCT '_ORI.nii']);
	cPathOri{2}   = fullfile( x.SUBJECTDIR, ['c2' x.P.STRUCT '_ORI.nii']);

	fprintf('%s','Saving segmentations before DARTEL (and/or cost function masking if Lesion_*.nii exist):   ')

	% First checking if they are the same size, or were resampled
	T1Path      = fullfile( x.SUBJECTDIR, [x.P.STRUCT '.nii']);
	T1im        = xASL_io_Nifti2Im(T1Path);

	%%% Second, we calculate the volumes (dimension, transformation matrix in the CAT12 space and original space)

	% The parameters of the CAT12 space are taken from the original space and modified
	volCAT = spm_vol(T1Path); % Loads the parameters of the original space
	volOriginalvoxelSize    = sqrt(sum(volCAT.mat(1:3,1:3).^2)); % calculates the voxelsize
	volOriginalvoxelSize    = round(volOriginalvoxelSize,2); % round to two digits after the decimal point

	% Take the size of the CAT12 space and adapt the resolution
	volCATvoxelSize   = round(volCAT.dim .* volOriginalvoxelSize./size(Ycls{1}));

	% Replicate what CAT12 does to recalculate the CAT12-space matrix
	volCAT        = rmfield(volCAT,'private');
	volCAT.descrip = '';
	imat      = spm_imatrix(volCAT.mat);
	volCAT.dim    = size(Ycls{1});
	imat(7:9) = volCATvoxelSize .* sign(imat(7:9));
	volCAT.mat    = spm_matrix(imat);

	% The parameters of the original T1 space are unchanged
 	volOriginal = spm_vol(T1Path);
	volOriginal = rmfield(volOriginal,'private');

	% Transfer CAT12 to original = volCAT -> volOriginal
	for ii=1:2
		xASL_TrackProgress(ii-1,2);

		% Save the image
		xASL_delete(cPathOri{ii});

		if  ~min(size(Ycls{1})==size(T1im)) % if sizes are not identical
			% Transform the CAT12 image to the original space
			imHigh2Save = xASL_im_ResampleIM(Ycls{ii},volCAT.mat,volOriginal.mat,volOriginal.dim,'linear');
			% And save the result
			xASL_io_SaveNifti(PathIn,cPathOri{ii},imHigh2Save,8);
		else
			xASL_io_SaveNifti(PathIn,cPathOri{ii},Ycls{ii},8);
		end

		xASL_TrackProgress(ii,2);
	end

	fprintf('\n');

	%%% ---------------------------------------------------------------------
	%%% Remove lesion from probability maps
	fprintf('Cost function masking before running DARTEL:   ');
	% Initialize output image
	if  length(Lesion_T1_list)>length(Lesion_FLAIR_list)
		LesionImOut     = false(size( xASL_io_Nifti2Im(Lesion_T1_list{1}) ));
	else
		LesionImOut     = false(size( xASL_io_Nifti2Im(Lesion_FLAIR_list{1}) ));
	end

	% Collect all lesions into output image
	TotalLength     = length(Lesion_T1_list) + length(rLesion_FLAIR_list);
	for iL=1:length(Lesion_T1_list)
		xASL_TrackProgress(iL,TotalLength);
		LesionIM    = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(Lesion_T1_list{iL}));
		LesionImOut(LesionIM)           = 1;
	end
	for iL=1:length(rLesion_FLAIR_list)
		xASL_TrackProgress(iL+length(Lesion_T1_list),TotalLength);
		LesionIM    = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(rLesion_FLAIR_list{iL}));
		LesionImOut(LesionIM)           = 1;
		delete(rLesion_FLAIR_list{iL}); % remove temporarily interpolated FLAIR lesion NIfTIs
	end

	% Resample output image to resampled CAT12 space, if necessary
	if  ~min(size(Ycls{1})==size(LesionImOut)) % if sizes are not identical
		% Resample from the original space to CAT12 space
		imHigh2Save = xASL_im_ResampleIM(LesionImOut,volOriginal.mat,volCAT.mat,volCAT.dim,'linear');
		% And save
		imHigh2Save = imHigh2Save > 0.5;
	else
		% Or save without transforming
		imHigh2Save = LesionImOut > 0.5;
	end

	% Remove lesion from CAT12 segmentations (cost function masking)
	for ii=1:length(Ycls)
		Ycls{ii}(logical(imHigh2Save))   = NaN;
	end
end
fprintf('\n');

return
