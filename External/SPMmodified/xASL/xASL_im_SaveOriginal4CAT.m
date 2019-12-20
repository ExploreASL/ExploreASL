function xASL_im_SaveOriginal4CAT( Ycls, PathIn )
%xASL_im_SaveOriginal4CAT Save the segmentation before lesion masking

%%%   If there are no lesions found, the images are untouched

[x.SUBJECTDIR, x.P.STRUCT , ~]     = xASL_fileparts(PathIn);
x.P.STRUCT    = 'T1';
x.P.FLAIR     = 'FLAIR';

Lesion_T1_list      = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.STRUCT '.*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);
Lesion_FLAIR_list   = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.FLAIR  '.*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);

if  isempty(Lesion_T1_list) && isempty(Lesion_FLAIR_list)
	fprintf('%s\n','No Lesion*.nii specified, not saving the segmentation before masking');
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
	volCATvoxelSize   = volCAT.dim .* volOriginalvoxelSize./size(Ycls{1});

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
end
fprintf('\n');

return
