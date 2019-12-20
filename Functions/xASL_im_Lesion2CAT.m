function LesionPathOut = xASL_im_Lesion2CAT(PathIn)
% xASL_im_Lesion2CAT For all lesion masks in the anatomical directory, load them, merge them and save them
% for the CAT segmentation

% PathIn - the name of the T1w file or another file in the same directory as the lesions

% LesionPathOut - the name and full path to the merged lesion file

%%%   If there are no lesions found, the images are untouched

% Parse the directory names
[x.SUBJECTDIR, x.P.STRUCT , ~]     = xASL_fileparts(PathIn);
x.P.STRUCT    = 'T1';
x.P.FLAIR     = 'FLAIR';

% Create the directory and specify the path to the merged lesion file
if ~exist(fullfile(x.SUBJECTDIR,'mri'),'dir'),mkdir(fullfile(x.SUBJECTDIR,'mri'));end
LesionPathOut = fullfile(x.SUBJECTDIR,'mri','LesionCAT.nii');

% Load the lesion names
Lesion_T1_list      = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.STRUCT '.*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);
Lesion_FLAIR_list   = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.FLAIR  '.*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);

% Initialize the rLesion list to have the same size as the original (or empty)
rLesion_FLAIR_list  = Lesion_FLAIR_list;

% First resample the lesion masks in FLAIR space to the T1 space
for iL=1:length(Lesion_FLAIR_list)
	[Fpath, Ffile, Fext]        = xASL_fileparts(Lesion_FLAIR_list{iL});
	rLesion_FLAIR_list{iL}      = fullfile(Fpath,['r' Ffile Fext]);
    xASL_spm_reslice( PathIn, Lesion_FLAIR_list{iL}, [], [], [], rLesion_FLAIR_list{iL}, 2);
end

%%% ---------------------------------------------------------------------
%%% Check if there are lesions

if  isempty(Lesion_T1_list) && isempty(rLesion_FLAIR_list)
	fprintf('%s\n','No Lesion*.nii specified, cost function masking skipped');
	LesionPathOut = '';
else

	%%% ---------------------------------------------------------------------
	fprintf('Cost function masking before running DARTEL:   ');
	% Initialize output image
	if  length(Lesion_T1_list)>length(rLesion_FLAIR_list)
		LesionImOut     = false(size(xASL_io_Nifti2Im(Lesion_T1_list{1})));
	else
		LesionImOut     = false(size(xASL_io_Nifti2Im(rLesion_FLAIR_list{1})));
	end

	% Collect all lesions into output image
	TotalLength     = length(Lesion_T1_list) + length(rLesion_FLAIR_list);
	for iL=1:length(Lesion_T1_list)
		xASL_TrackProgress(iL-1,TotalLength);
		LesionIM    = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(Lesion_T1_list{iL}));
		LesionImOut(LesionIM)           = 1;
	end
	for iL=1:length(rLesion_FLAIR_list)
		xASL_TrackProgress(iL+length(Lesion_T1_list)-1,TotalLength);
		LesionIM    = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(rLesion_FLAIR_list{iL}));
		LesionImOut(LesionIM)           = 1;
		xASL_delete(rLesion_FLAIR_list{iL}); % remove temporarily interpolated FLAIR lesion NIfTIs
	end
	LesionImOut = LesionImOut > 0.5;
	xASL_io_SaveNifti(PathIn,LesionPathOut,LesionImOut,16,0);xASL_TrackProgress(1,1);
end
fprintf('\n');

return
