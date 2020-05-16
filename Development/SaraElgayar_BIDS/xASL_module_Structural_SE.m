function [result x] = xASL_module_Structural(job, x)
% $Rev:: 72                    $:  Revision of last commit
% $Author:: hjmutsaerts        $:  Author of last commit
% $Date:: 2016-07-23 10:37:55 #$:  Date of last commit

%% ExploreASL T1 module

% 1    Rigid-body register T1 -> MNI (run for multiple time points in
%      case of longitudinal registration)
% 2    Longitudinal registration (if requested)
% 3    Segmentation T1
% 4    Get tissue volumes
% 5    Reslice to common space
% 6    Visual QC

%% Administration

x        = xASL_init_GenericMutexModules( x, 'T1' ); % starts mutex locking process to ensure that everything will run only once
x        = xASL_init_FileSystem(x);


% Find current subject index
    for subi=1:length(x.SUBJECTS)
        C = strsplit(x.SUBJECTS{subi},'-');
       x.SUBJECTS{subi} = C{2}; %C{1} = 'sub' C{2} = '001'
    end
       iSubj   = find(strcmp(x.SUBJECTS,x.P.SubjectID));

% 0.9 change working directory to make sure that unspecified output will go there...
oldFolder = cd(x.SUBJECTDIR);



%% Start with throwing a warning if the current subject does not have a structural image
if ~xASL_exist(x.P.Path_T1,'file') %&& ~xASL_exist(x.D.Path_T1_ORI,'file')
    warning('This subject didnt have a structural image, skipping...');
    result = true;
    return;
end

% deprecated - will pass the quality option directly to toolboxes through matlabbatch
%QName   = fullfile(x.SUBJECTDIR,'LowQ.status');

% Other
if  (~x.mutex.HasState('010_coreg_T12MNI') || ~x.mutex.HasState('020_coreg_FLAIR2T1w') || ~x.mutex.HasState('030_resample_FLAIR2T1w') || ~x.mutex.HasState('040_segment_FLAIR') || ~x.mutex.HasState('050_LesionFilling')) && ~x.mutex.HasState('999_ready')
    %xASL_adm_UnzipNativeFiles4Processing(x.D.SubDer,x); % if the script searches for *.nii
end

% Which WMH segmentation algorithm to use
if  xASL_exist(x.P.Path_WMH_PreSEGM,'file')
    x.WMHsegmAlg  = 'LPA';
	% force fast WMH segmentation, because already WMH_SEGM.nii existing
    % which we assume to be high quality

elseif ~isfield(x,'WMHsegmAlg')
        x.WMHsegmAlg  = 'LPA'; % default WMH segmentation algorithm, seems more robust & faster than LGA
end

%% -----------------------------------------------------------------------------
%% Lesion NIfTIs admin

Lesion_T1_list      = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.STRUCT '_\d*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);
Lesion_FLAIR_list   = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.FLAIR '_\d*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);

ROI_T1_list         = xASL_adm_GetFileList(x.SUBJECTDIR,['^ROI_' x.P.STRUCT '_\d*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);
ROI_FLAIR_list      = xASL_adm_GetFileList(x.SUBJECTDIR,['^ROI_' x.P.FLAIR '_\d*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);

if  ~isempty(Lesion_T1_list) || ~isempty(Lesion_FLAIR_list)
    CreateDir(x.D.LesionCheckDir);
end
if  ~isempty(ROI_T1_list) || ~isempty(ROI_FLAIR_list)
    CreateDir(x.D.ROICheckDir);
end


if      strcmp(x.WMHsegmAlg,'LPA')
        rWMHPath               = fullfile( x.SUBJECTDIR, ['ples_lpa_mr' x.P.FLAIR '.nii']);
elseif  strcmp(x.WMHsegmAlg,'LGA')
        rWMHPath                = fullfile( x.SUBJECTDIR, ['ples_lga_0.3_rmr' x.P.FLAIR '.nii']);% 0.3 an initial threshold.
else
        error('Unknown WMH segmentation option');
end
[fPath, File, ~]    = fileparts(rWMHPath);
T1_filledName       = fullfile(fPath,['T1_filled_' File(6:end) '.nii']);

if ~isfield(x,'T1wBasedFLAIRBiasfieldCorrection')
    x.T1wBasedFLAIRBiasfieldCorrection = false;
end

%% -----------------------------------------------------------------------------
%% Check whether we require certain files to be present, or otherwise to skip this subject
% By default we don't check this (because a population might not contain a
% FLAIR or M0, or this differs between subjects). This is useful when the
% data is incomplete, but one wants to start image processing nevertheless

if ~isfield(x,'SkipIfNoFlair');  x.SkipIfNoFlair    = 0;  end
if ~isfield(x,'SkipIfNoASL');    x.SkipIfNoASL      = 0;  end
if ~isfield(x,'SkipIfNoM0');     x.SkipIfNoM0       = 0;  end

bContinue = true;
if x.SkipIfNoFlair && ~xASL_exist(x.D.Path_FLAIR,'file')
   bContinue  = false;
end
% for iSess=1:x.nSessions % need to define various ASL4D_session files still
%     if ~xASL_exist(x.P.Path_ASL4D,'file') && x.SkipIfNoASL
%        bContinue  = false;
%     elseif ~xASL_exist(x.P.Path_M0,'file') && x.SkipIfNoM0
%        bContinue  = false;
%     end
%
% end


if ~bContinue
    fprintf('%s\n','Skipping this subject, because no FLAIR, ASL or M0 present whereas this was requested');
    return;
end



%% -----------------------------------------------------------------------------
%% Check if we need to reset the original T1w/FLAIR,
%% This is the case if any of the .status before 007 is missing

if ~x.mutex.HasState('010_coreg_T12MNI') || ~x.mutex.HasState('070_segment_T1w')

    DelList     = {x.P.Path_c1T1 x.P.Path_c2T1 x.P.Path_yT1 x.P.Path_PseudoCBF};
    for iD=1:length(DelList)
        if  xASL_exist(DelList{iD},'file'); xASL_delete(DelList{iD}); end
    end
end

%% Check for existing WMH_SEGM.nii, in which case we will perform low quality FLAIR segmentation,
%  & use the WMH_SEGM.nii instead
if  xASL_exist(x.P.Path_WMH_PreSEGM,'file') && ~x.mutex.HasState('999_ready')
    warning('Existing WMH_SEGM.nii detected, so we will perform a low quality FLAIR segmentation');
    fprintf('%s\n','& use the existing WMH_SEGM.nii');
    fprintf('%s\n','However, if this WMH_SEGM.nii was created in a previous run, consider removing it first!');
end


if  x.mutex.HasState('999_ready')
    bO  = 0; % no Output, as everything has been done already
else
    bO  = 1; % yes Output, not completely done, so we want to know what has or has not been done
end


% Deprecated - pass quality parameter directly
%%% Create low quality status file, for low quality processing hacking inside SPM toolboxes
%if  x.Quality==0
%    fclose(fopen(QName, 'w'));
%else
%    xASL_delete(QName);
%end



%% -----------------------------------------------------------------------------
%% 1    Register T1w -> MNI

if ~x.mutex.HasState('010_coreg_T12MNI') % tracks progress through lock/*.status files, & locks current run

    %% Restore the orientation matrix of all images, in case we perform a
    % re-run: but only when we don't have lesion maps

    if  length(Lesion_T1_list)<1  && ~xASL_exist(x.P.Path_WMH_PreSEGM,'file') && length(Lesion_FLAIR_list)<1
        % only if there are no lesion maps that should remain registered to
        % T1, restore orientation of the T1w image
        xASL_im_RestoreOrientation(x.P.Path_T1);
        xASL_im_RestoreOrientation(x.P.Path_FLAIR);

        % then also restore other functional images that should remain in
        % registration with T1w
        %for iSess = 1:x.nSessions
        %    xASL_im_RestoreOrientation( x.P.Path_ASL4D ); % this needs to be redone with multiple sessions
        %    xASL_im_RestoreOrientation( x.P.Path_M0 );
        %end
    end



    %% Register T1w linearly to the center of MNI space
    %% Align with ACPC
    % & the same transformations are applied to all other related scans
    % (ASL4D, M0, FLAIR, etc.)

    xASL_im_ClipExtremes( x.P.Path_T1,0.999,0); % first we clip high vascular intensities

    % saving ClipExtremes output just for testing
    % NewIM = xASL_im_ClipExtremes( x.D.Path_T1,0.999,0); % first we clip high vascular intensities
    %newIMpath = fullfile(x.D.SubDer, BIDS_FileName(x.P.SubjectID,x.P.SessionID, '','','ClipExtremes','T1', 'nii'));
    %xASL_save_nii( x.D.Path_T1, newIMpath, NewIM, [], 0 );
    %MoveFileIf(newIMpath, x.D.Path_T1);

    refPath = fullfile(x.D.MapsDir,'T1_template.nii'); % = SPM8 T1 template, to make sure it is always the same
    % This T1 reference image, is a blurred one. In SPM12 there are new reference images, but less blurry.
    % This one empirically works fine
    OtherList{1,1}                  = x.P.Path_FLAIR;
    OtherList{end+1,1}              = x.P.Path_WMH_PreSEGM;

    % Add lesion masks to the registration list
    for iS=1:length(Lesion_T1_list)
        OtherList{end+1,1}      = Lesion_T1_list{iS};
    end
    for iS=1:length(Lesion_FLAIR_list)
        OtherList{end+1,1}      = Lesion_FLAIR_list{iS};
    end

    % Add ASL images to the registration list
    % for iSess = 1:x.nSessions
    %      OtherList{end+1,1}      = x.P.Path_ASL4D;
    %    OtherList{end+1,1}      = x.P.Path_M0;
    %end

    % First, start with center of mass detection & realign with this
    xASL_im_CenterOfMass(x.P.Path_T1, OtherList);
    xASL_spm_coreg( refPath, x.P.Path_T1, OtherList, x, [4 2]);

    x.mutex.AddState('010_coreg_T12MNI');
   % xASL_adm_CompareDataSets([], [], x); % unit testing
    x.mutex.DelState('020_coreg_FLAIR2T1w');
elseif  bO; fprintf('%s\n','010_coreg_T12MNI has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 2    Register FLAIR -> T1w (if FLAIR.nii exists)
if ~x.mutex.HasState('020_coreg_FLAIR2T1w') % tracks progress through lock/ *.status files, & locks current run
    if  xASL_exist(x.P.Path_FLAIR,'file')

        fprintf('\n');
        fprintf('%s\n','---------------------------------');
        fprintf('%s\n','FLAIR.nii detected, processing...');

        OtherList           = {x.P.Path_WMH_PreSEGM };

        for iS=1:length(Lesion_FLAIR_list)
            OtherList{end+1,1}      = Lesion_FLAIR_list{iS};
        end

        xASL_im_ClipExtremes( x.P.Path_FLAIR,0.999999,0); % first we clip very high intensities
        %newIMpath = fullfile(x.D.SubDer, BIDS_FileName(x.P.SubjectID,x.P.SessionID, '','','ClipExtremes','FLAIR', 'nii'));
        %xASL_save_nii( x.D.Path_FLAIR, newIMpath, NewIM, [], 0 );
        %MoveFileIf(newIMpath, x.D.Path_FLAIR);

        xASL_im_CenterOfMass(x.P.Path_FLAIR, OtherList);
        xASL_spm_coreg( x.P.Path_T1, x.P.Path_FLAIR, OtherList, x);

        x.mutex.AddState('020_coreg_FLAIR2T1w');
       % xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState('025_BiasfieldCorrection');
    elseif  bO; fprintf('%s\n','020_coreg_FLAIR2T1w: No FLAIR data found, skipping...');
    end

elseif  bO; fprintf('%s\n','020_coreg_FLAIR2T1w has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 1.5  T1w-based FLAIR biasfield correction

% Here we acquire a biasfield from the T1 and apply it to the FLAIR, as with large lesions, the biasfield
% cannot be reliably estimated on the FLAIR. By default = disabled, enable this in case of large lesions.


if  xASL_exist(x.P.Path_WMH_PreSEGM,'file') || ~xASL_exist(x.P.Path_FLAIR,'file')
    % We skip now
    if  bO
        fprintf('%s\n','WMH segmentation already existing, or FLAIR.nii missing');
        fprintf('%s\n','Hence we skip the FLAIR biasfield correction');
    end
else

    if ~x.mutex.HasState('025_BiasfieldCorrection') && xASL_exist(x.P.Path_FLAIR,'file')

        if ~x.T1wBasedFLAIRBiasfieldCorrection
            fprintf('%s\n','025_BiasfieldCorrection not requested, skipping...');
        else
            xASL_spm_BiasfieldCorrection( x.P.Path_T1, x );

            %% apply the same biasfield to the FLAIR image
            if  xASL_exist(x.P.Path_FLAIR,'file')
                % Create biasfield image
                StructIm_BF0    = xASL_io_Nifti2Im(x.P.Path_T1);
                StructIm_BF1    = xASL_io_Nifti2Im(x.P.Path_mT1);
                BiasfieldIM     = StructIm_BF1./StructIm_BF0;

                % Extrapolating biasfield
                fprintf('%s\n','Extrapolating biasfield:   ');
                BiasfieldIM = xASL_im_ndnanfilter(BiasfieldIM,'gauss',[8 8 8]);
                xASL_save_nii(x.P.Path_mT1,x.P.Path_BiasField_T1,BiasfieldIM,[],0);

                % First resample to FLAIR space
                xASL_spm_reslice( x.P.Path_FLAIR, x.P.Path_BiasField_T1, [], [], x.Quality, x.P.Path_BiasField_FLAIR);

                % Open FLAIR & biasfield & multiply
                IM1     = xASL_io_Nifti2Im(x.P.Path_FLAIR);
                BF1     = xASL_io_Nifti2Im(x.P.Path_BiasField_FLAIR);
                IM2     = IM1.*BF1;
                %%SAve
                xASL_save_nii(x.P.Path_FLAIR,x.P.Path_mFLAIR,IM2,16,0);

                %MoveFileIf(x.P.Path_FLAIR, x.P.Path_FLAIR_ORI);
                %MoveFileIf(x.P.Path_mFLAIR, x.P.Path_FLAIR,1);
                xASL_delete(x.P.Path_BiasField_FLAIR);
            end

            MoveFileIf(x.P.Path_T1, x.P.Path_T1_ORI);
            MoveFileIf(x.P.Path_mT1,x.P.Path_T1,1);
            xASL_delete(x.P.Path_BiasField_T1);
            %%Save File

            xASL_io_ConvertNii16bit(x.P.Path_T1);
        end

        x.mutex.AddState('025_BiasfieldCorrection'); % tracks progress through lock/ *.status files, & locks current run
       % xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState('030_resample_FLAIR2T1w');
    elseif  bO; fprintf('%s\n','025_BiasfieldCorrection has already been performed, skipping...');
    end
end



%% -----------------------------------------------------------------------------
%% 3    Resample native FLAIR to native T1w space for lesion segmentation (if FLAIR.nii exists)

if ~x.mutex.HasState('030_resample_FLAIR2T1w')  % tracks progress through lock/ *.status files, & locks current run
    if  xASL_exist(x.P.Path_FLAIR,'file')

        xASL_io_ConvertNii16bit(x.P.Path_FLAIR);

        % save output in BIDS filename formate
        [Fpath, Ffile, Fext] = xASL_fileparts(x.P.Path_FLAIR);
        %get modality suffix
        C = strsplit(Ffile,'_');
        suffix = C{length(C)};
        NewName         = fullfile(Fpath, BIDS_FileName(x.P.SubjectID, x.P.SessionID,'','','','resample',suffix,'nii'));
        xASL_spm_reslice( x.P.Path_T1, x.P.Path_FLAIR, [], [], x.Quality, NewName);

        if  xASL_exist(x.P.Path_WMH_PreSEGM,'file')
            % save output in BIDS filename formate
            NewName         = fullfile(Fpath, BIDS_FileName(x.P.SubjectID, x.P.SessionID,'','WMH','','ResampleSegmented','','nii'));
            xASL_spm_reslice( x.P.Path_T1, x.P.Path_WMH_PreSEGM, [], [], x.Quality, NewName);
        end

        x.mutex.AddState('030_resample_FLAIR2T1w');
        xASL_adm_CompareDataSets([], [], x); % unit testing
    elseif  bO; fprintf('%s\n','030_resample_FLAIR2T1w: No FLAIR data found, skipping...');
    end

elseif  bO; fprintf('%s\n','030_resample_FLAIR2T1w has already been performed, skipping...');
end




%% -----------------------------------------------------------------------------
%% 4    Lesion Segmentation Toolbox: segment FLAIR.nii if available
if ~x.mutex.HasState('040_segment_FLAIR')  % tracks progress through lock/ *.status files, & locks current run
    if  xASL_exist(x.P.Path_rFLAIR,'file')

        if  strcmp(x.WMHsegmAlg,'LPA')
            fprintf('%s\n','WMH segmentation performed using LST LPA')
            % Lesion Prediction Algorithm (LPA)
             rWMHPath = xASL_wrp_LST_lpa(x);

        elseif strcmp(x.WMHsegmAlg,'LGA')
            fprintf('%s\n','WMH segmentation performed using LST LGA');
            % Lesion Growth Algorithm (LGA)
            rWMHPath= xASL_wrp_LST_lga(x);

        elseif ~strcmp(x.WMHsegmAlg,'LPA') || ~strcmp(x.WMHsegmAlg,'LGA')
                error('Wrong WMH segmentation defined -> x.WMHsegmAlg should be LPA or LGA');
        end

        fprintf('\n');
        xASL_io_ConvertNii16bit(rWMHPath);

        x.mutex.AddState('040_segment_FLAIR');
        %xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState('050_LesionFilling');

    elseif bO; fprintf('%s\n','040_segment_FLAIR: No FLAIR data found, skipping...');
    end
elseif  bO; fprintf('%s\n','040_segment_FLAIR has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 5    Lesion filling
if ~x.mutex.HasState('050_LesionFilling')  % tracks progress through lock/ *.status files, & locks current run

    if ~xASL_exist(rWMHPath,'file') && bO
        % check if the FLAIR segmentation exists from LST

        fprintf('%s\n','050_LesionFilling: No FLAIR data found, skipping...');
    elseif ~xASL_exist(rWMHPath,'file')
        % Just skip it without notifying

    else

        % Here we should choose between LGA and LPA. LGA seems more robust,
        % and also does a good job on removing the FLAIR biasfield.
        % LPA is newer, and quicker

        if  xASL_exist(x.P.Path_WMH_rSEGM,'file')
            % We assume that the pre-existing WMH_SEGM.nii contains a
            % better segmentation than the LST segmentation,
            % therefore replace the LST segmentation by the WMH_SEGM.nii

            FaultyIM                = xASL_io_ReadNifti(rWMHPath);
            CorrectIM               = xASL_io_Nifti2Im(x.P.Path_WMH_rSEGM);

            if  min(size(FaultyIM.dat)==size(CorrectIM))
                % first verify whether image sizes are identical
                xASL_save_nii(rWMHPath,rWMHPath,CorrectIM,[],0);

            end
        end


        fprintf('\n%s\n','----------------------------------------');
        fprintf('%s\n','Removing segmented WMH from T1w, and fill lesions with values interpolated from neighborhood');

        % Create a copy of the WMH segmentation
       if ~xASL_exist(x.P.Path_WMH_PreSEGM, 'file')
            xASL_Copy(rWMHPath, x.P.Path_WMH_PreSEGM);
        end

        % Clip the WMH segmentation used for lesion filling
        IM = xASL_io_Nifti2Im(rWMHPath);
        IM(IM>=0.5) = 1;
        IM(IM<0.5) = 0;
        SumIm = xASL_st_sumnan(IM(:));
        xASL_save_nii(rWMHPath, rWMHPath, IM, 16, 0);

        if SumIm>0
            %clear matlabbatch
            %matlabbatch{1}.spm.tools.LST.filling.data           = {x.P.Path_T1};
            %matlabbatch{1}.spm.tools.LST.filling.data_plm       = {x.P.rWMHPath};
            %matlabbatch{1}.spm.tools.LST.filling.html_report    = 0; % saves time
            %spm_jobman('run',matlabbatch); close all
            xASL_wrp_LST_lesfill( x );
        end

        xASL_delete(rWMHPath);

        x.mutex.AddState('050_LesionFilling');
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState('060_Get_WMH_vol');
        x.mutex.DelState('070_segment_T1w');
    end

elseif  bO; fprintf('%s\n','050_LesionFilling has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 6    Save WMH volume
if ~x.mutex.HasState('060_Get_WMH_vol')  % tracks progress through lock/ *.status files, & locks current run
    if  xASL_exist(x.P.Path_rFLAIR,'file')  && xASL_exist(rWMHPath,'file')

        % This function prints the WMH volume (mL) and number of WMH in a CSV file
        clear matlabbatch
        matlabbatch{1}.spm.tools.LST.tlv.data_lm            = {x.P.Path_WMH_rSEGM};
        BinThresh                                           = 0.5; % default setting
        matlabbatch{1}.spm.tools.LST.tlv.bin_thresh         = BinThresh;
        spm_jobman('run',matlabbatch); close all

        Fname       = GetFileList(x.SUBJECTDIR, ['^LST_tlv_' num2str(BinThresh) '.*\.csv$']);

        %oPath       = fullfile(x.D.TissueVolumeDir, ['WMH_LST_' x.WMHsegmAlg '_' x.P.SubjectID '.csv']);
        %saved in the subj folder instead of Tissue Volume Dir in population
        %folder
        oPath       = fullfile(SUBJECTDIR, BIDS_FileName(x.P.SubjectID,SessionID,'','WMH',x.WMHsegmAlg ,'','','cvs'));
        MoveFileIf(Fname{1},oPath,1 );
        xASL_bids_csv2tsvReadWrite(oPath,1); % convert to tsv per BIDS

        x.mutex.AddState('060_Get_WMH_vol');
        xASL_adm_CompareDataSets([], [], x); % unit testing
    elseif  bO; fprintf('%s\n','060_Get_WMH_vol: No FLAIR data found, skipping...');
    end

elseif  bO; fprintf('%s\n','060_Get_WMH_vol has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 6.5    File management -> move within separate scripts

% Here, the above WMH segmentation is only renamed to WMH_SEGM.nii if no
% other WMH_SEGM.nii already existed

% if LGA, clean up afterwards
%CleanUpDir      = GetFsList(x.SUBJECTDIR,'^LST_tmp_.*$',1,[],[],[0 Inf]);
%nList           = length(CleanUpDir);
for iD=1:nList
    DeleteFileList(fullfile(x.SUBJECTDIR,CleanUpDir{iD}), '^.*\.(nii|nii\.gz)$', [], [0 Inf]);
    DeleteFileList(fullfile(x.SUBJECTDIR,CleanUpDir{iD}), '^.*\.mat$', [], [0 Inf]);
    rmdir( fullfile(x.SUBJECTDIR,CleanUpDir{iD}),'s' );
end

if  ~xASL_exist(x.P.Path_WMH_PreSEGM,'file') && xASL_exist(x.P.Path_WMH_rSEGM,'file')
    % We assume here that any existing WMH_SEGM.nii is better than the one segmented by LST above.
    % Therefore, only if no WMH_SEGM.nii exists, this will be created from the LGA or LPA
    MoveFileIf(x.P.Path_WMH_rSEGM, x.P.Path_WMH_PreSEGM, 1);
end

% if  xASL_exist(x.P.Path_T1,'file') && ~xASL_exist(x.P.Path_T1_ORI)
%     % If T1w backup was not created through biasfield correction, create it now
%     CopyFileIf(x.P.Path_T1 , x.P.Path_T1_ORI , 1);
% end

if  xASL_exist(x.P.Path_T1,'file') && xASL_exist(T1_filledName,'file')
    MoveFileIf(T1_filledName, x.P.Path_T1,1);
end

WMHmat              = GetFileList(fPath,'^LST_.*FLAIR\.mat$','FPList',[0 Inf]);

if  x.DELETETEMP && ~x.bReproTesting
    xASL_delete(x.P.Path_WMH_rSEGM);
    % WMH segmentation, if already a WMH_SEGM.nii existed, we assume
    % this was better and delete the LST segmentation.
    % The LST segmentation is then only used for T1 lesion filling


    if ~isempty(WMHmat) % LST mat-file
        xASL_delete(WMHmat{1});
    end

    xASL_delete( x.P.Path_rFLAIR);  % LPA
    xASL_delete( x.P.Path_mrFLAIR); % LPA
    xASL_delete(x.P.Path_rmrFLAIR); % LGA
end










%% -----------------------------------------------------------------------------
%% 7    CAT12 segmentation
% CAT12 performs even better than SPM12. Therefore, always run CAT12, unless
% this crashes, then we can run SPM12

if ~isfield(x,'Segment_SPM12')
    x.Segment_SPM12 = 0;
end
% by default, use CAT12, not SPM12 for segmentation

if ~x.mutex.HasState('070_segment_T1w') || ~xASL_exist(x.P.Path_c1T1,'file')  || ~xASL_exist(x.P.Path_c2T1,'file')

    x = xASL_wrp_segment( x );

    x.mutex.AddState('070_segment_T1w');  % tracks progress through lock/ *.status files, & locks current run
    xASL_adm_CompareDataSets([], [], x); % unit testing
    x.mutex.DelState('080_TissueVolume')
    x.mutex.DelState('090_reslice2DARTEL');

elseif  bO; fprintf('%s\n','070_segment_T1w has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 8    Get tissue volumes
if ~x.mutex.HasState('080_TissueVolume')   % tracks progress through lock/ *.status files, & locks current run

    % PM: this should be corrected for any
    % Lesion_(T1|FLAIR)_(1|2|3|..)\.nii files

    if ~x.mutex.HasState('070_segment_T1w')
        fprintf('%s\n','070_segment_T1w not performed, 080_TissueVolume skipped');
    else
        % run only if segmentation has been run

        catVolFile  = fullfile(x.D.TissueVolumeDir,['cat_' x.P.STRUCT '_' x.P.SubjectID '.mat']);
        mat_file    = fullfile( x.SUBJECTDIR      , [x.P.STRUCT '_seg8.mat']);

        if  xASL_exist(catVolFile,'file') % for CAT12 segmentation
            catVol  = load(catVolFile);
            GMvol  = catVol.S.subjectmeasures.vol_abs_CGW(2)/1000; % volume was in cm^3, /1000 = dm^3 or Liter
            WMvol  = catVol.S.subjectmeasures.vol_abs_CGW(3)/1000;
            CSFvol = catVol.S.subjectmeasures.vol_abs_CGW(1)/1000;

            SaveFile= fullfile(x.D.TissueVolumeDir,['TissueVolume_' x.P.SubjectID '.csv']);
            %SaveFile= fullfile(x.D.TissueVolumeDir,BIDS_FileName(['TissueVolume_' x.P.SubjectID '.csv']));
            FileID  = fopen(SaveFile,'wt');
            fprintf(FileID,'%s\n', 'File,GM volume (L),WM volume (L),CSF volume (L)');
            fprintf(FileID,'%s', [catVolFile ',' num2str(GMvol) ',' num2str(WMvol) ',' num2str(CSFvol)]);
            fclose(FileID);


        elseif exist(mat_file,'file') % for SPM12 segmentation

            % CAVE: this only works in SPM12


            % Make sure that filename is correct, otherwise this crashes.
            % This is the case if the files have been run previously in
            % different folder structure
            load(mat_file,'-mat');
            for ii=1:6
                tpm(ii).fname = fullfile(x.SPMDIR, 'tpm', 'TPM.nii');
            end
            image.fname                 = fullfile( x.SUBJECTDIR, [x.P.STRUCT '.nii']);

            save(mat_file,'image','tpm','Affine','lkp','MT','Twarp','Tbias','mg','mn','vr','wp','ll');
            clear image tpm Affine lkp MT Twarp Tbias mg mn vr wp ll volumes

            job         = fullfile( x.D.FunctionsDir, 'xASL_spm_job_TissueVolume.m');
            mask_ICV    = GetImageList3D( fullfile(x.SPMDIR, 'tpm'), '^mask_ICV\.(nii|nii\.gz)$');
            SaveFile    = fullfile( x.D.TissueVolumeDir, ['TissueVolume_' x.P.SubjectID '.csv']);

            spm_jobman('serial', job, '', {mat_file}, mask_ICV, SaveFile);
        end

        if  exist('SaveFile','var')
            xASL_bids_csv2tsvReadWrite(SaveFile,1);
        end

        x.mutex.AddState('080_TissueVolume');
        xASL_adm_CompareDataSets([], [], x); % unit testing
    end
elseif  bO; fprintf('%s\n','080_TissueVolume has already been performed, skipping...');
end




%% -----------------------------------------------------------------------------
%% 9    Reslice to DARTEL folder
if ~x.mutex.HasState('090_reslice2DARTEL')  % tracks progress through lock/ *.status files, & locks current run
    if ~x.mutex.HasState('070_segment_T1w')
        fprintf('%s\n','070_segment_T1w not performed, 090_reslice2DARTEL skipped');
    else

        % Run only if has not been run, if segmentation has been run &
        % rigid-body registrations (for all volumes) have been run

        fprintf('%s\n','Reslice structural images & move to DARTEL folder');

        % Define pathname list of structural images to be resliced
        INname         = {x.P.Path_T1; x.P.Path_c1T1; x.P.Path_c2T1; x.P.Path_FLAIR};
        OUTname        = {x.P.Pop_Path_rT1; x.P.Pop_Path_rc1T1; x.P.Pop_Path_rc2T1; x.P.Pop_Path_rFLAIR};
        xASL_spm_deformations(x,INname,OUTname);

        % Lesion probability maps (do linear interpolation to avoid negative edge effects)
        xASL_spm_deformations(x,x.P.Path_WMH_PreSEGM,x.P.Pop_Path_rWMH_SEGM,1);

        [INname,OUTname]     = xASL_adm_LesionResliceList(x, Lesion_T1_list, Lesion_FLAIR_list,ROI_T1_list,ROI_FLAIR_list);

        if  ~isempty(INname) && ~isempty(OUTname)

            % First dilate ROIs, if they were e.g. used for annotation (single voxel only)
            % Do linear interpolation to avoid negative edge effects
            for iL=1:length(INname)
                xASL_im_dilateROI(INname{iL});
            end

            xASL_spm_deformations(x,INname,OUTname,1);
        end

        x.mutex.AddState('090_reslice2DARTEL');
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState('100_visualize');
    end

elseif  bO; fprintf('%s\n','090_reslice2DARTEL has already been performed, skipping...');
end




%% -----------------------------------------------------------------------------
%% 10    Visual QC
if ~x.mutex.HasState('100_visualize')  % tracks progress through lock/ *.status files, & locks current run
    if  ~x.mutex.HasState('070_segment_T1w') || ~x.mutex.HasState('090_reslice2DARTEL')
        fprintf('%s\n','070_segment_T1w or 090_reslice2DARTEL not performed, visual QC skipped');
    else

        xASLMAT  = fullfile(x.D.ROOT,'x.mat'); % later optimize this, do this in each xWrapper
        if  xASL_exist(xASLMAT,'file')
            Oldx          = load(xASLMAT);
            x.Output_im   = Oldx.x.Output_im;
			if isfield(Oldx.x,'Output')
				x.Output      = Oldx.x.Output;
			end
        end

        %% Run the SPM Univariate Plus (UP) QC for T1w
        anatQA = xASL_spmup_anatQA(x.P.Path_T1, x.P.Path_c1T1, x.P.Path_c2T1, x);
        x = xASL_adm_AddToQC(x, iSubj, anatQA, 'T1w');

        if  xASL_exist(x.P.Path_FLAIR,'file')
            % do the same for FLAIR
            xASL_spm_reslice(x.P.Path_FLAIR, x.P.Path_c1T1, [], [], x.Quality);
            xASL_spm_reslice(x.P.Path_FLAIR, x.P.Path_c2T1, [], [], x.Quality);

            anatQA = xASL_spmup_anatQA(x.P.Path_FLAIR, x.P.Path_rc1T1, x.P.Path_rc2T1, x);
            x = xASL_adm_AddToQC(x, iSubj, anatQA, 'FLAIR');

            xASL_delete(x.P.Path_rc1T1);
            xASL_delete(x.P.Path_rc2T1);
        end

        fprintf('%s\n','Saving QC images');

        % This function contains several visualizations, that are re-performed after longitudinal registration
        % So this is also called by longitudinal registration, therefore keep this function separate
        x = xASL_wrp_VisualCheckCollective_Structural(x);

        %% Visualize lesions
        xASL_wrp_VisualCheckLesionRemoval( x, Lesion_T1_list, Lesion_FLAIR_list);
        xASL_im_VisualizeROIs( x, ROI_T1_list, ROI_FLAIR_list);

        % Convert ROIs & lesions to specific masks
        for iS=1:length(Lesion_T1_list)
            xASL_im_Lesion2Mask( fullfile(x.D.PopDir, ['rLesion_' x.P.STRUCT '_' num2str(iS) '_' x.P.SubjectID '.nii']), [], [], [], x );
        end
        for iS=1:length(ROI_T1_list)
            xASL_im_Lesion2Mask( fullfile(x.D.PopDir, ['rROI_' x.P.STRUCT '_' num2str(iS) '_' x.P.SubjectID '.nii']), [], [], [], x );
        end
        for iS=1:length(Lesion_FLAIR_list)
            xASL_im_Lesion2Mask( fullfile(x.D.PopDir, ['rLesion_' x.P.FLAIR '_' num2str(iS) '_' x.P.SubjectID '.nii']), [], [], [], x );
        end
        for iS=1:length(ROI_FLAIR_list)
            xASL_im_Lesion2Mask( fullfile(x.D.PopDir, ['rROI_' x.P.FLAIR '_' num2str(iS) '_' x.P.SubjectID '.nii']), [], [], [], x );
        end

        % Show lesions individually
        for iS=1:length(Lesion_T1_list)
            xASL_im_CreateVisualFig( x, {x.P.Pop_Path_rT1, fullfile(x.D.PopDir,['rLesion_' x.P.STRUCT '_' num2str(iS) '_' x.P.SubjectID '.nii'])},x.D.LesionCheckDir,[0.8 1],[],[]);
        end
        for iS=1:length(Lesion_FLAIR_list)
            xASL_im_CreateVisualFig( x, {x.P.Pop_Path_rFLAIR, fullfile(x.D.PopDir,['rLesion_' x.P.FLAIR '_' num2str(iS) '_' x.P.SubjectID '.nii'])},x.D.LesionCheckDir,[0.8 1],[],[]);
        end
        % Visualize ROIs (these are manually added native space ROIs)
        for iS=1:length(ROI_T1_list)
            xASL_im_CreateVisualFig( x, {x.P.Pop_Path_rT1, fullfile(x.D.PopDir,['rROI_' x.P.STRUCT '_' num2str(iS) '_' x.P.SubjectID '.nii'])},x.D.ROICheckDir,[0.8 1],[],[]);
        end
        for iS=1:length(ROI_FLAIR_list)
            xASL_im_CreateVisualFig( x, {x.P.Pop_Path_rFLAIR, fullfile(x.D.PopDir,['rROI_' x.P.FLAIR '_' num2str(iS) '_' x.P.SubjectID '.nii'])},x.D.ROICheckDir,[0.8 1],[],[]);
        end

        x = xASL_qc_CollectParameters(x, iSubj, 0); % Quick & Dirty solution -> 0 indicates no ASL results
        save(xASLMAT,'x'); % future: do this in each xWrapper

        x.mutex.AddState('100_visualize');
        xASL_adm_CompareDataSets([], [], x); % unit testing
    end

    xASL_adm_CreateFileReport(x);
%   This function checks which files are existing/missing, providing overview of correctly imported files (native space)
%   & correctly processed files (standard space)

    xASL_qc_PrintOrientation(x.D.ROOT, ['^' x.P.STRUCT '\.(nii|nii\.gz)$'], x.D.ROOT, 'RigidRegT1');
    % This function summarizes the T1w orientation. Especially check the determinant, for left-right flips

    xASL_qc_CreateOutputPDF(x, iSubj);

elseif  bO; fprintf('%s\n','100_visualize has already been performed, skipping...');
end





%% -----------------------------------------------------------------------------
%% 11   Deface for anonymity
% SPM version. Later this can be improved by x version, as registration
% has already been performed, to only mask the T1, rather than having to
% register first, and with a smooth mask, rather than hard clipping at the
% edge, to make it more robust for slight misalignment.

%% This part should go to the import section

if ~isfield(x,'Deface')
    x.Deface=0; % by default, turn off (quick & dirty fix)
end

if ~x.mutex.HasState('110_Deface')  % tracks progress through lock/ *.status files, & locks current run
    if ~x.Deface
        if bO; fprintf('%s\n','Defacing skipped'); end
    else
        xASL_spm_deface( x.P.Path_T1); % this requires improvement, as sometimes too much face is removed. use gradient
        x.mutex.AddState('110_Deface');
        xASL_adm_CompareDataSets([], [], x); % unit testing
    end
elseif  bO; fprintf('%s\n','110_Deface has already been performed, skipping...');
end




%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
cd(oldFolder);
result = true;

% xASL_delete(QName);


end


function [x] = xASL_adm_AddToQC(x, iSubj, anatQA, Modality)
%xASL_adm_AddToQC % Add SPM U+ parameters to the QC list

    FN = fieldnames(anatQA);
    FNnew = cellfun(@(x) [Modality '_' x], FN, 'UniformOutput',false);

    for iL=1:length(FN)
        x.Output.Structural(iSubj).(FNnew{iL}) = anatQA.(FN{iL});
    end

end
