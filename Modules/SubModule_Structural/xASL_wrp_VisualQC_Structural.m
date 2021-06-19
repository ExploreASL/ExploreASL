function xASL_wrp_VisualQC_Structural(x)
%xASL_wrp_VisualQC_Structural Submodule of ExploreASL Structural Module, that performs several visualizations for QC
%
% FORMAT: xASL_wrp_VisualQC_Structural(x)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule performs several visualizations for visual & quantitative QC.
%
% 1. After initial admin
% 2. It starts with the SPM UP parameters (courtesy of Cyril Pernet, his SPM UP scripts were
%    made more robust & accurate by Jan & Henk, & are implemented here for T1w (& optionally FLAIR).
% 3. Then it performs a collection of visualizations
% 4. Also repeated specifically for lesions & manually provided ROIs
% 5. Finally, this contains a report of all missing raw & derivative files, in native & standard space,
%    printing the NIfTI orientation matrix content before (hdr.mat0) & after registrations (hdr.mat)
%    The determinant of these matrices should be the same, otherwise LeftRight has flipped. This should
%    also be the same across a group scanned at the same scanner
%    Then various other QC functions are called & all are summarized in a PDF report.
%
% EXAMPLE: xASL_wrp_VisualQC_Structural(x);
% __________________________________
% Copyright (C) 2015-2020 ExploreASL


%% -----------------------------------------------------------------------------------
%% 1) Admin
PathX = fullfile(x.dir.SUBJECTDIR,'x.mat');
iSubject = find(strcmp(x.SUBJECTS, x.P.SubjectID)); % Find current subject index
x = xASL_adm_LoadX(x, PathX, true); % assume x.mat is newer than x

% Clear the current QC images
if isfield(x,'Output_im') && isfield(x.Output_im,'Structural')
    x.Output_im = rmfield(x.Output_im,'Structural');
end
if isfield(x,'Output') && isfield(x.Output,'Structural')
    x.Output = rmfield(x.Output,'Structural');
end

%% -----------------------------------------------------------------------------------
%% 2) Run SPM Univariate Plus (UP) QC for T1w (& FLAIR if exists)
% These files need to be inversely transformed from MNI to native space.
% To speed up, we reslice them instead of transforming from MNI, if they
% already exist, which is the case when doing this for the FLAIR, after the
% T1w, below.

% First, delete any from a previous run
Path_NativeDeepWM = fullfile(x.dir.SUBJECTDIR,'CentralWM_QC.nii');
Path_NativeLRMask = fullfile(x.dir.SUBJECTDIR,'LeftRight.nii');
Path_MNI_DeepWM = fullfile(x.D.PopDir,['CentralWM_QC_' x.P.SubjectID '.nii']);
xASL_delete(Path_NativeDeepWM);
xASL_delete(Path_NativeLRMask);
xASL_delete(Path_MNI_DeepWM);

PCP_T1w = xASL_qc_PCPStructural(x.P.Path_T1, x.P.Path_c1T1, x.P.Path_c2T1, x, x.P.Pop_Path_rT1);
x = xASL_adm_AddToQC(x, PCP_T1w, 'T1w');

if xASL_exist(x.P.Path_FLAIR,'file') % do the same for FLAIR
    xASL_spm_reslice(x.P.Path_FLAIR, x.P.Path_c1T1, [], [], x.Quality, [], 1);
    xASL_spm_reslice(x.P.Path_FLAIR, x.P.Path_c2T1, [], [], x.Quality, [], 1);

    PCP_FLAIR = xASL_qc_PCPStructural(x.P.Path_FLAIR, x.P.Path_rc1T1, x.P.Path_rc2T1, x, x.P.Pop_Path_rFLAIR);
    x = xASL_adm_AddToQC(x, PCP_FLAIR, 'FLAIR');

    xASL_delete(x.P.Path_rc1T1);
    xASL_delete(x.P.Path_rc2T1);
end

fprintf('%s\n','Saving QC images');


%% -----------------------------------------------------------------------------------
%% 3) Perform several visualizations
x = xASL_wrp_VisualCheckCollective_Structural(x);



%% -----------------------------------------------------------------------------------
%% 4) Visualize lesions
Lesion_list = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['(?i)^Lesion_(' x.P.STRUCT '|' x.P.FLAIR ')_\d*\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);
ROI_list = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['^ROI_(' x.P.STRUCT '|' x.P.FLAIR ')_\d*\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);

xASL_adm_VisualCheckLesionRemoval(x, Lesion_list);
xASL_vis_VisualizeROIs(x, ROI_list);

% Show lesions individually
for iLesion=1:length(Lesion_list)
    [~, Ffile] = xASL_fileparts(Lesion_list{iLesion});
    LesionFile = fullfile(x.D.PopDir, ['r' Ffile '_' x.P.SubjectID '.nii']);
    xASL_im_Lesion2Mask(LesionFile, x); % Convert ROIs & lesions to specific masks
    xASL_vis_CreateVisualFig(x, {LesionFile}, x.D.LesionCheckDir,[0.8 1], 'Lesions'); % Show lesions individually
end
% Visualize ROIs (these are manually added native space ROIs)
for iROI=1:length(ROI_list)
    [~, Ffile] = xASL_fileparts(ROI_list{iROI});
    ROIFile = fullfile(x.D.PopDir, ['r' Ffile '_' x.P.SubjectID '.nii']);
    xASL_im_Lesion2Mask(ROIFile, x); % Convert ROIs & lesions to specific masks    
    xASL_vis_CreateVisualFig(x, {ROIFile}, x.D.ROICheckDir,[0.8 1], 'Lesions'); % Show lesions individually
end

%% -----------------------------------------------------------------------------------
%% 5) Final QCs
xASL_qc_PrintOrientation(x.dir.SUBJECTDIR, x.P.Path_T1, x.dir.SUBJECTDIR, 'RigidRegT1');
% This function summarizes the T1w orientation. Especially check the determinant, for left-right flips

x = xASL_qc_CollectParameters(x, iSubject, 'Structural');
xASL_delete(PathX);
xASL_adm_SaveX(x);

xASL_qc_CreatePDF(x, iSubject);


end




%% ========================================================================
%% ========================================================================
function [x] = xASL_wrp_VisualCheckCollective_Structural(x)
%xASL_wrp_VisualCheckCollective_Structural Runs a collection of visual QC
%functions

T.ModuleName = 'Structural';
Path_CentralWM_QC = fullfile(x.D.PopDir,['CentralWM_QC_' x.P.SubjectID '.nii']);
x = xASL_adm_ResetVisualizationSlices(x);

%% T1w before lesion filling
if xASL_exist(x.P.Path_T1_ORI, 'file') % if T1w was lesion-filled or bias-field corrected
    xASL_spm_deformations(x,x.P.Path_T1_ORI,x.P.Pop_Path_rT1_ORI);
end


%% Parameters for creating visual QC Figures
%  General parameters (empty cells are default)
T.IntScale    = {[] [1 1] [1 1] [] [] [1 1] [1 1] []};
T.Suffix      = {'Src' 'Reg' 'Seg' 'QC_ROI'};
T.Suffix(5:8) = cellfun(@(x) ['Cor_' x],T.Suffix(1:4),'UniformOutput',false);
T.Suffix(1:4) = cellfun(@(x) ['Tra_' x],T.Suffix(1:4),'UniformOutput',false);
T.ColorMapIs  = {[] [] [] {x.S.gray x.S.jet256} [] [] [] {x.S.gray x.S.jet256}};
T.TraSlices   = {[] [] [] [] 'n/a' 'n/a' 'n/a' 'n/a'}; % repmat(x.S.slicesLarge(8),[1 12])
T.CorSlices   = {'n/a' 'n/a' 'n/a' 'n/a' [] [] [] []}; % repmat(x.S.slicesLarge(8),[1 12])
T.bClip       = {[] [] [1 0] [1 0] [] [1 0] [] [1 0]};

if  xASL_exist(x.P.Path_T1,'file')
    % T1w parameters
    % Source T1w  , T1+WMH filling, T1w+WMmask  , T1w+WMrefMask;
    T.ImIn        = {{x.P.Pop_Path_rT1_ORI}     {x.P.Pop_Path_rT1}    {x.P.Pop_Path_rT1 x.P.Pop_Path_rc2T1} {x.P.Pop_Path_rT1 Path_CentralWM_QC}};
    T.ImIn(5:8)   = T.ImIn(1:4);
    T.DirOut(1:8) = {x.D.T1CheckDir};
    x = CreateVisual(x,T);
end

if  xASL_exist(x.P.Path_FLAIR,'file')
    % Same as above but then for FLAIR
    T.ImIn        = {{x.P.Pop_Path_rFLAIR} {x.P.Pop_Path_rFLAIR x.P.Pop_Path_rWMH_SEGM} {x.P.Pop_Path_rFLAIR x.P.Pop_Path_rc2T1} {x.P.Pop_Path_rFLAIR Path_CentralWM_QC}};
    T.ImIn(5:8)   = T.ImIn(1:4);
    T.DirOut(1:8) = {x.D.FLAIR_CheckDir};

    x = CreateVisual(x, T);
end

if  xASL_exist(x.P.Path_T2,'file')
	% It is created only if the T2 file exists
	xASL_adm_CreateDir(x.D.T2_CheckDir);
	
    % Same as above but then for T2
    T.ImIn        = {{x.P.Pop_Path_rT2} {x.P.Pop_Path_rT2 x.P.Pop_Path_rc2T1}};
    T.ImIn(3:4)   = T.ImIn(1:2);
    T.DirOut(1:4) = {x.D.T2_CheckDir};

    x = CreateVisual(x, T);
end

if  xASL_exist(x.P.Path_T1c,'file')
	% It is created only if the T1c file exists
	xASL_adm_CreateDir(x.D.T1c_CheckDir);
	
    % Same as above but then for T1c
    T.ImIn        = {{x.P.Pop_Path_rT1c} {x.P.Pop_Path_rT1c x.P.Pop_Path_rc2T1}};
    T.ImIn(3:4)   = T.ImIn(1:2);
    T.DirOut(1:4) = {x.D.T1c_CheckDir};

    x = CreateVisual(x, T);
end

xASL_delete(Path_CentralWM_QC); % householding

end


%% ========================================================================
%% ========================================================================
function x = CreateVisual(x, T)
%CreateVisual This part creates the figure

    for iM=1:length(T.ImIn)
        % First check file existence (otherwise fill with empty image?)
        FExistence   = 1;
        for iC=1:length(T.ImIn{iM})
            if ~xASL_exist(T.ImIn{iM}{iC})
                FExistence   = 0;
            end
        end

        %% Manage slices to show
        if FExistence
            % Sagittal
            x.S.SagSlices = []; % show no sagittal slices
            % Transversal
            TraSlicesDefault = x.S.slicesLarge;
            TraSlicesDefault(8) = 68; % 1 slice lower

            if isempty(T.TraSlices{iM})
                    x.S.TraSlices   = TraSlicesDefault;
            elseif  strcmp(T.TraSlices{iM},'n/a')
                    x.S.TraSlices = [];
            elseif  isnumeric(T.TraSlices{iM})
                    x.S.TraSlices = T.TraSlices{iM};
            else
                    warning('Wrong slice choice');
            end
            % Coronal
            CorSlicesDefault = x.S.slicesLarge+7;
            if isempty(T.CorSlices{iM})
                    x.S.CorSlices = CorSlicesDefault;
            elseif  strcmp(T.CorSlices{iM},'n/a')
                    x.S.CorSlices = [];
            elseif  isnumeric(T.CorSlices{iM})
                    x.S.CorSlices = T.CorSlices{iM};
            else
                    warning('Wrong slice choice');
            end

            % Create the image
            T.IM = xASL_vis_CreateVisualFig(x, T.ImIn{iM}, T.DirOut{iM}, T.IntScale{iM}, T.Suffix{iM}, T.ColorMapIs{iM}, T.bClip{iM});
            
            T.bCrop = false;
            X = size(T.IM,1); Y = size(T.IM,2);
            if strcmp(T.CorSlices{iM},'n/a')
                % This is a transversal slice
                T.IM = squeeze(T.IM(ceil(0.33*X)+2:floor(0.67*X)-1,ceil(Y*3/4+1):end,:)); % slice 8
            elseif strcmp(T.TraSlices{iM},'n/a')
                % This is a coronal slice
                T.IM = squeeze(T.IM(ceil(0.33*X)+2:floor(0.67*X)-1,ceil(Y/4+1):floor(Y/2),:)); % slice 6
            end

            x = xASL_vis_AddIM2QC(x,T);
        end
    end

end



%% ========================================================================
%% ========================================================================
function [x] = xASL_adm_AddToQC(x, anatQA, Modality)
%xASL_adm_AddToQC % Add SPM U+ parameters to the QC list

    FN = fieldnames(anatQA);
    FNnew = cellfun(@(x) [Modality '_' x], FN, 'UniformOutput',false);

    for iL=1:length(FN)
        x.Output.Structural.(FNnew{iL}) = anatQA.(FN{iL});
    end

end


%% ========================================================================
%% ========================================================================
function xASL_adm_VisualCheckLesionRemoval(x, Lesion_list)
% xASL_adm_VisualCheckLesionRemoval Creates for each subject  a JPEG image containing
% the segmented image with removed lesion, T1+lesion ROI, segmentation before removal

    if ~isempty(Lesion_list)

        % c1T1 & c2T1 before lesion masking
        if  xASL_exist(x.P.Path_c1T1_ORI,'file') && xASL_exist(x.P.Path_c2T1_ORI,'file') % if T1w was lesion-filled or bias-field corrected
            xASL_io_ReadNifti(x.P.Path_c1T1_ORI); % unzip if needed
            xASL_io_ReadNifti(x.P.Path_c2T1_ORI); % unzip if needed
            xASL_spm_deformations(x,{x.P.Path_c1T1_ORI;x.P.Path_c2T1_ORI},{x.P.Pop_Path_rc1T1_ORI;x.P.Pop_Path_rc2T1_ORI}); % no DARTEL yet, no LongReg yet
        end

        LesionIM    = zeros(121,145,121); % assuming 1.5 mm MNI

		% Call this for Lesions only and not ROIs
        [INname,OUTname]     = xASL_adm_LesionResliceList(x,1,1,0,0);

        if ~isempty(INname) && ~isempty(OUTname)
            % First dilate ROIs, if they were e.g. used for annotation (single voxel only)
            % Do linear interpolation to avoid negative edge effects
            for iO=1:length(OUTname)
                if ~xASL_exist(OUTname{iO})
                    for iLesion=1:length(INname)
                        xASL_im_dilateROI(INname{iLesion}, [], 40);
                    end
                    xASL_spm_deformations(x, INname, OUTname, 1);
                end
            end
        end

        for iLesion=1:length(Lesion_list)
            [~, Ffile] = xASL_fileparts(Lesion_list{iLesion});
            LesionFile = fullfile(x.D.PopDir, ['r' Ffile '_' x.P.SubjectID '.nii']);
            TempIM = xASL_io_Nifti2Im(LesionFile);
            LesionIM = min(1,LesionIM+xASL_im_ConvertMap2Mask(TempIM(:,:,:,1)));
        end
        
		ImOut1 = xASL_vis_CreateVisualFig( x, {x.P.Pop_Path_rc1T1_ORI x.P.Pop_Path_rc2T1_ORI}, [], [1 0.8], [], []);
		ImOut2 = xASL_vis_CreateVisualFig( x, {x.P.Pop_Path_rT1 LesionIM},[], [1 0.8], [], []);
		ImOut3 = xASL_vis_CreateVisualFig( x, {x.P.Pop_Path_rc1T1 x.P.Pop_Path_rc2T1},[], [1 0.8], [], []);

        IM = [ImOut1,ImOut2,ImOut3];
        xASL_adm_CreateDir(x.D.LesionCheckDir);
        xASL_vis_Imwrite((IM+eps)./max(IM(:)), fullfile(x.D.LesionCheckDir , ['Lesion_corr_' x.P.SubjectID '.jpg']));

    end
end
