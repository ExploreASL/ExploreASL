function xASL_wrp_VisualQC_ASL(x)
%xASL_wrp_VisualQC_ASL Submodule of ExploreASL ASL Module, that performs several visualizations for QC
%
% FORMAT: xASL_wrp_VisualQC_ASL(x)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule performs several visualizations for visual & quantitative QC.
%              1) After initial admin
%              2) It starts with making ASL NIfTIs ready for visualization & conversion to DICOM
%              3) Then it performs a collection of visualizations
%              4) Visualizes results of the TopUp geometric distortion correction
%              5) Visualization of slice gradient
%              6) Visualization & calculation of temporal QC parameters
%              7) Compute DICE overlap/intersection of ASL brain in FoV & T1w, to calculate coverage
%              8) Summarize orientation & check left-right flips
%              9) Collect several other parameters & store in PDF overview
%
% EXAMPLE: xASL_wrp_VisualQC_ASL(x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


%% -----------------------------------------------------------------------------------
%% 1) Admin
% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end
tempnii = xASL_io_ReadNifti(x.P.Path_despiked_ASL4D);
nVolumes = double(tempnii.hdr.dim(5));

fprintf('%s\n','print visual quality assurance checks');
Parms.ModuleName = 'ASL';
close all % close all Figures to avoid capturing & saving the wrong Figure
PathX = fullfile(x.SUBJECTDIR,'x.mat');
x = xASL_adm_LoadX(x, PathX, true); % assume x.mat is newer than x

% Clear any previous QC images
if isfield(x,'Output_im') && isfield(x.Output_im,'ASL')
   x.Output_im = rmfield(x.Output_im,'ASL');
end
if isfield(x,'Output') && isfield(x.Output,'ASL')
   x.Output = rmfield(x.Output,'ASL');
end

%% -----------------------------------------------------------------------------------
%% 2) Make ASL NIfTIs ready for visualization & conversion to DICOM
if xASL_exist(x.P.Path_T1_ORI, 'file')
    InputT1Oripath = x.P.Path_T1_ORI;
else
    InputT1Oripath = x.P.Path_T1;
end
   
xASL_io_MakeNifti4DICOM(x.P.Path_CBF, x, 'UINT16', 1, x.P.Path_T1); % Create uint16 NIfTI in 12 bit scale
% xASL_io_MakeNifti4DICOM(x.P.Pop_Path_qCBF, x);


%% -----------------------------------------------------------------------------------
%% 3) Perform several visualizations
x = xASL_wrp_VisualCheckCollective_ASL(x);

%% -----------------------------------------------------------------------------------
%% 4) Visualization TopUp results (quick & dirty)
PathPopB0 = fullfile(x.D.PopDir, ['rASL_B0_' x.P.SubjectID '_' x.P.SessionID '.nii']);
PathPopUnwarped = fullfile(x.D.PopDir, ['rASL_Unwarped_' x.P.SubjectID '_' x.P.SessionID '.nii']);

if xASL_exist(PathPopB0,'file') && xASL_exist(PathPopUnwarped,'file')% if we have TopUp results
    [Output1, Output2] = xASL_im_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, x.iSubject, x.D.ASLCheckDir);
    x.Output.ASL.MeanAI_PreTopUp_Perc = Output1;
    x.Output.ASL.MeanAI_PostTopUp_Perc = Output2;
    xASL_delete(PathPopB0);
    xASL_delete(PathPopUnwarped);
end

%% -----------------------------------------------------------------------------------
%% 5) Visualize SliceGradient
% Here we create a Figure that shows the ASL slices/FoV with respect to standard space/
% in standard space. Usually these are rotated forwards, for the ACPC angulation, and to
% accommodate the forwards-rotated shape of the cerebrum

xASL_adm_CreateDir(x.D.SliceCheckDir);
x.S.TraSlices = [];
x.S.CorSlices = 78;
x.S.SagSlices = 59;

% Create brainmask
WB = logical(xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir, 'brainmask.nii')));

% Create SliceGradient mask with T1 template as background
SliceMask = xASL_io_Nifti2Im(x.P.Pop_Path_SliceGradient)~=0;
T1template = fullfile(x.D.MapsSPMmodifiedDir, 'rT1.nii');

IM1 = xASL_im_CreateVisualFig(x, T1template, [], [], [], [], [], WB, true);
IM2 = xASL_im_CreateVisualFig(x, x.P.Pop_Path_SliceGradient, [], [], [], x.S.jet256, [], SliceMask, true);
IM3 = xASL_im_CreateVisualFig(x, {T1template x.P.Pop_Path_SliceGradient}, [],[0.65 0.6],[],{x.S.gray,x.S.jet256},[],{WB SliceMask},true);

Parms.IM = [IM1;IM2;IM3];

xASL_adm_CreateDir(x.D.SliceCheckDir);
OutputFile = fullfile(x.D.SliceCheckDir,['SliceGradient' x.P.SubjectID '_' x.P.SessionID '.jpg']);
xASL_delete(OutputFile);
if sum(isfinite(Parms.IM(:)))>100
    xASL_imwrite(Parms.IM,OutputFile);
end

%% -----------------------------------------------------------------------------------
%  6) Temporal QC parameters, based on SPM Univariate+ Toolbox
%  Run this part only if we have > 10 time points
if nVolumes>10
    if xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
        % If segmented files exist, then the PVgm should already be prepared previously

        if xASL_exist(x.P.Path_PWI4D, 'file')
            ExistPWI4D = true;
        else
            ExistPWI4D = false;
            [ControlIM, LabelIM] = xASL_quant_GetControlLabelOrder(xASL_io_Nifti2Im(x.P.Path_ASL4D));
            xASL_io_SaveNifti(x.P.Path_ASL4D, x.P.Path_PWI4D, ControlIM-LabelIM);
        end

        TempASL = xASL_qc_temporalSNR(x.P.Path_PWI4D,{x.P.Path_PVgm x.P.Path_PVwm});
        if ~isempty(TempASL)
            ASLFields = fields(TempASL);
            for iA=1:length(ASLFields)
                x.Output.ASL.(['ASL_' ASLFields{iA}]) = TempASL.(ASLFields{iA});
            end
        end

        if ~ExistPWI4D
            xASL_delete(x.P.Path_PWI4D);
        end
    else
        warning('Skipping SPM UP QC because structural files didnt exist');
    end
end


%% 7) Compute overlap (intersection)/DICE coefficient with T1w brainmask
if xASL_exist(x.P.Path_mean_control,'file')
    % if mean control doesnt exist, changes are we have a 3D file, then
    % BET doesnt really work, and coverage is usually fine
    x.Output.ASL.ASL_CoveragePerc = xASL_qc_ComputeFoVCoverage(x.P.Path_mean_control, x);
end


%% 8) Summarize ASL orientation & check for left-right flips
xASL_qc_PrintOrientation(x.SESSIONDIR, x.P.Path_ASL4D, x.SESSIONDIR, 'RigidRegASL');


%% 9) Collect several other parameters & store all in PDF overview
x = xASL_qc_CollectParameters(x, x.iSubject, 'ASL'); % Quick & Dirty solution, 0 == skip structural part
xASL_delete(PathX);
save(PathX,'x'); % future: do this in each xWrapper
xASL_qc_CreatePDF(x, x.iSubject);

end
