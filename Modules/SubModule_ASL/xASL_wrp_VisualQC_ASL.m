function xASL_wrp_VisualQC_ASL(x)
%xASL_wrp_VisualQC_ASL Submodule of ExploreASL ASL Module, that performs several visualizations for QC
%
% FORMAT: xASL_wrp_VisualQC_ASL(x)
%
% INPUT:
%   x                 - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P               - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   x.MakeNIfTI4DICOM - Boolean, true for resampling CBF native space to
%                       original T1w & ASL spaces, and other processing for use in
%                       DICOM image/server (OPTIONAL, DEFAULT=false);
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule performs several visualizations for visual & quantitative QC.
%
%              1. After initial admin
%              2. It starts with making ASL NIfTIs ready for visualization
%                 & conversion to DICOM (though skipped by default)
%              3. Then it performs a collection of visualizations
%              4. Visualizes results of the TopUp geometric distortion correction
%              5. Visualization of slice gradient
%              6. Visualization & calculation of temporal QC parameters
%              7. Compute DICE overlap/intersection of ASL brain in FoV & T1w, to calculate coverage
%              8. Summarize orientation & check left-right flips
%              9. Collect several other parameters & store in PDF overview
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

x = xASL_adm_LoadX(x, [], true); % assume x.mat is newer than x

%% -----------------------------------------------------------------------------------
%% 2) Make ASL NIfTIs ready for visualization & conversion to DICOM
if isfield(x,'MakeNIfTI4DICOM') && x.MakeNIfTI4DICOM
    if xASL_exist(x.P.Path_T1_ORI, 'file')
        InputT1Oripath = x.P.Path_T1_ORI;
    else
        InputT1Oripath = x.P.Path_T1;
    end

    if xASL_exist(x.P.Path_ASL4D_ORI, 'file')
        InputASLOripath = x.P.Path_ASL4D_ORI;
    else
        InputASLOripath = x.P.Path_ASL4D;
    end

    xASL_io_MakeNifti4DICOM(x.P.Path_CBF, x, 'UINT16', x.P.Path_T1, InputT1Oripath);
    xASL_io_MakeNifti4DICOM(x.P.Path_CBF, x, 'UINT16', x.P.Path_ASL4D, InputASLOripath); 
end

%% -----------------------------------------------------------------------------------
%% 3) Perform several visualizations
x = xASL_qc_VisualCheckCollective_ASL(x);

%% -----------------------------------------------------------------------------------
%% 4) Visualization TopUp results (quick & dirty)
PathPopB0 = fullfile(x.D.PopDir, ['rASL_B0_' x.P.SubjectID '_' x.P.SessionID '.nii']);
PathPopUnwarped = fullfile(x.D.PopDir, ['rASL_Unwarped_' x.P.SubjectID '_' x.P.SessionID '.nii']);

if xASL_exist(PathPopB0,'file') && xASL_exist(PathPopUnwarped,'file')% if we have TopUp results
    [Output1, Output2] = xASL_vis_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, x.iSubject, x.D.ASLCheckDir);
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

IM1 = xASL_vis_CreateVisualFig(x, T1template, [], [], [], [], [], WB, true);
IM2 = xASL_vis_CreateVisualFig(x, x.P.Pop_Path_SliceGradient, [], [], [], x.S.jet256, [], SliceMask, true);
IM3 = xASL_vis_CreateVisualFig(x, {T1template x.P.Pop_Path_SliceGradient}, [],[0.65 0.6],[],{x.S.gray,x.S.jet256},[],{WB SliceMask},true);

Parms.IM = [IM1;IM2;IM3];

xASL_adm_CreateDir(x.D.SliceCheckDir);
OutputFile = fullfile(x.D.SliceCheckDir,['SliceGradient' x.P.SubjectID '_' x.P.SessionID '.jpg']);
xASL_delete(OutputFile);
if sum(isfinite(Parms.IM(:)))>100
    xASL_vis_Imwrite(Parms.IM,OutputFile);
end

%% -----------------------------------------------------------------------------------
%  6) Temporal QC parameters, based on SPM Univariate+ Toolbox
%  Run this part only if we have > 10 time points
if nVolumes>10
    if xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
        % If segmented files exist, then the PVgm should already be prepared previously

        if ~xASL_exist(x.P.Path_PWI4D, 'file')
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

    else
        warning('Skipping SPM UP QC because structural files didnt exist');
    end
end

xASL_delete(x.P.Path_PWI4D);

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

xASL_adm_SaveX(x); % future: do this in each xWrapper
xASL_qc_CreatePDF(x, x.iSubject);

end




%% *****************************************************************************************************************

function x = xASL_qc_VisualCheckCollective_ASL(x)
%xASL_qc_VisualCheckCollective_ASL Runs a collection of visual QC functions


%% Get visualization settings
% Parameters for creating visual QC Figures:
% CBF, CBF with overlay c2T1, CBF with overlay c2T1
% MeanControl SD
% M0 NoSmoothM0 NoSmoothM0 with overlay c1T1
% TT TT with overlay c2T1

x               = xASL_adm_ResetVisualizationSlices(x);

T.ImIn          = {x.P.Pop_Path_qCBF  x.P.Pop_Path_SD {x.P.Pop_Path_qCBF x.P.Pop_Path_PV_pWM} x.P.Pop_Path_SNR};
T.ImIn(5:8)     = {x.P.Pop_Path_mean_control x.P.Pop_Path_noSmooth_M0 {x.P.Pop_Path_noSmooth_M0 x.P.Pop_Path_PV_pGM} x.P.Pop_Path_M0};
T.ImIn(9:10)    = {x.P.Pop_Path_TT  {x.P.Pop_Path_TT x.P.Pop_Path_PV_pWM}};

T.DirOut        = {x.D.ASLCheckDir x.D.SNRdir      x.D.ASLCheckDir       x.D.SNRdir};
T.DirOut(5:8)   = {x.D.RawDir      x.D.M0CheckDir  x.D.M0regASLdir       x.D.M0CheckDir};
T.DirOut(9:11)  = {x.D.TTCheckDir  x.D.TTCheckDir  x.D.ASLCheckDir};

T.IntScale(2)   = {[1 1]};
T.IntScale{8}   = [0.75 0.65];

T.ColorMapIs{10}= x.S.jet256;
T.ColorMapIs{11}= {x.S.jet256};

T.NameExt       = {[] [] 'Reg_pWM' []};
T.NameExt( 5:8) = {[] [] 'Reg_pGM' []};
T.NameExt(9:11) = {[] 'Reg_pWM' 'AnalysisMask'};

% Fill missing cells
Pars = {'ImIn' 'DirOut' 'ClipZero' 'IntScale' 'NameExt' 'ColorMapIs'}; % default pars
for iM=1:length(T.ImIn)
    for iP=1:length(Pars)
        if ~isfield(T,Pars{iP})
            T.(Pars{iP}) = [];
        elseif length(T.(Pars{iP}))<iM
            T.(Pars{iP}){iM} = [];
        end
    end
end

%% Take only NIfTIs that exist
for iL=1:length(T.ImIn)
    if  ischar(T.ImIn{iL})
        ExistInd(iL) = logical(xASL_exist(T.ImIn{iL},'file'));
    else
        ExistInd(iL) = min(cellfun(@(x) logical(xASL_exist(char(x),'file')),T.ImIn{iL}));
    end
end
        
for iP=1:length(Pars)
    T.(Pars{iP}) = T.(Pars{iP})(ExistInd);
end

%% Clone each row into a transversal & coronal row
nIms = length(T.(Pars{1}));
nRows = ceil( nIms/4);

for iN=1:nRows
    clear T2
    ImsI                    = (iN-1)*4+1:min(nIms,iN*4);
    nImsRow                 = length(ImsI);
    nRow1                   = 1:nImsRow;
    nRow2                   = nImsRow+1:2*nImsRow;
    T2.TraSlices(nRow1)     = {[]};
    T2.TraSlices(nRow2)     = {'n/a'};
    T2.CorSlices(nRow1)     = {'n/a'};
    T2.CorSlices(nRow2)     = {[]};

    T2.NameExt(nRow1)       = cellfun(@(x) ['Tra_' x], T.NameExt(ImsI), 'UniformOutput',false);
    T2.NameExt(nRow2)       = cellfun(@(x) ['Cor_' x], T.NameExt(ImsI), 'UniformOutput',false);
    T2.ImIn                 = [T.ImIn(ImsI) T.ImIn(ImsI)];
    T2.DirOut               = [T.DirOut(ImsI) T.DirOut(ImsI)];
    T2.IntScale             = [T.IntScale(ImsI) T.IntScale(ImsI)];
    T2.ColorMapIs           = [T.ColorMapIs(ImsI) T.ColorMapIs(ImsI)];
    T2.ModuleName           = 'ASL';

%%  Perform the visualization
% Perhaps at the end of the row we need to generate empty images, as transversal & coronal have different sizes
% they don't concatenate well horizontally, need to be concatenated vertically

    fprintf('%s','Printing images...  ');
    for iM=1:length(T2.ImIn)
        xASL_TrackProgress(iM,length(T2.ImIn)*nRows);

        % Manage slices to show
        % Sagittal
        x.S.SagSlices   = []; % show no sagittal slices
        % Transversal
        if      isempty(T2.TraSlices{iM})
                x.S.TraSlices   = x.S.slicesLarge;
        elseif  strcmp(T2.TraSlices{iM},'n/a')
                x.S.TraSlices   = [];
        else
                warning('Wrong slice choice');
        end
        % Coronal
        if      isempty(T2.CorSlices{iM})
                x.S.CorSlices   = x.S.slicesLarge+7;
        elseif  strcmp(T2.CorSlices{iM},'n/a')
                x.S.CorSlices   = [];
        else
                warning('Wrong slice choice');
        end    

        % Create the image
        T2.IM = xASL_vis_CreateVisualFig( x, T2.ImIn{iM}, T2.DirOut{iM}, T2.IntScale{iM}, T2.NameExt{iM}, T2.ColorMapIs{iM});
        % add single slice to QC collection
        if sum(~isnan(T2.IM(:)))>0 % if image is not empty
            x = xASL_vis_AddIM2QC(x,T2);
        end
    end
end

fprintf('\n');
   
end