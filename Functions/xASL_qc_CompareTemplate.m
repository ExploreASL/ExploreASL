function [QC] = xASL_qc_CompareTemplate(x, ScanTypePrefix, iSubjectSession)
%xASL_qc_CompareTemplate QC Function that computes several QC parameters
%
% FORMAT: [QC] = xASL_qc_CompareTemplate(x, ModPrefix, iSubjectSession)
%
% INPUT:
%   x               - structure containing fields with all information required to run this submodule (REQUIRED)
%   ScanTypePrefix  - states which ScanType we are analyzing (OPTIONAL, DEFAULT=raw EPI)
%   iSubjectSession - Index of current subject/session
%
% OUTPUT:
%    QC           - QC parameters computed in this function:
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function computes several advanced template-based QC parameters:
%              RMSE_Perc        - Root Mean Square Error between image and template (%)
%              nRMSE_Perc       - Same but then normalized
%              AI_Perc          - Asymmetry Index between image and template (%)
%              Mean_SSIM_Perc   - mean structural similarity index -> xASL_stat_MeanSSIM.m
%              PeakSNR_Ratio    - peak signal-to-noise ratio -> xASL_stat_PSNR.m
%
% EXAMPLE: [QC] = xASL_qc_CompareTemplate(x, ModPrefix, iSubjectSession)
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


%% ---------------------------------------------------------
%  Admin
if nargin<2 || isempty(ScanTypePrefix)
    ScanTypePrefix = '';
end

if strcmp(ScanTypePrefix,'qCBF')
    TemplateIM = xASL_io_Nifti2Im(fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_CBF.nii'));
else % if strcmp(ModPrefix,'mean_control') -> compare with EPI image
    TemplateIM = xASL_io_Nifti2Im(fullfile(x.D.TemplateDir,'Philips_2DEPI_noBsup_Control.nii'));
end

if nargin<3 || isempty(iSubjectSession)
    fprintf('Creating QC parameters:   ');
    iSubject = 1:x.nSubjects;
    iSession = 1:x.nSessions;
    iSubjectSession = ((iSubject-1)*x.nSessions)+iSession;
else % when running for all scans
    [iSubject(1), iSession(1)] = xASL_adm_ConvertSubjSess2Subj_Sess(x.nSessions, iSubjectSession);
end

TemplateIM(~x.skull) = NaN;
TemplateColumnAll = xASL_im_IM2Column(TemplateIM,x.skull);
% PM: we can add different vendors templates, or use the 4D template from spatial CoV


%% ---------------------------------------------------------
%  Start iteration

for iImage=1:length(iSubjectSession)
    if iImage>1 % only verbose when iterating
        xASL_TrackProgress(iImage, length(iSubjectSession));
    end

%     QC.RMSE_Perc = NaN;
%     QC.AI_Perc = NaN;
%     QC.Mean_SSIM = NaN;
%     QC.PSNR = NaN;
%     QC.nRMSE_Perc = NaN;

    % Load image
    QC.ID      = [x.SUBJECTS{iSubject(iImage)} '_' x.SESSIONS{iSession(iImage)}];
    PathIM     = fullfile(x.D.PopDir,[ScanTypePrefix '_' QC.ID '.nii']);
    if xASL_exist(PathIM,'file')

        ASLim                   = xASL_io_Nifti2Im(PathIM);
        ASLim(~x.skull)         = NaN;

        ASLim                   = xASL_im_ndnanfilter(ASLim,'gauss',[8 8 8]./[1.5 1.5 1.5]);
        ASLcolumn               = xASL_im_IM2Column(ASLim,x.skull);

        % Mask out non-finite values
        FiniteMask              = isfinite(ASLcolumn) & isfinite(TemplateColumnAll);
        ASLcolumn               = ASLcolumn(FiniteMask);
        TemplateColumn          = TemplateColumnAll(FiniteMask);

        %% Calculate linear normalization to minimize RMS
        X = ASLcolumn(:);
        X = [ones(length(X),1), X];
        Solution = pinv(X)*TemplateColumn(:);
        ASLcolumnNorm = ASLcolumn*Solution(2)+Solution(1);
        TemplateColumnNorm = TemplateColumn;

        DiffColumn = ASLcolumn-TemplateColumn;
        DiffColumnNorm = ASLcolumnNorm-TemplateColumnNorm;

        MeanASL = xASL_stat_MeanNan(ASLcolumn) + xASL_stat_MeanNan(TemplateColumn) /2;
        MeanASLNorm = xASL_stat_MeanNan(ASLcolumnNorm) + xASL_stat_MeanNan(TemplateColumnNorm) /2;

        QC.RMSE_Perc = xASL_round(100 * xASL_stat_MeanNan(DiffColumn.^2).^0.5 / MeanASL,4); % RootMeanSquare (%)
        QC.nRMSE_Perc = xASL_round(100 * xASL_stat_MeanNan(DiffColumnNorm.^2).^0.5 / MeanASLNorm,4); % normalized RootMeanSquare (%)

        QC.Mean_SSIM_Perc = xASL_round(xASL_stat_MeanSSIM(TemplateIM,ASLim),3);
        QC.PeakSNR_Ratio = xASL_round(100*xASL_stat_PSNR(TemplateIM,ASLim),3);

        %% Calculate asymmetry Index (AI)
        
        AI_Perc = xASL_qc_AsymmetryIndex(ASLim);
        QC.AI_Perc = AI_Perc;
    end
end



fprintf('\n');


end








