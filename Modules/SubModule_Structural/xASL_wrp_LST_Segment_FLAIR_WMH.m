function xASL_wrp_LST_Segment_FLAIR_WMH(x, rWMHPath, WMHsegmAlg)
%xASL_wrp_LST_Segment_FLAIR_WMH Submodule of ExploreASL Structural Module, that performs a biasfield correction on T1w
% & applies it on the FLAIR
%
% FORMAT: xASL_wrp_LST_Segment_FLAIR_WMH(x, rWMHPath[, WMHsegmAlg])
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   rWMHPath - path of the WMH segmentation (either performed by LST or a copy of x.P.Path_WMH_SEGM) (REQUIRED)
%   WMHsegmAlg - Choose the LST algorithm 'LGA'/'LPA' (OPTIONAL, DEFAULT = 'LPA')

% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule runs the LST WMH segmentation, either with LGA or LPA.
% LPA is the default choice, it outperforms LGA a bit, depending on the image quality. These algorithms
% perform optimally with 3T images, with good contrast. Generally, LPA oversegments whereas LGA undersegments.
% The LPA oversegmentation is corrected in a later submodule.
%
% EXAMPLE: xASL_wrp_LST_Segment_FLAIR_WMH(x, rWMHPath);
%
% REFERENCE:
% de Sitter A, Steenwijk MD, Ruet A, et al. Performance of five research-domain automated WM lesion segmentation methods in a multi-center MS study. Neuroimage. 2017;163(August):106?114.
% Schmidt P, Gaser C, Arsic M, et al. An automated tool for detection of FLAIR-hyperintense white-matter lesions in Multiple Sclerosis. Neuroimage. 2012;59(4):3774?3783.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-05-02 HJM


%% ADMIN

if nargin < 2 || isempty(rWMHPath)
	error('xASL_wrp_LST_Segment_FLAIR_WMH: Requires at least 2 input arguments.');
end

if nargin < 3 || isempty(WMHsegmAlg)
	WMHsegmAlg = 'LPA';
end



%% ----------------------------------------------------------
%% 1) Reslice FLAIR (& WMH_SEGM, if exists) to T1w
xASL_io_SaveNifti(x.P.Path_FLAIR, x.P.Path_FLAIR, xASL_io_Nifti2Im(x.P.Path_FLAIR), 16, false); % convert to 16 bit
xASL_spm_reslice(x.P.Path_T1, x.P.Path_FLAIR, [], [], x.Quality);

if xASL_exist(x.P.Path_WMH_SEGM, 'file')
    % first backup externally provided WMH_SEGM
    if ~xASL_exist(x.P.Path_WMH_SEGM_ORI, 'file')
        xASL_Copy(x.P.Path_WMH_SEGM, x.P.Path_WMH_SEGM_ORI);
    end
    
    % if an externally provided WMH_SEGM exists, resample it to the T1w space
    % use linear resampling, to avoid B-spline edge effects
    xASL_spm_reslice(x.P.Path_T1, x.P.Path_WMH_SEGM, [], [], x.Quality, x.P.Path_WMH_SEGM, 1);
end



%% ----------------------------------------------------------
%% 2) Define parameters for segmentation

fprintf('\n%s\n','----------------------------------------');

switch WMHsegmAlg
    case 'LPA'
        fprintf('%s\n','WMH segmentation performed using LST LPA');

        % Lesion Prediction Algorithm (LPA)
        matlabbatch{1}.spm.tools.LST.lpa.data_F2        = {x.P.Path_rFLAIR}; % FLAIR resampled to T1 native space
        matlabbatch{1}.spm.tools.LST.lpa.data_coreg     = {''};
        matlabbatch{1}.spm.tools.LST.lpa.html_report    = 0; % no HTML report output (takes a long time)
        
        if xASL_exist(x.P.Path_WMH_SEGM, 'file')
            matlabbatch{1}.spm.tools.LST.lpa.xasl_quality   = 2; % ultralow quality
        else
            matlabbatch{1}.spm.tools.LST.lpa.xasl_quality   = x.Quality;
        end

    case 'LGA'
        fprintf('%s\n','WMH segmentation performed using LST LGA');

        % Lesion Growth Algorithm (LGA)
        matlabbatch{1}.spm.tools.LST.lga.data_T1            = {x.P.Path_T1};
        matlabbatch{1}.spm.tools.LST.lga.data_F2            = {x.P.Path_rFLAIR};
        InitialParm                                         = 0.3; % default
        matlabbatch{1}.spm.tools.LST.lga.opts_lga.initial   = InitialParm;
        matlabbatch{1}.spm.tools.LST.lga.opts_lga.mrf       = 1; % default
        matlabbatch{1}.spm.tools.LST.lga.html_report        = 0; % no HTML report generation, takes too long
        matlabbatch{1}.spm.tools.LST.lga.xasl_quality       = x.Quality;

        if x.Quality
            matlabbatch{1}.spm.tools.LST.lga.opts_lga.maxiter = 100; % default=50
        else
            matlabbatch{1}.spm.tools.LST.lga.opts_lga.maxiter = 3;
		end
	otherwise
		error('Unknown or undefined segmentation algorithm.');
end



%% ----------------------------------------------------------
%% 3 Run the segmentation
fprintf('\n');
spm_jobman('run',matlabbatch); close all


%% ----------------------------------------------------------
%% 4) File management
xASL_io_SaveNifti(rWMHPath, rWMHPath, xASL_io_Nifti2Im(rWMHPath), 16, false); % convert to 16 bit

% Create a copy of the WMH segmentation, if no externally provided WMH_SEGM
if ~xASL_exist(x.P.Path_WMH_SEGM, 'file')
    xASL_Copy(rWMHPath, x.P.Path_WMH_SEGM);
end

% Delete the LST folder & its contents
CleanUpDir = xASL_adm_GetFileList(x.SUBJECTDIR,'^LST_tmp_.*$','List',[0 Inf], true);
nList = length(CleanUpDir);
for iD=1:nList
    xASL_adm_DeleteFileList(fullfile(x.SUBJECTDIR,CleanUpDir{iD}), '^.*\.(nii|nii\.gz)$', [], [0 Inf]);
    xASL_adm_DeleteFileList(fullfile(x.SUBJECTDIR,CleanUpDir{iD}), '^.*\.mat$', [], [0 Inf]);
    rmdir( fullfile(x.SUBJECTDIR,CleanUpDir{iD}),'s' );
end

if x.DELETETEMP
    xASL_delete(x.P.Path_mrFLAIR); % LPA
    xASL_delete(x.P.Path_rmrFLAIR); % LGA

    if ~x.bReproTesting
        xASL_delete(x.P.Path_rFLAIR);
    end
end

%% Remove NaNs from segmentations & fix image edges
if xASL_exist(x.P.Path_WMH_SEGM, 'file')
    tIM = xASL_io_Nifti2Im(x.P.Path_WMH_SEGM);
    tIM(isnan(tIM)) = 0; % set NaNs to zeros

    tIM = xASL_im_ZeroEdges(tIM); % set edges to zero
    xASL_io_SaveNifti(x.P.Path_WMH_SEGM, x.P.Path_WMH_SEGM, tIM, [], false);
end
if xASL_exist(rWMHPath, 'file')
    tIM = xASL_io_Nifti2Im(rWMHPath);
    tIM(isnan(tIM)) = 0; % set NaNs to zeros

    tIM = xASL_im_ZeroEdges(tIM); % set edges to zero
    xASL_io_SaveNifti(rWMHPath, rWMHPath, tIM, [], false);
end




end
