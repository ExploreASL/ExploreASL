function [SSE_WM_Mean] = xQC_SSE_WM_mean(SSEpath, WMpath)

% Reslice the WM mask into the SSE space
[ddir, dname, dext]=fileparts(SSEpath)
xASL_spm_reslice(SSEpath, WMpath, [], [], [], fullfile(ddir, 'tmp_resampled_GM.nii'), [])

%Read SSE and WM
SSE = xASL_io_Nifti2Im(SSEpath);
WM = xASL_io_Nifti2Im(fullfile(ddir, 'tmp_resampled_GM.nii'));
WM = WM>0.5;


SSE_WM = SSE(WM);

SSE_WM_Mean = xASL_stat_MeanNan(SSE_WM(:));

end