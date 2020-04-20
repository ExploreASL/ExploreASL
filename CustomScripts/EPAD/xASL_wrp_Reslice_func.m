function xASL_wrp_Reslice_func(x)
% xASL_wrp_Reslice_func
%
%
% 1    Registration ASL -> T1 (& M0 if there is a separate M0 or mean_control image without background suppression)
% 2    Create slice gradient image for quantification reference, in case of 2D ASL
% 3    Reslice ASL time series to MNI space (currently 1.5 mm^3)
% 4    Create mean control image, masking to 20% of max value if used as M0 (no background suppression)
% 5    Smart smoothing mean_control if used as M0
%


%% ------------------------------------------------------------------------------------------
%% Administration
tempnii             = xASL_io_ReadNifti(x.P.Path_func_bold);
% min_voxelsize       = double(min(tempnii.hdr.pixdim(2:4) )); % repmat, 1,3
nVol                = double(tempnii.hdr.dim(5));

if ~isfield(x,'SavePWI4D')
    x.SavePWI4D   = 0;
end












%% ------------------------------------------------------------------------------------------
%% 7    Compute time-series-related images & save maps

tPWI                                = xASL_io_Nifti2Im(x.P.Path_func_bold);
PWI                                 = xASL_stat_MeanNan(tPWI,4);

% SortInt                             = sort(PWI(isfinite(PWI)));
% Thresh                              = SortInt(round(0.75*length(SortInt)));
% MaskIM                              = PWI>Thresh;

% std0                                = std(PWI(MaskIM & isfinite(PWI)));
% mean0                               = mean(PWI(MaskIM & isfinite(PWI)));
% CoV0                                = (std0/mean0)*100;

if size(tPWI,4)>1

    new_SD                          = xASL_stat_StdNan( tPWI,0,4);
    new_SNR                         = PWI./new_SD;
    new_SNR(new_SNR<0)              = 0; % clip @ zero

    xASL_io_SaveNifti(x.P.Path_func_bold, x.P.Path_SD , new_SD , 32, 0);
    xASL_io_SaveNifti(x.P.Path_func_bold, x.P.Path_SNR, new_SNR, 32, 0);
    fprintf('%s\n','Standard space timeseries-related images saved (SD & SNR), & part1 & part2 for reproducibility');
else
    fprintf('%s\n',['Standard space SD & SNR maps were not created because of only ' num2str(size(tPWI,4)) ' frame(s)']);
end

xASL_io_SaveNifti(x.P.Path_func_bold, x.P.Path_mean_control, PWI, 32, false);

% Deform
InputPaths = {x.P.Path_mean_control x.P.Path_SD x.P.Path_SNR};
OutputPaths = {x.P.Pop_Path_mean_control x.P.Pop_Path_SD x.P.Pop_Path_SNR};

if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % BACKWARDS COMPATIBILITY, CAN BE REMOVED
	AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
else
	AffineTransfPath = [];
end

xASL_spm_deformations(x, InputPaths, OutputPaths, 2, [], AffineTransfPath, x.P.Path_y_ASL);

fprintf('%s\n','Standard space func maps saved');

end
