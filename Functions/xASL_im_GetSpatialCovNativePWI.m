function [spatCoV] = xASL_im_GetSpatialCovNativePWI(x)
%xASL_im_GetSpatialCovNativePWI Acquires spatial CoV from the native space ASL
%image, using registered mask

PathMaskTemplate = fullfile(x.SESSIONDIR,'Mask_Template.nii');
if ~xASL_exist(x.PathMask,'file')
    xASL_Copy(PathMaskTemplate, x.PathMask);
end
PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);

xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, PathMaskTemplate, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, x.PathMask,0);
MaskIM = xASL_io_Nifti2Im(x.PathMask);
   
% Ensure that the mask is binary
MaskIM = MaskIM > 0.5;

%% Simple intersection check
sortInt = sort(PWIim(:));
ThrInt = sortInt(round(0.75*length(sortInt)));
maskGM = PWIim>ThrInt;
maskBig = maskGM | MaskIM;
maskGM = maskGM & MaskIM;
JointMasks = sum(maskGM(:))/sum(maskBig(:));
fprintf('%s\n',['Joint brainmask for automatic spatial CoV detection was ' num2str(100*JointMasks,3) '%']);

if JointMasks<0.25
    warning('Registration off, spatial CoV detection unreliable');
end

%% Determine spatial CoV
spatCoV = xASL_stat_ComputeSpatialCoV(PWIim, MaskIM, 0);

spatCoV = spatCoV/1.5; % correction native space including WM to MNI spatial CoV, excluding WM

if spatCoV<0
    warning('Native space whole-brain spatial CoV was negative! (i.e. <0)');
    fprintf('Defaulting to spatial CoV of 40%\n');
    spatCoV = 0.4;
end

fprintf('%s\n',['Standard space whole-brain spatial CoV estimated as = ' num2str(100*spatCoV,3) '%']);

end
