function [spatCoV] = xASL_im_GetSpatialCovNativePWI(x)
%xASL_im_GetSpatialCovNativePWI Acquires spatial CoV from the native space ASL
%image, using registered mask

PathMask = fullfile(x.SESSIONDIR,'Mask_Template.nii');
PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);
MaskIM = xASL_io_Nifti2Im(PathMask);

if min(size(PWIim)~=size(MaskIM))
    xASL_spm_reslice(x.P.Path_mean_PWI_Clipped, Mask_Native, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality, Mask_Native,0);
    MaskIM = xASL_io_Nifti2Im(PathMask);
end
    
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

if  JointMasks<0.25
    warning('Registration off, spatial CoV detection unreliable');
end

%% Determine spatial CoV
spatCoV = xASL_stat_ComputeSpatialCoV(PWIim, MaskIM, 0);

spatCoV = spatCoV/1.5; % correction native space including WM to MNI spatial CoV, excluding WM

if spatCoV<0
    error('Native space whole-brain spatial CoV was negative! (i.e. <0)');
end

fprintf('%s\n',['Standard space whole-brain spatial CoV estimated as = ' num2str(100*spatCoV,3) '%']);

end
