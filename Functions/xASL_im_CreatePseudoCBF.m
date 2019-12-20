%% ----------------------------------------------------------------------------------------
%% ----------------------------------------------------------------------------------------
%% ----------------------------------------------------------------------------------------
function xASL_im_CreatePseudoCBF(x, spatCoV)
%xASL_im_CreatePseudoCBF This function creates a pseudo-CBF image from mean CBF template,
% arterial transit time (ATT) bias field & vascular artifacts, weighted through spatial CoV

% Check, spatCoV should be a single number
if  numel(spatCoV)~=1
    error('Wrong definition of spatial CoV, this should be a single number');
elseif spatCoV<0
    error('Native space whole-brain spatial CoV was negative! (i.e. <0)');
end

Mean_Native         = fullfile(x.SESSIONDIR,'Mean_CBF_Template.nii');
Vasc_Native         = fullfile(x.SESSIONDIR,'VascularArtifact_Template.nii');
Bias_Native         = fullfile(x.SESSIONDIR,'ATT_BiasField.nii');

Mean_IM             = xASL_io_Nifti2Im(Mean_Native);
Bias_IM             = xASL_io_Nifti2Im(Bias_Native);
PWIim               = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);

PseudoCBF_Native    = fullfile(x.SESSIONDIR,'PseudoCBF.nii');

%% PseudoCBF equation. For 3D spiral, we nearly see no vascular artifacts because of smoothing
if  ~strcmp(x.Sequence,'3D_spiral')
    Vasc_IM             = xASL_io_Nifti2Im(Vasc_Native);
    PseudoCBFim         = Mean_IM - Bias_IM.*5.*spatCoV.^2 + Vasc_IM./5.*spatCoV.^2;
else
    PseudoCBFim         = Mean_IM - Bias_IM.*5.*spatCoV.^2 ;
end

PseudoCBFim(PseudoCBFim<0)  = 0; % clip at 0

xASL_io_SaveNifti(Mean_Native, PseudoCBF_Native, PseudoCBFim, [], 0);

%% scale mean_PWI_Clipped to same values as PseudoCBF
pseudoIM            = xASL_io_Nifti2Im(PseudoCBF_Native);

minPWI              = min(PWIim(:));
maxPWI              = max(PWIim(:));
minPseu             = min(pseudoIM(:));
maxPseu             = max(pseudoIM(:));
DiffMin             = minPWI - minPseu;
PWIim               = PWIim-DiffMin;
RatioMax            = maxPseu/maxPWI;
PWIim               = PWIim.*RatioMax;

xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, PWIim, [], 0);
fprintf('%s\n',['PseudoCBF.nii created for spatial CoV=' num2str(100*spatCoV,3) '% & rescaled mean_PWI_Clipped.nii to this']);

end
