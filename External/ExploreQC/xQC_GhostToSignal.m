function [ Func_GSR ] = xQC_GhostToSignal(Im4DPath, PathToTemplate)
% xQC_GhostToSignal extract the Ghost to Signal ratio in a 4D image

% FORMAT: [[ Func_GSR ] = xQC_GhostToSignal(Im4DPath, PathToTemplate)
%
% INPUT:
%   Im4DPath   	     - Path to the Functional Image, this needs to be a
%                          path
%  PathToTemplate    - Path to the Template with Ghost Rois. This need to
%                      have 1 for Signal Roi, 2 for NonGhost Rois, 3 for Ghost Rois
%  
% OUTPUT:
%   Func_GSR          - Ghost to Signal Ratio
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 
% 
% EXAMPLE:  Func_GSR  = xQC_GhostToSignal(Path/to/Bold.nii, Path/To/GSR_rois.nii)
%
%


SubjFold = fileparts(Im4DPath);
% Resampled Rois Mask into 4D image space
xASL_spm_deformation(Im4DPath, PathToTemplate, [], [], [], 'Tmp_ghost_resampled', [])
xASL_spm_deformations([], PathToTemplate, fullfile(SubjFold, 'Tmp_ghost_resampled.nii'), [], Im4DPath, [], fullfile(SubjFold, 'y_ASL.nii'), [])


%Load Images
GhostTemplate = xASL_io_Nifti2Im('Tmp_ghost_resampled.nii');
Im4D = xASL_io_Nifti2Im(Im4DPath);
xASL_delete('Tmp_ghost_resampled.nii' );

%Create Rois Mask
SignalMask = GhostTemplate == 1;
NonGhostMask = GhostTemplate == 2; 
GhostMask = GhostTemplate == 3; 

%Extract Signal from Masks
Signal = bsxfun(@times,Im4D,SignalMask); 
NonGhost = bsxfun(@times,Im4D,NonGhostMask);
Ghost = bsxfun(@times,Im4D,GhostMask);

%Compute Ghost to Signal 
Func_GSR = (xASL_stat_MeanNan(Ghost(:))- xASL_stat_MeanNan(NonGhost(:)) ) / xASL_stat_MeanNan(Signal(:));

end



