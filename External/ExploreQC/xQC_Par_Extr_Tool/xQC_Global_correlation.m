function [GlobalCorr] = xQC_Global_correlation(InputBoldPath, InputGMPath, SubjDir)
%xASL_qc_Global_correlation extract Global Correlation, i.e. the mean
%correlation between each pairs of Gray Mattter voxels timeseries, from a functional MRI
%image
%
% FORMAT: [GlobalCorr] = xQC_Global_correlation(InputBoldPath, InputGMPath, SubjDir)
%
% INPUT:
%   InputBoldPath 	     - Path to the Functional Image, this needs to be a
%                          path
%   InputGMPath          - Path to the GM segmentation image, this needs to
%                          be a path
%   SubjDir              - Path to the subject directory, this is needed to
%                          store temporary files
%
% OUTPUT:
%   GlobalCorr           - Global Correlation coefficient
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 
% 
% EXAMPLE: GlobalCorr = xQC_Global_correlation('/examplepath/exampleSubj/func.nii, /examplepath/exampleSubj/c1T1.nii, /examplepath/exampleSubj/)
%
%

%Create Image from Path
BoldIm = xASL_io_Nifti2Im(InputBoldPath);
GMim = xASL_io_Nifti2Im(InputGMPath);       

% reslice GM into fMRI space
xASL_spm_reslice(InputBoldPath, InputGMPath, [], [], [], [SubjDir 'tmp_GM_mask.nii'], 0)

% load GM image in lower func resolution
GMim = xASL_io_Nifti2Im([SubjDir 'tmp_GM_mask.nii']);       
xASL_delete([SubjDir 'tmp_GM_mask.nii'])

% 1 mask pGM higher than 0.5 of the maximal pGM
GMim = GMim>0.5;

% make mask of functional infinite values (3D)
BoldImMask = isfinite(BoldIm);
BoldImMask = prod(BoldImMask, 4);

% 3 combine masks by multiplying GM mask with isfinite mask in a logical
% (==1) operator
Mask = ((GMim.*BoldImMask) ==1);

% Mask The bold Image, Reshape it into 2 dimensions (voxels TimeSeries),
% clean the rows that are zeros
BoldMasked = bsxfun(@times,BoldIm,Mask);     
Voxels_TS = reshape(BoldMasked, [size(BoldMasked,1)*size(BoldMasked,2)*size(BoldMasked,3), size(BoldMasked,4)]);
Voxels_TS_clean = Voxels_TS(any(Voxels_TS,2),:);     %


% perform the corrcoef, to obtain correlation matrix
% mean but take into account potential nans
CC_all = corrcoef(Voxels_TS_clean');
CC_vox_mean = mean(CC_all);
GlobalCorr = mean(CC_vox_mean);

end
