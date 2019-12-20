function [FA_Outliers_mL] = xASL_qc_FA_Outliers(InputFA)
%xASL_qc_FA_Outliers extract the number of FA outliers, i.e. values of FA
%above 1 or below 0, from a FA image. 
%
% FORMAT: [ FA_Outliers_mL ] = xASL_qc_FA_Outliers(InputFA)
%
% INPUT:
%   InputFA 	     - This can either be the path to the FA.nii file or an
%                      image itself
%  
%
% OUTPUT:
%   FA_Outliers_mL   - Number of Outlier expressed in mL
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% 
% 
% EXAMPLE: [FA_Outliers_mL] = xASL_qc_FA_Outliers('/examplepath/subjectdir/FA.nii');
%
%

% Read the header and get voxelsize
FAnii = xASL_io_ReadNifti(InputFA);
VoxelSize_mm = prod(FAnii.hdr.pixdim(2:4));

%Read image
FAimage = xASL_io_Nifti2Im(InputFA);

%Detect FA outliers
nFA_OutliersPos = xASL_stat_SumNan(FAimage(:)>1); % get number of positive FA outliers
nFA_OutliersNeg = xASL_stat_SumNan(FAimage(:)<0); % get number of negative FA outliers
nFA_Outliers = nFA_OutliersPos+nFA_OutliersNeg; % get total number of outliers

FA_Outliers_mm = nFA_Outliers*VoxelSize_mm;
FA_Outliers_mL = FA_Outliers_mm/1000;


end 
