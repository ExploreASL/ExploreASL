function [AI_perc] = xQC_AsymmetryIndex(ImageIn)
% xASL_qc_Asymmetry Extract voxel-wise asymmetry index for QC purposes
% 
%  FORMAT: [AI_perc] = xASL_qc_AsymmetryIndex(ImageIn)
%
% INPUT:
%   ImageIn 	     - Inut Image Matrix, an also be a path to a NIfTI image
%
% OUTPUT:
%   AI_perc          - Average Percentage Voxel-Wise asymmetry
% 
% EXAMPLE: AI_perc = xASL_qc_AsymmetryIndex('/data/StudyName/SubjectName/T1.nii');
%  Luigi Lorenzini
% ----------------------------------------------------------------------------------------------------------------------------------------------------

% if it's a path, make an image
ImageIn = xASL_io_Nifti2Im(ImageIn);

%flip it
NiiMirror = ImageIn(end:-1:1,:,:);  

%extract Asymmetry Index
Asymm = ImageIn - NiiMirror;
AI = Asymm ./ (0.5.*(ImageIn + NiiMirror));
AI_perc = xASL_round(100 * xASL_stat_MedianNan(abs(AI(:))),4); % convert to percentages

end 



