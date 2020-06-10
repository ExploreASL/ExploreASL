function xASL_qc_Default_Test_Unzip(x)
%xASL_qc_Default_Test_Unzip Script to unzip Nifti files for test runs.
%
% FORMAT:       x = xASL_qc_Default_Test_Unzip(x);
% 
% INPUT:        x structure
%
% OUTPUT:       Files
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Unzip T1 and FLAIR Nifti files for test runs.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES:     x = xASL_qc_Default_Test_Unzip(x);
% __________________________________
% Copyright 2015-2020 ExploreASL

% Unzip Nifti files (T1)
IM = xASL_io_Nifti2Im(x.P.Path_T1);
if size(IM,4)>1 || size(IM,5)>1 || size(IM,6)>1 || size(IM,7)>1
    warning(['Too many dims, using first: ' x.P.Path_T1]);
    xASL_Copy(x.P.Path_T1, x.P.Path_T1_ORI, true);
    xASL_io_SaveNifti(x.P.Path_T1, x.P.Path_T1, IM(:,:,:,1,1,1,1), [], false);
end

% Unzip Nifti files (FLAIR)
IM = xASL_io_Nifti2Im(x.P.Path_FLAIR);
if size(IM,4)>1 || size(IM,5)>1 || size(IM,6)>1 || size(IM,7)>1
    warning(['Too many dims, using first: ' x.P.Path_T1]);
    xASL_Move(x.P.Path_FLAIR, x.P.Path_FLAIR_ORI, true);
    xASL_io_SaveNifti(x.P.Path_FLAIR, x.P.Path_FLAIR, IM(:,:,:,1,1,1,1), [], false);
end







