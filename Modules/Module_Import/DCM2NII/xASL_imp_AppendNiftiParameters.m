function s = xASL_imp_AppendNiftiParameters(nii_files)
%xASL_imp_AppendNiftiParameters Append Nifti Parameters.
%
% FORMAT: s = xASL_imp_AppendNiftiParameters(nii_files)
% 
% INPUT:
%   nii_files  - List of NIFTI files (REQUIRED, CELL ARRAY)
%
% OUTPUT:
%   s          - String
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Append Nifti Parameters.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     s = xASL_imp_AppendNiftiParameters(nii_files);
% __________________________________
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


    %% Append Nifti Parameters

    % This function outputs s = [FileName voxel size XYZ matrix size XYZ]
    s = [];

    if ischar(nii_files)
        nii_files = {nii_files};
    end

    for iNii=1:length(nii_files)
        [~, Ffile, Fext] = fileparts(nii_files{iNii});
        s = sprintf(',"%s"', [Ffile Fext]); % filename

        tempnii = xASL_io_ReadNifti(nii_files{iNii});
        s = [s sprintf(',%g', tempnii.hdr.pixdim(2:5) )]; % voxel size XYZ
        s = [s sprintf(',%g', tempnii.hdr.dim(2:5) )]; % matrix size XYZ
    end
    
end



