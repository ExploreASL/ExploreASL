function xASL_bids_MergeNifti_Delete(NiftiPaths)
%xASL_bids_MergeNifti_Delete Delete NiftiPaths and associated JSONs
%
% FORMAT: xASL_bids_MergeNifti_Delete(NiftiPaths);
% 
% INPUT:
%   NiftiPaths - Nifti paths
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Delete NiftiPaths and associated JSONs.
%
% EXAMPLE:     ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Delete NiftiPaths and associated JSONs
    for iFile=1:length(NiftiPaths)
        [Fpath, Ffile] = xASL_fileparts(NiftiPaths{iFile});

        xASL_delete(NiftiPaths{iFile});

        pathJSON = fullfile(Fpath,[Ffile '.json']);
        % Delete JSONs
        xASL_delete(pathJSON); % already checks if exists before deleting
    end

end
