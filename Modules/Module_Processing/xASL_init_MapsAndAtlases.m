function x = xASL_init_MapsAndAtlases(x)
%xASL_init_MapsAndAtlases Add maps and atlases
%
% FORMAT: x = xASL_init_MapsAndAtlases(x)
%
% INPUT:
%   x      - Struct containing pipeline environment parameters, also useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x      - Struct containing pipeline environment parameters, also useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Add maps and atlases to the x struct. This function is called by xASL_init_Process.m.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_init_MapsAndAtlases(x);
%
% __________________________________
% Copyright (c) 2015-2024 ExploreASL


%% --------------------------------------------------------
%% 1. Define directories
if isfield(x.opts, 'MyPath')
    % ExploreASL license
    x.D.MapsDir             = fullfile(x.opts.MyPath, 'Maps');
    x.D.TemplateDir         = fullfile(x.opts.MyPath, 'Maps', 'Templates');    

    % SPM modified, for image processing
    x.D.MapsSPMmodifiedDir  = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded');
    x.D.ResliceRef          = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'rgrey.nii');
    x.D.IdentityTransfRef   = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'Identity_Deformation_y_T1.nii');

    % Atlases 4 ROIs
    % This includes:
    % /LicenseLimited: several atlases with a limited license (only free for non-commercial usage)
    % /LicensePermissive: several atlases with a permissive license (free for any usage)
    % Inside /LicensePermissive, we also have a subfolder with atlases from SPM or CAT, which are GNU GPL licensed
    x.D.AtlasDir            = fullfile(x.opts.MyPath, 'External', 'Atlases4ROIs');
else
    warning('MyPath field not defined...');
end


%% --------------------------------------------------------
%% 2. Add all atlas/ROI NIFTIs to the x.D.Atlas field

% Get all NIfTI files in atlas subfolders (recursively)
% Note that these atlases need to have the .nii.mat sidecar
filesInAtlasDir = xASL_adm_GetFileList(x.D.AtlasDir, '^.+\.nii\.mat$', 'FPListRec');

% Iterate over atlases
for iFile=1:numel(filesInAtlasDir)
    % Get current atlas
    [fPath, currentAtlas] = xASL_fileparts(filesInAtlasDir{iFile});
    x.D.Atlas.(currentAtlas) = fullfile(fPath, [currentAtlas '.nii']);
end


end