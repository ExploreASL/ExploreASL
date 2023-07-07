function x = xASL_init_MapsAndAtlases(x)
%xASL_init_MapsAndAtlases Add maps and atlases
%
% FORMAT: x = xASL_adm_CleanUpX(x)
%
% INPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Add maps and atlases.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_init_MapsAndAtlases(x);
%
% __________________________________
% Copyright (c) 2015-2022 ExploreASL

if isfield(x.opts, 'MyPath')
    x.D.MapsDir             = fullfile(x.opts.MyPath, 'Maps');
    x.D.MapsSPMmodifiedDir  = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded');
    x.D.ResliceRef          = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'rgrey.nii');
    x.D.IdentityTransfRef   = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'Identity_Deformation_y_T1.nii');
    x.D.TemplateDir         = fullfile(x.opts.MyPath, 'Maps', 'Templates');
    x.D.AtlasDir            = fullfile(x.opts.MyPath, 'External', 'Atlases');
else
    warning('MyPath field not defined...');
end

end


