function xASL_im_SkullStrip(InPath, PathMNIMask, x, OutPath)
%xASL_im_SkullStrip Creates skull-stripped T1w image based on
% MNI -> native space registration from segmentation
%
% FORMAT:       xASL_im_SkullStrip(InPath, PathMNIMask, x, OutPath)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Creates skull-stripped {{T1w}} image based on {{MNI}} -> native
%               space registration from segmentation.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL


    if nargin<4 || isempty(OutPath)
        [Fpath, Ffile, Fext] = xASL_fileparts(InPath);
        OutPath = fullfile(Fpath, ['mask_' Ffile Fext]);
    end

    if nargin<3
        x = [];
    end

    [Fpath, Ffile, Fext] = xASL_fileparts(OutPath);
    TempPath = fullfile(Fpath, ['Temp_' Ffile Fext]);

    %% Reslice MNI brainmask to native space
    xASL_spm_deformations(x, PathMNIMask, TempPath, 1, InPath); % trilinear interpolation & inverse

    %% Image arithmetics (masking)
    Image = double(xASL_io_Nifti2Im(InPath)).*(xASL_io_Nifti2Im(TempPath)>0);

    %% Save the masked image
    xASL_io_SaveNifti(InPath, OutPath, Image,[],0);
    xASL_delete(TempPath);
end
