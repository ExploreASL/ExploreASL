function [X Y Z] = xASL_adm_Voxel2RealWorldCoordinates(X,Y,Z,VoxelSize)
%xASL_adm_Voxel2RealWorldCoordinates ...
%
% FORMAT:       [X Y Z] = xASL_adm_Voxel2RealWorldCoordinates(X,Y,Z,VoxelSize)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Converts MNI coordinates from voxel coordinates/indices.
%               Assumes X Y Z = LR LeftRight AP AnteriorPosterior IS InferiorSuperior.
%               VoxelSize should be [1 3]-sized input.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL
    
    % World Coordinates:
    % X = -90:90  mm LR (x.S.SagSlices) -> this is mirrored, see below
    % Y = -126:90 mm AP (x.S.CorSlices)
    % Z = -72:108 mm IS (x.S.TraSlices)
    
    
    if ~exist('VoxelSize','var')
        VoxelSize   = [1.5 1.5 1.5]; % ExploreASL toolbox default
    end

    % X = slice 1 = -90 (right)
    % slice 121 = 90 (left)
    % 0 mm = midline (slice 60), every slice 1.5 mm
    % -90+120*1.5

    if ~exist('X','var')
        X   = [];
    else

        X           = -90+VoxelSize(1).*(X-1);
        X           = -X; % practical implementation to show exact same numbers
                          % as SPM, this left right may be different in reality
        % This is now fixed. with dip_image, left is under, right is above,
        % with MriCron, left is left & right is right (MNI == neurological orientation)
        % SPM loads NIfTIs voxels in radiological orientation
        % Harvard-Oxford atlas (subcortical) notation of left is identical to
        % SPM coordinates (X=-) and right (X=+)

        % below=Left 	= (61:end,:,:)
        % Up   =Right 	= ( 1: 60,:,:)

        % Y = slice 1 = -126 (posterior)
        % slice 145 = 90 (anterior)
        % 0 mm = just anterior to thalamus (slice 84), every slice 1.5 mm
        % -126+144*1.5 = 90
    end

    if ~exist('Y','var')
        Y   = [];
    else
        Y           = -126+VoxelSize(2).*(Y-1);

        % Z:
        % Slice 1   = -72 mm; inferior
        % Slice 121 = 108 mm; superior
        % 0 mm = ACPC line (slice 48), every slice 1.5 mm
        % -72+1.5*120 = 108;
    end

    if ~exist('Z','var')
        Z   = [];
    else    
        Z           = -72+VoxelSize(3).*(Z-1);
    end

end