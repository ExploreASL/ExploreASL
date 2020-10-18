function xASL_im_Upsample(PathOrig, PathDest, NewVoxelSize, LeaveEmpty, PaddingDim, Kernel)
%xASL_im_Upsample upsamples an ASL image, without changing the orientation
%matrix, which can be used e.g. for PVEc in higher resolution but same space
%
% FORMAT:       xASL_im_Upsample(PathOrig, PathDest, NewVoxelSize, LeaveEmpty, PaddingDim, Kernel)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Upsamples an ASL image, without changing the orientation
%               matrix, which can be used e.g. for PVEc in higher
%               resolution but same space.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

    %% Admin

    xASL_io_ReadNifti(PathOrig);
    volASL      = spm_vol(PathOrig);
    imASL       = spm_read_vols(volASL);

    matLow      = volASL.mat;
    dimLow      = size(imASL);
    imLow       = imASL;

    if length(dimLow)>3 || length(dimLow)>NewVoxelSize
        error([PathOrig ' has too many dimensions']);
    end        
    
    if nargin<4 || isempty(LeaveEmpty)
        LeaveEmpty = false;
        % by default we upsample the ASL image
        % if this parameter is set to true we simply create an empty (zeros)
        % matrix
    end

    if nargin<4 || isempty(PaddingDim)
        PaddingDim = [0 0 0];
        % This parameter divided by 2, will be padded zeros at either side
        % of the image matrix (in image space, not in Fourier/frequency domain)
    end

    if nargin<6 || isempty(Kernel)
        Kernel = 'linear';
    end

    %% Petr ASL upsample function
    % [imHighASL,tmp,matHigh] = aslUpsample(imASL,size(imASL),volASL.mat,[1 1 1]);

    % First voxel in NIFTII is [1,1,1] - this refers to the center of the voxel
    % imLow - low resolution image
    % dimLow - size of the low resolution image
    % matLow - the orientation matrix from nifti
    % resHigh - resolution of the upsampled image
    % imHigh - the upsampled image
    % dimHigh - size of the upsampled image in pixels
    % matHigh - the upsampled orientation matrix

    matHigh = matLow;
    dimHigh = dimLow;
    VoxDim  = {'X' 'Y' 'Z'};

    %% Calculate the new dimension and transformation matrix
    for iDim=1:3
        OriginalRes(iDim)   = norm(matLow(1:3,iDim));
        matHigh(1:3,iDim)   = matHigh(1:3,iDim)./OriginalRes(iDim)*NewVoxelSize(iDim);

        % Insert extra voxel to have a slightly larger FOV
        dimHigh(iDim)       = (dimHigh(iDim)*OriginalRes(iDim)/norm(matHigh(1:3,iDim)))+PaddingDim(iDim);

        % The new dimensions have to be rounded and ROUND or CEIL is
        % selected based on which will be more stable
        % I.e. ROUND for numbers close to integers as even for some variation ROUND(19.9), ROUND (20.2) will return 20.
        % and CEIL for others as CEIL (20.4) and CEIL(20.6) will return 21.
        % While ROUND(20.4) = 20 and ROUND(20.6) = 21.
        if abs(dimHigh(iDim) - round(dimHigh(iDim))) < 0.25
            dimHigh(iDim) = round(dimHigh(iDim));
        else
            dimHigh(iDim) = ceil(dimHigh(iDim));
        end;
    end

    fprintf('%s\n',['Voxelsize [' num2str(OriginalRes(1)) ' ' num2str(OriginalRes(2)) ' ' num2str(OriginalRes(3)) '] mm changed to [' num2str(NewVoxelSize(2)) ' ' num2str(NewVoxelSize(2)) ' ' num2str(NewVoxelSize(3)) '] mm']);

    % Change the offset to accommodate the extra pixels
    % This for cycle needs to be run after the full previous for-cycle is
    % finished
    for iDim=1:3
        centerOld = (dimLow(iDim) + 1) / 2 * matLow(1:3,iDim);
        centerNew = (dimHigh(iDim) + 1) / 2 * matHigh(1:3,iDim);
        matHigh(1:3,4) = matHigh(1:3,4) + centerOld - centerNew;
    end

    %% Calculate the old and new coordinates for center of the image
    if  LeaveEmpty
        imHigh              = ones(dimHigh);
        SaveBit             = 8; % low quality, doesn't matter
    else
        % Interpolate image
        SaveBit             = 32; % high quality
        Xold = 0:(dimLow(2)-1);
        Yold = 0:(dimLow(1)-1);
        Zold = 0:(dimLow(3)-1);
        [Xnew,Ynew,Znew] = meshgrid(0:(dimHigh(2)-1),0:(dimHigh(1)-1),0:(dimHigh(3)-1));

        Xnew = (Xnew*norm(matHigh(1:3,2)) - ((dimHigh(2)-1)/2)*norm(matHigh(1:3,2)) + ((dimLow(2)-1)/2)*norm(matLow(1:3,2)))/norm(matLow(1:3,2));
        Ynew = (Ynew*norm(matHigh(1:3,1)) - ((dimHigh(1)-1)/2)*norm(matHigh(1:3,1)) + ((dimLow(1)-1)/2)*norm(matLow(1:3,1)))/norm(matLow(1:3,1));
        Znew = (Znew*norm(matHigh(1:3,3)) - ((dimHigh(3)-1)/2)*norm(matHigh(1:3,3)) + ((dimLow(3)-1)/2)*norm(matLow(1:3,3)))/norm(matLow(1:3,3));

        imHigh = interp3(Xold,Yold,Zold,imLow,Xnew,Ynew,Znew,Kernel);
        imHigh(isnan(imHigh)) = 0;
    end

    volASL.dim                 = size(imHigh);
    volASL.dt                  = [16 0];
    volASL.mat                 = matHigh;
    volASL.pinfo               = [1;0;0];
    volASL.fname               = PathDest;
    volASL.private.mat         = volASL.mat;

	if xASL_exist(volASL.fname)
		xASL_delete(volASL.fname);
	end
    spm_write_vol(volASL,imHigh);

    if SaveBit==8
       xASL_io_SaveNifti(PathDest,PathDest,imHigh,size(imHigh,4),8);
    end


end
