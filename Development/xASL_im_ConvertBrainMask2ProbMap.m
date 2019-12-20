function xASL_im_ConvertBrainMaskProbMap( InPath, OutPath, SteepnessFactor, InitialDilations )
%xASL_im_ConvertBrainMaskProbMap Converts a brainmask to a probability map by
% repetitive dilation
% SteepnessRate defines how fast the mask falls down to zero
% Assume a 3D image
%
% This function is handy if we don't trust non-linear registration
% to create a BrainMask from MNI completely, but for let's say 80%.
% Then gradually decreasing image intensity still removes the stuff outside
% the skull, but in a gentle way, limiting the effect of slight misalignments

    % Admin

    % Calculate steepness
    % With ASL image matrix (108800 voxels) we want SteepnessRate=0.5
    % With T1w image matrix (18489600 voxels) we want SteepnessRate=0.05;




    %% Open NIfTI
    IM              = xASL_io_ReadNifti(InPath);
    IM              = IM.dat(:,:,:); % assume a 3D image
    IM(isnan(IM))   = 0; % deal with NaNs
    IM              = single(IM>(0.5.*max(IM(:))) ); % convert map to mask

    %% SteepnessFactor intialization
    if     ~exist('SteepnessFactor', 'var')
            SteepnessFactor      = 3;
    elseif  SteepnessFactor<0 % clip at 0
            SteepnessFactor     = 0;
    end

    SteepnessRate    = 10^6/numel(IM) .* SteepnessFactor;

    %% InitialDilations
    % This allows to dilate the mask a bit before gradually masking
    if ~exist('InitialDilations', 'var')
        InitialDilations    = 0;
    end

    % Do initial dilations
    for iD=1:InitialDilations
        IM         = xASL_im_DilateErodeSeparable(IM,'dilate',[1 1 1 1 1],[1 1 1 1 1],[1]);
    end


%% Start map creation

    ItN             = 1;
    MaskCheck       = IM;
    while sum(sum(sum(logical(MaskCheck)==0)))>0 % go until no voxels left
        ItN        = ItN-SteepnessRate;
        if  ItN<0
            ItN    = 0;
        end
        Im2         = logical(MaskCheck);
        Im2         = xASL_im_DilateErodeSeparable(Im2,'dilate',[1 1 1 1 1],[1 1 1 1 1],[1]);
        Im2         = shiftdim(Im2,1);
        Im2         = xASL_im_DilateErodeSeparable(Im2,'dilate',[1 1 1 1 1],[1 1 1 1 1],[1]);
        Im2         = shiftdim(Im2,2);
        Diff        = Im2-logical(MaskCheck);
        IM          = IM+(ItN.* Diff );
        MaskCheck   = MaskCheck+logical(Diff);
    end

    xASL_io_SaveNifti( InPath, OutPath, IM );


end
