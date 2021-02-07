function [NewIM] = xASL_im_ClipExtremes(InputIm, ThreshHigh, ThreshLow, bVerbose, bNormalize)
%xASL_im_ClipExtremes Clips image to threshold
%
% FORMAT:  [NewIM] = xASL_im_ClipExtremes(InputIm[, ThreshHigh, ThreshLow, bVerbose])
% 
% INPUT:
%  InputIm      - path to image or image matrix (REQUIRED)
%  ThreshHigh   - upper clipping threshold [0-1]
%                 e.g. set to 1 for no clipping
%                 (OPTIONAL, DEFAULT=0.999)
%  ThreshLow    - lower clipping threshold [0-1]
%                 e.g. set to 0 for no clipping
%                 (OPTIONAL, DEFAULT=1-ThreshHigh)
%  bVerbose     - boolean specifying verbosity
%                 (OPTIONAL, DEFAULT=true)
%  bNormalize   - normalizes maximum value to 4096 (12 bit max, 12^2)
%                 This can help with anatomical images that are incorrectly
%                 exported by the scanner with very high values.
%                 This is disabled by default, for quantitative images.
%                 (OPTIONAL, DEFAULT=false)
%
% OUTPUT:
%  NewIM        - clipped image (OPTIONAL)
% OUTPUT FILES:
%  NIfTI containing the clipped image, if InputIm was a path
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function clips an image to a given percentile. The percentile is found
%               using non-zeros sorted intensities, so both isfinite & non-zeros.
%               This function performs the following steps:
%
%               1. Constrain clippable intensities
%               2. Clip high intensities
%               3. Clip low intensities
%               4. Normalize to 4096 (12 bit, 12^2)
%               5. Save as NIfTI if the input was a NIfTI
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_im_ClipExtremes('MyStudy/anat/T1w.nii.gz');
% __________________________________
% Copyright 2015-2020 ExploreASL


    %% 0. Admin
    if nargin<5 || isempty(bNormalize)
        bNormalize = false;
    end
    
    if nargin<4 || isempty(bVerbose)
        bVerbose = true; % default is output clipped intensities
    end

    if nargin<2 || isempty(ThreshHigh)    
        % To disable high clipping, set to 1
        if bVerbose
            fprintf('No upper threshold set, set to 0.999\n');
        end
        ThreshHigh = 0.999;
    end

    if nargin<3 || isempty(ThreshLow)
        if bVerbose
            fprintf('No lower threshold set, set to 1-high threshold\n');
        end
        % By default, ThreshLow is 1-ThreshHigh (i.e. 0.001)
        % To disable, set ThreshLow to 0
           ThreshLow = 1-ThreshHigh;
    elseif ThreshLow == 1 % backward compatibility
           ThreshLow = 1-ThreshHigh;
    end    
    
    %% ------------------------------------------------------------------------------------------
    % 1. Constrain clippable intensities
    if ThreshLow~=0
        if ThreshHigh<0.5
            ThreshHigh=0.5;
        elseif ThreshHigh>1
            ThreshHigh=1;
        end
    elseif ThreshHigh<0
            ThreshHigh=0;
    elseif ThreshHigh>1
            ThreshHigh=1;
    end



    %% ------------------------------------------------------------------------------------------
    NewIM = xASL_io_Nifti2Im(InputIm);
    SortedInt = sort(NewIM(NewIM~=0 & isfinite(NewIM)));

    % 2. Clip high intensities
    IndexNumber = round(ThreshHigh*length(SortedInt));
    if IndexNumber<=length(SortedInt) && IndexNumber>0

        maxInt = max(NewIM(:));
        ClipInt = SortedInt(round(ThreshHigh*length(SortedInt)));
        NewIM(NewIM>ClipInt) = ClipInt;
        maxIntNew = max(NewIM(:));
        if bVerbose
            fprintf('%s\n',['Clipped maximal intensity from ' num2str(maxInt) ' to ' num2str(maxIntNew)]);
        end
    end



    %% ------------------------------------------------------------------------------------------
    if ThreshLow~=0
        % 3. Clip low intensities
        IndexNumber = round(ThreshLow*length(SortedInt));
        if IndexNumber>0
            minInt = min(NewIM(:));
            ClipIntLow = SortedInt(round(ThreshLow*length(SortedInt)));
            NewIM(NewIM<ClipIntLow) = ClipIntLow;
            minIntNew = min(NewIM(:));
            if bVerbose
                fprintf('%s\n',['Clipped minimal intensity from ' num2str(minInt) ' to ' num2str(minIntNew)]);
            end
        end
    end

    
    %% ------------------------------------------------------------------------------------------    
    % 4. Normalize to 4096
    if bNormalize
        NewIM = single(NewIM); % avoid clipping artifacts
        MaxIM = max(NewIM(:));
        NewIM = 4096.*(NewIM./MaxIM);
    end


    %% ------------------------------------------------------------------------------------------
    % 5. Save as NIfTI if the input was a NIfTI
    if ischar(InputIm) && xASL_exist(InputIm, 'file')
        [~, ~, Fext] = xASL_fileparts(InputIm);
        if strcmp(Fext, '.nii.gz')
            xASL_io_SaveNifti(InputIm, InputIm, NewIM, [], true);
        else
            xASL_io_SaveNifti(InputIm, InputIm, NewIM, [], false);
        end
    end

end
