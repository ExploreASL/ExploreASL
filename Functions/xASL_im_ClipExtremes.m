function [NewIM] = xASL_im_ClipExtremes(InputIm, ThreshHigh, ThreshLow, bVerbose)
% xASL_im_ClipExtremes Clips image to given percentile. The percentile is found using non-zeros sorted intensities
% , so both isfinite & non-zeros
%
% FORMAT:       [NewIM] = xASL_im_ClipExtremes(InputIm, ThreshHigh, ThreshLow, bVerbose)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Clips image to given percentile. The percentile is found
%               using non-zeros sorted intensities, so both isfinite & non-zeros.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

    %% Admin
    if nargin<4 || isempty(bVerbose)
        bVerbose = true; % default is output clipped intensities
    end

    if nargin<2 || isempty(ThreshHigh)    
        % To disable high clipping, set to 1
        fprintf('No upper threshold set, set to 0.999');
        ThreshHigh = 0.999;
    end

    if nargin<3 || isempty(ThreshLow)
        fprintf('No lower threshold set, set to 1-high threshold\n');
        % By default, ThreshLow is 1-ThreshHigh (i.e. 0.001)
        % To disable, set ThreshLow to 0
           ThreshLow = 1-ThreshHigh;
    elseif ThreshLow == 1 % backward compatibility
           ThreshLow = 1-ThreshHigh;
    end    
    
    %% ------------------------------------------------------------------------------------------
    % Constrain clippable intensities
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

    % Clip high intensities
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
        % Clip low intensities
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
    % Save as NIfTI if the input was a NIfTI
    if xASL_exist(InputIm,'file')
        xASL_io_SaveNifti(InputIm, InputIm, NewIM, [], false);
    end

end
