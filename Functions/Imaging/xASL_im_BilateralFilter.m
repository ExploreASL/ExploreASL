function [ovol] = xASL_im_BilateralFilter(volIM, mask, VoxelSize, x)
%xASL_im_BilateralFilter This function runs a spatial lowpass temporally
%highpass filter, and removes outliers within this signal, and adapts the
%time-series accordingly
%
% FORMAT:       [ovol] = xASL_im_BilateralFilter(volIM, mask, VoxelSize, x)
% 
% INPUT:        volIM      - Image volume (REQUIRED)
%               mask       - Mask (REQUIRED)
%               VoxelSize  - Voxel size (REQUIRED)
%               x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:       ovol       - Output volume
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function runs a spatial lowpass temporally
%               highpass filter, and removes outliers within this signal, and adapts the
%               time-series accordingly.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2021 ExploreASL

if  size(volIM,4)<10
    % skip filtering, for no time-series
    fprintf('%s\n',['Filtering skipped because of only ' num2str(size(volIM,4)) 'averages/dynamics/subtractions']);
    ovol    = volIM;

elseif  x.settings.BILAT_FILTER>0

        %% Get label order
        zerolabel   = zeros(1,size(volIM,4)/2);
        onelabel    = ones(1,size(volIM,4)/2);

        [~, ~, OrderContLabl]   = xASL_quant_GetControlLabelOrder(volIM);

        if ~OrderContLabl % control is higher than label
            label       = [zerolabel(1) onelabel(1)];

            for iI=2:size(volIM,4)/2
                label       = [label zerolabel(iI) onelabel(iI)];
            end
        else % label is higher than control
            label       = [onelabel(1) zerolabel(1)];

            for iI=2:size(volIM,4)/2
                label       = [label onelabel(iI) zerolabel(iI)];
            end
        end



        %% Run filter
        if ~deployed % skip this for compilation, to avoid DIP image issues
            if      x.settings.BILAT_FILTER==1
                    % original bilateral filter by Matthan Caan that filters raw
                    % images/timeseries
                    % Temporarily disabled
                    % ovol    = xASL_bilateralFilter(volIM,VoxelSize,double(squeeze(label')),double(mask) );
                    ovol = volIM;
                    % old version
                    %ovol    = bilateralFilterASL_HJM_LV(volIM,VoxelSize,double(squeeze(label')),double(mask) );
                    % Temporarily use this filter

            elseif  x.settings.BILAT_FILTER==2
                    % more recent bilateral filter by Matthan Caan that filters subtracted
                    % images/timeseries (==n/2)
                    % Temporarily disabled
                    % [ovol,varargout]    = extrapolateWMFilterASL_HJM(volIM,VoxelSize,double(squeeze(label')),double(mask) );
                    ovol = volIM;
            end
        else
            ovol = volIM;
        end



end

end
