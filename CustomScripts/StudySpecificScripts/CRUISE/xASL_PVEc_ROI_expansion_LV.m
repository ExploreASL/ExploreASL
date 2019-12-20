function [ ROIim ] = xASL_PVEc_ROI_expansion_LV( ROIim, GMim, WMim, WMdist)
%PVEc_ROI_expansion Dilates a GM ROI inwards to include more WM voxels
%   for partial volume correction
%   ExploreASL 2017

% NB: the ROI should be a mask, not a probability map

    %% Dilation loop
    % Here we dilate the ROI, until we have sufficient number of WM voxels
    % The dilation should grow inwards

    GMmask          = GMim>0.7;
    WMmask          = WMim>0.8;

    Iteration       = 1;
    MaxIt           = 100;
    
    GMmaskROIsum    = sum(sum(sum(GMmask(ROIim ))));
    WMmaskROIsum    = sum(sum(sum(WMmask(ROIim ))));
    GMmapROIsum     = sum(sum(sum(GMim(ROIim ))));
    WMmapROIsum     = sum(sum(sum(WMim(ROIim ))));
    
    if  xASL_stat_SumNan(ROIim(:))<296
        fprintf('%s','Empty or too small ROI, skipping PVEc expansion');
        % Sometimes there can be dummy ROIs, where we skip this part
    else
            %         fprintf('%s','PVEc ROI expansion, iteration   ');
            while   Iteration<MaxIt && GMmaskROIsum>(2*WMmaskROIsum) && ~(GMmaskROIsum<(4*WMmaskROIsum) && GMmapROIsum<WMmapROIsum)
                %                 xASL_TrackProgress(Iteration);
                 
                % Faster dilation that does it on cropped image
                tempMapAx = dilate_erode_full(ROIim,'dilate',[0 1 0; 1 1 1; 0 1 0]);
                tempMapSag = dilate_erode_full(ROIim,'dilate',cat(3,1,1,1));
                NewMask = (tempMapAx + tempMapSag) > 0;
                
                %
                NewMask                     = NewMask-ROIim; % select dilation edge only
                [IndVoxX,IndVoxY,IndVoxZ]   = ind2sub(size(NewMask),find(NewMask)); % find indexes
                
                % Check validity
                if  length(IndVoxX)~=sum(NewMask(:))
                    error('length(IndVoxX)~=sum(NewMask(:)), check whether ROI is binary mask');
                end
                
                % Every voxel has 6 neighbors
                NeiVoxX = repmat(IndVoxX,[1,6]);
                NeiVoxY = repmat(IndVoxY,[1,6]);
                NeiVoxZ = repmat(IndVoxZ,[1,6]);
                
                % Calculate the neighbors
                NeiVoxX = NeiVoxX+repmat([-1,0,1,0,0,0],[size(IndVoxX,1),1]);
                NeiVoxY = NeiVoxY+repmat([0,-1,0,1,0,0],[size(IndVoxY,1),1]);
                NeiVoxZ = NeiVoxZ+repmat([0,0,0,0,-1,1],[size(IndVoxZ,1),1]);
                
                % Check if they are still in the image
                % But this will probably not happen since the ROIs will never be close
                % to the image boundaries - especially when everything is in the MNI
                % space
                NeiVoxX(NeiVoxX<1) = 1;
                NeiVoxY(NeiVoxY<1) = 1;
                NeiVoxZ(NeiVoxZ<1) = 1;
                NeiVoxX(NeiVoxX>size(NewMask,1)) = size(NewMask,1);
                NeiVoxY(NeiVoxY>size(NewMask,2)) = size(NewMask,2);
                NeiVoxZ(NeiVoxZ>size(NewMask,3)) = size(NewMask,3);
                
                % loop across indexes (i.e. voxels from dilation)
                NewMask2Include             = zeros(size(NewMask));
                for iT=1:sum(NewMask(:))
                    % Get the distance of the voxel on the border
                    oldDistance = WMdist(IndVoxX(iT),IndVoxY(iT),IndVoxZ(iT));
                    
                    % The indexes of the neighbors - put them back to linear index
                    NeiInd = NeiVoxX(iT,:) + (NeiVoxY(iT,:)-1) * size(NewMask,1) + (NeiVoxZ(iT,:)-1) * size(NewMask,1) * size(NewMask,2);
                    
                    % Calculate their distances and if they were in the ROI
                    newDistance = WMdist(NeiInd).*ROIim(NeiInd);
                    % Exclude the zeros = those that were not in the ROI
                    newDistance = newDistance(newDistance>0);
                    NeighborDistance = mean(newDistance);
                    
                    % Compare distance of SingleVoxel and NeighborVoxel
                    % If distance of SingleVoxel is smaller than NeighborVoxel (i.e.
                    % the dilated voxel is closer to the skeleton of the WM) it will be
                    % included
                    
                    if  oldDistance<=NeighborDistance % only include voxel if it is not farther away from the WM skeleton as the original ROI
                        NewMask2Include(IndVoxX(iT),IndVoxY(iT),IndVoxZ(iT)) = 1;
                    end
                end
                ROIim   = ROIim+NewMask2Include;
                
                Iteration   = Iteration+1;
                GMmaskROIsum    = sum(sum(sum(GMmask(logical(ROIim )))));
                WMmaskROIsum    = sum(sum(sum(WMmask(logical(ROIim )))));
                GMmapROIsum     = sum(sum(sum(GMim(logical(ROIim )))));
                WMmapROIsum     = sum(sum(sum(WMim(logical(ROIim )))));
            end;% while Iteration
        
    end % if  xASL_stat_SumNan(ROIim(:))==0
%     fprintf('\n');
end
