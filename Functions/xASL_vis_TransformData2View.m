function FigureOut = xASL_vis_TransformData2View(ImagesIn, x)
%xASL_vis_TransformData2View Reshapes MNI image matrix data into visualization figure
%
% FORMAT: FigureOut = xASL_vis_TransformData2View(ImagesIn, x)
%
% INPUT:
%   ImagesIn - single or series of images that need to be reshaped (REQUIRED)
%              input format needs to be a matrix, with [121 145 121] (i.e. 1.5mm MNI)
%              images, can be concatenated over first or fourth dimension
%   x        - structure containing fields with generic information from ExploreASL (OPTIONAL, 
%              default tiles a few transversal slices)
%              ORIENTATION SETTINGS:
%               x.S.TraSlices - which transversal slices to show (n, orientation omitted if empty, DEFAULT = 20:7:97)
%               x.S.CorSlices - which coronal slices to show (n, orientation omitted if empty, DEFAULT = empty)
%               x.S.SagSlices - which sagittal slices to show (n, orientation omitted if empty, DEFAULT = empty)
%              CROP SETTINGS:
%               x.S.bCrop - number of voxels you want to crop from the default size, single value,
%                           ranging from 50 (pieces of brain cut out) to -Inf (creating space between images) (DEFAULT=0)
%              CONCATENATION SETTINGS:
%              x.S.ConcatSliceDims % whether orientations should be concatenated horizontally (1, DEFAULT) or vertically (0)
%              x.S.Square % whether we prefer to tile square (DEFAULT = true)
%
% OUTPUT:
%   FigureOut - reshaped image or series of images
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function changes the dimensionality and reshapes the input images
% in such a way that they are nicely tiled in a mosaic for visualization purposes.
% Reshaping a series of images with this function can be useful for
% visualization of SPM/voxel-based analyses.
%
% EXAMPLE: FigureOut = xASL_vis_TransformData2View(ImagesIn);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Admin

if nargin<2 || isempty(x)
    x = struct;
end
if ~isfield(x,'S')
    x.S = struct;
end

% Remove empty orientation fields (enabling to concatenate only requested slices)
OriNames = {'TraSlices' 'CorSlices' 'SagSlices'};
for iO=1:length(OriNames)
    if isfield(x.S, OriNames{iO}) && isempty(x.S.(OriNames{iO}))
        x.S = rmfield(x.S, OriNames{iO});
    end
end

% Define default settings
%   If there is no dimensionality specified at all, default is a
%   rectangular grid with transversal slices
if ~isfield(x.S,'TraSlices') && ~isfield(x.S,'CorSlices') && ~isfield(x.S,'SagSlices')
    x.S.TraSlices = 20:7:97; % by default show 9 transversal slices
end
if ~isfield(x.S,'ConcatSliceDims')
    x.S.ConcatSliceDims = true; % horizontal concatenation
end
if ~isfield(x.S,'Square')
    x.S.Square = true; % default is to try to create a square tile
end

% Manage ImagesIn for cell
if ~iscell(ImagesIn) % make cell if it is only a matrix, to accommodate both data types
    IsCell = false;
    TempData{1} = ImagesIn;
else
    TempData = ImagesIn;
    IsCell = true;
end

%% Cropping optimization
%  This assumes MNI. If this script is used with other space (e.g. native space) images,
%  best/easiest to disable cropping
if ~isfield(x.S,'bCrop') || isempty(x.S.bCrop)
    x.S.bCrop = 0; % crop zero voxels by default
end
x.S.TransCrop = [12 133 0 121 7 128] + x.S.bCrop .* [1 -1 1 -1 1 -1];

%% COMMENT OUT THE STUFF HERE IF THIS SCRIPT CRASHES
% ResultingSizeIs = [x.S.TransCrop(2)-x.S.TransCrop(1), x.S.TransCrop(4)-x.S.TransCrop(3), x.S.TransCrop(6)-x.S.TransCrop(5)];
% Here we crop a bit more, depending on the orientations
CroppingScheme = [0  0 10 -10 10 -10];

if isfield(x.S,'SagSlices') && ~x.S.ConcatSliceDims
       % vertical concatening with SagSlices
       % we can't crop LR but we can crop IS
       x.S.TransCrop = x.S.TransCrop + ([0 0 0 0 1 1].*CroppingScheme);
elseif isfield(x.S,'TraSlices') && x.S.ConcatSliceDims
       % horizontal concatening with TraSlices
       % we can't crop IS but we can crop LR
       x.S.TransCrop = x.S.TransCrop + ([0 0 1 1 0 0].*CroppingScheme);
       % otherwise we concatenate both LR & IS
else
       x.S.TransCrop = x.S.TransCrop + ([0 0 1 1 1 1].*CroppingScheme);
end


%% Which orientations to visualize
DIMS = {'Tra' 'Cor' 'Sag'};
DIMn(1:3) = 0;
for iDim=1:3
    DimField = [DIMS{iDim} 'Slices'];
    if isfield(x.S, DimField) && ~isempty(x.S.(DimField))
        DIMn(iDim)=1;
    end
end

% % % 
% % % 
% % % if nMultiSlice
% % %     % Number of slices need to be equal, to be able to concatenate
% % %     warning('Equal number of slices required in each dimension for concatenation, skipping...');
% % %     FigureOut{1} = NaN;
% % %     return;
% % % end    
% % % 



% Note dims
DimsTotal = 0;
NextN = 1;
for iDim=1:3
    if DIMn(iDim)
        DimsTotal(NextN) = iDim;
        NextN = NextN+1;
    end
end

%% =====================================================
%% Reshape data into visualization dimension/orientation
for iSet=1:length(TempData)
    DimData = size(TempData{iSet});
    % First accommodate compressed columns
    if  DimData(3)==1 
        TempData{iSet} = xASL_im_Column2IM(TempData{iSet}, x.S.masks.WBmask);
    end

    % Then accommodate both single & multiple subjects.
    % By default this function expects multiple subjects,
    % if there is only one subject, we create a dummy dimension for this
    if length(DimData)<4
        TempData2{iSet}(1,:,:,:) = TempData{iSet};
        TempData{iSet} = TempData2{iSet};
        DimShift = false;
    else
        CheckDim = DimData==[121 145 121 1];
        CheckDim2 = DimData==[1 121 145 121];
        if sum(CheckDim(1:3))==3
            DimShift = true;
            TempData{iSet} = shiftdim(TempData{iSet},3);
        elseif sum(CheckDim2(2:4))==3
            DimShift = false;
        else
            warning('Input image(s) should be [121 145 121] matrix, skipping...');
            continue;
        end
    end
    
    for iSubject=1:size(TempData{iSet},1)
        % Collect image slices for all dimensions

        % run cropping
        if DIMn(1)
            TempIm{1} = xASL_im_rotate(squeeze(TempData{iSet}(iSubject,:,:, x.S.TraSlices)),90);
            CropParms = x.S.TransCrop(1:4);
            TempIm{1} = xASL_vis_CropParmsApply(TempIm{1}, CropParms);
        end
        if  DIMn(2)
            TempIm{2} = xASL_im_rotate(xASL_vis_FlipOrientation(squeeze(TempData{iSet}(iSubject,:,x.S.CorSlices,:))),90);
            CropParms = [x.S.TransCrop(5),x.S.TransCrop(6),x.S.TransCrop(3),x.S.TransCrop(4)];
            TempIm{2} = xASL_vis_CropParmsApply(TempIm{2}, CropParms);
            % Pragmatic correction for strange flipping if nSlices==1
            if length(x.S.CorSlices)==1
                TempIm{2} = xASL_im_rotate(TempIm{2},180);
            end
        end
        if  DIMn(3)
            TempIm{3} = xASL_vis_FlipOrientation2(squeeze(TempData{iSet}(iSubject,x.S.SagSlices,:,:)));
            CropParms = [x.S.TransCrop(5),x.S.TransCrop(6),x.S.TransCrop(1),x.S.TransCrop(2)];
            TempIm{3} = xASL_vis_CropParmsApply(TempIm{3}, CropParms);
            % Pragmatic correction for strange flipping if nSlices==1
            if length(x.S.SagSlices)==1
                TempIm{3} = xASL_im_rotate(TempIm{3},180);
            end
        end

        % Determine columns of final image
        nColumnsTot = 0;
        for iDim=DimsTotal
            IndColumns = size(TempIm{iDim},3);
            nColumnsTot = nColumnsTot + size(TempIm{iDim},3);
        end

        % Concatenation has to do with the ratio between number of
        % individual slices and sqrt of number of total slices
        % So this will prefer the number of individual slices,
        % until the number of individual slices becomes too big
        % to make a square
        if x.S.ConcatSliceDims % horizontal concatination
            if  x.S.Square
                nColumns = floor(sqrt(nColumnsTot))/length(DimsTotal);
            else
                nColumns = nColumnsTot/length(DimsTotal);
            end
        else % vertical concatenation
            if x.S.Square
                nColumns = ceil(nColumnsTot^0.5); % make a square, bias width (prefer width above height)
            else
                nColumns = nColumnsTot;
            end
            nColumns = (ceil((nColumns/IndColumns)*2)/2) *IndColumns; % bias individual n slices, 
            % only deviate from individual number of slices when too many                
        end

        % Automatic row/colum size determination
        if length(DimsTotal)==1 % for single direction, keep it simple
            if x.S.Square % no preference
                nColumns = ceil(nColumnsTot^0.5);
            elseif x.S.ConcatSliceDims % prefer horizontal
                nColumns = nColumnsTot;
            else % prefer vertical
                nColumns = 1;
            end
        end

        if isfield(x.S,'nColumns')
            nColumns = x.S.nColumns; 
            % by specifying this we override the above automatic column creation
        end

        % Check for equal number of slices for requested dims
        UnequalSlices = false;
        if DIMn(1) && DIMn(2) && length(x.S.TraSlices)~=length(x.S.CorSlices)
               UnequalSlices = true;
        elseif DIMn(1) && DIMn(3) && length(x.S.TraSlices)~=length(x.S.SagSlices)
               UnequalSlices = true;
        elseif DIMn(2) && DIMn(3) && length(x.S.CorSlices)~=length(x.S.SagSlices)
               UnequalSlices = true;
        end
        
        % put data into "mosaic" form
        % if we have unequal number of slices between orientations, then 
        % simply stack them in the same axis as we stack the orientations
        % (i.e. we only have 1 row or column)
        % otherwise, try to stack them in the opposite direction (e.g. for
        % rows or columns with multiple slices from an orientation, where
        % rows or columns are the orientations)
        for iDim=DimsTotal
            if UnequalSlices && ~x.S.ConcatSliceDims
                TempIm{iDim} = xASL_vis_TileImages(TempIm{iDim}, 1);
            elseif UnequalSlices && x.S.ConcatSliceDims
                TempIm{iDim} = xASL_vis_TileImages(TempIm{iDim}, size(TempIm{iDim},3));
            else % keep default
                TempIm{iDim} = xASL_vis_TileImages(TempIm{iDim}, nColumns);
            end
        end
        
        % Concatenate them            
        NextN = 1;
        for iDim=DimsTotal
            if NextN==1 % start new image
                DIMs2Combine = TempIm{iDim};
            else % concatenate with existing image
                if x.S.ConcatSliceDims % concatenate horizontally
                    DIMs2Combine = [DIMs2Combine TempIm{iDim}];
                else % concatenate vertically
                     DIMs2Combine = [DIMs2Combine; TempIm{iDim}];
                end
            end
            NextN = NextN+1; % count                    
        end

        FigureOut{iSet}(:,:,iSubject) = DIMs2Combine;
        clear DIMs2Combine TempIm
        clear nColumns nColumnsTot IndColumns
    end % for iSubject
end

if ~exist('FigureOut', 'var')
    FigureOut{1} = NaN;
end

if ~IsCell % return cell to matrix dataform, to accommodate both data types
    FigureOut = FigureOut{1}; % default output for if function crashes
end


end




%% ===================================================================================================
%% ===================================================================================================



function image_out = xASL_vis_FlipOrientation(image_in)
%FlipOrientation This function flips the 3 dimensions from sagittal to
%transversal or tra to cor. Leaves other dimensions untouched.


image_out   = shiftdim(xASL_im_rotate(image_in,90),1);

% OLD CODE
% temp=single(xASL_im_rotate(image_in,90));
% dim=size(temp);
% image_out=zeros(dim(2),dim(3),dim(1),size(temp,4));
% 
% 
% 
% for l=1:size(temp,4)
%     for i=1:dim(2)
%         for j=1:dim(3)
%             for k=1:dim(1)
%                 image_out(i,j,k,l)=temp(k,i,j,l);
%             end
%         end
%     end
% end


%new_y=old_z dim(3)->dim(2)
%new_x=old_y dim(2)->dim(1)
%new_z=old_x dim(1)->dim(3)



end

%% ===================================================================================================
%% ===================================================================================================


function image_out = xASL_vis_FlipOrientation2(image_in)
%FlipOrientation This function flips the 3 dimensions from sagittal to
%cor or tra to sag. Leaves other dimensions untouched.


image_out               = xASL_im_rotate(shiftdim(xASL_im_rotate(image_in,90),2),180);
image_out               = image_out(:,:,size(image_out,3):-1:1); % flip


% 
% 
% temp=xASL_im_rotate(single(image_in),90);
% dim=size(temp);
% image_out=zeros(dim(3),dim(1),dim(2),size(temp,4));
% 
% for l=1:size(temp,4)
%     for i=1:dim(3)
%         for j=1:dim(1)
%             for k=1:dim(2)
%                 image_out(i,j,k,l)=temp(j,k,i,l); % new z-direction is magnified/blurred by 2. 2*i-1 = old i & 2*i is average of old i & old (i+1)
%             end
%         end
%     end
% end
% 
% image_out=xASL_im_rotate(image_out,180);
% 
%new_y=old_z dim(3)->dim(2)
%new_x=old_y dim(2)->dim(1)
%new_z=old_x dim(1)->dim(3)



end


