function imOutput = xASL_im_ResampleLinear(imInput, newSize)
%xASL_im_ResampleLinear Downsample or upsample a 1D/2D/3D image
%
% FORMAT:       imOutput = xASL_im_ResampleLinear(imInput, newSize)
%
% INPUT:        imInput     - Image matrix (REQUIRED, DOUBLE, SINGLE or INT)
%               newSize     - Size of ouput image (REQUIRED, INTEGER ARRAY)
%
% OUTPUT:       imOutput    - Resampled image matrix (SINGLE)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Downsample or upsample an image from its old to a new resolution.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      image = xASL_io_Nifti2Im('.\M0.nii.gz');
%               output = xASL_im_ResampleLinear(image, [2,2,2]);
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Input check
    if nargin<1 || isempty(imInput)
        error('Missing input image...');
    end
    if numel(size(imInput))>3
        warning('Input image is not one, two or three dimensional...');
        imOutput = imInput;
        return
    end
    if nargin<2 || isempty(newSize)
        warning('Missing new size to resample image...');
        imOutput = imInput;
        return
    end
    if numel(newSize)>3
        warning('New image size is not one, two or three dimensional...');
        imOutput = imInput;
        return
    end
    if numel(newSize)>1 && (numel(newSize)~=numel(size(imInput)))
        error('Mismatch of new size and image dimension...');
    end
    

    %% Calculation

    % Check for dimension == 1
    [dimEqual1,numDim,multi1] = xASL_Resample_Check1D(newSize);
    
    % Remove "empty" dimension
    [newSize,imInput] = xASL_Resample_Fix1D(dimEqual1,numDim,newSize,imInput,multi1);
    
    if numel(newSize)==1
        
        if newSize(1)==0
            imOutput = imInput;
        else
            % Define the new grid
            [Xq] = ndgrid(1:newSize(1));
            
            % Adapt the new grid to perfectly match the original grip - i.e. stretch or squeeze to the original grid coordinates
            Xq = (Xq-1)/(newSize(1)-1)*(size(imInput,1)-1) + 1;
            imOutput = interp1(imInput, Xq);
        end
        
    elseif numel(newSize)==2
        
        % Define the new grid
        [Xq, Yq] = ndgrid(1:newSize(1), 1:newSize(2));
        
        % Adapt the new grid to perfectly match the original grip - i.e. stretch or squeeze to the original grid coordinates
        Xq = (Xq-1)/(newSize(1)-1)*(size(imInput,1)-1) + 1;
        Yq = (Yq-1)/(newSize(2)-1)*(size(imInput,2)-1) + 1;
        imOutput = interpn(imInput, Xq, Yq);
        
    elseif numel(newSize)==3
        
        % Define the new grid
        [Xq, Yq, Zq] = ndgrid(1:newSize(1), 1:newSize(2), 1:newSize(3));
        
        % Adapt the new grid to perfectly match the original grip - i.e. stretch or squeeze to the original grid coordinates
        Xq = (Xq-1)/(newSize(1)-1)*(size(imInput,1)-1) + 1;
        Yq = (Yq-1)/(newSize(2)-1)*(size(imInput,2)-1) + 1;
        Zq = (Zq-1)/(newSize(3)-1)*(size(imInput,3)-1) + 1;
        imOutput = interpn(imInput, Xq, Yq, Zq);
        
    end
    

end


%% Check if one dimension is equal to 1 (does not work if multiple dimensions are set to 1)
function [dimEqual1,numDim,multi1] = xASL_Resample_Check1D(newSize)

    dimEqual1 = false;
    numDim = 0;
    multi1 = false;
    for iDim = 1:numel(newSize)
        if newSize(iDim)==1
            if ~dimEqual1
                dimEqual1 = true;
            else
                multi1 = true;
            end
            numDim = iDim;
        end
    end

end


%% Remove "empty" dimension
function [newSize,imInput] = xASL_Resample_Fix1D(dimEqual1,numDim,newSize,imInput,multi1)

    % Are multiple dimensions 1?
    if multi1
        if numel(newSize)==2
            % Convert 2D -> 1x1
            newSize = 0;
            imInput = mean(imInput(:)');
            return
        elseif numel(newSize)==3
            % Convert 3D -> 1x1
            if newSize(1)==1 && newSize(2)==1 && newSize(3)==1
                newSize = 0;
                imInput = mean(imInput(:)');
                return
            end
            % Convert 3D -> 1D
            if newSize(1)==1 && newSize(2)==1
                imInput = imInput(1,1,:);
                newSize = newSize(3);
            elseif newSize(1)==1 && newSize(3)==1
                imInput = imInput(1,:,1);
                newSize = newSize(2);
            elseif newSize(2)==1 && newSize(3)==1
                imInput = imInput(:,1,1);
                newSize = newSize(1);
            end
            return
        end
    end

    if dimEqual1
        if numel(newSize)==1
            % Convert a 1D image/signal to a single value -> just return the mean
            newSize = 0;
            imInput = mean(imInput(:)');
        elseif numel(newSize)==2
            % Convert a 2D image to a 1D array
            if numDim==1
                newSize = newSize(2);
            else
                newSize = newSize(1);
            end
        elseif numel(newSize)==3
            % Convert a 3D image to a 2D image
            if numDim==1
                newSize = [newSize(2), newSize(3)];
                imInput = imInput(1,:,:);
            elseif numDim==2
                newSize = [newSize(1), newSize(3)];
                imInput = imInput(:,1,:);
            else
                newSize = [newSize(1), newSize(2)];
                imInput = imInput(:,:,1);
            end
        end
    end

end







