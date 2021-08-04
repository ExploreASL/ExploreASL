function [LR_flip_YesNo] = xASL_im_DetermineFlip(PathOrientationResults)
%xASL_im_DetermineFlip Detect potential left-right flip
%
% FORMAT: [LR_flip_YesNo] = xASL_im_DetermineFlip(PathOrientationResults)
% 
% INPUT:        
% PathOrientationResults    - path to TSV file containing the orientation
%                             parameters
%
% OUTPUT:
% LR_flip_YesNo             - for a single image:  true if left-right flip is found
%                             for multiple images: indices for images with
%                                                  left-right flip
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This functions check determinants before and after image processing
%               (nii.mat0 vs nii.mat, respectively) to find any potential
%               left-right processing. This function performs the following
%               steps:
%               1. Determine correct row, differs between Matlab versions
%               2. If units are printed as second row, the data starts on the third row
%               3. Determine column indices
%               4. Find left-right flips
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: LR_flip_YesNo = xASL_im_DetermineFlip(PathOrientationResults);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% ============================================================
    %% Admin
    if nargin<1 || isempty(PathOrientationResults)
        error('PathOrientationResults missing');
    elseif ~exist(PathOrientationResults, 'file')
        error([PathOrientationResults ' non-existing, skipping']);
        % PM: this error may be a warning
    end

    LR_flip_YesNo = NaN; % default

    [~, CellTSV] = xASL_bids_csv2tsvReadWrite(PathOrientationResults);

    %% ============================================================    
    %% 1. Determine correct row, differs between Matlab versions
    if size(CellTSV,1)<2 || size(CellTSV,2)<26
        warning(['Missing data in: ' PathOrientationResults]);
        return;
    end
    
    %% ============================================================    
    %% 2. If units are printed as second row, the data starts on the third row
    if isempty(CellTSV{2,13})
        firstRow = 3;
    else
        firstRow = 2;
    end

    %% ============================================================
    %% 3. Determine column indices
    columnDetOri = find(strcmp(CellTSV(1,:), 'DetOrigin'));
    columnDetNew = find(strcmp(CellTSV(1,:), 'DetCurrent'));
    
    if isempty(columnDetOri) || ~isnumeric(columnDetOri)
        warning('Information missing for original determinant, skipping');
        return;
    elseif isempty(columnDetNew) || ~isnumeric(columnDetNew)
        warning('Information missing for new determinant, skipping');
        return;
    end
    
    DeterminantOri = xASL_str2num(CellTSV(firstRow:end, columnDetOri));
    DeterminantNew = xASL_str2num(CellTSV(firstRow:end, columnDetNew));
    
    %% ============================================================
    %% 4. Find left-right flips
    if numel(DeterminantOri)==1
        % standard ExploreASL behavior is to do this per subject in
        % xASL_wrp_VisualQC* -> xASL_qc_CollectQC*        
        LR_flip_YesNo = max((DeterminantOri.*DeterminantNew)<0);
    else
        LR_flip_YesNo = find((DeterminantOri.*DeterminantNew)<0);
    end
    
    
end