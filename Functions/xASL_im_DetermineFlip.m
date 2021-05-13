function [LR_flip_YesNo] = xASL_im_DetermineFlip(x, iSubject, PathOrientationResults)
%xASL_im_DetermineFlip Check determinants, should be the same
% before & after registration, otherwise a left-right flip is applied
% This is not visible, but detrimental for image analysis/stats
%
% FORMAT:       [LR_flip_YesNo] = xASL_im_DetermineFlip(x, iS, PathOrientationResults)
% 
% INPUT:        iSubject -> OPTIONAL, DEFAULT = check all
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check determinants, should be the same
%               before & after registration, otherwise a left-right flip is applied
%               This is not visible, but detrimental for image analysis/stats.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL


if nargin<3 || isempty(PathOrientationResults)
    error('PathOrientationResults missing');
end
if nargin<2 || isempty(iSubject)
    bSingleNifti = false;
else
    bSingleNifti = true;
end

    LR_flip_YesNo = NaN; % default

    if ~exist(PathOrientationResults, 'file')
        warning(['File missing: ' PathOrientationResults]);
        return;
    end        
        
    [~, CellTSV] = xASL_bids_csv2tsvReadWrite(PathOrientationResults);

    % Determine correct row, differs between Matlab versions
    if size(CellTSV,1)<3 || size(CellTSV,2)<26
            warning(['Missing data in: ' PathOrientationResults]);
            return;
    end
    
    % if units are printed as second row, the data starts on the third row
    if isempty(CellTSV{2,13})
        firstRow = 3;
    else
        firstRow = 2;
    end

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
    
    if bSingleNifti
        % standard ExploreASL behavior is to do this per subject in
        % xASL_wrp_VisualQC* -> xASL_qc_CollectQC*        
        LR_flip_YesNo = max((DeterminantOri.*DeterminantNew)<0);
        
        if LR_flip_YesNo>0
            fprintf(['LR flip found for ' x.SUBJECTS{iSubject}]);
            LR_flip_YesNo = 1;
        end        
    else
        LR_flip_YesNo = find((DeterminantOri.*DeterminantNew)<0);
    end
    
end