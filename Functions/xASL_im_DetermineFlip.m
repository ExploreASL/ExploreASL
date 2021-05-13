function [LR_flip_YesNo] = xASL_im_DetermineFlip(x, iSubject, PathOrientationResults)
%xASL_im_DetermineFlip Check determinants, should be the same
% before & after registration, otherwise a left-right flip is applied
% This is not visible, but detrimental for image analysis/stats
%
% FORMAT:       [LR_flip_YesNo] = xASL_im_DetermineFlip(x, iS, PathOrientationResults)
% 
% INPUT:        x                      - ExploreASL x structure (STRUCT, REQUIRED)
%               iS                     - Subject number (INTEGER, REQUIRED)
%               PathOrientationResults - Path orientation results (CHAR ARRAY, REQUIRED)
%               QCstruct               - QC struct (STRUCT, REQUIRED)
%
% OUTPUT:       QCstruct               - QC struct (STRUCT)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check determinants, should be the same
%               before & after registration, otherwise a left-right flip is applied
%               This is not visible, but detrimental for image analysis/stats.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2021 ExploreASL

if nargin<3 || isempty(PathOrientationResults)
    error('PathOrientationResults missing');
end
if nargin<2 || isempty(iSubject)
    bSingleNifti = false;
else
    bSingleNifti = true;
end

    LR_flip_YesNo = NaN; % default

    if exist(PathOrientationResults,'file')
        [~, CellTSV] = xASL_bids_csv2tsvReadWrite(PathOrientationResults);

        % Determine correct row, differs between Matlab versions
        
        if size(CellTSV,1)<3 || size(CellTSV,2)<26
                warning(['Missing data in: ' PathOrientationResults]);
                return;
        elseif ~isempty(CellTSV{2,13})
                iRow    = 2;
        elseif ~isempty(CellTSV{3,13})
                iRow    = 3;
        end
    else
        warning(['File missing: ' PathOrientationResults]);
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