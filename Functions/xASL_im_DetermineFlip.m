function [QCstruct] = xASL_im_DetermineFlip(x,iS,PathOrientationResults,QCstruct)
%xASL_im_DetermineFlip Check determinants, should be the same
% before & after registration, otherwise a left-right flip is applied
% This is not visible, but detrimental for image analysis/stats
%
% FORMAT:       [QCstruct] = xASL_im_DetermineFlip(x,iS,PathOrientationResults,QCstruct)
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

    QCstruct.LR_flip_YesNo = NaN; % default

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

    Det_Ori = xASL_str2num(CellTSV{iRow,13});
    Det_New = xASL_str2num(CellTSV{iRow,26});
    QCstruct.LR_flip_YesNo = uint8(~min((Det_Ori.*Det_New)>0));

    if QCstruct.LR_flip_YesNo>0
        fprintf(['LR flip found for ' x.SUBJECTS{iS}]);
        QCstruct.LR_flip_YesNo = 1;
    end
    
end