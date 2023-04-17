function [parms, pathDcmDict] = xASL_imp_DCM2NII_Subject_StoreJSON(imPar, SavePathJSON, first_match, bUseDCMTK, pathDcmDict)
%xASL_imp_DCM2NII_Subject_StoreJSON Store JSON.
%
% FORMAT: [parms, pathDcmDict] = xASL_imp_DCM2NII_Subject_StoreJSON(imPar, SavePathJSON, first_match, bUseDCMTK, pathDcmDict)
% 
% INPUT:
%   imPar        - Structure with import parameters (REQUIRED, STRUCT)
%   SavePathJSON - Save path(s) JSON (CELL ARRAY, REQUIRED)
%   first_match  - First match (CHAR ARRAY, PATH, REQUIRED)
%   bUseDCMTK    - Use DCMTK (BOOLEAN, REQUIRED)
%   pathDcmDict  - Path to Dicom dictionary (CHAR ARRAY, PATH, REQUIRED)
%
% OUTPUT:
%   parms       - Parameters cell array
%   pathDcmDict - Path to Dicom dictionary
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Store JSON.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2023 ExploreASL

    %% Store JSON
    if exist(SavePathJSON{1}, 'file') && ~isempty(first_match)
        [~, ~, fext] = fileparts(first_match);
        if  strcmpi(fext,'.PAR')
            if length(SavePathJSON) > 1
                warning('xASL_bids_Par2JSON only works with a single file');
            end
            parms = xASL_bids_Par2JSON(first_match, SavePathJSON{1});
        elseif strcmpi(fext,'.nii')
            parms = [];
        elseif imPar.bMatchDirectories
            Fpath  = fileparts(first_match);
            [parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, Fpath, SavePathJSON, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
            clear Fpath Ffile Fext
        else
            [parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, first_match, SavePathJSON, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
        end
    end

    % Fallback
    if ~exist('parms','var')
        parms = cell(1,1);
    end

end


