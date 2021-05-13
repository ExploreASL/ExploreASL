function imageType = xASL_bids_determineImageTypeGE(jsonPar)
%xASL_bids_determineImageTypeGE Determine the image type of a GE DICOM.
%
% FORMAT: imageType = xASL_bids_determineImageTypeGE(jsonPar)
% 
% INPUT:
%   jsonPar    - Header of DICOM file that is usually stored as a JSON (STRUCT, REQUIRED)
%
% OUTPUT:
%   imageType  - (CHAR ARRAY)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Determine the image type of a GE DICOM.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     iHeader = xASL_io_DcmtkRead(iFile);
%              imageType = xASL_bids_determineImageTypeGE(iHeader);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Starts looking for the correct image type
	imageType = '';
    
    % Check ImageType field
    if isfield(jsonPar, 'ImageType')
        if iscell(jsonPar.ImageType) && length(jsonPar.ImageType)>1
            warning('jsonPar.ImageType had a wrong format of multiple cells');
        elseif iscell(jsonPar.ImageType) && length(jsonPar.ImageType)==0
            warning('jsonPar.ImageType was an empty cell');
        elseif iscell(jsonPar.ImageType)
            jsonPar.ImageType = jsonPar.ImageType{1};
        end
            
        if ~iscell(jsonPar.ImageType)
            jsonPar.ImageType = strsplit(jsonPar.ImageType,'\')';
        end
    else
        return;
    end

    % ["ImageType": ["DERIVED", "PRIMARY", "ASL", "PERFUSION", "ASL"] - deltaM
    if length(jsonPar.ImageType) == 5 && strcmpi(jsonPar.ImageType{1},'DERIVED') && strcmpi(jsonPar.ImageType{2},'PRIMARY') &&...
            strcmpi(jsonPar.ImageType{3},'ASL') && strcmpi(jsonPar.ImageType{4},'PERFUSION') && strcmpi(jsonPar.ImageType{5},'ASL')
        imageType = 'deltam';
    end

    % ["DERIVED", "PRIMARY", "ASL", "PERFUSION_ASL"] - deltaM
    if length(jsonPar.ImageType) == 4 && strcmpi(jsonPar.ImageType{1},'DERIVED') && strcmpi(jsonPar.ImageType{2},'PRIMARY') &&...
            strcmpi(jsonPar.ImageType{3},'ASL') && strcmpi(jsonPar.ImageType{4},'PERFUSION_ASL')
        imageType = 'deltam';
    end

    % ["ORIGINAL", "PRIMARY", "ASL"] - M0
    if length(jsonPar.ImageType) == 3 && strcmpi(jsonPar.ImageType{1},'ORIGINAL') && strcmpi(jsonPar.ImageType{2},'PRIMARY') &&...
            strcmpi(jsonPar.ImageType{3},'ASL')
        imageType = 'm0scan';
    end

    % ["DERIVED", "PRIMARY", "CBF", "CBF"] - CBF
    if length(jsonPar.ImageType) == 4 && strcmpi(jsonPar.ImageType{1},'DERIVED') && strcmpi(jsonPar.ImageType{2},'PRIMARY') &&...
            strcmpi(jsonPar.ImageType{3},'CBF') && strcmpi(jsonPar.ImageType{4},'CBF')
        imageType = 'cbf';
    end

end