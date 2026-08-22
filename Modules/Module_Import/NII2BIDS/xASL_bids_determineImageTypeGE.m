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
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



    %% Starts looking for the correct image type
	imageType = '';
    
    % Check ImageType field
    if isfield(jsonPar, 'ImageType')
        if iscell(jsonPar.ImageType) && length(jsonPar.ImageType)==0
            warning('jsonPar.ImageType was an empty cell');
        elseif iscell(jsonPar.ImageType) && length(jsonPar.ImageType) == 1
            jsonPar.ImageType = jsonPar.ImageType{1};
        end
            
        if ~iscell(jsonPar.ImageType)
            jsonPar.ImageType = strsplit(jsonPar.ImageType,'\')';
        end
    else
        return;
    end

    % ["ImageType": ["DERIVED", "PRIMARY", "ASL", "PERFUSION", "ASL"] - deltaM
    if length(jsonPar.ImageType) == 5 && ~isempty(regexpi(jsonPar.ImageType{1},'DERIVED','once')) && ~isempty(regexpi(jsonPar.ImageType{2},'PRIMARY','once')) &&...
            ~isempty(regexpi(jsonPar.ImageType{3},'^ASL','once')) && ~isempty(regexpi(jsonPar.ImageType{4},'PERFUSION','once')) && ~isempty(regexpi(jsonPar.ImageType{5},'^ASL','once'))
        imageType = 'deltam';
	end
	
	% ["ImageType": ["DERIVED", "PRIMARY", "ASL", "PERFUSION", "ASL", "REAL"] - deltaM
    if length(jsonPar.ImageType) == 6 && ~isempty(regexpi(jsonPar.ImageType{1},'DERIVED','once')) && ~isempty(regexpi(jsonPar.ImageType{2},'PRIMARY','once')) &&...
            ~isempty(regexpi(jsonPar.ImageType{3},'^ASL','once')) && ~isempty(regexpi(jsonPar.ImageType{4},'PERFUSION','once')) &&...
			~isempty(regexpi(jsonPar.ImageType{5},'^ASL','once')) && ~isempty(regexpi(jsonPar.ImageType{6},'REAL','once'))
        imageType = 'deltam';
    end

    % ["DERIVED", "PRIMARY", "ASL", "PERFUSION_ASL"] - deltaM
    if length(jsonPar.ImageType) == 4 && ~isempty(regexpi(jsonPar.ImageType{1},'DERIVED','once')) && ~isempty(regexpi(jsonPar.ImageType{2},'PRIMARY','once')) &&...
            ~isempty(regexpi(jsonPar.ImageType{3},'^ASL','once')) && ~isempty(regexpi(jsonPar.ImageType{4},'PERFUSION_ASL','once'))
        imageType = 'deltam';
	end

	% ["ORIGINAL", "PRIMARY", "ASL", "REAL"] - M0
    if length(jsonPar.ImageType) == 4 && ~isempty(regexpi(jsonPar.ImageType{1},'ORIGINAL','once')) && ~isempty(regexpi(jsonPar.ImageType{2},'PRIMARY','once')) &&...
            ~isempty(regexpi(jsonPar.ImageType{3},'^ASL','once')) && ~isempty(regexpi(jsonPar.ImageType{4},'REAL','once'))
        imageType = 'm0scan';
	end
	
    % ["ORIGINAL", "PRIMARY", "ASL"] - M0
    if length(jsonPar.ImageType) == 3 && ~isempty(regexpi(jsonPar.ImageType{1},'ORIGINAL','once')) && ~isempty(regexpi(jsonPar.ImageType{2},'PRIMARY','once')) &&...
            ~isempty(regexpi(jsonPar.ImageType{3},'^ASL','once'))
        imageType = 'm0scan';
    end

    % ["DERIVED", "PRIMARY", "CBF", "CBF"] - CBF
    if length(jsonPar.ImageType) == 4 && ~isempty(regexpi(jsonPar.ImageType{1},'DERIVED','once')) && ~isempty(regexpi(jsonPar.ImageType{2},'PRIMARY','once')) &&...
            ~isempty(regexpi(jsonPar.ImageType{3},'CBF','once')) && ~isempty(regexpi(jsonPar.ImageType{4},'CBF','once'))
        imageType = 'cbf';
	end

	% ["DERIVED", "PRIMARY", "CBF", "CBF", "REAL"] - CBF
    if length(jsonPar.ImageType) == 5 && ~isempty(regexpi(jsonPar.ImageType{1},'DERIVED','once')) && ~isempty(regexpi(jsonPar.ImageType{2},'PRIMARY','once')) &&...
            ~isempty(regexpi(jsonPar.ImageType{3},'CBF','once')) && ~isempty(regexpi(jsonPar.ImageType{4},'CBF','once')) && ~isempty(regexpi(jsonPar.ImageType{5},'REAL','once'))
        imageType = 'cbf';
    end

end