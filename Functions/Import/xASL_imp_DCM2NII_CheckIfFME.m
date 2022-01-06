function [resultJSON, bTimeEncoded, bTimeEncodedFME] = xASL_imp_DCM2NII_CheckIfFME(nii_files, bTimeEncoded, bTimeEncodedFME)
%xASL_imp_DCM2NII_CheckIfFME Check if the current sequence is a FME (Fraunhofer Mevis) time encoded sequence
%
% FORMAT: [resultJSON, bTimeEncoded, bTimeEncodedFME] = xASL_imp_DCM2NII_CheckIfFME(nii_files, bTimeEncoded, bTimeEncodedFME)
% 
% INPUT:
%  nii_files          - List of nifti files (REQUIRED)
%  bTimeEncoded       - Time encoded sequence or not (BOOLEAN, REQUIRED)
%  bTimeEncodedFME    - Time encoded FME sequence or not (BOOLEAN,  REQUIRED)
%
% OUTPUT:
%   resultJSON         - JSON structure (STRUCT)
%   bTimeEncoded       - Time encoded sequence or not (BOOLEAN)
%   bTimeEncodedFME    - Time encoded FME sequence or not (BOOLEAN)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check if the current sequence is a FME (Fraunhofer Mevis) time encoded sequence.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Initialize JSON struct
    resultJSON = struct;

    if numel(nii_files)>=1
        [resultPath, resultFile] = xASL_fileparts(nii_files{1});
        % Check if we have the corresponding JSON file
        if exist(fullfile(resultPath, [resultFile '.json']), 'file')
            % Load the JSON
            resultJSON = spm_jsonread(fullfile(resultPath, [resultFile '.json']));
            % Check if we have the SeriesDescription field
            if isfield(resultJSON,'SeriesDescription')
                % Determine if we have the specific FME Hadamard sequence from Bremen
                bTimeEncodedFME = ~isempty(regexp(resultJSON.SeriesDescription,'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
                % If the FME sequence was detected we can always set the general bTimeEncoded to true as well
                if bTimeEncodedFME
                    bTimeEncoded = true;
                end
            end
        end
    end

end


