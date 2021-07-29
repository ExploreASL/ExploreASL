function jsonIn = xASL_bids_MergeStudyPar(jsonIn,studyPar,bidsModality)
%xASL_bids_MergeStudyPar Check if required fields exist in studyPar but not in jsonIn 
% or if we can find them out in other ways
%
% FORMAT: jsonIn = xASL_bids_MergeStudyPar(jsonIn,studyPar,bidsModality);
%
% INPUT:
%   jsonIn       - JSON with the input fields - from DICOMs (REQUIRED)
%   studyPar     - Manually defined parameters (REQUIRED)
%   bidsModality - Modality (CHAR ARRAY, OPTIONAL, DEFAULT = 'asl')
%
% OUTPUT: 
%   jsonIn    - ordered and checked JSON structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check if required fields exist in studyPar but not in jsonIn or if we can find them in other ways.
%              
%              The BIDSification of JSON metadata requires at least some basic fields. If dcm2niix can't extract
%              fields like Manufacturer from the DICOM data (strict anonymization), we need to be able to read them
%              from the studyPar JSON (manually inserted). Alternatively we can check other DICOM tags for information.
%              
%              This function can potentially be enhanced in future release to fix other fields besides the Manufacturer as well.
%              To enable this functionality for different modalities, we introduced the bidsModality parameter.
%              
%              This function is called by:
%              
%              - `xASL_bids_BIDSifyM0`
%              - `xASL_bids_BIDSifyASLJSON`
%              - `xASL_bids_BIDSifyAnatJSON`
%              
% EXAMPLE:     jsonIn = xASL_bids_MergeStudyPar(jsonIn,studyPar,'asl');
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Default
    if nargin < 3
        bidsModality = 'asl';
    end

    % List of required fields
    switch bidsModality
        case 'asl'
            fieldList = {'Manufacturer', 'TotalAcquiredPairs', 'FlipAngle'};
        case 'm0'
            fieldList = {'Manufacturer', 'FlipAngle'};
        otherwise
            fieldList = {'Manufacturer'};
    end
    
    % Iterate over fields
    for iField = 1:numel(fieldList)
        if isfield(studyPar,fieldList{iField}) && ~isfield(jsonIn,fieldList{iField})
            jsonIn.(fieldList{iField}) = studyPar.(fieldList{iField});
        end
    end

    % Check if manufacturer field does not exist, but we know the model name
    if ~isfield(jsonIn,'Manufacturer') && isfield(jsonIn,'ManufacturersModelName')
        modelList = xASL_bids_BIDSifyFixBasicFields_GetModelList();
        for iModel = 1:size(modelList,1)
            % Determine Manufacturer based on ManufacturersModelName
            if ~isempty(regexpi(jsonIn.ManufacturersModelName,modelList{iModel,1}))
                jsonIn.Manufacturer = modelList{iModel,2};
            end
        end
    end

end


%% ManufacturerModel list: if the Manufacturer DICOM tag was deleted (strict anonymization), we can try to check the ManufacturersModelName tag instead.
function modelList = xASL_bids_BIDSifyFixBasicFields_GetModelList()

    modelList = {
        'Signa',             'GE'; ...
        'Discovery',         'GE'; ...
        'Achieva',           'Philips'; ...
        'Ingenuity',         'Philips'; ...
        'Intera',            'Philips'; ...
        'Ingenia',           'Philips'; ...
        'Skyra',             'Siemens'; ...
        'TrioTim',           'Siemens'; ...
        'Verio',             'Siemens'; ...
        'Prisma',            'Siemens'
    };


end

