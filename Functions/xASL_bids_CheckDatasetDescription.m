function [bImportedExploreASL, bImportedSameVersion, versionExploreASLBIDS, bImportedBETA] = xASL_bids_CheckDatasetDescription(datasetDescription, versionExploreASL)
%xASL_bids_CheckDatasetDescription Check the dataset_description.json field
% contents with special regard to the import version
%
% FORMAT: [bImportedExploreASL, bImportedSameVersion, versionExploreASLBIDS, bImportedBETA] = xASL_bids_CheckDatasetDescription(datasetDescription, versionExploreASL)
% 
% INPUT:
%   datasetDescription - Structure which defines the dataset_description fields (STRUCT, REQUIRED)
%   versionExploreASL  - Version of the current ExploreASL (STRING, REQUIRED)
%
% OUTPUT:
%   bImportedExploreASL   - Was the dataset imported with ExploreASL or not (BOOLEAN)
%   bImportedSameVersion  - True if the version matches exactly the current version (full match with BETA) (BOOLEAN)
%   versionExploreASLBIDS - Empty (if above == false), or string with version (excluding the BETA) (STRING)
%   bImportedBETA         - True if BETA version was used -> see above, otherwise false (BOOLEAN)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check the dataset_description.json field contents with special regard to the import version.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      [bImportedExploreASL, bImportedSameVersion, versionExploreASLBIDS, bImportedBETA] = xASL_bids_CheckDatasetDescription(datasetDescription, x.Version);
%               
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Input check
    if nargin<1 || isempty(datasetDescription)
        error('Please provide a struct based on the dataset_description.json BIDS file...');
    end
    if nargin<2 || isempty(versionExploreASL)
        error('Please provide the version of your ExploreASL installation...');
    end
    
    
    %% Defaults
    bImportedExploreASL = false;
    bImportedSameVersion = false;
    versionExploreASLBIDS = [];
    bImportedBETA = false;
    
    
    %% Check structure
    
    % Was the dataset imported with ExploreASL or not
    if isfield(datasetDescription,'Acknowledgements') && ...
       ~isempty(regexp(datasetDescription.Acknowledgements,'Imported with ExploreASL', 'once'))
        bImportedExploreASL = true;
    end
    
    % True if the version matches exactly the current version
    if bImportedExploreASL
        indexVersionString = regexp(datasetDescription.Acknowledgements,'Imported with ExploreASL', 'once')+length('Imported with ExploreASL');
        fileVersion = datasetDescription.Acknowledgements(indexVersionString+1:end);
        if strcmp(fileVersion,versionExploreASL)
            bImportedSameVersion = true;
        end
    end
    
    % Empty (if above == false), or string with version
    if bImportedExploreASL
        indexUnderscore = regexp(fileVersion,'_');
        if ~isempty(indexUnderscore)
            versionExploreASLBIDS = fileVersion(1:indexUnderscore-1);
        else
            versionExploreASLBIDS = fileVersion;
        end
    end
    
    % True if BETA version was used -> see above, otherwise false
    if bImportedExploreASL
        if ~isempty(regexp(fileVersion,'BETA', 'once'))
            bImportedBETA = true;
        end
    end










end

