function flavorList = xASL_test_Flavors_Compare(testConfig,folderA,folderB)
%xASL_test_Flavors_Compare Compare rawdata or derivatives with reference
%
% FORMAT: flavorList = xASL_test_Flavors_Compare(testConfig,folderA,folderB)
%
% INPUT:        testConfig - Struct containing test infos (STRUCT, REQUIRED)
%               folderA    - Directory A (STRING, REQUIRED)
%               folderB    - Directory A (STRING, REQUIRED)
%
% OUTPUT:       flavorList - Flavor list
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Compare rawdata or derivatives with reference.
%
% EXAMPLE:      flavorList = xASL_test_Flavors_Compare(testConfig,folderA,folderB);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

    % List all studies in the import directory
    flavorList = testConfig.flavorList;
    
    % Iterate over flavors
    for iCompare = 1:length(flavorList)
        
        % Compare the imported data in the 'rawdata' subdirectory with the counterpart
        fprintf('%s\n', ['Dataset: '  flavorList{iCompare}]);
        convertedData = fullfile(testConfig.pathFlavorDatabase, flavorList{iCompare}, folderA);
        referenceData = fullfile(testConfig.pathFlavorDatabase, flavorList{iCompare}, folderB);
        
        % Run comparison script
        [flavorList{iCompare,2},flavorList{iCompare,3}] = ...
            xASL_bids_CompareStructures(referenceData,convertedData,[],[],0,1);
    end

end
