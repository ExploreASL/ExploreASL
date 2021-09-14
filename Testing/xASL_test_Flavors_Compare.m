%% Compare conversion
function flavorList = xASL_test_Flavors_Compare(testConfig,folderA,folderB)

    % List all studies in the import directory
    flavorList = xASL_adm_GetFileList(testConfig.pathFlavorDatabase, '^.+$', 'List', [], true);
    
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


