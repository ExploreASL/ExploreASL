function xASL_test_Flavors_RemoveExistingTestData(testConfig)
%xASL_test_Flavors_RemoveExistingTestData Remove existing test data
%
% FORMAT: xASL_test_Flavors_RemoveExistingTestData(testConfig)
%
% INPUT:        testConfig - Struct containing test infos (REQUIRED)
%
% OUTPUT:       n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Remove existing test data.
%
% EXAMPLE:      xASL_test_Flavors_RemoveExistingTestData(testConfig);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


    fprintf('Remove existing test data...\n');

    % Iterate over the flavors
    for iFlavor = 1:length(testConfig.flavorList)
        % Get current flavor directory
        thisFlavorDir = fullfile(testConfig.pathFlavorDatabase,testConfig.flavorList{iFlavor});

        % Check for temp folder
        if xASL_exist(fullfile(thisFlavorDir,'temp'), 'dir')==7
            xASL_delete(fullfile(thisFlavorDir,'temp'),1);
        end
        
        % Check for rawdata folder
        if xASL_exist(fullfile(thisFlavorDir,'rawdata'), 'dir')==7
            xASL_delete(fullfile(thisFlavorDir,'rawdata'),1);
        end

        % Check for derivatives folder
        if xASL_exist(fullfile(thisFlavorDir,'derivatives'), 'dir')==7
            xASL_delete(fullfile(thisFlavorDir,'derivatives'),1);
        end

        % Check for log and other files
        allFilesInThisDir = xASL_adm_GetFileList(thisFlavorDir,'^.+$');

        % Do not delete sourceStructure.json, studyPar.json and dataPar.json
        for iFile=1:numel(allFilesInThisDir)
            [~, thisFile, thisExtension] = fileparts(allFilesInThisDir{iFile});
            if ~strcmp([thisFile thisExtension],'sourceStructure.json') && ...
                    ~strcmp([thisFile thisExtension],'studyPar.json') && ...
					~strcmp([thisFile thisExtension],'README.md') && ...
                    ~strcmp([thisFile thisExtension],'dataPar.json')
                xASL_delete(allFilesInThisDir{iFile},1);
            end
        end
    end


end


