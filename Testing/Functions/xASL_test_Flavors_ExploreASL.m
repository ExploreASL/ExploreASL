function loggingTable = xASL_test_Flavors_ExploreASL(testConfig,loggingTable)
%xASL_test_Flavors_ExploreASL Run ExploreASL on all Legacy-converted data
%
% FORMAT: loggingTable = xASL_test_Flavors_ExploreASL(testConfig,loggingTable)
%
% INPUT:        testConfig   - Struct containing test infos (STRUCT, REQUIRED)
%               loggingTable - Logging table (REQUIRED)
%
% OUTPUT:       loggingTable - Logging table
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Run ExploreASL on all Legacy-converted data.
%
% EXAMPLE:      loggingTable = xASL_test_Flavors_ExploreASL(testConfig,loggingTable);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL
    
    % Iterate over flavors
    for iList=1:numel(testConfig.flavorList)
        currentFlavor = fullfile(testConfig.pathFlavorDatabase,testConfig.flavorList{iList});
        pathDerivatives = fullfile(currentFlavor,'derivatives');
        % Process data that were converted to derivatives
        if exist(pathDerivatives,'dir')
            pathExploreASLflavor = fullfile(pathDerivatives,'ExploreASL');
            if exist(pathExploreASLflavor,'dir')
                % Don't run population module
                xFlavor = ExploreASL(currentFlavor, 0, 0, [1 1 0], 0);
                if isfield(xFlavor,'logging')
                    loggingTable = xASL_test_AddLoggingEntryToTable(testConfig.flavorList{iList},loggingTable,xFlavor.logging);
                end
            end
        end
    end

end


