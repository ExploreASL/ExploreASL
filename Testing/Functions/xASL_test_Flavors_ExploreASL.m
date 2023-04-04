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
% Copyright (c) 2015-2022 ExploreASL
    
    % Iterate over flavors
    for iList=1:numel(testConfig.flavorList)
        % Get the current flavor
        currentFlavor = fullfile(testConfig.pathFlavorDatabase,testConfig.flavorList{iList});
        pathDerivatives = fullfile(currentFlavor,'derivatives');
        % Process data that were converted to derivatives
        if exist(pathDerivatives,'dir')
            % Determine the ExploreASL directory
            pathExploreASLflavor = fullfile(pathDerivatives,'ExploreASL');
            if exist(pathExploreASLflavor,'dir')
                % Run ExploreASL structural and ASL module (don't run the population module)
                try
                    xFlavor = ExploreASL(currentFlavor, 0, [1 1 0], 0);
                catch ME
                    xFlavor.logging.message = ME.message;
                    xFlavor.logging.name = testConfig.flavorList{iList};
                    xFlavor.logging.stack = ME.stack;
                end 
                
                if isfield(xFlavor,'logging')
                    loggingTable = xASL_test_AddLoggingEntryToTable(testConfig.flavorList{iList},loggingTable,xFlavor.logging);
                end
            end
        end
    end

end


