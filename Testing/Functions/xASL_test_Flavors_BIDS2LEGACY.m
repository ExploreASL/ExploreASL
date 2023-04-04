function loggingTable = xASL_test_Flavors_BIDS2LEGACY(testConfig, loggingTable)
%xASL_test_Flavors_BIDS2LEGACY Test the BIDS to Legacy conversion
%
% FORMAT: xASL_test_Flavors_BIDS2LEGACY(testConfig)
%
% INPUT:        testConfig - Struct containing test infos (STRUCT, REQUIRED)
%
% OUTPUT:       n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Run the BIDS to Legacy conversion.
%
% EXAMPLE:      xASL_test_Flavors_BIDS2LEGACY(testConfig);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

    % Go through all studies
    for iList=1:numel(testConfig.flavorList)
        currentFlavor = fullfile(testConfig.pathFlavorDatabase,testConfig.flavorList{iList});
        % Convert only those containing raw data
        if exist(fullfile(currentFlavor,'rawdata'),'dir')
            % ensure that no file is locked
            diary('off');
            fclose('all');
            % Run the legacy conversion
            try
                xFlavor = ExploreASL(currentFlavor);
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


