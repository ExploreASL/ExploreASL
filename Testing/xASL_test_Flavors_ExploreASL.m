%% Run ExploreASL on all Legacy-converted data
function loggingTable = xASL_test_Flavors_ExploreASL(testConfig,loggingTable)
    
    % Iterate over flavors
    for iList=1:numel(testConfig.flavorList)
        pathDerivatives = fullfile(testConfig.flavorList{iList},'derivatives');
        % Process data that were converted to derivatives
        if exist(pathDerivatives,'dir')
            pathDerivatives = fullfile(pathDerivatives,'ExploreASL');
            if exist(pathDerivatives,'dir')
                % Don't run population module
                xFlavor = ExploreASL(testConfig.flavorList{iList}, 0, [1 1 0], 0);
                if isfield(xFlavor,'logging')
                    loggingTable = xASL_test_AddLoggingEntryToTable(testConfig.flavorList{iList},loggingTable,xFlavor.logging);
                end
            end
        end
    end

end


