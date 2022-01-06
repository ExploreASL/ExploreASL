function xASL_test_Flavors_BIDS2LEGACY(testConfig)
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
            ExploreASL(currentFlavor, 0, [1 0 0 0], 0);
        end
    end

end


