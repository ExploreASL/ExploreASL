%% Run the BIDS to Legacy conversion
function xASL_test_Flavors_BIDS2LEGACY(testConfig)

    % Default dataPar.json for the testing that is fast to run
    
    % Default regular expression
    defaultDataPar.x.dataset.subjectRegexp = '^sub-.*$';
    
    % When there is no structural data, use ASL-MNI registration
    defaultDataPar.x.modules.asl.bUseMNIasDummyStructural = 1;
    
    % Speed up testing and delete temporary files
    defaultDataPar.x.settings.Quality = 0; 
    defaultDataPar.x.settings.DELETETEMP = 1;

    % Go through all studies
    for iList=1:numel(testConfig.flavorList)
        % Convert only those containing raw data
        if exist(fullfile(testConfig.flavorList{iList},'rawdata'),'dir')

            % Currently, we clean the old data for unix only
            if isunix
                if exist(fullfile(testConfig.flavorList{iList}, 'derivatives'), 'dir')
                    diary('off');
                    fclose('all'); % ensure that no file is locked
                    system(['rm -rf ' fullfile(testConfig.flavorList{iList}, 'derivatives')]);
                end
            else
                % Use xASL_delete on windows
                diary('off');
                fclose('all'); % ensure that no file is locked
                xASL_delete(fullfile(testConfig.flavorList{iList}, 'derivatives'),true);
            end

            % Run the legacy conversion
            % Check if a dataPar is provided, otherwise use the defaults
            fListDataPar = xASL_adm_GetFileList(testConfig.flavorList{iList},'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
            if length(fListDataPar) < 1
                % Fill the dataPars with default parameters
                pathDefaultDataPar = fullfile(testConfig.flavorList{iList},'dataPar.json');
                spm_jsonwrite(pathDefaultDataPar,defaultDataPar);
                ExploreASL(testConfig.flavorList{iList}, [0 0 0 1], 0, 0);
                xASL_delete(pathDefaultDataPar);
            else
                % Fill the dataPars with the provided parameters
                ExploreASL(testConfig.flavorList{iList}, [0 0 0 1], 0, 0);
            end

        end
    end

end


