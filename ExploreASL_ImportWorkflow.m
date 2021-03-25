function [x] = ExploreASL_ImportWorkflow(x)
%ExploreASL_ImportWorkflow Multi-step import workflow for DCM2NII, NII2BIDS & BIDS2LEGACY.
%
% FORMAT: [x] = ExploreASL_ImportWorkflow(x)
%
% INPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Multi-step import workflow for DCM2NII, NII2BIDS & BIDS2LEGACY.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Import Workflow
    
    % We expect x.DataParPath to be the sourceStructure.json file, so the root directory should be the study directory
    [x.StudyRoot,~,~] = xASL_fileparts(x.DataParPath);
    
    % Check if at least one of the three steps should be performed
    if sum(x.ImportArray)>0
        if ~isempty(x.DataParPath)
            if exist(x.DataParPath,'file')
                % DICOM TO NII
                if x.ImportArray(1)==1
                    ExploreASL_ImportBIDS(x.StudyRoot, x.DataParPath, [], [1 0 0], false, true, false, false, x);
                end
                % NII TO BIDS
                if x.ImportArray(2)==1
                    ExploreASL_ImportBIDS(x.StudyRoot, x.DataParPath, [], [0 1 0], false, true, false, false, x);
                end
                % ANONYMIZE
                if x.ImportArray(3)==1
                    ExploreASL_ImportBIDS(x.StudyRoot, x.DataParPath, [], [0 0 1], false, true, false, false, x);
                end
                % BIDS TO LEGACY
                if x.ImportArray(4)==1
                    x = xASL_Import_BIDS2LEGACY(x);
                end
            else
                warning('ImportArray was set to 1, but the sourceStructure file does not exist');
            end
        else
            warning('ImportArray was set to 1, but there is no DataParPath, Import will not be executed');
        end
    end
    
    
end

%% xASL_Import_BIDS2LEGACY
function [x] = xASL_Import_BIDS2LEGACY(x)

    % Go through all studies
    ListFolders = xASL_adm_GetFileList(x.StudyRoot, '^.+$', 'FPListRec', [0 Inf], 1);
    for iList=1:numel(ListFolders)
        % Convert only those containing raw data
        [thisRootFolder,thisFolderName,~] = xASL_fileparts(ListFolders{iList});
        if strcmp(thisFolderName,'rawdata') && exist(ListFolders{iList},'dir')
            
            % Clean the old data
            if isunix
                if exist(fullfile(ListFolders{iList}, 'derivatives'), 'dir')
                    diary('off');
                    fclose('all'); % ensure that no file is locked
                    system(['rm -rf ' fullfile(ListFolders{iList}, 'derivatives')]);
                end
            else
            	if exist(fullfile(ListFolders{iList}, 'derivatives'), 'dir')
	                fprintf('Delete existing derivatives folders...\n');
	                diary('off');
	                fclose('all'); % ensure that no file is locked
	                xASL_delete(fullfile(ListFolders{iList}, 'derivatives'),true);
	            end
            end
            
            % Default dataPar.json for the testing that is fast to run
            dataPar.x = x;

            % Run the legacy conversion: Check if a dataPar is provided, otherwise use the defaults
            fListDataPar = xASL_adm_GetFileList(ListFolders{iList},'(?i)(^dataPar.*\.json$)', 'FPList', [], 0);
            if length(fListDataPar) < 1
                % Fill the dataPars with default parameters
                xASL_bids_BIDS2Legacy(thisRootFolder, 1, dataPar);
            else
                % Fill the dataPars with the provided parameters
                dataPar = spm_jsonread(fListDataPar{1});
                xASL_bids_BIDS2Legacy(thisRootFolder, 1, dataPar);
            end
        end
    end
    
    % Overwrite DataParPath in x structure
    fprintf('Overwriting x.DataParPath...\n');
    ListDataParFiles = xASL_adm_GetFileList(x.StudyRoot, '^DataPar.+$', 'FPListRec', [0 Inf]);
    x.DataParPath = ListDataParFiles{1,1};
    
end




