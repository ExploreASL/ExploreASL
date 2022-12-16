function xASL_qc_bids_validator(rootDir, bReference)
    % xASL_qc_bids_validator script to run bids validation over a studies
    %
    % FORMAT: [summaryStruct] = xASL_qc_bids_validator(['Path/To/Root/Dir')
    %
    % INPUT:        rootDir       - Location of the database you want to validate
    %               bReference    - bool to handle if we're checking for rawdataReference (OPTIONAL, default = false)
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION:  Function to run bids-validator over a specific directory
    %
    % EXAMPLE:      xASL_qc_bids_validator('../StudiesRootDir');
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % Copyright 2015-2022 ExploreASL
    
    %% Admin
    if nargin < 1 || isempty(rootDir)
        warning('Input must contain a folder containing studies data');
        return;
    end

    if nargin < 2 || isempty(bReference)
        bReference = false;
    end

    % Dependency check, bids validator needs to be installed
    [sts, ~] = xASL_system('bids-validator --version', false);
    if sts
      warning('Requires bids-validator from https://github.com/bids-standard/bids-validator');
      return
    end

    % Generate all jsons 
    if bReference
        iFolder = fullfile(rootDir, 'rawdataReference');
    else 
        iFolder = fullfile(rootDir, 'rawdata');
    end

    xASL_adm_CreateDir(rootDir, 'derivatives');
    iFile   = fullfile(rootDir, 'derivatives', 'bids-validation.json');
    input   = ['bids-validator ', iFolder , ' --no-color --json > ', iFile];
    xASL_system(input, true);
end