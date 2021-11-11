function [x] = xASL_imp_CheckDirectoriesAndPermissions(x)
%xASL_imp_CheckDirectoriesAndPermissions Check directories and permissions
%
% FORMAT: [x] = xASL_imp_CheckDirectoriesAndPermissions(x)
%
% INPUT:
%   x            - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Check directories and permissions.
%
% - Check if the RawRoot exists
%     - Search for temp, derivatives, source and sourcedata
%     - Raise error if not a single directory exists
% - Check the access rights for temp and rawdata
% - DCMTK & DicomInfo realted settings
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    % Create the basic folder structure for sourcedata & derivative data
    if x.opts.ImportModules(1)
        if ~exist(x.modules.import.imPar.RawRoot, 'dir')
            warning(['Could not find ' x.modules.import.imPar.RawRoot ', trying to find a different folder instead...']);

            % Find any folder except for temp, sourcedata, rawdata, derivatives
            % xASL_adm_GetFileList uses regular expressions, to create a nice list of foldernames,
            % with/without FullPath (FPList), with/without recursive (FPListRec)
            % very powerful once you know how these work
            FolderNames = xASL_adm_GetFileList(fullfile(x.modules.import.imPar.RawRoot, x.modules.import.imPar.studyID), ...
                '^(?!(temp|derivatives|source|sourcedata)).*$', 'FPList', [0 Inf], true);

            if length(FolderNames)==1
                x.modules.import.imPar.RawRoot = FolderNames{1};
                fprintf('%s\n', ['Found ' x.modules.import.imPar.RawRoot ' as sourcedata folder']);
            else
                error('Couldnt find a sourcedata folder, please rename one, or move other folders...');
            end
        end
    end
    
    % Check access rights of temp and rawdata directories
    if x.modules.import.settings.bCheckPermissions
        if x.opts.ImportModules(1)
            % don't need execution permisions
            xASL_adm_CheckPermissions(x.modules.import.imPar.RawRoot, false);
        end
        if x.opts.ImportModules(2)
            % don't need execution permisions
            xASL_adm_CheckPermissions(x.modules.import.imPar.TempRoot, false);
        end
    end

    % Path to the dictionary to initialize - we need to keep track if the dictionary has been set, 
    % because Dicominfo can be used despite bUSEDCMTK==1 when DCMTK fails
    if x.opts.ImportModules(1) || x.opts.ImportModules(2)
        x.modules.import.pathDcmDict = fullfile(x.opts.MyPath,'External','xASL_DICOMLibrary.txt');
        if ~x.modules.import.settings.bUseDCMTK
            % Initialize dicom dictionary by appending private philips stuff to a temporary copy
            dicomdict('set', x.modules.import.pathDcmDict);
        end
    end


end



