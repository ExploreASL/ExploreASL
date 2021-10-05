function [x,imPar] = xASL_imp_CheckDirectoriesAndPermissions(x,imPar)
%xASL_imp_CheckDirectoriesAndPermissions Check directories and permissions
%
% FORMAT: [x,imPar] = xASL_imp_CheckDirectoriesAndPermissions(x,imPar)
%
% INPUT:
%   x            - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   imPar        - JSON file with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   imPar        - JSON file with structure with import parameters
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Check directories and permissions.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    % Create the basic folder structure for sourcedata & derivative data
    if ~exist(imPar.RawRoot, 'dir')
        warning(['Could not find ' imPar.RawRoot ', trying to find a different folder instead...']);
        
        % find any folder except for temp, sourcedata, rawdata, derivatives
        % xASL_adm_GetFileList uses regular expressions, to create a nice list of foldernames,
        % with/without FullPath (FPList), with/without recursive (FPListRec)
        % very powerful once you know how these work
        FolderNames = xASL_adm_GetFileList(fullfile(imPar.RawRoot, imPar.studyID), ...
            '^(?!(temp|derivatives|source|sourcedata)).*$', 'FPList', [0 Inf], true);
        
        if length(FolderNames)==1
            imPar.RawRoot = FolderNames{1};
            fprintf('%s\n', ['Found ' imPar.RawRoot ' as sourcedata folder']);
        else
            error('Couldnt find a sourcedata folder, please rename one, or move other folders');
        end
    end
    
    % Check access rights of temp and rawdata directories
    if x.modules.import.settings.bCheckPermissions
        xASL_adm_CheckPermissions(imPar.RawRoot, false); % don't need execution permisions
        xASL_adm_CheckPermissions(imPar.TempRoot, false);  % don't need execution permisions
    end

    % Path to the dictionary to initialize - we need to keep track if the dictionary has been set, 
    % because Dicominfo can be used despite bUSEDCMTK==1 when DCMTK fails
    x.modules.import.pathDcmDict = fullfile(x.opts.MyPath,'External','xASL_DICOMLibrary.txt');
    if ~x.modules.import.settings.bUseDCMTK
        % Initialize dicom dictionary by appending private philips stuff to a temporary copy
        dicomdict('set', x.modules.import.pathDcmDict);
    end


end



