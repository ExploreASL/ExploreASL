function x = xASL_imp_BasicParameterChecks(x)
%xASL_imp_BasicParameterChecks Basic parameter checks for the import pipeline
%
% FORMAT: x = xASL_imp_BasicParameterChecks(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Basic parameter checks for the import pipeline.
%
%    - check that x.dir.DatasetRoot is correct
%    - determine x.dir.sourceStructure
%    - determine x.dir.studyPar
%    - check x.opts.ImportModules
%    - set the defaults for ...
%       - x.modules.import.settings.bCopySingleDicoms
%       - x.modules.import.settings.bUseDCMTK
%       - x.modules.import.settings.bCheckPermissions
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2021 ExploreASL


    %% Basic checks

    % Fix x.dir.DatasetRoot if the last character is a forward or backward slash
    if strcmp(x.dir.DatasetRoot(end),'\') || strcmp(x.dir.DatasetRoot(end),'/')
        x.dir.DatasetRoot = x.dir.DatasetRoot(1:end-1);
    end

    % Check the imagePar input file
    if isempty(x.dir.sourceStructure) && x.opts.ImportModules(1)
        % If the path is empty, then try to find sourceStructure.json or sourcestruct.json
        fListImPar = xASL_adm_GetFileList(x.dir.DatasetRoot,'(?i)^source(struct(ure|)\.json$', 'List', [], 0);
        if length(fListImPar) < 1
            error('Could not find the sourceStructure.json file...');
        end
        x.dir.sourceStructure = fullfile(x.dir.DatasetRoot,fListImPar{1});
    elseif isempty(x.dir.sourceStructure) && (x.opts.ImportModules(2) || x.opts.Deface)
        % For BIDS2Legacy we do not need the imPar struct (and we should not need it for NII2BIDS or DEFACE)
        x.dir.sourceStructure = [];
    else
        fpath = fileparts(x.dir.sourceStructure);
        if isempty(fpath)
            x.dir.sourceStructure = fullfile(x.dir.DatasetRoot,x.dir.sourceStructure);
        end
    end

    % Find the studyPar input file
    if isempty(x.dir.studyPar)
        % If the path is empty, then try to find studyPar.json
        fListStudyPar = xASL_adm_GetFileList(x.dir.DatasetRoot,'(?i)^studypar\.json$', 'List', [], 0);
        if length(fListStudyPar) < 1
            warning('Could not find the studyPar.json file...');
        else
            x.dir.studyPar = fullfile(x.dir.DatasetRoot,fListStudyPar{1});
        end
    else
        fpath = fileparts(x.dir.studyPar);
        if isempty(fpath)
            x.dir.studyPar = fullfile(x.dir.DatasetRoot,x.dir.studyPar);
        end
    end

    if isempty(x.opts.ImportModules)
        x.opts.ImportModules = [1 1];
    else
        if length(x.opts.ImportModules) ~= 2
            error('x.opts.ImportModules must have length 2...');
        end
    end

    % By default don't copy DICOMs for anonymization reasons
    if isempty(x.modules.import.settings.bCopySingleDicoms)
        x.modules.import.settings.bCopySingleDicoms = false;
    end

    if isempty(x.modules.import.settings.bUseDCMTK)
        % Default set to using DCM-TK
        x.modules.import.settings.bUseDCMTK = true;
    elseif ~x.modules.import.settings.bUseDCMTK && isempty(which('dicomdict'))
        error('Dicomdict missing, image processing probably not installed, try DCMTK instead');
    end

    if isempty(x.modules.import.settings.bCheckPermissions)
        if isunix
            x.modules.import.settings.bCheckPermissions = true;
        else
            x.modules.import.settings.bCheckPermissions = false;
        end
    end




end



