function [x] = xASL_imp_CheckImportSettings(x)
%xASL_imp_CheckImportSettings Basic import checks before execution
%
% FORMAT: [x] = xASL_imp_CheckImportSettings(x)
%
% INPUT:
%   x            - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Basic import checks before execution.
%
% - Check permissions for DCM2NII
% - Get correct DCMNII version
% - Define DCM Extension Filter
% - Set default for skip subjects option
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    % Check bCheckPermissions
    if x.modules.import.settings.bCheckPermissions
        dcm2niiDir = fullfile(x.opts.MyPath, 'External', 'MRIcron');
        xASL_adm_CheckPermissions(dcm2niiDir, true); % dcm2nii needs to be executable
    end
    if ~isfield(x.modules.import.imPar,'dcm2nii_version') || isempty(x.modules.import.imPar.dcm2nii_version)
        % OR for PARREC x.modules.import.imPar.dcm2nii_version = '20101105'; THIS IS AUTOMATED BELOW
        x.modules.import.imPar.dcm2nii_version = '20190902';
    end
    if ~isfield(x.modules.import.imPar,'dcmExtFilter') || isempty(x.modules.import.imPar.dcmExtFilter)
        % dcmExtFilter: the last one is because some convertors save files without extension, 
        % but there would be a dot/period before a bunch of numbers
        x.modules.import.imPar.dcmExtFilter = '^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$';
    end
    
    % Check SkipSubjectIfExists
    if ~isfield(x.modules.import.imPar,'SkipSubjectIfExists') || isempty(x.modules.import.imPar.SkipSubjectIfExists)
        % allows to skip existing subject folders in the temp folder, when this is set to true,
        % avoiding partly re-importing/converting dcm2niiX when processing has been partly done
        x.modules.import.imPar.SkipSubjectIfExists = false;
    else
        warning('Skipping existing subjects in temp folder...');
        fprintf('If you want to overwrite, first remove the full subject folder...');
    end

end



