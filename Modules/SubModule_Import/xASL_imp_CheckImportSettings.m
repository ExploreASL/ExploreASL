function [x, imPar, fid_summary] = xASL_imp_CheckImportSettings(x, imPar)
%xASL_imp_CheckImportSettings Basic import checks before execution
%
% FORMAT: [x, imPar, fid_summary] = xASL_imp_CheckImportSettings(x, imPar)
%
% INPUT:
%   x            - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   imPar        - JSON file with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   imPar        - JSON file with structure with import parameters
%   fid_summary  - Variable to catch errors and close if valid
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Basic import checks before execution.
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
    if ~isfield(imPar,'dcm2nii_version') || isempty(imPar.dcm2nii_version)
        % OR for PARREC imPar.dcm2nii_version = '20101105'; THIS IS AUTOMATED BELOW
        imPar.dcm2nii_version = '20190902';
    end
    if ~isfield(imPar,'dcmExtFilter') || isempty(imPar.dcmExtFilter)
        % dcmExtFilter: the last one is because some convertors save files without extension, 
        % but there would be a dot/period before a bunch of numbers
        imPar.dcmExtFilter = '^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$';
    end
    
    % Check SkipSubjectIfExists
    if ~isfield(imPar,'SkipSubjectIfExists') || isempty(imPar.SkipSubjectIfExists)
        % allows to skip existing subject folders in the temp folder, when this is set to true,
        % avoiding partly re-importing/converting dcm2niiX when processing has been partly done
        imPar.SkipSubjectIfExists = false;
    else
        warning('Skipping existing subjects in temp folder...');
        fprintf('If you want to overwrite, first remove the full subject folder...');
    end

    % Initialize to be able to catch errors and close if valid
    fid_summary = -1;

end



