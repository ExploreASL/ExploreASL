function xASL_imp_DCM2NII(x, imPar)
%xASL_imp_DCM2NII Run the dcm2nii part of the import.
%
% FORMAT: xASL_imp_DCM2NII(imPar, x)
% 
% INPUT:
%   x                  - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar              - JSON file with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run the dcm2nii part of the import.
%
% - Initialize defaults of dcm2nii
% - Create the basic folder structure for sourcedata & derivative data
% - Here we try to fix backwards compatibility, but this may break
% - Redirect output to a log file
% - Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
% - Sanity check for missing elements
% - Import subject by subject, visit by visit, session by session, scan by scan
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_DCM2NII(imPar, x);
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    
    %% Initialize defaults of dcm2nii
    
    % We may need to restart the logging
    diary(fullfile(x.dir.DatasetRoot,'xASL_module_Import.log'));
    
    % Print feedback
    fprintf('================================== DICOM to NIFTI CONVERSION =================================\n');
    
    % Print matching files
    if isfield(imPar,'bVerbose') && imPar.bVerbose
        fprintf('\nMatching files (#=%g):\n',length(x.modules.import.matches));
        for iMatch=1:size(x.modules.import.matches,1)
            fprintf('%s\n', x.modules.import.matches{iMatch,1});
        end
    end
    
    % Create the temp directory for DCM2NII
    xASL_adm_CreateDir(imPar.TempRoot);
    
    % Initialize to be able to catch errors and close if valid
    fid_summary = -1;
    
    %% Import subject by subject, visit by visit, session by session, scan by scan
    fprintf('\nRunning DCM2NIIX...\n');

    % Initialization of an empty catched errors struct
    dcm2niiCatchedErrors = struct;
    
    % Iterate over subjects
    for iSubject=1:x.modules.import.numOf.nSubjects
        [x, imPar, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_imp_DCM2NII_Subject(x, imPar, iSubject, x.modules.import.matches, dcm2niiCatchedErrors);
    end
    
    % Create summary file
    xASL_imp_CreateSummaryFile(imPar, PrintDICOMFields, x, fid_summary);
    
    % Clean-Up
    xASL_imp_DCM2NII_CleanUp(x, imPar, dcm2niiCatchedErrors);
    
    
end



%% Clean-Up
function xASL_imp_DCM2NII_CleanUp(x,imPar,dcm2niiCatchedErrors)


    if ~x.modules.import.settings.bUseDCMTK || isempty(x.modules.import.pathDcmDict)
        dicomdict('factory');
    end
    diary('off');
    
    if ~isempty(fields(dcm2niiCatchedErrors))
        fclose all;
        SavePath = fullfile(imPar.TempRoot, 'dcm2niiCatchedErrors.mat');
        SaveJSON = fullfile(imPar.TempRoot, 'dcm2niiCatchedErrors.json');
        xASL_delete(SavePath);
        xASL_delete(SaveJSON);
        save(SavePath,'dcm2niiCatchedErrors');
        spm_jsonwrite(SaveJSON, dcm2niiCatchedErrors);
    end
    
    fprintf('\n');


end



