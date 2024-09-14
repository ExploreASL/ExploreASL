function x = xASL_wrp_DCM2NII(x)
%xASL_wrp_DCM2NII Run the dcm2nii part of the import.
%
% FORMAT: x = xASL_wrp_DCM2NII(x)
% 
% INPUT:
%   x                      - ExploreASL x structure (REQUIRED, STRUCT)
%   x.modules.import.imPar - A substructure with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   x      - ExploreASL x structure
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run the dcm2nii part of the import.
%
% 1. Initialize defaults of dcm2nii
% 2. Create the temp directory for DCM2NII
% 3. Import visit by visit, session by session, scan by scan for current subject
% 4. Create summary file
% 5. Clean-up
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     x = xASL_wrp_DCM2NII(x);
%
% __________________________________
% Copyright (c) 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    
    %% 1. Initialize defaults of dcm2nii
    
    % Make sure that logging is still active
    diary(x.dir.diaryFile);
    
    % Print feedback
    xASL_adm_BreakString('DICOM to NIFTI CONVERSION');
    
    %% 2. Create the temp directory for DCM2NII
    xASL_adm_CreateDir(x.modules.import.imPar.TempRoot);
    
    % Initialization of an empty catched errors struct
    dcm2niiCatchedErrors = struct;
    %% 3. Import visit by visit, session by session, scan by scan for current subject
    fprintf('\nRunning DCM2NIIX...\n');
    [x, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_wrp_DCM2NII_Subject(x, x.modules.import.matches, dcm2niiCatchedErrors);
    
    % We do not iterate over subjects anymore, since this is done in xASL_Iteration now
    iSubject = strcmp(x.SUBJECT,x.SUBJECTS);
    
    % In case there were illegal characters in the subject ID, we get rid of them within x.SUBJECT & x.SUBJECTS here
    x.SUBJECT = xASL_adm_CorrectName(x.SUBJECTS{iSubject},2);
    x.SUBJECTS{iSubject} = xASL_adm_CorrectName(x.SUBJECTS{iSubject},2);
    
    %% 4. Create summary file
    overviewSubjects = fieldnames(x.importOverview);
    thisSubject = x.importOverview.(overviewSubjects{iSubject});
    xASL_imp_CreateSummaryFile(thisSubject, PrintDICOMFields, x);
    
    %% 5. Clean-Up
    xASL_imp_DCM2NII_CleanUp(x, dcm2niiCatchedErrors);
end
%% Clean-Up
function xASL_imp_DCM2NII_CleanUp(x, dcm2niiCatchedErrors)
    if ~x.modules.import.settings.bUseDCMTK || isempty(x.modules.import.pathDcmDict)
        dicomdict('factory');
    end
    diary('off');
    
    if ~isempty(fields(dcm2niiCatchedErrors))
        fclose all;
        SavePath = fullfile(x.modules.import.imPar.TempRoot, 'dcm2niiCatchedErrors.mat');
        SaveJSON = fullfile(x.modules.import.imPar.TempRoot, 'dcm2niiCatchedErrors.json');
        xASL_delete(SavePath);
        xASL_delete(SaveJSON);
        save(SavePath,'dcm2niiCatchedErrors');
        xASL_io_WriteJson(SaveJSON, dcm2niiCatchedErrors);
    end
    
    fprintf('\n');
end
