function x = xASL_imp_DCM2NII(x, imPar)
%xASL_imp_DCM2NII Run the dcm2nii part of the import.
%
% FORMAT: x = xASL_imp_DCM2NII(x, imPar)
% 
% INPUT:
%   x      - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar  - JSON file with structure with import parameters (REQUIRED, STRUCT)
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
% EXAMPLE:     x = xASL_imp_DCM2NII(x, imPar);
%
% __________________________________
% Copyright (c) 2015-2021 ExploreASL

    
    %% 1. Initialize defaults of dcm2nii
    
    % Make sure that logging is still active
    diary(x.dir.diaryFile);
    
    % Print feedback
    xASL_adm_BreakString('DICOM to NIFTI CONVERSION');
    
    %% 2. Create the temp directory for DCM2NII
    xASL_adm_CreateDir(imPar.TempRoot);
    
    % Initialization of an empty catched errors struct
    dcm2niiCatchedErrors = struct;

    %% 3. Import visit by visit, session by session, scan by scan for current subject
    fprintf('\nRunning DCM2NIIX...\n');
    [x, imPar, PrintDICOMFields, dcm2niiCatchedErrors] = xASL_imp_DCM2NII_Subject(x, imPar, x.modules.import.matches, dcm2niiCatchedErrors);
    
    % We do not iterate over subjects anymore, since this is done in xASL_Iteration now
    iSubject = strcmp(x.SUBJECT,x.SUBJECTS);
    
    % In case there were illegal characters in the subject ID, we get rid of them within x.SUBJECT & x.SUBJECTS here
    x.SUBJECT = xASL_adm_CorrectName(x.SUBJECTS{iSubject},2);
    x.SUBJECTS{iSubject} = xASL_adm_CorrectName(x.SUBJECTS{iSubject},2);
    
    %% 4. Create summary file
    overviewSubjects = fieldnames(x.overview);
    thisSubject = x.overview.(overviewSubjects{iSubject});
    xASL_imp_CreateSummaryFile(thisSubject, imPar, PrintDICOMFields, x);
    
    %% 5. Clean-Up
    xASL_imp_DCM2NII_CleanUp(x, thisSubject, imPar, dcm2niiCatchedErrors);
    
    
end



%% Clean-Up
function xASL_imp_DCM2NII_CleanUp(x, thisSubject, imPar, dcm2niiCatchedErrors)


    % Clean up DCMTK and DCMDICT
    if ~x.modules.import.settings.bUseDCMTK || isempty(x.modules.import.pathDcmDict)
        dicomdict('factory');
    end
    
    % Make sure diary is off
    diary('off');
    
    % Check catched errors
    if ~isempty(fields(dcm2niiCatchedErrors))
        fclose all;
        SavePath = fullfile(imPar.TempRoot, 'dcm2niiCatchedErrors.mat');
        SaveJSON = fullfile(imPar.TempRoot, 'dcm2niiCatchedErrors.json');
        xASL_delete(SavePath);
        xASL_delete(SaveJSON);
        save(SavePath,'dcm2niiCatchedErrors');
        spm_jsonwrite(SaveJSON, dcm2niiCatchedErrors);
    end
    
    % If zipped sourcedata was imported, we can remove tempSourcedata/sourcedata/[subject]
    tempSourceSubject = fullfile(x.dir.xASLDerivatives,'tempSourcedata','sourcedata',thisSubject.name);
    if xASL_exist(tempSourceSubject,'dir')
        xASL_delete(tempSourceSubject,1);
    end
    
    % If the tempSourcedata directory exists and is empty, we delete it completely
    tempSourcedataFiles = xASL_adm_GetFileList(fullfile(x.dir.xASLDerivatives,'tempSourcedata'),'.','FPListRec');
    if numel(tempSourcedataFiles)<1 && xASL_exist(fullfile(x.dir.xASLDerivatives,'tempSourcedata'),'dir')
        xASL_delete(fullfile(x.dir.xASLDerivatives,'tempSourcedata'),1);
    end
    
    % Finish printing
    fprintf('\n');

end



