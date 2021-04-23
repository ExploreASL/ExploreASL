function xASL_imp_ANONYMIZE(imPar)
%xASL_imp_ANONYMIZE Run defacing.
%
% FORMAT: xASL_imp_ANONYMIZE(imPar)
% 
% INPUT:
%   imPar      - JSON file with structure with import parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run defacing.
%
% 1. Iterate over list of subjects
% 2. Get subject labels
% 3. Process all anatomical files (`xASL_spm_deface`)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_imp_ANONYMIZE(imPar);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Anonymize
    
    % 1. Iterate over list of subjects
    listSubjects = xASL_adm_GetFileList(imPar.AnalysisRoot,[],false,[],true);
    for iSubject = 1:length(listSubjects)

        % 2. Get subject labels
        subjectLabel = xASL_adm_CorrectName(listSubjects{iSubject},2);

        % Check if the anatomical directory exists
        if exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'),'dir')
            % 3. Process all anatomical files
            fAnat = xASL_adm_GetFileList(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'),'^.+\.nii',false,[]);
            for iAnat = 1:length(fAnat)
                %Unzip the file for SPM
                pathUnzipped = xASL_adm_UnzipNifti(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat',fAnat{iAnat}));
                % Remove the face
                xASL_spm_deface(pathUnzipped,true);
                % Zip again
                gzip(pathUnzipped);
                delete(pathUnzipped);
            end
        end
    end

end



