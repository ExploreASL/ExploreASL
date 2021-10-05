function xASL_imp_Deface(x,imPar)
%xASL_imp_Deface Run defacing.
%
% FORMAT: xASL_imp_Deface(imPar)
% 
% INPUT:
%   x          - ExploreASL x structure (REQUIRED, STRUCT)
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
% EXAMPLE:     xASL_imp_Deface(imPar);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% We may need to restart the logging
    diary(fullfile(x.dir.DatasetRoot,'xASL_module_Import.log'));

    %% 1. Iterate over list of subjects
    listSubjects = xASL_adm_GetFileList(imPar.BidsRoot,[],false,[],true);
    for iSubject = 1:length(listSubjects)

        %% 2. Get subject labels
        subjectLabel = listSubjects{iSubject};

        % Check if the anatomical directory exists
        if exist(fullfile(imPar.BidsRoot,subjectLabel,'anat'),'dir')
            %% 3. Process all anatomical files
            fAnat = xASL_adm_GetFileList(fullfile(imPar.BidsRoot,subjectLabel,'anat'),'^.+\.nii',false,[]);
            for iAnat = 1:length(fAnat)
                %Unzip the file for SPM
                pathUnzipped = xASL_adm_UnzipNifti(fullfile(imPar.BidsRoot,subjectLabel,'anat',fAnat{iAnat}));
                % Remove the face
                xASL_spm_deface(pathUnzipped,true);
                % Zip again
                gzip(pathUnzipped);
                delete(pathUnzipped);
            end
        end
    end

end



