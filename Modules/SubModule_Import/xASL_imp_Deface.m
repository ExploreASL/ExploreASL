function xASL_imp_Deface(x,imPar)
%xASL_imp_Deface Run defacing.
%
% FORMAT: xASL_imp_Deface(x,imPar)
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
% EXAMPLE:     xASL_imp_Deface(x,imPar);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Initialize
    if nargin<1
        error('Please provide an x and an imPar struct...');
    end
    if nargin<2
        error('Please provide an imPar struct...');
    end
    
    % Make sure that logging is still active
    diary(x.dir.diaryFile);

    % Print feedback
    fprintf('\n[\b========================================== DEFACING ==========================================]\b\n');
    
    
    % We do not iterate over subjects anymore, since this is done in xASL_Iteration now
    iSubject = strcmp(x.SUBJECT,x.SUBJECTS);
    subjectName = x.SUBJECTS{iSubject};
    
    
    %% 1. Iterate over list of subjects
    listSubjects = xASL_adm_GetFileList(imPar.BidsRoot,[],false,[],true);
    for iSubject = 1:length(listSubjects)
        
        % Only run it for the current subject (maybe we can do this more elegantly in the future)
        if ~isempty(regexpi(listSubjects{iSubject},subjectName))

            %% 2. Get subject labels
            subjectLabel = listSubjects{iSubject};

            %% 3. Process all anatomical files
            if exist(fullfile(imPar.BidsRoot,subjectLabel,'anat'),'dir') 
                % Single-session
                fAnat = xASL_adm_GetFileList(fullfile(imPar.BidsRoot,subjectLabel,'anat'),'^.+\.nii',false,[]);
                xASL_imp_RunDeface(imPar,fAnat,subjectLabel,[]);
            else
                % Multi-session
                sessionDirs = xASL_adm_GetFileList(fullfile(imPar.BidsRoot,subjectLabel),[],false,[],true);
                for iSession = 1:numel(sessionDirs)
                    if exist(fullfile(imPar.BidsRoot,subjectLabel,sessionDirs{iSession},'anat'),'dir')
                        sessionName = sessionDirs{iSession};
                        fAnat = xASL_adm_GetFileList(fullfile(imPar.BidsRoot,subjectLabel,sessionName,'anat'),'^.+\.nii',false,[]);
                        xASL_imp_RunDeface(imPar,fAnat,subjectLabel,sessionName);
                    end
                end
            end
            
        end
    end

end



% Actual defacing
function xASL_imp_RunDeface(imPar,fAnat,subjectLabel,sessionName)

    % Check that list is not empty
    if ~isempty(fAnat)
        for iAnat = 1:length(fAnat)
            % Get filename
            if nargin<4 || isempty(sessionName)
                fileName = fullfile(imPar.BidsRoot,subjectLabel,'anat',fAnat{iAnat});
            else
                fileName = fullfile(imPar.BidsRoot,subjectLabel,sessionName,'anat',fAnat{iAnat});
            end
            
            % Print feedback
            fprintf('\nDeface %s...\n',fAnat{iAnat});
            
            % Check if file exists
            if ~xASL_exist(fileName)
                % Print warning
                warning('Defacing was not run, file %s not found...', fAnat{iAnat});
            else
                % Unzip the file for SPM
                pathUnzipped = xASL_adm_UnzipNifti(fileName);
                % Remove the face
                xASL_spm_deface(pathUnzipped,true);
                % Zip again
                gzip(pathUnzipped);
                delete(pathUnzipped);
            end
            
        end
    end

end


