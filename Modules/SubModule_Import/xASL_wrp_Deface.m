function xASL_wrp_Deface(x)
%xASL_wrp_Deface Run defacing.
%
% FORMAT: xASL_wrp_Deface(x)
% 
% INPUT:
%   x          - ExploreASL x structure (REQUIRED, STRUCT)
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
% EXAMPLE:     xASL_wrp_Deface(x);
% __________________________________
% Copyright 2015-2022 ExploreASL


    %% Initialize
    if nargin<1
        error('Please provide an x struct...');
    end

    if ~isfield(x.modules.import,'BidsRoot')
        error('It seems as if the initialization of the deface module failed...');
    end

    % Print feedback
    xASL_adm_BreakString('DEFACING');
    
    % We do not iterate over subjects anymore, since this is done in xASL_Iteration now
    iSubject = strcmp(x.SUBJECT,x.SUBJECTS);
    subjectName = x.SUBJECTS{iSubject};
    
    
    %% 1. Iterate over list of subjects
    listSubjects = xASL_adm_GetFileList(x.modules.import.BidsRoot, [], false, [], true);
    for iSubject = 1:length(listSubjects)
        
        % Only run it for the current subject (maybe we can do this more elegantly in the future)
        if ~isempty(regexpi(listSubjects{iSubject}, subjectName))

            %% 2. Get subject labels
            subjectLabel = listSubjects{iSubject};

            %% 3. Process all anatomical files
            if exist(fullfile(x.modules.import.BidsRoot, subjectLabel, 'anat'), 'dir') 
                % Single-session
                fileAnat = xASL_adm_GetFileList(fullfile(x.modules.import.BidsRoot, subjectLabel, 'anat'), '^.+\.nii', false, []);
                xASL_imp_RunDeface(x, fileAnat, subjectLabel, []);
            else
                % Multi-session
                sessionDirs = xASL_adm_GetFileList(fullfile(x.modules.import.BidsRoot, subjectLabel), [], false, [], true);
                for iSession = 1:numel(sessionDirs)
                    if exist(fullfile(x.modules.import.BidsRoot, subjectLabel, sessionDirs{iSession}, 'anat'), 'dir')
                        sessionName = sessionDirs{iSession};
                        fileAnat = xASL_adm_GetFileList(fullfile(x.modules.import.BidsRoot, subjectLabel, sessionName, 'anat'), '^.+\.nii', false, []);
                        xASL_imp_RunDeface(x, fileAnat, subjectLabel, sessionName);
                    end
                end
            end
            
        end
    end

end



% Actual defacing
function xASL_imp_RunDeface(x, fileAnat, subjectLabel, sessionName)

    % Check that list is not empty
    if ~isempty(fileAnat)
        for iAnat = 1:length(fileAnat)
            % Get filename
            if nargin<4 || isempty(sessionName)
                fileName = fullfile(x.modules.import.BidsRoot, subjectLabel, 'anat', fileAnat{iAnat});
            else
                fileName = fullfile(x.modules.import.BidsRoot, subjectLabel, sessionName, 'anat', fileAnat{iAnat});
            end
            
            % Print feedback
            fprintf('\nDeface %s...\n', fileAnat{iAnat});
            
            % Check if file exists
            if ~xASL_exist(fileName)
                % Print warning
                warning('Defacing was not run, file %s not found...', fileAnat{iAnat});
            else
                % Unzip the file for SPM
                pathUnzipped = xASL_adm_UnzipNifti(fileName);
                % Remove the face
                xASL_spm_deface(pathUnzipped, true);
                % Zip again
                gzip(pathUnzipped);
                delete(pathUnzipped);
            end
            
        end
    end

end


