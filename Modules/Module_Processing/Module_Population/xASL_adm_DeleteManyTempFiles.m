function xASL_adm_DeleteManyTempFiles(x)
%xASL_adm_DeleteManyTempFiles This function removes as many files as possible
%
% FORMAT:       xASL_adm_DeleteManyTempFiles(x)
% 
% INPUT:        x          - x structure
%
% OUTPUT:       n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function removes as many files as possible.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_adm_DeleteManyTempFiles(x);
% __________________________________
% Copyright 2015-2024 ExploreASL

Files2DelNativeStructural = {}; % (currently no files here)
Files2DelNativeASL = {'ATT_BiasField.nii' 'Mask_Template.nii' 'Mean_CBF_Template.nii' 'PseudoCBF.nii' 'RawTemplate.nii' 'VascularArtifact_Template.nii' 'mean_PWI_Clipped.nii' 'SliceGradient_extrapolated.nii' 'FoV.nii'};
Files2DelStandardStructural = {'rT1_ORI_' 'rT1_' 'rFLAIR_' '(m|)rc\dT1_'};
Files2DelStandardASL = {'noSmooth_M0_' 'mean_control_' 'PWI_' 'SliceGradient_' 'SNR'};

%% First, we create a list of subjects and sessions of which the structural and ASL modules, respectively, are completed
[listSubjects, listSessions] = xASL_adm_GetListCompletedModules(x);

if isfield(x,'D') && isfield(x.dir,'xASLDerivatives')

    %% 1. Native space — structural module files (currently no files here)
    % fprintf('Deleting temporary structural module files in native space folders:    ');
    % 
    % for iP=1:length(Files2DelNativeStructural)
    %     xASL_TrackProgress(iP, length(Files2DelNativeStructural));
    % 
    %     for iSubject=1:length(listSubjects)
    %         % Delete this file, if it exists
    %         xASL_delete(fullfile(x.dir.xASLDerivatives, listSubjects{iSubject}, Files2DelNativeStructural{iP}));
    %     end
    % end
    % fprintf('\n');
    
    %% 2. Native space — ASL module files
    fprintf('Deleting temporary ASL files in native space folders:    ');

    for iP=1:length(Files2DelNativeASL)
        xASL_TrackProgress(iP, length(Files2DelNativeASL));
        
        for iSession=1:length(listSessions)
            sessionIndex = regexp(listSessions{iSession}, 'ASL_\d*');
            sessionIs = listSessions{iSession}(sessionIndex:end);
            subjectIs = listSessions{iSession}(1:sessionIndex-2);
            
            % Delete this file, if it exists
            xASL_delete(fullfile(x.dir.xASLDerivatives, subjectIs, sessionIs, Files2DelNativeASL{iP}));
        end
    end
    fprintf('\n');
else
    fprintf('Missing root directory, deleting temporary files failed...\n');
end


%% 3. Standard space — structural module files
fprintf('Deleting temporary structural files in population folder:    ');

for iPath=1:length(Files2DelStandardStructural)
    xASL_TrackProgress(iP,length(Files2DelStandardStructural));

    for iSubject=1:length(listSubjects)
        % Delete any existing files
        xASL_adm_DeleteFileList(x.D.PopDir, ['^' Files2DelStandardStructural{iPath} listSubjects{iSubject} '.*$'], 1, [0 Inf]);
    end
end
fprintf('\n');

%% 4. Standard space — ASL module files
%% Now we do the same for the standard space files in the population folder
fprintf('Deleting temporary ASL files in population folder:    ');

for iPath=1:length(Files2DelStandardASL)
    xASL_TrackProgress(iP,length(Files2DelStandardASL));
    
    for iSession=1:length(listSessions)
        % Delete any existing files
        xASL_adm_DeleteFileList(x.D.PopDir, ['^' Files2DelStandardASL{iPath} listSessions{iSession} '.*$'], 1, [0 Inf]);
    end
end
fprintf('\n');


end