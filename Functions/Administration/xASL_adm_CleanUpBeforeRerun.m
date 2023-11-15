function xASL_adm_CleanUpBeforeRerun(AnalysisDir, iModule, bRemoveWMH, bAllSubjects, SubjectID, SessionID)
%xASL_adm_CleanUpBeforeRerun Delete previous results for fresh rerun
%
% FORMAT: xASL_adm_CleanupBeforeCompleteRerun(AnalysisDir, iModule, bRemoveWMH, bAllSubjects, SubjectID)
%
% INPUT:
%   AnalysisDir - path to analysis folder containing derivatives in
%                 ExploreASL data structure (REQUIRED)
%   iModule     - integer value for the module from which to clean the
%                 derivatives, 1 = structural, 2 = ASL and 3 = population module, e.g. [1
%                 3] will clean the structural and population module, but not the ASL
%                 module (REQUIRED when bAllSubjects==false)
%   bRemoveWMH  - remove WMH segmentations (WMH_SEGM.nii) 1 = removal, 0 = keep
%                 them. (OPTIONAL, DEFAULT = false)
%   bAllSubjects  - true for removing all derivatives from all subjects
%                 (OPTIONAL, DEFAULT = false)
%   SubjectID   - string containing subject name/ID from which the
%                 derivatives should be removed (REQUIRED if
%                 bAllSubjects==false)
%   SessionID   - string containing session ID (e.g. 'ASL_1', 'ASL_2', etc... ASL_n') 
%                 from which the derivatives should be removed (OPTIONAL,
%                 DEFAULT = detect ASL sessions automatically & remove all
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function (partly) reverts previous ExploreASL runs,
%              deleting derivatives, while keeping raw data intact.
%              if bAllSubjects==true, then all subjects and all module
%              derivatives will be removed. This function performs the
%              following steps:
%   
%              1. If a Population folder doesn't exist yet but dartel does, rename it
%              2. Remove whole-study data files in AnalysisDir if bAllSubjects
%              3. Remove lock files/folders for reprocessing
%              4. Restore backupped _ORI (original) files
%              5. Delete native space CAT12 temporary folders (always, independent of iModule)
%              6. Remove native space files for iModule
%              7. Remove standard space files for iModule
%              8. Remove population module files
%              9. Remove or clean up stored x-struct & QC file -> THIS HAS NO SESSION SUPPORT YET
%
% NB: still need to add xASL_module_func & xASL_module_dwi for EPAD
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
% - For rerunning full study: xASL_adm_CleanupBeforeCompleteRerun('/PathToMyStudyAnalysisDir', [], 0, 1);
% - For rerunning ASL & population modules of single subject: xASL_adm_CleanupBeforeCompleteRerun('/PathToMyStudyAnalysisDir', [2 3], 0, 0, 'Sub-001')
%
% __________________________________
% Copyright 2015-2020 ExploreASL


try

    % ===========================================================================================
    %% 0) Admin
    if nargin<4 || isempty(bAllSubjects) || ~bAllSubjects
        bAllSubjects = false;
    elseif bAllSubjects && (nargin<2 || isempty(iModule))
        iModule = [1 2 3];
    elseif bAllSubjects
        % continue
    else
        error(['Invalid bAllSubjects option' xASL_num2str(bAllSubjects)]);
    end
    if nargin<3 || isempty(bRemoveWMH)
        bRemoveWMH = false;
    end
    if isempty(AnalysisDir)
        warning('Incomplete input arguments');
        return;
    elseif ~exist(AnalysisDir, 'dir')
        warning('Invalid folder input argument');
        return;
    end

    if bAllSubjects
        warning('Cleaning data from all subjects!');
        SubjectID = '.*';
        SessionID = {'ASL_\d+'};
        SubjectDir = AnalysisDir;
        SessionDir = AnalysisDir;
        nSessions = 1;
    elseif nargin<5 || isempty(SubjectID)
        fprintf('No subjectID specified, skipping...\n');
        return;
    else
        fprintf(['Cleaning data to rerun subject ' SubjectID '\n']);
        SubjectDir = fullfile(AnalysisDir, SubjectID);
    end
    if (nargin<6 || isempty(SessionID)) && ~bAllSubjects
        % Obtain ASL session dirs as well
        fprintf('No SessionID provided, finding session dirs automatically:\n');
        SessionID = xASL_adm_GetFileList(SubjectDir, '^ASL_\d+$', 'List', [0 Inf], true);
    else
        if ~iscell(SessionID)
            SessionID = {SessionID};
        end
    end
    if ~bAllSubjects && ~isempty(SessionID)
        SessionDir = cellfun(@(y) fullfile(SubjectDir, y), SessionID, 'UniformOutput', 0);
        nSessions = length(SessionDir);
    end

    if any(iModule>3)
        error('Invalid iModule input argument');
    end

    PopulationDir = fullfile(AnalysisDir, 'Population');

    % ===========================================================================================
    %% 1) If a Population folder doesn't exist yet but dartel does, rename it
    DARTELdir = fullfile(AnalysisDir,'dartel');
    if exist(DARTELdir,'dir') && ~exist(PopulationDir, 'dir')
        fprintf('Renaming folders: dartel to Population\n');
        xASL_Move(DARTELdir, PopulationDir);
    elseif exist(DARTELdir,'dir') && exist(PopulationDir, 'dir')
        error('Please first merge or delete old DARTEL dir with PopulationDir, before continuing');
	else
        xASL_adm_CreateDir(PopulationDir);
    end


    % ===========================================================================================
    %% 2) Remove whole-study data files in AnalysisDir if bAllSubjects
    if bAllSubjects
            fprintf('Removing root data files:   ');
        RootFiles = {'AcquisitionTime' 'CheckOrientation_RigidRegASL' 'FileReportSummary' '.*Motion' 'GM_ICV' 'GMWM' '.*_vol' 'import_log' 'QC_collection\.json' '.*_count' 'xASL\.mat'};
        for iFile=1:length(RootFiles)
            xASL_TrackProgress(iFile, length(RootFiles));
            xASL_adm_DeleteFileList(AnalysisDir, [RootFiles{iFile} '.*'], false, [0 Inf]);
        end
    end
    fprintf('n');


    % ===========================================================================================
    %% 3) Remove lock files/folders for reprocessing
    LockDir = fullfile(AnalysisDir, 'lock');
    LockDirs = {'xASL_module_Structural' 'xASL_module_ASL' 'xASL_module_Population'};
    LockDirs2 = fullfile(LockDir, LockDirs);
    
    fprintf('Removing lock folders:   ');

    for iDir=iModule
        CurrentDir = {};
        % remove population lock files
        if ~isempty(regexp(LockDirs2{iDir},'xASL_module_Population'))
            xASL_adm_DeleteFileList(LockDirs2{iDir}, '.*', true, [0 Inf]);
        
        elseif bAllSubjects % and not population folder
            % remove for all subjects
            xASL_adm_DeleteFileList(LockDirs2{iDir}, '.*', true, [0 Inf]);
            if isempty(xASL_adm_GetFileList(LockDirs2{iDir}, '.*', 'FPListRec', [0 Inf]))
                % then remove subfolders
                DirList = xASL_adm_GetFileList(LockDirs2{iDir}, '.*', 'FPList', [0 Inf], true);
                for iList=1:length(DirList)
                    xASL_TrackProgress(iList, length(DirList));
                    xASL_delete(DirList{iList});
                end
            else
                warning(['Couldnt delete contents of ' LockDir']);
            end
            
        else % remove for specific subject/session
            if ~isempty(regexp(LockDirs2{iDir},'xASL_module_Structural'))
                % specify subject ID
                CurrentDir(end+1) = xASL_adm_GetFileList(LockDirs2{iDir}, SubjectID, 'FPList', [0 Inf], true);         
            elseif ~isempty(regexp(LockDirs2{iDir},'xASL_module_ASL'))
                % specify session ID
                if exist(LockDirs2{iDir}, 'dir')
                    for iSession=1:nSessions
                        CurrentDir(end+1) = xASL_adm_GetFileList(fullfile(LockDirs2{iDir}, SubjectID), ['^xASL_module_ASL_' SessionID{iSession} '$'], 'FPList', [0 Inf], true);
                    end
                end
            end
        end
        
        if ~isempty(CurrentDir)
            for iCurrent=1:length(CurrentDir)
                xASL_TrackProgress(iCurrent, length(CurrentDir));
                if ~isempty(CurrentDir{iCurrent}) % bugfix, otherwise this deletes any files in the current folder (e.g. FLAIR/T1)
                    xASL_adm_DeleteFileList(CurrentDir{iCurrent}, '.*', true, [0 Inf]);
                    LockedDir = xASL_adm_GetFileList(CurrentDir{iCurrent}, 'locked', 'FPListRec', [0 Inf], true);
                    if isempty(LockedDir) % keep locked dir for mutex
                        xASL_delete(CurrentDir{iCurrent});
                    end
                end
            end
        end
    end
    fprintf('\n');


    % ===========================================================================================
    %% 4) Restore backupped _ORI (original) files
    % Here we search recursively for ORI files, if any found, and the original as well, replace original by ORI
    OriList = {}; % initiate list
    fprintf('Restoring backupped _ORI (original) files:   ');

    if bAllSubjects
        OriList = [OriList;xASL_adm_GetFileList(SubjectDir, '(?i).*_ORI\.nii$', 'FPListRec', [0 Inf])]; % for all subjects/sessions
    elseif ~isempty(find(iModule==1)) % if we remove the structural data
        OriList = [OriList;xASL_adm_GetFileList(SubjectDir, '(?i).*_ORI\.nii$', 'FPList', [0 Inf])]; % within the native space SubjectDir
        OriList = [OriList;xASL_adm_GetFileList(PopulationDir, ['(?i).*_ORI_' SubjectID '(?!_ASL_\d+)\.nii$'], 'FPListRec', [0 Inf])]; % within the standard space PopulationDir
    elseif ~isempty(find(iModule==2)) % if we remove the ASL data
        for iSession=1:nSessions
            OriList = [OriList;xASL_adm_GetFileList(SessionDir{iSession}, '(?i).*_ORI\.nii$', 'FPList', [0 Inf])]; % within the native space SessionDir
            OriList = [OriList;xASL_adm_GetFileList(PopulationDir, ['(?i).*_ORI_' SubjectID '_' SessionID{iSession} '\.nii$'], 'FPListRec', [0 Inf])]; % within the standard space PopulationDir
        end
    end

    for iList=1:length(OriList)
        xASL_TrackProgress(iList, length(OriList));
        NonOriPath = strrep(OriList{iList},'_ORI','');
        xASL_Move(OriList{iList}, NonOriPath, 1, 0); % if the non-ori file existed, overwrite it
    end
    fprintf('\n');


    % ===========================================================================================
    %% 5) Delete native space CAT12 temporary folders (always, independent of iModule)
    fprintf('Deleting native space CAT12 temporary folders:   ');

    CAT12Folders = {'label' 'mri' 'report'};
    for iFolder=1:length(CAT12Folders)
        FolderList = xASL_adm_GetFileList(SubjectDir, CAT12Folders{iFolder}, 'FPListRec', [0 Inf], true);
        for iDir=1:length(FolderList)
            xASL_TrackProgress(iDir, length(FolderList));
            xASL_adm_DeleteFileList(FolderList{iDir}, '.*', true, [0 Inf]);
            xASL_delete(FolderList{iDir});
        end
    end
    fprintf('\n');

    % ===========================================================================================
    %% 6) Remove native space files for iModule
    % To save computation/harddisk burden/time, combine as many regexps here as
    % possible. But need to take care not to delete e.g. _ORI or _Backup
    NativeSpaceFiles{1} = {'.*\.ps$' 'xASL_Report.*' '(r|)c(1|2|3)T1(\.mat|\.nii|\.nii\.gz)$' 'WADQC.*\.json' ...
        '(ples.*|WMH_SEGM_CleanUp|j_T1|rT1|y_T1|CentralWM_QC|LeftRight)(\.mat|\.nii|\.nii\.gz)$' 'catreport_T1\.pdf$'...
        'T1_seg8\.mat$'};

    NativeSpaceFiles{2} = {'.*\.topup.*log' 'WADQC.*\.json' '.*(rep_|BeforeSpikeExclusion|despiked).*' '(VascularArtifact|xASL_qc|qcdc).*' 'rp_ASL4D\.txt' '(ASL4D|M0|.*_sn)\.mat'...
        '(MaskASL|Mean_CBF_Template|mean_PWI_Clipped|PVgm|PVwm|CBFgm|CBFwm|PVcsf|PVwmh|FoV|M0_biasfield|(m|)mean_control|(mean|)PWI|SD|SNR|TopUp.*|(Mask_|Raw)Template|ATT_bias|B0|(Mean_|)(q|)CBF)(\.mat|\.nii|\.nii\.gz)$'...
        '(SliceGradient(_extrapolated|)|slice_gradient|rgrey|PWI4D|y_ASL|Pseudo(CBF|Tissue)|rM0|Field|Unwarped|rASL4D)(\.mat|\.nii|\.nii\.gz)$'...
        '(.*beforeMoCo|CBF(|_Visual2DICOM)|(r|)temp.*|Mask(|Vascular)|mask)(\.mat|\.nii|\.nii\.gz)$'};

    if bRemoveWMH
        NativeSpaceFiles{1}{end+1} = 'WMH_SEGM(\.mat|\.nii|\.nii\.gz)$';
    end

    if bAllSubjects % only remove log files when running this for the full population
        NativeSpaceFiles{1}{end+1} = '.*\.log$';
        NativeSpaceFiles{2}{end+1} = '.*\.log$';
    end

    for iList=iModule
        if iList~=3 % skip population module here
            fprintf(['Deleting native space files for module ' num2str(iList) ':   ']);
            % We do Structural & ASL modules only, population module doesnt have native space files
            if iList==1 % structural module
                Dir2Check = SubjectDir;
            elseif iList==2 && ~isempty(SessionID)
                Dir2Check = SessionDir;
            else
                Dir2Check = '';
            end
            if ~iscell(Dir2Check)
                Dir2Check = {Dir2Check};
            end

            for iFile=1:length(NativeSpaceFiles{iList})
                xASL_TrackProgress(iFile, length(NativeSpaceFiles{iList}));
                for iDir=1:length(Dir2Check)
                    if bAllSubjects
                        % do this recursively throughout all subject/session folders
                        xASL_adm_DeleteFileList(Dir2Check{iDir}, NativeSpaceFiles{iList}{iFile}, true, [0 Inf]);
                    else
                        % don't do this recursive, so that we can remove
                        % individual sessions separately
                        xASL_adm_DeleteFileList(Dir2Check{iDir}, NativeSpaceFiles{iList}{iFile}, false, [0 Inf]);
                    end
                end
            end
            fprintf('\n');
        end
    end

    % ===========================================================================================
    %% 7) Remove standard space files for iModule
    if ~isempty(find(iModule==1)) % if we remove the structural data
        fprintf('Deleting all standard space files for module 1\n');
        xASL_adm_DeleteFileList(PopulationDir, ['.*' SubjectID '(?!_ASL_\d+)'], true, [0 Inf]);
    end
    if ~isempty(find(iModule==2)) && ~isempty(SessionID) % if we remove the ASL data
        fprintf('Deleting all standard space files for module 2\n');
        for iSession=1:nSessions
            xASL_adm_DeleteFileList(PopulationDir, ['.*' SubjectID '_' SessionID{iSession}], true, [0 Inf]);
        end
    end

    % ===========================================================================================
    %% 8) Remove population module files
    if ~isempty(find(iModule==3))
        PopulationFiles = {'(xASL_module|)Population(_module|)\.log' '.*quantification_parameters\.csv' 'Overview.*\.jpg'};
        fprintf('Deleting all standard space files for module 3:   ');

        for iFile=1:length(PopulationFiles)
            xASL_TrackProgress(iFile, length(PopulationFiles));
            xASL_adm_DeleteFileList(PopulationDir, PopulationFiles{iFile}, true, [0 Inf]);
        end
    end
    fprintf('\n');

    % ===========================================================================================
    %% 9) Remove or clean up stored x-struct & QC file
    % THIS HAS NO SESSION SUPPORT YET
    if ~isempty(find(iModule==1)) && ~isempty(find(iModule==2))
        % if we remove both modules Structural & ASL, remove x.mat & QC_collection.*.json completely
        fprintf('Removing x-struct and QC-file\n');
        xASL_adm_DeleteFileList(SubjectDir, '^x\.mat$', true, [0 Inf]);
        xASL_adm_DeleteFileList(SubjectDir, ['^QC_collection_' SubjectID '\.json$'], true, [0 Inf]);
    else
        ListXmat = xASL_adm_GetFileList(SubjectDir, '^x\.mat$', 'FPListRec', [0 Inf]);
        ListQCjson = xASL_adm_GetFileList(SubjectDir, ['(?i)^QC_collection_' SubjectID '\.json$'], 'FPListRec', [0 Inf]);

        RemoveFields = {};
        if ~isempty(find(iModule==1)) % if we remove the structural provenance
            RemoveFields{end+1} = 'Structural';
        end
        if ~isempty(find(iModule==2)) % if we remove the ASL provenance
            RemoveFields{end+1} = 'ASL';
        end        
        
        % keep the x.mat file but remove parts of it        
        if ~isempty(ListXmat)
            fprintf('Cleaning up x-struct:   ');

            for iX=1:length(ListXmat)
                xASL_TrackProgress(iX, length(ListXmat));
                if exist(ListXmat{iX}, 'file')
                    % Load the x structure
                    xStruct = load(ListXmat{iX},'-mat','x'); % We should use xASL_adm_LoadX here!
                    % delete any mutex folder that was accidentally created
                    if isfield(xStruct.x,'mutex')
                        NewMutexFolder = fullfile(pwd,xStruct.x.mutex.Root);
                        while exist(NewMutexFolder,'dir') && length(NewMutexFolder)>length(pwd) && isempty(xASL_adm_GetFileList(NewMutexFolder,'.*'))
                            xASL_delete(NewMutexFolder);
                            NewMutexFolder = fileparts(NewMutexFolder);
                        end
                        xStruct.x = rmfield(xStruct.x,'mutex'); % remove the mutex field
                    end
                    x = xStruct.x;

                    % Optional: Remove Output and Output_im fields
                    for iMod=iModule
                        if iMod~=3 % if not a Population module
                            if isfield(x, 'Output')
                                if isfield(x.Output, RemoveFields{iMod})
									if strcmpi(RemoveFields{iMod}, 'structural')
										x.Output = rmfield(x.Output, RemoveFields{iMod});
									else
										 for iSession=1:nSessions
											 if isfield(x.Output.(RemoveFields{iMod}), SessionID{iSession})
												 x.Output.(RemoveFields{iMod}) = rmfield(x.Output.(RemoveFields{iMod}), SessionID{iSession});
											 end
										 end
									end
                                end
                            end
                            if isfield(x,'Output_im')
                                if isfield(x.Output_im, RemoveFields{iMod})
									% If the module is structural or if the subfield, e.g., x.Output_im.ASL is a cell, then it contains the outdated list of images and should be delete completely
									if strcmpi(RemoveFields{iMod}, 'structural') || iscell(x.Output_im.(RemoveFields{iMod}))
										x.Output_im = rmfield(x.Output_im, RemoveFields{iMod});
									else
										% Only in case we have non-structural subfields for sessions, then we remove the one for the current session
										for iSession=1:nSessions
											 if isfield(x.Output_im.(RemoveFields{iMod}), SessionID{iSession})
												 x.Output_im.(RemoveFields{iMod}) = rmfield(x.Output_im.(RemoveFields{iMod}), SessionID{iSession});
											 end
										 end
									end
                                end
                            end
                        end
                    end

                    % Always remove previous paths, irrespective of iModule
                    PathFields = {'P', 'D', 'STRUCT_TEMPLATE_IM', 'LOCKDIR', 'PathPop_MaskSusceptibility' 'StudyAtlasDir'};

                    for iField=1:length(PathFields)
                        if isfield(x, PathFields{iField})
                            x = rmfield(x, PathFields{iField});
                        end
                    end
                    
                    % Remove symbolic path of SUBJECTDIR & MyPath
                    if isfield(x.dir,'SUBJECTDIR')
                        x.dir = rmfield(x.dir, 'SUBJECTDIR');
                    end
                    if isfield(x.dir,'SESSIONDIR')
                        x.dir = rmfield(x.dir, 'SESSIONDIR');
                    end
                    if isfield(x.opts,'MyPath')
                        x.opts = rmfield(x.opts, 'MyPath');
                    end
                    
                    % Store the new x.mat
                    xASL_delete(ListXmat{iX});
                    save(ListXmat{iX}, 'x');
                    clear x
                end
            end
            fprintf('\n');
        end
        
        if ~isempty(ListQCjson)
            % Do the same for the QC_Collection json
            fprintf('Cleaning up QC_Collection json:   ');
            for iQC=1:length(ListQCjson)
                xASL_TrackProgress(iQC, length(ListQCjson));
                if exist(ListQCjson{iQC}, 'file')
                    QCmat = xASL_io_ReadJson(ListQCjson{iQC});
                    for iMod=iModule % SAME CODE AS ABOVE
                        if iMod~=3 && iMod~=2 % if not a Population or ASL module
                            if isfield(QCmat, RemoveFields{iMod})
                                QCmat = rmfield(QCmat, RemoveFields{iMod});
                            end
                        end
                    end
                    xASL_delete(ListQCjson{iQC});
                    xASL_io_WriteJson(ListQCjson{iQC}, QCmat);
                end
            end
            fprintf('\n');
        end
    end

    % ===========================================================================================
    %% 10) Done
    if bAllSubjects
        fprintf('Done cleaning up!\n');
        fprintf('Remember to remove the derivative .mat & .csv files in the AnalysisDir\n');
        fprintf('But keep those that were manually created!\n');
    end
    
catch ME
   warning(ME.message);
    
end

end
