function xASL_adm_CleanUpBeforeRerun(AnalysisDir, iModule, bRemoveWMH, bFullRerun, SubjectID)
%xASL_adm_CleanUpBeforeRerun Delete previous results for fresh rerun
%
% FORMAT: xASL_adm_CleanupBeforeCompleteRerun(AnalysisDir, iModule, bRemoveWMH, bFullRerun, SubjectID)
%
% INPUT:
%   AnalysisDir - path to analysis folder containing derivatives in
%                 ExploreASL data structure (REQUIRED)
%   iModule     - integer value for the module from which to clean the
%                 derivatives, 1 = structural, 2 = ASL and 3 = population module, e.g. [1
%                 3] will clean the structural and population module, but not the ASL
%                 module (REQUIRED when bFullRerun==false)
%   bRemoveWMH  - remove WMH segmentations (WMH_SEGM.nii) 1 = removal, 0 = keep
%                 them. (OPTIONAL, DEFAULT = false)
%   bFullRerun  - true for removing all derivatives from all subjects
%                 (OPTIONAL, DEFAULT = false)
%   SubjectID   - string containing subject name/ID from which the
%                 derivatives should be removed (REQUIRED if
%                 bFullRerun==false)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function (partly) reverts previous ExploreASL runs,
%              deleting derivatives, while keeping raw data intact.
%              if bFullRerun==true, then all subjects and all module
%              derivatives will be removed. This function performs the
%              following steps:
%   

%
% NB: still need to add xASL_module_func & xASL_module_dwi for EPAD
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE for rerunning full study: xASL_adm_CleanupBeforeCompleteRerun('/PathToMyStudyAnalysisDir', [], 0, 1);
% EXAMPLE for rerunning ASL & population modules of single subject: xASL_adm_CleanupBeforeCompleteRerun('/PathToMyStudyAnalysisDir', [2 3], 0, 0, 'Sub-001')
% __________________________________
% Copyright 2015-2019 ExploreASL



% ===========================================================================================
%% 0) Admin
if nargin<4 || isempty(bFullRerun)
    bFullRerun = false;
elseif bFullRerun
    iModule = [1 2 3];
elseif ~bFullRerun
    % continue
else
    error(['Invalid bFullRerun option' xASL_num2str(bFullRerun)]);
end
if nargin<3 || isempty(bRemoveWMH)
    bRemoveWMH = false;
end
if nargin<2 || isempty(iModule) || isempty(AnalysisDir)
    warning('Incomplete input arguments');
    return;
end
    
if bFullRerun
    warning('Cleaning data from all subjects!');
    SubjectID = '.*';
    SubjectDir = AnalysisDir;
elseif nargin<5 || isempty(SubjectID)
    fprintf('No subjectID specified, skipping...\n');
    return;
else
    warning(['Cleaning data to rerun subject ' SubjectID]);
    SubjectDir = fullfile(AnalysisDir, SubjectID);
    % Obtain ASL session dirs as well
    SessionDirs = xASL_adm_GetFileList(SubjectDir, '^ASL_\d+$', 'FPList', [0 Inf], true);
    nSessions = length(SessionDirs);
end
if any(iModule>3)
    error('Invalid iModule input argument');
end

PopulationDir = fullfile(AnalysisDir, 'Population');

% To save computation/harddisk burden/time, combine as many regexps here as
% possible. But need to take care not to delete e.g. _ORI or _Backup
NativeSpaceFiles{1} = {'(xASL_module|)Structural(_module|)\.log$' '.*\.ps$' 'xASL_Report.*' '(r|)c(1|2|3)T1(\.mat|\.nii|\.nii\.gz)$'...
    '(ples.*|WMH_SEGM_CleanUp|j_T1|rT1|y_T1|CentralWM_QC|LeftRight)(\.mat|\.nii|\.nii\.gz)$' 'catreport_T1\.pdf$'...
    'T1_seg8\.mat$'};

NativeSpaceFiles{2} = {'(xASL_module_|)ASL(_module|)\.log' '(VascularArtifact|xASL_qc).*' 'rp_(despiked_|)ASL4D\.txt'...
    '(PVgm|PVwm|PVcsf|FoV|M0_biasfield|(m|)mean_control|(mean|)PWI|SD|SNR|TopUp|RawTemplate|ATT_bias|B0|qCBF)(\.mat|\.nii|\.nii\.gz)$'...
    '(SliceGradient(_extrapolated|)|slice_gradient|rgrey|rASL4D|y_ASL|Pseudo(CBF|Tissue)|rM0|Field|Unwarped|rp_ASL4D)(\.mat|\.nii|\.nii\.gz)$'...
    '(CBF(|_Visual2DICOM)|rtemp|temp|Mask(|Vascular)|mask|r|despiked_ASL4D)(\.mat|\.nii|\.nii\.gz)$'};

if bRemoveWMH
    NativeSpaceFiles{1}{end+1} = 'WMH_SEGM(\.mat|\.nii|\.nii\.gz)$';
end


%% 1) If a Population folder doesn't exist yet but dartel does, rename it
DARTELdir = fullfile(AnalysisDir,'dartel');
if exist(DARTELdir,'dir') && ~exist(PopulationDir, 'dir')
    xASL_Move(DARTELdir, PopulationDir);
elseif exist(DARTELdir,'dir') && exist(PopulationDir, 'dir')
    error('Please first merge or delete old DARTEL dir with PopulationDir, before continuing');
elseif ~exist(PopulationDir,'dir')
    xASL_adm_CreateDir(PopulationDir);
end


% ===========================================================================================
%% 2) Remove whole-study data files in AnalysisDir if bFullRerun
if bFullRerun
    RootFiles = {'AcquisitionTime' 'CheckOrientation_RigidRegASL' 'FileReportSummary' '.*Motion' 'GM_ICV' 'GMWM' '.*_vol' 'import_log' 'QC_collection\.json' '.*_count' 'xASL\.mat'};
    for iFile=1:length(RootFiles)
        xASL_adm_DeleteFileList(AnalysisDir, [RootFiles{iFile} '.*'], false, [0 Inf]);
    end
end


% ===========================================================================================
%% 3) Remove lock files/folders for reprocessing
LockDir = fullfile(AnalysisDir, 'lock');

if bFullRerun % remove all lock dirs/files
    xASL_adm_DeleteFileList(LockDir, '.*', true, [0 Inf]);
    if isempty(xASL_adm_GetFileList(LockDir, '.*', 'FPListRec', [0 Inf]))
        % then remove subfolders
        DirList = xASL_adm_GetFileList(LockDir, '.*', 'FPList', [0 Inf], true);
        for iList=1:length(DirList)
            xASL_delete(DirList{iList});
        end
    else
        warning(['Couldnt delete contents of ' LockDir']);
    end
else
    LockDirs = {'xASL_module_Structural' 'xASL_module_ASL' 'xASL_module_Population'};
    LockDirs2 = fullfile(LockDir, LockDirs);

    for iDir=iModule
        if ~isempty(regexp(LockDirs2{iDir},'xASL_module_Population'))
            CurrentDir = {LockDirs2{iDir}}; % specify for population module, no subject ID
        else % specify subject ID
            CurrentDir = xASL_adm_GetFileList(LockDirs2{iDir}, SubjectID, 'FPList', [0 Inf], true);
        end

        if ~isempty(CurrentDir)
            for iCurrent=1:length(CurrentDir)
                xASL_adm_DeleteFileList(CurrentDir{iCurrent}, '.*', true, [0 Inf]);
                LockedDir = xASL_adm_GetFileList(CurrentDir{iCurrent}, 'locked', 'FPListRec', [0 Inf], true);
                if isempty(LockedDir) % keep locked dir for mutex
                    xASL_delete(CurrentDir{iCurrent});
                end
            end
        end
    end
end

% ===========================================================================================
%% 4) Restore backupped _ORI (original) files
% Here we search recursively for ORI files, if any found, and the original as well, replace original by ORI
OriList = ''; % initiate list
if ~isempty(find(iModule==1)) % if we remove the structural data
    OriList{end+1} = xASL_adm_GetFileList(SubjectDir, '.*_ORI\.nii$', 'FPList', [0 Inf]); % within the native space SubjectDir
    OriList{end+1} = xASL_adm_GetFileList(PopulationDir, ['.*_ORI_' SubjectID '(?!_ASL_\d+)\.nii$'], 'FPListRec', [0 Inf]); % within the standard space PopulationDir
elseif ~isempty(find(iModule==2)) % if we remove the ASL data
    for iSession=1:nSessions
        OriList{end+1} = xASL_adm_GetFileList(SessionDirs{iSession}, '.*_ORI\.nii$', 'FPList', [0 Inf]); % within the native space SessionDir
    end
    OriList{end+1} = xASL_adm_GetFileList(PopulationDir, ['.*_ORI_' SubjectID '_ASL_\d+\.nii$'], 'FPListRec', [0 Inf]); % within the standard space PopulationDir
end
% Merge lists
OriListTotal = '';
for iList=1:length(OriList)
    OriListTotal = [OriListTotal OriList{iList}]; % combine the lists
end
    
for iList=1:length(OriListTotal)
    NonOriPath = strrep(OriListTotal{iList},'_ORI','');
    xASL_delete(NonOriPath); % if the non-ori file existed, overwrite it
    xASL_Move(OriListTotal{iList}, NonOriPath);
end
    

% ===========================================================================================
%% 5) Delete native space CAT12 temporary folders (always, independent of iModule)
CAT12Folders = {'label' 'mri' 'report'};
for iFolder=1:length(CAT12Folders)
    FolderList = xASL_adm_GetFileList(SubjectDir, CAT12Folders{iFolder}, 'FPListRec', [0 Inf], true);
    for iDir=1:length(FolderList)
        xASL_adm_DeleteFileList(FolderList{iDir}, '.*', true, [0 Inf]);
        xASL_delete(FolderList{iDir});
    end
end


% ===========================================================================================
%% 6) Remove native space files for iModule
for iList=iModule
    if iList<3 % structural & ASL modules only, population module doesnt have native space files
        for iFile=1:length(NativeSpaceFiles{iList})
            % NB: here we delete files recursively!
            xASL_adm_DeleteFileList(SubjectDir, NativeSpaceFiles{iList}{iFile}, true, [0 Inf]);
        end
    end
end

% ===========================================================================================
%% 7) Remove standard space files for iModule
if bFullRerun
    % Remove the full population folder contents
    xASL_adm_DeleteFileList(PopulationDir, '.*', true, [0 Inf]);
    % Remove all empty folders within the population folder, for overview
    FolderList = xASL_adm_GetFileList(PopulationDir, '.*', 'FPListRec', [0 Inf], true);
    for iFolder=1:length(FolderList)
        if isempty(xASL_adm_GetFileList(FolderList{iFolder}, '.*', 'FPListRec', [0 Inf]))
            xASL_delete(FolderList{iFolder});
        else
            warning(['Couldnt delete contents of ' FolderList{iFolder}]);
        end
    end
end
if ~isempty(find(iModule==1)) % if we remove the structural data
    xASL_adm_DeleteFileList(PopulationDir, ['.*' SubjectID '(?!_ASL_\d+)'], true, [0 Inf]);
end
if ~isempty(find(iModule==2)) % if we remove the ASL data
    xASL_adm_DeleteFileList(PopulationDir, ['.*' SubjectID '_ASL_\d+'], true, [0 Inf]);
end

% ===========================================================================================
%% 8) Remove population module files
if ~isempty(find(iModule==3))
    PopulationFiles = {'(xASL_module|)Population(_module|)\.log' '.*quantification_parameters\.csv' 'Overview.*\.jpg'};

    for iFile=1:length(PopulationFiles)
        xASL_adm_DeleteFileList(PopulationDir, PopulationFiles{iFile}, true, [0 Inf]);
    end
end


% ===========================================================================================
%% 9) Remove or clean up stored x-struct & QC file
PathX = fullfile(SubjectDir, 'x.mat');
PathQC = fullfile(SubjectDir, ['QC_collection_' SubjectID '.json']);
if ~isempty(find(iModule==1)) && ~isempty(find(iModule==2))
    % if we remove both modules Structural & ASL, remove x.mat & QC_collection.*.json completely
    xASL_delete(PathX);
    xASL_delete(PathQC);
else
    RemoveFields = {'Structural' 'ASL'};
    if ~isempty(find(iModule==1)) % if we remove the structural provenance
        RemoveFields{end+1} = 'Structural';
    end
    if ~isempty(find(iModule==2)) % if we remove the ASL provenance
        RemoveFields{end+1} = 'ASL';
    end    
    
    % keep the x.mat file but remove parts of it
    if exist(PathX, 'file')
        % Load the x structure
        xStruct = load(PathX,'-mat','x');
        % delete any mutex folder that was accidentally created
        if isfield(xStruct.x,'mutex')
            NewMutexFolder = xStruct.x.mutex.Root;
            while exist(NewMutexFolder,'dir')
                xASL_delete(NewMutexFolder);
                NewMutexFolder = fileparts(NewMutexFolder);
            end
            xStruct.x = rmfield(xStruct.x,'mutex'); % remove the mutex field
        end
        x = xStruct; 

        for iMod=iModule
            if iMod~=3 % if not a Population module
                if isfield(x,'Output')
                    if isfield(x.Output, RemoveFields{iMod})
                        x.Output = rmfield(x.Output, RemoveFields{iMod});
                    end
                end
                if isfield(x,'Output_im')
                    if isfield(x.Output_im, RemoveFields{iMod})
                        x.Output_im = rmfield(x.Output_im, RemoveFields{iMod});
                    end
                end
            end
        end

        % Always remove previous paths, irrespective of iModule
        PathFields = {'MyPath' 'SpaghettiDir', 'HistogramDir', 'P', 'D', 'StatsMaps',...
            'SPMDIR', 'SPMpath', 'STRUCT_TEMPLATE_IM', 'LOCKDIR', 'SUBJECTDIR',...
            'SESSIONDIR', 'PathPop_MaskSusceptibility' 'StudyAtlasDir'};

        for iField=1:length(PathFields)
            if isfield(x, PathFields{iField})
                x = rmfield(x, PathFields{iField});
            end
        end
        % Store the new x.mat
        xASL_delete(PathX);
        save(PathX, 'x');
        clear x
    end
    % Do the same for the QC_Collection json
    if exist(PathQC, 'file')    
        QCmat = xASL_import_json(PathQC);
        for iMod=iModule % SAME CODE AS ABOVE
            if iMod~=3 % if not a Population module
                if isfield(QCmat, RemoveFields{iMod})
                    QCmat = rmfield(QCmat, RemoveFields{iMod});
                end
            end
        end    
        xASL_delete(PathQC);
        xASL_adm_SaveJSON(QCmat, PathQC);
    end
end


end
