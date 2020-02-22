%% Data curation CICERO
% This script moves the data from Nolan's data format into ExploreASL's
% data format

if isempty(which('ExploreASL_Master'))
    cd ../../..;
    if isempty(which('ExploreASL_Master'))
        error('Please start this script in ExploreASL folder');
    end
end

x = ExploreASL_Master('',0); % initialization

if ispc
    Rdir = 'C:\BackupWork\ASL\CICERO_Nolan';
else
    Rdir = '/Volumes/C/BackupWork/ASL/CICERO_Nolan';
    Rdir = '/Users/henk/ExploreASL/CICERO_Nolan';
end
    
xASL_delete(fullfile(Rdir,'.DS_Store'));
RawDir = fullfile(Rdir, 'raw');
AnalysisDir = fullfile(Rdir, 'analysis');
PathCohort = fullfile(AnalysisDir, 'Cohort.mat');
PathScanDate = fullfile(AnalysisDir, 'ScanDate.mat');
xASL_adm_CreateDir(RawDir);
xASL_adm_CreateDir(AnalysisDir);

CohortDirs = {'HEALTHY' 'OCCL-ASYM' 'OCCL-SYMM' 'STEN-ASYM' 'STEN-SYMM'}; % here we specify the cohort folders
% Here we specify which files to look for to copy to the analysis &
% subfolders, throughout the iterations/loops below
FileIn = {'CBF.nii' 'EPI.nii' 'pcasl.nii' 'vepcasl.nii' 'AP1.nii' 'AP2.nii' 'RL.nii'...
    'vepcasl_1_CBF' 'vepcasl_2_RL' 'vepcasl_3_AP1.nii' 'vepcasl_4_AP2.nii' 'ASL.nii'};
FileOut = {'CBF.nii' 'mean_control.nii' 'ASL4D.nii' 'VE_ASL4D.nii' 'AP1.nii' 'AP2.nii' 'RL.nii'...
    'CBF.nii' 'RL.nii' 'AP1.nii' 'AP2.nii' 'ASL4D.nii'};

clear ScanDate Cohort

NextN = 1;
for iCohort=1:length(CohortDirs)
    fprintf(['\nProcessing ' CohortDirs{iCohort} ':   ']);
    MainDir = fullfile(Rdir,CohortDirs{iCohort});
    SubjectList = xASL_adm_GetFileList(MainDir,'^C3T-S\d{3}-\d*$', 'List',[0 Inf], true);
    
    for iSubject=1:length(SubjectList)
        xASL_TrackProgress(iSubject,length(SubjectList));
        SubjectDir = fullfile(MainDir, SubjectList{iSubject});
        
        %% Get SubjectID, manage multiple visits
        if iSubject>1 && strcmp(SubjectList{iSubject}(1:8), SubjectList{iSubject-1}(1:8))
            CurrentSubject = str2num(SubjectID(end))+1;
        else
            CurrentSubject = 1;
        end
            
        SubjectID = [SubjectList{iSubject}(1:8) '_' num2str(CurrentSubject)];
        
        %% Get ScanDate
        ScanDateN = SubjectList{iSubject}(10:end);
        ScanDate{NextN,1} = SubjectID;
        ScanDate{NextN,2} = ScanDateN;

        %% Get Cohort
        Cohort{NextN,1} = SubjectID;
        Cohort{NextN,2} = CohortDirs{iCohort};
        
        NextN = NextN + 1;
        AnalysisSubjectDir = fullfile(AnalysisDir, SubjectID);
        ASLdir = fullfile(AnalysisSubjectDir,'ASL_1');
        
        %% Move raw data
        RawSubjectDir = fullfile(RawDir, SubjectID);
        xASL_adm_CreateDir(RawSubjectDir);
        FileList = xASL_adm_GetFileList(SubjectDir,['^' SubjectID(1:end-2) '.*\.(PAR|REC)$'], 'FPList');
        for iFile=1:length(FileList)
            [~, Ffile, Fext] = fileparts(FileList{iFile});
            DestFile = fullfile(RawSubjectDir,[Ffile, Fext]);
            xASL_Move(FileList{iFile}, DestFile);
        end
        % Delete unnecessary files
        xASL_delete(fullfile(SubjectDir,'.DS_Store'));
        if ~isempty(regexp(SubjectDir,'C3T-S048-20080909'))
            xASL_delete(fullfile(SubjectDir,'leesmij.txt'));
        elseif ~isempty(regexp(SubjectDir,'C3T-S056-20080910'))
            xASL_delete(fullfile(SubjectDir,'C3T-S056-20080910.lnk'));
        elseif ~isempty(regexp(SubjectDir,'C3T-S085-20090113'))
            xASL_delete(fullfile(SubjectDir,'C3T-S085-20090113_13_1_T2_FFE.hdr'));
            xASL_delete(fullfile(SubjectDir,'C3T-S085-20090113_13_1_T2_FFE.img'));
        end
        if ~isempty(xASL_adm_GetFileList(SubjectDir,'^.*$'))
            warning(['Still files in ' SubjectDir]);
        end
       
        %% AnalysisSubjectDir
        xASL_adm_CreateDir(AnalysisSubjectDir);
        xASL_adm_CreateDir(ASLdir);
        
        %% Move ANA - dwi
        DirANA = fullfile(SubjectDir,'ANA');
        Files2Copy = {'.*(ADC|adc).*\.nii' '.*(DWI|dwi).*\.nii'};
        DestDir = fullfile(AnalysisSubjectDir, 'dwi');
        xASL_adm_CreateDir(DestDir);
        for iFile=1:length(Files2Copy)
            ListIs = xASL_adm_GetFileList(DirANA,Files2Copy{iFile}, 'FPList');
            for iList=1:length(ListIs)
                [~, Ffile, Fext] = xASL_fileparts(ListIs{iList});
                xASL_Move(ListIs{iList}, fullfile(DestDir,[Ffile, Fext]));
            end
        end
        %% Move ANA - FLAIR
        PathOrig = fullfile(DirANA,'FLAIR.nii');
        PathDest = fullfile(AnalysisSubjectDir, 'FLAIR.nii');
        xASL_delete(PathOrig);
%         if xASL_exist(PathOrig,'file') % too low FLAIR resolution to be
%         useful for image processing
%             xASL_Move(PathOrig, PathDest);
%         end
        
        % Delete unnecessary files
        xASL_adm_DeleteFileList(DirANA, '.*\.(png|jpg)', true, [0 Inf]);
        xASL_delete(fullfile(DirANA,'FIG'));
        
        DirsAre = xASL_adm_GetFileList(DirANA,'^ACQ\d*','FPList',[0 Inf], true);
        if length(DirsAre)>2
            fprintf('%s\n',['Too many folders in ' DirANA]);
        end
        for iDir=1:length(DirsAre)
            xASL_delete(fullfile(DirsAre{iDir},'FIG'));
            xASL_adm_DeleteFileList(DirsAre{iDir}, '^(r|)(x|)(DWI|ADC|FLAIR|t1).*\.nii$', false, [0 Inf]);
        end
        if ~isempty(xASL_adm_GetFileList(DirANA,'^.*$','FPListRec'))
            warning(['Still files in ' DirANA]);
        else
            xASL_delete(DirANA);
        end

        %% Move M0 data
        DirM0 = fullfile(SubjectDir,'M0');
        
        FileOrig = {'T1.nii' 'brain.nii' 'T1_gm.nii' 'T1_wm.nii' 'T1_csf.nii' 'T1_inv_sn.mat' 'T1_sn.mat'};
        FileDest = {'T1.nii' 'T1_mask.nii' 'Nolan_c1T1.nii'  'Nolan_c2T1.nii' 'Nolan_c3T1.nii' 'Nolan_T1_inv_sn.mat' 'Nolan_T1_sn.mat'};
        for iPath=1:length(FileOrig)
            PathOrig = fullfile(DirM0,FileOrig{iPath});
            PathDest = fullfile(AnalysisSubjectDir,FileDest{iPath});
            if xASL_exist(PathOrig,'file')
                xASL_Move(PathOrig, PathDest);
            end
        end
        
        PathM0 = fullfile(DirM0,'m0.nii');
        PathM0Dest = fullfile(ASLdir,'M0.nii');
        PathT1mapping_Dest = fullfile(ASLdir,'M0_T1mapping.nii');
        
        if xASL_exist(PathM0,'file')
            tIM = xASL_io_Nifti2Im(PathM0);
            if size(tIM,4)~=10
                warning(['M0 size differed for ' PathM0]);
            end
            xASL_io_SaveNifti(PathM0, PathM0Dest, tIM(:,:,:,end));
            clear tIM
            xASL_Move(PathM0, PathT1mapping_Dest);
        end
        
        % Delete unnecessary files
        xASL_delete(fullfile(DirM0,'T1_seg.nii'));
        xASL_delete(fullfile(DirM0,'brain.roi'));
        xASL_adm_DeleteFileList(DirM0, '^.*\.(png|jpg)$', false, [0 Inf]);
        
        if ~isempty(regexp(DirM0,'C3T-S006-20080205'))
            xASL_delete(fullfile(DirM0,'gm.nii'));
        elseif ~isempty(regexp(DirM0,'C3T-S036-20080709'))
            xASL_adm_DeleteFileList(DirM0, '^(pre|post)_(ASL|CBF|EPI).*\.nii$', false, [0 Inf]);
        end
            
        if ~isempty(xASL_adm_GetFileList(DirM0,'^.*$','FPListRec'))
            warning(['Still files in ' DirM0]);
        else
            xASL_delete(DirM0);
        end
            
        %% Move MRA MIP pictures
        DirMRA = fullfile(SubjectDir,'MRA');
        if exist(DirMRA,'dir')
            ScantypeRegExp = {'^.*MRA_(COR|SAG|TRA)\.(png|jpg)$' '^.*2DPC_(RL|AP)\.(png|jpg)$'};
            DestDirs = {fullfile(AnalysisSubjectDir,'MRA') fullfile(AnalysisSubjectDir,'PC')};
            IndexIs = [6 4];
            
            for iScantype=1:length(ScantypeRegExp)
                FileList = xASL_adm_GetFileList(DirMRA,ScantypeRegExp{iScantype},'List');
                % Get list of sequences
                if ~isempty(FileList)
                    SequenceList = cellfun(@(y) y(1:5), FileList, 'UniformOutput',false);
                    SequenceList = unique(SequenceList);

                    for iFile=1:length(FileList)
                        iSequence = find(cellfun(@(y) strcmp(y,FileList{iFile}(1:5)), SequenceList));
                        PathOrig = fullfile(DirMRA, FileList{iFile});
                        [~, Ffile, Fext] = fileparts(FileList{iFile});
                        DestDir = [DestDirs{iScantype} '_' num2str(iSequence)];
                        xASL_adm_CreateDir(DestDir);
                        
                        PathDest = fullfile(DestDir, [Ffile(end-IndexIs(iScantype):end), Fext]);
                        xASL_Move(PathOrig, PathDest);
                    end
                end
            end
                
            if ~isempty(xASL_adm_GetFileList(DirMRA,'^.*$','FPListRec'))
                warning(['Still files in ' DirMRA]);
            else
                xASL_delete(DirMRA);
            end
        end

        %% Move PCASL files
        DirACZ = fullfile(SubjectDir,'ACZ');
        if exist(DirACZ,'dir')
            DirsASL = xASL_adm_GetFileList(DirACZ, '^ACQ\d*$','FPList',[0 Inf], true);
            % First remove Figure subdir
            DirFig = fullfile(DirACZ,'FIG');
            if exist(DirFig,'dir')
                xASL_adm_DeleteFileList(DirFig, '^.*\.(png|jpg)$', true, [0 Inf]);
                if ~isempty(xASL_adm_GetFileList(DirFig,'^.*$','FPListRec'))
                    warning(['Still files in ' DirFig]);
                else
                    xASL_delete(DirFig);
                end
            end

            for iASL=1:length(DirsASL)
                for iFile=1:length(FileIn)
                    DirDest = [ASLdir(1:end-1) num2str(iASL)];
                    PathOrig = fullfile(DirsASL{iASL},FileIn{iFile});
                    PathDest = fullfile(DirDest,FileOut{iFile});
                    
                    if xASL_exist(PathOrig) && ~xASL_exist(PathDest)
                        xASL_adm_CreateDir(DirDest);
                        xASL_Move(PathOrig,PathDest);
                    elseif xASL_exist(PathOrig) && xASL_exist(PathDest)
                        xASL_delete(PathOrig);                        
                    end
                end
            end
            % Delete unnecessary files
            xASL_adm_DeleteFileList(DirACZ, '^.*\.(png|jpg)$', true, [0 Inf]);
            xASL_adm_DeleteFileList(DirACZ, '^(m|r|)(pre|post|)(_|)CBF(_|)(gm|wm|csf|seg|sn|inv_sn|).*\.(nii|nii\.gz|mat)$', true, [0 Inf]);
            xASL_adm_DeleteFileList(DirACZ, '^brain(|_2)\.nii$', false, [0 Inf]);
            DirsASL{end+1} = DirACZ;
            for iDir=1:length(DirsASL)
                if ~isempty(xASL_adm_GetFileList(DirsASL{iDir},'^.*$','FPListRec'))
                    warning(['Still files in ' DirsASL{iDir}]);
                else
                    xASL_delete(DirsASL{iDir});
                end
            end
        end

        %% Move VE-PCASL (vessel-encoded) files
        DirVEPCASL = fullfile(SubjectDir,'VEPCASL');
        if xASL_exist(DirVEPCASL,'dir')
            % Delete unnecessary files
            xASL_adm_DeleteFileList(DirVEPCASL, '^x.*\.(nii|nii\.gz|png|jpg|roi|mat)$', true, [0 Inf]); % delete mirrored files
            xASL_adm_DeleteFileList(DirVEPCASL, '^CBF_(gm|wm|csf|sn|inv_sn|seg)\.(mat|nii|nii\.gz|png|jpg)$', true, [0 Inf]); % delete CBF segmentations
            xASL_adm_DeleteFileList(DirVEPCASL, '^(brain|mask)(_2|)\.nii$', true, [0 Inf]); % delete brainmask
            
            DirsASL = xASL_adm_GetFileList(DirVEPCASL, '^ACQ\d*$','FPList',[0 Inf], true);
            
            for iASL=1:length(DirsASL)
                VEdir = fullfile(AnalysisSubjectDir,['VE_ASL_' num2str(iASL)]);
                xASL_adm_CreateDir(VEdir);
                xASL_adm_DeleteFileList(DirsASL{iASL}, '^.*\.(png|jpg)$', false, [0 Inf]); % delete ASL figures, but keep those from 4V/K3 reconstructions
                vepcaslFiles = xASL_adm_GetFileList(DirsASL{iASL}, '^vepcasl_\d_.*\.nii$','List');
                for iFile=1:length(vepcaslFiles) % remove double files
                    [~, Ffile, Fext] = xASL_fileparts(vepcaslFiles{iFile});
                    DoublePath = fullfile(DirsASL{iASL}, [Ffile(11:end) '.nii']);
                    if xASL_exist(DoublePath,'file')
                        xASL_delete(vepcaslFiles{iFile});
                    end
                end
                % Move filetypes of which we change the name
                for iFile=1:length(FileIn)
                    PathOrig = fullfile(DirsASL{iASL}, FileIn{iFile});
                    PathDest = fullfile(VEdir, FileOut{iFile});
                    if xASL_exist(PathOrig) && ~xASL_exist(PathDest)
                        xASL_Move(PathOrig, PathDest)
                    elseif xASL_exist(PathOrig) && xASL_exist(PathDest)
                        xASL_delete(PathOrig);
                    end
                end
                % Move other folders
                FoldersAre = {'4V' 'K3'};
                for iFolder=1:length(FoldersAre)
                    FolderIn = fullfile(DirsASL{iASL}, FoldersAre{iFolder});
                    FolderOut = fullfile(VEdir, FoldersAre{iFolder});
                    if exist(FolderIn,'dir')
                        xASL_Move(FolderIn, FolderOut);
                    end
                end
                % Delete residual files
                xASL_adm_DeleteFileList(DirsASL{iASL}, '.*vepcasl.*(CBF|RL|AP).*\.nii', false, [0 Inf]);
                
                % Delete unnecessary figure files
                DirFig = fullfile(DirsASL{iASL},'FIG');
                if exist(DirFig,'dir')
                    xASL_adm_DeleteFileList(DirFig, '^.*\.(png|jpg)$', true, [0 Inf]);
                    if ~isempty(xASL_adm_GetFileList(DirFig,'^.*$','FPListRec'))
                        warning(['Still files in ' DirFig]);
                    else
                        xASL_delete(DirFig);
                    end
                end                
                
                if ~isempty(xASL_adm_GetFileList(DirsASL{iASL},'^.*$','FPListRec'))
                    warning(['Still files in ' DirsASL{iASL}]);
                else
                    xASL_delete(DirsASL{iASL});
                end
            end

            % Delete unnecessary figure files
            DirFig = fullfile(DirVEPCASL,'FIG');
            if exist(DirFig,'dir')
                xASL_adm_DeleteFileList(DirFig, '^.*\.(png|jpg)$', true, [0 Inf]);
                if ~isempty(xASL_adm_GetFileList(DirFig,'^.*$','FPListRec'))
                    warning(['Still files in ' DirFig]);
                else
                    xASL_delete(DirFig);
                end
            end
            DirArchive = fullfile(DirVEPCASL,'archive');
            if exist(DirArchive,'dir')
                xASL_adm_DeleteFileList(DirArchive, '.*$', true, [0 Inf]);
                if ~isempty(xASL_adm_GetFileList(DirArchive,'^.*$','FPListRec'))
                    warning(['Still files in ' DirArchive]);
                else
                    xASL_delete(DirArchive);
                end
            end
            if ~isempty(regexp(DirVEPCASL,'C3T-S021-20080602'))
                xASL_delete(fullfile(DirVEPCASL,'j.roi'));
            end
            
            % Delete residual acquisition non-specific files
            xASL_adm_DeleteFileList(DirVEPCASL, '.*', false, [0 Inf]);
            
            if ~isempty(xASL_adm_GetFileList(DirVEPCASL,'^.*$','FPListRec'))
                warning(['Still files in ' DirVEPCASL]);
            else
                xASL_delete(DirVEPCASL);
            end
            
        %% Move DARTEL files
%         DirDARTEL = fullfile(SubjectDir,'DARTEL');
%         
%         
%         if exist(DirDARTEL,'dir')
%             xASL_adm_DeleteFileList(DirDARTEL, '^x.*\.(mat|nii|nii\.gz|png|jpg|roi)$', true, [0 Inf]); % delete x-files, are mirrored copies
%             DirCoreg = fullfile(DirDARTEL,'COREG');
%             if exist(DirCoreg,'dir')
%                 xASL_adm_DeleteFileList(DirCoreg, '^.*\.(mat|nii|nii\.gz|png|jpg|roi)$', false, [0 Inf]); % delete coregistration intermediate files
%                 if ~isempty(xASL_adm_GetFileList(DirCoreg,'^.*$','FPListRec'))
%                     warning(['Still files in ' DirCoreg]);
%                 else
%                     xASL_delete(DirCoreg);
%                 end
%             end
%             
%             DirSubject = fullfile(DirDARTEL,'SUBJECT');
%             if exist(DirSubject,'dir')
%                 xASL_adm_DeleteFileList(DirCoreg, '^(r|)(c4|c5|c6)T1\.nii$', false, [0 Inf]); % delete redundant intermediate files
%                 xASL_adm_DeleteFileList(DirCoreg, '^MSK\.nii$', false, [0 Inf]); % delete redundant intermediate files
%                 
%                 % Keep the minimum number of files to have correct registration
%                 % 
%             
            
        end % VEPCASL
    end % Subjects
end % Cohorts
fprintf('\n');

save(PathCohort,'Cohort');
save(PathScanDate,'ScanDate');

%% Delete subjects without M0/T1:
Subjects2Delete = {'C3T-S034_1' 'C3T-S086_1' 'C3T-S087_1'};
for iSubject=1:length(Subjects2Delete)
    SubjectDir = fullfile(AnalysisDir, Subjects2Delete{iSubject});
    xASL_adm_DeleteFileList(SubjectDir, '.*', true, [0 Inf]);
    xASL_delete(SubjectDir);
end

%% Replace ASL4D by CBF
% CBF is already quantified, saves lot of hassle, and more reproducible
% with previous paper (using same quantification)
for iSubject=1:x.nSubjects
    xASL_TrackProgress(iSubject,x.nSubjects);
    SubjectDir = fullfile(x.D.ROOT, x.SUBJECTS{iSubject});
    ASLdirs = xASL_adm_GetFileList(SubjectDir,'^ASL_\d$','FPList',[0 Inf],true);
    % Get ASL folders
    for iA=1:length(ASLdirs)
        PathASL4D = fullfile(ASLdirs{iA}, 'ASL4D.nii');
        
        xASL_delete(fullfile(ASLdirs{iA}, 'ASL4D.mat'));
        
        PathCBF = fullfile(ASLdirs{iA}, 'CBF.nii');
        if ~xASL_exist(PathCBF,'file')
            warning(['Non existing:' PathCBF]);
            fprintf('\n\n\n');
        else
            % Move CBF.nii -> ASL4D
            xASL_Move(PathCBF, PathASL4D, true, false);
        end
        PathM0 = fullfile(ASLdirs{iA}, 'M0.nii');
        PathM01 = fullfile(ASLdirs{1}, 'M0.nii');
        % same for M0 but copying
        if ~xASL_exist(PathM01,'file')
            warning(['Non existing:' PathM01]);
            fprintf('\n\n\n');
        else
            % Move CBF.nii -> ASL4D
            xASL_Copy(PathM01, PathM0, true, false);
        end        
    end
end

% -> Disable quantification in DataPar.json!

%% ExploreASL pipeline

% Can run CAT12 with low resolution T1, after upsampling it (set x.bFixResolution to true)
% CAT12 segmentation works nicely, but DARTEL to existing high resolution
% templates doesnt work -> Disable CAT12 DARTEL

OrigPath = fullfile(Rdir, 'DataParameters_HiQ.json');
DestPath = fullfile(AnalysisDir, 'DataParameters_HiQ.json');
xASL_Copy(OrigPath, DestPath, true);
cd(fileparts(which('ExploreASL_Master')));
% ExploreASL_Master(DestPath, 1, 1, 1, 1, 1); % but need to run DARTEL afterwards




% diamox failed at C3T-S048_
% warnings "too many folders in \\ANA" seem that M0/T1mapping were cloned
% or processed twice, while this was acquired once

% x = mirrored, just remove
% Processing skipped the second timepoint (which was post-CEA)




% -> CICERO1 seems to have failed scan first time point -> Nolan need to
% confirm
% Examples of weird registration standard space 069, 066, 053_2 (compare
% with 053_1 which worked)
% Lower resolution DARTEL didnt help