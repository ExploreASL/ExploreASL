function x = xASL_adm_CreateFileReport(x, bHasFLAIR, bHasMoCo, bHasM0, bHasLongitudinal)
%xASL_adm_CreateFileReport Prints a summary of created files or the
%individual modules (i.e. Structural, Longiutudinal & ASL modules)
% Provides a quick check to see what has been skipped, an whether all files
% are present.
%
% FORMAT:       x = xASL_adm_CreateFileReport(x, bHasFLAIR, bHasMoCo, bHasM0, bHasLongitudinal)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Prints a summary of created files or the individual modules
% (i.e. Structural, Longiutudinal & ASL modules). Provides a quick check to
% see what has been skipped, an whether all files are present.
%
% This script iterates across:
% Native space 1) subject and 2) session files,
% Resampled 3) subject and 4) session files,
% 5) Lock files and 6) QC Figure files.
%
% For all we perform a:
%
% - A) Count of the files present, summarized in FileReportSummary.csv
% - B) List of the missing files in "Missing*.csv" files
%
% PM: Simplify/optimize this code, to make filename variable changing,
% search within subject-directories, etc. Combine the parts searching for
% missing & summarizing count.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL





fclose all;
if nargin<1
    error('x missing from input, please input x table');
end
% Following booleans are to determine which files/locks need to be checked
% By default, we assume the minimal amount of image processing,
% to avoid confusion when all data are nicely processed but
% "missing_files.csv" are still showing up

fprintf('\n%s\n','====================================================================================');
fprintf('%s\n','Printing file reports... Check this for a summary of missing files');
fprintf('Detected');

if nargin<2 || isempty(bHasFLAIR)
    bHasFLAIR = false;
    iSubject = 1;
    while iSubject<=x.nSubjects && ~bHasFLAIR
        if xASL_exist(fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'FLAIR.nii'))
            bHasFLAIR = true;
        else
            iSubject = iSubject+1;
        end
    end
    if bHasFLAIR
        fprintf(': FLAIR ');
    else
        fprintf(': no FLAIR ');
    end
end
if nargin<3 || isempty(bHasMoCo)
    bHasMoCo = false;
    iSubject = 1;
    while iSubject<=x.nSubjects && ~bHasMoCo % speeds up instead of for-loop
        CurrentSubject = iSubject;
        for iSession=1:x.dataset.nSessions
            if xASL_exist(fullfile(x.D.ROOT, x.SUBJECTS{CurrentSubject}, x.SESSIONS{iSession}, 'ASL4D.mat'))
                bHasMoCo = true;
            else
                iSubject = iSubject+1;
            end
        end
    end
    if bHasMoCo
        fprintf(': MoCo results ');
    else
        fprintf(': no MoCo results ');
    end
end
if nargin<4 || isempty(bHasM0)
    bHasM0 = false;
    iSubject = 1;
    while iSubject<=x.nSubjects && ~bHasM0 % speeds up instead of for-loop
        CurrentSubject = iSubject;        
        for iSession=1:x.dataset.nSessions
            if xASL_exist(fullfile(x.D.ROOT, x.SUBJECTS{CurrentSubject}, x.SESSIONS{iSession}, 'M0.nii'))
                bHasM0 = true;
            else
                iSubject = iSubject+1;
            end
        end
    end
    if bHasM0
        fprintf(': M0 ');
    else
        fprintf(': no M0 ');
    end
end
if nargin<5 || isempty(bHasLongitudinal)
    bHasLongitudinal = false; %% THIS LONGITUDINAL PART HASNT BEEN AUTOMATED YET
end

fprintf('\n');

%% Create list of longitudinal registration subjects
LongRegSubj = '';
if ~isfield(x.S, 'SetsName')
    warning('x.S.SetsName missing, cannot define LongitudinalTimepoints');
    LongRegSubj = x.SUBJECTS;
else
    for iSet=1:length(x.S.SetsName)
        if strcmp(x.S.SetsName{iSet},'LongitudinalTimePoint')  % This is the longitudinal registration set
            for iS=1:x.nSubjects-1
                if x.S.SetsID(iS,iSet)==1 && x.S.SetsID(iS+1,iSet)==2
                    % This volume/TimePoint is the first out of multiple volumes/TimePoints
                    LongRegSubj{end+1} = x.SUBJECTS{iS};
                end
            end
        end
    end
end

%% -----------------------------------------------------------------------------
%% Administration SaveFiles
FileTypes = {'Native_subject' 'Native_session' 'MNI_subject' 'MNI_session' 'Lock' 'Figure'};

if bHasLongitudinal
    FileTypes{end+1} = 'LongReg';
end

LockDir = fullfile(x.D.ROOT, 'lock');
ReportName = fullfile(x.D.ROOT, 'FileReportSummary.csv');
xASL_delete(ReportName);
SummaryFid = fopen(ReportName,'wt');

xASL_adm_DeleteFileList(x.D.ROOT, '^Missing_.*_files\.csv$', false, [0 Inf]);
for iScanType=1:length(FileTypes)
    FileMissing{iScanType} = fullfile(x.D.ROOT, ['Missing_' FileTypes{iScanType} '_files.csv']); 
    SummaryFid_{iScanType} = fopen(FileMissing{iScanType},'wt');
end

CountMissing = [0 0 0 0 0 0 0];

%% -----------------------------------------------------------------------------
%% Define where to search for
% 1) Native space subject level
NativeRegExp_SubjectLevel = {[x.P.STRUCT '.nii'],['c1' x.P.STRUCT '.nii'],['c2' x.P.STRUCT '.nii']};
% 2) Native space session level
NativeRegExp_SessionLevel = {[x.P.ASL4D '.nii']};
% 3) MNI subject level
MNI_subject_prefix = {['r' x.P.STRUCT],['rc1' x.P.STRUCT],['rc2' x.P.STRUCT]};
% 4) MNI session level
MNI_session_prefix = {['q' x.P.CBF] 'PWI' 'SliceGradient'};
% 5) LongReg


% 5) Lock dirs
lockDIRS = {'xASL_module_Structural', 'xASL_module_ASL'};
lockPrefix{1} = {'010_LinearReg_T1w2MNI', '060_Segment_T1w', '080_Resample2StandardSpace', '090_GetVolumetrics', '100_VisualQC_Structural'}; % struct % struct
lockPrefix{2} = {'030_RegisterASL', '040_ResampleASL', '050_PreparePV', '070_CreateAnalysisMask', '080_Quantification', '090_VisualQC_ASL'}; % ASL


% Add dirs for specific image processing
if bHasFLAIR
    NativeRegExp_SubjectLevel(end+1:end+2) = {[x.P.FLAIR '.nii'],[x.P.WMH_SEGM '.nii']};
    MNI_subject_prefix(end+1:end+2) = {['r' x.P.FLAIR],['r' x.P.WMH_SEGM]};
    lockPrefix{1}(end+1:end+4) = {'020_LinearReg_FLAIR2T1w', '040_LST_Segment_FLAIR_WMH', '050_LST_T1w_LesionFilling_WMH','070_CleanUpWMH_SEGM'};
    % '030_FLAIR_BiasfieldCorrection' was taken out, for when external
    % WMH_SEGM.nii were used
end
if bHasMoCo
    lockPrefix{2}{end+1} = '020_RealignASL';
end
if bHasM0
    lockPrefix{2}{end+1} = '060_ProcessM0';
end
if bHasLongitudinal
    lockDIRS{end+1} = ['xASL_module_LongReg_' x.P.STRUCT];
    RegExp_LongReg = {['LongReg_y_' x.P.STRUCT '.nii']};
    lockPrefix{3} = {'001_LongitudinalRegistration','002_visualize','003_Resample','004_Re-visualize_Structural'}; % LongReg
end

%% -----------------------------------------------------------------------------
%% 6) Define figure files, where to search for
% FigureDIR{1}        = {x.D.T1CheckDir x.D.TissueVolumeDir x.D.CoregDir x.D.FLAIR_CheckDir x.D.FLAIR_REGDIR x.D.FlowFieldCheck};
% 
% FigPrefix{1}{1}     = {['r' x.P.STRUCT] ['r' x.P.STRUCT]};
% FigPrefix{1}{2}     = {['cat_' x.P.STRUCT] ['catROI_' x.P.STRUCT] 'Tissue_volume'};
% FigPrefix{1}{3}     = {['r' x.P.STRUCT '_ORI']};
% FigPrefix{1}{4}     = {['r' x.P.FLAIR]};
% FigPrefix{1}{5}     = {['r' x.P.FLAIR] 'segm_corr'};
% FigPrefix{1}{6}     = {['y_' x.P.STRUCT]};
% 
% FigEndfix{1}{1}     = {'.jpg' '_reg.jpg'}; % stackoverflow.com/questions/406230/regular-expression-to-match-a-line-that-doesnt-contain-a-word  ->  (?!reg) 
% FigEndfix{1}{2}     = {'.mat' '.mat' '.csv'};
% FigEndfix{1}{3}     = {'.jpg'};
% FigEndfix{1}{4}     = {'.jpg' '_reg.jpg'}; 
% FigEndfix{1}{5}     = {'_reg\.jpg' '.jpg'};
% FigEndfix{1}{6}     = {'.niiColor.jpg'};
% 
% FigureDIR{2}        = {x.D.ASLCheckDir x.D.T1_ASLREGDIR x.D.SliceCheckDir};
% 
% FigPrefix{2}{1}    = {['q' x.P.CBF '_'] ['q' x.P.CBF '_untreated_']};
% FigPrefix{2}{2}    = {'PWI_'};
% FigPrefix{2}{3}    = {'SliceGradient_'};
% 
% FigEndfix{2}{1}    = {'.jpg' '.jpg'};
% FigEndfix{2}{2}    = {'_reg.jpg'};
% FigEndfix{2}{3}    = {'_reg.jpg'};
% 
% if bHasLongitudinal
%     FigureDIR{3}        = {x.D.LongRegCheckDir};
%     FigPrefix{3}{1}     = {['LongReg_y_' x.P.STRUCT '_']};
%     FigEndfix{3}{1}     = {'.jpg'};
% end
    
%% Facultative ASL/session additions
% if time series exists (e.g. 3D GRASE and 2D EPI), then process them.
% FigureDIR{2}{end+1}         = x.SNRdir;
% FigPrefix{2}{end+1}         = {'SD_map' 'SNR_map'};
% FigEndfix{2}{end+1}         = {'.jpg' '.jpg'};
% MNI_session_prefix{end+1}   = 'SD';
% MNI_session_prefix{end+1}   = 'SNR';

% FigureDIR{2}{end+1}         = x.MotionDir;
% FigPrefix{2}{end+1}         = {'motion_correction_NDV' 'rp'};
% FigEndfix{2}{end+1}         = {'.mat' '_motion.jpg'};

% FigureDIR{2}{end+1}         = x.D.M0CheckDir;
% FigPrefix{2}{end+1}         = {x.P.M0};
% FigEndfix{2}{end+1}         = {'.jpg'};
% FigureDIR{2}{end+1}         = x.D.M0regASLdir;
% FigPrefix{2}{end+1}         = {x.P.M0};
% FigEndfix{2}{end+1}         = {'_reg.jpg'};
    
% MNI_session_prefix{end+1}           = x.P.M0;
% NativeRegExp_SessionLevel{end+1}    = [x.P.M0 '.nii'];
% NativeRegExp_SessionLevel{end+1}    = [x.P.M0 '_parms.mat'];
    
% FigureDIR{2}{end+1}     = x.D.TTCheckDir;
% FigPrefix{2}{end+1}     = {'TT'};
% FigEndfix{2}{end+1}     = {'\.jpg'};
            
fprintf('Checking native space: ');

%% -----------------------------------------------------------------------------
%% 1) Native space subject files (e.g. T1w, FLAIR)
fprintf('anat files:   ');
for iExp=1:length(NativeRegExp_SubjectLevel)
    xASL_TrackProgress(iExp,length(NativeRegExp_SubjectLevel));
    for iSubject=1:x.nSubjects
        FilePathNii    = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, NativeRegExp_SubjectLevel{iExp});
        if ~xASL_exist(FilePathNii,'file')
            [Fpath, Ffile] = xASL_fileparts(FilePathNii);
            RegExpBIDS = ['(.*' Ffile '.*run.*|.*run.*' Ffile '.*)\.(nii|nii\.gz)'];
            RegExpBIDS = xASL_adm_GetFileList(Fpath, RegExpBIDS, 'FPList', [0 Inf], false);
            
            if isempty(RegExpBIDS)
                fprintf(SummaryFid_{1},'%s\n', FilePathNii );
                CountMissing(1) = CountMissing(1)+1;
                NativeSubjectExists(iSubject,1) = false;
            else
                NativeSubjectExists(iSubject,1) = true;
            end
        else
            NativeSubjectExists(iSubject,1) = true;
        end
    end
    fprintf(SummaryFid,'%s\n', [num2str(CountMissing(1)) ' ' NativeRegExp_SubjectLevel{iExp} ' missing']);
end
fprintf(', ');


%% -----------------------------------------------------------------------------
%% 2) Native space session files (e.g. ASL4D.nii, M0)
fprintf('ASL files:   ');
for iExp=1:length(NativeRegExp_SessionLevel)
    xASL_TrackProgress(iExp,length(NativeRegExp_SessionLevel));
    for iSubject=1:x.nSubjects
        for iSession=1:x.dataset.nSessions
            iSubjectSession = (iSubject-1)*x.dataset.nSessions+iSession;
            FilePathNii    = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, NativeRegExp_SessionLevel{iExp});

            if ~xASL_exist(FilePathNii,'file')
                [Fpath, Ffile] = xASL_fileparts(FilePathNii);
                RegExpBIDS = ['(.*' Ffile '.*run.*|.*run.*' Ffile '.*)\.(nii|nii\.gz)'];
                RegExpBIDS = xASL_adm_GetFileList(Fpath, RegExpBIDS, 'FPList', [0 Inf], false);

                if isempty(RegExpBIDS)                
                    fprintf(SummaryFid_{2},'%s\n', FilePathNii );
                    CountMissing(2)     = CountMissing(2)+1;
                    NativeSessionExist(iSubjectSession,1) = false;
                else
                    NativeSessionExist(iSubjectSession,1) = true;
                end
            else
                NativeSessionExist(iSubjectSession,1) = true;
            end     
        end
    end
    fprintf(SummaryFid,'%s\n', [num2str(CountMissing(2)) ' ' NativeRegExp_SessionLevel{iExp} ' missing']);
end

fprintf('\nChecking common space: ');

%% -----------------------------------------------------------------------------
%% 3) Common space subject files (e.g. rT1, rFLAIR)
fprintf('anat files:   ');
for iPrefix=1:length(MNI_subject_prefix)
    xASL_TrackProgress(iPrefix,length(MNI_subject_prefix));
    for iSubject=1:x.nSubjects
        if NativeSubjectExists(iSubject,1)
            FilePathNii    = fullfile(x.D.PopDir, [MNI_subject_prefix{iPrefix} '_' x.SUBJECTS{iSubject} '.nii']);

            if ~xASL_exist(FilePathNii,'file')
                fprintf(SummaryFid_{3},'%s\n', FilePathNii );
                CountMissing(3) = CountMissing(3)+1;
            end
        end
    end
    fprintf(SummaryFid,'%s\n', [num2str(CountMissing(3)) ' ' MNI_subject_prefix{iExp} ' missing in MNI']);
end
fprintf(', ');

%% -----------------------------------------------------------------------------
%% 4) Common space session files (e.g. qCBF)
fprintf('ASL files:   ');
for iPrefix=1:length(MNI_session_prefix)
    xASL_TrackProgress(iPrefix,length(MNI_session_prefix));
    for iSubject=1:x.nSubjects
        for iSession=1:x.dataset.nSessions
            iSubjectSession = (iSubject-1)*x.dataset.nSessions+iSession;
            if NativeSessionExist(iSubjectSession,1)
                FilePathNii    = fullfile(x.D.PopDir, [MNI_session_prefix{iPrefix} '_' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession} '.nii']);

                if ~xASL_exist(FilePathNii,'file')
                    fprintf(SummaryFid_{4},'%s\n', FilePathNii );
                    CountMissing(4)     = CountMissing(4)+1;
                end
            end
        end
    end
    fprintf(SummaryFid,'%s\n', [num2str(CountMissing(4)) ' ' MNI_session_prefix{iExp} ' missing in MNI']);
end
fprintf(', checking ');

%% -----------------------------------------------------------------------------
%% 5) Lock files
fprintf('lock files:   ');
for iDIR=1:length(lockDIRS)
    xASL_TrackProgress(iDIR,length(lockDIRS));
    for iP=1:length(lockPrefix{iDIR})
        if strcmp(lockDIRS{iDIR},['xASL_module_DARTEL_' x.P.STRUCT]) || strcmp(lockDIRS{iDIR},'Population')
            % only 1 status file
                FilePath    = fullfile(LockDir, lockDIRS{iDIR}, [lockDIRS{iDIR}],[lockPrefix{iDIR}{iP} '.status']);
                if ~exist(FilePath,'file')
                    fprintf(SummaryFid_{5},'%s\n', FilePath);
                    CountMissing(5)     = CountMissing(5)+1;
                end
            
        elseif strcmp(lockDIRS{iDIR},['xASL_module_LongReg_' x.P.STRUCT]) % status files for Longitudinal Registration
                for iSubj=1:length(LongRegSubj)
                    FilePath = fullfile(LockDir, lockDIRS{iDIR}, LongRegSubj{iSubj}, 'xASL_module_LongReg',[lockPrefix{iDIR}{iP} '.status']);
                    if ~exist(FilePath,'file') % check subjects
                        fprintf(SummaryFid_{5},'%s\n', FilePath);
                        CountMissing(5) = CountMissing(5)+1;
                    end                    
                end
        else % status file for each Subject/Session
            for iSubject=1:x.nSubjects
                if  strcmp(lockDIRS{iDIR},'xASL_module_ASL')
                    for iSession=1:x.dataset.nSessions % check sessions
                        FilePath = fullfile(LockDir, lockDIRS{iDIR}, x.SUBJECTS{iSubject}, [lockDIRS{iDIR} '_' x.SESSIONS{iSession}],[lockPrefix{iDIR}{iP} '.status']);
                        if ~exist(FilePath,'file')
                            fprintf(SummaryFid_{5},'%s\n', FilePath);
                            CountMissing(5) = CountMissing(5)+1;
                        end
                    end

                else
                    FilePath = fullfile(LockDir, lockDIRS{iDIR}, x.SUBJECTS{iSubject}, [lockDIRS{iDIR}],[lockPrefix{iDIR}{iP} '.status']);
                    if ~exist(FilePath,'file') % check subjects
                        fprintf(SummaryFid_{5},'%s\n', FilePath);
                        CountMissing(5) = CountMissing(5)+1;
                    end
                end
            end
        end
    end
end
fprintf(SummaryFid,'%s\n', [num2str(CountMissing) ' status/lock files missing']);

fprintf('\n');

%% -----------------------------------------------------------------------------
%% 6) Figure files (total & missing files)
% fprintf('%s\n','Checking Figure files...  ');
% for iT=1:3 % struct, asl (subjects & sessions) & LongReg
%     for iDir=1:length(FigureDIR{iT})
%         xASL_TrackProgress(iDir,length(FigureDIR{iT}));
%         for iP=1:length(FigPrefix{iT}{iDir})
% 
%             if  strcmp(FigureDIR{3},x.D.LongRegCheckDir) % LongReg check data
%                 for iSubj=1:length(LongRegSubj)
%                     for iS=1:x.nSubjects
%                         if  strcmp(LongRegSubj{iSubj},x.SUBJECTS{iS}) && iS<x.nSubjects
%                             LRCompareName   = [LongRegSubj{iSubj} '_vs_' x.SUBJECTS{iS+1}];
%                         end
%                     end                    
%                     
%                     FilePath    = fullfile(FigureDIR{iT}{iDir}, [FigPrefix{iT}{iDir}{iP} LRCompareName FigEndfix{iT}{iDir}{iP}]);
%                     if ~exist(FilePath,'file')
%                         fprintf(SummaryFid_{6},'%s\n', FilePath );
%                         CountMissing(6)     = CountMissing(6)+1;
%                     end   
%                 end
%             else
% 
%                 for iSubject=1:x.nSubjects
%                     for iSession=1:x.dataset.nSessions
% 
%                         if  iT==2 % if it concerns sessions
%                             SearchFile  = [x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession}];
%                         else
%                             SearchFile  = x.SUBJECTS{iSubject};
%                         end
% 
%                         if ~(iT==1 && iSession>1) % skip if subject & session other than 1
% 
%                             % B) list missing files
%                             FileSearch  = [FigPrefix{iT}{iDir}{iP} '_' SearchFile FigEndfix{iT}{iDir}{iP}];
%                             FilePath = fullfile(FigureDIR{iT}{iDir}, FileSearch);
%                             if ~exist(FilePath,'file')
%                                 fprintf(SummaryFid_{6},'%s\n', FilePath);
%                                 CountMissing(6)     = CountMissing(6)+1;
%                             end
%                         end
% 
%                     end
%                 end
%             end % if  strcmp(FigureDIR{3},x.D.LongRegCheckDir) % LongReg check data
%         end % for iP=1:length(FigPrefix{iT}{iDir})
%     end % for iDir=1:length(FigureDIR{iT})
% end
% fprintf('\n');

%% -----------------------------------------------------------------------------
%% 7) LongReg files (i.e. LongReg_y_T1.nii)
if bHasLongitudinal
    fprintf('%s\n','Counting longitudinal registration files:  ');
    for iExp=1:length(RegExp_LongReg)
        xASL_TrackProgress(iExp,length(RegExp_LongReg));
        for iSubj=1:length(LongRegSubj) % this will skip if there are no longitudinal scans
            FilePath        = fullfile(x.D.ROOT, LongRegSubj{iSubj}, RegExp_LongReg{iExp});
            if ~exist(FilePath,'file')
                fprintf(SummaryFid_{7},'%s\n', FilePath );
                CountMissing(7)     = CountMissing(7)+1;
            end
            fprintf(SummaryFid,'%s\n', [num2str(CountMissing(7)) ' ' RegExp_LongReg{iExp} ' missing']);
        end
    end
    fprintf('\n');
end

%% -----------------------------------------------------------------------------
%% Housekeeping

fclose(SummaryFid);
for iScanType=1:length(SummaryFid_)
    fclose(SummaryFid_{iScanType});
    if  CountMissing(iScanType)==0; delete(FileMissing{iScanType}); % fprintf('%s\n',['No missing ' FileTypes{iM} ' files']);
    else
        fprintf('%s\n',['Please check missing ' FileTypes{iScanType} ' files: ' FileMissing{iScanType}])
    end
end

fprintf('%s\n','Done generating file reports');
fprintf('%s\n','====================================================================================');

end