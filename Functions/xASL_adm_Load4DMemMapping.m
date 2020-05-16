function LoadFile = xASL_adm_Load4DMemMapping( x, WhichModality )
%xASL_adm_Load4DMemMapping Part of ExploreASL analysis module
% Loads data & maps it to memory mapping file on disc, if not done before

%% Admin, specify for each modality how to treat it
ModalitiesList  = {'ASL' 'ASL_notScaled' 'ASL_untreated' 'ASL_PVEc'  'ASL_HctCohort'    'ASL_HctCorrInd' 'PWI' 'M0' 'mean_control' 'SD' 'SNR' 'TT' 'PV_pGM' 'PV_pWM'  'FLAIR'   'T1'  'WMH_SEGM'  'c1T1' 'R1'};
PrefixList      = {''    'notScaled'     'untreated'     'GM_PVEC'   'HctCorr_cohort'   'HctCorrInd'     'PWI' 'M0' 'mean_control' 'SD' 'SNR' 'TT' 'PV_pGM' 'PV_pWM' 'rFLAIR' 'rT1' 'rWMH_SEGM'  'rc1T1' 'R1'};
qCBFprefix      = [1     1               1               1           1                  1                0     0    0              0    0     0    0        0        0        0      0           0       0];
NoSession       = [0     0               0               0           0                  0                0     0    0              0    0     1    1        1        1        1      1           1       1];

%% Which modality are we currently processing
for iC=1:length(ModalitiesList)
    if  strcmp(ModalitiesList{iC},WhichModality)
        iM=iC;
    end
end

%% What files should we load
if  qCBFprefix(iM)
    PreFix      = [x.P.CBF_Resliced '_' PrefixList{iM}];
else
    PreFix      = PrefixList{iM};
end

if  strcmp(PreFix(end),'_')
    PreFix      = PreFix(1:end-1);
end

AssumedSizeColumns  = sum(x.WBmask(:));

%% Load ASL, create 4D memory mapping file
% Create single 4D-file for memory mapping, if doesn't exist yet
LoadFile  = fullfile( x.D.PopDir, [ModalitiesList{iM} '_' num2str(x.nSubjects) '.dat']);

%% check if size of file matches; otherwise recreate it
if  exist( LoadFile ,'file' )
    x.S = dir(LoadFile);
    if  NoSession(iM)
        RequiredN   = x.nSubjects;
    else
        RequiredN   = x.nSubjectsSessions;
    end

    if  x.S.bytes ~= RequiredN * AssumedSizeColumns * 4
        warning('file size doesn''t match; recreating %s',LoadFile);
        delete(LoadFile);
    end
end

%% Start memory mapping
% NB: memory mapping fails if exclusions are not explicitly defined, all inclusions should be available
if  exist( LoadFile ,'file' )
    fprintf('%s\n',['Memory mapping ' WhichModality ' skipped, file already existed']);
    fprintf('%s\n','If re-mapping is desired, please delete this .dat-file first');

else
    if      NoSession(iM)
            nNeeded         = x.nSubjects;
            nSessions       = 1;
    elseif  NoSession(iM)==0
            nNeeded         = x.nSubjectsSessions;
            nSessions       = x.nSessions;
    else    error('Wrong NoSession definition');
    end


    if  length(xASL_adm_GetFileList(x.D.PopDir,['^' PreFix '_.*\.(nii|nii\.gz)$'],'FPList',[0 Inf]))==0
        fprintf('%s\n',['Memory mapping ' WhichModality ' skipped, no images exist']);
    else
        fprintf('\n');
        fprintf('%s\n','Pre-allocating memory: NB this can fail with large datasets and limited memory');
        tempData    = zeros( nNeeded, AssumedSizeColumns,'single'); % Pre-allocate memory
        fprintf(['Memory mapping ' WhichModality ' images...  ']);
        for iSubject=1:x.nSubjects
            xASL_TrackProgress(iSubject,x.nSubjects);

            for iSession=1:nSessions

                if  NoSession(iM)
                    SubjectSessionFile      = [PreFix '_' x.SUBJECTS{iSubject} '.nii'];
                else
                    SubjectSessionFile      = [PreFix '_' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession} '.nii'];
                end

                load_file                   = fullfile(x.D.PopDir, SubjectSessionFile);
                if  exist(load_file,'file') | exist([load_file '.gz'],'file')

                    iSubjSess                   = ((iSubject-1)*nSessions)+iSession;
                    tDat                        = xASL_io_Nifti2Im( load_file );
                    tempData(iSubjSess,:,:,:)   = xASL_im_IM2Column(tDat(:,:,:,1),x.WBmask);
                end
            end
        end

        fileID = fopen(LoadFile,'w');
        fwrite( fileID,tempData,'single');
        fclose(fileID);
    end
end


end
