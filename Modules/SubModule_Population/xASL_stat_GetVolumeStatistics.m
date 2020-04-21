function xASL_stat_GetVolumeStatistics( x, GetStats)
%xASL_stat_GetVolumeStatistics Summarize volume values
% INPUT
% x = from ExploreASL
% masks   = physically loaded masks for each subject
% ASL     = the data to be analyzed. Could be ASL, or e.g. SD or SNR masks
%
%
% By HJMM Mutsaerts, ExploreASL 2019
%



%% -----------------------------------------------------------------------------------------------
%% Load volume files

if ~exist('GetStats','var')
    GetStats    = 0;
end

x.S.NamesROI                              = {'GM_vol' 'WM_vol' 'CSF_vol' 'WMH_vol' 'WMH_count'}; % 'Bone (L)' 'Soft tissue (L)' 'Air (L)'
x.S.output_ID                             = 'volume';
x.S.unit                                  = 'L';

fprintf('%s\n',['Printing csv-files with ' x.S.output_ID ' statistics...  ']);




%% -----------------------------------------------------------------------------------------------
for iSubject=1:x.nSubjects
    xASL_TrackProgress(iSubject,x.nSubjects);

    x.S.SUBJECTID{iSubject,1}   = x.SUBJECTS{iSubject};

    csv_load                    = fullfile( x.D.TissueVolumeDir, ['TissueVolume_' x.SUBJECTS{iSubject} '.csv']);
    if ~exist(csv_load,'file')
        csv_load                = fullfile( x.D.TissueVolumeDir, ['TissueVolume_' x.SUBJECTS{iSubject} '.tsv']);
    end
    
    GM_vol{iSubject,1}          = x.SUBJECTS{iSubject};
    WM_vol{iSubject,1}          = x.SUBJECTS{iSubject};
    CSF_vol{iSubject,1}         = x.SUBJECTS{iSubject};
    GM_ICVRatio{iSubject,1}     = x.SUBJECTS{iSubject};
    GMWM_ICVRatio{iSubject,1}   = x.SUBJECTS{iSubject};
    WMH_vol{iSubject,1}         = x.SUBJECTS{iSubject};
    WMH_count{iSubject,1}       = x.SUBJECTS{iSubject};      

    GM_vol{iSubject,2}          = NaN;
    WM_vol{iSubject,2}          = NaN;
    CSF_vol{iSubject,2}         = NaN;
    GM_ICVRatio{iSubject,2}     = NaN;
    GMWM_ICVRatio{iSubject,2}   = NaN;
    WMH_vol{iSubject,2}         = NaN;
    WMH_count{iSubject,2}       = NaN;          

    DidExist = 0;
    if exist(csv_load,'file')
        DidExist = 1;
        try
            [~, TempCell] = xASL_adm_csv2tsv(csv_load);
            for iMeas=1:size(TempCell,2)-1
                vol(iSubject,iMeas) = str2num(TempCell{2,iMeas+1});
            end

            GM_vol{iSubject,2} = vol(iSubject,1);
            WM_vol{iSubject,2} = vol(iSubject,2);
            CSF_vol{iSubject,2} = vol(iSubject,3);
            GM_ICVRatio{iSubject,2} = vol(iSubject,1)/sum(vol(iSubject,1:3));
            GMWM_ICVRatio{iSubject,2} = sum(vol(iSubject,1:2))/sum(vol(iSubject,1:3));
            DidExist = 2;
        catch ME
            fprintf('%s\n',ME.message);
        end
    end
    if DidExist==0
        fprintf('%s\n',['Tissue volume for subject ' x.SUBJECTS{iSubject} ' was not found   ']);
        vol(iSubject,1:3) = NaN;
    elseif DidExist==1
        fprintf('%s\n',['Tissue volume for subject ' x.SUBJECTS{iSubject} ' could not be read   ']);
        vol(iSubject,1:3) = NaN;
    end

    %% WMH
    csv_load = xASL_adm_GetFileList( x.D.TissueVolumeDir, ['^WMH_LST_(LGA|LPA)_' x.SUBJECTS{iSubject} '(\.csv|\.tsv)$'], 'FPList', [0 Inf]);
    
    DidExist = 0;
    if ~isempty(csv_load)
        DidExist = 1;
        try
            [~, TempCell] = xASL_adm_csv2tsv(csv_load{1});

            WMH_vol{iSubject,2} = TempCell{2,4};
            vol(iSubject,4) = str2num(WMH_vol{iSubject,2});
            WMH_count{iSubject,2} = TempCell{2,5};
            vol(iSubject,5) = str2num(WMH_count{iSubject,2});
            DidExist = 2;
        catch ME
            fprintf('%s\n',ME.message);
        end
    end
    if DidExist==0
        fprintf('%s\n',['WMH volume for subject ' x.SUBJECTS{iSubject} ' was not found   ']);
        vol(iSubject,4:5) = NaN;
    elseif DidExist==1
        fprintf('%s\n',['WMH volume for subject ' x.SUBJECTS{iSubject} ' could not be read   ']);
        vol(iSubject,4:5) = NaN;
    end    
end

fprintf('\n');



%% -----------------------------------------------------------------------------------------------
% Save mat-file for stats
Var2Save                            = {'GM_vol' 'WM_vol' 'CSF_vol' 'GM_ICVRatio' 'GMWM_ICVRatio'};

if  exist('WMH_vol','var') && exist('WMH_count','var')
    Var2Save{end+1}     = 'WMH_vol';
    Var2Save{end+1}     = 'WMH_count';
end

for iV=1:length(Var2Save)
    MatFile{iV}                     = fullfile( x.D.ROOT, [Var2Save{iV} '.mat']);
    save(MatFile{iV},Var2Save{iV});
end






%% -----------------------------------------------------------------------------------------------
if  GetStats
    % Restructure DATA per ROI/measurement
    x.S.DAT       = vol;

    % Restructure x.S.SetsID into single session only
    if  isfield(x.S,'SetsID')
        if  numel(x.S.SetsID)>0
            TempID  = x.S.SetsID;
            x.S       = rmfield(x.S,'SetsID');
            for iSubject=1:x.nSubjects
                x.S.SetsID(iSubject,:)    = TempID( ((iSubject-1)*x.nSessions)+1,:);
            end
        end
    end

    % Initiation parameters for xASL_wrp_PermuteSets1

    x.nSessions   = 1;
    x.S.KISS              = 1; % keeps it simple
    x.S.StatsDir          = x.S.StatsDir;
    x.S.function2call     = @xASL_stat_PrintBasicStats;

    xASL_wrp_PermuteSets1( x );
end
    
    
    
end
