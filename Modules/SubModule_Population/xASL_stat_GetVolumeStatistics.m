function xASL_stat_GetVolumeStatistics(x)
%xASL_stat_GetVolumeStatistics Summarize volume values
%
% FORMAT: xASL_stat_GetVolumeStatistics(x)
% 
% INPUT:
%   x - struct containing pipeline environment parameters (REQUIRED)
%
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects motion stats, with the following
% steps:
% 1) Collect structural volume data
% 2) Collect WMH data
% 3) Add stats in participants.tsv
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_stat_GetVolumeStatistics(x);
% __________________________________
% Copyright 2016-2020 ExploreASL


%% -----------------------------------------------------------------------------------------------
%% 1) Collect structural volume data
x.S.NamesROI = {'GM_vol' 'WM_vol' 'CSF_vol' 'WMH_vol' 'WMH_count'}; % 'Bone (L)' 'Soft tissue (L)' 'Air (L)'
x.S.output_ID = 'volume';
x.S.unit = 'L';

fprintf('%s\n',['Printing csv-files with ' x.S.output_ID ' statistics...  ']);


for iSubject=1:x.nSubjects
    xASL_TrackProgress(iSubject,x.nSubjects);

    PathCSV = fullfile( x.D.TissueVolumeDir, ['TissueVolume_' x.SUBJECTS{iSubject} '.csv']);
    if ~exist(PathCSV,'file')
        PathCSV = fullfile( x.D.TissueVolumeDir, ['TissueVolume_' x.SUBJECTS{iSubject} '.tsv']);
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
    if exist(PathCSV,'file')
        DidExist = 1;
        try
            [~, CellArray] = xASL_bids_csv2tsvReadWrite(PathCSV);
            for iMeas=1:size(CellArray,2)-1
                vol(iSubject,iMeas) = str2num(CellArray{2,iMeas+1});
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
        fprintf('%s\n',['Tissue volume for subject ' x.SUBJECTS{iSubject} ' was not found   \n']);
        vol(iSubject,1:3) = NaN;
    elseif DidExist==1
        fprintf('%s\n',['Tissue volume for subject ' x.SUBJECTS{iSubject} ' could not be read   \n']);
        vol(iSubject,1:3) = NaN;
    end

    %% -----------------------------------------------------------------------------------------------
    %% 2) Collect WMH data
    PathCSV = xASL_adm_GetFileList(x.D.TissueVolumeDir, ['^WMH_LST_(LGA|LPA)_' x.SUBJECTS{iSubject} '(\.csv|\.tsv)$'], 'FPList', [0 Inf]);
    
    DidExist = 0;
    if ~isempty(PathCSV)
        DidExist = 1;
        try
            [~, CellArray] = xASL_bids_csv2tsvReadWrite(PathCSV{1});

            WMH_vol{iSubject,2} = CellArray{2,4};
            vol(iSubject,4) = str2num(WMH_vol{iSubject,2});
            WMH_count{iSubject,2} = CellArray{2,5};
            vol(iSubject,5) = str2num(WMH_count{iSubject,2});
            DidExist = 2;
        catch ME
            fprintf('%s\n',ME.message);
        end
    end
    if DidExist==0
        fprintf('%s\n',['WMH volume for subject ' x.SUBJECTS{iSubject} ' was not found']);
        vol(iSubject,4:5) = NaN;
    elseif DidExist==1
        fprintf('%s\n',['WMH volume for subject ' x.SUBJECTS{iSubject} ' could not be read']);
        vol(iSubject,4:5) = NaN;
    end    
end

fprintf('\n');



%% -----------------------------------------------------------------------------------------------
%% 3) Add stats in participants.tsv
VarName = {'GM_vol' 'WM_vol' 'CSF_vol' 'GM_ICVRatio' 'GMWM_ICVRatio' 'WMH_vol' 'WMH_count'};
VarData = {GM_vol WM_vol CSF_vol GM_ICVRatio GMWM_ICVRatio WMH_vol WMH_count};

for iData=1:length(VarData)
    if exist(VarName{iData}, 'var')
        xASL_bids_Add2ParticipantsTSV(VarData{iData}, VarName{iData}, x);
    end
end


    
end
