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
% DESCRIPTION: This functions collects motion stats, with the following. Steps:
%
% 1. Collect structural volume data
% 2. Collect WMH data
% 3. Add stats in participants.tsv
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_stat_GetVolumeStatistics(x);
% __________________________________
% Copyright 2016-2021 ExploreASL


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
    
    % In case this subject had a FLAIR (in native space) then we default it
    % by NaNs, and issue a warning if no WMH volumetric statistics were
    % found
    Path_FLAIR = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'FLAIR.nii');
    Path_WMH = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'WMH_SEGM.nii');
    if xASL_exist(Path_FLAIR, 'file') || xASL_exist(Path_WMH, 'file')
        HasWMH = true;
    else
        % skip FLAIR volumetrics for this subject
        HasWMH = false;
    end
    
    % SUBJECT & NaN (for absent/empty NIfTIs) definitions
    GM_vol{iSubject,1}          = x.SUBJECTS{iSubject};
    WM_vol{iSubject,1}          = x.SUBJECTS{iSubject};
    CSF_vol{iSubject,1}         = x.SUBJECTS{iSubject};
    GM_ICVRatio{iSubject,1}     = x.SUBJECTS{iSubject};
    GMWM_ICVRatio{iSubject,1}   = x.SUBJECTS{iSubject};

    GM_vol{iSubject,2}          = NaN;
    WM_vol{iSubject,2}          = NaN;
    CSF_vol{iSubject,2}         = NaN;
    GM_ICVRatio{iSubject,2}     = NaN;
    GMWM_ICVRatio{iSubject,2}   = NaN;    
    
    if HasWMH
        WMH_vol{iSubject,1}         = x.SUBJECTS{iSubject};
        WMH_count{iSubject,1}       = x.SUBJECTS{iSubject};

        WMH_vol{iSubject,2}         = NaN;
        WMH_count{iSubject,2}       = NaN;
    end

    DidExist = 0;
    if exist(PathCSV,'file')
        DidExist = 1;
        try
            [~, CellArray] = xASL_bids_csv2tsvReadWrite(PathCSV);
            for iMeas=1:size(CellArray,2)-1
                vol(iSubject,iMeas) = xASL_str2num(CellArray{2,iMeas+1});
            end

            GM_vol{iSubject,2} = vol(iSubject,1);
            WM_vol{iSubject,2} = vol(iSubject,2);
            CSF_vol{iSubject,2} = vol(iSubject,3);
            GM_ICVRatio{iSubject,2} = vol(iSubject,1)/sum(vol(iSubject,1:3));
            GMWM_ICVRatio{iSubject,2} = sum(vol(iSubject,1:2))/sum(vol(iSubject,1:3));
            DidExist = 2;
        catch ME
            fprintf('%s\n', ME.message);
        end
    end
    if DidExist==0
        fprintf('%s\n', ['Tissue volume for subject ' x.SUBJECTS{iSubject} ' was not found']);
        vol(iSubject,1:3) = NaN;
    elseif DidExist==1
        fprintf('%s\n', ['Tissue volume for subject ' x.SUBJECTS{iSubject} ' could not be read']);
        vol(iSubject,1:3) = NaN;
    end

    %% -----------------------------------------------------------------------------------------------
    %% 2) Collect WMH data
    PathCSV = xASL_adm_GetFileList(x.D.TissueVolumeDir, ['(?i)^WMH_LST_(LGA|LPA)_' x.SUBJECTS{iSubject} '(\.csv|\.tsv)$'], 'FPList', [0 Inf]);
    
    DidExist = 0;
    if ~isempty(PathCSV)
        DidExist = 1;
        try
            [~, CellArray] = xASL_bids_csv2tsvReadWrite(PathCSV{1});
            WMH_vol{iSubject,1} = x.SUBJECTS{iSubject};
            WMH_vol{iSubject,2} = CellArray{2,4};
            vol(iSubject,4) = xASL_str2num(WMH_vol{iSubject,2});
            
            WMH_count{iSubject,1} = x.SUBJECTS{iSubject};
            WMH_count{iSubject,2} = CellArray{2,5};
            vol(iSubject,5) = xASL_str2num(WMH_count{iSubject,2});
            DidExist = 2;
            
            HasWMH = true; % even if no WMH existed in native space
        catch ME
            fprintf('%s\n', ME.message);
        end
    end
    if HasWMH && DidExist==0
        fprintf('%s\n', ['WMH volume for subject ' x.SUBJECTS{iSubject} ' was not found']);
        vol(iSubject,4:5) = NaN;
    elseif HasWMH && DidExist==1
        fprintf('%s\n', ['WMH volume for subject ' x.SUBJECTS{iSubject} ' could not be read']);
        vol(iSubject,4:5) = NaN;
    end    
end

fprintf('\n');



%% -----------------------------------------------------------------------------------------------
%% 3) Add stats in participants.tsv
VarName = {'GM_vol' 'WM_vol' 'CSF_vol' 'GM_ICVRatio' 'GMWM_ICVRatio'};
VarData = {GM_vol WM_vol CSF_vol GM_ICVRatio GMWM_ICVRatio};

if HasWMH
    VarName(end+1:end+2) = {'WMH_vol' 'WMH_count'};
    VarData(end+1:end+2) = {WMH_vol WMH_count};
end
    
for iData=1:length(VarData)
    if exist(VarName{iData}, 'var')
        xASL_bids_Add2ParticipantsTSV(VarData{iData}, VarName{iData}, x);
    end
end


    
end
