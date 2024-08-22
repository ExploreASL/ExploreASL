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

fprintf('%s\n',['Adding ' x.S.output_ID ' statistics to participants.tsv...  ']);

bAnyWMHFound = false; % by default we don't add WMH volumes to participants.tsv below,
% only if these are found

for iSubject=1:x.dataset.nSubjects
    xASL_TrackProgress(iSubject,x.dataset.nSubjects);

    PathTSV = fullfile( x.D.TissueVolumeDir, ['TissueVolume_' x.SUBJECTS{iSubject} '.tsv']);
    if ~exist(PathTSV, 'file')
        % backward compatibility
        PathTSV = fullfile( x.D.TissueVolumeDir, ['TissueVolume_' x.SUBJECTS{iSubject} '.csv']);
    end
    
    % In case this subject had a FLAIR (in native space) then we default it
    % by NaNs, and issue a warning if no WMH volumetric statistics were
    % found
    pathFLAIR = fullfile(x.dir.xASLDerivatives, x.SUBJECTS{iSubject}, 'FLAIR.nii');
    pathWMH = fullfile(x.dir.xASLDerivatives, x.SUBJECTS{iSubject}, 'WMH_SEGM.nii');
    bHasWMH = xASL_exist(pathFLAIR, 'file') || xASL_exist(pathWMH, 'file');
    
    % SUBJECT & NaN (for absent/empty NIfTIs) definitions
    GM_vol{iSubject,1}          = x.SUBJECTS{iSubject};
    WM_vol{iSubject,1}          = x.SUBJECTS{iSubject};
    CSF_vol{iSubject,1}         = x.SUBJECTS{iSubject};
    GM_ICVRatio{iSubject,1}       = x.SUBJECTS{iSubject};
    GMWM_ICVRatio{iSubject,1}     = x.SUBJECTS{iSubject};
    WMH_vol{iSubject,1}         = x.SUBJECTS{iSubject};
    WMH_count{iSubject,1}          = x.SUBJECTS{iSubject};

    GM_vol{iSubject,2}          = NaN;
    WM_vol{iSubject,2}          = NaN;
    CSF_vol{iSubject,2}         = NaN;
    GM_ICVRatio{iSubject,2}       = NaN;
    GMWM_ICVRatio{iSubject,2}     = NaN;    
    WMH_vol{iSubject,2}         = NaN;
    WMH_count{iSubject,2}          = NaN;

    DidExist = 0;
    if exist(PathTSV, 'file')
        DidExist = 1;
        try
            [~, CellArray] = xASL_bids_csv2tsvReadWrite(PathTSV);
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
    PathTSV = xASL_adm_GetFileList(x.D.TissueVolumeDir, ['(?i)^WMH_LST_(LGA|LPA)_' x.SUBJECTS{iSubject} '(\.csv|\.tsv)$'], 'FPList', [0 Inf]);
    
    DidExist = 0;

    if ~isempty(PathTSV)
        DidExist = 1;
        try
            [~, CellArray] = xASL_bids_csv2tsvReadWrite(PathTSV{1});
            WMH_vol{iSubject,1} = x.SUBJECTS{iSubject};
            WMH_vol{iSubject,2} = CellArray{2,3};
            vol(iSubject,4) = xASL_str2num(WMH_vol{iSubject,2});
            
            WMH_count{iSubject,1} = x.SUBJECTS{iSubject};
            WMH_count{iSubject,2} = CellArray{2,4};
            vol(iSubject,5) = xASL_str2num(WMH_count{iSubject,2});
            DidExist = 2;
            
            bHasWMH = true; % even if no WMH existed in native space
            bAnyWMHFound = true; % to keep track that we want to add this to participants.tsv below
        catch ME
            fprintf('%s\n', ME.message);
        end
    end
    if bHasWMH && ~DidExist
        % when a native space FLAIR exists, but a standard space WMH does
        % not exist, we want to default and issue this message
        fprintf('\n%s\n', ['WMH volume for subject ' x.SUBJECTS{iSubject} ' was not found']);
        vol(iSubject,4:5) = NaN;
    elseif ~bHasWMH
        % when a native space FLAIR does not exists, we want to default
        % only
        vol(iSubject,4:5) = NaN;
    end    
end

fprintf('\n');



%% -----------------------------------------------------------------------------------------------
%% 3) Add stats in participants.tsv
VarName = {'GM_vol' 'WM_vol' 'CSF_vol' 'GM_ICVRatio' 'GMWM_ICVRatio'};
VarData = {GM_vol WM_vol CSF_vol GM_ICVRatio GMWM_ICVRatio};

if bAnyWMHFound
    VarName(end+1:end+2) = {'WMH_vol' 'WMH_count'};
    VarData(end+1:end+2) = {WMH_vol WMH_count};
else
    fprintf('%s\n', 'No WMH detected, so not reporting this in participants.tsv');
end
    
for iData=1:length(VarData)
    if exist(VarName{iData}, 'var')
        xASL_bids_Add2ParticipantsTSV(VarData{iData}, VarName{iData}, x);
    end
end


    
end