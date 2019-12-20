function EPAD_CopyFLAIR_WMH_Carole(ROOT, CaroleDir)
%EPAD_ASL_parmsPrepare Prepare ASL parameter files for different EPAD vendors/sites
%
% FORMAT: EPAD_CopyFLAIR_WMH_Carole(ROOT, CaroleDir)
% 
% INPUT:
%   ROOT        - path to root folder containing the /analysis folder with NIfTI/BIDS data (REQUIRED)
%   CaroleDir   - path to root folder of Carole' segmentations & FLAIRs
%                 (OPTIONAL, trying a previous location if not provided)
%   
%
% OUTPUT: n/a
% OUTPUT FILE:
%    /ROOT/analysis/NotCopiedList.json - contains list of missing WMH
%                                        segmentations (where we couldnt find both the FLAIR/WMH_SEGM from
%                                        Carole)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function searches for Carole' FLAIRs & WMH_SEGM & copies them into the ExploreASL/BIDS compliant
%              directory structure, for further use in EPAD analyses & QC.
%              We copy the FLAIR in addition to the WMH, as the WMH
%              segmentation will be aligned with Carole' FLAIR but not
%              necessarily with ours. For ExploreASL, these files are
%              copied to:
%              /ROOT/analysis/SubjectName/FLAIR.nii[.gz]
%              /ROOT/analysis/SubjectName/WMH_SEGM.nii[.gz]
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_CopyFLAIR_WMH_Carole('data/RAD/share/EPAD500/analysis','data/RAD/share/epad-carole/clean');
% __________________________________
% Copyright 2015-2019 ExploreASL


if nargin<2
    cd(ROOT);
    try
        cd('../epad-carole/clean');
        CaroleDir=pwd;
    catch
        warning('Didnt find folder for Caroles scans, skipping');
        return;
    end
end

if ~exist(CaroleDir,'dir')
    warning('Didnt find folder for Caroles scans, skipping');
    return;
end

AnalysisDir = fullfile(ROOT, 'analysis');

SubjList = xASL_adm_GetFileList(AnalysisDir, '^\d{3}EPAD\d*$', 'FPList', [0 Inf], true);

NotCopiedList = struct;

fprintf('Copying FLAIR & WMH segmentations Carole:   ');
for iS=1:length(SubjList)
    xASL_TrackProgress(iS, length(SubjList));
    [~, NameNew] = fileparts(SubjList{iS});
    % here we check if the scanner prefix is a site prefix (e.g. 010 instead of 110)
    % & also if there is a '-' instead of 'EPAD'
    Name1 = ['0' NameNew(2:3) '-' NameNew(8:end)];
    Name2 = NameNew(8:end);
    
    FLAIRnew = fullfile(SubjList{iS}, 'FLAIR.nii');
    WMHnew = fullfile(SubjList{iS}, 'WMH_SEGM.nii');
    
    FLAIRdir1 = fullfile(CaroleDir, Name1, 'FLAIR2T1');
    FLAIRdir2 = fullfile(CaroleDir, Name2, 'FLAIR2T1');
    WMHdir1 = fullfile(CaroleDir, Name1, 'WMH', 'Lesion');
    WMHdir2 = fullfile(CaroleDir, Name2, 'WMH', 'Lesion');
    
    if exist(FLAIRdir1,'dir') || exist(FLAIRdir2, 'dir')
        
        FLAIR1 = xASL_adm_GetFileList(FLAIRdir1, '^FLAIR.*\.nii$', 'FPList', [0 Inf]);
        FLAIR2 = xASL_adm_GetFileList(FLAIRdir2, '^FLAIR.*\.nii$', 'FPList', [0 Inf]);
        WMH1 = xASL_adm_GetFileList(WMHdir1, '^.*Lesion.*\.nii$', 'FPList', [0 Inf]);
        WMH2 = xASL_adm_GetFileList(WMHdir2, '^.*Lesion.*\.nii$', 'FPList', [0 Inf]);

        % check if complete, to proceed copying
        bCopy = false;
        if ~isempty(FLAIR1) && ~isempty(WMH1)
            FLAIRold = FLAIR1{1};
            WMHold = WMH1{1};
            bCopy = true;
        elseif ~isempty(FLAIR2) && ~isempty(WMH2)
            FLAIRold = FLAIR2{1};
            WMHold = WMH2{1};
            bCopy = true;
        end

        % let's copy. we need to overwrite here, as we trust Carole' FLAIR
        % more than ours (we assume that her WMH segmentation is aligned
        % with her FLAIR, not necessarily with ours)
        if bCopy
            % first we remove all pre-existing FLAIRs/WMHs to make sure we wont have doubles
            xASL_adm_DeleteFileList(SubjList{iS}, '^.*(FLAIR|WMH).*$', false, [0 Inf]); % non-recursively
            
            xASL_Copy(FLAIRold, FLAIRnew, true); % need to overwrite here
            xASL_Copy(WMHold, WMHnew, true); % need to overwrite here
        else
            NotCopiedList.(['n' num2str(length(fields(NotCopiedList)))]) = NameNew;
        end
    else
        NotCopiedList.(['n' num2str(length(fields(NotCopiedList)))]) = NameNew;
    end
end

%% SAVE THE MISSING LIST HERE IN JSON FILES

if length(fields(NotCopiedList))>0
    warning('We couldnt find WMH segmentations for all subjects, please check //analysis/NotCopiedWMH_List.json');
    SavePath = fullfile(AnalysisDir, 'NotCopiedWMH_List.json');
    xASL_delete(SavePath);
    xASL_adm_SaveJSON(NotCopiedList, SavePath);
end

end
