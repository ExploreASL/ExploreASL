function [x] = xASL_init_FileSystem(x)
%xASL_init_FileSystem Define path/names used throughout study
%
% FORMAT: [x] = xASL_init_FileSystem(x)
%
% INPUT:
%   DataParPath - path to data parameter file (OPTIONAL, required when ProcessData=true, will then be prompted if omitted)
%   ProcessData - true if running pipeline on data, false if only initializing ExploreASL (e.g. path management etc, REQUIRED, will be prompted if omitted)
%   iWorker     - allows parallelization when called externally. iWorker defines which of the parallel ExploreASL calls we are (OPTIONAL, DEFAULT=1)
%   nWorkers    - allows parallelization when called externally. nWorkers defines how many ExploreASL calls are made in parallel (OPTIONAL, DEFAULT=1)
%
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function initializes the file system used throughout ExploreASL, for processing a single dataset/scan.
% It is repeated for each scan, and runs the following parts:
%
% 1. Create folders
% 2. Subject/session definitions
% 3. Add prefixes & suffixes
% 4. Add Subject-specific prefixes
% 5. Add sidecars
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_init_FileSystem(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% ------------------------------------------------------------------------------------
%% Admin
x.P = struct; % (re-)initiate P (otherwise the looping below gets too much)




%% 1) Create folders
xASL_adm_CreateDir(x.D.T1CheckDir); 
% Put this dir creation in separate scripts
% Create derivative dir as subfolder in each subject folder, where to put output,
% keeping original files
% /dartel -> /Population folder, 
% standardize naming creation

xASL_adm_CreateDir(x.D.TissueVolumeDir);
% xASL_adm_CreateDir(x.D.CoregDir);
% xASL_adm_CreateDir(x.D.FlowFieldCheck);

xASL_adm_CreateDir(x.D.FLAIR_CheckDir);
% xASL_adm_CreateDir(x.D.FLAIR_REGDIR);

xASL_adm_CreateDir(x.D.ASLCheckDir);
xASL_adm_CreateDir(x.D.MotionDir);
% xASL_adm_CreateDir(x.D.RawEPIdir);

xASL_adm_CreateDir(x.D.M0CheckDir);  
xASL_adm_CreateDir(x.D.M0regASLdir);


%% ------------------------------------------------------------------------------------
%% 2) Subject/session definitions
if isfield(x,'SUBJECTDIR')
    [~, x.P.SubjectID] = xASL_fileparts(x.SUBJECTDIR);
    x.iSubject = find(cellfun(@(y) strcmp(y, x.P.SubjectID), x.SUBJECTS));
end
if isfield(x,'SESSIONDIR')
    [~, x.P.SessionID] = xASL_fileparts(x.SESSIONDIR);
    x.iSession = find(cellfun(@(y) strcmp(y, x.P.SessionID), x.SESSIONS));
end
if isfield(x,'iSubject') && isfield(x,'iSession')
    x.iSubjectSession = (x.iSubject-1)*x.nSessions + x.iSession; % It goes Sub1Sess1, Sub1Sess2, Sub2Sess1, Sub2Sess2, Sub3Sess1...
end


%% ------------------------------------------------------------------------------------
%% 3) Add prefixes & suffixes
x.P.STRUCT = 'T1';

% FileTypes in SUBJECTDIR
FileDef{1} = {'FLAIR' 'T1' 'T1_filled' 'c1T1' 'c2T1' 'c3T1' 'j_T1' 'y_T1' 'WMH_SEGM' 'R1' 'PV_pGM' 'PV_pWM' 'PV_WMH_SEGM'};

% FileTypes in SESSIONDIR
FileDef{2} = {'y_ASL' 'ASL4D' 'ASL4D_RevPE'...
	           'CBF' 'qCBF' 'qCBF4D' 'qCBF_untreated' 'qCBF_masked' 'despiked_ASL4D' ...
			   'PseudoCBF' 'PWI' 'PWI4D' 'mean_PWI_Clipped' 'mean_PWI_Clipped_DCT' 'M0' 'mean_control' 'SD' 'SNR' 'SD_control'...
			   'SNR_control' 'SliceGradient' 'SliceGradient_extrapolated' 'FoV' 'TT' 'PVgm' 'PVwm' 'PVcsf' 'PVwmh' 'CBFgm' 'CBFwm' 'MaskSusceptibilityPop' 'TotalGMPop' 'DeepWMPop' 'HammersPop' 'MNIStructuralPop' 'LeftRightPop' 'HOcort_CONNPop' 'HOsub_CONNPop'}; 

Prefix = {'r' 'm' 's' 'mr' 'rmr' 'rr' 'temp_' 'rtemp_' 'mask_' 'BiasField_' 'noSmooth_'}; % r=resample m=modulate s=smooth w=warp q=quantified p=probability % USE t for TEMP? replace w by r
Suffix = {'_backup' '_ORI'};
% use "b" for backup & "o" for original, to reduce the number of suffixes
% need to define various ASL4D_session files still

for iFD=1:length(FileDef)
    for iD=1:length(FileDef{iFD})
        x.P.(FileDef{iFD}{iD}) = FileDef{iFD}{iD};
    end
end
% Add suffixes
for iSuffix=1:length(Suffix)
    x.P.([FileDef{iFD}{iD} Suffix{iSuffix}]) = [FileDef{iFD}{iD} Suffix{iSuffix}];
end     
% Add prefixes
for iPref=1:length(Prefix)
    x.P.([Prefix{iPref} FileDef{iFD}{iD}]) = [Prefix{iPref} FileDef{iFD}{iD}];

    % Add prefixes & suffixes
    for iSuffix=1:length(Suffix)
        x.P.([Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}]) = [Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}];
    end            
end        
        
        
%% ------------------------------------------------------------------------------------
%% 4) Add Subject-specific prefixes
if isfield(x.P,'SubjectID')
    
    for iSess=1:length(x.SESSIONS)
        x.P.SessionDir{iSess} = fullfile(x.D.ROOT,x.P.SubjectID,x.SESSIONS{iSess});
    end

    if ~isfield(x.P,'SessionID')
        x.P.SessionID = 'ASL_1';
    end

    %% ------------------------------------------------------------------------------------
    %% File definitions
    Path{1} = x.SUBJECTDIR;
    Path{2} = fullfile(x.D.ROOT,x.P.SubjectID,x.P.SessionID);

    for iFD=1:length(FileDef)
        if iFD==1
               Pop_suffix = [x.P.SubjectID]; % pop == population analysis, which is in common/standard space
        elseif iFD==2
               Pop_suffix = [x.P.SubjectID '_' x.P.SessionID]; % cave multiple sessions
        end     

        for iD=1:length(FileDef{iFD})
            % Create file & path
            x.P.(['File_' FileDef{iFD}{iD}]) = [FileDef{iFD}{iD} '.nii'];
            x.P.(['Path_' FileDef{iFD}{iD}]) = fullfile(Path{iFD},[FileDef{iFD}{iD} '.nii']);
            x.P.(['Pop_Path_' FileDef{iFD}{iD}]) = fullfile(x.D.PopDir,[FileDef{iFD}{iD} '_' Pop_suffix '.nii']);

            % Add suffixes
            for iSuffix=1:length(Suffix)
                x.P.(['File_' FileDef{iFD}{iD} Suffix{iSuffix}]) = [FileDef{iFD}{iD} Suffix{iSuffix} '.nii'];
                x.P.(['Path_' FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(Path{iFD},[FileDef{iFD}{iD} Suffix{iSuffix} '.nii']);
                x.P.(['Pop_Path_' FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(x.D.PopDir,[FileDef{iFD}{iD} Suffix{iSuffix} '_' Pop_suffix '.nii']);
            end     

            % Add prefixes
            for iPref=1:length(Prefix)
                x.P.(['File_' Prefix{iPref} FileDef{iFD}{iD}]) = [Prefix{iPref} FileDef{iFD}{iD} '.nii'];
                x.P.(['Path_' Prefix{iPref} FileDef{iFD}{iD}]) = fullfile(Path{iFD},[Prefix{iPref} FileDef{iFD}{iD} '.nii']);
                x.P.(['Pop_Path_' Prefix{iPref} FileDef{iFD}{iD}]) = fullfile(x.D.PopDir,[Prefix{iPref} FileDef{iFD}{iD} '_' Pop_suffix '.nii']);

                % Add prefixes & suffixes
                for iSuffix=1:length(Suffix)
                    x.P.(['File_' Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}]) = [Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix} '.nii'];
                    x.P.(['Path_' Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(Path{iFD},[Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix} '.nii']);
                    x.P.(['Pop_Path_' Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(x.D.PopDir,[Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix} '_' Pop_suffix '.nii']);
                end            
            end
        end
    end
    
    %% ------------------------------------------------------------------------------------------
    %% Add custom cases
    if isfield(x, 'SESSIONDIR')
        x.P.Path_MaskVascular = fullfile(x.SESSIONDIR, 'MaskVascular.nii');
    end
    x.P.Pop_Path_MaskVascular = fullfile(x.D.PopDir, ['MaskVascular_' x.P.SubjectID '_' x.P.SessionID '.nii']);
    x.P.Path_Pop_MaskSusceptibility = fullfile(x.D.PopDir, ['rMaskSusceptibility_' x.P.SubjectID '_' x.P.SessionID '.nii']);    
end


%% ------------------------------------------------------------------------------------------
%% 5) Add sidecars
FieldList = fieldnames(x.P);
for iL=12:length(FieldList) % here we should remove the prefixes
    
    % Add .mat sidecars
    CarSuf1 = {'_parms_mat' '_sn_mat' '_mat' '_json'}; % suffixes for sidecars (dots replaced by underscores)
    CarSuf2 = {'_parms.mat' '_sn.mat' '.mat' '.json'}; % suffixes for sidecars
    for iS=1:length(CarSuf1)
        x.P.([FieldList{iL} CarSuf1{iS}]) = [x.P.(FieldList{iL})(1:end-4) CarSuf2{iS}];
    end
end    


end