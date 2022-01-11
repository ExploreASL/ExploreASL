function [x] = xASL_init_FileSystem(x)
%xASL_init_FileSystem Define path/names used throughout study
%
% FORMAT: [x] = xASL_init_FileSystem(x)
%
% INPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
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
% 6. Add atlas paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_init_FileSystem(x);
% __________________________________
% Copyright 2015-2021 ExploreASL

%% ------------------------------------------------------------------------------------
%% Admin
x.P = struct; % (re-)initiate P (otherwise the looping below gets too much)


%% 1) Create folders
[fixedPopDir,~] = xASL_adm_GetPopulationDirs(x);

% Create fixed directories
for iDir = 1:numel(fixedPopDir)
    xASL_adm_CreateDir(fixedPopDir{iDir});
end

%% ------------------------------------------------------------------------------------
%% 2) Subject/session definitions
if isfield(x.dir,'SUBJECTDIR')
    [~, x.P.SubjectID] = xASL_fileparts(x.dir.SUBJECTDIR);
    x.iSubject = find(cellfun(@(y) strcmp(y, x.P.SubjectID), x.SUBJECTS));
end
if isfield(x.dir,'SESSIONDIR')
    [~, x.P.SessionID] = xASL_fileparts(x.dir.SESSIONDIR);
    x.iSession = find(cellfun(@(y) strcmp(y, x.P.SessionID), x.SESSIONS));
end
if isfield(x,'iSubject') && isfield(x,'iSession')
    x.iSubjectSession = (x.iSubject-1)*x.dataset.nSessions + x.iSession; % It goes Sub1Sess1, Sub1Sess2, Sub2Sess1, Sub2Sess2, Sub3Sess1...
end


%% ------------------------------------------------------------------------------------
%% 3) Add prefixes & suffixes
x.P.STRUCT = 'T1';

% FileTypes in SUBJECTDIR
FileDef{1} = {'FLAIR' 'T1' 'T1c' 'T2' 'T1_filled' 'c1T1' 'c2T1' 'c3T1' 'j_T1' 'y_T1' 'WMH_SEGM' 'R1' 'PV_pGM' 'PV_pWM' 'PV_WMH_SEGM'};

% FileTypes in SESSIONDIR
FileDef{2} = {'y_ASL' 'ASL4D' 'ASL4D_RevPE' 'M0_RevPE'...
	           'CBF' 'qCBF' 'qCBF4D' 'qCBF_untreated' 'qCBF_masked' 'despiked_ASL4D' ...
			   'PseudoCBF' 'PWI' 'PWI4D' 'mean_PWI_Clipped' 'mean_PWI_Clipped_DCT' 'M0' 'mean_control' 'SD' 'SNR' 'SD_control'...
			   'SNR_control' 'SliceGradient' 'SliceGradient_extrapolated' 'FoV' 'TT' 'PVgm' 'PVwm' 'PVcsf' 'PVwmh' 'CBFgm' 'CBFwm' 'MaskSusceptibilityPop' 'TotalGMPop' 'TotalWMPop' 'DeepWMPop' 'HammersPop' 'MNIStructuralPop' 'LeftRightPop' 'HOcort_CONNPop' 'HOsub_CONNPop'}; 

% Use "b" for backup & "o" for original, to reduce the number of suffixes (need to define various ASL4D_session files still)
% r=resample m=modulate s=smooth w=warp q=quantified p=probability % USE t for TEMP? replace w by r
Prefix = {'r' 'm' 's' 'mr' 'rmr' 'rr' 'temp_' 'rtemp_' 'mask_' 'BiasField_' 'noSmooth_'};
Suffix = {'_backup' '_ORI'};

% Define population directories and corresponding subfolder names
PopTable = {'qCBF', 'perf';...
            'PWI', 'perf';...
            'PWI_part1', 'perf';...
            'PWI_part2', 'perf';...
            'M0', 'perf';...
            'noSmooth_M0', 'perf';...
            'mean_control', 'perf';...
            'SD', 'perf';...
            'SNR', 'perf';...
            'rT1_ORI', 'anat';...
            'rT1', 'T1w';...
            'mrc1T1', 'T1w';...
            'mrc2T1', 'T1w';...
            'mrc3T1', 'T1w';...
            'rc1T1', 'T1w';...
            'rc2T1', 'T1w';...
            'rc3T1', 'T1w';...
            'PV_pGM', 'T1w';...
            'PV_pWM', 'T1w';...
            'rFLAIR', 'FLAIR';...
            'WMH_SEGM', 'FLAIR';...
            'rWMH_SEGM', 'FLAIR';...
            'mrWMH_SEGM', 'FLAIR';...
            'PV_WMH_SEGM', 'anat';...
            'rT2', 'FLAIR';...
            'SliceGradient', 'SliceGradient';...
            'SliceGradient_extrapolated', 'SliceGradient'...
            };

% Iterate over file definitions
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
    Path{1} = x.dir.SUBJECTDIR;
    Path{2} = fullfile(x.D.ROOT,x.P.SubjectID,x.P.SessionID);
    x = xASL_init_FileSystemDefinitions(x,FileDef,Path,Prefix,Suffix,PopTable);
    
    %% ------------------------------------------------------------------------------------------
    %% Add custom cases
    if isfield(x.dir, 'SESSIONDIR')
        x.P.Path_MaskVascular = fullfile(x.dir.SESSIONDIR, 'MaskVascular.nii');
    end
    x.P.Pop_Path_MaskVascular = fullfile(x.D.PopDir, ['MaskVascular_' x.P.SubjectID '_' x.P.SessionID '.nii']);
    x.P.Path_Pop_MaskSusceptibility = fullfile(x.D.PopDir, ['rMaskSusceptibility_' x.P.SubjectID '_' x.P.SessionID '.nii']);    
	
	x.P.Path_ASL4Dcontext = fullfile(Path{2},'ASL4Dcontext.tsv');
	x.P.Path_ASL4Dcontext_Source = fullfile(Path{2},'ASL4Dcontext_Source.tsv');
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


%% ------------------------------------------------------------------------------------------
%% 6) Add atlas paths
x.P.Atlas.TotalGM                           = fullfile(x.D.MapsSPMmodifiedDir, 'TotalGM.nii');
x.P.Atlas.TotalWM                           = fullfile(x.D.MapsSPMmodifiedDir, 'TotalWM.nii');
x.P.Atlas.DeepWM                            = fullfile(x.D.MapsSPMmodifiedDir, 'DeepWM.nii');
x.P.Atlas.WholeBrain                        = fullfile(x.D.MapsSPMmodifiedDir, 'WholeBrain.nii');
x.P.Atlas.Tatu_ACA_MCA_PCA                  = fullfile(x.D.MapsSPMmodifiedDir, 'VascularTerritories', 'CortVascTerritoriesTatu.nii.gz');
x.P.Atlas.Tatu_ICA_PCA                      = fullfile(x.D.MapsSPMmodifiedDir, 'VascularTerritories', 'TatuICA_PCA.nii.gz');
x.P.Atlas.Tatu_ICA_L_ICA_R_PCA              = fullfile(x.D.MapsSPMmodifiedDir, 'VascularTerritories', 'LabelingTerritories.nii.gz');
x.P.Atlas.Tatu_ACA_MCA_PCA_Prox_Med_Dist    = fullfile(x.D.MapsSPMmodifiedDir, 'VascularTerritories', 'ATTbasedFlowTerritories.nii.gz');

% Add atlases from atlasDir
x = xASL_init_AtlasList(x,x.D.AtlasDir);

end

% Add all NIFTIs from the atlasDir to the x.P.Atlas field
function [x] = xASL_init_AtlasList(x,atlasDir)

% Get all files in atlas directory
filesInAtlasDir = xASL_adm_GetFileList(atlasDir,'^.+\.nii$');

% Iterate over atlases
for iFile=1:numel(filesInAtlasDir)
    % Get current atlas
    [~,currentAtlas] = xASL_fileparts(filesInAtlasDir{iFile});
    x.P.Atlas.(currentAtlas) = filesInAtlasDir{iFile};
end

end


%% xASL_init_GetSpecificPopDir
function specificPopDir = xASL_init_GetSpecificPopDir(x,currentName,PopulationTokens,PopulationDirs)

    specificPopDir = fullfile(x.D.PopDir);
    popDirIndex = find(ismember(PopulationTokens, currentName), 1);
    if ~isempty(popDirIndex) && isfield(x.D,PopulationDirs{popDirIndex})
        specificPopDir = x.D.(PopulationDirs{popDirIndex});
    end

end


%% xASL_adm_CreatePopulationDirs
function [fixedPopDir,conditionalPopDir] = xASL_adm_GetPopulationDirs(x)

    % Always create these population directories
    fixedPopDir =         { x.D.PopDir,...
                            x.D.perf,...
                            x.D.anat,...
                            x.D.stats,...
                            x.D.T1w,...
                            x.D.FLAIR,...
                            x.D.ASL,...
                            x.S.StatsDir,...
                            x.D.StatsMaps,...
                            x.D.HistogramDir,...
                            x.D.ASLCheckDir,...
                            x.D.M0CheckDir,...
                            x.D.TTCheckDir,...
                            x.D.FLAIR_CheckDir, ...
                            x.D.MotionDir,...
                            x.D.M0regASLdir, ...
                            x.D.TissueVolumeDir, ...
                            x.D.SliceGradient ...
                            };
    
    % Conditionally create these population directories
    conditionalPopDir =   { x.D.T1CheckDir, ...
                            x.D.CoregDir, ...
                            x.D.T1c_CheckDir, ...
                            x.D.T2_CheckDir, ...
                            x.D.FLAIR_REGDIR, ...
                            x.D.FlowFieldCheck, ...
                            x.D.LongRegCheckDir, ...
                            x.D.LesionCheckDir, ...
                            x.D.ROICheckDir, ...
                            x.D.ExclusionDir, ...
                            x.D.DICOMparameterDir, ...
                            x.D.SNRdir, ...
                            x.D.SliceCheckDir, ...
                            x.D.RawDir, ...
                            x.D.RawEPIdir, ...
                            x.D.T1_ASLREGDIR, ...
                            x.D.TemplatesStudyDir, ...
                            x.D.SpaghettiDir};

end


%% xASL_init_FileSystemDefinitions
function x = xASL_init_FileSystemDefinitions(x,FileDef,Path,Prefix,Suffix,PopTable)

    PopulationTokens = PopTable(:,1)';
    PopulationDirs = PopTable(:,2)';

    for iFD=1:length(FileDef)
        if iFD==1
            % pop == population analysis, which is in common/standard space
            Pop_suffix = [x.P.SubjectID];
        elseif iFD==2
            % cave multiple sessions
            Pop_suffix = [x.P.SubjectID '_' x.P.SessionID];
        end

        for iD=1:length(FileDef{iFD})
            % Create file & path
            x.P.(['File_' FileDef{iFD}{iD}]) = [FileDef{iFD}{iD} '.nii'];
            x.P.(['Path_' FileDef{iFD}{iD}]) = fullfile(Path{iFD},[FileDef{iFD}{iD} '.nii']);

            % Put files in specific population directory sub-directories (default to main directory)
            specificPopDir = xASL_init_GetSpecificPopDir(x,FileDef{iFD}{iD},PopulationTokens,PopulationDirs);
            x.P.(['Pop_Path_' FileDef{iFD}{iD}]) = fullfile(specificPopDir,[FileDef{iFD}{iD} '_' Pop_suffix '.nii']);

            % Add suffixes
            for iSuffix=1:length(Suffix)
                x.P.(['File_' FileDef{iFD}{iD} Suffix{iSuffix}]) = [FileDef{iFD}{iD} Suffix{iSuffix} '.nii'];
                x.P.(['Path_' FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(Path{iFD},[FileDef{iFD}{iD} Suffix{iSuffix} '.nii']);

                % Put files in specific population directory sub-directories (default to main directory)
                specificPopDir = xASL_init_GetSpecificPopDir(x,[FileDef{iFD}{iD} Suffix{iSuffix}],PopulationTokens,PopulationDirs);
                x.P.(['Pop_Path_' FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(specificPopDir,[FileDef{iFD}{iD} Suffix{iSuffix} '_' Pop_suffix '.nii']);
            end

            % Add prefixes
            for iPref=1:length(Prefix)
                x.P.(['File_' Prefix{iPref} FileDef{iFD}{iD}]) = [Prefix{iPref} FileDef{iFD}{iD} '.nii'];
                x.P.(['Path_' Prefix{iPref} FileDef{iFD}{iD}]) = fullfile(Path{iFD},[Prefix{iPref} FileDef{iFD}{iD} '.nii']);

                % Put files in specific population directory sub-directories (default to main directory)
                specificPopDir = xASL_init_GetSpecificPopDir(x,[Prefix{iPref} FileDef{iFD}{iD}],PopulationTokens,PopulationDirs);
                x.P.(['Pop_Path_' Prefix{iPref} FileDef{iFD}{iD}]) = fullfile(specificPopDir,[Prefix{iPref} FileDef{iFD}{iD} '_' Pop_suffix '.nii']);

                % Add prefixes & suffixes
                for iSuffix=1:length(Suffix)
                    x.P.(['File_' Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}]) = [Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix} '.nii'];
                    x.P.(['Path_' Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(Path{iFD},[Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix} '.nii']);

                    % Put files in specific population directory sub-directories (default to main directory)
                    specificPopDir = xASL_init_GetSpecificPopDir(x,[Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}],PopulationTokens,PopulationDirs);
                    x.P.(['Pop_Path_' Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix}]) = fullfile(specificPopDir,[Prefix{iPref} FileDef{iFD}{iD} Suffix{iSuffix} '_' Pop_suffix '.nii']);
                end
            end
        end
    end

end

