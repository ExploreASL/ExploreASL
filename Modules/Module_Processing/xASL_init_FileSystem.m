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
% So note that ExploreASL-wide parameters, that are not specific to a single scan, are not defined here.
% E.g., the atlas paths are defined in xASL_init_MapsAndAtlases.m
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
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% ------------------------------------------------------------------------------------
%% Admin
x.P = struct; % (re-)initiate P (otherwise the looping below gets too much)


%% 1) Create folders
% Put this dir creation in separate scripts
% Create derivative dir as subfolder in each subject folder, where to put output,
% keeping original files
% /dartel -> /Population folder, standardize naming creation

xASL_adm_CreateDir(x.D.PopDir);


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
FileDef{1} = {'FLAIR' 'T1' 'T1c' 'T2' 'T1_filled' 'c1T1' 'c2T1' 'c3T1' 'j_T1' 'y_T1' 'WMH_SEGM' 'R1' 'PV_pGM' 'PV_pWM' 'PV_pCSF' 'PV_WMH_SEGM'};

% FileTypes in SESSIONDIR
FileDef{2} = {'y_ASL' 'ASL4D' 'ASL4D_RevPE' 'M0_RevPE'...
	           'CBF' 'qCBF' 'CBF4D' 'qCBF4D' 'qCBF_untreated' 'qCBF_masked' 'despiked_ASL4D' ...
			   'PseudoCBF' 'PWI' 'PWI3D' 'PWI4D' 'mean_PWI_Clipped' 'mean_PWI_Clipped_DCT' 'M0' 'mean_control' 'SD' 'SNR' 'SD_control'...
			   'SNR_control' 'SliceGradient' 'SliceGradient_extrapolated' 'FoV' 'TT' 'ATT' 'Tex' 'ABV' 'ITT' 'PVgm' 'PVwm' 'PVcsf' 'PVt1' 'PVwmh' 'CBFgm' 'CBFwm'}; 

Prefix = {'r' 'm' 'mr' 'rmr' 'rr' 'temp_' 'rtemp_' 'mask_' 'BiasField_' 'noSmooth_'}; % r=resample m=modulate s=smooth w=warp q=quantified p=probability % USE t for TEMP? replace w by r

ListFieldsSuffixBackup{1} = {};
ListFieldsSuffixBackup{2} = {'M0'};
ListFieldsSuffixORI{1} = {'T1' 'T2' 'T1c' 'c1T1' 'c2T1' 'FLAIR' 'WMH_SEGM'};
ListFieldsSuffixORI{2} = {'mean_PWI_Clipped_ORI' 'ASL4D'};
ListFieldsSuffixORI_r_Pop{1} = {'T1' 'c1T1' 'c2T1' 'rc1T1' 'rc2T1'};
ListFieldsSuffixORI_r_Pop{2} = {};

% List of prefixes 
ListFieldsPrefixBiasField{1} = {'FALIR','T1'};ListFieldsPrefixBiasField{2} = {};
ListFieldsPrefixMask{1} = {'T1'};ListFieldsPrefixMask{2} = {'M0'};
ListFieldsPrefixRtemp{1} = {};ListFieldsPrefixRtemp{2} = {'despiked_ASL4D'};
ListFieldsPrefixTemp{1} = {};ListFieldsPrefixTemp{2} = {'despiked_ASL4D'};
ListFieldsPrefixRR{1} = {};ListFieldsPrefixRR{2} = {'M0'};
ListFieldsPrefixMR{1} = {'FLAIR'};ListFieldsPrefixMR{2} = {};
ListFieldsPrefixRMR{1} = {'FLAIR'};ListFieldsPrefixRMR{2} = {};
ListFieldsPrefixM{1} = {'FLAIR' 'T1'};ListFieldsPrefixM{2} = {};
ListFieldsPrefixR{1} = {'c1T1' 'c2T1' 'c3T1' 'WMH_SEGM' 'T1'};ListFieldsPrefixR{2} = {'despiked_ASL4D' 'M0' 'PWI' 'mean_PWI_Clipped'};

ListFieldsPrefixNosmoothPop{1} = {};ListFieldsPrefixNosmoothPop{2} = {'M0'};
ListFieldsPrefixMaskPop{1} = {};ListFieldsPrefixMaskPop{2} = {'M0'};
ListFieldsPrefixMRPop{1} = {'c1T1' 'c2T1' 'c3T1' 'WMH_SEGM'};ListFieldsPrefixMRPop{2} = {};
ListFieldsPrefixRPop{1} = {'c1T1' 'c2T1' 'c3T1' 'WMH_SEGM' 'T1' 'T2' 'T1c'};ListFieldsPrefixRPop{2} = {};

% use "b" for backup & "o" for original, to reduce the number of suffixes
% need to define various ASL4D_session files still

for iFD=1:length(FileDef)
    for iD=1:length(FileDef{iFD})
        x.P.(FileDef{iFD}{iD}) = FileDef{iFD}{iD};
    end
end

%% Add prefixes
for iPref=1:length(Prefix)
    x.P.([Prefix{iPref} FileDef{iFD}{iD}]) = [Prefix{iPref} FileDef{iFD}{iD}];
end        
        
        
%% ------------------------------------------------------------------------------------
%% 4) Add Subject-specific prefixes
if isfield(x.P,'SubjectID')
    
    for iSess=1:length(x.SESSIONS)
        x.P.SessionDir{iSess} = fullfile(x.dir.xASLDerivatives,x.P.SubjectID,x.SESSIONS{iSess});
    end

    if ~isfield(x.P,'SessionID')
        % No sessions found, defaulting to a single ASL_1 session
        % This is the case for non-ASL modules
        x.P.SessionID = 'ASL_1';
    end

    %% ------------------------------------------------------------------------------------
    %% File definitions. Make it ready for use without further using the time-consuming function fullfile
    Path{1} = fullfile(x.dir.SUBJECTDIR,' ');
	Path{1} = Path{1}(1:(end-1));
    Path{2} = fullfile(x.dir.xASLDerivatives,x.P.SubjectID,x.P.SessionID,' ');
	Path{2} = Path{2}(1:(end-1));
	PathPop = fullfile(x.D.PopDir, ' ');
	PathPop = PathPop(1:(end-1));

	% Custom 'File' prefix
	x.P.File_ASL4D = 'File_ASL4D.nii';
	x.P.File_despiked_ASL4D = 'File_despiked_ASL4D.nii';

    for iFD=1:length(FileDef)
        if iFD==1
               Pop_suffix = [x.P.SubjectID]; % pop == population analysis, which is in common/standard space
        elseif iFD==2
               Pop_suffix = [x.P.SubjectID '_' x.P.SessionID]; % cave multiple sessions
        end     

		% Add suffix BACKUP
		for iSuffix = 1:length(ListFieldsSuffixBackup{iFD})
			x.P.(['Path_' ListFieldsSuffixBackup{iFD}{iSuffix} '_backup']) = [Path{iFD} ListFieldsSuffixBackup{iFD}{iSuffix} '_backup.nii'];
		end

		% Add suffix ORI
		for iSuffix = 1:length(ListFieldsSuffixORI{iFD})
			x.P.(['Path_' ListFieldsSuffixORI{iFD}{iSuffix} '_ORI']) = [Path{iFD} ListFieldsSuffixORI{iFD}{iSuffix} '_ORI.nii'];
		end
		% Add suffix ORI with prefix R and in Population folder
		for iSuffix = 1:length(ListFieldsSuffixORI_r_Pop{iFD})
			x.P.(['Pop_Path_r' ListFieldsSuffixORI_r_Pop{iFD}{iSuffix} '_ORI']) = [PathPop 'r' ListFieldsSuffixORI_r_Pop{iFD}{iSuffix} '_ORI.nii'];
		end

		% Add prefixes BiasField
		for iPref=1:length(ListFieldsPrefixBiasField{iFD})
			x.P.(['Path_BiasField_' ListFieldsPrefixBiasField{iFD}{iPref}]) = [Path{iFD} 'BiasField_' ListFieldsPrefixBiasField{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixMask{iFD})
			x.P.(['Path_mask_' ListFieldsPrefixMask{iFD}{iPref}]) = [Path{iFD} 'mask_' ListFieldsPrefixMask{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixRtemp{iFD})
			x.P.(['Path_rtemp_' ListFieldsPrefixRtemp{iFD}{iPref}]) = [Path{iFD} 'rtemp_' ListFieldsPrefixRtemp{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixTemp{iFD})
			x.P.(['Path_temp_' ListFieldsPrefixTemp{iFD}{iPref}]) = [Path{iFD} 'temp_' ListFieldsPrefixTemp{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixRR{iFD})
			x.P.(['Path_rr' ListFieldsPrefixRR{iFD}{iPref}]) = [Path{iFD} 'rr' ListFieldsPrefixRR{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixMR{iFD})
			x.P.(['Path_mr' ListFieldsPrefixMR{iFD}{iPref}]) = [Path{iFD} 'mr' ListFieldsPrefixMR{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixRMR{iFD})
			x.P.(['Path_rmr' ListFieldsPrefixRMR{iFD}{iPref}]) = [Path{iFD} 'rmr' ListFieldsPrefixRMR{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixM{iFD})
			x.P.(['Path_m' ListFieldsPrefixM{iFD}{iPref}]) = [Path{iFD} 'm' ListFieldsPrefixM{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixR{iFD})
			x.P.(['Path_r' ListFieldsPrefixR{iFD}{iPref}]) = [Path{iFD} 'r' ListFieldsPrefixR{iFD}{iPref} '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixNosmoothPop{iFD})
			x.P.(['Pop_Path_noSmooth_' ListFieldsPrefixNosmoothPop{iFD}{iPref}]) = [PathPop 'noSmooth_' ListFieldsPrefixNosmoothPop{iFD}{iPref} '_' Pop_suffix '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixMaskPop{iFD})
			x.P.(['Pop_Path_mask_' ListFieldsPrefixMaskPop{iFD}{iPref}]) = [PathPop 'mask_' ListFieldsPrefixMaskPop{iFD}{iPref} '_' Pop_suffix '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixMRPop{iFD})
			x.P.(['Pop_Path_mr' ListFieldsPrefixMRPop{iFD}{iPref}]) = [PathPop 'mr' ListFieldsPrefixMRPop{iFD}{iPref} '_' Pop_suffix '.nii'];
		end

		for iPref=1:length(ListFieldsPrefixRPop{iFD})
			x.P.(['Pop_Path_r' ListFieldsPrefixRPop{iFD}{iPref}]) = [PathPop 'r' ListFieldsPrefixRPop{iFD}{iPref} '_' Pop_suffix '.nii'];
		end

        for iD=1:length(FileDef{iFD})
            % Create file & path
            x.P.(['Path_' FileDef{iFD}{iD}]) = [Path{iFD} FileDef{iFD}{iD} '.nii'];
            x.P.(['Pop_Path_' FileDef{iFD}{iD}]) = [PathPop FileDef{iFD}{iD} '_' Pop_suffix '.nii'];
        end
    end
    
    %% ------------------------------------------------------------------------------------------
    %% Add custom cases
    % Here we add ASL-specific files, that we only need when processing ASL
    % data
    if isfield(x.dir, 'SESSIONDIR')
        x.P.Path_MaskVascular = fullfile(x.dir.SESSIONDIR, 'MaskVascular.nii');
        x.P.Path_BrainMaskProcessing = fullfile(x.dir.SESSIONDIR, 'BrainMaskProcessing.nii');
    end
    x.P.Pop_Path_MaskVascular = [PathPop 'MaskVascular_' x.P.SubjectID '_' x.P.SessionID '.nii'];
    x.P.Pop_Path_BrainMaskProcessing = [PathPop 'BrainMaskProcessing_' x.P.SubjectID '_' x.P.SessionID '.nii'];
    x.P.Pop_Path_MaskSusceptibility = [PathPop 'rMaskSusceptibility_' x.P.SubjectID '_' x.P.SessionID '.nii'];    
	
	x.P.Path_ASL4Dcontext = [Path{2} 'ASL4Dcontext.tsv'];
	x.P.Path_ASL4Dcontext_Source = [Path{2} 'ASL4Dcontext_Source.tsv'];

	%% ------------------------------------------------------------------------------------------
	%% 5) Add sidecars
	% Add custom sidecars
	ListFields = {'Path_ASL4D' 'Path_M0' 'Path_M0_backup'};
	for iL = 1:length(ListFields)
		x.P.([ListFields{iL} '_parms_mat']) = [x.P.(ListFields{iL})(1:end-4) '_parms.mat'];
	end
	ListFields = {'Path_mean_PWI_Clipped'};
	for iL = 1:length(ListFields)
		x.P.([ListFields{iL} '_sn_mat']) = [x.P.(ListFields{iL})(1:end-4) '_sn.mat'];
	end

	ListFields = {'Path_ASL4D' 'Path_despiked_ASL4D' 'File_ASL4D' 'File_despiked_ASL4D'};
	for iL = 1:length(ListFields)
		x.P.([ListFields{iL} '_mat']) = [x.P.(ListFields{iL})(1:end-4) '.mat'];
	end

	ListFields = {'Path_ASL4D' 'Path_despiked_ASL4D' 'Path_rdespiked_ASL4D' 'Path_rtemp_despiked_ASL4D' 'Path_mean_PWI_Clipped' 'Path_mean_control' 'Path_M0' 'Path_M0_backup' 'Pop_Path_PWI4D'};
	for iL = 1:length(ListFields)
		x.P.([ListFields{iL} '_mat']) = [x.P.(ListFields{iL})(1:end-4) '.mat'];
	end
end

end
