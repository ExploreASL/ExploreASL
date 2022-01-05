function xASL_adm_UpdateStructure(datasetRoot)
%ExploreASL_UpdateStructure Updates the structure of the temp & lock folder from old formats to the current format
%
% FORMAT: ExploreASL_UpdateStructure(datasetRoot)
%
% INPUT: (either datasetRoot or x is REQUIRED, the other OPTIONAL)
%   datasetRoot   - root directory of your BIDS dataset
%
% OUTPUT:
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function allows backward compatibility between different ExploreASL WIP versions
%              However, run with caution. If after running this function
%              you still get unexpected errors, try rerunning the full
%              pipeline (or only for a single subject), to see if any other
%              changes were not backwards compatible
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: ExploreASL_UpdateStructure([], x);
% __________________________________
% Copyright 2015-2021 ExploreASL



%% Check if we need to initialize ExploreASL
if ~isempty(datasetRoot) && exist(datasetRoot,'dir')
    x = ExploreASL_Initialize(datasetRoot);
else
    error(['Dataset root ' datasetRoot ' is invalid/doesnt exist']);
end

%% Rename the directory dartel to Population
if exist(fullfile(x.D.ROOT, 'dartel'),'dir')
	xASL_delete(fullfile(x.D.ROOT,'Population'));
	xASL_Move(fullfile(x.D.ROOT, 'dartel'),fullfile(x.D.ROOT, 'Population'))
end

%% Rename the structural lock directories
if exist(fullfile(x.D.ROOT, 'lock', 'T1'),'dir')
	xASL_Move(fullfile(x.D.ROOT, 'lock', 'T1'),fullfile(x.D.ROOT, 'lock', 'xASL_module_Structural'));
end
pathStruct = fullfile(x.D.ROOT, 'lock', 'xASL_module_Structural');

% List all the patients
fList = xASL_adm_GetFileList(pathStruct,'^.+$',[],[],1);
for iL=1:length(fList)
    xASL_TrackProgress(iL,length(fList));
    % Rename the module directory
    OldDir = fullfile(fList{iL},'T1_module');
    pathLoc = fullfile(fList{iL},'xASL_module_Structural');
    if exist(OldDir,'dir')
        xASL_Move(OldDir, pathLoc);
    end

    % Rename the lock-files
    OldNames{1} = {'001_coreg_T12MNI'      '002_segment_T1'         '005_reslice2DARTEL'         '003_tissue_volume'         '006_visualize'};
    OldNames{2} = {'010_LinearReg_T1w-MNI' '070_CAT12_Segment_T1w'  '090_Resample-StandardSpace' '110_Get_TissueVolume'      '120_Visual_QC'...
        '020_LinearReg_FLAIR-T1w' '030_FLAIR_BiasfieldCorrect'      '050_LST_Segment_FLAIR_WMH'  '060_LST_T1w_LesionFilling_WMH' '080_CleanUp_WMH_SEGM'};
    NewNames    = {'010_LinearReg_T1w2MNI' '060_Segment_T1w'        '080_Resample2StandardSpace' '090_GetVolumetrics'        '100_VisualQC'...
        '020_LinearReg_FLAIR2T1w' '030_FLAIR_BiasfieldCorrection'   '040_LST_Segment_FLAIR_WMH'  '050_LST_T1w_LesionFilling_WMH' '070_CleanUpWMH_SEGM'};

	DeleteNames = {'100_visualize' '100_VisualQC' '999_ready'};

	for iGen=1:length(OldNames)
		if length(OldNames{iGen})>length(NewNames)
			error('Something wrong with the lock-file updating');
		end
		for iLock=1:length(OldNames{iGen})
			OldPath = fullfile(pathLoc, [OldNames{iGen}{iLock} '.status']);
			NewPath = fullfile(pathLoc, [NewNames{iLock} '.status']);
			if exist(OldPath, 'file')
				xASL_Move(OldPath, NewPath, true);
			end
		end
	end

	for iLock = 1:length(DeleteNames)
		NewPath = fullfile(pathLoc, [DeleteNames{iLock} '.status']);
		xASL_delete(NewPath);
	end
end


%% Rename the ASL lock directories
if exist(fullfile(x.D.ROOT, 'lock', 'ASL'),'dir')
	xASL_Move(fullfile(x.D.ROOT, 'lock', 'ASL'),fullfile(x.D.ROOT, 'lock', 'xASL_module_ASL'));
end
pathASL = fullfile(x.D.ROOT, 'lock', 'xASL_module_ASL');

% List all the patients
fList = xASL_adm_GetFileList(pathASL,'^.+$',[],[],1);
for iL=1:length(fList)
    xASL_TrackProgress(iL,length(fList));
    % List all the session
    fSess = xASL_adm_GetFileList(fList{iL},'^.+$',false,[],1);
    for iS = 1:length(fSess)
        % Rename the session directory
        sessionName = ['ASL_' num2str(iS)];
        DirOld = fullfile(fList{iL},['ASL_module_' sessionName]);
        pathLoc = fullfile(fList{iL},['xASL_module_ASL_' sessionName]);
        if exist(DirOld,'dir')
            xASL_Move(DirOld, pathLoc);
        end

        % Rename the lock-files (OldNames1 & 2 are different generations)
        OldNames{1} = {'002_realign_ASL' '0025_register_ASL' '003_reslice_ASL' '0035_PreparePV' '004_realign_reslice_M0' '005_quantification' '006_visualize'};
        OldNames{2} = {'020_realign_ASL' '025_register_ASL'  '030_reslice_ASL' '035_PreparePV'  '040_realign_reslice_M0' '050_quantification' '060_VisualQC'};
        NewNames    = {'020_MoCoASL'     '030_RegisterASL'   '040_ResliceASL'  '050_PreparePV'  '060_ProcessM0'          '070_Quantification' '090_VisualQC'};

		DeleteNames = {'090_VisualQC'};

		for iGen=1:length(OldNames) % different generations
			if length(OldNames{iGen})~=length(NewNames)
				error('Something wrong with the lock-file updating');
			end
			for iLock=1:length(OldNames{iGen})
				OldPath = fullfile(pathLoc,[OldNames{iGen}{iLock} '.status']);
				NewPath = fullfile(pathLoc,[NewNames{iLock} '.status']);
				if xASL_exist(OldPath, 'file')
					xASL_Move(OldPath, NewPath, true);
				end
			end
		end

		for iLock = 1:length(DeleteNames)
			NewPath = fullfile(pathLoc, [DeleteNames{iLock} '.status']);
			xASL_delete(NewPath);
		end
    end
end

%% Rename the population lock directories
if exist(fullfile(x.D.ROOT, 'lock', 'QC'),'dir')
	xASL_Move(fullfile(x.D.ROOT, 'lock', 'QC'),fullfile(x.D.ROOT, 'lock', 'xASL_module_Population'));
	xASL_Move(fullfile(x.D.ROOT, 'lock', 'xASL_module_Population','QC_module'),fullfile(x.D.ROOT, 'lock', 'xASL_module_Population', 'xASL_module_Population'));

	pathLoc = fullfile(x.D.ROOT, 'lock', 'xASL_module_Population', 'xASL_module_Population');

	if exist(fullfile(pathLoc,'999_ready.status'),'file'), xASL_delete(fullfile(pathLoc,'999_ready.status'));end
	if exist(fullfile(pathLoc,'001_QA_checkImreg.status'),'file'), xASL_delete(fullfile(pathLoc,'001_QA_checkImreg.status'));end

	if exist(fullfile(x.D.ROOT,'lock','NORMALIZE','NORMALIZE_module','002_mean_warp_templates.status'),'file')
		xASL_Move(fullfile(x.D.ROOT,'lock','NORMALIZE','NORMALIZE_module','002_mean_warp_templates.status'),fullfile(pathLoc,'010_mean_warp_templates.status'))

		xASL_delete(fullfile(x.D.ROOT,'lock','NORMALIZE','NORMALIZE_module','999_ready.status'));
		xASL_delete(fullfile(x.D.ROOT,'lock','NORMALIZE'));
	end

    % Rename the lock-files
    OldNames = {'002_QA_checkDICOMvalues' '003_QA_volume_stats' '004_QA_motion_stats' '005_ROI_analysis'};
    NewNames = {'050_QA_checkDICOMvalues' '060_QA_volume_stats' '070_QA_motion_stats' '080_ROI_analysis'};

    if length(OldNames)~=length(NewNames)
        error('Something wrong with the lock-file updating');
    end
    for iLock=1:length(OldNames)
        OldPath = fullfile(pathLoc,OldNames{iLock});
        NewPath = fullfile(pathLoc,NewNames{iLock});
        if exist(OldPath, 'file')
            xASL_Move(OldPath, NewPath, true);
        end
    end

end

%% Delete the population lock files
pathLoc = fullfile(x.D.ROOT, 'lock', 'xASL_module_Population', 'xASL_module_Population');
xASL_delete(fullfile(pathLoc,'010_mean_warp_templates.status'));
xASL_delete(fullfile(pathLoc,'050_QA_checkDICOMvalues.status'));
xASL_delete(fullfile(pathLoc,'060_QA_volume_stats.status'));
xASL_delete(fullfile(pathLoc,'070_QA_motion_stats.status'));
xASL_delete(fullfile(pathLoc,'080_ROI_analysis.status'));
xASL_delete(fullfile(pathLoc,'999_ready.status'));

%% Delete the .mat files
xASL_delete(fullfile(x.D.ROOT,'x.mat'));

% And for each subject
fList = xASL_adm_GetFileList(x.D.ROOT,'^.+$',[],[],1);
for iL=1:length(fList)
	xASL_delete(fullfile(fList{iL},'x.mat'));
end
end
