function [result, x] = xASL_module_ASL(x)
%xASL_module_ASL ExploreASL module for ASL processing
%
% FORMAT: [result, x] = xASL_module_ASL(x)
%
% INPUT:
%   x  - x structure containing all input parameters (REQUIRED)
%   x.dir.SUBJECTDIR  -  anatomical directory, containing the derivatives of anatomical images (REQUIRED)
%   x.dir.SESSIONDIR  -  ASL directory, containing the derivatives of perfusion images (REQUIRED)
%
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters (REQUIRED)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL module processes the ASL
% images, i.e. ASL4D, M0, etc (if present), on an individual (i.e. session-to-session, where session indicates BIDS "run") basis.
% Both 2D and 3D options are automatially chosen, as well as processing of time-series (if available), such as motion correction and outlier
% exclusion. This module performs the following initialization/admin steps:
%
% A. Check existence and validity of ASL
% B. Manage mutex state processing step
% C. Cleanup before rerunning
% D. ASL processing & quantification parameters
% E. Backward and forward compatibility of filenames
% F. Split ASL and M0 within the ASL time series
% G. DeltaM and CBF parsing - check if all/some volumes are deltams or CBFs
%
% This module has the following submodules/wrappers:
%
% - `010_TopUpASL`           - FSL TopUp geometric distortion correction (if M0 images with reversed phase-encoding are present)
% - `020_RealignASL`         - If time-series are present, motion correction and outlier exclusion (ENABLE)
% - `030_RegisterASL`        - Registration of ASL to T1w anatomical images (if lacking, to MNI images)
% - `040_ResampleASL`        - Resample ASL images to standard space
% - `050_PreparePV`          - Create partial volume images in ASL space with ASL resolution
% - `060_ProcessM0`          - M0 image processing
% - `070_CreateAnalysisMask` - Create mask using FoV, vascular outliers & susceptibility atlas
% - `080_Quantification`     - CBF quantification
% - `090_VisualQC_ASL`       - Generate QC parameters & images
% - `100_WADQC`              - QC for WAD-QC DICOM server (OPTIONAL)
%
% EXAMPLE: [~, x] = xASL_module_ASL(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2024 ExploreASL


%% Admin
[x] = xASL_init_SubStructs(x);
x = xASL_init_InitializeMutex(x, 'ASL'); % starts mutex locking process to ensure that everything will run only once

if x.mutex.bAnyModuleLocked
    % If any module is locked for this subject, we skip this module for
    % this subject    
    result = true;
    return;
end


%% A. Check existence and validity of ASL
[x, bSkip] = xASL_module_ASL_CheckASL(x);
if bSkip; result = true; return; end
x = xASL_init_FileSystem(x); % Do this only here, to save time when skipping this module
oldFolder = cd(x.dir.SESSIONDIR); % Change working directory to make sure that unspecified output will go there...

tempASL = xASL_io_Nifti2Im(x.P.Path_ASL4D);
if isempty(tempASL) || max(tempASL(:))==0 || numel(unique(tempASL(:)))==1
    error('Invalid ASL image');
else
    nVolumes = size(tempASL, 4);
end


%% B. Manage mutex state processing step
StateName{ 1} = '010_TopUp';
StateName{ 2} = '020_RealignASL';
StateName{ 3} = '030_RegisterASL';
StateName{ 4} = '040_ResampleASL';
StateName{ 5} = '050_PreparePV';
StateName{ 6} = '060_ProcessM0';
StateName{ 7} = '070_CreateAnalysisMask';
StateName{ 8} = '080_Quantification';
StateName{ 9} = '090_VisualQC_ASL';
StateName{10} = '100_WADQC';

if ~x.mutex.HasState('999_ready')
    bO = true; % generate output, some processing has and some has not been yet done
else
    bO = false; % skip output, as all processing has been performed
end


%% C. Cleanup before rerunning
bCompleteRerun = false;
if ~x.mutex.HasState(StateName{3}) || ~x.mutex.HasState(StateName{4})
    % If we rerun the ASL module, then clean it fully for proper rerunning
    % This function cleans all ASL sessions, so only run this (once) for the first session
    xASL_adm_CleanUpBeforeRerun(x.dir.xASLDerivatives, 2, false, false, x.P.SubjectID, x.P.SessionID);
	
	% If we rerun completely, then we do not the reload QC for that session, otherwise we load it from QC.json if present.
	bCompleteRerun = true;
end

x = xASL_adm_LoadX(x, [], true); % assume x.mat is newer than x

% Check and remove all outdated QC fields that are not used anymore
x = xASL_qc_CleanOldQC(x, bCompleteRerun);



%% D. ASL processing & quantification parameters
x = xASL_module_ASL_ParseParameters(x, bO, nVolumes);


%% E. Backward and forward compatibility of filenames
if ~xASL_exist(x.P.Path_M0, 'file')
    % First try to find one with a more BIDS-compatible name & rename it (QUICK & DIRTY FIX)
    FileList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(.*M0.*run.*|.*run.*M0.*)\.nii$','FPList',[0 Inf]);
    if ~isempty(FileList)
        xASL_Move(FileList{1}, x.P.Path_M0);
    end
end
% Do the same for the ancillary files
FileList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(.*M0.*run.*|.*run.*M0.*)_parms\.mat$','FPList',[0 Inf]);
if ~isempty(FileList)
    xASL_Move(FileList{1}, x.P.Path_M0_parms_mat);
end
FileList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(.*M0.*run.*|.*run.*M0.*)\.json$','FPList',[0 Inf]);
if ~isempty(FileList)
    xASL_Move(FileList{1}, fullfile(x.dir.SESSIONDIR, 'M0.json'));
end

% Backward compatibility
if xASL_exist(x.P.Pop_Path_qCBF_untreated,'file')
    xASL_Move(x.P.Pop_Path_qCBF_untreated, x.P.Pop_Path_qCBF, true);
end


%% F. Split ASL and M0 within the ASL time series
% Run this when the data hasn't been touched yet
% The first three states are here, because the first two are run only conditionally

% NOTE THAT THIS PART IS BEING PHASED OUT, AND NOW PROCESSED IN
% `xASL_bids_parseM0` (BIDS->Legacy conversion). We keep this here only for
% backward compatibility, for data imported a long time ago
if ~x.mutex.HasState(StateName{1}) && ~x.mutex.HasState(StateName{2}) && ~x.mutex.HasState(StateName{3})
	% Split the M0 and dummy scans from the ASL time-series
	xASL_io_SplitASL(x.P.Path_ASL4D, x.modules.asl.M0PositionInASL4D, x.modules.asl.DummyScanPositionInASL4D);
	
	% Do the same for the ancillary files
	FileList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(.*ASL4D.*run.*|.*run.*ASL4D.*)_parms\.mat$','FPList',[0 Inf]);
	if ~isempty(FileList)
		xASL_Move(FileList{1}, x.P.Path_ASL4D_parms_mat);
	end
	FileList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(.*ASL4D.*run.*|.*run.*ASL4D.*)\.json$','FPList',[0 Inf]);
	if ~isempty(FileList)
		xASL_Move(FileList{1}, fullfile(x.dir.SESSIONDIR, 'ASL4D.json'));
	end
end

%% G. DeltaM and CBF parsing - check if all/some volumes are deltams or CBFs
% If TSV file exist
% We don't have a subtraction image by default
x.modules.asl.bContainsSubtracted = false;
% Load TSV file
if xASL_exist(x.P.Path_ASL4Dcontext, 'file')
	aslContext = xASL_tsvRead(x.P.Path_ASL4Dcontext);
	bidsPar = xASL_bids_Config;
	% Check for presence of deltaM subtraction volumes
	if ~isempty(regexpi(strjoin(aslContext(2:end)),bidsPar.stringDeltaM, 'once'))
		x.modules.asl.bContainsSubtracted = true;
	end

	% Check for presence of CBF volumes
	if ~isempty(regexpi(strjoin(aslContext(2:end)),bidsPar.stringCbf, 'once'))
		x.modules.asl.bContainsSubtracted = true;
	end
else
	% In case the ASL-context file is missing we set containsDeltaM to true for all NIfTIs with a single volume only
	if xASL_exist(x.P.Path_ASL4D, 'file')
		niftiASL = xASL_io_ReadNifti(x.P.Path_ASL4D);
		if size(niftiASL.dat,4) == 1
			warning('ASL4Dcontext.tsv is missing, assuming that this is a single deltaM volume');
			x.modules.asl.bContainsSubtracted = true;
		end
	end
end



%% ========================================================================================================================
%% 1 TopUp (WIP, only supported if FSL installed)
Path_RevPE = xASL_adm_GetFileList(x.dir.SESSIONDIR, '^(ASL4D|M0).*RevPE\.nii$', 'FPList', [0 Inf]);

iState = 1;
if xASL_exist(x.P.Path_M0,'file') && ~isempty(Path_RevPE)
    if ~x.mutex.HasState(StateName{iState}) || ~xASL_exist(fullfile(x.dir.SESSIONDIR, 'TopUp_fieldcoef.nii'),'file')

        xASL_adm_DeleteFileList(x.dir.SESSIONDIR,'^(B0|Field|TopUp|Unwarped).*$',[],[0 Inf]); % delete previous TopUp stuff first
        bSuccess = xASL_fsl_TopUp(x.dir.SESSIONDIR, 'asl', x, x.P.Path_ASL4D);

        if bSuccess
            x.mutex.AddState(StateName{iState});
            x.mutex.DelState(StateName{iState+1});
        else
            warning('TopUp failed, may affect results rest of pipeline');
        end
    elseif bO; fprintf('%s\n',[StateName{iState} 'has already been performed, skipping...']);
    end
elseif bO; fprintf('%s\n','No TopUp scans available, skipping...');
end


%% ========================================================================================================================
%% 2    Motion correction ASL (& center of mass registration)
iState = 2;
if ~x.modules.asl.motionCorrection
    if bO; fprintf('%s\n','Motion correction was disabled, skipping'); end
    x.mutex.AddState(StateName{iState});
elseif ~x.mutex.HasState(StateName{iState})

        % Remove previous files
        DelList = {x.P.Path_mean_PWI_Clipped_sn_mat, x.P.File_ASL4D_mat, x.P.File_despiked_ASL4D, x.P.File_despiked_ASL4D_mat, 'rp_ASL4D.txt'};
        for iD=1:length(DelList)
            FileName = fullfile(x.dir.SESSIONDIR, DelList{iD});
            xASL_delete(FileName);
        end

%         % First, solve dimensionality (in case there are empty dims, that need restructuring)
        
        % Then, check matrix size: throw error if 2D data with 3 dimensions only

        if nVolumes>1 && ~x.modules.asl.bContainsSubtracted && mod(nVolumes, 2) ~= 0
            error('Uneven number of control-label frames, either incomplete pairs or M0 image in time-series!');
        end

        % Before motion correction, we align the images with ACPC
        PathB0 = fullfile(x.dir.SESSIONDIR, 'B0.nii');
        PathRevPE = fullfile(x.dir.SESSIONDIR, 'ASL4D_RevPE.nii');
        PathField = fullfile(x.dir.SESSIONDIR, 'Field.nii');
        PathFieldCoeff = fullfile(x.dir.SESSIONDIR, 'TopUp_fieldcoef.nii');
        PathUnwarped = fullfile(x.dir.SESSIONDIR, 'Unwarped.nii');
        
        OtherList = {x.P.Path_M0, PathB0, PathRevPE, PathField, PathFieldCoeff, PathUnwarped, x.P.Path_ASL4D_ORI}; % all other files will be created
        if x.settings.bAutoACPC
            xASL_im_CenterOfMass(x.P.Path_ASL4D, OtherList, 10); % set CenterOfMass to lower accepted distance for when rerunning wrong registration
        end
        
        if nVolumes>1 % skip realignment if there are too few volumes
            xASL_wrp_RealignASL(x);
        end

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState(StateName{iState+1});
else
    xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
    if bO; fprintf('%s\n',[StateName{iState} ' session has already been performed, skipping...']); end
end


%% ========================================================================================================================
%% 3    Registration ASL sessions to T1w
iState = 3;
if ~x.mutex.HasState(StateName{iState})

    x = xASL_wrp_RegisterASL(x);

	% And saved the cleaned up version
	xASL_adm_SaveX(x);

    x.mutex.AddState(StateName{iState});
    xASL_adm_CompareDataSets([], [], x); % unit testing
    x.mutex.DelState(StateName{iState+1});
    x.mutex.DelState(StateName{iState+2}); % also dependent on ASL registration
else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if  bO
		fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']);
	end
end


%% ========================================================================================================================
%% 4    Resample ASL data
iState = 4;
if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-1})

    xASL_wrp_ResampleASL(x);

    x.mutex.AddState(StateName{iState});
    xASL_adm_CompareDataSets([], [], x); % unit testing
    x.mutex.DelState(StateName{iState+1});
    x.mutex.DelState(StateName{iState+2});
    x.mutex.DelState(StateName{iState+3});
    x.mutex.DelState(StateName{iState+4});
    x.mutex.DelState(StateName{iState+5});
else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if  bO
		fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']);
	end
end


%% ========================================================================================================================
%% 5    Resolution estimation & prepare partial volume maps
iState = 5;
% We currently only do this once per subject, as we do not anticipate resolution differences between ASL sessions
% We skip this in case there are no structural segmentations
if xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
    % Or if this is the first session, redo this (if ExploreASL is rerun)
    if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-1})

        bStandardSpace = strcmp(x.P.SessionID,x.SESSIONS{1});
        % in standard space we only have 1 set of PV maps in ASL resolution
        % per T1w, therefore only create these for the first ASL session
        x = xASL_wrp_PreparePV(x, bStandardSpace);

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing

	else
		xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
		if  bO; fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']); end
    end
elseif  bO; fprintf('%s\n',['there were no pGM/pWM, skipping ' StateName{iState} '...']);
    x = xASL_adm_DefineASLResolution(x);
end


%% ========================================================================================================================
%% 6    Process M0
iState = 6;
if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-3})
        if strcmp(x.Q.M0, 'Absent')
			if xASL_exist(x.P.Path_M0,'file')
				warning('M0 NIfTI detected, but skipping M0 processing because x.Q.M0="Absent"');
			else
				fprintf('%s\n', 'x.Q.M0="Absent", skipping M0 processing');
			end
        elseif xASL_exist(x.P.Path_M0,'file')

            xASL_wrp_ProcessM0(x);

            x.mutex.AddState(StateName{iState});
            xASL_adm_CompareDataSets([], [], x); % unit testing
            x.mutex.DelState(StateName{iState+1});
			x.mutex.DelState(StateName{iState+2});
        elseif  bO; fprintf('%s\n',[StateName{iState} ' skipped, because no M0 available']);
        end
elseif  xASL_exist(x.P.Path_M0,'file')
		xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
elseif  bO; fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']);
end


%% ========================================================================================================================
%% 7    Create analysis mask
iState = 7;
if ~x.mutex.HasState(StateName{iState})
    if xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
        if size(xASL_io_Nifti2Im(x.P.Path_PWI),4)>1
            warning('Skipped vascular masking because we had too many images');
		else
            xASL_wrp_CreateIndividualMask(x);

            x.mutex.AddState(StateName{iState});
            xASL_adm_CompareDataSets([], [], x); % unit testing
        end
    else
        if  bO;fprintf('%s\n',['T1w-related images missing, skipping ' StateName{iState}]);end
    end
else
    xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
    if  bO;fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']);end
end


%% ========================================================================================================================
%% 8    Quantification

%% #1543 Here we call xASL_im_MergePWI4D as first function before any quantification (so the user can move it somewhere else/it can be a separate step in a GUI)
%% This is not part of the quantification, just as HadamardDecoding is also not part of the quantification, it is a resorting/reconcatenation of image & parameter matrices


% Quantification is performed here according to ASL consensus paper (Alsop, MRM 2016)
% Including PVC
iState = 8;
if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-4})

    fprintf('%s\n','Quantifying ASL:   ');
    % If BASIL quantification will be performed, only native space analysis is possible
    if isfield(x.Q, 'bUseBasilQuantification') && x.Q.bUseBasilQuantification
        % Quantification in native space only for BASIL

        xASL_wrp_Quantify(x, x.P.Path_PWI4D, x.P.Path_CBF, x.P.Path_rM0, x.P.Path_SliceGradient);
        % Transform the output volumes of BASIL to standard space (running
        % it in standard space takes too much time and the
        % noise-distributions may differ due to the B-spline interpolation
        InputPaths = {x.P.Path_CBF, x.P.Path_ATT, x.P.Path_ABV};
        OutputPaths = {x.P.Pop_Path_qCBF, x.P.Pop_Path_ATT, x.P.Pop_Path_ABV}; 
        xASL_spm_deformations(x, InputPaths, OutputPaths, 4, [], [], x.P.Path_y_ASL);
    else
        % Quantification in standard space:
        xASL_wrp_Quantify(x);
        % Quantification in native space:
        xASL_wrp_Quantify(x, x.P.Path_PWI4D, x.P.Path_CBF, x.P.Path_rM0, x.P.Path_SliceGradient);
    end
    
	% allow 4D quantification as well, storing CBF4D. This is currently not implemented for multi-PLD/multi-TE (BASIL/FABBER)
	if isfield(x.Q, 'SaveCBF4D') && x.Q.SaveCBF4D
        if x.Q.nUniqueInitial_PLD>1 || x.modules.asl.bQuantifyMultiTE
            warning('Saving CBF4D was requested but not implemented yet for multi-PLD or multi-TE, skipping');
		else
            if size(xASL_io_Nifti2Im(x.P.Path_ASL4D), 4) == 1
                warning('x.Q.SaveCBF4D was requested but only one volume exists, skipping');
            else
                fprintf('%s\n','Quantifying CBF4D in native space');
                xASL_wrp_Quantify(x, x.P.Path_PWI4D, x.P.Path_qCBF4D, x.P.Path_rM0, x.P.Path_SliceGradient);
    
                fprintf('%s\n','Quantifying CBF4D in standard space');
                xASL_wrp_Quantify(x, x.P.Pop_Path_PWI4D, x.P.Pop_Path_qCBF4D);
            end
        end
	end
	
	if x.modules.asl.bPVCNativeSpace
		fprintf('%s\n','Partial volume correcting ASL in native space:   ');
		if xASL_exist(x.P.Path_PVgm,'file') && xASL_exist(x.P.Path_PVwm,'file') && xASL_exist(x.P.Path_CBF,'file')
			xASL_wrp_PVC(x);
		else
			warning(['Skipped PV: ' x.P.Path_PVgm ', ' x.P.PathPVwm ', or ' x.P.Path_CBF ' missing']);
		end
	end

    x.mutex.AddState(StateName{iState});
    xASL_adm_CompareDataSets([], [], x); % unit testing
    x.mutex.DelState(StateName{iState+1});
    x.mutex.DelState(StateName{iState+2});
else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if  bO; fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']); end
end


%% ========================================================================================================================
%% 9    Visual QC
iState = 9;
if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-2})
	xASL_wrp_VisualQC_ASL(x);
	x.mutex.AddState(StateName{iState});
else
	if  bO; fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']); end
end


%% ========================================================================================================================
%% 10    WAD-QC
iState = 10;
if ~x.mutex.HasState(StateName{iState}) && x.DoWADQCDC
    xASL_qc_WADQCDC(x, x.iSubject, 'ASL');
    x.mutex.AddState(StateName{iState});
elseif x.mutex.HasState(StateName{iState}) && bO
    fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);
end


%% ========================================================================================================================
%% 999 Ready
x.mutex.AddState('999_ready');
cd(oldFolder);

x.mutex.Unlock();
x.result = true;
result = true;
close all;

end


%% ========================================================================================================================
%% ========================================================================================================================
function [x, bSkip] = xASL_module_ASL_CheckASL(x)
%% xASL_module_ASL_CheckASL Check if ASL exists, otherwise skip this module

    bSkip = false;
    x.P.Path_ASL4D = fullfile(x.dir.SESSIONDIR, 'ASL4D.nii');
    x.P.Path_ASL4D_json = fullfile(x.dir.SESSIONDIR, 'ASL4D.json');
    x.P.Path_ASL4D_parms_mat = fullfile(x.dir.SESSIONDIR, 'ASL4D_parms.mat');

    if ~xASL_exist(x.P.Path_ASL4D, 'file')
        % First try to find one with a more BIDS-compatible name & rename it (QUICK & DIRTY FIX)
        FileList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(?i)ASL4D.*\.nii$');

        if ~isempty(FileList)
            xASL_Move(FileList{1}, x.P.Path_ASL4D);
            [Fpath, Ffile] = xASL_fileparts(x.P.Path_ASL4D);
            jsonPath = fullfile(Fpath, [Ffile '.json']);
            parmsPath = fullfile(Fpath, [Ffile '_parms.mat']);
            if exist(jsonPath,'file')
                xASL_Move(jsonPath, x.P.Path_ASL4D_json);
            end
            if exist(parmsPath,'file')
                xASL_Move(parmsPath, x.P.Path_ASL4D_parms_mat);
            end
        else
            fprintf('%s\n',['No ASL found, skipping: ' x.dir.SESSIONDIR]);
            bSkip = true;
        end
    end

end






%% ========================================================================================================================
%% =============================================================================================
function [x] = xASL_module_ASL_ParseParameters(x, bO, nVolumes)
%% xASL_module_ASL_ParseParameters Manage ASL-specific quantification and processing parameters
% 1. Load ASL parameters (inheritance principle)
% 2. Default ASL processing settings in the x.modules.asl field
% 3. Multi-parameter parsing TE-PLD-LD
% 4. TimeEncoded parsing
% 5. Session merging
% 6. Default quantification parameters in the Q field
% 7. Define sequence (educated guess based on the Q field)


%% 1. Load ASL parameters (inheritance principle)
[~, x] = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x, bO);


%% 2. Default ASL processing settings in the x.modules.asl field
if ~isfield(x.modules.asl,'bPVCNativeSpace') || isempty(x.modules.asl.bPVCNativeSpace)
	x.modules.asl.bPVCNativeSpace = 0;
end

% Initialize the DummyScan and M0 position fields - by default empty
if ~isfield(x.modules.asl,'DummyScanPositionInASL4D') 
	x.modules.asl.DummyScanPositionInASL4D = [];
end
if ~isfield(x.modules.asl,'M0PositionInASL4D') 
	x.modules.asl.M0PositionInASL4D = [];
end
if ~isfield(x.modules.asl,'motionCorrection')
    x.modules.asl.motionCorrection = 1;
end
if ~isfield(x,'DoWADQCDC')
    x.DoWADQCDC = false; % default skip WAD-QC stuff
end


%% 3. Multi-parameter parsing
parNames = {'EchoTime' 'Initial_PLD' 'LabelingDuration'};
parAbbreviation = {'TE' 'PLD' 'LD'};
parTolerance = {0.001 0 0};
parLowestValue = {0 0 0};

for iPar=1:length(parNames)

    if ~isfield(x.Q, parNames{iPar})
        % if the field is completely missing
        warning(['Missing field x.Q.' parNames{iPar} ', this may be needed for quantification']);
        % After the warning, we default to single-parameter processing
        fprintf('%s\n', ['Defaulting to single-' parAbbreviation{iPar} ' ASL processing']);
	    x.Q.(['nUnique' parNames{iPar}]) = 1;
        x.modules.asl.(['bQuantifyMulti' parAbbreviation{iPar}]) = false;
    
    elseif isempty(x.Q.(parNames{iPar})) || ~isnumeric(x.Q.(parNames{iPar}))
        % if the field exists but is illegal
        warning(['Illegal field x.Q.' parNames{iPar} ', this should not be empty and should contain numerical values']);
        % After the warning, we default to single-parameter processing
        fprintf('%s\n', ['Defaulting to single-' parAbbreviation{iPar} ' ASL processing']);
	    x.Q.(['nUnique' parNames{iPar}]) = 1;
        x.modules.asl.(['bQuantifyMulti' parAbbreviation{iPar}]) = false;
    else
        % Check, with allowed tolerance (0 is without tolerance) what the unique parameters are
        x.Q.(['unique' parNames{iPar}]) = uniquetol(x.Q.(parNames{iPar}), parTolerance{iPar});
    
        x.Q.(['nUnique' parNames{iPar}]) = length(x.Q.(['unique' parNames{iPar}]));
        % Obtain the number of unique parameters
    
        if sum(x.Q.(parNames{iPar})<=parLowestValue{iPar})>0
            warning(['x.Q.' parNames{iPar} ' should be larger than ' num2str(parLowestValue{iPar})]);
        end
    
        % Handle different number of parameters (== potentially different ASL sequences)
        if x.Q.(['nUnique' parNames{iPar}]) == 1
            fprintf('%s\n', ['Single-' parAbbreviation{iPar} ' ASL detected']);
        elseif x.Q.(['nUnique' parNames{iPar}]) == 2
            fprintf('%s\n', ['Dual-' parAbbreviation{iPar} ' ASL detected']);
            
            if strcmp(parAbbreviation{iPar}, 'TE')
                warning('Dual-echo ASL detected, dual-echo ASL processing not yet implemented');
                fprintf('%s\n', 'We will process these ASL images as if they have a single echo, using the first echo only');
                fprintf('%s\n', 'Assuming that this is a dual-echo ASL+BOLD sequence');
    
                x.Q.EchoTime = min(x.Q.EchoTime);
		        % Update the other fields as well
		         x.Q.uniqueEchoTime = x.Q.EchoTime;
		         x.Q.nUniqueEchoTime = 1;
            end
    
        elseif x.Q.(['nUnique' parNames{iPar}])>2
            fprintf('%s\n', ['Multiple ' parNames{iPar} 's detected, processing this as multi-' parAbbreviation{iPar} ' ASL']);
            fprintf('%s\n', 'Note that this feature is still under development');
    
            if strcmp(parAbbreviation{iPar}, 'PLD')
                if isfield(x.modules.asl, 'bQuantifyMultiPLD') && ~x.modules.asl.bQuantifyMultiPLD
                    fprintf('%s\n', 'Multi PLDs detected, but multi-PLD quantification set OFF.');
                else
                % In case it wasn't defined or it was defined as true, we can set it to true
                x.modules.asl.bQuantifyMultiPLD = true;
                end
            elseif strcmp(parAbbreviation{iPar}, 'LD')
                if x.modules.asl.bQuantifyMultiPLD
                    fprintf('%s\n', 'Multi labeling durations detected, taking this into account for multi-PLD ASL processing');
                else
                    warning('Multi labeling durations detected but multiPLD-quantification is turned off, is this correct?');
                end
            end
        end
    
        % Now we print to the screen which unique parameter-values we detected
        fprintf('%s\n', ['Detected the following unique ' parNames{iPar} 's (ms):'])
        for iVol = 1:x.Q.(['nUnique' parNames{iPar}])
            fprintf('%.2f, ', round(x.Q.(['unique' parNames{iPar}])(iVol), 4));
        end
        fprintf('\n');
    end
    
    % If multi-parameter data is detected, we switch on multi-parameter quantification by default. Unless this has been deactived in dataPar.json
    if x.Q.(['nUnique' parNames{iPar}])>1
	    if ~isfield(x.modules.asl, ['bQuantifyMulti' parAbbreviation{iPar}]) || isempty(x.modules.asl.(['bQuantifyMulti' parAbbreviation{iPar}]))
		    % By default, initialize as Multi-par quantification as true
		    x.modules.asl.(['bQuantifyMulti' parAbbreviation{iPar}]) = true;
	    end
    else
	    x.modules.asl.(['bQuantifyMulti' parAbbreviation{iPar}]) = false;
    end
    
    % Deal with too short vectors (e.g., repetitions in the case of single-parameter)
    nPar = length(x.Q.([parNames{iPar}]));
    if mod(nVolumes, nPar)~=0
        error(['Number of ' parNames{iPar} 's (n=' xASL_num2str(nPar) ') does not fit in nVolumes (n=' xASL_num2str(nVolumes) ')']);
        % If we don't issue an error here, repmat will crash
    end
    
    factorPar = nVolumes/nPar;
    x.Q.(parNames{iPar}) = repmat(x.Q.(parNames{iPar})(:), [factorPar 1]);
    nPar = length(x.Q.(parNames{iPar}));
    
    % Number of echoes should now be equal to the number of ASL volumes
    if nPar~=nVolumes
        warning(['Number of ' parNames{iPar} ' (n=' xASL_num2str(nPar) ') should be equal to number of ASL volumes (n=' xASL_num2str(nVolumes) ')']);
    end

end



%% 4. TimeEncoded parsing
% Check if TimeEncoded is defined
if isfield(x.Q, 'TimeEncodedMatrixType') || isfield(x.Q, 'TimeEncodedMatrixSize') || isfield(x.Q, 'TimeEncodedMatrix')
    fprintf(2,'Time encoded data detected, we will process this but this is a new feature that is still under development\n');
    x.modules.asl.bTimeEncoded = true;
	
	% Initialize as empty if missing
	if ~isfield(x.Q, 'TimeEncodedMatrixType')
		x.Q.TimeEncodedMatrixType = [];
	end

	if ~isfield(x.Q, 'TimeEncodedMatrixSize')
		x.Q.TimeEncodedMatrixSize = [];
	end

	if ~isfield(x.Q, 'TimeEncodedMatrix')
		x.Q.TimeEncodedMatrix = [];
	end

	% Check for duplicate definitions
	if isempty(x.Q.TimeEncodedMatrixType) && isempty(x.Q.TimeEncodedMatrix)
		error('Neither TimeEncodedMatrixType nor TimeEncodedMatrix were provided for a TimeEncoded acquisition.')
	end
	
	if isempty(x.Q.TimeEncodedMatrixSize)
		if isempty(x.Q.TimeEncodedMatrix)
			error('TimeEncodedMatrixSize not defined and TimeEncodedMatrix not provided.')
		else
			x.Q.TimeEncodedMatrixSize = size(x.Q.TimeEncodedMatrix,1);
		end
	else
		if (x.Q.TimeEncodedMatrixSize < 2)
			fprintf(2,'TimeEncodedMatrixSize field missing (should be a multiple of 4)...\n');
		end
		if ~isempty(x.Q.TimeEncodedMatrix) && (x.Q.TimeEncodedMatrixSize ~= size(x.Q.TimeEncodedMatrix, 1))
			error('Mismatch between TimeEncodedMatrix size and TimeEncodedMatrixSize parameters.');
		end
	end
else
	x.modules.asl.bTimeEncoded = false;
end

% Check if there is Decoding Matrix as input
%(some datasets will have a decoding matrix that we can use directly in the decoding part)

%% 5. Session merging
% Initialization
nLists = 0;

% SessionMergingList == total list of sessions that need concatenating (can be multiple lists for multiple concatenations)
% Read the list of sessions to merge, defaulting to empty
if ~isfield(x.modules.asl, 'SessionMergingList')
	x.modules.asl.SessionMergingList = {};
else
	% The parameter should be a list of lists, but if a single list is provided, it converts it
	if ~isempty(x.modules.asl.SessionMergingList) && iscell(x.modules.asl.SessionMergingList) && ~iscell(x.modules.asl.SessionMergingList{1})
		x.modules.asl.SessionMergingList = {x.modules.asl.SessionMergingList};
	end
    nLists = numel(x.modules.asl.SessionMergingList);
    fprintf('\n%s', ['-> Detected ' xASL_num2str(nLists) ' list(s) for concatenating sessions']);
end

x.modules.asl.bMergingSessions = 0; % By default, we do not merge any sessions
% sessionsToMerge = session list that we will concatenate for this specific session & run of xASL_module_ASL
% e.g., if we want to concatenate ASL_1 & ASL_2, this will only contain both sessions for the xASL_module_ASL iteration of ASL_2
x.modules.asl.sessionsToMerge = {};


% Find the lists containing the current session
% Note that these can be multiple lists
bSessionListed = false; % It is possible that a single session is in multiple list. That is not illegal
for iList=1:nLists
	if ~isempty(x.modules.asl.SessionMergingList{iList}) && sum(ismember(x.SESSION, x.modules.asl.SessionMergingList{iList}))
		if bSessionListed
			fprintf(['Session ' x.SESSION ' appears in more than one list of sessions to merge.\n']);
		else
			x.modules.asl.sessionsToMerge = x.modules.asl.SessionMergingList{iList};
			bSessionListed = true;
		end

        if strcmp(x.SESSION, x.modules.asl.SessionMergingList{iList}{end})
            % If the current session is the last of the list then we set the merging to TRUE. Otherwise, we keep merging to later

			% If merging is already set to true, it means that there was a previous list to be merged. We report this as a warning and we keep the first fitting list to be merged
			if x.modules.asl.bMergingSessions
				warning(['Session ' x.SESSION ' is at the end of more than one list of sessions to merge. We ignore all such lists but the first.'])
			else
				% We have to assign the list again, because it is possible that sessionsToMerge was already initialized to a list that contained the current session, but that was not ending it
				x.modules.asl.sessionsToMerge = x.modules.asl.SessionMergingList{iList}; 				                                                            
				x.modules.asl.bMergingSessions = 1;
			end
        end
	end
end

if x.modules.asl.bMergingSessions
    if ~isempty(x.modules.asl.sessionsToMerge)
        fprintf('%s', ' and will now concatenate the following sessions:')
        for iSession = 1:numel(x.modules.asl.sessionsToMerge)
            fprintf([' ' x.modules.asl.sessionsToMerge{iSession}]);
        end
        fprintf('\n');
    else
        fprintf('%s\n\n', [' but will not concatenate them here (xASL_module_ASL: ' x.SESSION ')']);
    end
end


%% 6. Default quantification parameters in the Q field
if ~isfield(x.modules.asl,'ApplyQuantification') || isempty(x.modules.asl.ApplyQuantification)
    x.modules.asl.ApplyQuantification = [1 1 1 1 1 1]; % by default we perform scaling/quantification in all steps
elseif length(x.modules.asl.ApplyQuantification)>6
    warning('x.modules.asl.ApplyQuantification had too many parameters');
    x.modules.asl.ApplyQuantification = x.modules.asl.ApplyQuantification(1:6);
elseif length(x.modules.asl.ApplyQuantification)<6
    warning('x.modules.asl.ApplyQuantification had too few parameters, using default 1');
    x.modules.asl.ApplyQuantification(length(x.modules.asl.ApplyQuantification)+1:6) = 1;
end

if ~isfield(x.Q,'BackgroundSuppressionNumberPulses') && isfield(x,'BackgroundSuppressionNumberPulses')
    % Temporary backwards compatibility that needs to go
    x.Q.BackgroundSuppressionNumberPulses = x.BackgroundSuppressionNumberPulses;
end

if ~isfield(x.Q, 'bUseBasilQuantification') || isempty(x.Q.bUseBasilQuantification)
	x.Q.bUseBasilQuantification = false;
    
    if x.modules.asl.bQuantifyMultiPLD || x.modules.asl.bQuantifyMultiTE
        x.Q.bUseBasilQuantification = true;
    end
end

% Manage absent M0
if ~x.modules.asl.ApplyQuantification(5) && ~xASL_exist(x.P.Path_M0) && ~strcmp(x.Q.M0, 'Absent')
    warning('M0 division was disabled & M0 missing, setting M0 to "absent"');
    x.Q.M0 = 'Absent';
end

if strcmp(x.Q.M0, 'Absent')
    fprintf('%s\n', 'x.Q.M0="Absent" so disabling M0 processing');
    x.modules.asl.ApplyQuantification([2, 4, 5]) = 0;
    
    if xASL_exist(x.P.Path_M0)
        warning('Ignoring existing M0 file because of x.Q.M0="Absent"');
    end
end


%% 7. Define sequence (educated guess based on the Q field)
x = xASL_adm_DefineASLSequence(x, bO);


end