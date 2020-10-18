function [result, x] = xASL_module_ASL(x)
%xASL_module_ASL ExploreASL module for ASL processing
%
% FORMAT: [result, x] = xASL_module_ASL(x)
%
% INPUT:
%   x  - x structure containing all input parameters (REQUIRED)
%   x.SUBJECTDIR  -  anatomical directory, containing the derivatives of anatomical images (REQUIRED)
%   x.SESSIONDIR  -  ASL directory, containing the derivatives of perfusion images (REQUIRED)
%
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters (REQUIRED)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL module processes the ASL
% images, i.e. ASL4D, M0, etc (if present), on an individual (i.e. session-to-session, where session indicates BIDS "run") basis.
% Both 2D and 3D options are automatially chosen, as well as processing of time-series (if available), such as motion correction and outlier
% exclusion. This module has the following submodules/wrappers:
%
% - 010_TopUpASL          - FSL TopUp geometric distortion correction (if M0 images with reversed phase-encoding are present)
% - 020_RealignASL        - If time-series are present, motion correction and outlier exclusion (ENABLE)
% - 030_RegisterASL       - Registration of ASL to T1w anatomical images (if lacking, to MNI images)
% - 040_ResliceASL        - Resample ASL images to standard space
% - 050_PreparePV         - Create partial volume images in ASL space with ASL resolution
% - 060_ProcessM0         - M0 image processing
% - 070_CreateAnalysisMask- Create mask using FoV, vascular outliers & susceptibility atlas
% - 080_Quantification    - CBF quantification
% - 090_VisualQC_ASL      - Generate QC parameters & images
% - 100_WADQC             - QC for WAD-QC DICOM server (OPTIONAL)
%
% EXAMPLE: [~, x] = xASL_module_ASL(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL


%% Admin
x = xASL_init_InitializeMutex(x, 'ASL'); % starts mutex locking process to ensure that everything will run only once
result = false;

if ~isfield(x,'ApplyQuantification') || isempty(x.ApplyQuantification)
    x.ApplyQuantification = [1 1 1 1 1]; % by default we perform scaling/quantification in all steps
elseif length(x.ApplyQuantification)>5
    warning('x.ApplyQuantification had too many parameters');
    x.ApplyQuantification = x.ApplyQuantification(1:5);
elseif length(x.ApplyQuantification)<5
    warning('x.ApplyQuantification had too few parameters, using default 1');
    x.ApplyQuantification(length(x.ApplyQuantification)+1:5) = 1;
end

% Only continue if ASL exists
x.P.Path_ASL4D = fullfile(x.SESSIONDIR, 'ASL4D.nii');
x.P.Path_ASL4D_json = fullfile(x.SESSIONDIR, 'ASL4D.json');
x.P.Path_ASL4D_parms_mat = fullfile(x.SESSIONDIR, 'ASL4D_parms.mat');

if ~xASL_exist(x.P.Path_ASL4D, 'file')
    % First try to find one with a more BIDS-compatible name & rename it (QUICK & DIRTY FIX)
    FileList = xASL_adm_GetFileList(x.SESSIONDIR, 'ASL4D.*\.(nii|nii\.gz)');

    if ~isempty(FileList) && isfield(x,'M0PositionInASL4D')
        % skip, managed below
    elseif ~isempty(FileList)
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
        fprintf('%s\n',['No ASL found, skipping: ' x.SESSIONDIR]);
        result = true;
        return;
    end
end

%% Manage an M0 within time series
if isfield(x,'M0PositionInASL4D')
    x.M0PositionInASL4D = xASL_str2num(x.M0PositionInASL4D); % make sure it is numeric when coming from JSON string
    if ~isnumeric(x.M0PositionInASL4D)
        warning('Something wrong with x.M0PositionInASL4D');
    end
    % make sure to split the M0 from the ASL time-series if needed
    xASL_io_SplitASL_M0(x.P.Path_ASL4D, x.M0PositionInASL4D);
end


% Do the same for the ancillary files
FileList = xASL_adm_GetFileList(x.SESSIONDIR, '(.*ASL4D.*run.*|.*run.*ASL4D.*)_parms\.mat','FPList',[0 Inf]);
if ~isempty(FileList)
    xASL_Move(FileList{1}, x.P.Path_ASL4D_parms_mat);
end
FileList = xASL_adm_GetFileList(x.SESSIONDIR, '(.*ASL4D.*run.*|.*run.*ASL4D.*)\.json','FPList',[0 Inf]);
if ~isempty(FileList)
    xASL_Move(FileList{1}, fullfile(x.SESSIONDIR, 'ASL4D.json'));
end

%% Input

x = xASL_init_FileSystem(x); % do this only here, to save time when skipping this module

if ~isfield(x,'SavePWI4D')
    x.SavePWI4D=0;
end
if ~isfield(x,'Q')
    x.Q = struct;
end
if ~isfield(x,'DoWADQCDC')
    x.DoWADQCDC = false; % default skip WAD-QC stuff
end
if ~isfield(x.Q,'BackGrSupprPulses') && isfield(x,'BackGrSupprPulses')
    % Temporary backwards compatibility that needs to go
    x.Q.BackGrSupprPulses = x.BackGrSupprPulses;
end


%% Change working directory to make sure that unspecified output will go there...
oldFolder = cd(x.SESSIONDIR);


%% Import fixes

if ~xASL_exist(x.P.Path_M0, 'file')
    % First try to find one with a more BIDS-compatible name & rename it (QUICK & DIRTY FIX)
    FileList = xASL_adm_GetFileList(x.SESSIONDIR, '(.*M0.*run.*|.*run.*M0.*)\.(nii|nii\.gz)','FPList',[0 Inf]);
    if ~isempty(FileList)
        xASL_Move(FileList{1}, x.P.Path_M0);
    end
end
% Do the same for the ancillary files
FileList = xASL_adm_GetFileList(x.SESSIONDIR, '(.*M0.*run.*|.*run.*M0.*)_parms\.mat','FPList',[0 Inf]);
if ~isempty(FileList)
    xASL_Move(FileList{1}, x.P.Path_M0_parms_mat);
end
FileList = xASL_adm_GetFileList(x.SESSIONDIR, '(.*M0.*run.*|.*run.*M0.*)\.json','FPList',[0 Inf]);
if ~isempty(FileList)
    xASL_Move(FileList{1}, fullfile(x.SESSIONDIR, 'M0.json'));
end


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


%% Delete old files
if ~x.mutex.HasState(StateName{3}) || ~x.mutex.HasState(StateName{4})
    % If we rerun the ASL module, then clean it fully for proper rerunning
    % This function cleans all ASL sessions, so only run this (once) for the first session
    xASL_adm_CleanUpBeforeRerun(x.D.ROOT, 2, false, false, x.P.SubjectID, x.P.SessionID);
	
	% Also clean the QC output files
	x = xASL_adm_LoadX(x, [], true); % assume x.mat is newer than x

	% Clear any previous QC images
	if isfield(x,'Output_im') && isfield(x.Output_im,'ASL')
		x.Output_im = rmfield(x.Output_im,'ASL');
	end

	if isfield(x,'Output') && isfield(x.Output,'ASL')
		x.Output = rmfield(x.Output,'ASL');
	end
	
	% And saved the cleaned up version
	xASL_adm_SaveX(x);
end

if ~isfield(x,'motion_correction')
    x.motion_correction = 1;
end

% Backward compatibility
if xASL_exist(x.P.Pop_Path_qCBF_untreated,'file')
    xASL_Move(x.P.Pop_Path_qCBF_untreated, x.P.Pop_Path_qCBF, true);
end

%% If ASL_4D_parms.mat (ASL parameters) file exist, load and overwrite existing parameters in the xASL file (inheritance principle)
if ~x.mutex.HasState('999_ready')
    bO = true; % generate output, some processing has and some has not been yet done
else
    bO = false; % skip output, as all processing has been performed
end

[~, x] = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x, bO);

if ~isfield(x,'Q')
    x.Q = struct;
    warning('x.Q didnt exist, are quantification parameters lacking?');
end

%% Define sequence (educated guess)
x = xASL_adm_DefineASLSequence(x);

%% -----------------------------------------------------------------------------
%% 1 TopUp (WIP, only supported if FSL installed)
Path_RevPE = xASL_adm_GetFileList(x.SESSIONDIR, '^(ASL4D|M0).*RevPE\.nii$', 'FPList', [0 Inf]);

iState = 1;
if xASL_exist(x.P.Path_M0,'file') && ~isempty(Path_RevPE)
    if ~x.mutex.HasState(StateName{iState}) || ~xASL_exist(fullfile(x.SESSIONDIR, 'TopUp_fieldcoef.nii'),'file')

        xASL_adm_DeleteFileList(x.SESSIONDIR,'^(B0|Field|TopUp|Unwarped).*$',[],[0 Inf]); % delete previous TopUp stuff first
        bSuccess = xASL_fsl_TopUp(x.SESSIONDIR, 'asl', x, x.P.Path_ASL4D);

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

%% -----------------------------------------------------------------------------
%% 2    Motion correction ASL (& center of mass registration)
iState = 2;
if ~x.motion_correction
    if bO; fprintf('%s\n','Motion correction was disabled, skipping'); end
    x.mutex.AddState(StateName{iState});
elseif ~x.mutex.HasState(StateName{iState})

        % Remove previous files
        DelList = {x.P.Path_mean_PWI_Clipped_sn_mat, x.P.File_ASL4D_mat, x.P.File_despiked_ASL4D, x.P.File_despiked_ASL4D_mat, 'rp_ASL4D.txt'};
        for iD=1:length(DelList)
            FileName = fullfile(x.SESSIONDIR, DelList{iD});
            xASL_delete(FileName);
        end

%         % First, solve dimensionality (in case there are empty dims, that need restructuring)
        nVolumes = size(xASL_io_Nifti2Im(x.P.Path_ASL4D),4);

        % Then, check matrix size: throw error if 2D data with 3 dimensions only

        if nVolumes>1 && nVolumes/2~=round(nVolumes/2)
            error('Uneven number of control-label frames, either incomplete pairs or M0 image in time-series!');
        end

        % Before motion correction, we align the images with ACPC
        PathB0 = fullfile(x.SESSIONDIR, 'B0.nii');
        PathRevPE = fullfile(x.SESSIONDIR, 'ASL4D_RevPE.nii');
        PathField = fullfile(x.SESSIONDIR, 'Field.nii');
        PathFieldCoeff = fullfile(x.SESSIONDIR, 'TopUp_fieldcoef.nii');
        PathUnwarped = fullfile(x.SESSIONDIR, 'Unwarped.nii');
        
        OtherList = {x.P.Path_M0, PathB0, PathRevPE, PathField, PathFieldCoeff, PathUnwarped, x.P.Path_ASL4D_ORI}; % all other files will be created
        if x.bAutoACPC
            xASL_im_CenterOfMass(x.P.Path_ASL4D, OtherList, 10); % set CenterOfMass to lower accepted distance for when rerunning wrong registration
        end

        if  nVolumes>1
            % Run motion Correction
            xASL_wrp_RealignASL(x);
        else
            fprintf('%s\n',['Skipping motion correction for ' x.P.SubjectID '_' x.P.SessionID ' because it had only ' num2str(nVolumes) ' 3D frames.']);
        end

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState(StateName{iState+1});
else
    xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
    if bO; fprintf('%s\n',[StateName{iState} ' session has already been performed, skipping...']); end
end


%% -----------------------------------------------------------------------------
%% 3    Registration ASL sessions to T1w
iState = 3;
if ~x.mutex.HasState(StateName{iState})

	% Load the previously saved QC Output
	x = xASL_adm_LoadX(x, [], true); % assume x.mat is newer than x

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




%% -----------------------------------------------------------------------------
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




%% -----------------------------------------------------------------------------
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



%% -----------------------------------------------------------------------------
%% 6    Process M0
iState = 6;
if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-3})
        if xASL_exist(x.P.Path_M0,'file')

            xASL_wrp_ProcessM0(x);

            x.mutex.AddState(StateName{iState});
            xASL_adm_CompareDataSets([], [], x); % unit testing
            x.mutex.DelState(StateName{iState+1});
        elseif  bO; fprintf('%s\n',[StateName{iState} ' skipped, because no M0 available']);
        end
elseif  xASL_exist(x.P.Path_M0,'file')
		xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
elseif  bO; fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']);
end




%% -----------------------------------------------------------------------------
%% 7    Create analysis mask
iState = 7;
if ~x.mutex.HasState(StateName{iState})
    if xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
        if size(xASL_io_Nifti2Im(x.P.Path_PWI),4)>1
            warning('Skipped vascular masking because we had too many images');
		else
            xASL_wrp_CreateAnalysisMask(x);

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



%% -----------------------------------------------------------------------------
%% 8    Quantification
% Quantification is performed here according to ASL consensus paper (Alsop, MRM 2016)
iState = 8;
if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-4})

    fprintf('%s\n','Quantifying ASL:   ');

    % Quantification in standard space:
    if ~xASL_exist(x.P.Pop_Path_PWI,'file')
        warning(['Skipped standard space quantification: ' x.P.Pop_Path_PWI ' missing']);
    else
        xASL_wrp_Quantify(x);
    end
    % Quantification in native space:
    if ~xASL_exist(x.P.Path_PWI,'file')
        warning('Skipped native space quantification: files missing');
    else
        xASL_wrp_Quantify(x, x.P.Path_PWI, x.P.Path_CBF, x.P.Path_rM0, x.P.Path_SliceGradient);
    end

    % allow 4D quantification as well
    if x.SavePWI4D
        if ~xASL_exist(x.P.Path_PWI4D,'file')
            fprintf('%s\n','Skipping, native space PWI4D missing');
        else
            fprintf('%s\n','Quantifying ASL timeseries in native space');
            xASL_wrp_Quantify(x, x.P.Path_PWI4D, x.P.Path_qCBF4D, x.P.Path_rM0, x.P.Path_SliceGradient);
        end
        if ~xASL_exist(x.P.Pop_Path_PWI4D,'file')
            fprintf('%s\n','Skipping, standard space PWI4D missing');
        else
            fprintf('%s\n','Quantifying ASL timeseries in standard space');
            xASL_wrp_Quantify(x, x.P.Pop_Path_PWI4D, x.P.Pop_Path_qCBF4D);
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


%% -----------------------------------------------------------------------------
%% 9    Visual QC
iState = 9;
if ~x.mutex.HasState(StateName{iState}) && x.mutex.HasState(StateName{iState-2})

    xASL_wrp_VisualQC_ASL(x);
    x.mutex.AddState(StateName{iState});

else
	if  bO; fprintf('%s\n',[StateName{iState} ' has already been performed, skipping...']); end
end

%% -----------------------------------------------------------------------------
%% 10    WAD-QC
iState = 10;
if ~x.mutex.HasState(StateName{iState}) && x.DoWADQCDC
    xASL_qc_WADQCDC(x, x.iSubject, 'ASL');
    x.mutex.AddState(StateName{iState});
elseif x.mutex.HasState(StateName{iState}) && bO
    fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);
end

%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
cd(oldFolder);

x.mutex.Unlock();
x.result = true;
result = true;
close all;

end
