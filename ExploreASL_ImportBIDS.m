function ExploreASL_ImportBIDS(studyPath, imParPath, studyParPath, bRunSubmodules, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source, x)
%ExploreASL_ImportBIDS Imports the DICOM or PAR/REC raw data to NIFTIs in ASL-BIDS format
%
% FORMAT: ExploreASL_Import(studyPath, imParPath, studyParPath, bRunSubmodules, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source, x)
%
% INPUT:
%   studyPath           - path to the study directory containing the 'raw' directory with the DICOM files (REQUIRED)
%   imParPath           - path to the JSON file with structure with import parameters, output of ExploreASL_ImportConfig.m (originally)
%                         All other input parameters are configured within this function. (OPTIONAL)
%                         The path is optional, but the file has to be there. Either provided as a full path or a filename in the path,
%                         or default names (case-insensitive) ImPar.json, ImagePar.json are seeked
%   studyParPath        - path to the JSON file with the BIDS parameters relevant for the whole study. These parameters are used
%                         if they cannot be extracted from the DICOMs automatically. (OPTIONAL)
%                         Looking automatically for file studyPar.json
%   bRunSubmodules      - Specify which of the parts should be run (OPTIONAL, DEFAULT [1 1 0])
%                         [1 0 0] - Run the DICOM to NIFTI conversion
%                         [0 1 0] - Run the NIFTI transformation to the proper ASL-BIDS
%                         [0 0 1] - Run the defacing and full anonymization
%   bCopySingleDicoms   - if true, copies a single DICOM with each NIfTI
%                         dataset/ScanType, that can be used to retrieve missing parameters from
%                         the DICOM header, or as dummy DICOM to dump embed data into (e.g. WAD-QC) (DEFAULT=false)
%   bUseDCMTK           - if true, then use DCMTK, otherwise use DICOMINFO from Matlab (DEFAULT=false)
%   bCheckPermissions   - if true, check whether data permissions are set correctly, before trying to read/copy the files (DEFAULT=false)
%   Clone2Source        - if true, then makes a copy of everything it converted to NIfTI.
%                         Can be useful to have a separate source BIDS structure to store
%                         all source NIfTIs, and to keep the derivatives in the
%                         analysisfolder (OPTIONAL, DEFAULT=false)
%   x                   - if x is provided, initialization of ExploreASL is skipped
%
%
% OUTPUT: n/a
%
% OUTPUT FILES:
%   //AnalysisDir/dcm2niiCatchedErrors.(mat|json) - overview of catched dcm2nii errors, or other errors in this function
%   //AnalysisDir/import_log_StudyID_yyyymmdd_hhmmss.txt - diary log of this function
%   //AnalysisDir/import_summary.csv - hence the name
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% Import batch T1, FLAIR, DWI, fMRI, M0, ASL data from dicom 2 NIfTI in ASL-BIDS format and structure.
% Uses dcm2niiX for the conversion, and additionally collects important DICOM header data
% and puts them in .json sidecars to be used with the ExploreASL pipeline.
% This function takes any folder input, but the folder input should be
% specified in the imPar definition. Follow the steps below, for study "MyStudy" located on "//MyDisk":
%
% 1) Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
%    Create a root folder with study ID name, and put the DICOMs in any structure in the raw folder within the study ID root folder
%    Example:
%    imPar.StudyID: MyStudy
%    StudyRoot folder: //MyDisk/MyStudy
%    Raw folder containing DICOMs: //MyDisk/MyStudy/raw
% 2) Make sure that your DICOM data has any structure that can be retrieved
%    from the folder and/or file names. This function doesn't yet read the DICOM headers
%    For a quick and dirty (but actually slow) function that converts a
%    DICOM folder/file structure into readable format, first run
%    ConvertDicomFolderStructure_CarefulSlow.m. This will read each DICOM
%    individually, and put it in a folder with the name identical to the
%    DICOMs SeriesName/ProtocolName.
% 3) Once you have all DICOMs in folderstructure with identifyable names
%    inside //MyDisk/MyStudy/raw, set up the folderstructure in
%    ExploreASL_ImportConfig.m. This setup uses the SPM form of regular
%    expressions, which can be daunting at first, but are very flexible.
%    Easiest is to study other examples, before creating your own.
%    For this example, let's say we have //MyDisk/MyStudy/raw/ScanType/SubjectName
%    because we downloaded our data from XNAT, ordered per ScanType first,
%    and then per subject.
%
%    BRIEF EXPLANATION:
%    Let's suppose we don't have sessions (only a single structural and functional scan per subject)
%    The names of our scans comes out of XNAT as '3D_FLAIR_eyesClosed', 'T1w_MPRAGE' and 'PCASL_10_min'
%
%    and the subject names are 'MyStudy001' .. 'MyStudy002' .. etc.
%
%    imPar.folderHierarchy   - contains a a cell array of regular expressions, with each cell specifying a directory layer/level
%                              the parts within brackets () tell the script that this is a token (i.e. subject, session, ScanType)
%                              Example:
%                              imPar.folderHierarchy = {'^(3D_FLAIR|T1w|PCASL).*', '^(Sub-\d{3})$'};
%                              here we say that there are two folder layers '', separated by comma ,
%                              where the names between brackets are used to define what is what.
%                              ^ means that the foldername has to start with the following, $ means that the previous has to be the end of the foldername
%                              .* means anything, anylength, \d{3} means three digits
%    imPar.tokenOrdering     - defines which tokens are captured by the brackets () in imPar.folderHierarchy: position 1==subject, 2==visit, 3==session, 4==ScanType
%                              Example:
%                              imPar.tokenOrdering = [2 3 0 1]; stating that subject is the 2nd token, visit is the 3rd token, session has no token (i.e. no session) and ScanType is the 1st token
%    imPar.tokenVisitAliases - cell array that defines the aliases for the Visits, i.e. it tells the script which scans are which timepoint/visit.
%                              Similar as explained below for ScanAliases.
%                              First column contains the names that are
%                              recognized in raw DICOM folders for visits,
%                              second column how it is named in NIfTI
%                              structure (should be _1 _2 _3 etc).
%                              Example:
%                              imPar.tokenVisitAliases = {'Screening','_1'; 'Month_12','_2'; 'Month_24','_3'; 'Month_36','_4'; 'Month_48','_5'};
%                              Note that if you specify tokenVisitAliases, the folders will receive
%                              the indices (e.g. _1 _2 _3), or even _1 only with a single Visit). If you don't specify
%                              them, they will not get this postfix.
%    imPar.tokenScanAliases  - cell array that defines the aliases for the ScanTypes, i.e. it tells the script which scans are which ScanType.
%                              First column should contain regular expression corresponding with the matching criteria in imPar.folderHierarchy
%                              whereas the second column contains the
%                              alias. Following valid aliases exist:
%                              'T1' 'FLAIR' 'ASL4D' 'M0' 'ASL4D_RevPE' 'func' 'func_NormPE' 'func_RevPE' 'dwi' 'dwi_RevPE' 'DSC4D'
%                              Example:
%                              imPar.tokenScanAliases = {'^3D_FLAIR$', 'FLAIR'; '^T1w$', 'T1'; '^PCASL$', 'ASL4D'};
%    imPar.tokenSessionAliases-same as tokenScanAliases but for sessions
%                              Example:
%                              imPar.tokenSessionAliases = {}; as we don't have sessions
%    imPar.bMatchDirectories - true if the last layer is a folder, false if the last layer is a filename (as e.g. with PAR/REC, enhanced DICOMs)
%
% EXAMPLE: ExploreASL_ImportBIDS('//MyDisk/MyStudy');
%          ExploreASL_ImportBIDS('//MyDisk/MyStudy','imPar.json','studyHiQ.json');
% __________________________________
% Copyright 2015-2020 ExploreASL

%% 1. Initialize the parameters

% First do the basic parameter admin and initialize the default values
if nargin < 1 ||  isempty(studyPath)
	error('The studyPath needs to be defined');
end

if strcmp(studyPath(end),'\') || strcmp(studyPath(end),'/')
    studyPath = studyPath(1:end-1); % bugfix
end

% Check the imagePar input file
if nargin < 2 || isempty(imParPath)
	% If the path is empty, then try to find impar.json or imagepar.json
	fListImPar = xASL_adm_GetFileList(studyPath,'^(im|image|Im|Image)(Par|par).json$', 'List', [], 0);
	if length(fListImPar) < 1
		error('Could not find the ImPar.json file');
	end
	imParPath = fullfile(studyPath,fListImPar{1});
else
	[fpath, ~, ~] = fileparts(imParPath);
	if isempty(fpath)
		imParPath = fullfile(studyPath,imParPath);
	end
end

% Find the studyPar input file
if nargin < 3 || isempty(studyParPath)
	% If the path is empty, then try to find impar.json or imagepar.json
	fListStudyPar = xASL_adm_GetFileList(studyPath,'^(study|Study)(Par|par).json$', 'List', [], 0);
	if length(fListStudyPar) < 1
		warning('Could not find the StudyPar.json file');
	else
		studyParPath = fullfile(studyPath,fListStudyPar{1});
	end
else
	[fpath, ~, ~] = fileparts(studyParPath);
	if isempty(fpath)
		studyParPath = fullfile(studyPath,studyParPath);
	end
end

if nargin<4 || isempty(bRunSubmodules)
	bRunSubmodules = [1 1 0];
else
	if length(bRunSubmodules) ~= 3
		error('bRunSubmodules must have length 3');
	end
end

if nargin<5 || isempty(bCopySingleDicoms)
    bCopySingleDicoms = false; % by default don't copy dicoms for anonymization reasons
end

if nargin<6 || isempty(bUseDCMTK)
    bUseDCMTK = true; % default set to using DCM-TK
elseif ~bUseDCMTK && isempty(which('dicomdict'))
    error('Dicomdict missing, image processing probably not installed, try DCMTK instead');
end

if nargin<7 || isempty(bCheckPermissions)
    if isunix
        bCheckPermissions = true;
    else
        bCheckPermissions = false;
    end
end

if nargin<8 || isempty(bClone2Source)
    bClone2Source = false;
end

if nargin<9 || isempty(x)
    x = ExploreASL_Initialize('',0); % only initialize ExploreASL if this wasnt initialized before
end

%% 2. Initialize the setup of the dicom2nii conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the imPar field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fpath, fname, fext] = fileparts(studyPath);

% Load the imPar from the file
imPar = spm_jsonread(imParPath);

% Specify paths
if ~isfield(imPar,'studyID') || isempty(imPar.studyID)
	imPar.studyID = [fname fext];
end
if ~isfield(imPar,'AnalysisRoot') || isempty(imPar.AnalysisRoot) 
	imPar.AnalysisRoot = fpath;
end
if ~isfield(imPar,'RawRoot') || isempty(imPar.RawRoot)
	imPar.RawRoot = fpath;
end
if ~isfield(imPar,'BidsRoot') || isempty(imPar.BidsRoot)
	imPar.BidsRoot = fpath;
end

% Finalize the directories
imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID, 'raw'); % default name
imPar.AnalysisRoot = fullfile(imPar.AnalysisRoot,imPar.studyID,'analysis');
imPar.BidsRoot = fullfile(imPar.BidsRoot,imPar.studyID,'bids');

% Specify the tokens
if ~isfield(imPar,'folderHierarchy')
	imPar.folderHierarchy = {}; % must define this per study; use a cell array of regular expressions. One cell per directory level.
end
if ~isfield(imPar,'tokenOrdering')
	imPar.tokenOrdering = []; % must match imPar.folderHierarchy: 1==subject, 2=visit, 3==session, 4==scan (if visit or session are omitted, they will be skipped)
else
	imPar.tokenOrdering = imPar.tokenOrdering(:)';
end
if ~isfield(imPar,'tokenScanAliases')
	imPar.tokenScanAliases = [];
else
	tokenScanAliasesOld = imPar.tokenScanAliases;
	imPar.tokenScanAliases = tokenScanAliasesOld(1:2:end);
	imPar.tokenScanAliases(:,2) = tokenScanAliasesOld(2:2:end);
end
if ~isfield(imPar,'tokenVisitAliases')
	imPar.tokenVisitAliases = [];
else
	tokenVisitAliasesOld = imPar.tokenVisitAliases;
	imPar.tokenVisitAliases = tokenVisitAliasesOld(1:2:end);
	imPar.tokenVisitAliases(:,2) = tokenVisitAliasesOld(2:2:end);
end
if ~isfield(imPar,'tokenSessionAliases')
	imPar.tokenSessionAliases = [];
else
	tokenSessionAliasesOld = imPar.tokenSessionAliases;
	imPar.tokenSessionAliases =tokenSessionAliasesOld(1:2:end);
	imPar.tokenSessionAliases(:,2) = tokenSessionAliasesOld(2:2:end);
end
if ~isfield(imPar,'bMatchDirectories')
	imPar.bMatchDirectories  = false;
end

% Specify the additional details of the conversion
if ~isfield(imPar,'bVerbose') || isempty(imPar.bVerbose)
    imPar.bVerbose = true;
end
if ~isfield(imPar,'bOverwrite') || isempty(imPar.bOverwrite)
    imPar.bOverwrite  = false; % NB, the summary file will be recreated anyway and dicom conversion in temp is always done, even if dest. exists
end
if ~isfield(imPar,'visitNames') || isempty(imPar.visitNames)
    imPar.visitNames = {};
end
if ~isfield(imPar,'nMaxVisits') || isempty(imPar.nMaxVisits)
    imPar.nMaxVisits = 0;
end
if ~isfield(imPar,'sessionNames') || isempty(imPar.sessionNames)
    imPar.sessionNames = {};
end
if ~isfield(imPar,'nMaxSessions') || isempty(imPar.nMaxSessions)
    imPar.nMaxSessions = 0;
end

%% 3. Run the DCM2NIIX
if bRunSubmodules(1)
	ExploreASL_ImportBIDS_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x);
end

%% 4. Run the NIIX to ASL-BIDS
if bRunSubmodules(2)
	% Loads the general configuration necessary for the conversion and BIDS saving
	bidsPar = xASL_bids_Config();

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load the study parameters + dataset description
	if ~exist(studyParPath,'file')
		warning('Study-par file is not provided.');
		studyPar = struct;
	else
		%studyParPath = xASL_init_ConvertM2JSON(studyParPath); % convert .m to .json - outdated import from an m-file
		studyPar = xASL_import_json(studyParPath);
	end
	
	% The Name has to be always assigned
	if ~isfield(studyPar,'Name')
		studyPar.Name = imPar.studyID;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Create the study description output and verify that all is there
	datasetDescription = ExploreASL_ImportBIDS_datasetDescription(studyPar,bidsPar.BIDSVersion,bidsPar.datasetDescription);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Make the output directory and save the description
		
	if ~exist(fullfile(imPar.BidsRoot),'dir')
		mkdir(fullfile(imPar.BidsRoot));
	end
	
	spm_jsonwrite(fullfile(imPar.BidsRoot,[bidsPar.datasetDescription.filename '.json']),datasetDescription);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Go through all subjects and check all the M0 and ASLs and modify the JSONs
	% This step should be completely automatic, just taking the info filled above and using it to convert to full BIDS.
	
	% Go through all subjects
	listSubjects = xASL_adm_GetFileList(imPar.AnalysisRoot,[],false,[],true);
	for iSubject = 1:length(listSubjects)
		
		subLabel = xASL_adm_CorrectName(listSubjects{iSubject},2);
		
		% Make a subject directory
		if ~exist(fullfile(imPar.BidsRoot,['sub-' subLabel]),'dir')
			mkdir(fullfile(imPar.BidsRoot,['sub-' subLabel]));
		end
	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Process all the anatomical files
		% Go throught the list of anat files
		for iAnatType = bidsPar.listAnatTypes
			
			% Check if it exists
			anatPath = '';
			if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},[iAnatType{1},'.nii']),'file')
				anatPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},iAnatType{1});
			end
			
			if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},[iAnatType{1} '_1'],[iAnatType{1},'.nii']),'file')
				anatPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},[iAnatType{1} '_1'],iAnatType{1});
			end
			
			% If anatomical file of this type exist, then BIDSify its structure
			if ~isempty(anatPath)
				
				% Create the anatomical directory
				if ~exist(fullfile(imPar.BidsRoot,['sub-' subLabel],'anat'),'dir')
					mkdir(fullfile(imPar.BidsRoot,['sub-' subLabel],'anat'));
				end
				
				% Copy the NiFTI file
				xASL_Copy([anatPath '.nii'],fullfile(imPar.BidsRoot,['sub-' subLabel],'anat',...
					['sub-' subLabel '_' iAnatType{1} '.nii.gz']));
				
				% Load the JSON
				jsonAnat = spm_jsonread([anatPath,'.json']);
				
				% Save the JSON
				jsonAnat = ExploreASL_ImportBIDS_jsonCheck(jsonAnat,bidsPar,0);
				spm_jsonwrite(fullfile(imPar.BidsRoot,['sub-' subLabel],'anat',['sub-' subLabel '_' iAnatType{1} '.json']),jsonAnat);
			end
		end
	
%% ERRORS FIXED	
% Dataset:		Philips_PCASL_3DGRASE_R5
% /sub-Sub1/fmap/sub-Sub1_dir-pa_m0scan.json:
%				Extra field: EffectiveEchoSpacing
%				Different value: TotalReadoutTime (0.029286 vs 0.013600)	

% Siemens_PASL_2DEPI_noBsup_AD - strange different value that breaks formatting

% Philips_PCASL_3DGRASE_functional
%			/sub-Patient1/perf/sub-Patient1_asl.json:
%				Different value: RepetitionTime (0.008077 vs 4.280000)
% %Siemens_PASL_2DEPI_noBsup_AD and the same for 2
% %xtra field: SiemensSliceTime
%% ERRORS TO FIX
%Philips_PCASL_2DEPI_GBM + strange error
%							Extra field: EffectiveEchoSpacing
%				Different value: SliceTiming (0.000000 vs 43.764700)

% Philips_PCASL_2DEPI_GBM	- NII
% Philips_PCASL_2DEPI_Ingenia_volunteer
% Philips_PCASL_2DEPI_Intera_volunteer
% Philips_PCASL_2DEPI_dummyLL
% Philips_PCASL_2DEPI_dummyQUASAR	
% Philips_PCASL_2DEPI_dummyMultiPLD
% Philips_PCASL_2DEPI_pharma - NII
% Philips_PCASL_2DEPI_pharma2 - NII - also M0
% Philips_PCASL_2DEPI_volunteer3 - NII
% Philips_PCASL_2DEPI_volunteer_1 - NII
% Philips_PCASL_2DEPI_volunteer_2 - NII
% Philips_PCASL_3DGRASE_functional - NII

% Philips_PCASL_2DEPI_volunteer3
%Extra field: EffectiveEchoSpacing
%				Different value: RepetitionTime (4.571680 vs 4.571681)
% Philips_PCASL_2DEPI_dummyQUASAR	
% Different value: VascularCrushingVenc (0.000000 vs 0.000000) + something other and strange
%Siemens_PASL_2DEPI_noBsup_AD and the same for 2
				%Different value: BolusCutOffDelayTime (0.800000 vs 1.400000)
%1.400000e+00				Different value: 				Different value: EffectiveEchoSpacing (0.000058 vs 0.000348)
%				Different value: SliceTiming (0.000000 vs 0.045000)
%Siemens_PASL_singleTI, Siemens_PCASL_2DEPI_AD
%Different value: BolusCutOffDelayTime (0.700000 vs 1.600000)
%Different value: EffectiveEchoSpacing (0.000128 vs 0.000510)
% Different value: SliceTiming (0.000000 vs 0.042500)
%Different value: TotalReadoutTime (0.032513 vs 0.032130)

	
	% For M0 in aslcontext, filling in PLD and labdur as zero

	% Test if all works equally well with dicominfo only
	
	% Multi TI - recognized as multiple TIs (and for similar fields as well) - give priority to the study-par file if the dicom value is single...
	
	% After final comparison, need to remove the random study descriptions and rename it correctly

	% Flip angle priority of studypar over data
	
	% Automatic reading of more ASL parameters
	
	% Automatic control/label order extraction
	
	% Q2TIPS has two timing entries
	
% 	if ~isfield(importStr{ii}.par,'ReadoutSegments') && isfield(importStr{ii}.x,'NumberSegments')
%		importStr{ii}.par.NumberSegments = importStr{ii}.x.NumberSegments;
%	end	
	
% Slice readouttiming for 3D can be set to 0

% M0 field can be numerical provided. Otherwise - true, separate, false, can be either copied as true or false. Or we simply assign true/false/file based on the fact if you provide M0 or if within series or if contol without BS
% Missing separate file+Bsup or missing separate+nocontrol should be evaluated as an error

% Update to new BIDS version

% Run the defacing module
% Do the anonymization

	% This can all go, will be in files
% 	% Labeling delays and durations
% 	if strcmpi(importStr{ii}.par.LabelingType,'PASL')
% 		%importStr{ii}.par.LabelingDuration = 0;% importStr{ii}.x.LabelingDuration           = 1.800;  % for PASL this is TI1
% 		importStr{ii}.par.PostLabelingDelay = importStr{ii}.x.InitialPostLabelDelay;
% 		if importStr{ii}.par.BolusCutOffFlag
% 			importStr{ii}.par.BolusCutOffDelayTime = importStr{ii}.par.BolusCutOffDelayTime + importStr{ii}.x.LabelingDuration;
% 		end
% 	else
% 		importStr{ii}.par.LabelingDuration = importStr{ii}.x.LabelingDuration;
% 		importStr{ii}.par.PostLabelingDelay = importStr{ii}.x.InitialPostLabelDelay;
% 	end

	
	
	fSes = xASL_adm_GetFileList(fullfile(imPar.AnalysisRoot,listSubjects{iSubject}),'^ASL.+$',false,[],true);
	
	% Go through all sessions
	for kk = 1:length(fSes)
		
		% Make a subject directory
		if length(fSes)>1
			sesLabel = ['ses-' fSes{kk}(5:end)];
			sesLabelUnd = ['_' sesLabel];
			if ~exist(fullfile(imPar.BidsRoot,['sub-' subLabel],sesLabel),'dir')
				mkdir(fullfile(imPar.BidsRoot,['sub-' subLabel],sesLabel));
				mkdir(fullfile(imPar.BidsRoot,['sub-' subLabel],sesLabel,'asl'));
			end
			inSesPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk});
			outSesPath = fullfile(imPar.BidsRoot,['sub-' subLabel],sesLabel);
		else
			sesLabel = '';
			sesLabelUnd = '';
			
			% Only one session - no session labeling
			if ~exist(fullfile(imPar.BidsRoot,['sub-' subLabel]),'dir')
				mkdir(fullfile(imPar.BidsRoot,['sub-' subLabel]));
				mkdir(fullfile(imPar.BidsRoot,['sub-' subLabel],'perf'));
			end
			inSesPath = fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk});
			outSesPath = fullfile(imPar.BidsRoot,['sub-' subLabel]);
		end
		
		% Check if there are multiple runs per session
		fRuns = xASL_adm_GetFileList(inSesPath,'^ASL4D_\d.nii+$',false,[],false);
		nSes = length(fRuns);
		
		for mm = 1:(max(nSes,1))
			if nSes
				aslLabel = ['ASL4D_' num2str(mm)];
				aslOutLabel = fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd '_run-' num2str(mm)]);
				aslOutLabelRelative = fullfile('perf',['sub-' subLabel sesLabelUnd '_run-' num2str(mm)]);
			else
				aslLabel = 'ASL4D';
				aslOutLabel = fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd]);
				aslOutLabelRelative = fullfile('perf',['sub-' subLabel sesLabelUnd]);
			end
			
			% Load the JSON
			jsonDicom = spm_jsonread(fullfile(inSesPath,[aslLabel '.json']));
			imNii = xASL_io_Nifti2Im(fullfile(inSesPath,[aslLabel '.nii']));
						
			if ~isempty(regexpi(jsonDicom.Manufacturer,'Philips'))
				scaleFactor = xASL_adm_GetPhilipsScaling(jsonDicom,xASL_io_ReadNifti(fullfile(inSesPath,[aslLabel '.nii'])));
			else
				scaleFactor = 0;
			end
			
			if scaleFactor
				imNii = imNii .* scaleFactor;
				xASL_io_SaveNifti(fullfile(inSesPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],imNii,[],1,[]);
			elseif size(imNii,4) == 1
				% The fourth dimension is 1, so we have to write the file again, to make sure the
				xASL_io_SaveNifti(fullfile(inSesPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],imNii,[],1,[]);
			else
				% Copy the ASL
				xASL_Copy(fullfile(inSesPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz']);
			end
			
			% Copy the basic ones
			jsonLocal = studyPar;
			
			% Copy all dicom ones
			for fn = fieldnames(jsonDicom)'
				if isfield(jsonLocal,fn{1}) && ~isequal(jsonLocal.(fn{1}),jsonDicom.(fn{1}))
					warning('User defined and DICOM field %s differ.\n',fn{1});
				end
				jsonLocal.(fn{1}) = jsonDicom.(fn{1});
			end
			
			% Check if BolusDuration field is present and not in conflict with the BolusCutoffDelayTime
			if isfield(jsonDicom,'BolusDuration')
				if ~isfield(studyPar,'BolusCutOffDelayTime')
					warning('Bolus duration obtained from DICOM, but not correctly redefined.');
				elseif ~isequal(jsonDicom.BolusDuration,studyPar.BolusCutOffDelayTime(1))
					warning('Bolus duration obtained from DICOM and the manuall defined one differ.');
				end
			end
			
			if isfield(jsonLocal,'GELabelingDuration')
				if isfield(studyPar,'LabelingDuration') && jsonLocal.GELabelingDuration ~= studyPar.LabelingDuration
					warning('Labeling duration mismatch with GE private field.');
				end
				if isfield(jsonDicom,'LabelingDuration') && jsonLocal.GELabelingDuration ~= studyPar.LabelingDuration
					warning('Labeling duration mismatch with GE private field.');
				end
				jsonLocal.LabelingDuration = jsonLocal.GELabelingDuration;
			end
			
			% BSup sanity check
			%if par.BackgroundSuppression = false
			%	then null and remove 'BackgroundSuppressionPulseTime and BackgroundSuppressionNumberPulses'
			%	verify that importStr{ii}.par.BackgroundSuppressionNumberPulses = length(importStr{ii}.par.BackgroundSuppressionPulseTime);
			
			% BSup sanity check - remove pulsenumbers and timings if BSup is OFF
			if jsonLocal.BackgroundSuppression == false
				if isfield(jsonLocal,'BackgroundSuppressionNumberPulses')
					jsonLocal = rmfield(jsonLocal,'BackgroundSuppressionNumberPulses');
				end
				if isfield(jsonLocal,'BackgroundSuppressionPulseTime')
					jsonLocal = rmfield(jsonLocal,'BackgroundSuppressionPulseTime');
				end
			else
				% If times are given, but not the number of pulses, then assign the length
				if isfield(jsonLocal,'BackgroundSuppressionPulseTime')
					if isfield(jsonLocal,'BackgroundSuppressionNumberPulses')
						if jsonLocal.BackgroundSuppressionNumberPulses ~= length(jsonLocal.BackgroundSuppressionPulseTime)
							fprintf('Warning: Number of pulses and their timings do not match.');
						end
					else
						jsonLocal.BackgroundSuppressionNumberPulses = length(jsonLocal.BackgroundSuppressionPulseTime);
					end
				end
			end
			
			% Free info about the sequence, now just the scanner type+software
			if isfield(jsonDicom,'ManufacturersModelName')
				jsonLocal.PulseSequenceDetails = jsonDicom.ManufacturersModelName;
				%else
				%	jsonLocal.PulseSequenceDetails = '';
			end
			if isfield(jsonDicom,'SoftwareVersions')
				if ~isempty(jsonLocal.PulseSequenceDetails)
					jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails '-'];
				end
				jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails jsonDicom.SoftwareVersions];
				%jsonLocal.MRSoftwareVersion = jsonDicom.SoftwareVersions;
			end
			
			% Process all the data and automatically fill in the missing parameters
			if strcmpi(jsonLocal.MRAcquisitionType,'2D')
				jsonLocal.PulseSequenceType = '2D_EPI';
			else
				if strcmpi(jsonLocal.Manufacturer,'GE') || strcmpi(jsonLocal.Manufacturer,'GE_WIP') || strcmpi(jsonLocal.Manufacturer,'GE_product')
					jsonLocal.PulseSequenceType = '3D_spiral';
				else
					jsonLocal.PulseSequenceType = '3D_GRASE';
				end
			end
			
			if isfield (studyPar,'SliceReadoutTime')
				if isfield(studyPar,'SliceTiming') && ~isequal(studyPar.SliceTiming,studyPar.SliceReadoutTime)
					sprintf('Warning, difference in SliceTiming and SliceReadoutTime');
				end
				studyPar.SliceTiming = studyPar.SliceReadoutTime;
			end
			
			% Fill in extra parameters based on the JSON from the data
			if jsonLocal.PulseSequenceType(1) == '2'
				% Take the studyPar as a prior source of SliceTiming since this is mostly wrong in DICOM otherwise
				if isfield(studyPar,'SliceTiming')
					jsonLocal.SliceTiming = studyPar.SliceTiming;
				end
				
				% The siemens field is rather reliable though
				if isfield(jsonLocal,'SiemensSliceTime') && ~isempty(jsonLocal.SiemensSliceTime)
					jsonLocal.SliceTiming = jsonLocal.SiemensSliceTime;
				end
				
				% If the length of SliceTiming fits to the number of slices, do nothing
				if length(jsonLocal.SliceTiming) ~= size(imNii,3)
					% if the length of studyPar.sliceTiming is higher than 1 and the difference non-zero then use this
					if length(jsonLocal.SliceTiming) > 1 && abs(jsonLocal.SliceTiming(2)-jsonLocal.SliceTiming(1)) > 0
						jsonLocal.SliceTiming = jsonLocal.SliceTiming(2)-jsonLocal.SliceTiming(1);
					end
					
					if abs(jsonLocal.SliceTiming) > 0
						jsonLocal.SliceTiming = ((0:(size(imNii,3)-1))')*jsonLocal.SliceTiming;
					end
				end
			else
				% 3D sequences should not have a SliceTiming or have it defined as zero
				if isfield(jsonLocal,'SliceTiming')
					jsonLocal = rmfield(jsonLocal,'SliceTiming');
				end
			end
			
			if isfield(jsonLocal,'EffectiveEchoSpacing')
				if jsonLocal.EffectiveEchoSpacing == 0
					jsonLocal = rmfield(jsonLocal,'EffectiveEchoSpacing');
				else
					jsonLocal.EffectiveEchoSpacing = abs(jsonLocal.EffectiveEchoSpacing);
				end
			end
			
			if isfield(jsonLocal,'TotalReadoutTime')
				if jsonLocal.TotalReadoutTime == 0
					jsonLocal = rmfield(jsonLocal,'TotalReadoutTime');
				else
					jsonLocal.TotalReadoutTime = abs(jsonLocal.TotalReadoutTime);
				end
			end
			
			%if isfield(studyPar,'RepetitionTime')
			%	jsonLocal.RepetitionTime = studyPar.RepetitionTime;
			%end
			
			% Check if TR is a vector - replace by the maximum then
			if length(jsonLocal.RepetitionTime) > 1
				jsonLocal.RepetitionTime = max(jsonLocal.RepetitionTime);
				warning('TR was a vector. Taking the maximum only.');
			end
			
			% Import the number of averages
			if isfield(jsonLocal,'NumberOfAverages') && (max(jsonLocal.NumberOfAverages) > 1)
				if isfield(studyPar,'TotalAcquiredVolumes')
					if max(jsonLocal.NumberOfAverages) ~= studyPar.TotalAcquiredVolumes
						warning('Discrepancy in the number of averages');
					end
				end
			end
			
			%m0 true - check asl context (copied as is from the input)
			%m0 false - mean control (or UseControlAsM0, no_background_suppression) (copied as is from the input)
			%m0 value - use value numerical (copied as is from the input)
			%m0 separate_file, or empty or anything - use that
			%m0 missing, empty or separate_scan - assign path, true or false based on the data
			
			% Type of an M0 image
			bJsonLocalM0isFile = 0;
			if ~isfield(studyPar,'M0') || isempty(studyPar.M0) || strcmpi(studyPar.M0,'separate_scan')
				if isfield(studyPar,'M0PositionInASL4D') && (max(studyPar.M0PositionInASL4D(:))>0)
					jsonLocal.M0 = true;
				elseif xASL_exist(fullfile(inSesPath,'M0.nii'))
					if length(fSes)>1
						jsonLocal.M0 = fullfile('perf',['sub-' subLabel sesLabelUnd]);
						bJsonLocalM0isFile = 1;
					else
						jsonLocal.M0 = fullfile('perf',['sub-' subLabel sesLabelUnd]);
						bJsonLocalM0isFile = 1;
					end
				else
					if ~isempty(strfind(studyPar.ASLContext,bidsPar.strM0scan))
						jsonLocal.M0 = true;
					else
						jsonLocal.M0 = false;
					end
				end
			else
				if strcmpi(studyPar.M0,'UseControlAsM0')
					jsonLocal.M0 = false;
				else
					if strcmpi(studyPar.M0,'no_background_suppression')
						jsonLocal.M0 = false;
					else
						jsonLocal.M0 = studyPar.M0;
					end
				end
			end
			
			% If Post-labeling delay or labeling duration is longer than 1, but shorten then number of volumes
			% then repeat it
			listFieldsRepeat = {'PostLabelingDelay', 'LabelingDuration','VascularCrushingVenc','FlipAngle','RepetitionTime'};
			for iRepeat = 1:length(listFieldsRepeat)
				if isfield(jsonLocal,(listFieldsRepeat{iRepeat})) && (length(jsonLocal.(listFieldsRepeat{iRepeat})) > 1) && (size(imNii,4) ~= length(jsonLocal.(listFieldsRepeat{iRepeat})))
					if mod(size(imNii,4),length(jsonLocal.(listFieldsRepeat{iRepeat})))
						error('Cannot find a match between the %s and the 4th dimension of the NIFTI.\n',listFieldsRepeat{iRepeat});
					else
						jsonLocal.(listFieldsRepeat{iRepeat}) = repmat(jsonLocal.(listFieldsRepeat{iRepeat}),[1 size(imNii,4)/length(jsonLocal.(listFieldsRepeat{iRepeat}))]);
					end
				end
			end
			
			if isfield(studyPar,'RepetitionTime') && studyPar.RepetitionTime ~= jsonLocal.RepetitionTime
				warning('User defined repetition time differs from DICOM, using user defined: %d vs %d \n',studyPar.RepetitionTime,jsonLocal.RepetitionTime);
			end
			
			% Reformat ASLcontext field
			% Remove ',' and ';' at the 
			if (jsonLocal.ASLContext(end) == ';') || (jsonLocal.ASLContext(end) == ',')
				jsonLocal.ASLContext = jsonLocal.ASLContext(1:(end-1));
			end
			
			% Replace all ',' and ';' by \n
			jsonLocal.ASLContext = strrep(jsonLocal.ASLContext, ' ','');
			jsonLocal.ASLContext = strrep(jsonLocal.ASLContext, ';',',');
			lengthASLContext = sum(jsonLocal.ASLContext == ',')+1;
			jsonLocal.ASLContext = strrep(jsonLocal.ASLContext, ',',sprintf('\n'));
			
			
			% Check if the length is the same
			if size(imNii,4) ~= lengthASLContext
				% Check if we can simply repeat it 
				if mod(size(imNii,4),lengthASLContext)
					error('Cannot find a match between the ASLContext and the 4th dimension of the NIFTI');
				else
					numRepeat = size(imNii,4)/lengthASLContext;
					tmpStr = jsonLocal.ASLContext;
					for iRepeat = 2:numRepeat
						jsonLocal.ASLContext = sprintf('%s\n%s',jsonLocal.ASLContext,tmpStr);
					end
				end
			end
			jsonLocal.ASLContext = sprintf('%s\n',jsonLocal.ASLContext);
			
			% Remove the AslContext field and save it as a separate file
			fContext = fopen([aslOutLabel '_' bidsPar.strAslContext '.tsv'],'w+');
			fwrite(fContext,sprintf('volume_type\n'));
			fwrite(fContext,jsonLocal.ASLContext);
			fclose(fContext);
			
			jsonLocal = rmfield(jsonLocal,'ASLContext');
			
			if mm == 1
				for nn = 1:2
					if nn == 1
						nnStrIn = '';
						if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk},'M0PERev.nii'))
							nnStrOut = '_dir-ap';
							
							tagPhaseEncodingDirection = 'j-';
							jsonLocal.PhaseEncodingDirection = 'j-';
							tagIntendedFor = [];
							tagTotalReadoutTime = studyPar.TotalReadoutTime;
							
							if bJsonLocalM0isFile
								jsonLocal.M0 = [jsonLocal.M0 nnStrOut '_' bidsPar.strM0scan '.nii.gz'];
							end
						else
							if bJsonLocalM0isFile
								jsonLocal.M0 = [jsonLocal.M0 '_' bidsPar.strM0scan '.nii.gz'];
							end
							nnStrOut = '';
							tagPhaseEncodingDirection = [];
							tagIntendedFor = [];
							tagTotalReadoutTime = [];
						end
					else
						nnStrIn = 'PERev';
						nnStrOut = '_dir-pa';
						tagPhaseEncodingDirection = 'j';
						tagIntendedFor = fullfile('perf',['sub-' subLabel sesLabelUnd '_dir-ap' '_' bidsPar.strM0scan '.nii.gz']);
						
						if isfield(studyPar,'TotalReadoutTime')
							tagTotalReadoutTime = studyPar.TotalReadoutTime;
						else
							tagTotalReadoutTime = [];
						end
					end
					% If M0, then copy M0 and add ASL path to the IntendedFor
					if xASL_exist(fullfile(imPar.AnalysisRoot,listSubjects{iSubject},fSes{kk},['M0' nnStrIn '.nii']))
						jsonM0 = spm_jsonread(fullfile(inSesPath,['M0' nnStrIn '.json']));
						imM0   = xASL_io_Nifti2Im(fullfile(inSesPath,['M0' nnStrIn '.json']));
						
								
						if ~isempty(regexpi(jsonDicom.Manufacturer,'Philips'))
							scaleFactor = xASL_adm_GetPhilipsScaling(jsonDicom,xASL_io_ReadNifti(fullfile(inSesPath,['M0' nnStrIn '.nii'])));
						else
							scaleFactor = 0;
						end
						
						if scaleFactor
							imM0 = imM0 .* scaleFactor;
						end
						
						jsonM0Write = jsonM0;
						
						if isfield(jsonLocal,'SliceTiming')
							% Issue a warning if the SliceTiming was already existing for M0, but still overwrite with ASL one
							if isfield(jsonM0Write,'SliceTiming')
								warning('SliceTiming already existed for M0, overwriting with ASL');
							end
							
							if size(imNii,3) == size(imM0,3)
								% Either copy if the save number of slices in M0 as in ASL
								jsonM0Write.SliceTiming = jsonLocal.SliceTiming;
							else
								% Or recalculate for M0 if the number of slices differ
								jsonM0Write.SliceTiming = ((0:(size(imM0,3)-1))')*(jsonLocal.SliceTiming(2)-jsonLocal.SliceTiming(1));
							end
						else
							if isfield(jsonM0Write,'SliceTiming')
								jsonM0Write = rmfield(jsonM0Write,'SliceTiming');
								warning('Removing pre-existing SliceTiming from M0, as there was no SliceTiming for ASL');
							end
						end
						
						if isfield(studyPar,'RepetitionTime')
							jsonM0Write.RepetitionTime = studyPar.RepetitionTime;
						else
							jsonM0Write.RepetitionTime = jsonM0.RepetitionTime;
						end
						
						jsonM0Write.IntendedFor = [aslOutLabelRelative '_asl.nii.gz'];
						
						if ~isempty(tagPhaseEncodingDirection)
							jsonM0Write.PhaseEncodingDirection = tagPhaseEncodingDirection;
						end
						
						if isfield(jsonM0Write,'EffectiveEchoSpacing')
							if jsonM0Write.EffectiveEchoSpacing == 0
								jsonM0Write = rmfield(jsonM0Write,'EffectiveEchoSpacing');
							else
								jsonM0Write.EffectiveEchoSpacing = abs(jsonM0Write.EffectiveEchoSpacing);
							end
						end
						
						if isfield(jsonM0Write,'TotalReadoutTime')
							if jsonM0Write.TotalReadoutTime == 0
								jsonM0Write = rmfield(jsonM0Write,'TotalReadoutTime');
							else
								jsonM0Write.TotalReadoutTime = abs(jsonM0Write.TotalReadoutTime);
							end
						end
						
						if ~isempty(tagIntendedFor)
							jsonM0Write.IntendedFor = tagIntendedFor;
						end
						
						if ~isempty(tagTotalReadoutTime)
							jsonM0Write.TotalReadoutTime = tagTotalReadoutTime;
						end
						
						if nn == 2 && ~exist(fullfile(outSesPath,'fmap'),'dir')
							mkdir(fullfile(outSesPath,'fmap'));
						end
						
						% if scaling modified then save instead of copy
						if scaleFactor || size(imM0,4) == 1
							if nn == 1
								xASL_io_SaveNifti(fullfile(inSesPath,['M0' nnStrIn '.nii']),fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
							else
								xASL_io_SaveNifti(fullfile(inSesPath,['M0' nnStrIn '.nii']),fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
							end
						else
							% Copy the M0
							if nn == 1
								xASL_Copy(fullfile(inSesPath,['M0' nnStrIn '.nii']),...
									fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']));
							else
								xASL_Copy(fullfile(inSesPath,['M0' nnStrIn '.nii']),...
									fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']));
							end
						end
						% Save JSON to new dir
						jsonM0Write = ExploreASL_ImportBIDS_jsonCheck(jsonM0Write,bidsPar,1);
						if nn == 1
							spm_jsonwrite(fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
						else
							spm_jsonwrite(fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
						end
					end
				end
			else
				if bJsonLocalM0isFile
					jsonLocal.M0 = [jsonLocal.M0 '.nii.gz'];
				end
			end
			% Save JSON to new dir
			jsonLocal = ExploreASL_ImportBIDS_jsonCheck(jsonLocal,bidsPar,1);
			spm_jsonwrite([aslOutLabel '_asl.json'],jsonLocal);
			
		end
	end
end









%%%%%%%%%%%% end todo
end


%% 5. Run defacing
if bRunSubmodules(3)
end

% end of import
end

% -----------------------------------------------------------------------------
%
% -----------------------------------------------------------------------------
function s = AppendNiftiParameters(nii_files)
% This function outputs s=[FileName voxel size XYZ matrix size XYZ]
s = [];

if ischar(nii_files)
    nii_files = {nii_files};
end

for iNii=1:length(nii_files)
    [~, Ffile, Fext] = fileparts(nii_files{iNii});
    s = sprintf(',"%s"', [Ffile Fext]); % filename

    tempnii = xASL_io_ReadNifti(nii_files{iNii});
    s = [s sprintf(',%g', tempnii.hdr.pixdim(2:5) )]; % voxel size XYZ
    s = [s sprintf(',%g', tempnii.hdr.dim(2:5) )]; % matrix size XYZ
end
end

% -----------------------------------------------------------------------------
%
% -----------------------------------------------------------------------------
function [s, FieldNames] = AppendParmsParameters(parms)
% This function outputs s=fields of _parms.mat
s = [];

FieldNames = {'RepetitionTime' 'EchoTime' 'NumberOfAverages' 'RescaleSlope' 'RescaleSlopeOriginal'...
    'MRScaleSlope' 'RescaleIntercept' 'AcquisitionTime' 'AcquisitionMatrix' 'TotalReadoutTime'...
    'EffectiveEchoSpacing'};

if ~isempty(parms)
    for iField=1:length(FieldNames)
        if isfield(parms,FieldNames{iField})
            s = [s ',' xASL_num2str(parms.(FieldNames{iField}))];
        else
            s = [s ',n/a'];
        end
    end
end
    
end

% -----------------------------------------------------------------------------
%
% -----------------------------------------------------------------------------
function [dcm2niiCatchedErrors] = CatchErrors(WarningID, WarningMessage, WarningLine, WarningFileName, WarningPath, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, StackIn)
% Catch reported warnings/errors, print them if verbose, & add them to a structure of warnings/errors to be stored for later QC

    if imPar.bVerbose % print warning if we want verbose
        warning(WarningMessage);
    end

    % Find index of the warning to store
    if isempty(fields(dcm2niiCatchedErrors))
        IndexN = 1;
    else
        IndexN = length(dcm2niiCatchedErrors)+1;
    end

    % store the warning/error
    dcm2niiCatchedErrors(IndexN).scan_name = scan_name;
    dcm2niiCatchedErrors(IndexN).scanpath = scanpath;
    dcm2niiCatchedErrors(IndexN).destdir = destdir;
    dcm2niiCatchedErrors(IndexN).identifier = WarningID;
    dcm2niiCatchedErrors(IndexN).message = WarningMessage;
    dcm2niiCatchedErrors(IndexN).cause = 'n/a';

    if exist('StackIn', 'var')
        dcm2niiCatchedErrors(IndexN).stack = StackIn;
    else
        dcm2niiCatchedErrors(IndexN).stack.file = fullfile(WarningPath, [WarningFileName '.m']);
        dcm2niiCatchedErrors(IndexN).stack.name = WarningFileName;
        dcm2niiCatchedErrors(IndexN).stack.line = WarningLine(end).line;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the import dcm2nii part of the import
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExploreASL_ImportBIDS_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize defaults of dcm2nii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dcm2niiCatchedErrors = struct; % initialization
if bCheckPermissions
	dcm2niiDir = fullfile(x.MyPath, 'External', 'MRIcron');
	xASL_adm_CheckPermissions(dcm2niiDir, true); % dcm2nii needs to be executable
end
if ~isfield(imPar,'dcm2nii_version') || isempty(imPar.dcm2nii_version)
    imPar.dcm2nii_version = '20190902'; % OR for PARREC imPar.dcm2nii_version = '20101105'; THIS IS AUTOMATED BELOW
end
if ~isfield(imPar,'dcmExtFilter') || isempty(imPar.dcmExtFilter)
    imPar.dcmExtFilter = '^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$'; % the last one is because some convertors save files without extension, but there would be a dot/period before a bunch of numbers
end
	
	if ~isfield(imPar,'SkipSubjectIfExists') || isempty(imPar.SkipSubjectIfExists)
		% allows to skip existing subject folders in the analysis folder, when this is set to true,
		% avoiding partly re-importing/converting dcm2niiX when processing has been partly done
		imPar.SkipSubjectIfExists = false;
	else
		warning('Skipping existing subjects in analysis folder');
		fprintf('If you want to overwrite, first remove the full subject folder');
	end
	
	% Create the basic folder structure for raw & derivative data
	if ~exist(imPar.RawRoot, 'dir')
		warning(['Couldnt find ' imPar.RawRoot ', trying to find a different folder instead...']);
		
		% find any folder except for analysis, source, raw, derivatives
		% xASL_adm_GetFileList uses regular expressions, to create a nice list of foldernames,
		% with/without FullPath (FPList), with/without recursive (FPListRec)
		% very powerful once you know how these work
		FolderNames = xASL_adm_GetFileList(fullfile(imPar.RawRoot, imPar.studyID), '^(?!(analysis|derivatives|source|raw)).*$', 'FPList', [0 Inf], true);
		
		if length(FolderNames)==1
			imPar.RawRoot = FolderNames{1};
			fprintf('%s\n', ['Found ' imPar.RawRoot ' as raw folder']);
		else
			error('Couldnt find a raw folder, please rename one, or move other folders');
		end
	end
	
	xASL_adm_CreateDir(imPar.AnalysisRoot);
	
	if bCheckPermissions
		xASL_adm_CheckPermissions(imPar.RawRoot, false); % don"t need execution permisions
		xASL_adm_CheckPermissions(imPar.AnalysisRoot, false);  % don"t need execution permisions
	end
	
	% here we try to fix backwards compatibility, but this may break
	if length(imPar.tokenOrdering)==3 % backwards compatibility Visits
		imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 imPar.tokenOrdering(2:3)]; % insert Visits(none)
	elseif length(imPar.tokenOrdering)==2 % backwards compatibility Visits & sessions
		imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 0 imPar.tokenOrdering(2)]; % insert sessions & visits(none)
	end
	
	% Path to the dictionary to initialize - we need to keep track if the dictionary has been set, because
	% Dicominfo can be used despite bUSEDCMTK==1 when DCMTK fails
	pathDcmDict = fullfile(x.MyPath,'External','xASL_DICOMLibrary.txt');
	if ~bUseDCMTK
		% -----------------------------------------------------------------------------
		% Initialize dicom dictionary by appending private philips stuff to a temporary copy
		% -----------------------------------------------------------------------------
		dicomdict('set', pathDcmDict);
	end
	
	fid_summary = -1; % initialize to be able to catch errors and close if valid
	
	% change dcmnii_version for PARREC if needed
	if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
		imPar.dcm2nii_version = '20101105';
	end
	
	% redirect output to a log file
	diary_filepath = fullfile(imPar.AnalysisRoot, ['import_log_' imPar.studyID '_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
	diary(diary_filepath);
	
	if imPar.bMatchDirectories
		strLookFor = 'Directories';
	else
		strLookFor = 'Files';
	end
	
	%
	% -----------------------------------------------------------------------------
	% Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
	% -----------------------------------------------------------------------------
	
	% Recursively scan the directory tree using regular exspressions at each directory level. Use ()-brackets to extract tokens
	% that will be used to identify subjects, sessions and scans. In the loop below it is possible to translate the tokens
	% to more convenient strings before using them in destination paths.
	[matches, tokens] = xASL_adm_FindByRegExp(imPar.RawRoot, imPar.folderHierarchy, 'StripRoot', true, 'Match', strLookFor,'IgnoreCase',true);
	if isempty(matches)
		warning('No matching files, skipping');
		return;
	elseif imPar.bVerbose
		fprintf('\nMatching files:\n');
		disp(matches);
		fprintf('#=%g\n',length(matches));
	end
	
	% Copy the columns into named vectors. This construction allows for arbitrary directory hierarchies.
	% Make sure to select the columns that correspond to the folder ordering defined using the regular expressions above.
	% Define Subjects
	
	% SUBJECTS
	vSubjectIDs = tokens(:,imPar.tokenOrdering(1)); % cell vector with extracted subject IDs (for all sessions and scans)
	
	% VISITS
	if imPar.tokenOrdering(2)==0
		% a zero means: no visits applicable
		bUseVisits = false;
		vVisitIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject has a single visit
		imPar.tokenVisitAliases = {'^1$', '_1'};
	else
		bUseVisits = true;
		vVisitIDs = tokens(:,imPar.tokenOrdering(2)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	% SESSIONS
	if imPar.tokenOrdering(3)==0
		% a zero means: no sessions applicable
		bUseSessions = false;
		vSessionIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject-visit has a single session
		imPar.tokenSessionAliases = {'^1$', 'ASL_1'};
	else
		bUseSessions = true;
		vSessionIDs = tokens(:,imPar.tokenOrdering(3)); % cell vector with extracted session IDs (for all subjects and scans)
	end
	
	% SCANTYPES
	vScanIDs = tokens(:,imPar.tokenOrdering(4)); % cell vector with extracted scan IDs (for all subjects and sessions)
	
	% Convert the vectors to unique & sort sets by the output aliases
	subjectIDs  = sort(unique(vSubjectIDs));
	nSubjects = length(subjectIDs);
	
	visitIDs  = unique(vVisitIDs);
	% sort by output
	if length(visitIDs)>1
		for iV=1:length(visitIDs)
			IDrow(iV) = find(cellfun(@(y) strcmp(y,visitIDs{iV}), imPar.tokenVisitAliases(:,1)));
		end
		visitIDs = visitIDs(IDrow);
	end
	nVisits = length(visitIDs);
	
	sessionIDs  = sort(unique(vSessionIDs));
	nSessions = length(sessionIDs);
	
	scanIDs = sort(unique(lower(vScanIDs)));
	nScans = length(scanIDs);
	
	% VISIT NAMES
	if isempty(imPar.visitNames)
		if isempty(visitIDs)
			imPar.visitNames = cell(nVisits,1);
			for kk=1:nVisits
				imPar.visitNames{kk}=sprintf('ASL_%g',kk);
			end
		else
			imPar.visitNames = visitIDs;
		end
	end
	
	% SESSION NAMES
	% optionaly we can have human readble session names; by default they are the same as the original tokens in the path
	if isempty(imPar.sessionNames)
		if isempty(sessionIDs)
			imPar.sessionNames = cell(nSessions,1);
			for kk=1:nSessions
				imPar.sessionNames{kk}=sprintf('ASL_%g',kk);
			end
		else
			imPar.sessionNames = sessionIDs;
		end
	end
	
	% SCAN NAMES
	scanNames = scanIDs;
	
	% sanity check for missing elements
	if nSubjects==0
		error('No subjects')
	end
	if nVisits==0
		error('No visits')
	end
	if nSessions==0
		error('No sessions')
	end
	if nScans==0
		error('No scans')
	end
	
	% preallocate space for (global) counts
	converted_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	skipped_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	missing_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
	
	% define a cell array for storing info for parameter summary file
	xASL_adm_CreateDir(imPar.AnalysisRoot);
	summary_lines = cell(nSubjects,nVisits,nSessions,nScans);
	
	% -----------------------------------------------------------------------------
	% import subject by subject, visit by visit, session by session, scan by scan
	% -----------------------------------------------------------------------------
	fprintf('%s\n', 'Running import (i.e. dcm2niiX)');
	separatorline = repmat(char('+'),1,80);
	for iSubject=1:nSubjects
		subjectID = subjectIDs{iSubject};
		
		for iVisit=1:nVisits
			visitID = visitIDs{iVisit};
			
			% convert visit ID to a suitable name
			if size(imPar.tokenVisitAliases,2)==2
				iAlias = find(~cellfun(@isempty,regexp(visitIDs{iVisit},imPar.tokenVisitAliases(:,1),'once')));
				if ~isempty(iAlias)
					imPar.visitNames{iVisit} = imPar.tokenVisitAliases{iAlias,2};
				end
			end
			
			if bUseVisits % only pad VisitID _1 _2 _3 etc if there are visits specified
				% Multiple visits is defined by the tokenVisitAliases.
				% If this is non-existing, it is set to 1, and if it does exist,
				% it will put the _1 _2 _3 etc in the folder
				% this fix allows to import a single visit from a range of
				% specified visits
				SubjDir = fullfile(imPar.AnalysisRoot, [subjectID imPar.visitNames{iVisit}]);
				% if strcmp(imPar.visitNames{iVisit},'_1') % only pad the visitID _1 _2 _3 etc if there are multiple visits
			else
				SubjDir = fullfile(imPar.AnalysisRoot, subjectID);
			end
			
			if imPar.SkipSubjectIfExists && exist(SubjDir, 'dir')
				continue; % we found the subject dir (i.e. SubjectVisit), so we skip it
				% this is ignored when imPar.SkipSubjectIfExists is set to
				% false (default)
			end
			
			fprintf('%s\nImporting subject=%s:   \n',separatorline,[subjectID imPar.visitNames{iVisit}]); % display subject-visit ID
			
			% loop through all sessions
			for iSession=1:nSessions
				sessionID = sessionIDs{iSession};
				
				% convert session ID to a suitable name
				if size(imPar.tokenSessionAliases,2)==2
					iAlias = find(~cellfun(@isempty,regexp(sessionID,imPar.tokenSessionAliases(:,1),'once')));
					if ~isempty(iAlias)
						imPar.sessionNames{iSession} = imPar.tokenSessionAliases{iAlias,2};
					end
				end
				
				for iScan=1:nScans
					scanID = scanIDs{iScan};
					summary_line = [];
					first_match = [];
					summary_lines{iSubject,iVisit,iSession,iScan} = 'n/a';
					
					if ~imPar.bVerbose % if not verbose, track % progress
						CounterN = (iSession-1)*nScans+iScan;
						CounterT = nSessions*nScans;
						xASL_TrackProgress(CounterN, CounterT);
					end
					
					% convert scan ID to a suitable name and set scan-specific parameters
					if size(imPar.tokenScanAliases,2)==2
						iAlias = find(~cellfun(@isempty,regexpi(scanID,imPar.tokenScanAliases(:,1),'once')));
						if ~isempty(iAlias)
							scanNames{iScan} = imPar.tokenScanAliases{iAlias,2};
						else
							% keep the original name
							WarningMessage = ['ExploreASL_Import: Unknown scan ID ' scanID ' found, don"t know what this is'];
							dcm2niiCatchedErrors = CatchErrors('isempty(iAlias)', WarningMessage, dbstack, mfilename, pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
						end
					end
					scan_name = scanNames{iScan};
					
					% minimalistic feedback of where we are
					if imPar.bVerbose; fprintf('>>> Subject=%s, visit=%s, session=%s, scan=%s\n',subjectID, visitID, num2str(iSession), scan_name); end
					
					bOneScanIsEnough = false; % default
					bPutInSessionFolder = true; % by default put in session folder
					switch scan_name
						case {'ASL4D', 'M0', 'ASL4D_RevPE', 'func_bold'}
							bPutInSessionFolder = true;
						case {'T1', 'WMH_SEGM', 'FLAIR'}
							bPutInSessionFolder = false;
					end
					
					if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
						imPar.dcm2nii_version = '20101105';
					end
					
					% now pick the matching one from the folder list
					iMatch = find(strcmp(vSubjectIDs,subjectID) & strcmp(vVisitIDs, xASL_adm_CorrectName(visitID,2,'_')) & strcmp(vSessionIDs,sessionID) & strcmpi(vScanIDs,scanID) ); % only get the matching session
					if isempty(iMatch)
						% only report as missing if we need a scan for each session (i.e. ASL)
						if sum(converted_scans(iSubject,iVisit,:,iScan))==0
							WarningMessage = ['Missing scan: ' [subjectID imPar.visitNames{iVisit}] ', ' num2str(iSession) ', ' scan_name];
							if imPar.bVerbose; warning(WarningMessage); end
							missing_scans(iSubject, iVisit, iSession, iScan) = 1;
						end
						
						summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
						continue;
						warning('Dont forget to comment continue here for debugging');
					end
					
					% determine input and output paths
					bSkipThisOne = false;
					branch = matches{iMatch};
					scanpath = fullfile(imPar.RawRoot,branch);
					
					if ~isempty(strfind(scanNames{iScan}, 'ASL4D')) || ~isempty(strfind(scanNames{iScan}, 'M0'))
						session_name = ['ASL_' num2str(iSession)];
					elseif ~isempty(strfind(scanNames{iScan}, 'DSC4D'))
						session_name = ['DSC_' num2str(iSession)];
					else
						session_name = [scanNames{iScan} '_' num2str(iSession)]; % Allow multiple ScanTypes for sessions
					end
					
					if bPutInSessionFolder
						destdir = fullfile(SubjDir, session_name);
					else % put in subject folder instead of session folder
						destdir = SubjDir;
					end
					
					if bOneScanIsEnough && sum(converted_scans(iSubject,iVisit,:,iScan))~=0
						% one scan is enough, so skip this one if there was already a scan converted of this type (i.e. T1)
						if imPar.bVerbose; fprintf('Skipping scan: %s, %s, %s\n',[subjectID imPar.visitNames{iVisit}],session_name,scan_name); end
						bSkipThisOne = true;
						destdir = []; % just in case
					end
					
					% start the conversion if this scan should not be skipped
					if bSkipThisOne
						summary_line = sprintf(',"skipped",,,,,,,,');
						skipped_scans(iSubject, iVisit, iSession, iScan) = 1;
					else
						nii_files = {};
						xASL_adm_CreateDir(destdir);
						
						% check if we have a nii(gz) file, or something that needs to be converted (parrec/dicom)
						if ~exist(scanpath, 'dir') && ~isempty(regexpi(scanpath,'(\.nii|\.nii\.gz)$'))
							% we found a NIfTI file
							% check if output exists
							first_match = fullfile(destdir, [scan_name '.nii']);
							if imPar.bOverwrite || ~xASL_exist(first_match,'file')
								[~, fname, fext] = fileparts(scanpath);
								destfile = fullfile(destdir, [fname fext]); % will be renamed later
								xASL_Copy(scanpath, destfile, imPar.bOverwrite, imPar.bVerbose);
								% gunzip if required
								destfile = xASL_adm_UnzipNifti(destfile);
								xASL_Move(destfile, first_match, imPar.bOverwrite, imPar.bVerbose);
							end
							nii_files{1} = first_match;
						else % we found dicom files
							% -----------------------------------------------------------------------------
							% start the conversion. Note that the dicom filter is only in effect when a directory is specified as input.
							% -----------------------------------------------------------------------------
							try
								[nii_files, scan_name, first_match, MsgDcm2nii] = xASL_io_dcm2nii(scanpath, destdir, scan_name, 'DicomFilter', imPar.dcmExtFilter, 'Verbose', imPar.bVerbose, 'Overwrite', imPar.bOverwrite, 'Version', imPar.dcm2nii_version, 'x', x);
								
								% If dcm2nii produced a warning or error, catch this & store it
								if ~isempty(MsgDcm2nii) && ~isempty(regexpi(MsgDcm2nii,'.*(warning|error).*')) % if it contains a warning/error
									dcm2niiCatchedErrors = CatchErrors('xASL_io_dcm2nii', MsgDcm2nii, dbstack, ['dcm2nii_' imPar.dcm2nii_version], pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
								end
								
							catch ME
								dcm2niiCatchedErrors = CatchErrors(ME.identifier, ME.message, [], [], [], scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, ME.stack);
								
								if imPar.bVerbose; warning(['dcm2nii ' scanpath ' crashed, skipping']); end
								if imPar.bVerbose; warning('Check whether the scan is complete'); end
								first_match = xASL_adm_GetFileList(scanpath, ['.*' imPar.dcmExtFilter],'FPList',[0 Inf]);
								if  ~isempty(first_match); first_match = first_match{1}; end
							end
						end
						
						% Merge NIfTIs if there are multiples
						if length(nii_files)>1
							1;
						end
						% For ASL or M0, merge multiple files
						if length(nii_files)>1
							if ~isempty(strfind(scan_name,'ASL4D'))
								nii_files = xASL_bids_MergeNifti(nii_files,'ASL');
							elseif  ~isempty(strfind(scan_name,'M0'))
								nii_files = xASL_bids_MergeNifti(nii_files,'M0');
							end
						end
						
						% Extract relevant parameters from nifti header and append to summary file
						summary_line = AppendNiftiParameters(nii_files);
						converted_scans(iSubject,iSession,iScan) = 1;
					end
					
					% extract relevant parameters from dicom header, if not
					% already exists
					% Find JSONpath that is there already
					SavePathJSON = {};
					SavePathJSON{1} = fullfile(destdir, [scan_name '.json']);
					SavePathJSON{2} = fullfile(destdir, [session_name '.json']);
					for iPath=1:length(nii_files)
						% now we add the path only if it didnt exist already in this list
						tmpNewPath = [nii_files{iPath}(1:end-4) '.json'];
						if ~max(cellfun(@(y) strcmp(y, tmpNewPath), SavePathJSON))
							SavePathJSON{end+1} = tmpNewPath;
						end
					end
					
					for iPath=1:length(SavePathJSON)
						if exist(SavePathJSON{iPath}, 'file') && ~isempty(first_match)
							[~, ~, fext] = fileparts(first_match);
							if  strcmpi(fext,'.PAR')
								parms = xASL_adm_Par2JSON(first_match, SavePathJSON{iPath});
							elseif strcmpi(fext,'.nii')
								parms = [];
							elseif imPar.bMatchDirectories
								Fpath  = fileparts(first_match);
								[parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, Fpath, SavePathJSON{iPath}, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
								clear Fpath Ffile Fext
							else
								[parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, first_match, SavePathJSON{iPath}, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
							end
						end
					end
					
					% correct nifti rescale slope if parms.RescaleSlopeOriginal =~1
					% but nii.dat.scl_slope==1 (this can happen in case of
					% hidden scale slopes in private Philips header,
					% that is dealt with by xASL_bids_Dicom2JSON but not by
					% dcm2niiX
					
					if ~isempty(nii_files) && exist('parms','var')
						[TempLine, PrintDICOMFields] = AppendParmsParameters(parms);
						summary_line = [summary_line TempLine];
					end
					
					if bClone2Source % make a copy of analysisdir in sourcedir
						if ~isempty(nii_files)
							for iFile=1:length(nii_files)
								% replace 'analysis' by 'source'
								[iStart, iEnd] = regexp(nii_files{iFile}, 'analysis');
								DestPath = [nii_files{iFile}(1:iStart-1) 'source' nii_files{iFile}(iEnd+1:end)];
								xASL_Copy(nii_files{iFile}, DestPath, true);
								% do the same for other extensions
								Extensions = {'.json' '_parms.json'};
								for iExt=1:length(Extensions)
									[Fpath, Ffile] = xASL_fileparts(nii_files{iFile});
									CopyPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
									[Fpath, Ffile] = xASL_fileparts(DestPath);
									DestPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
									if xASL_exist(CopyPath)
										xASL_Copy(CopyPath, DestPath, true);
									end
								end
							end
						end
					end
					
					% Copy single dicom as QC placeholder
					if bCopySingleDicoms && ~isempty(first_match)
						xASL_Copy(first_match, fullfile(destdir, ['DummyDicom_' scan_name '.dcm']), imPar.bOverwrite, imPar.bVerbose);
					end
					
					% store the summary info so it can be sorted and printed below
					summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
				end % scansIDs
			end % sessionIDs
		end % visitIDs
	end % subjectIDs
	
	
	
	
	% -----------------------------------------------------------------------------
	% create summary file
	% -----------------------------------------------------------------------------
	summary_filepath = fullfile(imPar.AnalysisRoot, 'import_summary.csv');
	fid_summary = fopen(summary_filepath,'wt');
	% Print headers for parameters obtained from NIfTI file
	fprintf(fid_summary,'subject,visit,session,scan,filename,dx,dy,dz,dt,nx,ny,nz,nt');
	% Print headers for parameters obtained from DICOM file
	if exist('PrintDICOMFields','var')
		for iField=1:length(PrintDICOMFields)
			fprintf(fid_summary,[',' PrintDICOMFields{iField}]);
		end
	end
	fprintf(fid_summary,'\n');
	
	for iScan=1:nScans
		for iSubject=1:nSubjects
			for iVisit=1:nVisits
				for iSession=1:nSessions
					if converted_scans(iSubject, iVisit, iSession, iScan) || skipped_scans(iSubject, iVisit, iSession, iScan) || missing_scans(iSubject, iVisit, iSession, iScan)
						fprintf(fid_summary,'"%s","%s","%s","%s"%s,\n', subjectIDs{iSubject}, visitIDs{iVisit}, imPar.sessionNames{iSession}, scanNames{iScan}, summary_lines{iSubject, iVisit, iSession, iScan});
					end
				end
			end
		end
	end
	fprintf(fid_summary,'\n');
	
	nMissing = sum(missing_scans(:));
	nSkipped = sum(skipped_scans(:));
	
	% report totals
	% header first
	fprintf(fid_summary,'\n');
	fprintf(fid_summary,'\nSubject,nConverted');
	fprintf(fid_summary,[',MissingScans (n=' num2str(nMissing) ')']);
	fprintf(fid_summary,[',SkippedScans (n=' num2str(nSkipped) ')\n']);
	% then subjects row-by-row
	for iSubject=1:nSubjects
		for iVisit=1:nVisits
			fprintf(fid_summary,'"%s"', [subjectIDs{iSubject} visitIDs{iVisit}]);
			fprintf(fid_summary,',%d',sum(converted_scans(iSubject,:,:,:)));
			
			for iSession=1:nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',scanNames{logical(missing_scans(iSubject, iVisit, iSession,:))});
				fprintf(fid_summary,'"');
			end
			
			for iSession=1:nSessions
				fprintf(fid_summary,',"');
				fprintf(fid_summary,'%s ',scanNames{logical(skipped_scans(iSubject, iVisit, iSession,:))});
				fprintf(fid_summary,'"');
			end
			
			fprintf(fid_summary,'\n');
		end
	end
	
	% and a grand total of missing and skipped
	if nMissing>0
		fprintf(2,'Number of missing scans: %d\n',nMissing);
	end
	if nSkipped>0
		fprintf(2,'Number of skipped scans: %d\n',nSkipped);
	end
	fclose(fid_summary);
	
	% cleanup
	if ~bUseDCMTK || isempty(pathDcmDict)
		dicomdict('factory');
	end
	diary('off');
	
	if ~isempty(fields(dcm2niiCatchedErrors))
		fclose all;
		SavePath = fullfile(imPar.AnalysisRoot, 'dcm2niiCatchedErrors.mat');
		SaveJSON = fullfile(imPar.AnalysisRoot, 'dcm2niiCatchedErrors.json');
		xASL_delete(SavePath);
		xASL_delete(SaveJSON);
		save(SavePath,'dcm2niiCatchedErrors');
		xASL_adm_SaveJSON(dcm2niiCatchedErrors, SaveJSON);
	end
	
	fprintf('\n');
end

%-----------------------------------------------------------------
% Create the study description output and verify that all is there
%-----------------------------------------------------------------
function datasetDescription = ExploreASL_ImportBIDS_datasetDescription(studyPar,BIDSVersion,datasetDescriptionConfig)

% Check and assign the REQUIRED fields
for iList = 1:length(datasetDescriptionConfig.Required)
	% Skip name and BIDSVersion as these were treated previously
	if isfield(studyPar,datasetDescriptionConfig.Required{iList})
		datasetDescription.(datasetDescriptionConfig.Required{iList}) = studyPar.(datasetDescriptionConfig.Required{iList});
		if strcmp(datasetDescriptionConfig.Required{iList},'BIDSVersion') && ~strcmp(BIDSVersion,studyPar.BIDSVersion)
			warning('Provided BIDS version is not consistent with the current BIDS-version-config.');
		end
	else
		% The BIDS version can be provided separately
		if strcmp(datasetDescriptionConfig.Required{iList},'BIDSVersion')
			datasetDescription.BIDSVersion = BIDSVersion;
		else
			error(['Missing a required field for datasetDescription: ' datasetDescriptionConfig.Required{iList}]);
		end
	end
end
		
% Browse through the recommended fields and report missing ones
listMissingFiles = '';
for iList = 1:length(datasetDescriptionConfig.Recommended)
	if isfield(studyPar,datasetDescriptionConfig.Recommended{iList})
		datasetDescription.(datasetDescriptionConfig.Recommended{iList}) = studyPar.(datasetDescriptionConfig.Recommended{iList});
	else
		if length(listMissingFiles)>1
			listMissingFiles = [listMissingFiles ', '];
		end
		listMissingFiles = [listMissingFiles datasetDescriptionConfig.Recommended{iList}];
	end
end
if length(listMissingFiles)>1
	% Report the missing fields
	fprintf('dataset_description.json is missing the following RECOMMENDED fields: \n%s\n',listMissingFiles);
end
	
% Browse through the optional fields and report missing ones
listMissingFiles = '';
for iList = 1:length(datasetDescriptionConfig.Optional)
	if isfield(studyPar,datasetDescriptionConfig.Optional{iList})
		datasetDescription.(datasetDescriptionConfig.Optional{iList}) = studyPar.(datasetDescriptionConfig.Optional{iList});
	else
		if length(listMissingFiles)>1
			listMissingFiles = [listMissingFiles ', '];
		end
		listMissingFiles = [listMissingFiles datasetDescriptionConfig.Optional{iList}];
	end
end

if length(listMissingFiles)>1
	% Report the missing fields
	fprintf('dataset_description.json is missing the following OPTIONAL fields: \n%s\n',listMissingFiles);
end

end

%%
% Final checking of the JSON structure, renaming and sorting last fields and checking last conditions
function jsonOut = ExploreASL_ImportBIDS_jsonCheck(jsonIn,bidsPar,bIsASL)
% Create an empty output structure and a structure with fields to delete
jsonOut = struct;
jsonRemove = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rename certain fields from the vendor-name to BIDS name

% Rename the coil names for different manufacturers
if isfield(jsonIn,'CoilString')
	switch (jsonIn.Manufacturer)
		case 'Philips'
			jsonOut.ReceiveCoilName = jsonIn.CoilString;
		case 'GE'
			jsonOut.ReceiveCoilName = jsonIn.CoilString;
		case 'Siemens'
			jsonOut.ReceiveCoilActiveElements = jsonIn.CoilString;
		otherwise
			error('Unknown manufacturer')
	end
	jsonRemove.CoilString = '';
end

% Rename the fields with number of segments
jsonRemove.NumberOfAverages = '';
if isfield(jsonIn,'NumberSegments')
	jsonOut.NumberShots = jsonIn.NumberSegments;
	jsonRemove.NumberSegments = '';
end

% Rename the phase encoding directions fields
if isfield(jsonIn,'PhaseEncodingAxis')
	if ~isfield(jsonIn,'PhaseEncodingDirection')
		jsonOut.PhaseEncodingDirection = jsonIn.PhaseEncodingAxis;
	end
	jsonRemove.PhaseEncodingAxis = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove fields not belonging to BIDS

% Remove certain empty fields
for iField = 1:length(bidsPar.listRemoveIfEmpty)
	if isfield(jsonIn,bidsPar.listRemoveIfEmpty{iField})
		if isempty(jsonIn.(bidsPar.listRemoveIfEmpty{iField}))
			jsonRemove.(bidsPar.listRemoveIfEmpty{iField}) = '';
		end
	end
end

% Remove non-BIDS fields
for iField = 1:length(bidsPar.listFieldsRemoveGeneral)
	if isfield(jsonIn,bidsPar.listFieldsRemoveGeneral{iField})
		jsonRemove.(bidsPar.listFieldsRemoveGeneral{iField}) = '';
	end
end

% Remove non-ASL-BIDS fields
if bIsASL
	for iField = 1:length(bidsPar.listFieldsRemoveASL)
		if isfield(jsonIn,bidsPar.listFieldsRemoveASL{iField})
			jsonRemove.(bidsPar.listFieldsRemoveASL{iField}) = '';
		end
	end
end

% Remove fields belonging to dataset_description
for iField = 1:length(bidsPar.datasetDescription.Required)
	if isfield(jsonIn,bidsPar.datasetDescription.Required{iField})
		jsonRemove.(bidsPar.datasetDescription.Required{iField}) = '';
	end
end
for iField = 1:length(bidsPar.datasetDescription.Recommended)
	if isfield(jsonIn,bidsPar.datasetDescription.Recommended{iField})
		jsonRemove.(bidsPar.datasetDescription.Recommended{iField}) = '';
	end
end
for iField = 1:length(bidsPar.datasetDescription.Optional)
	if isfield(jsonIn,bidsPar.datasetDescription.Optional{iField})
		jsonRemove.(bidsPar.datasetDescription.Optional{iField}) = '';
	end
end

% Go through all input fields and copy to output, but skip those in jsonRemove
for nameField = fieldnames(jsonIn)'
	if ~isfield(jsonRemove, nameField{1})
		jsonOut.(nameField{1}) = jsonIn.(nameField{1});
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check field requirements and dependencies

% Check required ASL fields
if bIsASL
	strReport = '';
	for iField = 1:length(bidsPar.ASLfields.Required)
		if ~isfield(jsonOut,bidsPar.ASLfields.Required{iField})
			if isempty(strReport)
				strReport = bidsPar.ASLfields.Required{iField};
			else
				strReport = [strReport ', ' bidsPar.ASLfields.Required{iField}];
			end
		end
	end
	if ~isempty(strReport)
		fprintf('%s\n\n',['Missing required ASL fields: ' strReport]);
	end
	
	strReport = '';
	for iField = 1:length(bidsPar.ASLfields.Recommended)
		if ~isfield(jsonOut,bidsPar.ASLfields.Recommended{iField})
			if isempty(strReport)
				strReport = bidsPar.ASLfields.Recommended{iField};
			else
				strReport = [strReport ', ' bidsPar.ASLfields.Recommended{iField}];
			end
		end
	end
	if ~isempty(strReport)
		fprintf('%s\n\n',['Missing Recommended ASL fields: ' strReport]);
	end
end

% Check ASL dependencies
if bIsASL
	for iCond = 1:length(bidsPar.ASLCondition)
		% First check if the field is present
		if isfield(jsonOut,bidsPar.ASLCondition{iCond}.field)
			
			% Checking if the condition is met, assuming no
			bCond = 0;
			if isempty(bidsPar.ASLCondition{iCond}.value)
				% Empty value means the field is only present
				bCond = 1;
			elseif ischar(bidsPar.ASLCondition{iCond}.value)
				% strings are checked with regexpi
				if regexpi(jsonOut.(bidsPar.ASLCondition{iCond}.field),bidsPar.ASLCondition{iCond}.value)
					bCond = 1;
				end
			elseif isequal(jsonOut.(bidsPar.ASLCondition{iCond}.field),bidsPar.ASLCondition{iCond}.value)
				% Logical and numbers are checked with isequal
				bCond = 1;
			end
			
			% Conditions are met, now check for dependencies
			if bCond
				% Check the required filled fields
				strReportFilled = '';
				for iField = 1:length(bidsPar.ASLCondition{iCond}.RequiredFilled)
					if ~isfield(jsonOut,bidsPar.ASLCondition{iCond}.RequiredFilled{iField}) || isempty(jsonOut.(bidsPar.ASLCondition{iCond}.RequiredFilled{iField}))
						if isempty(strReportFilled)
							strReportFilled = bidsPar.ASLCondition{iCond}.RequiredFilled{iField};
						else
							strReportFilled = [strReportFilled ', ' bidsPar.ASLCondition{iCond}.RequiredFilled{iField}];
						end
					end
				end
				
				% Check the required empty fields
				strReportEmpty = '';
				for iField = 1:length(bidsPar.ASLCondition{iCond}.RequiredEmpty)
					if isfield(jsonOut,bidsPar.ASLCondition{iCond}.RequiredEmpty{iField}) && ~isempty(jsonOut.(bidsPar.ASLCondition{iCond}.RequiredEmpty{iField}))
						if isempty(strReportEmpty)
							strReportEmpty = bidsPar.ASLCondition{iCond}.RequiredEmpty{iField};
						else
							strReportEmpty = [strReportEmpty ', ' bidsPar.ASLCondition{iCond}.RequiredEmpty{iField}];
						end
					end
				end
				
				% Check the recommended filled fields
				strReportRecommended = '';
				for iField = 1:length(bidsPar.ASLCondition{iCond}.RecommendedFilled)
					if ~isfield(jsonOut,bidsPar.ASLCondition{iCond}.RecommendedFilled{iField}) || isempty(jsonOut.(bidsPar.ASLCondition{iCond}.RecommendedFilled{iField}))
						if isempty(strReportRecommended)
							strReportRecommended = bidsPar.ASLCondition{iCond}.RecommendedFilled{iField};
						else
							strReportRecommended = [strReportRecommended ', ' bidsPar.ASLCondition{iCond}.RecommendedFilled{iField}];
						end
					end
				end
				
				% One of the dependencies was not fulfilled
				if ~isempty(strReportFilled) || ~isempty(strReportEmpty) || ~isempty(strReportRecommended)
					% Report the conditional field
					if isempty(bidsPar.ASLCondition{iCond}.value)
						fprintf('The field %s is empty, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field);
					elseif islogical(bidsPar.ASLCondition{iCond}.value)
						if bidsPar.ASLCondition{iCond}.value
							fprintf('The field %s is true, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field);
						else
							fprintf('The field %s is false, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field);
						end
					else
						fprintf('The field %s is %s, please check the dependencies below:\n',bidsPar.ASLCondition{iCond}.field,bidsPar.ASLCondition{iCond}.value);
					end
					
					% Report the incorrect dependencies
					if ~isempty(strReportFilled)
						fprintf('The required fields are missing: %s\n',strReportFilled);
					end
					if ~isempty(strReportEmpty)
						frpintf('The following fields should be empty: %s\n',strReportEmpty);
					end
					if ~isempty(strReportRecommended)
						fprintf('The recommended fields are missing: %s\n',strReportRecommended);
					end
					fprintf('\n');
				end
			end
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort fields in a predefined order

% Create the structure with the correct field order
fieldOrderStruct = [];
for iField=1:length(bidsPar.listFieldOrder)
	fieldOrderStruct.(bidsPar.listFieldOrder{iField}) = '';
end

% And sort the fields
jsonOut = xASL_adm_OrderFields(jsonOut,fieldOrderStruct);
end
