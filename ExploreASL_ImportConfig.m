function imPar = ExploreASL_ImportConfig(StudyRoot)
%ExploreASL_ImportConfig Configures the import parameters used by ExploreASL_Import
%
% FORMAT: imPar = ExploreASL_ImportConfig(StudyRoot)
% 
% INPUT: root of study folder containing DICOMs, e.g. '//MyDisk/MyStudy'
% OUTPUT: imPar which is input to ExploreASL_Import
%
% Please read the help of ExploreASL_Import for more information
% __________________________________
% Copyright 2015-2019 ExploreASL

% Just configuring, no need to check if the directory exists yet
%if nargin<1 || ~exist(StudyRoot,'dir')
%    error('Please provide path to StudyRoot folder as input argument');
%end

if strcmp(StudyRoot(end),'\') || strcmp(StudyRoot(end),'/')
    StudyRoot = StudyRoot(1:end-1); % bugfix
end

[fpath, fname, fext] = fileparts(StudyRoot);

imPar.studyID = [fname fext];
imPar.AnalysisRoot = fpath;
imPar.RawRoot = fpath;

% Does not need to check for the RAW dir at this moment
%RawDir = fullfile(StudyRoot,'raw');

%if ~exist(RawDir,'dir')
%    error(['Couldnt find raw DICOM data in ' RawDir]);
%end

% ExploreASL does not need to be initialized here
% ExploreASL_Master('',0);

% fprintf('\n%s\n\n',['Found ' RawDir ' folder...']);

%% -----------------------------------------------------------------------------
%% Initialize the study-specific parameters
%% -----------------------------------------------------------------------------
imPar.folderHierarchy = {}; % must define this per study; use a cell array of regular expressions. One cell per directory level.
imPar.tokenOrdering = []; % must match imPar.folderHierarchy: 1==subject, 2=visit, 3==session, 4==scan (if visit or session are omitted, they will be skipped)
imPar.tokenScanAliases = [];
imPar.tokenSessionAliases = [];
imPar.bMatchDirectories  = false;

% -----------------------------------------------------------------------------
% Study specific parameters
% -----------------------------------------------------------------------------
switch imPar.studyID

	case 'SABRE'
		imPar.folderHierarchy = {'^(\d.*)$', '.*', '.*', 'DICOM', '(T1_3D|pCASL_main|calibration)'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^pCASL_main$', 'ASL4D'; '^T1_3D$', 'T1'; '^calibration$', 'M0'};
		imPar.bMatchDirectories = true;    
    
	case 'SydneyMS_Controls'
		imPar.folderHierarchy = {'^(\d{3}.*\d{1})$' '^ASL$' '^(ASL4D)$'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^ASL4D$', 'ASL4D'; '^T1$', 'T1'};
		imPar.bMatchDirectories = true;      
    
	case 'SydneyMS_Controls'
		imPar.folderHierarchy = {'^(\d{3}.*\d{1})$' '^(T1)\.nii$'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^asl$', 'ASL4D'; '^T1$', 'T1'};
		imPar.bMatchDirectories = false;     
    
	case 'MS_Sidney'
		imPar.folderHierarchy = {'^(\d{3}.*\d{1})$' '^(T1)\.nii\.gz$'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^asl$', 'ASL4D'; '^T1$', 'T1'};
		imPar.bMatchDirectories = false; 

% 	case 'NOMARED'
% 		imPar.folderHierarchy = {'^(NMRD026)$' '^.*(MPRAGE|FLAIR3D).*$'};
% 		imPar.tokenOrdering = [1 0 2];
% 		imPar.tokenSessionAliases = {};
% 		imPar.tokenScanAliases = {'^MPRAGE$', 'T1'; '^FLAIR3D$', 'FLAIR'};
% 		imPar.bMatchDirectories = true;        
%         
	case 'NOMARED'
		imPar.folderHierarchy = {'^(NMRD026)$' '^.*(dASL_restSS)-(ns|ss)$'};
		imPar.tokenOrdering = [1 3 2];
		imPar.tokenSessionAliases = {'^ns$','ASL_1';'^ss$','ASL_2'};
		imPar.tokenScanAliases = {'^dASL_restSS$','ASL4D'};
		imPar.bMatchDirectories = true;
        
	case 'Leiden'
		imPar.folderHierarchy = {'^(sciad\d{3})$' '(structural|M0|asl).*'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^asl$', 'ASL4D'; '^structural$', 'T1'; '^M0$', 'M0'};
		imPar.bMatchDirectories = false; 

	case 'WSU'
		imPar.folderHierarchy = {'^(wsu-\w{3}-\d{4}-01)$' '^(T1_dicoms|asl_dicoms)'}; % {'^(wsu-\w{3}-\d{4}-01)$' '^(structural|asl)'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^asl_dicoms', 'ASL4D'; '^T1_dicoms', 'T1'};
		imPar.bMatchDirectories = true;

	case 'session_by_date' % example config when using session dates instead of IDs
		% Define a folder hierarchy by using a cell array, with regular expressions for each directory level. Use ()-brackets
		% to extract relevant subject, session and scan ID's (tokens).
		imPar.folderHierarchy = { '^(\d{3}).*', '^(\d+)$','(PSEUDO_10_min|T1W3D|M0).*\.PAR$'};
		% Define subject,session,scan order of tokens in folder hierarchy. A negative session value will ignore the actual
		% session name and replace it with a sequential ID (i.e. to replace a date with ASL_1 etc) or with the successive name
		% found in sessionNames{}
		imPar.tokenOrdering = [ 1 -2 3]; % must match imPar.folderHierarchy: entry 1==subject index, entry 2==session index, entry 3==scan index
		% Define aliases for sessions and scans by using a a cell array of size {N,2}. First column contains regular expressions
		% that should correspond with the matching criteria in imPar.folderHierarchy. The second column contains the alias.
		imPar.tokenScanAliases = { '^PSEUDO_10_min$', 'ASL4D'; '^T1W3D$', 'T1' }; % NB: The IDs 'T1' and 'ASL4D' are used below, so cannot be changed here
		imPar.bMatchDirectories = false;
		% need this if number of sessions cannot be determined from unique session names (i.e. dates)
		% sessionNames = { 'baseline', 'followup'}; < optional, if ASL_# is not applicable
		imPar.nMaxSessions = 2;

    case 'Obesitas_Nijmegen'
        imPar.folderHierarchy = {'(.*)' '.*' '.*(t1|t2_FLAIR|asl_2d)'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^asl_2d$', 'ASL4D'; '^t1$', 'T1'; '^t2_FLAIR$', 'FLAIR'};
        imPar.bMatchDirectories = true;        
        
    case 'RADAR_Sheffield'
        imPar.folderHierarchy = {'^(\d{7}_\d)$' '^(T1|ASL|FLAIR|M0)$'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^ASL$', 'ASL4D'; '^T1$', 'T1'; '^FLAIR$', 'FLAIR';  '^M0$', 'M0'};
        imPar.bMatchDirectories = true;
        
    case 'RADAR_Glasgow'
        imPar.folderHierarchy = {'^(\d{7}_\d)$' '^(T1|ASL|FLAIR)$'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^ASL$', 'ASL4D'; '^T1$', 'T1'; '^FLAIR$', 'FLAIR'};
        imPar.bMatchDirectories = true;
        
    case 'RADAR_Bristol'
        imPar.folderHierarchy = {'^(\d{7}_\d)$' '^(T1|ASL|FLAIR)$'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^ASL$', 'ASL4D'; '^T1$', 'T1'; '^FLAIR$', 'FLAIR'};
        imPar.bMatchDirectories = true;        
        
    case 'NoseWorthy'
        imPar.folderHierarchy = {'^(Sub-\d{3})$' '^(3DT1|asl)$'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^3DT1$', 'ASL4D'; '^asl$', 'T1'};
        imPar.bMatchDirectories = true;  
        
    case 'InsomniaHC'
        imPar.folderHierarchy = {'(.*_\d{4})$' '.*(T1_MPR|pcasl|M0).*\.PAR$'};
        imPar.tokenOrdering = [1 0 2]; % = [1 3 2]; for later, when we include the multiPhase
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^T1_MPR$', 'T1'; '^pcasl$', 'ASL4D'; '^M0$', 'M0'};
        imPar.bMatchDirectories = false;         
        
    case 'QSMSC'
        imPar.folderHierarchy = {'^.*(Sub-\d{3}).*$' '.*(3D_FLAIR|3D_T1|SOURCE.*pCASL)'};
        imPar.tokenOrdering = [1 0 2]; % = [1 3 2]; for later, when we include the multiPhase
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^3D_FLAIR$', 'FLAIR'; '^3D_T1$', 'T1'; '^SOURCE.*pCASL$', 'ASL4D'};
        imPar.bMatchDirectories = true;   
        
    case 'QSMSC_MultiPhase' % multiPhase, we will run this later
        imPar.folderHierarchy = {'^.*(Sub-\d{3}).*$' '.*(3D_FLAIR|3D_T1|SOURCE).*(pCASL|MultiPhase)'};
        imPar.tokenOrdering = [1 3 2];
        imPar.tokenSessionAliases = {'^pCASL$', 'ASL_1'; '^MultiPhase$', 'ASL_2'};
        imPar.tokenScanAliases = {'^3D_FLAIR$', 'FLAIR'; '^3D_T1$', 'T1'; '^SOURCE$', 'ASL4D'};
        imPar.bMatchDirectories = true;    
        
    case 'HongerWinter2'
        imPar.folderHierarchy = {'^.*(Sub-\d{3}).*$' '.*(MPRAGE|SOURCE|M0).*'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^MPRAGE$', 'T1'; '^SOURCE$', 'ASL4D'; '^M0$', 'M0'};
        imPar.bMatchDirectories = true;        
        
    case 'Dent'
        imPar.folderHierarchy = {'^.*(MCI-\d{4}).*$' '.*(FLAIR|T1W|SOURCE.*pCASL).*'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^FLAIR$', 'FLAIR'; '^T1W$', 'T1'; '^SOURCE.*pCASL$', 'ASL4D'};
        imPar.bMatchDirectories = true;
		
    case 'Miami_GE'
        imPar.folderHierarchy = {'^(RMCI_HIM_\d{3}).*$' '.*(FLAIR|T1w|ASL).*'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^FLAIR$', 'FLAIR'; '^T1w$', 'T1'; '^ASL$', 'ASL4D'};
        imPar.bMatchDirectories = true;
        
    case 'Miami_Siemens_Rockland'
        imPar.folderHierarchy = {'^(Sub-\d{3})$' '.*(MPRAGE|pCASL).*'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^MPRAGE$', 'T1'; '^pCASL$', 'ASL4D'};
        imPar.bMatchDirectories = true;
        
    case 'SPIR_artifact_Koen'
        imPar.folderHierarchy = {'^(Sub-001)$' '^.*(SOURCE_-_ASL-PSEUDO)(_1800|_1800_FrOff_135Hz|_1800_FrOff_220Hz)$'};
        imPar.tokenOrdering = [ 1 0 2];
        imPar.tokenSessionAliases = {'^_1800$', 'ASL_1'; '^_1800_FrOff_135Hz$', 'ASL_2'; '^_1800_FrOff_220Hz$', 'ASL_3'};
        imPar.tokenScanAliases = {'^SOURCE_-_ASL-PSEUDO$', 'ASL4D'};
        imPar.bMatchDirectories = true;
		
	case 'ExploreASLtest' % Test for Dicom import
		imPar.folderHierarchy = { '^(\d{3}).*', '^Session0([12])$','^(EP2D|M0|T1).*$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^1$', 'ASL_1'; '^2$', 'ASL_2' };
		imPar.tokenScanAliases = { '^EP2D*$', 'ASL4D'; '^M0$', 'M0'; '^T1_*$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'Trial_patient' % Maria Ivanova
		imPar.folderHierarchy = { '^(Subj\d{3}).*', '^Session([12])$','^(ASL|M|T1|FLAIR).*$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^1$', 'ASL_1'; '^2$', 'ASL_2' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'; '^M*$', 'M0'; '^T1_*$', 'T1';'^FLAIR$', 'FLAIR'};
		imPar.bMatchDirectories = true;
		
	case 'GENFI_Test'
		imPar.folderHierarchy = {'^(.*_C9ORF.*)$' '^.*(ASL|M0).*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^ASL$', 'ASL4D'; '^T1$', 'T1'; '^T2$', 'FLAIR'; '^M0$', 'M0' };
		imPar.bMatchDirectories = true;
		
	case 'RunDMC_inTENse'
		imPar.folderHierarchy = {'^(Sub-001)$' '^(ASL|T1)$' '^.*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^ASL$', 'ASL4D'; '^T1$', 'T1'; '^T2$', 'FLAIR'; '^M0$', 'M0' };
		imPar.bMatchDirectories = true;
		
	case 'MAGNIMS_Verona'
		imPar.folderHierarchy = {'^(Sub-001)$' '^(ASL|T1|T2|M0)$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^ASL$', 'ASL4D'; '^T1$', 'T1'; '^T2$', 'FLAIR'; '^M0$', 'M0' };
		imPar.bMatchDirectories = true;
		
	case 'Frontier'
		imPar.folderHierarchy = {'^Raw_perfustion_data_DICOM$' '^(P\d{2})$' '^(ASL|DSC)$' '^1.*$' '^1.*$' '^1.*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^ASL$', 'ASL4D'};
		imPar.bMatchDirectories = true;
	case 'FRONTIER'
		imPar.folderHierarchy = {'^(P\d{2})$' '^(ASL|DSC)$' '^1.*$' '^1.*$' '^1.*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^ASL$', 'ASL4D';'^DSC$','DSC4D'};
		imPar.bMatchDirectories = true;
	case 'BrnoEpilepsy'
		imPar.folderHierarchy = {'^(\d{4}A)$' '.*(?:PCASL|Anat).*' '.*(FLAIR|mprage|f.nii).*'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = {'^FLAIR$', 'FLAIR';'^mprage$','T1';'^f.nii$','ASL4D'};
		imPar.bMatchDirectories = false;		
	case 'BoleStudien'
		imPar.folderHierarchy = { '^(\d{3})$' '^.*(FLAIR)$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = {};
		imPar.tokenScanAliases = { '^FLAIR$', 'FLAIR'; '^M0$', 'M0'; '^T1$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'Sleep2'
		imPar.folderHierarchy = { '^(\d{3})_(1|2|3)$' '^(ASL)$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^1$', 'ASL_1'; '^2$', 'ASL_2'; '^3$', 'ASL_3' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'; '^M0$', 'M0'; '^T1$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'HIFU_ASL_Vera'
		imPar.folderHierarchy = {'^(Sub-\d{3})$' '^(ASL|T1w)$'};
		imPar.tokenOrdering = [ 1 0 2]; % In this example, there is no session (2nd token). The first brackets in imPar.folderHierarchy contain subject (1st token), and the second brackets in imPar.folderHierarchy contain the scan ID (3rd token).
		imPar.tokenSessionAliases = {''};
		imPar.tokenScanAliases = {'^T1w$','T1';'^ASL$', 'ASL4D'};
		imPar.bMatchDirectories = true;
		
	case 'Sleep_2018'
		imPar.folderHierarchy = {'^(Sub-\d{5}_\d)$' '^(ASL|MEMPRAGE RMS)$'};
		imPar.tokenOrdering = [ 1 0 2]; % In this example, there is no session (2nd token). The first brackets in imPar.folderHierarchy contain subject (1st token), and the second brackets in imPar.folderHierarchy contain the scan ID (3rd token).
		imPar.tokenSessionAliases = {''};
		imPar.tokenScanAliases = {'^MEMPRAGE RMS$','T1';'^ASL$', 'ASL4D'};
		imPar.bMatchDirectories = true;
		
	case 'Sleep_2018_2'
		imPar.folderHierarchy = { '^(\d{3})_(1|2)$' '^(ASL|T1|M0)$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^1$', 'ASL_1'; '^2$', 'ASL_2' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'; '^M0$', 'M0'; '^T1$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'CP_Tavi'
		imPar.folderHierarchy = {'^.*_\d{2}$' '^(.*_\d{2}_\d)$' '^scans$' '^\d*1$' '^.*(ASL|MPRAGE_ADNI|3DFLAIR|M0).*$'};
		imPar.tokenOrdering = [ 1 0 2]; % In this example, there is no session (2nd token). The first brackets in imPar.folderHierarchy contain subject (1st token), and the second brackets in imPar.folderHierarchy contain the scan ID (3rd token).
		imPar.tokenSessionAliases = {''};
		imPar.tokenScanAliases = {'^MPRAGE_ADNI$','T1';'^ASL$', 'ASL4D';'^3DFLAIR$', 'FLAIR';'^M0$', 'M0'};
		imPar.bMatchDirectories = true;
		
	case 'Iris4Jan_LongASL'
		imPar.folderHierarchy = { '^(Subj-\d{3})$' '^.*(3DT1|WIP_SOURCE).*(1000PLD|1800PLD|2200PLD).*\.PAR$'}; % '^.*$'
		imPar.tokenOrdering = [ 1 3 2]; % In this example, there is no session (2nd token). The first brackets in imPar.folderHierarchy contain subject (1st token), and the second brackets in imPar.folderHierarchy contain the scan ID (3rd token).
		imPar.tokenSessionAliases = {'^1000PLD$','ASL_1';'^1800PLD$','ASL_2';'^2200PLD$','ASL_3'};
		imPar.tokenScanAliases = {'^3DT1$', 'T1';'^WIP_SOURCE$', 'ASL4D'};
		imPar.bMatchDirectories = false;
		
	case 'Maarten_Lequin' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(Subj-\d{3})' '^(ASL|T1|M0|FLAIR)$' };
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1$', 'T1';'^M0$', 'M0';'^FLAIR$', 'FLAIR'};
		imPar.bMatchDirectories = true;
		
	case 'Hardy'
		imPar.folderHierarchy = { '^(HD\d{3}_\d)$' '^.*(PERFUSION_WEIGHTED)_(1|2)$'}; % '^.*$'
		imPar.tokenOrdering = [ 1 3 2]; % In this example, there is no session (2nd token). The first brackets in imPar.folderHierarchy contain subject (1st token), and the second brackets in imPar.folderHierarchy contain the scan ID (3rd token).
		imPar.tokenSessionAliases = {'1' 'ASL_1';'2' 'ASL_2'};
		imPar.tokenScanAliases = {'^FLAIR$', 'FLAIR';'^SAG_MPRAGE$', 'T1';'^PERFUSION_WEIGHTED$', 'ASL4D'};
		imPar.bMatchDirectories = true;
		
	case 'Hardy2'
		imPar.folderHierarchy = { '^(HD\d{3}_\d)$' '^.*(FLAIR|SAG_MPRAGE|ASL)_(1|2|3)$'}; % '^.*$'
		imPar.tokenOrdering = [ 1 3 2]; % In this example, there is no session (2nd token). The first brackets in imPar.folderHierarchy contain subject (1st token), and the second brackets in imPar.folderHierarchy contain the scan ID (3rd token).
		imPar.tokenSessionAliases = {'^1$' 'ASL_1';'^2$' 'ASL_2';'^3$' 'ASL_3'};
		imPar.tokenScanAliases = {'^FLAIR$', 'FLAIR';'^SAG_MPRAGE$', 'T1';'^ASL$', 'ASL4D'};
		imPar.bMatchDirectories = true;
		
	case 'Chile'
		imPar.folderHierarchy = { '^(Sub-\d{4})$' '^(ASL)$' '^DICOM' 'pCASL'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'};
		imPar.bMatchDirectories = true;
		%imPar.folderHierarchy = { '^(Sub-\d{4})$' '^(Session_1|Session_2)$'  '^(ASL|T1)$' '^DICOM' };
		%imPar.tokenOrdering = [ 1 2 3];
		%imPar.tokenSessionAliases = { 'Session_1', 'ASL_1' ;'Session_2' , 'ASL_2' };
		%imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1$', 'T1'};
		
	case 'APGEM'
		imPar.folderHierarchy = {'^(A80\d{3}_1)$' '^(PASL|T1|T2_FLAIR|T2_tse)$'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = {''};
		imPar.tokenScanAliases = {'^T2_FLAIR$', 'FLAIR';'^T1$', 'T1';'^PASL$', 'ASL4D';'^T2_tse$', 'T2'};
		imPar.bMatchDirectories = true;
		
	case 'Divers_Bonn' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(Patient\d{1})$' '^(ASL)_(Session1|Session2|Session3|Session4|Session5|Session6|Session7)$' '^.*$' '^.*$'};
		imPar.tokenOrdering = [ 1 3 2];
		imPar.tokenSessionAliases = {'Session1', 'ASL_1' ; 'Session2', 'ASL_2' ; 'Session3', 'ASL_3' ; 'Session4', 'ASL_4' ; 'Session5', 'ASL_5' ; 'Session6', 'ASL_6' ; 'Session7', 'ASL_7'};
		imPar.tokenScanAliases = {'^ASL$', 'ASL4D';'^t1$', 'T1'};
		imPar.bMatchDirectories = true;
		imPar.dcmwildcard = '*.';
		
	case 'GSP_perfusion_phantom' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(Scan8)$' '^(ASL3D_4sh_multiTI350-2600_3dyn_.*|t1_fl2d_ax_2mm_.*)$'}; % '^(1|2|3|4)$'
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '1', 'ASL_1';'2', 'ASL_2';'3', 'ASL_3';'4', 'ASL_4'};
		imPar.tokenScanAliases = { '^ASL3D_4sh_multiTI350$', 'ASL4D';'^t1_fl2d_ax_2mm_$', 'T1'};
		imPar.bMatchDirectories = true;
		imPar.dcmwildcard = '*.';
		
	case 'Nijmegen_RunDMC_Serial_Imaging' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(RUNDMCSI_\d{4}).*$' '^(asl|t1|m0)$' '^(ASL_1|ASL_2|ASL_3|ASL_4)$'};
		imPar.tokenOrdering = [ 1 3 2];
		imPar.tokenSessionAliases = {'^ASL_1$','ASL_1';'^ASL_2$','ASL_2';'^ASL_3$','ASL_3';'^ASL_4$','ASL_4'};
		imPar.tokenScanAliases = { '^asl$', 'ASL4D';'^t1$', 'T1';'^m0$','M0'};
		imPar.bMatchDirectories = true;
		imPar.dcmwildcard = '*.IMA';
		
	case 'Nijmegen_RunDMC_Serial_Imaging2' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(Test1)$' '^(LONGPLDLONGBD)$' '^(1|2|3|4)$' };
		imPar.tokenOrdering = [ 1 3 2];
		imPar.tokenSessionAliases = { '1', 'ASL_1';'2', 'ASL_2';'3', 'ASL_3';'4', 'ASL_4'};
		imPar.tokenScanAliases = { '^LONGPLDLONGBD$', 'ASL4D';'^003_t1_mprage$', 'T1';'^M0$','M0'};
		imPar.bMatchDirectories = true;
		imPar.dcmwildcard = '*.IMA';
		
	case 'Michiel_Utrecht' % 2D EPI with other direction phase
		imPar.folderHierarchy = { 'DICOM','^(PT_00002)$','^(ASL|FLAIR|T1w|M0)\.dcm$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'; '^FLAIR$', 'FLAIR'; '^T1w$', 'T1'; '^M0$', 'M0'};
		imPar.bMatchDirectories = false;
		imPar.tokenSessionAliases = { '', ''};
		imPar.dcm2nii_version = '20101105';
		
	case 'ACTL' % John Wood % CTL SCD
		% two-shot 3D grase. 2 background suppression pulses, Ldur 1600, Ldelay 2000. & MO
		imPar.folderHierarchy = { '^(ACTL\d{3})$','^(3D bs ASL Mo|3D bs ASL turboz|T1W_3d|PhaseImageASL).*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenScanAliases = { '^3D bs ASL Mo$', 'M0'; '^3D bs ASL turboz$', 'ASL4D'; '^T1W_3d$', 'T1';'PhaseImageASL','Phase4D'};
		imPar.bMatchDirectories = true;
		imPar.tokenSessionAliases = { '', ''};
		
	case 'Atle_WIP_Siemens_3DGRASE' % SI Atle WIP
		imPar.folderHierarchy = { '^(ASL_WIP_Atle)$','^(1)$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenScanAliases = { '^1$', 'ASL4D' };
		imPar.bMatchDirectories = true;
		imPar.tokenSessionAliases = { '', ''};
		imPar.dcmwildcard = '*.IMA';
    
    case 'MEDIRI' % TEST
		imPar.folderHierarchy = { '^(\d{3}).*', '^(\d{3}).*', '^(ASL|M0|MPRAGE).*$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^1$', 'ASL_1'; '^2$', 'ASL_2' };
		imPar.tokenScanAliases = { '^EP2D*$', 'ASL4D'; '^M0$', 'M0'; '^T1_*$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'Iris_unilateral_sclerosis'
		imPar.folderHierarchy = { '^PAR_REC_(\d{3})$' '^BFC_\d{3}.?_.*(O)(FF1|FF2|FF3|FF4|N1|N2|N3|N4).*\.PAR$'};
		imPar.tokenOrdering = [ 1 3 2];
		imPar.tokenSessionAliases = {'N1', 'ASL_1'; 'FF1', 'ASL_2'; 'N2', 'ASL_3'; 'FF2', 'ASL_4'; 'N3', 'ASL_5'; 'FF3', 'ASL_6'; 'N4', 'ASL_7'; 'FF4', 'ASL_8'};
		imPar.tokenScanAliases = { '^FLAIR$', 'FLAIR'; '^O$', 'ASL4D';'^MPRAGE$', 'T1';'^M0$','M0'};
		imPar.bMatchDirectories = false;
		
	case 'Iris_unilateral_sclerosis2'
		imPar.folderHierarchy = { '^PAR_REC_(\d{3})$' '^BFC_\d{3}.?_WIP_(MPRAGE|AXIAL_FLAIR).*\.PAR$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^AXIAL_FLAIR$', 'FLAIR'; '^O$', 'ASL4D';'^MPRAGE$', 'T1';'^Axial_T2-Star$','M0'};
		imPar.bMatchDirectories = false;
		
	case 'WRAP'
		imPar.folderHierarchy = { '^.*$' '^(wrap_00\d{2})_WIP_.*(pcasl|T1|M0).*\.PAR$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^pcasl$', 'ASL4D';'^T1$', 'T1';'^M0$','M0'};
		imPar.bMatchDirectories = false;
		
	case 'VICI'
		imPar.folderHierarchy = { '^.*$' '^(.*)_WIP_.*(pcasl|T1|M0).*\.PAR$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^pcasl$', 'ASL4D';'^T1$', 'T1';'^M0$','M0'};
		imPar.bMatchDirectories = false;

% 	case 'BioCog_Repro' % mosaic files from PRISMA, changed *.IMA into *.dcm
% 		imPar.folderHierarchy = { '^BioCog_(MRU\d{3}).*$' '^mri$' '^WINTERER.*$' '^.*(ep2d)(_pcasl_3x3x7_M0|_pcasl_3x3x7).*$'};
% 		imPar.tokenOrdering = [ 1 2 3];
% 		imPar.tokenSessionAliases = {'^ep2d$' 'ASL_2'};
% 		imPar.tokenScanAliases = { '^_pcasl_3x3x7_M0$', 'M0'; '^_pcasl_3x3x7$', 'ASL4D' }; % '^t1_mprage$', 'T1'; '^t2_spc_da-fl_irprep$', 'FLAIR';   NB: The IDs 'T1', 'M0' and 'ASL4D' are used below, so cannot be changed here
% 		imPar.bMatchDirectories = true;
% 		
% 	case 'BioCog' % mosaic files from PRISMA, changed *.IMA into *.dcm
% 		imPar.folderHierarchy = { '^(BICON\d{3}_\d)$' '^BioCog.*$' '^mri$' '^WINTERER.*$' '^.*(ep2d_pcasl_3x3x7_M0|ep2d_pcasl_3x3x7|t1_mprage|t2_spc_da-fl_irprep).*$'};
% 		imPar.tokenOrdering = [ 1 0 2];
% 		imPar.tokenSessionAliases = { '', ''};
% 		imPar.tokenScanAliases = { '^ep2d_pcasl_3x3x7_M0$', 'M0'; '^t1_mprage$', 'T1'; '^t2_spc_da-fl_irprep$', 'FLAIR'; '^ep2d_pcasl_3x3x7$', 'ASL4D' };
% 		imPar.bMatchDirectories = true;

	case 'BioCog_Repro' % Utrecht files
		imPar.folderHierarchy = { '^(MRU\d{3})_MRI$' '^(ASL|T1w|M0|FLAIR)\.dcm$' };
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1$', 'T1';'^M0$', 'M0';'^FLAIR$', 'FLAIR'};
		imPar.bMatchDirectories = false;

% 	case 'BioCog' % mosaic files from PRISMA, changed *.IMA into *.dcm
% 		imPar.folderHierarchy = { '^(CCC\d{3}_\d)$' '^data$' '^(ASL|T1|M0)\.PAR$' };
% 		imPar.tokenOrdering = [ 1 0 2];
% 		imPar.tokenSessionAliases = { '', ''};
% 		imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1$', 'T1';'^M0$', 'M0'};
% 		imPar.bMatchDirectories = false;

	case '22q11' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(11)$' '^(ASL|T1w)$' };
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1w$', 'T1'};
		imPar.bMatchDirectories = true;
		
	case 'BioCog' % enhanced dicom files BioCog study Utrecht
		imPar.folderHierarchy = { '^(BCU\d{3}_\d)$' '^.*(T1W|FLAIR|pCASL|M0_17).*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { ''};
		imPar.tokenScanAliases = { '^T1W$', 'T1'; '^FLAIR$', 'FLAIR'; '^pCASL$', 'ASL4D'; '^M0_17$','M0'};
		imPar.bMatchDirectories = true;
		
	case 'BioCogSiemens' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^BioCog_(BIM\d{3}_pre_\d{6})$' '^.*$' '^MRI$' '^BIM\d{3}_.*$' '^(003_t1_mprage|_.*|014_ep2d_pcasl_3x3x7|012_ep2d_pcasl_3x3x7_M0).*$' };
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^014_ep2d_pcasl_3x3x7$', 'ASL4D';'^003_t1_mprage$', 'T1';'^012_ep2d_pcasl_3x3x7_M0$','M0'};
		imPar.bMatchDirectories = true;
		
	case 'XNAT_BioCog' % enhanced dicom files BioCog study Utrecht
		imPar.folderHierarchy = { '^(BIM\d{3}_\d{1})$' '^scans$' '^(ASL|M0|FLAIR)$' '^DICOM$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { ''};
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^FLAIR$', 'FLAIR';'^T1$', 'T1';'^M0$','M0'};
		imPar.bMatchDirectories = true;
		
	case 'DDI_IVS' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(D.*-1)_.*\d{8}$','^\D.*(T1|PCASL|M0).*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^PCASL$', 'ASL4D';'^T1$', 'T1';'^M0$','M0'};
		imPar.bMatchDirectories = true;
		
	case 'Hongerwinter'
		imPar.folderHierarchy = { '^(SOE_\d{4})$','^.*$','^MR$','^\d*_(MPRAGE_ADNI|ASL-PSEUDO|3D_BV_FLAIR.*)$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^ASL-PSEUDO$', 'ASL4D';'^MPRAGE_ADNI$', 'T1';'^3D_BV_FLAIR.*$','FLAIR'};
		imPar.bMatchDirectories = true;
		
	case 'CPC_NEURO_NPH'
		imPar.folderHierarchy = { '^(S\d{5})$','^S(8|9)010$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^8$', 'ASL4D';'^9$', 'M0'};
		imPar.bMatchDirectories = true;
		imPar.dcmwildcard         = '*.';
		
	case 'Siemens_PASL_DCM' % mosaic files from PRISMA, changed *.IMA into *.dcm
		imPar.folderHierarchy = { '^(Subject1)_(ASL)$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^t1$', 'T1'};
		imPar.bMatchDirectories = true;
		
	case 'GENFI_T1MatrixTrial'
		imPar.folderHierarchy = { '^(GRN018|GRN037|GRN078)$','^GENFI_(t1)$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '1', 'ASL_1'};
		imPar.tokenScanAliases = { '^t1$', 'T1'};
		imPar.bMatchDirectories = true;
		
	case 'CRUISE_pCASL_artefact'
		imPar.folderHierarchy = { '^(01303)_(ASL)$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { ''};
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'};
		imPar.bMatchDirectories = true;
           
	case 'GENFI_DF12' % SI Trio NSA5 ASL
		imPar.folderHierarchy = { '^(SI_Trio_NSA5_C9ORF037_2)$','^(GENFI_asl)_(tag|control|tag_2|control_2|tag_3|control_3|tag_4|control_4|tag_5|control_5)$'};
		imPar.tokenOrdering = [ 1 3 2];
		imPar.tokenSessionAliases = { '^tag$', 'ASL_1'; '^control$', 'ASL_2'; '^tag_2$', 'ASL_3'; '^control_2$', 'ASL_4'; '^tag_3$', 'ASL_5'; '^control_3$', 'ASL_6'; '^tag_4$','ASL_7'; '^control_4$', 'ASL_8'; '^tag_5$','ASL_9'; '^control_5$', 'ASL_10' };
		imPar.tokenScanAliases = { '^GENFI_asl$', 'ASL4D' };
		imPar.bMatchDirectories = true;
	
% 	case 'GIFMI_STUMOOD' % SI ASL Gent ASL sessions
% 		for ii=1:38
% 			if      ii==1
% 				NUMBERseries='01';
% 				imPar.tokenSessionAliases{ii,1}='^01$';
% 				imPar.tokenSessionAliases{ii,2}='ASL_1';
% 			elseif  ii<10
% 				NUMBERseries=[NUMBERseries '|0'   num2str(ii)];
% 				imPar.tokenSessionAliases{ii,1}=['^0'   num2str(ii) '$'];
% 				imPar.tokenSessionAliases{ii,2}=['ASL_' num2str(ii)];
% 			else    NUMBERseries=[NUMBERseries '|'    num2str(ii)];
% 				imPar.tokenSessionAliases{ii,1}=['^'    num2str(ii) '$'];
% 				imPar.tokenSessionAliases{ii,2}=['ASL_' num2str(ii)];
% 			end
% 		end
% 		imPar.folderHierarchy = { '^(SUB\d{3})$','^SUB\d{3}$',['^(TE00_TI1700)_(' NUMBERseries ')$']};
% 		imPar.tokenOrdering = [ 1 3 2];
% 		imPar.tokenScanAliases = { '^TE00_TI1700$', 'ASL4D' };
% 		imPar.bMatchDirectories = true;
% 		imPar.dcmwildcard = '*.IMA';

	case 'GIFMI_STUMOOD' % SI ASL Gent T1/M0
		imPar.folderHierarchy = { '^(SUB\d{3})$','^SUB\d{3}$','^(T1_MPRAGE|SS_TE00_TI5000)_.*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenScanAliases = { '^T1_MPRAGE$', 'T1'; '^SS_TE00_TI5000$', 'M0'};
		imPar.tokenSessionAliases     = {};
		imPar.bMatchDirectories = true;
		imPar.dcmwildcard = '*.IMA';
		 		
% 	case 'Retest6mnd'
% 		% sessions 1=no crushed, 2=crushed
% 		%                 imPar.RawRoot = 'C:\Users\lvaclavu\prediva_testdataset\raw';
% 		%imPar.folderHierarchy = { '^\d{8}_(\d{6})$' '(flair|ADNI|pseudo_crush|crush_m0)$'};
% 		imPar.folderHierarchy = { '^\d{8}_(\d{6})$', '(T1)$'};
% 		imPar.tokenOrdering = [1 0 2];
% 		imPar.tokenSessionAliases = { '^nocrush$','ASL_1' };
% 		%imPar.tokenSessionAliases = { '^pseudo_crush$','ASL_2' };
% 		imPar.tokenScanAliases = {  '^T1$', 'T1' }; % NB: The IDs 'T1' and 'ASL4D' are used below, so cannot be changed here
% 		imPar.bMatchDirectories = true;
                
	case 'FollowUp'
		% sessions 1=no crushed, 2=crushed
		% imPar.RawRoot = 'C:\Users\lvaclavu\prediva_testdataset\raw';
		% imPar.folderHierarchy = { '^\d{8}_(\d{6})$' '(flair|ADNI|pseudo_crush|crush_m0)$'};
		imPar.folderHierarchy = { '^\d{8}_(\d{6})$', '^(noC|C)(rush)$'};
		imPar.tokenOrdering = [1 2 3];
		imPar.tokenSessionAliases = { '^noC$','ASL_1';'^C$','ASL_2' };
		%imPar.tokenSessionAliases = { '^pseudo_crush$','ASL_2' };
		imPar.tokenScanAliases = {  '^rush$','ASL4D' }; % NB: The IDs 'T1' and 'ASL4D' are used below, so cannot be changed here
		imPar.bMatchDirectories = true;
		
	case 'FInd'
		imPar.folderHierarchy = {'^(\d{3})_1186_FIND$', '^20$','^MR$' '^(\d{4})_(pseudo)_(nocrush)$'};
		imPar.tokenOrdering   = [1 3 2];
		imPar.tokenSessionAliases = { '^nocrush$','ASL_1' };
		imPar.tokenScanAliases = { '^pseudo$', 'ASL4D'}; % NB: The IDs 'T1' and 'ASL4D' are used below, so cannot be changed here
		imPar.bMatchDirectories = true;
                
% 	case 'CRUISE'
% 		imPar.folderHierarchy = { '^(\d{3})', '^(before|after|ramp)$','^(PSEUDO|T1|M0|3d_flair)$'};
% 		imPar.tokenOrdering = [1 2 3];
% 		imPar.tokenSessionAliases = { '^before$', 'before'; '^ramp$', 'ramp';'^after$','after'};
% 		imPar.tokenScanAliases = { '3d_flair$', 'FLAIR';'^PSEUDO$', 'ASL4D'; '^M0$', 'M0'; '^T1$', 'T1' };
% 		imPar.bMatchDirectories = true;

	case 'CRUISE'
		imPar.folderHierarchy = { '^CRUISE_(\d{3})$' '^20' '^(MR)$' '(PSEUDO_1800|T2W_TSE|flair|M0)$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^MR$','ASL_1' };
		imPar.tokenScanAliases = { 'flair$', 'FLAIR'; '^PSEUDO_1800$', 'ASL4D'; '^M0$', 'M0'; '^T2W_TSE$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case '3D_Antipsychotics'
		imPar.folderHierarchy = { '^(\d{3})_(0|S)$','^(3D structural|Mzero|PCASL)$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^0$', 'ASL_1'; '^S$', 'ASL_2' };
		imPar.tokenScanAliases = { '^PCASL$', 'ASL4D'; '^Mzero$', 'M0'; '^3D structural$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case {'Trial_ePOD_Paul',}
		imPar.folderHierarchy = { 'ASL|T1', '^(\d{3}).*', '^sessie([12])$','(PSEUDO_10_min|T1W3D|M0).*\.PAR$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^1$', 'ASL_1'; '^2$', 'ASL_2' };
		imPar.tokenScanAliases = { '^PSEUDO_10_min$', 'ASL4D'; '^T1W3D$', 'T1' }; 
		imPar.bMatchDirectories = false;
		
	case 'ASL_epod_MPH'
		imPar.RawRoot  = '/mnt/l-disk/Projects/epod/ePOD-MPH/Post-treatment/MRI/Original';
		imPar.AnalysisRoot = '/scratch/pfgroot/ASL_epod_MPH/pt';
		imPar.folderHierarchy = { '^(\d{3}).*', '^sessie([12])$','(PSEUDO_10_min|T1W3D|M0).*_1\.PAR$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^1$', 'ASL_1'; '^2$', 'ASL_2' };
		imPar.tokenScanAliases = { '^PSEUDO_10_min$', 'ASL4D'; '^T1W3D$', 'T1' }; 
		imPar.bMatchDirectories = false;
		imPar.dcm2nii_version = '20101105';
		
	case 'Score'
		imPar.folderHierarchy = { '^(\d{3}_\d)$', '^.*(MPRAGE|PSEUDO|FLAIR).*$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '' };
		imPar.tokenScanAliases = { '^PSEUDO$', 'ASL4D'; '^MPRAGE$', 'T1'; '^FLAIR$', 'FLAIR' }; 
		imPar.bMatchDirectories = true;
                
	case 'Score_2'  % .../ASL/SCORE_000_1/scans/601/DICOM
		imPar.folderHierarchy = { '^(MPRAGE)$', '^SCORE_(\d{3})$', '^scans$','^\d*$','^DICOM$'};
		imPar.tokenOrdering = [ 2 0 1];
		imPar.tokenSessionAliases = { '' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'; '^MPRAGE$', 'T1' };
		imPar.bMatchDirectories = true;
		% imPar.dcm2nii_version = '20101105';
		
	case 'ePOD_AMPH'
		imPar.folderHierarchy = { '^(session1)$', '^(ASL)$', '^AMPHMRI_(\d{3})_WIP.*PCASL.*\.PAR$'};
		imPar.tokenOrdering = [ 3 1 2];
		imPar.tokenSessionAliases = { '^session1$', 'ASL_1' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D' };
		imPar.bMatchDirectories = false;
		
    case 'Novice'
        imPar.folderHierarchy = {'^(NOV\d{3}_\d)$' '.*(MPRAGE|FLAIR|SOURCE.*ASL).*' '.*'}; % '^MR$'
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {};
        imPar.tokenScanAliases = {'^MPRAGE$', 'T1'; '^SOURCE.*ASL$', 'ASL4D'; '^FLAIR$', 'FLAIR'};
        imPar.bMatchDirectories = true;        
        
	case 'NOVICE'
		imPar.folderHierarchy = { '^(session1)$' '^(NOV\d{3})$' '^(ASL)$' '^\d+$'};
		imPar.tokenOrdering = [ 2 1 3];
		imPar.tokenSessionAliases = { '^session1$', 'ASL_1' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D' };
		imPar.bMatchDirectories = true;
		
	case 'NOVICE2'
		imPar.folderHierarchy = { '^session1$' '^(T1)$' '^(NOV\d{3})$' '^scans$' '^\d+$' '^DICOM$'};
		imPar.tokenOrdering = [3 0 2];
		imPar.tokenSessionAliases = { '', '' };
		imPar.tokenScanAliases = { '^T1$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'Antipsychotics'
		imPar.folderHierarchy = { '^Antipsychotics_Sorted_for_Analysis$','^(\d{3})_([0|S])','^(ASL|M0|T1)$'};
		imPar.tokenOrdering = [ 1 2 3];
		imPar.tokenSessionAliases = { '^0$', 'ASL_1'; '^S$', 'ASL_2' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D'; '^M0$', 'M0'; '^T1$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'preDIVA_followUp1'
		imPar.folderHierarchy = { '^\d{8}_(\d{6})$' '^(no)(Crush)$'};
		imPar.tokenOrdering = [1 2 3];
		imPar.tokenSessionAliases = { '^no$', 'ASL_1' };
		imPar.tokenScanAliases = { '^Crush$', 'ASL4D' };
		imPar.bMatchDirectories = true;
		
	case 'preDIVA_followUp2'
		imPar.folderHierarchy = { '^\d{8}_(\d{6})$' '^(C)(rush)$'};
		imPar.tokenOrdering = [1 2 3];
		imPar.tokenSessionAliases = { '^C$', 'ASL_2' };
		imPar.tokenScanAliases = { '^rush$', 'ASL4D' };
		imPar.bMatchDirectories = true;
		
	case 'preDIVA_followUp3'
		imPar.folderHierarchy = { '^\d{8}_(\d{6})$' '^(no)(CrushM0)$'};
		imPar.tokenOrdering = [1 2 3];
		imPar.tokenSessionAliases = { '^no$', 'ASL_1' };
		imPar.tokenScanAliases = { '^CrushM0$', 'M0' };
		imPar.bMatchDirectories = true;
		
	case 'preDIVA_followUp4'
		imPar.folderHierarchy = { '^\d{8}_(\d{6})$' '^(Crush)(M0)$'};
		imPar.tokenOrdering = [1 2 3];
		imPar.tokenSessionAliases = { '^Crush$', 'ASL_2' };
		imPar.tokenScanAliases = { '^M0$', 'M0' };
		imPar.bMatchDirectories = true;
		
	case 'preDIVA_baselinePipeline'
		imPar.RawRoot = 'L:\basic\divi\Projects\prediva\DICOM';
		imPar.folderHierarchy = { '^\d{8}_(\d{6})$' '(flair|ADNI)'};
		imPar.tokenOrdering = [1 0 2];
		imPar.tokenSessionAliases = { '', '' };
		imPar.tokenScanAliases = { '^flair$', 'FLAIR' ; '^ADNI$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'preDIVA_baselinePipeline2'
		imPar.RawRoot = 'L:\basic\divi\Projects\prediva\DICOM';
		imPar.folderHierarchy = { '^(\d{8}_\d{6})$' '^(noC|C)(rush)$'};
		imPar.tokenOrdering = [1 2 3];
		imPar.tokenSessionAliases = { 'noC', 'ASL_1';'C', 'ASL_2' };
		imPar.tokenScanAliases = { '^rush$', 'ASL4D' };
		imPar.bMatchDirectories = true;
		
	case 'GE_trial'
		imPar.folderHierarchy = { '^(HART\d{2})(A)$' '^(T1|ASL|M0)$'};
		imPar.tokenOrdering = [1 2 3];
		imPar.tokenSessionAliases = { 'A', 'ASL_1' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D' ; '^T1$', 'T1' ; '^M0$', 'M0' };
		imPar.bMatchDirectories = true;
		
	case 'COBRA_ICL_BL'
		imPar.RawRoot = fullfile('/home/pfcgroot/lood_storage/divi/Projects/ageiv/MRI/raw/COBRA/ICL_baseline/DICOM');
		imPar.folderHierarchy = {'^COB_I(\d{3})_v1$', '^scans$', '^.*(MPRAGE|pasl).*$', '^resources$', '^DICOM$', '^files$' };
		imPar.tokenOrdering   = [1 0 2];
		imPar.tokenSessionAliases = { '', '' };
		imPar.tokenScanAliases = { '^pasl$', 'ASL4D' ; '^MPRAGE$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'COBRA_ICL_FU'
		imPar.RawRoot = fullfile('/home/pfcgroot/lood_storage/divi/Projects/ageiv/MRI/raw/COBRA/ICL_followup/DICOM');
		imPar.folderHierarchy = {'^COB_I(\d{3})_v2$', '^scans$', '^.*(MPRAGE|pasl).*$', '^resources$', '^DICOM$', '^files$' };
		imPar.tokenOrdering   = [1 0 2];
		imPar.tokenSessionAliases = { '', '' };
		imPar.tokenScanAliases = { '^pasl$', 'ASL4D' ; '^MPRAGE$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'COBRA_AMC'
		imPar.RawRoot = '/scratch/mwcaan/AgeIV/asl';
		imPar.folderHierarchy = { 'CNS_(\d{3}_..._\d)$', '^scans$', '^\d+_(T1|ASL)$', '^DICOM$'};
		imPar.tokenOrdering   = [1 0 2];
		imPar.tokenSessionAliases = { '', '' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D' ; '^T1$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'Parelsnoer_ASL_HJ'
		imPar.folderHierarchy = {'^\d{3}$' '^(\d{3}.*)$' '^.*(ASL|T1)$'};
		imPar.tokenOrdering   = [1 0 2];
		imPar.tokenSessionAliases = { '', '' };
		imPar.tokenScanAliases = { '^ASL$', 'ASL4D' ; '^T1$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'INOX'
		imPar.folderHierarchy = {'^INOX(\d{5})$' '^(000)(2|3)$'}; %
		imPar.tokenOrdering   = [1 3 2];
		imPar.tokenSessionAliases = { '2', 'ASL_1';'3', 'ASL_2' };
		imPar.tokenScanAliases = { '^000$', 'ASL4D' };
		imPar.bMatchDirectories = true;
		
	case 'INOX2'
		imPar.folderHierarchy = {'^INOX(\d{5})$' '^(0007|0008|0015)$'}; %
		imPar.tokenOrdering   = [1 0 2];
		imPar.tokenSessionAliases = { '' };
		imPar.tokenScanAliases = { '^0007$', 'T1';'^0007$', 'T1';'^0015$', 'T1' };
		imPar.bMatchDirectories = true;
		
	case 'FIND_STUDIE'
		imPar.folderHierarchy = {'^(\d{3})_1186_FIND$' '^\d{8}_$' '^.*$' '^MR$' '^0\d{3}1_(pseudo)_(nocrush)$'}; %
		imPar.tokenOrdering   = [1 3 2];
		imPar.tokenSessionAliases = { '^nocrush$','ASL_1' };
		imPar.tokenScanAliases = { '^pseudo$', 'ASL4D'};
		imPar.bMatchDirectories = true;
		
	case {'WMH_AGEIV',}
		imPar.folderHierarchy = { '(CNS_\d{3}_...).*$', '^(wnu|lesprob_thr_90_wnu).*\.nii\.gz$'};
		imPar.tokenOrdering = [ 1 0 2];
		imPar.tokenSessionAliases = { '', ''};
		imPar.tokenScanAliases = { '^wnu$', 'FLAIR'; '^lesprob_thr_90_wnu$', 'WMH_SEGM' };
		imPar.bMatchDirectories = false;
		
    case 'Sagittal_Sinus'
        imPar.folderHierarchy = {'^(\d{3})_1186_FIND$', '^\d{3}_1186_FIND$', '^20.*$', '^MR$', '^\d{5}_(pseudo_nocrush|T2_TRA$)'};
        imPar.tokenOrdering = [1 0 2];
        imPar.tokenSessionAliases = {'',''};
        imPar.tokenScanAliases = {'^pseudo_nocrush$','ASL4D';'^T2_TRA$', 'T1'}; 
        imPar.bMatchDirectories = true;

	otherwise % by default, append study ID and raw (or analysis) to the roots
		warning('Unknown study: %s', imPar.studyID);
        fprintf('%s\n','Could not define import parameters, please make sure that they are added to ExploreASL_ImportConfig.m');
        fprintf('%s\n','Also make sure that the root foldername is equal to the StudyID in ExploreASL_ImportConfig.m');
end


end
