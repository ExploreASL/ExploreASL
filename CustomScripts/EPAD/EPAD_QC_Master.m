function EPAD_QC_Master(ROOT, bUseDCMTK, bCheckPermissions, bRunExploreASL)
%EPAD_QC_Master Runs EPAD QC pipeline
%
% FORMAT: EPAD_QC_Master(ROOT[, bUseDCMTK, bCheckPermissions, bRunExploreASL])
% 
% INPUT:
%   ROOT              - path to root folder containing the data. The raw DICOM data should be in //ROOT/raw, then ExploreASL_Import will create
%                       //ROOT/source (BIDS NIfTI data) and //ROOT/analysis (NIfTI data with derivatives) (REQUIRED)
%   bUseDCMTK         - true for replacing Matlab's DICOMINFO by Jan Petr's compilation of DCM Tk, which goes much faster,
%                       but may not be as extensive as DICOMINFO (i.e. exceptions need to be defined by Jan Petr) (OPTIONAL, DEFAULT=false)
%   bCheckPermissions - Check whether file permissions are set correctly.
%                       Checks each dicom individually, hence very accurate but very slow. Can be useful when wanting to run this script without SUDO (OPTIONAL, DEFAULT=false)
%                       This can be useful to notify if a script crashes because of permission issues. Data, including the ROOT folder) needs to be readable and
%                       writable for users and groups (i.e. at least 660). For ExploreASL executables/data folders we want to read, write and execute for users and groups (i.e. at least 770)
%                       xASL_adm_CheckPermissions will try resetting this, and throw a warning if it cannot.
%   bRunExploreASL    - set to true to run ExploreASL (optional, default = TRUE)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: 
% This script contains the EPAD QC pipeline, which is based on ExploreASL, Matlab, SPM, SPM U+ and WAD-QC
% ExploreASL generic scripts are prefixed with xASL_*, EPAD QC specific scripts with EPAD_*
% 
% This pipeline will run the following steps:
% 0) Admin: manages paths, permissions, etc. If data are found inside //ROOT they are copied to //ROOT/raw
%    NB: NEED TO IMPLEMENT: cp --preserve=mode,ownership
%
% 1) EPAD_dcmCuration: this function will curate, clean, rename and sort the raw DICOMs in //ROOT/raw
%    dicoms are typically downloaded from IXICO and structured by Viktor
%    using IXICO tags (type help EPAD_dcmCuration for more information)
%    Curation is performed according to ScanType_LabelsConfig.csv, which is
%    overwritten from //ExploreASL/CustomScripts/EPAD.
%    RESULTS IN: DCM_StructureList.csv & DCM_StructureList.mat
%    SKIPPED IF: //ROOT/raw/DCM_StructureList.mat exists
%
% 2) dcm2niiX: this part runs ExploreASL_Import for converting dicoms in
%    //ROOT/raw (curated DICOM data) to //ROOT/source (BIDS-compliant NIfTI data) and //ROOT/analysis (derivatives)
%    Additionally, this function will copy single dummy DICOMs for each ScanType to //ROOT/analysis to embed QC results in for WAD-QC
%    RESULTS IN: import_summary.csv
%    SKIPPED IF: import_summary.csv exists
%    SUBJECTS SKIPPED IF: //ROOT/analysis/SubjectFolder already exists
%
% 3) Additional modifications
%    Here, we do some additional fixes in //ROOT/analysis (PM: THESE MODIFICATIONS TO BE CLONED TO //ROOT/source AS WELL
%    EPAD_BIDS_Fix_SWI           - Fix dcm2nii SWI conversion for NIfTI with both Magnitude&Phase contrast
%    EPAD_BIDS_Fix_DTI           - fix ADC conversion & create NormPE NIfTI
%    EPAD_BIDS_Fix_PE            - Combine separate NormPE/RevPE dirs into single dir per contrast
%    EPAD_Manage_PEPolar_Parms   - put the correct PEPolar parameters in the JSON files
%    EPAD_BIDS_Fix_ASL           - Fix dcm2nii ASL conversion errors
%    EPAD_ASL_parmsPrepare       - Prepare different ASL parameter files for different vendors
%    EPAD_CreateASLJSONPars      - List ASL sequence parameters & populate the JSONs for quantification
%    EPAD_CopyFLAIR_WMH_Carole   - Copy WMHs segmented by Carole, making sure to also copy the FLAIRs that are in alignment with the WMH
%    EPAD_ReportMissingFiles     - Creates report of missing data, according to ScanType_LabelsConfig.csv
%    PM: let xASL_adm_CreateFileReport provide option for
%        ScanType_LabelsConfig table, and only provide derivative results for files of which the source files exist
%
% 4) ExploreASL_Master_EPAD This function is an EPAD extension of ExploreASL_Master,
%    by including QC modules xASL_module_func, xASL_module_dwi.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_QC_Master('/data/RAD/share/EPAD500_new');
% __________________________________
% Copyright 2017-2019 Amsterdam UMC
% Henk-Jan Mutsaerts, Silvia Ingala, Viktor Wottschel, Alle Meije Wink, Joost Kuijer


%% Deal with input arguments
if nargin < 2 || isempty(bUseDCMTK)
	bUseDCMTK = false;
end

if nargin < 3 || isempty(bCheckPermissions)
    if isunix || ismac
        bCheckPermissions = false; % set this to TRUE LATER!!!!!!
    else
        bCheckPermissions = false;
    end
end
if nargin<4 || isempty(bRunExploreASL)
    bRunExploreASL = true;
end

% Ensure we start in the EPAD folder
if exist(fullfile(pwd,'ExploreASL_Master.m')) % we are in the ExploreASL root folder
    cd(fullfile('CustomScripts','EPAD'));
end
if ~exist(fullfile(pwd,'EPAD_QC_Master.m'))
    error('Please start this script in the EPAD CustomScripts folder');
end

%% ========================================================
%% 0) Admin & copy data
RawDir = fullfile(ROOT, 'raw');
AnalysisDir = fullfile(ROOT, 'analysis');
SourceDir = fullfile(ROOT, 'source');

x = EPAD_InitxASL; % call ExploreASL for initialization only

ScanType_ConfigPath = 'ScanType_LabelsConfig.csv';
ScanType_ConfigPath2 = fullfile(RawDir, 'ScanType_LabelsConfig.csv');
xASL_Copy(ScanType_ConfigPath, ScanType_ConfigPath2, true);

bCopy = false;
if ~exist(RawDir, 'dir' )
    xASL_adm_CreateDir(RawDir); % to copy DICOMs from //ROOT to //ROOT/raw
end

DirList = xASL_adm_GetFsList(ROOT, '^\d{3}$', 1, 1, [], [0 Inf]); % get folders in //ROOT
DirListRaw = xASL_adm_GetFsList(RawDir, '^\d{3}$', 1, 1, [], [0 Inf]);  % get folders in //ROOT/raw

if isempty(DirListRaw) && ~isempty(DirList)
    bCopy = true; % to copy DICOMs to //ROOT/raw
elseif isempty(DirList) && isempty(DirListRaw)
    error('No data found'); % give error if no DICOM data found at all
end

if bCopy % copy DICOMs from //ROOT to //ROOT/raw
    if bCheckPermissions; xASL_adm_CheckPermissions(RawDir); end
    for iD=1:length(DirList)
        DirPath = fullfile(ROOT, DirList{iD});
        DirPathNew = fullfile(RawDir, DirList{iD});
        if bCheckPermissions; xASL_adm_CheckPermissions(DirPath, 0); end
        fprintf('%s\n', 'Copying data:...');
        EPAD_CopyDirsTrackProgress(DirPath, DirPathNew); % wrapper around xASL_Copy to track progress
    end
end


% Add this folder to path & go to ExploreASL root
addpath(genpath(pwd));
cd ../..;

%% ========================================================
%% 1) Run the curation script, to sort the DICOMs
fprintf('\n\n%s\n', '----------------------------------------------');
fprintf('%s\n', '1) EPAD_dcmCuration');
StatusPath = fullfile(RawDir, 'DCM_StructureList.mat'); % ->>>> DELETE THIS FILE TO REPEAT SORTING
if ~exist(StatusPath, 'file')
    EPAD_dcmCuration(RawDir, bUseDCMTK, bCheckPermissions);
end


%% ========================================================
%% 2) run DICOM2NII/BIDS
fprintf('\n%s\n', '----------------------------------------------');
fprintf('%s\n', '2) DICOM2NII/BIDS');

StatusPath = fullfile(AnalysisDir, 'import_summary.csv'); % ->>>> DELETE THIS FILE TO REPEAT DCM2NII CONVERSION
if ~exist(StatusPath, 'file')
    % now start reading the parameters from ScanType_LabelsConfig.csv
    % & putting them in import parameter struct imPar 
    [ScanTypeConfig] = xASL_csvRead(fullfile(RawDir, 'ScanType_LabelsConfig.csv'));
    for iL=1:size(ScanTypeConfig,1)
        BIDSlist{iL,1} = [ScanTypeConfig{iL,2} '_' ScanTypeConfig{iL,3}];
    end
    % Deal with empty cells
    xASLList = ScanTypeConfig(:,4);
    NonEmptyBIDS = cellfun(@(x) ~isempty(x), BIDSlist);
    NonEmptyXASL = cellfun(@(x) ~isempty(x), xASLList);
    NonEmptyList = NonEmptyBIDS & NonEmptyXASL;
    BIDSlist = BIDSlist(NonEmptyList);
    xASLList = xASLList(NonEmptyList);
    ConcatBIDS = '.*(';
    for iL=2:length(BIDSlist)
        ConcatBIDS = [ConcatBIDS BIDSlist{iL}];
        TokenScan{iL-1,1} = ['^' BIDSlist{iL} '$'];
        TokenScan{iL-1,2} = xASLList{iL};
        if iL~=length(BIDSlist)
            ConcatBIDS = [ConcatBIDS '|']; % "|" does the either or concatenation
        end
    end
    ConcatBIDS = [ConcatBIDS ').*']; % BIDS name for scantypes

    % Define ExploreASL_Import settings
    [Fpath, Ffile, Fext] = fileparts(ROOT);
    imPar.studyID = [Ffile Fext];
    imPar.AnalysisRoot = Fpath;
    imPar.RawRoot = Fpath;
    imPar.tokenOrdering = [1 2 0 3]; % Subject Visit Session Scan
    imPar.tokenSessionAliases = {};
    imPar.tokenVisitAliases = {'Screening','_1'; 'Month_12','_2'; 'Month_24','_3'; 'Month_36','_4'; 'Month_48','_5'};
    imPar.tokenScanAliases = TokenScan;
    imPar.bMatchDirectories = true;
    imPar.bVerbose = false;
    imPar.bOverwrite = false; % if we want to redo stuff, best to simply remove the analysisDir
    imPar.SkipSubjectIfExists = true;
    bCopySingleDicoms = true; % copies also a dummy DICOM copy for each NIfTI ScanType, to allow retracing
    % original DICOM parameters, and as dummy DICOM to embed QC parameters in for WAD-QC
    bUseDCMTK = false; % speed up process by using DCMTK compilation instead of Matlab's dicominfo
    bClone2Source = false; % copies NIfTI analysis folders (derived scans) also to source folder
    % that can go to XNAT, before continuing processing below. But let's do
    % this now only after full import & JSON wrappers were completed
    
    % Now throw error or warning depending on pre-existing TimePoints
    if ~isempty(xASL_adm_GetFileList(AnalysisDir, '^\d{3}EPAD\d*$','FPList',[0 Inf],true))
        error('Pre-existing analysis-subject folders without TimePoint suffix');
    elseif ~isempty(xASL_adm_GetFileList(AnalysisDir, '^\d{3}EPAD\d*_\d*$','FPList',[0 Inf],true))
        warning('Make sure that all existing analysis-subject folders have the correct suffix');
    end
    
    % Define FolderHierarchy
    imPar.folderHierarchy = {'^\d{3}$' '^(\d{3}EPAD\d*)$' '^(Screening|Month_12|Month_24|Month_36|Month_48)$' ConcatBIDS '^(?!inactive).*(?!recon)$'};
    % Note that this folder layout states that we ignore DICOM folders containing the text "Inactive" and "recon".
    % This avoids converting derivative DICOMs to NIfTI
    % If only inactive exists that is complete (as was set by IXICO),
    % then this is set as active by the curation script above.

    % run import script
    ExploreASL_Import(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, [], bClone2Source, x);
end


%% ========================================================
%% 3) Some extra modifications
fprintf('\n%s\n', '----------------------------------------------');
fprintf('%s\n', '3) Some modifications: fixes & preparations');

if bCheckPermissions; xASL_adm_CheckPermissions(AnalysisDir, false); end % to be sure

EPAD_BIDS_Fix_SWI(AnalysisDir); % Fix dcm2nii SWI conversion for NIfTI with both Magnitude&Phase contrast
EPAD_BIDS_Fix_DTI(AnalysisDir); % fix ADC conversion & create NormPE NIfTI
EPAD_BIDS_Fix_PE(AnalysisDir); % Combine separate NormPE/RevPE dirs into single dir per contrast
EPAD_Manage_PEPolar_Parms(AnalysisDir); % put the correct PEPolar parameters in the JSON files
EPAD_BIDS_Fix_ASL(AnalysisDir); % Fix dcm2nii ASL conversion errors
% EPAD_ASL_parmsPrepare(AnalysisDir); % REPLACE THIS FUNCTION BY EPAD_CREATEASLJSONPARS.M
EPAD_CreateASLJSONPars(AnalysisDir); % List ASL sequence parameters & populate the JSONs for quantification
EPAD_CopyFLAIR_WMH_Carole(AnalysisDir, '/radshare/EPAD/Carole/EPAD_All_CaroleOutput'); % Copy WMHs segmented by Carole, making sure to also copy the FLAIRs that are in alignment with the WMH

% Check availability files
EPAD_ReportMissingFiles(ROOT, false); % set latter to true for removing incomplete subjects for re-import

%% ========================================================
%% Here we should copy the analysis folder for storing in XNAT!

if bRunExploreASL

    % 1) FOR ALL NormPE/RevPE -> KEEP ONLY FIRST, REMOVE OTHER TIMEPOINTS
    % 4) MOVE T1 & FLAIR TO SUBFOLDERS PER BIDS
    % SubjList = xASL_adm_GetFsList(SourceDir, '^\d{3}-\d{4}$', 1);
    % for iL=1:length(SubjList)
    %     SubDir = fullfile(SourceDir, SubjList{iL});
    %     AnatDir = fullfile(SubDir, 'anat');
    %     xASL_delete(fullfile(SubDir, 'FLAIR_parms.mat'));
    %     xASL_delete(fullfile(SubDir, 'T1_parms.mat'));
    %     xASL_adm_CreateDir(AnatDir);
    %     FileList = {'FLAIR.json' 'FLAIR.nii' 'T1.json' 'T1.nii'};
    %     for iI=1:length(FileList)
    %         OriPath = fullfile(SubDir, FileList{iI});
    %         DestPath = fullfile(AnatDir, FileList{iI});
    %         if xASL_exist(OriPath, 'file')
    %             xASL_Move(OriPath, DestPath, 1);
    %         end
    %     end
    % end
    %  ZIP ALL FILES
    %  xASL_adm_GzipAllFiles(SourceDir);



    %% ========================================================
    %% 4) Run ExploreASL for total dataset
    %  Prepare parameter file
    if isempty(which('DataParametersEPAD_HighQuality'))
        error('ParmsFile missing!'); % checks if Matlab knows where the file is
    else
        OriPath = which('DataParametersEPAD_HighQuality'); %% Can change this into HIGH/LOW QUALITY
        DestPath = fullfile(AnalysisDir, 'DataParametersEPAD_HighQuality.m'); %% Can change this into HIGH QUALITY
        xASL_Copy(OriPath, DestPath, true);
    end

    % -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> ->
    % NOW CALL EXPLOREASL PARALLEL
    cd(fileparts(fileparts(fileparts(which('DataParametersEPAD_HighQuality'))))); % go to ExploreASL folder
    ExploreASL_Master_EPAD(DestPath, true, true);

end

end
