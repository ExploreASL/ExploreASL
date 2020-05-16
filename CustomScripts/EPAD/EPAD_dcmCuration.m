function [CurationTable] = EPAD_dcmCuration(RawDir, bUseDCMTK, bCheckPermissions)
%EPAD_dcmCuration This script manages & corrects the DICOM structure of raw data
%
% FORMAT: [CurationTable] = EPAD_dcmCuration(RawDir[, bUseDCMTK, bCheckPermissions])
%
% INPUT:
%   RawDir      - folder containing the DICOMs that need curating.
%                 folder structure should be (per Viktor's IXICO download):
%                 //raw/Site/Subject/Visit/ScanType/ScanName, making for:
%
%                 RawDir = //raw
%                 SiteDir = //raw/Site
%                 SubjectDir = //raw/Site/Subject
%                 VisitDir = //raw/Site/Subject/Visit (e.g. 'screen')
%                 ScanTypeDir = //raw/Site/Subject/Visit/ScanType
%                 ScanDir = //raw/Site/Subject/Visit/ScanType/ScanName
%
%                 Likewise, iterations/looping occurs with
%                 iSite, iSubject, iVisit, iScanType, iScan
%                 & lists: SiteList, SubjectList, VisitList, ScanTypeList, ScanList
% 
%                 EXAMPLE:
%                 Scanner - 010, 020, 022, 030, etc. Similar to site IDs from
%                      IXICO. Difference is that IXICO Sites always have a leading
%                      '0', which is changed into 1 2 3 etc for different scanners at
%                      the same site. E.g.: IXICO site 010 has scanners 010 and 110.
%            
%                 Subject - should be structured \d{3}EPAD\d* (i.e. 3 digits, EPAD, then any number of digits)
%                         the first 3 digits are the scanner ID, the last digits are not
%                         specified as this differs between IXICO and ARIDHIA
%                         IDs, allowing this script to serve both
%                 Visit - Screening, Month_12, Month_24, Month_36 etc
%                 ScanType - can be anything, but will be sorted by this
%                           function into: anat_T1w anat_FLAIR anat_T2star anat_T2w
%                           swi_part_mag func_RevPE func_NormPE func_bold dwi_RevPE
%                           dwi_dwi asl_RevPE asl_FitM0AndT1 asl_M0 asl_asl
%                           Other_PhasContrast Other_Survey Other_Reconstructions
%                           Other_MoCo Other_Localizer Other_ExamCard Other_ReferenceScan
%                           Other_B1Calibration (see ScanType_LabelsConfig.csv for
%                           up-to-date list)
%                 ScanName - can be anything. First part by IXICO will be
%                           kept (e.g. 030-00019_Screening_MRI_T2W_MR-) and the part after
%                           MR- will be replaced by the SeriesDescription/ProtocolName from
%                           the DICOM header (e.g. Axial_T2W_TSE_with_Fat_Sat_CLEAR)
%                           'Inactive' will be prepended, 'recon' postpended, 'run-2'
%                           postpended, indicating that resp. this sequence should
%                           preferably not be used, this is a reconstruction, or that this
%                           is the 2nd complete sequence of this ScanType, for this
%                           subject (located in
%                           //raw/Site/Subject/Visit/ScanType/ScanName)
%   bUseDCMTK         - true for replacing Matlab's DICOMINFO by Jan Petr's compilation of DCM Tk, which goes much faster,
%                       but may not be as extensive as DICOMINFO (i.e. exceptions need to be defined by Jan Petr) (OPTIONAL, DEFAULT=false)
%   bCheckPermissions - Check whether file permissions are set correctly.
%                       Checks each dicom individually, hence very accurate but very slow. Can be useful when wanting to run this script without SUDO (OPTIONAL, DEFAULT=false)
%                       This can be useful to notify if a script crashes because of permission issues. Data, including the ROOT folder) needs to be readable and
%                       writable for users and groups (i.e. at least 660). For ExploreASL executables/data folders we want to read, write and execute for users and groups (i.e. at least 770)
%                       xASL_adm_CheckPermissions will try resetting this, and throw a warning if it cannot.
%
% INPUT FILES:
%   ScanType_LabelsConfig.csv
% 
% OUTPUT:
%   CurationTable - This Table contains the folder structure as found and curated.
%                   It is also stored as DCM_StructureList.csv & DCM_StructureList.mat files
%                   The structure of this list is slightly different from the folder structure:
%                   CurationTable{Sites}{ScanTypes}{1+SubjectVisit,6+Sequence}
%                   the first row contains the headers with the names of the
%                   columns, and the first five columns are not-sequence specific
%                   but subject & ScanType specific:
%                   EXAMPLE:
%                   SubjectName	    VisitName        nFoundDuplicates nDeletedDuplicates nScansFound nScansUnknownName nSlices/volumes       nSlices/volumes nSlices/volumes
%                   022EPAD00004	Screening        0                0                  1           []            '170_Complete_run_1'  n/a             n/a
%                   This shows that for subject 022EPAD00004 no duplicate sequences were found, no duplicate sequences were deleted, 1 true scan was
%                   found, no scans with unknown ScanTypeID, the first sequence had DICOM files, and there were no DICOMs for the 2nd or 3rd sequences
%
% OUTPUT FILES:
% //RawDir/DICOM_Sort_Log_yyyymmdd_HHMMSS.txt - contains the log (screendump) of this function call
% //RawDir/DifferentStudyUIDList.json         - Table of subjects/visits with different StudyInstanceUIDs
% //RawDir/DeleteActiveList.json              - Table of deleted active duplicate scans where the inactive duplicate scans were more complete
% //RawDir/ActiveMissingScan.json             - Table of inactive complete scans
% //RawDir/DCM_StructureList.(csv|mat)        - export of CurationTable (described above)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: 
% This script performs an initial DICOM QC, for DICOM sorting, curation,
% deletion of duplicates etc. etc. to ensure that a complete, unambiguous
% and correctly named datasets arrive on XNAT and to the analysis
% pipelines. This function is based on Matlab, SPM, ExploreASL.
% Throughout this script, the CurationTable, ScanList etc is kept for
% reference, and updated according to the changes in the foldernames.
% The foldernames should always be leading (as opposed to the CurationTable
% names)
% It will run the following steps:
% 0) Admin
%    Here we manage input arguments, start the diary log, load the
%    ScanType_LabelsConfig.csv and manage ranges for number of slices that
%    we call complete (more info in subfunction header):
%    - PrepareSliceRange - Specify the number of slices that we expect, for this script to call a scan/sequence "Complete"
% 1) Renaming Subjects
%    Here we make sure that the subject IDs have site as a prefix, and 'EPAD'
%    instead of '-'. The IXICO & ARIDHIA IDs have different numbers of digits,
%    which we can ignore here (simply by using \d* instead of specifying the
%    number of digits)
% 2) - ReorderFolders First part of the curation/sorting is to reorder the
%      folder structure (info in subfunction header)
% 3) Iterate over ScanTypes & manage them (info in subfunction headers):
%    - GetScanInfo Obtain and parse scan information from the DICOM header 
%    - CorrectDICOMFolderDepth Move DICOMs to correct folder layer (i.e. //RawDir/SiteDir/SubjectDir/VisitDir/ScanTypeDir/ScanDir/DICOMs)
%    - GetScanInfo Obtain and parse scan information from the DICOM header
%    - AssignComplete Assign "complete" scans in the foldername (ScanDir) & the CurationTable
%    - AssignRun Assign the 'run' of each scan for scan repetitions within MRI examination (i.e. run-1 run-2 etc)
%    - AssignActiveMissing Assigns an inactive complete scan as active if active complete scan is missing
%    - VerifyIdenticalStudyUID Check that all scans belong to the same MRI examination visit
% 4) Store the CurationTable & other tables/lists in CSV/JSON files
%
% Note that for each processed site, we store the empty file
% 'Processed.status' in the VisitDir. This will skip this Subject for the next
% time this script will be run (e.g. if the server crashed & we need to
% repeat this script). Please delete these files if you wish to rerun this
% script!!!
%
% The following variables are declared in this function:
% ActiveMissingList     - cell structure containing names of scans where an active complete was missing, but inactive complete existed & was used instead
% bUseDCMTK             - true for replacing Matlab's DICOMINFO by Jan Petr's compilation of DCM Tk
% CurationTable         - This Table contains the folder structure as found and curated. It is also stored as DCM_StructureList.csv & DCM_StructureList.mat files
% CurrentSubject        - name of current subject, e.g. '010EPAD00001' (SubjList{iS})
% dcmList               - list of paths to DICOM files we want to put in correct folder layer
% DeleteActiveList      - cell structure containing names of scans where an active complete was missing, but inactive complete existed & was used instead
% DicomPath             - path to the DICOM we aim to read
% DifferentStudyUIDList - list of subjects/visits with scans from different MRI examinations/visits/studies/scan sessions
% EchoTime              - matrix containing EchoTimes, structured as (iSite,iSubject,iVisit,iScanType,iScan)
% Info                  - struct containing DICOM fields and their information
% iScan                 - index for current Scan
% iScanType             - index of current ScanType
% IsExamcard            - true if DICOM is an examcard
% iSite                 - index of current site
% iSubject              - index of current Subject
% iVisit                - index of current visit
% iSubjectVisit         - combined index/counter for subject/visit rows
% IXICORegExpList ?     - regular expression list of names IXICO/sites used for ScanTypes, taken from ScanType_LabelsConfig.csv
% nSliceNrTable         - table with number of slices we require each scan to have, for each scanner, taken from ScanType_LabelsConfig.csv
% nSubjects             - number of subjects length(SubjList)
% RemovedDuplicateList  - cell structure containing names of deleted scans
% RepetitionTime        - matrix containing RepetitionTimes, structured as (iSite,iSubject,iVisit,iScanType,iScan)
% ScanDir               - the folder where this scan/sequence is located (i.e. //raw/Site/Subject/Visit/ScanType/ScanName)
% ScanList              - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
% ScanTypeDir           - the folder where the ScanType is located (i.e. //raw/Site/Subject/Visit/ScanType)
% ScanTypeName          - composite name of scantype, taken from BIDS: i.e. the group name _ contrast name (e.g. anat_T1w, dwi_dwi, dwi_RevPE)
% SeriesDate            - cell structure containing scanning date, structured as {iSite}{iSubject}{iVisit}{iScanType}{iScan}
% SeriesDescription     - DICOM field containing the name of the sequence (e.g. Axial_T2W_TSE_with_Fat_Sat_CLEAR).
% SeriesInstanceUID     - Unique Identifier of current scan/sequence, consisting of organization UID etc
% SeriesTime            - cell structure containing scanning time, structured as {iSite}{iSubject}{iVisit}{iScanType}{iScan}
% SiteDir               - path to SiteFolder, e.g. '/EPAD500_new/010'
% SiteNrString          - string containing the number of the site/scanner
% SliceRange            - cell structure containing number of slices we expect to call a scan/sequence complete, structured as SliceRange{iSite}{iScanType}
% StudyInstanceUID      - Unique Identifier of current MRI examination/visit/study/scan session, consisting of organization UID etc
% SubjectList           - list of subjects
% UIDWarningReported    - true if a warning was thrown already for this subject/visit, to avoid throwing multiple warnings for the same MRI examination/visit
% VisitDir              - path to VisitFolder, e.g. '/EPAD500_new/010/010EPAD00001/screening'
%
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_dcmCuration('/data/RAD/share/EPAD500_new/raw');
% __________________________________
% Copyright 2017-2019 Amsterdam UMC
% Luigi Lorenzini, Silvia Ingala, Viktor Wottschel, Alle Meije Wink, Joost Kuijer, Henk-Jan Mutsaerts


%% Admin

%% Manage input arguments
if nargin<2 || isempty(bUseDCMTK)
	bUseDCMTK = false;
end

if nargin<3 || isempty(bCheckPermissions)
    if isunix || ismac
        bCheckPermissions = true;
    else
        bCheckPermissions = false; % with windows we usually have less permission issues
    end
end

if nargin<1 || isempty(RawDir)
    error('Please provide the input path to RawDir');
end
if bCheckPermissions
    xASL_adm_CheckPermissions(RawDir, false);
end


%% Start diary
DiaryPath = fullfile(RawDir, ['DICOM_Sort_Log_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
diary(DiaryPath);

%% Admin: here specify the specifics of this study
xASL_adm_CreateDir(RawDir);

% SiteName     = {'010GEHDxT'   '010SIPrisma' '011SiPrisma' '012SiVerio' '015SiPrisma'  '020PhAchieva' '021Si' '022PhIngenia' '024SiArea'  '030PhIngenia' '031SiTrio' '032SiVida'  '040PhIngenuity' '050PhIngenia' '060Skyra'};

%% Load the ScanType Labels Configuration
[~, ScanTypeConfig] = xASL_bids_csvRead(fullfile(RawDir, 'ScanType_LabelsConfig.csv'));

SiteNrString = ScanTypeConfig(1,5:end); % NB: THIS ASSUMES THAT VENDORS START AT THE 5TH COLUMN
SiteNrString = cellfun(@(x) x(5:end), SiteNrString, 'UniformOutput',false);

IXICORegExpList = ScanTypeConfig(2:end,1); % what regexp names can the scans have
%  Make sure that the order is such, that the most specific scan types are first:
%  e.g. FLAIR before T2, such that if a FLAIR is found, it doesn't specify this scan also as a T2
for iB=2:size(ScanTypeConfig,1)
    ScanTypeName{iB-1,1} = [ScanTypeConfig{iB,2} '_' ScanTypeConfig{iB,3}];
end
nSliceNrTable = ScanTypeConfig(:,5:end); % HERE WE ASSUME THAT THE 5TH COLUMN CONTAINS THE FIRST SITE

%% Load the number of slices that we expect, for the script to call them "Complete"
[SliceRange] = PrepareSliceRange(SiteNrString, ScanTypeName, nSliceNrTable, IXICORegExpList);


%% Initiate lists
DifferentStudyUIDList = struct;
ActiveMissingList = struct;
StudyInstanceUID{1}{1}{1}{1} = '';
SeriesInstanceUID{1}{1}{1}{1} = '';
SeriesDescription{1}{1}{1}{1} = '';
SeriesDate{1}{1}{1}{1} = [];
SeriesTime{1}{1}{1}{1} = [];
EchoTime = 0;
RepetitionTime = 0;
RemovedDuplicateList = {}; % removed duplicates
DeleteActiveList = struct; % these were removed while they were active, keeping the inactive (because the inactive was more complete)


%% 1) We rename the folders here, such that the SubjectNames contain the Site in a proper way
% Here we make sure that the subject IDs have site as a prefix, and 'EPAD'
% instead of '-'. The IXICO & ARIDHIA IDs have different numbers of digits,
% which we can ignore here (simply by using \d* instead of specifying the
% number of digits)

fprintf('%s\n','Renaming subjects to have site as prefix');
for iSite=1:length(SiteNrString)
    SiteDir = fullfile(RawDir, SiteNrString{iSite});
    if exist(SiteDir, 'dir')
        % For ARIDHIA code (originally had subject ID only)
        SubjectList = xASL_adm_GetFsList(SiteDir,'^\d*$', 1, [], [], [0 Inf]);
        for iSubject=1:length(SubjectList)
            xASL_TrackProgress(iSubject, length(SubjectList));
            OldDir = fullfile(RawDir, SiteNrString{iSite}, SubjectList{iSubject});
            NewDir = fullfile(RawDir, SiteNrString{iSite}, [SiteNrString{iSite} 'EPAD' SubjectList{iSubject}]);
            xASL_MoveDirMerge(OldDir, NewDir, true);
        end
        % For IXICO code (originally had site-subjectID)
        SubjectList = xASL_adm_GetFsList(SiteDir,['^' SiteNrString{iSite} '-\d*$'], true, [], [], [0 Inf]);
        for iSubject=1:length(SubjectList)
            xASL_TrackProgress(iSubject, length(SubjectList));            
            OldDir = fullfile(RawDir, SiteNrString{iSite}, SubjectList{iSubject});
            NewDir = fullfile(RawDir, SiteNrString{iSite}, [SiteNrString{iSite} 'EPAD' SubjectList{iSubject}(5:end)]);
            xASL_MoveDirMerge(OldDir, NewDir, true);
        end        
    end
end



%% Start running script
for iSite=1:length(SiteNrString) % iSite loops over sites (i.e. scanners)
    
    SiteDir = fullfile(RawDir, SiteNrString{iSite});
    
    if ~exist(SiteDir, 'dir')
        fprintf('%s\n',['Scanner ' SiteNrString{iSite} ' didnt exist, skipping...']);
        continue;
    end

    fprintf('%s',['Processing scanner ' SiteNrString{iSite} ':   '])
    SubjectList = xASL_adm_GetFsList(SiteDir,'^\d{3}EPAD\d*$', true, [], [], [0 Inf]);

    % Now we process it
    for iScanType=1:length(ScanTypeName) % Initialize the CurationTable
        CurationTable{iSite}{iScanType}(1,1:6) = {'SubjectName' 'VisitName' 'nFoundDuplicates' 'nDeletedDuplicates' 'nScansFound' 'nScansUnknownName'};
    end

    for iSubject=1:length(SubjectList) % iSubject loops over subjects within an iSite
        UIDWarningReported = false;
        SubjectDir = fullfile(SiteDir, SubjectList{iSubject});
        
        VisitList = xASL_adm_GetFileList(SubjectDir,'^(Screening|Month_\d*)$','List',[0 Inf],true);
        if length(VisitList)>1
            TempVisitList=VisitList;
            VisitList{1} = TempVisitList{end};
            VisitList(2:length(TempVisitList)) = TempVisitList(1:end-1);
            clear TempVisitList;
        end
        for iVisit=1:length(VisitList)
            VisitDir = fullfile(SubjectDir, VisitList{iVisit});

            PathStatus = fullfile(VisitDir,'Processed.status');
            if exist(PathStatus, 'file')
                fprintf('%s\n',[VisitDir ' was already processed, skipping...']);
                continue;
            end            
            
            %% 2) First part of the sorting: reorder folders
           ReorderFolders(VisitDir, IXICORegExpList, ScanTypeName, length(SubjectList), iSubject, bUseDCMTK);

           iSubjectVisit = (iSubject-1)*length(VisitDir)+iVisit;
           
            %% 3) Second part of the sorting: manage ScanTypes
            %  Here we will delete them, check them, etc & add them to the list
            for iScanType=1:length(ScanTypeName) % iM loops over "modalities"
                TotalCount = length(SubjectList)*length(VisitList)+length(ScanTypeName);
                CurrentCount = (iSubject-1)*length(VisitList)*length(ScanTypeName)+iVisit*length(ScanTypeName)+iScanType;
                xASL_TrackProgress(CurrentCount, TotalCount);                    

                ScanTypeDir = fullfile(VisitDir, ScanTypeName{iScanType});

                % Count number of Scans
                if ~exist(ScanTypeDir,'dir')
                    ScanList = [];
                else
                    ScanList = xASL_adm_GetFsList(ScanTypeDir,'.*.$',true,[],[],[0 Inf]);
                end

                % Check existing scans
                for iK=7:9
                    CurationTable{iSite}{iScanType}{   1,iK} = 'nSlices/volumes';
                    CurationTable{iSite}{iScanType}{iSubjectVisit+1,iK} = 'n/a';
                end

                % First remove empty dirs
                for iScan=1:length(ScanList)
                    ScanDir = fullfile(ScanTypeDir, ScanList{iScan});
                    dcmList = xASL_adm_GetFileList(ScanDir,'^.*[^.(csv|tsv|gz)]$','FPListRec',[0 Inf]);
                    if isempty(dcmList)
                        rmdir(ScanDir, 's');
                    end
                end

                % Now update the ScanList
                if exist(ScanTypeDir,'dir')
                    ScanList = xASL_adm_GetFsList(ScanTypeDir,'.*.$',true,[],[],[0 Inf]);
                end

                for iScan=1:length(ScanList)
                    ScanDir = fullfile(ScanTypeDir, ScanList{iScan});            

                    %% Here we obtain for each individual scan its SeriesInstanceUID, StudyInstanceUID, and 
                    %  number of slices (the latter goes into the List)
                    [CurationTable, ScanList, StudyInstanceUID, SeriesInstanceUID, dcmList, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime] = GetScanInfo(CurationTable, ScanList, ScanDir, StudyInstanceUID, SeriesInstanceUID, SiteNrString, iSite, iSubject, iVisit, iSubjectVisit, iScanType, iScan, bUseDCMTK, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime);

                    % Move DICOMs to correct folder layer, if they are "too deep"
                    CorrectDICOMFolderDepth(dcmList, ScanDir);
                end

                nRemovedDuplicates = 0; % default
                if ~isempty(ScanList)
                    %% Search for true duplicate scans & delete them
                    % NB: in this submodule we also change the lists based on deleted scans, this could induce errors!
                    PreviousNDuplicates = length(RemovedDuplicateList);
                    [CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, RemovedDuplicateList, DeleteActiveList] = RemoveDuplicates(CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, iSite, iSubject, iVisit, iSubjectVisit, iScanType, RemovedDuplicateList, ScanTypeDir, DeleteActiveList);

                    %% Assign "complete" scans in the list & to the foldername
                    % If the scan has correct number of slices, call it 'Complete'
                    [CurationTable, ScanList] = AssignComplete(ScanList, CurationTable, ScanTypeName, ScanTypeDir, SliceRange, iSite, iScanType, iSubjectVisit, SiteNrString);

                    %% Assign the 'run' for complete scans (preferably actives)
                    [CurationTable, ScanList] = AssignRun(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubject, iVisit, iSubjectVisit, SeriesTime);

                    %% If we only have an inactive scan, but it is complete, then replace "inactive" by "ActiveMissing"
                    %  this way, it will be used by the dcm2nii import
                    [ScanList, ActiveMissingList] = AssignActiveMissing(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubjectVisit, ActiveMissingList);

                    nRemovedDuplicates = length(RemovedDuplicateList)-PreviousNDuplicates;
                end

                % Save in Table
                CurationTable{iSite}{iScanType}(iSubjectVisit+1,1:5) = {SubjectList{iSubject} VisitList{iVisit} nRemovedDuplicates nRemovedDuplicates length(ScanList)}; % SubjectName VisitName nRemovedDuplicates nRemovedDuplicates number of scans found
            end % loop over ScanTypes

            %% 11) Check for identical StudyInstanceUIDs:
            % Check that all scans belong to the same subject/visit, otherwise throw warning
            [UIDWarningReported, DifferentStudyUIDList] = VerifyIdenticalStudyUID(UIDWarningReported, StudyInstanceUID, iSite, iSubject, SubjectList, iVisit, VisitList, DifferentStudyUIDList);

            fclose(fopen(PathStatus, 'wt')); % create processed file
        end % loop over visits
    end % loop over subjects
    fprintf('\n');
end % loop over sites


%% Here we store the lists in CSV file
if exist('CurationTable','var')
    SavePath = fullfile(RawDir, 'DCM_StructureList.mat');
    save(SavePath, 'CurationTable');

    SavePath = fullfile(RawDir, 'DCM_StructureList.csv');
    fclose all;
    xASL_delete(SavePath);
    FID = fopen(SavePath, 'wt');
    % Print header
    fprintf(FID, '%s\n', 'ScanType,Subject,Visit,nFoundDuplicates,nDeletedDuplicates,nScansFound,nScansUnknownName,nSlices/volumes,nSlices/volumes,nSlices/volumes');

    % Print all the lists
    for iSite=1:length(CurationTable)
        fprintf(FID, '\n%s', SiteNrString{iSite});
        for iScanType=1:length(CurationTable{iSite})
            if size(CurationTable{iSite}{iScanType},1)>1
                fprintf(FID, '\n,');
                % print cell content:
                for iX=2:size(CurationTable{iSite}{iScanType},1) % omit the first row
                    fprintf(FID, '\n%s,', ScanTypeName{iScanType});
                    for iY=1:size(CurationTable{iSite}{iScanType},2)
                        PrintNow = CurationTable{iSite}{iScanType}{iX,iY};
                        if isnumeric(PrintNow)
                            fprintf(FID, '%s,', num2str(PrintNow));
                        elseif strcmp(PrintNow,'n/a')
                            fprintf(FID, '%s,', '0');
                        else
                            fprintf(FID, '%s,', PrintNow);
                        end
                    end
                end
            end
        end
        fprintf(FID, '\n');
    end
    fclose(FID);
end
diary('off');

% Store other cell structures/lists
if length(fields(DifferentStudyUIDList))>0
    warning('Some subjects had scans with different Study IDs, check this in the list, please check ./DifferentStudyUIDList.json');
    SavePath = fullfile(RawDir, 'DifferentStudyUIDList.json');
    xASL_adm_SaveJSON(DifferentStudyUIDList, SavePath);
end
if length(fields(DeleteActiveList))>0
    warning('We deleted active scans as the inactive scans were more complete, please check ./DeleteActiveList.json');
    SavePath = fullfile(RawDir, 'DeleteActiveList.json');
    xASL_adm_SaveJSON(DeleteActiveList, SavePath);
end
if length(fields(ActiveMissingList))>0
    warning('Some subjects only had inactive complete scans, check this in the list, please check ./ActiveMissingList.json');
    SavePath = fullfile(RawDir, 'ActiveMissingList.json');
    xASL_adm_SaveJSON(ActiveMissingList, SavePath);
end

end


















%% SUBFUNCTIONS
                
%% =============================================================================================================================                
%% =============================================================================================================================
function ReorderFolders(VisitDir, IXICORegExpList, ScanTypeName, nSubjects, iSubject, bUseDCMTK)
%ReorderFolders First part of the curation/sorting is to reorder the folder structure
%
% FORMAT: ReorderFolders(VisitDir, IXICORegExpList, ScanTypeName, nSubjects, iSubject, bUseDCMTK)
% 
% INPUT:
%   VisitDir          - path to VisitFolder, e.g. '/EPAD500_new/010/010EPAD00001/screening'
%   IXICORegExpList   - regular expression list of names IXICO/sites used for ScanTypes, taken from ScanType_LabelsConfig.csv
%   ScanTypeName      - composite name of scantypes, taken from BIDS: i.e. the group name _ contrast name
%                       (e.g. anat_T1w, dwi_dwi, dwi_RevPE)
%   nSubjects         - number of subjects length(SubjList)
%   iSubject          - index of current subject
%   bUseDCMTK         - true for replacing Matlab's DICOMINFO by Jan Petr's compilation of DCM Tk, which goes much faster,
%                       but may not be as extensive as DICOMINFO (i.e. exceptions need to be defined by Jan Petr) (OPTIONAL, DEFAULT=false)
%
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: First part of the sorting: reorder folders:
% a) Manage Scan-prefixes & suffixes: here we prefix "inactive" if the scan was labeled as inactive (they are moved to their respective ScanType
%    folder, and the inactive folder is deleted). Also IXICO's suffixes 'null' 'examcard' and 'version' are removed
% b) File management (unzipping, fix dcm extensions, remove dicom0 files)
% c) xASL_io_ReadTheDicom & GetSeriesDescription: reads a DICOM within the ScanDir
%    & extracts its SeriesDescription/ProtocolName
% d) RenameScanFolderWithSeriesDescription: Rename scan folders using SeriesDescription
% e) Deal with ExamCards: change DICOM extensions from .dcm into .ExamCard
     % change ScanDir suffix into '_MR-ExamCard'
% e) MoveScanToScantypeFolder: Move scans to the correct ScanType folder
% f) Move sequence to correct ScanTypeFolder
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: ReorderFolders(VisitDir, IXICORegExpList, ScanTypeName, nSubjects, iSubject, bUseDCMTK)
% __________________________________
% Copyright 2015-2019 ExploreASL


% Create list of found scans within all subdirs (e.g. "inactive", "other", & scan names)
ScanTypeList = xASL_adm_GetFsList(VisitDir,'.*', true, [], [], [0 Inf]);

% First move scans to the correct ScanType folder
for iScanType=1:length(ScanTypeList)

    ScanTypeDir = fullfile(VisitDir, ScanTypeList{iScanType});
    ScanList = xASL_adm_GetFsList(ScanTypeDir,'.*.$', true, [], [], [0 Inf])';

    for iScan=1:length(ScanList)
%         TotalCount = nSubjects*length(ScanTypeList)*length(ScanList);
%         CurrentCount = (iSubject-1)*length(ScanTypeList)*length(ScanList)+(iScanType-1)*length(ScanList)+iScan;
%         xASL_TrackProgress(CurrentCount, TotalCount);

        %% a) Manage Scan-prefixes & suffixes
        ScanDir = fullfile(ScanTypeDir,ScanList{iScan});

        % If inactive folder, prefix 'Inactive_'
        if strcmp(lower(ScanTypeList{iScanType}),'inactive') && exist(ScanDir,'dir')
            ScanList{iScan} = RemoveExistingRegExp(ScanList{iScan}, 'inactive(_|)');
            ScanList{iScan} = ['Inactive_' ScanList{iScan}];
            NewDir = fullfile(ScanTypeDir, ScanList{iScan});
            if ~strcmp(ScanDir, NewDir)
                while exist(NewDir,'dir')
                    NewDir = [NewDir '_2'];
                end
                xASL_Move(ScanDir, NewDir);
                ScanDir = NewDir;
            end
        end

        % Remove suffixes
        Suffexes = {'version' 'null' 'examcard'};
        for iSuf=1:length(Suffexes)
            if ~isempty(strfind(lower(ScanList{iScan}),Suffexes{iSuf})) && exist(ScanDir,'dir')
                ScanList{iScan} = RemoveExistingRegExp(ScanList{iScan}, [Suffexes{iSuf} '\d*']);
                NewDir = fullfile(ScanTypeDir, ScanList{iScan});
                if ~strcmp(ScanDir, NewDir) % if the new directory differs from the old
                    while exist(NewDir,'dir')
                        NewDir = [NewDir '_2'];
                    end
                    xASL_Move(ScanDir, NewDir);
                    ScanDir = NewDir;
                end
            end
        end
        
        %% b) File management (unzipping, fix dcm extensions, remove dicom0 files)
        % Unzip dicoms
        if ~isempty(xASL_adm_GetFileList(ScanDir, '^.*\.dcm\.gz$', 'FPList', [0 Inf]))
            gunzip(fullfile(ScanDir, '*.gz'));
            xASL_adm_DeleteFileList(ScanDir, '^.*\.dcm\.gz$',0,[0 Inf]);
        end
        % Fixing dcm suffixes with non-.dcm extension (e.g. IMG (GE), IMA (Siemens), empty extension)
        SuffList = xASL_adm_GetFileList(ScanDir, '(.*\.img|.*\.dic|.*\.IMA|^([^.]*))$', 'FPList', [0 Inf]);
        for iSuff=1:length(SuffList)
            [Fpath, Ffile] = fileparts(SuffList{iSuff});
            xASL_Move(SuffList{iSuff}, fullfile(Fpath, [Ffile '.dcm']));
        end
        % Removing odd dicom0 files
        if ~isempty(xASL_adm_GetFileList(ScanDir, '^.*\.dcm0$', 'FPList', [0 Inf]))
            xASL_adm_DeleteFileList(ScanDir, '^.*\.dcm0$',0,[0 Inf]);
        end        
        
        dcmList = xASL_adm_GetFileList(ScanDir, '^.*[^.(csv|tsv|gz)]$', 'FPListRec', [0 Inf]); % \.dcm|.*\.img|.*\.IMA|
        if ~isempty(dcmList)
            %% c) xASL_io_ReadTheDicom & GetSeriesDescription
            % First we check the folder name given by EPAD:
            IndexDCM = ceil(0.5*length(dcmList));
            Info = xASL_io_ReadTheDicom(bUseDCMTK, dcmList{IndexDCM});
            [SeriesDescription, IsExamcard] = GetSeriesDescription(Info);

            %% d) Rename ScanFolder to SeriesDescription, if needed
            if ~isempty(SeriesDescription)
                ScanList{iScan} = RenameScanFolderWithSeriesDescription(ScanList{iScan}, SeriesDescription, ScanDir, IXICORegExpList, ScanTypeDir);
            end
            
            %% e) Deal with ExamCards
            if IsExamcard && length(dcmList)==1
                % a) change file extension from dcm into ExamCard
                [Dpath, Dfile] = xASL_fileparts(dcmList{1});
                NewFile = fullfile(Dpath, [Dfile '.ExamCard']);
                if ~strcmp(dcmList{1}, NewFile)
                    xASL_Move(dcmList{1}, NewFile);
                end
                
                % b) add "ExamCard" to folder name
                [Dpath, Dfile] = fileparts(Dpath);
                Dfile = xASL_adm_CorrectName(Dfile, [], '*');
                if strcmp(Dfile(end-2:end),'MR_') || strcmp(Dfile(end-2:end),'MR-')
                    Dfile = [Dfile(1:end-3) 'MR-ExamCard'];
                elseif strcmp(Dfile(end-2:end),'_MR') || strcmp(Dfile(end-2:end),'-MR')
                    Dfile = [Dfile(1:end-3) '_MR-ExamCard'];
                else
                    Dfile = [Dfile '_ExamCard'];
                end
                
                xASL_Move(fileparts(dcmList{1}), fullfile(Dpath, Dfile));
                ScanList{iScan} = Dfile;
            end
            
            %% f) Move sequence to correct ScanTypeFolder
            % Get the AcquisitionName from the FolderName (FolderName was adapted above)
            MoveScanToScantypeFolder(ScanList{iScan}, IXICORegExpList, ScanTypeName, ScanTypeDir, VisitDir);            
        end
    end
    %% Remove ScanTypeDir if isempty
    if isempty(xASL_adm_GetFileList(ScanTypeDir, '.*', 'FPListRec', [0 Inf], false))
        xASL_delete(ScanTypeDir);
    end
end

fprintf('\n');

end



%% =============================================================================================================================                
%% =============================================================================================================================
function [CurationTable, ScanList] = AssignComplete(ScanList, CurationTable, ScanTypeName, ScanTypeDir, SliceRange, iSite, iScanType, iSubjectVisit, SiteNrString)
%AssignComplete Assign "complete" scans in the foldername (ScanDir) & the CurationTable
%
% FORMAT: [CurationTable, ScanList] = AssignComplete(ScanList, CurationTable, ScanTypeName, ScanTypeDir, SliceRange, iSite, iScanType, iSubjectVisit, SiteNrString)
% 
% INPUT:
%   ScanList            - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%                         this list is adapted/updated with the changed folder names
%   CurationTable       - This Table contains the folder structure as found and curated.
%                         It is also stored as DCM_StructureList.csv & DCM_StructureList.mat files
%                         The structure of this list is slightly different from the folder structure:
%                         CurationTable{Sites}{ScanTypes}{1+Subject,5+Sequence}
%   ScanTypeName        - name of the scantype (e.g. anat_FLAIR, dwi_dwi, dwi_RevPE)
%   ScanTypeDir         - the folder where the ScanType is located (i.e. //raw/Site/Subject/Visit/ScanType)
%   SliceRange          - cell structure containing number of slices we
%                         expect to call a scan/sequence complete, structured as SliceRange{iSite}{iScanType}
%   iSite               - index of current site
%   iScanType           - index of current ScanType
%   iSubjectVisit       - combined index/counter for subject/visit rows
%   SiteNrString        - string containing the number of the site/scanner
% 
% OUTPUT:
%   CurationTable       - as input
%   ScanList            - as input
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function adds the suffix "Complete" to a ScanDir (and likewise in the CurationTable) when a ScanDir contains the number of DICOMs
%              that lies within the SliceRange
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [CurationTable, ScanList] = AssignComplete(ScanList, CurationTable, ScanTypeName, ScanTypeDir, SliceRange, iSite, iScanType, iSubjectVisit, SiteNrString)
% __________________________________
% Copyright 2015-2019 ExploreASL

for iScan=1:length(ScanList)
    IsComplete = false;
    CurrentScan = xASL_adm_CatchNumbersFromString(CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+iScan});

    % Here we need a numeric value (i.e. nSlices) & a valid slice range
    if isempty(SliceRange{iSite}{iScanType})
        warning(['SliceRange unknown for site ' SiteNrString{iSite} ' ScanType ' ScanTypeName{iScanType}]);
    elseif ~isempty(CurrentScan) % we allow for the exact expected number of slices, or a few more
        if  max(CurrentScan==SliceRange{iSite}{iScanType}) % if sufficient nSlices
            IsComplete = true; % set to complete
        end
    end

    if IsComplete
        TempList = RemoveExistingRegExp(CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+iScan}, '_complete'); % remove the "complete" without capital
        CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+iScan} = [TempList '_Complete'];
        
        ScanDir = fullfile(ScanTypeDir, ScanList{iScan});
        ScanList{iScan} = RemoveExistingRegExp(ScanList{iScan}, '_complete'); % remove the "complete" without capital
        ScanList{iScan} = [ScanList{iScan} '_Complete']; % add "complete" with capital
        NewPath = fullfile(ScanTypeDir, ScanList{iScan});

        if ~strcmp(ScanDir, NewPath)
            while exist(NewPath,'dir')
                NewPath = [NewPath '_2'];
            end            
            xASL_Move(ScanDir, NewPath);
        end
    end
end


end



%% =============================================================================================================================                
%% =============================================================================================================================
function [SliceRange] = PrepareSliceRange(SiteNrString, ScanTypeName, nSliceNrTable, IXICORegExpList)
%PrepareSliceRange Specify the number of slices that we expect, for this script to call a scan/sequence "Complete"
%
% FORMAT: [SliceRange] = PrepareSliceRange(SiteNrString, ScanTypeName, nSliceNrTable, IXICORegExpList)
% 
% INPUT:
%   SiteNrString     - string containing the number of the site/scanner
%   ScanTypeName     - name of the scantype (e.g. anat_FLAIR, dwi_dwi, dwi_RevPE)
%   IXICORegExpList  - regular expression list of names IXICO/sites used for ScanTypes, taken from ScanType_LabelsConfig.csv
%   nSliceNrTable    - table with number of slices we require each scan to have, for each scanner, taken from ScanType_LabelsConfig.csv
% 
% OUTPUT:
%   SliceRange       - cell structure containing number of slices we
%                      expect to call a scan/sequence complete, structured as SliceRange{iSite}{iScanType}
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads the number of slices that we expect, for the script to call them "Complete"
%              This can be more than 1 number. Will be called 'NaN' if we don't know it (or if it is not necessary)
%              we allow for 3 more or less slices than the
%              specified number (e.g. 170 becomes a range of 167-173 slices). This
%              speeds up the inspection of incomplete scans in case of additional
%              descriptive DICOMs, or the case where 1 or 2 slices are missing, which may not be too bad.
%              The significance of having a few slices more or less than the
%              specified protocol range becomes apparent in the NIfTI QC stage.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [SliceRange] = PrepareSliceRange(SiteNrString, ScanTypeName, nSliceNrTable, IXICORegExpList)
% __________________________________
% Copyright 2015-2019 ExploreASL

% 1) Create empty Table for all sites:
for iY=1:length(SiteNrString)
    for iX=1:length(ScanTypeName)
        nExpectedSlices{iY,iX} = NaN;
    end
end

% 2) Fill with known numbers
for iX=1:size(nSliceNrTable,2)
    % find the site number
    for iSite=1:length(SiteNrString)
        if ~isempty(findstr(SiteNrString{iSite}, nSliceNrTable{1,iX}))
            CurrSiteN = iSite;
        end
    end
    for iY=2:size(nSliceNrTable,1)
        nExpectedSlices{iY-1, CurrSiteN} = sort(xASL_adm_CatchNumbersFromString(nSliceNrTable{iY,iX}));
    end
end   
        
% 3) Create slice ranges to validate complete scans
RR  = 3; % allowed buffer/range above/below
for iSite=1:length(SiteNrString)
    for iScanType=1:length(IXICORegExpList)
        for iN=1:length(nExpectedSlices{iScanType,iSite})
            if ~isfinite(nExpectedSlices{iScanType,iSite}(iN))
                SliceRange{iSite}{iScanType} = NaN;
            elseif nExpectedSlices{iScanType,iSite}(iN)==0
                SliceRange{iSite}{iScanType} = 0;
            else
                
                RangeIs = [nExpectedSlices{iScanType,iSite}(iN)-RR:nExpectedSlices{iScanType,iSite}(iN)+RR];
                if iN==1
                    SliceRange{iSite}{iScanType} = RangeIs;
                else 
                    SliceRange{iSite}{iScanType}(end+1:end+length(RangeIs)) = RangeIs;
                end
            end
        end
    end
end


end



%% =============================================================================================================================                
%% =============================================================================================================================
function [ScanList] = RenameScanFolderWithSeriesDescription(ScanList, SeriesDescription, ScanDir, IXICORegExpList, ScanTypeDir)
%RenameScanFolderWithSeriesDescription Rename sequence-folder to SeriesDescription, if needed
%
% FORMAT: [ScanList] = RenameScanFolderWithSeriesDescription(ScanList, SeriesDescription, ScanDir, IXICORegExpList, ScanTypeDir)
% 
% INPUT:
%   ScanList            - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%                         this list is adapted/updated with the changed folder names
%   SeriesDescription   - DICOM field containing the name of the sequence (e.g. Axial_T2W_TSE_with_Fat_Sat_CLEAR).
%                         If DICOM field SeriesDescription doesnt exist, it is taken from DICOM field ProtocolName.
%                         SeriesDescription is used to identify the ScanType (e.g. the word T2W in this example)
%   ScanDir             - the folder where this scan/sequence is located (i.e. //raw/Site/Subject/Visit/ScanType/ScanName)
%   IXICORegExpList     - regular expression list of names IXICO/sites used for ScanTypes, taken from ScanType_LabelsConfig.csv
%   ScanTypeDir         - the folder where the ScanType is located (i.e. //raw/Site/Subject/Visit/ScanType)
% 
% OUTPUT:
%   ScanList            - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%                         this list is adapted/updated with the changed folder names
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function concerns the name of the ScanDir after the part 'MR-', where IXICO printed the SeriesDescription (sometimes
% incorrectly). The part before this is left untouched. If the part after 'MR-' doesnot contain a recognizable ScanType identification (e.g. T2W),
% then we replace it by the contents of the SeriesDescription 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ScanList] = RenameScanFolderWithSeriesDescription(ScanList, SeriesDescription, ScanDir, IXICORegExpList, ScanTypeDir)
% __________________________________
% Copyright 2015-2019 ExploreASL


% First make sure that we ignore everything before MR- (because we will do this below as well)
% We assume that IXICO printed the SeriesDescription after MR- (anything before, we ignore)
ScanListSeriesName = strfind(ScanList,'MR-');
if isempty(ScanListSeriesName)
    ScanListSeriesName = ScanList;
elseif ScanListSeriesName==length(ScanList)-2
    ScanListSeriesName = '';
else 
    ScanListSeriesName = ScanList(ScanListSeriesName+3:end);
end


CorrAcqName = xASL_adm_CorrectName(lower(ScanListSeriesName), 2, '*'); % The asterisk '*' is left in for T2*
CorrDcmName = xASL_adm_CorrectName(lower(SeriesDescription), 2, '*');
if isempty(strfind(CorrAcqName, CorrDcmName))
    % if there is a name discrepancy between dicom field & folder
    % (the SeriesDescription should be inside the FolderName)
    % Check first whether the dicom field contains a known ScanType
    % Now search if this FolderName fits with any of the known modalities
    ScanTypeMatch = find(cellfun(@(x) ~isempty(regexp(lower(CorrDcmName),['.*' x '.*'])) && ~isempty(x), lower(IXICORegExpList)));
    
    IndN   = strfind(ScanList,'MR-');
    PreFix = ScanList(1:IndN+2); % to comply with the search criterion '-MR' below    
    Suffix = ScanList(IndN+3:end); % to comply with the search criterion '-MR' below    
    
    if ~isempty(ScanTypeMatch) && ~strcmp(lower(Suffix),lower(SeriesDescription))
        % we only replace what IXICO did by the SeriesDescription, if this makes sense (i.e. equals with what we found)
        % & if this differed from what it currently is
        
        % Then rename folder
        ScanList = xASL_adm_CorrectName(SeriesDescription, 1, '*');

        if ~isempty(PreFix)
            ScanList = [PreFix ScanList];
        end

        % if the folder already exist, add count-number suffix
        while exist(fullfile(ScanTypeDir, ScanList),'dir')
              ScanList = [ScanList '_2'];
        end
        xASL_Move(ScanDir, fullfile(ScanTypeDir, ScanList));
    end
end

end





%% =============================================================================================================================                
%% =============================================================================================================================
function MoveScanToScantypeFolder(ScanList, IXICORegExpList, ScanTypeName, ScanTypeDir, VisitDir)
%MoveScanToScantypeFolder Move scan/sequence to correct ScanTypeFolder
%
% FORMAT: MoveScanToScantypeFolder(ScanList, IXICORegExpList, ScanTypeName, ScanTypeDir, VisitDir)
% 
% INPUT:
%   ScanList            - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%                         this list is adapted/updated with the changed folder names
%   IXICORegExpList     - regular expression list of names IXICO/sites used for ScanTypes, taken from ScanType_LabelsConfig.csv
%   ScanTypeName        - name of the scantype (e.g. anat_FLAIR, dwi_dwi, dwi_RevPE)
%   SiteDir             - path to the folder where the Site is located (i.e. //raw/Site)
%   ScanTypeDir         - path to the folder where the ScanType is located (i.e. //raw/Site/Subject/Visit/ScanType)
%   VisitDir            - path to VisitFolder, e.g. '/EPAD500_new/010/010EPAD00001/screening'
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function move a scan/sequence to correct
% ScanTypeFolder. It recognizes the ScanType by the foldername, which
% should contain the ScanType identification as we ensured this in RenameScanFolderWithSeriesDescription above.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: MoveScanToScantypeFolder(ScanList, IXICORegExpList, ScanTypeName, ScanTypeDir, VisitDir)
% __________________________________
% Copyright 2015-2019 ExploreASL


IndN        = strfind(ScanList,'MR-');
IndV        = strfind(ScanList,'_version');
if  isempty(IndN)
    IndN    = 1;
end
if ~isempty(IndV)
    IndV    = IndV-1;
else
    IndV    = length(ScanList);
end
AcqName     = ScanList(IndN+3:IndV);

% Remove empty regular expressions:
NotEmptyList = cellfun(@(x) ~isempty(x), IXICORegExpList);
IXICORegExpList = IXICORegExpList(NotEmptyList);
ScanTypeName = ScanTypeName(NotEmptyList);

% Sort this sequence to the correct ScanType folder
% First search if this FolderName fits with any of the known modalities

ScanTypeMatch = find(cellfun(@(x) ~isempty(regexp(lower(AcqName),['.*' x '.*'])) && ~isempty(x), lower(IXICORegExpList)));

% If it does, put it in the ScanType folder
if ~isempty(ScanTypeMatch)
    % always take the first ScanType match (note the hierarchy specified above)
    ScanTypeMatch = ScanTypeName{ScanTypeMatch(1)};
    ScanTypeMatch = xASL_adm_CorrectName(ScanTypeMatch, [], '*');
    ScanDir = fullfile(ScanTypeDir, ScanList);
    NewDir = fullfile(VisitDir, ScanTypeMatch);
    NewScanDir = fullfile(NewDir, ScanList);

    % But only if it is not already in the correct ScanType folder
    % (Allowing Upper/Lower-case differences)
    if ~strcmp(lower(ScanTypeDir), lower(NewDir))
        % if the folder already exist, add count-number suffix
        while exist(NewScanDir,'dir')
            NewScanDir = [NewScanDir '_2'];
        end
        if ~strcmp(ScanDir, NewScanDir)
			xASL_Move(ScanDir, NewScanDir);
        end
    end
end


end



%% =============================================================================================================================
%% =============================================================================================================================
function [CurationTable, ScanList, StudyInstanceUID, SeriesInstanceUID, dcmList, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime] = GetScanInfo(CurationTable, ScanList, ScanDir, StudyInstanceUID, SeriesInstanceUID, SiteNrString, iSite, iSubject, iVisit, iSubjectVisit, iScanType, iScan, bUseDCMTK, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime)
%GetScanInfo Obtain and parse scan information from the DICOM header
%
% FORMAT: [CurationTable, ScanList, StudyInstanceUID, SeriesInstanceUID, dcmList, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime] = GetScanInfo(CurationTable, ScanList, ScanDir, StudyInstanceUID, SeriesInstanceUID, SiteNrString, iSite, iSubject, iVisit, iSubjectVisit, iScanType, iScan, bUseDCMTK, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime)

% CurationTable
% ScanList              - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%                         this list is adapted/updated with the changed folder names
% ScanDir               - path to the folder where this scan/sequence is located (i.e. //raw/Site/Subject/Visit/ScanType/ScanName)
% StudyInstanceUID      - Unique Identifier of current MRI examination/visit/study/scan session, consisting of organization UID etc
% SeriesInstanceUID     - Unique Identifier of current scan/sequence, consisting of organization UID etc
% SiteNrString          - string containing the number of the site/scanner
% iSite                 - index for current site
% iSubject              - index for current subject
% iVisit                - index for current visit
% iSubjectVisit         - combined index/counter for subject/visit rows
% iScanType             - index for current ScanType
% iScan                 - index for current Scan
% bUseDCMTK             - true for replacing Matlab's DICOMINFO by Jan Petr's compilation of DCM TK
% EchoTime              - matrix containing EchoTimes, structured as (iSite,iSubject,iScanType,iScan)
% RepetitionTime        - matrix containing RepetitionTimes, structured as (iSite,iSubject,iScanType,iScan)
% SeriesDescription     - String containing the SeriesDescription (empty if not found)
% SeriesDate            - cell structure containing scanning date, structured as {iSite}{iSubject}{iScanType}{iScan}
% SeriesTime            - cell structure containing scanning time, structured as {iSite}{iSubject}{iScanType}{iScan}
%
% OUTPUT: n/a
%   dcmList             - list of paths to DICOM files we want to put in correct folder layer
%   Otherwise same as input: CurationTable, ScanList, StudyInstanceUID, SeriesInstanceUID, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function obtains and parses scan information from the header of the
%               average DICOM image (i.e. half of the list). Then reads the DICOM
%               (xASL_io_ReadTheDicom) & obtains the EchoTime, RepetitionTime, SeriesDate/AcquisitionDate, SeriesTime/AcquisitionTime
%               Then it checks for field ImageType, which can show that a scan is a
%               reconstruction by the words derived|projection, in which case it will add
%               the suffix '_Recon' to the ScanDir & the item in the CurationTable
%               Finally, it acquires the SeriesDescription (GetSeriesDescription) & 
%               StudyInstanceUID & SeriesInstanceUID. For site 040, the latter two are
%               replaced by StudyDate & SeriesTime, as this site had some scans
%               anonimized locally, which changed the UIDs.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [CurationTable, ScanList, StudyInstanceUID, SeriesInstanceUID, dcmList, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime] = GetScanInfo(CurationTable, ScanList, ScanDir, StudyInstanceUID, SeriesInstanceUID, SiteNrString, iSite, iSubject, iVisit, iSubjectVisit, iScanType, iScan, bUseDCMTK, EchoTime, RepetitionTime, SeriesDescription, SeriesDate, SeriesTime)
% __________________________________
% Copyright 2015-2019 ExploreASL



if ~exist(ScanDir,'dir')
    nSlices  = 0;
    dcmList = [];
else
    dcmList = xASL_adm_GetFileList(ScanDir,'^.*[^.(csv|tsv|gz)]$','FPListRec',[0 Inf]); % \.dcm|.*\.img|.*\.IMA|
    nSlices = length(dcmList);
end

CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+iScan} = num2str(nSlices); % number of files for this scan

% Obtain DICOM identifying data from the average DICOM image
if ~isempty(dcmList)
    IndexDCM = ceil(0.5*length(dcmList));

    Info = xASL_io_ReadTheDicom(bUseDCMTK, dcmList{IndexDCM});

    if isfield(Info,'EchoTime')
        EchoTime(iSite,iSubject,iVisit,iScanType,iScan) = xASL_str2num(Info.EchoTime);
    end
    if isfield(Info,'RepetitionTime')
        RepetitionTime(iSite,iSubject,iVisit,iScanType,iScan) = xASL_str2num(Info.RepetitionTime);
    end
    
    if isfield(Info,'SeriesDate')
        SeriesDate{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.SeriesDate;
    elseif isfield(Info,'AcquisitionDate')
        SeriesDate{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.AcquisitionDate;
    else
        SeriesDate{iSite}{iSubject}{iVisit}{iScanType}{iScan} = [];
    end

    if isfield(Info,'SeriesTime')
        SeriesTime{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.SeriesTime;
    elseif isfield(Info,'AcquisitionTime')
        SeriesTime{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.AcquisitionTime;
    else
        SeriesTime{iSite}{iSubject}{iVisit}{iScanType}{iScan} = [];
    end    
    
    % ImageType is a field that can designate reconstruction scans
    if isfield(Info,'ImageType')
        if ~isempty(Info.ImageType)
            ImageType   = Info.ImageType;
            if ~isempty(regexp(lower(ImageType), '(derived|projection)', 'once'))
                % add to list
                CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+iScan} = [CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+iScan} '_Recon'];
                
                % add to foldername
                ScanList{iScan} = [RemoveExistingRegExp(ScanList{iScan}, '_recon') '_Recon'];
                NewPath = fullfile(fileparts(ScanDir), ScanList{iScan});
                
                if ~strcmp(ScanDir, NewPath)
                    while exist(NewPath,'dir')
                        NewPath = [NewPath '_2'];
                    end                    
                    xASL_Move(ScanDir, NewPath);
                end
            end
        end
    end        

    % Get the Series Description as well
    SeriesDescription{iSite}{iSubject}{iVisit}{iScanType}{iScan} = GetSeriesDescription(Info);

    if  strcmp(SiteNrString{iSite},'040') % VUmc
        % For VUmc the UIDs are always different, because scans were anonymized locally
        % So we cannot use them to check for different UIDs for the same
        % scan session/MRI examination, as we do in function VerifyIdenticalStudyUID
        StudyInstanceUID{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.StudyDate;
        SeriesInstanceUID{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.SeriesTime;
    else
        StudyInstanceUID{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.StudyInstanceUID;
        SeriesInstanceUID{iSite}{iSubject}{iVisit}{iScanType}{iScan} = Info.SeriesInstanceUID;
    end
end


end






%% =============================================================================================================================
%% =============================================================================================================================
function [dcmList] = CorrectDICOMFolderDepth(dcmList, ScanDir)
%CorrectDICOMFolderDepth Move DICOMs to correct folder layer (i.e. //RawDir/SiteDir/SubjectDir/VisitDir/ScanTypeDir/ScanDir/DICOMs)
% 
%
% FORMAT: [dcmList] = CorrectDICOMFolderDepth(dcmList, ScanDir)
% 
% INPUT:
%   dcmList     - list of paths to DICOM files we want to put in correct folder layer
%   ScanDir     - path to folder for one scan/sequence -> //RawDir/SiteDir/SubjectDir/VisitDir/ScanTypeDir/ScanDir
%                             
% OUTPUT:
%   dcmList     - as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions move DICOMs to correct folder layer, if they are "too deep".
%              Sometimes, for whatever reason, the DICOMs are not directly stored in the ScanDir,
%              but there are some subfolders in between. These are confusing for the import script.
%              Here we check whether the "in-between layers" are indeed empty, and will then solve this
%              by moving the DICOMs. If the "in-between layers" are not empty, we might have multiple DICOM series
%              stored into a single ScanDir, which will issue a warning.
%              
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [dcmList] = CorrectDICOMFolderDepth(dcmList, ScanDir)
% __________________________________
% Copyright 2015-2019 ExploreASL

if ~isempty(dcmList) && ~isempty(ScanDir)
    DICOMpath = xASL_fileparts(dcmList{1});
    if ~strcmp(ScanDir, DICOMpath)
        
        AllList = xASL_adm_GetFileList(ScanDir,'.*','FPListRec',[0 Inf], false);
        
        fprintf('Moving DICOMs to correct directory level:\n');
        fprintf('%s\n', ['From ' DICOMpath]);
        fprintf('%s', ['to ' ScanDir '   ']);

        % Move files & delete dirs
        if exist(ScanDir,'dir')
            for iD=1:length(AllList)
                OldPath = AllList{iD};
                [~, Ffile, Fext] = xASL_fileparts(AllList{iD});
                NewPath = fullfile(ScanDir, [Ffile Fext]);
                if ~strcmp(OldPath,NewPath) % They shouldnot be the same
                    while exist(NewPath, 'file') % make sure not to overwrite files, but add a counter suffix instead
                        NewPath = fullfile(ScanDir, [Ffile '_2' Fext]);
                    end
                    xASL_Move(AllList{iD}, NewPath);
                end
            end
        end

        Dirs2Del = xASL_adm_GetFileList(ScanDir,'.*','FPListRec',[0 Inf], true);
        for iD=[length(Dirs2Del):-1:1]
            if isempty(xASL_adm_GetFileList(Dirs2Del{iD},'.*','FPListRec',[0 Inf]))
                rmdir(Dirs2Del{iD});
            end
        end

        fprintf('\n'); % now recreate the list to be sure
        dcmList = xASL_adm_GetFileList(ScanDir,'^.*[^.(csv|tsv|gz)]$','FPListRec',[0 Inf]); % \.dcm|.*\.img|.*\.IMA|
    end
end

end





%% =============================================================================================================================                
%% =============================================================================================================================
function [UIDWarningReported, DifferentStudyUIDList] = VerifyIdenticalStudyUID(UIDWarningReported, StudyInstanceUID, iSite, iSubject, SubjectList, iVisit, VisitList, DifferentStudyUIDList)
%VerifyIdenticalStudyUID Check that all scans belong to the same MRI examination visit
%
% FORMAT: [UIDWarningReported, DifferentStudyUIDList] = VerifyIdenticalStudyUID(UIDWarningReported, StudyInstanceUID, iSite, iSubject ,SubjectList, iVisit, VisitList, DifferentStudyUIDList)
% 
% INPUT:
%   UIDWarningReported      - true if a warning was thrown already for this subject/visit
%                             to avoid throwing multiple warnings for the same MRI examination/visit
%   StudyInstanceUID        - Unique Identifier of current MRI examination/visit/study/scan session
%                             This is built up out of organization-unique values, and the latest value is an ID specified by
%                             the organization (image center, scanner), which could be e.g. DateTime
%   iSite                   - current site index
%   iSubject                - current subject index
%   SubjectList             - list of subjects
%   iVisit                  - index for current visit
%   VisitList               - list of visits
%   DifferentStudyUIDList   - list of subjects/visits with scans from different MRI examinations/visits/studies/scan sessions
%                             
% OUTPUT:
%   UIDWarningReported      - as input
%   DifferentStudyUIDList   - as input
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function compares StudyInstanceUIDs between scans from the same subject/visit/scan session
%              and throws a warning if they differ, keeping track that we only throw a warning per subject/visit once,
%              and that we store this in the DifferentStudyUIDList.
%              
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [UIDWarningReported, DifferentStudyUIDList] = VerifyIdenticalStudyUID(UIDWarningReported, StudyInstanceUID, iSite, iSubject, SubjectList, iVisit, VisitList, DifferentStudyUIDList)
% __________________________________
% Copyright 2015-2019 ExploreASL


if length(StudyInstanceUID)>=iSite && length(StudyInstanceUID{iSite})>=iSubject && length(StudyInstanceUID{iSite}{iSubject})>=iVisit
    for iScan1=1:length(StudyInstanceUID{iSite}{iSubject}{iVisit})
        if ~isempty(StudyInstanceUID{iSite}{iSubject}{iVisit}{iScan1})
            for iScan2=1:length(StudyInstanceUID{iSite}{iSubject}{iVisit})
                if ~isempty(StudyInstanceUID{iSite}{iSubject}{iVisit}{iScan2})
                    for iSequence=1:length(StudyInstanceUID{iSite}{iSubject}{iVisit}{iScan2})
                        if ~isempty(StudyInstanceUID{iSite}{iSubject}{iVisit}{iScan2}{iSequence}) && ~strcmp(StudyInstanceUID{iSite}{iSubject}{iVisit}{iScan1}{1},StudyInstanceUID{iSite}{iSubject}{iVisit}{iScan2}{iSequence})
                            if ~UIDWarningReported % issue warning only once per subject
                                DifferentStudyUIDList.(['n' num2str(length(fields(DifferentStudyUIDList)))]) = [SubjectList{iSubject} '_' VisitList{iVisit}];
%                               warning([SubjList{iS} ' had scans with different StudyInstanceUIDs']);
                                UIDWarningReported = true;
                            end
                        end
                    end
                end
            end
        end
    end
end


end




%% =============================================================================================================================
%% =============================================================================================================================
function [CurationTable, ScanList] = AssignRun(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubject, iVisit, iSubjectVisit, SeriesTime)
%AssignRun Assign the 'run' of each scan for scan repetitions within MRI examination (i.e. run-1 run-2 etc)
%
% FORMAT: [CurationTable, ScanList] = AssignRun(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubject, iVisit, iSubjectVisit, SeriesTime)
% 
% INPUT:
%   CurationTable   - This Table contains the folder structure as found and curated. It is also stored as DCM_StructureList.csv & DCM_StructureList.mat files
%   ScanList        - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%   ScanTypeDir ? ? - the folder where the ScanType is located (i.e. //raw/Site/Subject/Visit/ScanType)
%   iSite ? ? ? ? ? - index of current site
%   iScanType ? ? ? - index of current ScanType
%   iSubject? ? ? ? - index of current Subject
%   iVisit          - index for current visit
%   iSubjectVisit   - combined index/counter for subject/visit rows
%   SeriesTime      - cell structure containing scanning time, structured as {iSite}{iSubject}{iScanType}{iScan}
%
%                             
% OUTPUT:
%   CurationTable   - as input
%   ScanList        - as input
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function assigns 'run' numbers to complete scans.
%              This run number is added as suffix (e.g. '_run-1', '_run-2', etc) to the ScanDir & the item in the CurationTable.
%              Complete scans are defined as having the desired number of DICOMs (this was defined as complete previously)
%              & not being a 'reconstruction' or 'inactive'. If a complete active scan doesn't exist, we take the complete inactive
%              scan as '_run-1'. Otherwise, the 'runs' are sorted by SeriesTime (if this exists)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [CurationTable, ScanList] = AssignRun(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubject, iVisit, iSubjectVisit, SeriesTime)
% __________________________________
% Copyright 2015-2019 ExploreASL


CompleteList = cellfun(@(x) ~isempty(strfind(lower(x), 'complete')), CurationTable{iSite}{iScanType}(iSubjectVisit+1,6+1:end));
ReconList = cellfun(@(x) ~isempty(strfind(lower(x), 'recon')), ScanList);
InactiveList = cellfun(@(x) ~isempty(strfind(lower(x), 'inactive')), ScanList);

% Equalize length lists, zerofilling
MaxL = max([length(CompleteList), length(ReconList), length(InactiveList)]);
CompleteList(length(CompleteList)+1:MaxL) = 0;
ReconList(length(ReconList)+1:MaxL) = 0;
InactiveList(length(InactiveList)+1:MaxL) = 0;

% Find complete scans
CompleteList = CompleteList & ~ReconList;
IsRun = find(CompleteList & ~InactiveList); % try first to find active complete scans
if isempty(IsRun)
    IsRun = find(CompleteList & InactiveList); % otherwise get inactive complete scans
end

if ~isempty(IsRun) % check if there are complete scans
    % First try to see if we have timing data
    SortTime = 0;
    for iR=1:length(IsRun)
        SortTime(iR,1) = iR;
        try
            SortTime(iR,2) = SeriesTime{iSite}{iScanType}{iSubject}{iVisit}{IsRun(iR)};
            TimeExist = true;
        catch
            % If time was incomplete
            TimeExist = false;
        end
    end
    
    % Check if the timing data is valid
    if TimeExist
        if sum(isfinite(SortTime(:,2)) & SortTime(:,2)~=0)==length(SortTime)
            SortTime = sortrows(SortTime,2); % if yes, then sort IsRun according to it
            IsRun = SortTime(:,1)';
        end % otherwise, keep IsRun as it was (alphabetically sorted on foldername)
    end
    
    for iR=1:length(IsRun) % Add the run-n suffixes to the foldernames
        RunString = ['run-' num2str(iR)];

        TempList = CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+IsRun(iR)};
        TempList = RemoveExistingRegExp(TempList, '_run-\d*');
        CurationTable{iSite}{iScanType}{iSubjectVisit+1,6+IsRun(iR)} = [TempList '_' RunString];
        
        ScanPath = fullfile(ScanTypeDir, ScanList{IsRun(iR)});
        ScanList{IsRun(iR)} = RemoveExistingRegExp(ScanList{IsRun(iR)}, '_run-\d*');
        ScanList{IsRun(iR)} = RemoveExistingRegExp(ScanList{IsRun(iR)}, '_complete');
        ScanList{IsRun(iR)} = [ScanList{IsRun(iR)} '_' RunString];
        NewPath = fullfile(ScanTypeDir, ScanList{IsRun(iR)});

        if ~strcmp(ScanPath, NewPath)
            while exist(NewPath,'dir')
                NewPath = [NewPath '_2'];
            end            
            xASL_Move(ScanPath, NewPath);
        end
    end
end


end


%% =============================================================================================================================
%% =============================================================================================================================
function [ScanList, ActiveMissingList] = AssignActiveMissing(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubjectVisit, ActiveMissingList)
%AssignActiveMissing Assigns an inactive complete scan as active if active complete scan is missing
%
% FORMAT: [ScanList, ActiveMissingList] = AssignActiveMissing(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubjectVisit, ActiveMissingList)
% 
% INPUT:
%   CurationTable     - This Table contains the folder structure as found and curated. It is also stored as DCM_StructureList.csv & DCM_StructureList.mat files
%   ScanList          - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%   ScanTypeDir ? ?   - the folder where the ScanType is located (i.e. //raw/Site/Subject/Visit/ScanType)
%   iSite ? ? ? ? ?   - index of current site
%   iScanType ? ? ?   - index of current ScanType
%   iSubjectVisit     - combined index/counter for subject/visit rows
%   ActiveMissingList - cell structure containing names of scans where an active complete was missing, but inactive complete existed & was used instead
%                             
% OUTPUT:
%   ScanList           - as input
%   ActiveMissingList  - as input
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function assigns an inactive complete scan as active if an active complete scan from the same scan session for the same ScanType is missing
%              This removes the "inactive_" prefix of the ScanDir, allowing dcm2nii to import this scan (i.e. not skipping it)
%              These instances are tracked in the ActiveMissingList, which
%              can be checked in a JSON file
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ScanList, ActiveMissingList] = AssignActiveMissing(CurationTable, ScanList, ScanTypeDir, iSite, iScanType, iSubjectVisit, ActiveMissingList)
% __________________________________
% Copyright 2015-2019 ExploreASL

CompleteList = cellfun(@(x) ~isempty(strfind(lower(x), '_run')), CurationTable{iSite}{iScanType}(iSubjectVisit+1,6+1:end));
InactiveList = cellfun(@(x) ~isempty(strfind(lower(x), 'inactive')), ScanList);

% Equalize length lists, zerofilling
MaxL = max([length(CompleteList), length(InactiveList)]);
CompleteList(length(CompleteList)+1:MaxL) = 0;
InactiveList(length(InactiveList)+1:MaxL) = 0;

ActiveMissing = CompleteList & ~InactiveList;
InactivePresent = find(CompleteList & InactiveList);

if ~sum(ActiveMissing) && ~isempty(InactivePresent) % check if there are only inactive complete scans
    for iP=1:length(InactivePresent)
        
        ScanPath = fullfile(ScanTypeDir, ScanList{InactivePresent(iP)});
        ScanList{InactivePresent(iP)} = RemoveExistingRegExp(ScanList{InactivePresent(iP)}, 'Inactive_');
        ActiveMissingList.(['n' num2str(length(fields(ActiveMissingList)))]) = ScanList{InactivePresent(iP)}; % add new name to list
        NewPath = fullfile(ScanTypeDir, ScanList{InactivePresent(iP)}); % assign new foldername without 'inactive'

        if ~strcmp(ScanPath, NewPath) % change the foldername
            while exist(NewPath,'dir')
                NewPath = [NewPath '_2'];
            end            
            xASL_Move(ScanPath, NewPath);
        end
    end
end

end


%% =============================================================================================================================
%% =============================================================================================================================
function [CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, RemovedDuplicateList, DeleteActiveList] = RemoveDuplicates(CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, iSite, iSubject, iVisit, iSubjectVisit, iScanType, RemovedDuplicateList, ScanTypeDir, DeleteActiveList)
%RemoveDuplicates Find true duplicate scans & delete them
%
% FORMAT: [CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, RemovedDuplicateList, DeleteActiveList] = RemoveDuplicates(CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, iSite, iSubject, iVisit, iSubjectVisit, iScanType, RemovedDuplicateList, ScanTypeDir, DeleteActiveList)
% 
% INPUT:
%   CurationTable        - This Table contains the folder structure as found and curated. It is also stored as DCM_StructureList.csv & DCM_StructureList.mat files
%   ScanList             - list of scannames/ScanNames (from directory names), as in //raw/Site/Subject/Visit/ScanType/ScanName
%   SeriesInstanceUID    - Unique Identifier of current scan/sequence, consisting of organization UID etc
%   StudyInstanceUID     - Unique Identifier of current MRI examination/visit/study/scan session, consisting of organization UID etc
%   SeriesDescription    - DICOM field containing the name of the sequence (e.g. Axial_T2W_TSE_with_Fat_Sat_CLEAR).
%   EchoTime             - matrix containing EchoTimes, structured as (iSite,iSubject,iScanType,iScan)
%   RepetitionTime       - matrix containing RepetitionTimes, structured as (iSite,iSubject,iScanType,iScan)
%   iSite ? ? ? ? ?      - index of current site
%   iSubject? ? ? ?      - index of current Subject
%   iVisit                - index for current visit
%   iSubjectVisit         - combined index/counter for subject/visit rows
%   iScanType ? ? ?      - index of current ScanType
%   RemovedDuplicateList - cell structure containing names of deleted scans
%   ScanTypeDir ? ?   - the folder where the ScanType is located (i.e. //raw/Site/Subject/Visit/ScanType)
%   DeleteActiveList  - cell structure containing names of scans where an active complete was missing, but inactive complete existed & was used instead
%                             
% OUTPUT: all same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function searches for true duplicate scans & deletes them.
%              If a duplicate is found, it is removed and added to the cell structure RemovedDuplicateList.
%              Missing active scans are added to the cell structure ActiveMissingList.
%              The input lists/tables are adapted accordingly (i.e. deleted scans removed)
%              
%              Duplicates are defined as having equal:
%              SeriesInstanceUID
%              StudyInstanceUID
%              SeriesDescription (or ProtocolName)
%              EchoTime
%              RepetitionTime
%              number of slices/DICOMs (in the CurationTable)
% 
%              If a duplicate is found, this doesnt mean that the scans are
%              equally stored. The can have the active designation from
%              IXICO, and they can be incompletely downloaded.
%              Here the preference is: 
%              a) keep scan with most slices (assuming that the more slices, the more complete a scan is)
%              b) to delete inactive & keep active scan
%              c) to delete the second & keep the first scan
%              If this leads to deleting active scans & keeping inactive
%              ones, this is tracked in DeleteActiveList
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, RemovedDuplicateList, DeleteActiveList] = RemoveDuplicates(CurationTable, ScanList, SeriesInstanceUID, StudyInstanceUID, SeriesDescription, EchoTime, RepetitionTime, iSite, iSubject, iVisit, iSubjectVisit, iScanType, RemovedDuplicateList, ScanTypeDir, DeleteActiveList)
% __________________________________
% Copyright 2015-2019 ExploreASL

                        
HasUIDs = false;
if  length(SeriesInstanceUID)>=iSite && length(StudyInstanceUID)>=iSite && length(SeriesDescription)>=iSite && size(EchoTime,1)>=iSite && size(RepetitionTime,1)>=iSite
    if length(SeriesInstanceUID{iSite})>=iSubject && length(StudyInstanceUID{iSite})>=iSubject && length(SeriesDescription{iSite})>=iSubject && size(EchoTime,2)>=iSubject && size(RepetitionTime,2)>=iSubject
        if length(SeriesInstanceUID{iSite}{iSubject})>=iVisit && length(StudyInstanceUID{iSite}{iSubject})>=iVisit && length(SeriesDescription{iSite}{iSubject})>=iVisit && size(EchoTime,3)>=iVisit && size(RepetitionTime,3)>=iVisit        
            if length(SeriesInstanceUID{iSite}{iSubject}{iVisit})>=iScanType && length(StudyInstanceUID{iSite}{iSubject}{iVisit})>=iScanType && length(SeriesDescription{iSite}{iSubject}{iVisit})>=iScanType && size(EchoTime,4)>=iScanType && size(RepetitionTime,4)>=iScanType
                HasUIDs = true;
            end
        end
    end
end

if ~HasUIDs
    return;
end

CompareSeriesUID = SeriesInstanceUID{iSite}{iSubject}{iVisit}{iScanType};
CompareStudyUID = StudyInstanceUID{iSite}{iSubject}{iVisit}{iScanType};

CompareSlicesN = cellfun(@(x) xASL_adm_CatchNumbersFromString(x), FixList(CurationTable{iSite}{iScanType}(iSubjectVisit+1,6+1:end)));
CompareEchoTime = FixList(EchoTime(iSite,iSubject,iVisit,iScanType,:));
CompareRepetitionTime = FixList(RepetitionTime(iSite,iSubject,iVisit,iScanType,:));

LengthN(1) = length(CompareSeriesUID);
LengthN(2) = length(CompareStudyUID);
LengthN(3) = length(CompareEchoTime);
LengthN(4) = length(CompareRepetitionTime);
LengthN(5) = length(CompareSlicesN);
MaxL = min(LengthN); % bugfix

ScanNDeleted = false(MaxL,1);

PerformedList = [];

for iL1=1:MaxL
    for iL2=1:MaxL
        IsProcessed = false;
        for iP=1:size(PerformedList,1)
            if (min(PerformedList(iP,:)==[iL1 iL2]))
                IsProcessed = true;
            end
        end
        
        if iL1~=iL2 && ~IsProcessed
            SameSeriesUID = strcmp(CompareSeriesUID{iL1}, CompareSeriesUID{iL2});
            SameStudyUID = strcmp(CompareStudyUID{iL1}, CompareStudyUID{iL2});
            SameEchoTime = CompareEchoTime(iL1) == CompareEchoTime(iL2);
            SameRepetitionTime = CompareRepetitionTime(iL1) == CompareRepetitionTime(iL2);

            if SameSeriesUID && SameStudyUID && SameEchoTime && SameRepetitionTime

                % Now we are reasonably convinced that this is a true duplicate
                % Now check which one to delete, according to preference described in this header

                if CompareSlicesN(iL1)>1.05*CompareSlicesN(iL2)
                    IndexN = iL2; % second has significantly (5%) fewer slices
                    KeepN = iL1;
                elseif CompareSlicesN(iL2)>1.05*CompareSlicesN(iL1)
                    IndexN = iL1; % first has significantly fewer (5%) slices
                    KeepN = iL2;
                elseif isempty(findstr(lower(ScanList{iL1}), 'inactive')) && ~isempty(strfind(lower(ScanList{iL2}), 'inactive'))
                    % the second is inactive, delete the second
                    IndexN = iL2;
                    KeepN = iL1;
                elseif isempty(findstr(lower(ScanList{iL2}), 'inactive')) && ~isempty(strfind(lower(ScanList{iL1}), 'inactive'))
                    % the first is inactive, delete the first
                    IndexN = iL1;
                    KeepN = iL2;
                else % either both are inactive, or both are active, then delete the second
                    IndexN = iL2;
                    KeepN = iL1;
                end
                
                Dir2Delete = fullfile(ScanTypeDir, ScanList{IndexN});
                
                if exist(Dir2Delete, 'dir')
                    rmdir(Dir2Delete, 's');
                    ScanNDeleted(IndexN) = true;
                    RemovedDuplicateList{end+1,1} = Dir2Delete;
                end
                if isempty(findstr(lower(ScanList{IndexN}), 'inactive')) && ~isempty(findstr(lower(ScanList{KeepN}), 'inactive'))
                    % if we keep the inactive one and delete the active one, which must be because the inactive is more complete
                    DeleteActiveList.(['n' num2str(length(fields(DeleteActiveList)))]) = Dir2Delete;
                end
            end
        end
        PerformedList(end+1,:) = [iL1 iL2];
        PerformedList(end+1,:) = [iL2 iL1];
    end
end

if sum(ScanNDeleted)>0
    
    % Now remove the deleted scans from the lists
    % 1) grab what we want to insert
    TempCell = SeriesInstanceUID{iSite}{iSubject}{iVisit}{iScanType}(~ScanNDeleted);
    % 2) empty the cell
    SeriesInstanceUID{iSite}{iSubject}{iVisit}{iScanType}(:) = {[]};
    % 3) insert the contents
    SeriesInstanceUID{iSite}{iSubject}{iVisit}{iScanType}(1:length(TempCell)) = TempCell;

    % Repeat this for other cells
    TempCell = StudyInstanceUID{iSite}{iSubject}{iVisit}{iScanType}(~ScanNDeleted);
    StudyInstanceUID{iSite}{iSubject}{iVisit}{iScanType}(:) = {[]};
    StudyInstanceUID{iSite}{iSubject}{iVisit}{iScanType}(1:length(TempCell)) = TempCell;

    TempCell = SeriesDescription{iSite}{iSubject}{iVisit}{iScanType}(~ScanNDeleted);
    SeriesDescription{iSite}{iSubject}{iVisit}{iScanType}(:) = {[]};
    SeriesDescription{iSite}{iSubject}{iVisit}{iScanType}(1:length(TempCell)) = TempCell;
    
    % Now for numbers
    TempCell = squeeze(EchoTime(iSite,iSubject,iVisit,iScanType,~ScanNDeleted))';
    EchoTime(iSite,iSubject,iVisit,iScanType,:) = false(1,size(EchoTime,5));
    EchoTime(iSite,iSubject,iVisit,iScanType,1:length(TempCell)) = TempCell;

    TempCell = squeeze(RepetitionTime(iSite,iSubject,iVisit,iScanType,~ScanNDeleted))';
    RepetitionTime(iSite,iSubject,iVisit,iScanType,:) = false(1,size(RepetitionTime,5));
    RepetitionTime(iSite,iSubject,iVisit,iScanType,1:length(TempCell)) = TempCell;

    % Now for part of cellist
    TempCell = CurationTable{iSite}{iScanType}(iSubjectVisit+1,6+find(~ScanNDeleted)');
    CurationTable{iSite}{iScanType}(iSubjectVisit+1,6+1:end) = {[]};

    CurationTable{iSite}{iScanType}(iSubjectVisit+1,6+1:6+length(TempCell)) = TempCell;

    % now for ScanList
    ScanList = ScanList(~ScanNDeleted);
end


end



%% =============================================================================================================================
%% =============================================================================================================================
function [SeriesDescription, IsExamcard] = GetSeriesDescription(Info)
%GetSeriesDescription Extracts SeriesDescription from DICOM header
%
% FORMAT: [SeriesDescription, IsExamcard] = GetSeriesDescription(Info)
% 
% INPUT:
%   Info              - DICOM header as struct, containing relevant fields
%
% OUTPUT: 
%   SeriesDescription - String containing the SeriesDescription (empty if not found)
%   IsExamcard        - true if DICOM is an examcard
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function tries to extract the SeriesDescription from a DICOM header, as follows:
%              A) if Info.SeriesDescription but no Info.ProtocolName -> use Info.SeriesDescription
%              B) if Info.ProtocolName but no Info.SeriesDescription -> use Info.ProtocolName
%              C) if exist(both), concatenate both if they are unequal
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [SeriesDescription, IsExamcard] = GetSeriesDescription(Info)
% __________________________________
% Copyright 2015-2019 ExploreASL


SeriesDescription = '';
if isfield(Info,'SeriesDescription') && ~isempty(Info.SeriesDescription)
   SeriesDescription = Info.SeriesDescription;
end

if isfield(Info,'ProtocolName') && ~isempty(Info.ProtocolName)
    if isempty(SeriesDescription)
        SeriesDescription = Info.ProtocolName;
    elseif ~strcmp(SeriesDescription,Info.ProtocolName)
        SeriesDescription = [SeriesDescription '_' Info.ProtocolName];
    end
end

if isempty(SeriesDescription) && isfield(Info,'ImageType') && ~isempty(Info.ImageType)
    SeriesDescription = xASL_adm_CorrectName(Info.ImageType, [], '*'); % only add ImageType if nothing else works
end    

if strcmp(lower(SeriesDescription), 'examcard')
    IsExamcard = true;
    SeriesDescription = '';
else
    IsExamcard = false;
end
    
end



%% =============================================================================================================================
%% =============================================================================================================================
function [ThroughCell] = FixList(ThroughCell)
%FixList Cleans a cell structure (e.g. list/table) by removing empty cells
%
% FORMAT: [ThroughCell] = FixList(ThroughCell)
% 
% INPUT/OUTPUT:
%   ThroughCell - cell structure we want to clean
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: There can be inequal lengths of list rows, because of differences between scans between ScanTypes,
%              Subjects, Visits, etc. Checking this in scripts provides errors, this script removes the empty cells
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ThroughCell] = FixList(ThroughCell)
% __________________________________
% Copyright 2015-2019 ExploreASL

if isnumeric(ThroughCell)
    for iFix=length(ThroughCell):-1:1
        if ThroughCell(iFix)==0
            ThroughCell = ThroughCell(1:iFix-1);
        end
    end
else
    for iFix=length(ThroughCell):-1:1
        if isempty(ThroughCell{iFix})
            ThroughCell = ThroughCell(1:iFix-1);
        end
    end
end

end



%% =============================================================================================================================
%% =============================================================================================================================
function [ThroughputString] = RemoveExistingRegExp(ThroughputString, InputRegExp)
%RemoveExistingRegExp Remove requested regular expression match from string
%
% FORMAT: [ThroughputString] = RemoveExistingRegExp(ThroughputString, InputRegExp)
% 
% INPUT:
%   ThroughputString - string we want to adapt
%   InputRegExp      - regular expression to match string we want to remove
%                      from the ThroughputString
%
% OUTPUT:
%   ThroughputString - string we want to adapt
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function removes part of string from a string, using a regular expression.
%              This is useful when we want to add a suffix or prefix to a folder or CurationTable name, but make
%              sure that we are not getting multiple suffixes or prefixes
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ThroughputString] = RemoveExistingRegExp(ThroughputString, InputRegExp)
% __________________________________
% Copyright 2015-2019 ExploreASL

[Ind1, Ind2] = regexp(lower(ThroughputString),lower(InputRegExp));

while ~isempty(Ind1) && ~isempty(Ind2)
    ThroughputString = [ThroughputString(1:Ind1(1)-1) ThroughputString(Ind2(1)+1:end)];
    [Ind1, Ind2] = regexp(lower(ThroughputString),lower(InputRegExp));
end

end
