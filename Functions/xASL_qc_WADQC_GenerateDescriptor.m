function xASL_qc_WADQC_GenerateDescriptor(x, iSubject, ScanTypeIs)
%xASL_qc_WADQC_GenerateDescriptor QC function to generate WAD-QC descriptor JSON
%
% FORMAT: xASL_qc_WADQC_GenerateDescriptor(x, iSubject)
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject    - index of current subject (REQUIRED)
%
% OUTPUT: n/a
% OUTPUT FILE:
%    WAD_Path = fullfile(x.D.ROOT, SubPath{iScan}, ['qcdc_' x.SUBJECTS{iSubject} '_' ScanTypes{iScan} '.json'])
%               This file is created for each subject/scantype, and contains the
%               descriptor for QCDC, such that it knows where it can find QC
%               information & images to later embed into the dummy DICOM, and later
%               send to the WAD-QC server
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This QC function generates a JSON descriptor for Gaspare'
%              QCDC script, by the following steps:
%
%              - a) include information about where to find the dummy DICOM (i.e. placeholder DICOM)
%              - b) For ExploreASL' QC fields (as passed through in
%                   x.Output), here we note all these QC fields for each
%                   ScanType, as the x.Output should have been collected
%                   equally in the QC file 'QC_collection_SubjectName.json'
%                   by function xASL_qc_CollectParameters
%              - c) Subfunction xASL_qc_WADQC_images - Includes visual standard space QC
%                   images, by searching them on prescribed paths within the
%                   Population folder (where currently all derivatives reside)
%              - d) Insert the PDF report; this PDF report is
%                   subject-specific, not scan-specific. For completeness it
%                   is added to each QCDC descriptor
%              - e) Add WAD-QC server details (i.e. IP address etc)
%              - f) Save the Descriptor JSON file.
% 
% EXAMPLE: xASL_qc_WADQC_GenerateDescriptor(x, iSubject)
% For more information about WAD-QC please visit: https://github.com/wadqc/WAD_Documentatie/wiki
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

 
    
if ~x.DoWADQCDC
    return;
else
    fprintf('Creating WAD-QC descriptor\n');
end

%% Define Scantypes
ScanTypes = {'Structural', 'ASL', 'func', 'dwi'};
SubPath = {x.SUBJECTS{iSubject}, fullfile(x.SUBJECTS{iSubject},...
    'ASL_1'), fullfile(x.SUBJECTS{iSubject}, 'func'), fullfile(x.SUBJECTS{iSubject}, 'dwi')};
NIfTIname = {'T1', 'ASL4D', 'func', 'dwi'};

%% Define scanType to process
if nargin>2 && ~isempty(ScanTypeIs)
    Indices = find(cellfun(@(x) strcmp(x,ScanTypeIs), ScanTypes));
else
    Indices = 1:length(ScanTypes)
end

% Define folders to check in population folder for QC images
FolderList{1} = {'FLAIRCheck' 'T1Check' 'TissueVolume'};
FolderList{2} = {'ASLCheck' 'M0Check' 'M0Reg_ASL' 'MotionASL' 'RawSourceIMCheck' 'SD_SNR' 'SliceGradientCheck'};
FolderList{3} = {'FuncCheck'};
FolderList{4} = {'dwiCheck'};

for iScan=Indices
    
    if ~exist(fullfile(x.D.ROOT, SubPath{iScan}), 'dir')
        continue; % skip this ScanType
    elseif isempty(xASL_adm_GetFileList(fullfile(x.D.ROOT, SubPath{iScan}), [NIfTIname{iScan} '.*\.nii'], [], [0 Inf]))
        continue; % skip this ScanType
    end
    
    wadqc = struct;
    %% a) DICOM placeholder
    wadqc.dicom_meta.dicom_wrapper.sub_path = SubPath{iScan};
    DummyFile = ['DummyDicom_' NIfTIname{iScan} '*.dcm'];
    wadqc.dicom_meta.dicom_wrapper.filename = DummyFile;
    wadqc.dicom_meta.dicom_wrapper.filepath_wadqc_placeholder = 'path_placeholder/';

    %% -----------------------------------------------------------------------
    %% b) QC parameters
    if isfield(x.Output,ScanTypes{iScan})
        Fields2 = fieldnames(x.Output.(ScanTypes{iScan}));
        for iL=1:length(Fields2)
            CurrFieldName = [ScanTypes{iScan} '_' Fields2{iL}];
            wadqc.qc_items.(CurrFieldName).type = 'json';
            wadqc.qc_items.(CurrFieldName).sub_path = x.SUBJECTS{iSubject};
            wadqc.qc_items.(CurrFieldName).filename = ['QC_collection_' x.SUBJECTS{iSubject} '.json'];
            wadqc.qc_items.(CurrFieldName).child = [ScanTypes{iScan} '/' Fields2{iL}];

            FieldContents = x.Output.(ScanTypes{iScan}).(Fields2{iL});
            if  isnumeric(FieldContents)
                wadqc.qc_items.(CurrFieldName).category = 'float';
            elseif ischar(FieldContents)
                wadqc.qc_items.(CurrFieldName).category = 'string';
            end
        end
    end



    %% -----------------------------------------------------------------------
    %% c) Include visual QC .jpg images
    wadqc = xASL_qc_WADQC_images(x, wadqc, iSubject, iScan, FolderList);


    %% -----------------------------------------------------------------------
    %% d) PDF report
    PDFFileName = ['xASL_Report_' x.SUBJECTS{iSubject} '.pdf'];

    itemName = ['PDF_' x.SUBJECTS{iSubject}];

    % Make sure that all illegal characters are replaced
    itemName = xASL_adm_CorrectName(itemName, 1);

    wadqc.qc_items.(itemName).type = 'file';  
    wadqc.qc_items.(itemName).sub_path = x.SUBJECTS{iSubject}; % this would change with multiple subjects/scans
    wadqc.qc_items.(itemName).filename = PDFFileName;
    wadqc.qc_items.(itemName).category = 'object';




    %% -----------------------------------------------------------------------
    %% e) WAD-QC Server
    wadqc.wad_qc_server.ip_address = 'wad-qc'; % 145.121.126.53/ localhost
    wadqc.wad_qc_server.port = '11112';
    wadqc.wad_qc_server.ae_title = 'dummy';

    %% f) Save
    WAD_Path = fullfile(x.D.ROOT, SubPath{iScan}, ['qcdc_' x.SUBJECTS{iSubject} '_' ScanTypes{iScan} '.json']);
    fclose all;
    xASL_delete(WAD_Path);
    spm_jsonwrite(WAD_Path, wadqc);
        
end

end


%% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
function [wadqc] = xASL_qc_WADQC_images(x, wadqc, iSubject, iScan, FolderList)
%xASL_qc_WADQC_images Include visual QC .jpg images
% % This script will look for all images, just to make it for 1 subject only, with the current run,
% let it search for a single subject ID (which exists in the filename of both the structural & ASL images)


    for iFolder=1:length(FolderList{iScan})
        SubFolder = fullfile(x.D.PopDir, FolderList{iScan}{iFolder});
        
        JPGlist = xASL_adm_GetFileList(SubFolder, ['^.*' x.SUBJECTS{iSubject} '.*\.(jpg|png|tiff|bmp)$'],'FPListRec',[0 Inf]);
        for iJ=1:length(JPGlist)
            [Fpath, Ffile, Fext] = fileparts(JPGlist{iJ});
            [StartInd, EndInd] = regexp(Ffile, x.SUBJECTS{iSubject});
            if StartInd==1
                NameWithoutID = Ffile(EndInd+2:end);
            else
                NameWithoutID = [Ffile(1:StartInd-2) Ffile(EndInd+2:end)];
            end

            if  isempty(NameWithoutID)
                NameWithoutID = Ffile;
            end

            NameWithoutID = xASL_adm_CorrectName(NameWithoutID,1);
            if NameWithoutID(1) == '_'
                % Tests if first character is illegal
                NameWithoutID = NameWithoutID(2:end);
            end

            SubFolder = xASL_adm_GetSubFolder(Fpath, x.D.ROOT);

            wadqc.qc_items.(NameWithoutID).type = 'file';
            wadqc.qc_items.(NameWithoutID).sub_path = SubFolder;
            wadqc.qc_items.(NameWithoutID).filename = [Ffile Fext];
            wadqc.qc_items.(NameWithoutID).category = 'object';
        end
    end

end




%% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
function [FolderPath] = xASL_adm_GetSubFolder(Fpath, RootIn)
%xASL_adm_GetSubFolder Obtain folder name by subtracting root folder & forcing forward slash

    FolderPath = Fpath(length(RootIn)+2:end);
    IndexSlash = strfind(FolderPath,'\');
    if ~isempty(IndexSlash)
        FolderPath(IndexSlash) = '/';
    end

end

