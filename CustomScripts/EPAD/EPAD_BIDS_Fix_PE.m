function EPAD_BIDS_Fix_PE(AnalysisDir)
%EPAD_BIDS_Fix_PE Manage the PE NIfTI files
%
% FORMAT: EPAD_BIDS_Fix_PE(AnalysisDir)
% 
% INPUT:
%   AnalysisDir - path to folder containing the NIfTI/BIDS data (e.g. /data/RAD/share/EPAD/analysis) (REQUIRED)
%
% OUTPUT: n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function manages the phase encoding (PE) NifTI files,
%              which are used for TopUp geometric distortion correction. It
%              does the following:
%              1) Move all files into a single folder per ScanType: e.g. 
%                 'func_bold_1' 'func_NormPE_1' 'func_RevPE_1' -> 'func'
%                 'dwi_1' 'dwi_RevPE' -> 'dwi'
%                 Note that ASL is not here, as its M0 & RevPE files
%                 already go to the correct folder.
%              2) File renaming:
%                 If there are .mat files, rename them to _parms.mat
%                 If there are .dcm files, rename them to DummyDicom
%              3) Delete Dummy DICOMs from the ancillary (e.g. ADC, NormPE, RevPE) scans
%              This function should be run after EPAD_BIDS_Fix_DTI
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: EPAD_BIDS_Fix_PE(AnalysisDir);
% __________________________________
% Copyright 2015-2019 ExploreASL


SubjectList = xASL_adm_GetFsList(AnalysisDir,'^\d{3}EPAD\d*$', true, [], [], [0 Inf]);

fprintf('%s','Adjusting DTI BIDS:  0%');

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    
    % Move the //dwi_1 content to //dwi
    InputDir = fullfile(AnalysisDir, SubjectList{iSubject}, 'dwi_1');
    OutputDir = fullfile(AnalysisDir, SubjectList{iSubject}, 'dwi');
    FileList = xASL_adm_GetFileList(InputDir, 'dwi.*(\.bval|\.bvec|\.json|\.dcm|_parms\.mat|\.(nii|nii\.gz))$', 'FPList', [0 Inf]);
    for iFile=1:length(FileList)
        PathOri = FileList{iFile};
        [~, FileOri, ExtOri] = xASL_fileparts(PathOri);
        
        % Get BIDS run index % SAME AS ABOVE
        [Ind1, Ind2] = regexp(FileOri,'run(-|_)\d');
        RunString = FileOri(Ind1:Ind2);
        if isempty(RunString)
            RunString = 'run-1'; % default
        end
        RunString(4) = '-';
    
        [Ind1, Ind2] = regexp(FileOri,'(RevPE|NormPE|ADC)');
        if ~isempty(Ind1) && ~isempty(Ind2)
            ContrastString = FileOri(Ind1:Ind2);
        else
            ContrastString = 'dwi';
        end
        
        FileDest = ['dwi_' RunString '_' ContrastString];
        if strcmp(ExtOri,'.mat')
            PathDest = fullfile(OutputDir, [FileDest '_parms.mat']);
        elseif strcmp(ExtOri,'.dcm')
            PathDest = fullfile(OutputDir, ['DummyDicom_' FileDest ExtOri]);
        else
            PathDest = fullfile(OutputDir, [FileDest ExtOri]);
        end
        
        xASL_adm_CreateDir(OutputDir);
        xASL_Move(PathOri, PathDest);
    end
        
    % Delete DICOMs from ancillary scans
    DICOMlist = xASL_adm_GetFileList(fullfile(AnalysisDir, SubjectList{iSubject}), '^DummyDicom.*(RevPE|NormPE|ADC).*\.dcm$', 'FPListRec', [0 Inf]);
    for iFile=1:length(DICOMlist); xASL_delete(DICOMlist{iFile}); end
    
    if exist(InputDir,'dir') && isempty(xASL_adm_GetFileList(InputDir, '.*', 'FPListRec', [0 Inf]))
        rmdir(InputDir);
    end
    clear PathOri FileOri ExtOri Ind1 Ind2 RunString ContrastString FileDest PathDest
    clear FileList
end

% The part below is very similar to the part above, can be merged later

Extensions = {'\.(nii|nii\.gz)' '\.json' '_parms\.mat' '\.dcm'};

ScanType = {'func' 'dwi'};

ScanDirs{1} = {'func_bold_1' 'func_NormPE_1' 'func_RevPE_1'};
ScanDirs{2} = {'dwi_RevPE_1'};

ContrastString{1} = {'bold' 'NormPE' 'RevPE'};
ContrastString{2} = {'RevPE'}; % 'dwi' 'NormPE' 

RegExp = {'func(|_bold|_NormPE|_RevPE)' 'dwi_RevPE'}; % (|_dwi|_NormPE| )

fprintf('\n');
fprintf('%s','Adjusting func & DTI BIDS:  0%');

for iSubject=1:length(SubjectList)
    xASL_TrackProgress(iSubject,length(SubjectList));
    
    for iScanType=1:length(ScanType)
        for iD=1:length(ScanDirs{iScanType})
            ScanDir{iD} = fullfile(AnalysisDir, SubjectList{iSubject}, ScanDirs{iScanType}{iD});
            for iExt=1:length(Extensions)
                FileList{iD}{iExt} = xASL_adm_GetFileList(ScanDir{iD}, [RegExp{iScanType} '.*' Extensions{iExt} '$'], 'FPList', [0 Inf]);
                 
                if ~isempty(FileList{iD}{iExt})
                    [DirOri, FileOri, ExtOri] = xASL_fileparts(FileList{iD}{iExt}{1});

                    % Get BIDS run index
                    [Ind1, Ind2] = regexp(FileOri,'run(-|_)\d');
                    RunString = FileOri(Ind1:Ind2);
                    if isempty(RunString)
                        RunString = 'run-1'; % default
                    end
                    RunString(4) = '-';

                    DirNew = fullfile(AnalysisDir, SubjectList{iSubject}, ScanType{iScanType});
                    DestName = [ScanType{iScanType} '_' RunString '_' ContrastString{iScanType}{iD}];

                    PathOri = fullfile(DirOri, [FileOri ExtOri]);
                    if strcmp(ExtOri,'.mat')
                        DestName = [DestName '_parms'];
                    elseif strcmp(ExtOri,'.dcm')
                        DestName = ['DummyDicom_' DestName];
                    else
                        
                    end
                    PathDest = fullfile(DirNew, [DestName ExtOri]);
                    
                    if exist(PathOri,'file')
                        xASL_adm_CreateDir(DirNew);
                        xASL_Move(PathOri, PathDest, true);
                    end
                end
            end
            
            % Delete DICOMs from ancillary scans
            DICOMlist = xASL_adm_GetFileList(fullfile(AnalysisDir, SubjectList{iSubject}), '^DummyDicom.*(RevPE|NormPE|ADC).*\.dcm$', 'FPListRec', [0 Inf]);
            for iFile=1:length(DICOMlist); xASL_delete(DICOMlist{iFile}); end            
            
            if exist(ScanDir{iD},'dir')
                % delete any residual empty folders
                if isempty(xASL_adm_GetFileList(ScanDir{iD}, '.*', 'FPListRec', [0 Inf], false)) % if no files
                    FolderList = xASL_adm_GetFileList(ScanDir{iD}, '.*', 'FPListRec', [0 Inf], true); % then delete subfolders
                    for iFile=1:length(FolderList)
                        rmdir(FolderList{iFile});
                    end
                    rmdir(ScanDir{iD});
                end
            end
        end
    end
end

fprintf('\n');

end