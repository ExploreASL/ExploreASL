function xASL_bids_BIDS2Legacy(Dir_RawData, bOverwrite)
%xASL_bids_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: xASL_bids_BIDS2xASL(Dir_RawData[, bOverwrite])
% 
% INPUT:
%   Dir_RawData - path to the folder containing the raw BIDS data (REQUIRED)
%   bOverwrite  - boolean, true for overwriting files (OPTIONAL, DEFAULT = true)
%   
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts BIDS rawdata (e.g. /StudyName/rawdata/) 
% to xASL legacy derivative format (e.g. /StudyName/derivatives/ExploreASL/)
%
% Can be updated step-by-step when ExploreASL's derivative structure moves to BIDS
% NB: ask how Visits/session layer is defined in bids-matlab (should be
% separate layer within subjects, but now isn't?)
%
% This function performs the following steps:
% 1. Parse a folder using bids-matlab
% 2. Define Subject
% 3. Define SubjectVisit
% 4. Parse modality
% 5. Parse scantype
% 6. Compile paths for copying
% 7. Manage sidecars to copy
% 8. Copy files
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_bids_BIDS2xASL('MyStudy/rawdata')
% __________________________________
% Copyright 2015-2021 ExploreASL



%% ------------------------------------------------------------------------------------
%% Configuration
% PM: This configuration section could be taken outside this script, e.g. to a JSON file
FolderNameConfiguration =  {'type'   'modality' 'foldernames' 'filenames'      'run_locations'   'run_1_index';...
                            'T1w'    'anat'     ''            'T1'             'file'            false;...
                            'FLAIR'  'anat'     ''            'FLAIR'          'file'            false;...
                            'asl'    'perf'     'ASL'         'ASL4D'          'folder'          true;...
                            'm0scan' 'perf'     'ASL'         'M0'             'folder'          true;...
                            'm0scan' 'fmap'     'ASL'         'M0_RevPE'       'folder'          true;};

types = FolderNameConfiguration(2:end, 1); % scantypes
modalities = FolderNameConfiguration(2:end, 2); % the BIDS domains of scantypes
foldernames = FolderNameConfiguration(2:end, 3); % xASL legacy subfolder name (empty means no subfolder)
filenames = FolderNameConfiguration(2:end, 4); % xASL legacy filename
run_locations = FolderNameConfiguration(2:end, 5); % xASL legacy location of run specification (e.g. T1_2.nii.gz vs ASL_2/ASL4D.nii.gz for file vs folder location)
run_1_index = FolderNameConfiguration(2:end, 6); % if xASL legacy requires to specify the first run (e.g. T1.nii.gz vs ASL_1/ASL4D.nii.gz)

% Sidecars definition
Sidecars = {'.json' '_aslcontext.tsv' '_labeling.jpg'};
SidecarRequired =[1 0 0];
SidecarTypeSpecific = {'no' 'asl' 'asl'};
SidecarSuffixType = [1 0 0]; % specifies if the sidecar suffix keeps the scantype (e.g. yes for *_asl.json, not for *_aslcontext


%% ------------------------------------------------------------------------------------
%% Admin
if nargin<2 || isempty(bOverwrite)
    bOverwrite = 1;
end

[~, Ffile] = xASL_fileparts(Dir_RawData);
if ~strcmp(Ffile, 'rawdata')
    warning('Invalid folder selected, needs to be rawdata folder');
    return;
end

Dir_xASL = fullfile(fileparts(Dir_RawData), 'derivatives', 'ExploreASL');
if exist(Dir_xASL, 'dir') && bOverwrite
    warning([Dir_xASL ' already existed, overwriting']);
elseif exist(Dir_xASL, 'dir')
    fprintf('%s\n', [Dir_xASL ' existed, merging']);
else
    xASL_adm_CreateDir(Dir_xASL);
end


%% ------------------------------------------------------------------------------------
%% 1. Parse a folder using bids-matlab
BIDS = bids.layout(Dir_RawData);

fprintf('%s\n', ['Converting ' Dir_RawData]);
fprintf('%s', ' -> ../derivatives/ExploreASL:   ');


%% -----------------------------------------------
%% 2. Define Subject
nSubjects = length(BIDS.subjects);
for iSubject=1:nSubjects % iterate over subjects
    xASL_TrackProgress(iSubject, nSubjects);
    SubjectID = BIDS.subjects(iSubject).name;
    % Currently, ExploreASL concatenates subject_visit/timepoint in the
    % same folder layer, so we only use SubjectSession

    %% -----------------------------------------------
    %% 3. Define SubjectVisit
    nVisits = max([1 length(BIDS.subjects(iSubject).session)]); % minimal 1 visit
    for iVisit=1:nVisits % iterate visits in this Subject
        % ExploreASL uses visit as a number (e.g. _1 _2 _3 etc)
        if nVisits==1
            Dir_xASL_SubjectVisit = fullfile(Dir_xASL, SubjectID);
            VisitString = '';
        else
            Dir_xASL_SubjectVisit = fullfile(Dir_xASL, [SubjectID '_' xASL_num2str(iVisit)]);
            VisitString = [' visit ' BIDS.subjects(iSubject).session];
        end
        SubjectVisit = [SubjectID VisitString];
        xASL_adm_CreateDir(Dir_xASL_SubjectVisit);

        %% -----------------------------------------------
        %% 4. Parse modality
        ModalitiesUnique = unique(modalities);
        nModalities = length(ModalitiesUnique);
        for iModality=1:nModalities % iterate modalities in this Subject/Visit
            ModalityIs = ModalitiesUnique{iModality};
            if isfield(BIDS.subjects(iSubject), ModalityIs) && ~isempty(BIDS.subjects(iSubject).(ModalityIs))
                ModalityFields = BIDS.subjects(iSubject).(ModalityIs);
                nScans = length(ModalityFields);

                % Parse fields of this modality for combinations ScanType & Run, in a reference table
                Reference = {'index' 'ScanType' 'run'};
                for iScan=1:nScans % iterate NIfTIs in this Subject/Visit/Modality
                    Reference{iScan+1, 1} = iScan; % index
                    Reference{iScan+1, 2} = ModalityFields(iScan).type; % ScanType
                    Reference{iScan+1, 3} = xASL_str2num(ModalityFields(iScan).run); % run
                    if isempty(Reference{iScan+1, 3}) || isnan(Reference{iScan+1, 3})
                        Reference{iScan+1, 3} = 1; % default to 1st run
                    end
                end
                Reference(2:end,:) = sortrows(Reference(2:end,:), [2, 3]); % first sort for ScanType then run

                RunsAre = cellfun(@(y) y, Reference(2:end, 3));
                RunsUnique = unique(RunsAre);
                
                %% -----------------------------------------------
                %% 5. Parse scantype
                for iType=1:length(types) % iterate scantypes in this Subject/Visit/Modality
                    TypeIs = types{iType};
                    TypeIndices = cellfun(@(y) strcmp(y, TypeIs), Reference(2:end, 2)); % this are the indices for this ScanType

                    if ~isempty(TypeIndices)
                        for iRun=1:length(RunsUnique) % iterate runs in this Subject/Visit/Modality
                            RunIs = RunsUnique(iRun);
                            RunIndices = RunsAre==RunsUnique(iRun);
                            TypeRunIndex = find(RunIndices & TypeIndices);
                            if length(TypeRunIndex)>1
                                warning(['Multiple NIfTIs found for ' SubjectVisit '_run-' xASL_num2str(RunIs) '_' TypeIs ', using first only']);
                                TypeRunIndex = TypeRunIndex(1);
                            end

                            %% -----------------------------------------------
                            %% 6. Compile paths for copying                          
                            if length(TypeRunIndex)==1 % if this scantype-run combination exists

                                % ModalityIs = current modality (e.g. 'anat' 'perf')
                                % TypeIs = current scantype, e.g. 'asl' 'm0' 't1w'
                                % RunIs = current run (e.g. 1, 2, 3)
                                % TypeRunIndex = index of current scantype & run inside the above created Reference Table
                                
                                % Define folder & filename
                                ConfigIndex = cellfun(@(y) strcmp(TypeIs, y), types); % find scantype index
                                ConfigIndex2 = cellfun(@(y) strcmp(ModalityIs, y), modalities); % find modality index
                                ConfigIndex = find(ConfigIndex & ConfigIndex2); % Combine them & convert to index
                                
                                FolderIs = foldernames{ConfigIndex};
                                FileIs = filenames{ConfigIndex};
                                
                                % Manage runs
                                PrintRun = false;
                                if RunIs==1 && run_1_index{ConfigIndex}==1
                                    % if we need to print first run
                                    PrintRun = true;
                                elseif RunIs>1
                                    PrintRun = true;
                                end
                                
                                if PrintRun==1
                                    if strcmp(run_locations{ConfigIndex}, 'file')
                                        FileIs = [FileIs '_' xASL_num2str(RunIs)];
                                    elseif strcmp(run_locations{ConfigIndex}, 'folder')
                                        FolderIs = [FolderIs '_' xASL_num2str(RunIs)];
                                    end
                                end
                                
                                Path_Orig{1} = fullfile(BIDS.subjects(iSubject).path, ModalityIs, ModalityFields(TypeRunIndex).filename);
                                [~, ~, Fext] = xASL_fileparts(ModalityFields(TypeRunIndex).filename);
                                Path_Dest{1} = fullfile(Dir_xASL_SubjectVisit, FolderIs, [FileIs Fext]);
                                
                                
                                %% -----------------------------------------------
                                %% 7. Manage sidecars to copy
                                
                                % Assuming that each .nii has a .json
                                % sidecar, do the same for .json (and for
                                % other sidecars only if they exist per
                                % SidecarRequired)
                                
                                iCount = 1;
                                for iCar=1:length(Sidecars)
                                    [Fpath, Ffile] = xASL_fileparts(Path_Orig{1});
                                    
                                    if ~SidecarSuffixType(iCar)
                                        Ffile = Ffile(1:end-length(TypeIs)-1);
                                    end
                                    TempSidecar = fullfile(Fpath, [Ffile Sidecars{iCar}]);
                                    
                                    if ~strcmp(SidecarTypeSpecific{iCar}, 'no') && ~strcmp(SidecarTypeSpecific{iCar}, TypeIs)
                                        % skip this sidecar (e.g. some
                                        % asl-specific sidecars for non-asl NIfTIs)
                                    elseif ~exist(TempSidecar, 'file') && SidecarRequired(iCar)
                                            warning([TempSidecar ' missing']);
                                    elseif exist(TempSidecar, 'file')
                                        Path_Orig{iCount+1} = TempSidecar;

                                        [Fpath, Ffile] = xASL_fileparts(Path_Dest{1});
                                        Path_Dest{iCount+1} = fullfile(Fpath, [Ffile Sidecars{iCar}]);
                                        iCount = iCount+1;
                                    end
                                end
                                    
                                %% -----------------------------------------------
                                %% 8. Copy files
                                for iFile=1:length(Path_Orig)
                                    xASL_bids_BIDS2xASL_CopyFile(Path_Orig{iFile}, Path_Dest{iFile}, bOverwrite);
                                end

                            end
                        end % iterate runs in this Subject/Visit/Modality
                    end
                end
            end
        end
    end
end

fprintf('\n');

%% Parse M0
ListASL4D = xASL_adm_GetFileList(Dir_xASL, '^ASL4D\.nii$', 'FPListRec');
if ~isempty(ListASL4D)
    for iList=1:numel(ListASL4D)
        xASL_bids_parseM0(ListASL4D{iList});
    end
    fprintf('%s\n', ['M0 parsed for ' ListASL4D{iList}]);
else
    warning(['No ASL4D file found in ' Dir_xASL]);
end

%% Create DataPar.json
PathDataPar = fullfile(Dir_xASL, 'DataPar.json');
JSON.x.subject_regexp = '^sub-.*$';
JSON.x.Quality = 0;
JSON.x.DELETETEMP = 1;
spm_jsonwrite(PathDataPar, JSON');


end


%% ===========================================================================
function xASL_bids_BIDS2xASL_CopyFile(Path_Orig, Path_Dest, bOverwrite)
%xASL_bids_BIDS2xASL_CopyFile

    % Create folder(s) if didnt exist
    xASL_adm_CreateDir(fileparts(Path_Dest));

    % check for existance & overwriting
    if ~xASL_exist(Path_Orig)
        warning(['Couldnt find ' Path_Orig, ' skipping']);
    elseif xASL_exist(Path_Dest) && ~bOverwrite
        warning([Path_Dest ' already existed, skipping']);
    else % if didnt already exist, or bOverwrite is true
        xASL_Copy(Path_Orig, Path_Dest, 1);
    end
    
    
end