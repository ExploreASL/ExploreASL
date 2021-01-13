function xASL_bids_BIDS_xASL(Dir_RawData, bOverwrite)
%xASL_bids_BIDS_xASL Convert BIDS rawdata to derivatives/ExploreASL
% Can be updated step-by-step when ExploreASL's derivative structure moves to BIDS
% NB: ask how Visits/session layer is defined in bids-matlab (should be
% separate layer within subjects, but now isn't?)
%
% PM: If multiple runs, this is fine. If multiple NIfTIs per run, issue a
% warning
% PM: could move folder/name configuration to JSON file
% PM: Add TopUp compatibility

FolderNameConfiguration =  {'type'   'modality' 'foldernames' 'filenames'      'run_locations'   'run_1_index';...
                            'T1w'    'anat'     ''            'T1'             'file'            false;...
                            'FLAIR'  'anat'     ''            'FLAIR'          'file'            false;...
                            'asl'    'perf'     'ASL'         'ASL4D'          'folder'          true;...
                            'm0scan' 'perf'     'ASL'         'M0'             'folder'          true;};

types = FolderNameConfiguration(2:end, 1);
modalities = FolderNameConfiguration(2:end, 2);
foldernames = FolderNameConfiguration(2:end, 3);
filenames = FolderNameConfiguration(2:end, 4);
run_locations = FolderNameConfiguration(2:end, 5);
run_1_index = FolderNameConfiguration(2:end, 6);

% Sidecars definition
Sidecars = {'.json' '_aslcontext.tsv' '_labeling.jpg'};
SidecarRequired =[1 0 0];
SidecarSuffixType = [1 0 0]; % specifies if the sidecar suffix keeps the scantype (e.g. yes for *_asl.json, not for *_aslcontext

%% ------------------------------------------------------------------------------------
%% Admin
if nargin<2 || isempty(bOverwrite)
    bOverwrite = 1;
end

Dir_RawData = '/Users/henk/ExploreASL/ASL/TestBIDS/GE_PCASL_3Dspiral/rawdata';

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
%% 1. Parse folder using bids-matlab
BIDS = bids.layout(Dir_RawData);

%% ------------------------------------------------------------------------------------
%% 2. Clone /rawdata as parsed by BIDS to /derivatives/ExploreASL 


% -----------------------------------------------
% 2a. Define Subject
nSubjects = length(BIDS.subjects);
for iSubject=1:nSubjects % iterate over subjects
    SubjectID = BIDS.subjects(iSubject).name;
    % Currently, ExploreASL concatenates subject_visit/timepoint in the
    % same folder layer, so we only use SubjectSession

    % -----------------------------------------------
    % 2b. Define SubjectVisit    
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

        % -----------------------------------------------
        % 2c. Parse modality
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
                
                % -----------------------------------------------
                % 2d. Parse scantype
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

                            % -----------------------------------------------
                            % 2e. Compile paths for copying                            
                            if length(TypeRunIndex)==1 % if this scantype-run combination exists

                                % ModalityIs = current modality (e.g. 'anat' 'perf')
                                % TypeIs = current scantype, e.g. 'asl' 'm0' 't1w'
                                % RunIs = current run (e.g. 1, 2, 3)
                                % TypeRunIndex = index of current scantype & run inside the above created Reference Table
                                
                                % Define folder & filename
                                ConfigIndex = find(cellfun(@(y) strcmp(TypeIs, y), types));
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
                                Path_Dest{1} = fullfile(Dir_xASL_SubjectVisit, FolderIs, [FileIs ModalityFields(TypeRunIndex).ext]);
                                
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
                                    
                                    if ~exist(TempSidecar, 'file') && SidecarRequired(iCar)
                                        warning([TempSidecar ' missing']);
                                    elseif exist(TempSidecar, 'file')
                                        Path_Orig{iCount+1} = TempSidecar;

                                        [Fpath, Ffile] = xASL_fileparts(Path_Dest{1});
                                        Path_Dest{iCount+1} = fullfile(Fpath, [Ffile Sidecars{iCar}]);
                                        iCount = iCount+1;
                                    end
                                end
                                    
                                % -----------------------------------------------
                                % 2f. Copy files
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