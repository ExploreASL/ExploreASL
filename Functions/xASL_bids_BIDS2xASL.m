function [ output_args ] = xASL_bids_BIDS_xASL(Dir_RawData, bOverwrite)
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
                            'T1w'    'anat'     ''            'T1'             'file'            'no';...
                            'FLAIR'  'anat'     ''            'FLAIR'          'file'            'no';...
                            'asl'    'perf'     'ASL'         'ASL4D'          'folder'          'yes';...
                            'm0scan' 'perf'     'ASL'         'M0'             'folder'          'yes';};

types = FolderNameConfiguration(2:end, 1);
modalities = FolderNameConfiguration(2:end, 2);
foldernames = FolderNameConfiguration(2:end, 3);
filenames = FolderNameConfiguration(2:end, 4);
run_locations = FolderNameConfiguration(2:end, 5);
run_1_index = FolderNameConfiguration(2:end, 6);
                        
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
                            TypeRunIndices = find(RunIndices & TypeIndices);
                            if length(TypeRunIndices)>1
                                warning(['Multiple NIfTIs found for ' SubjectVisit '_run-' xASL_num2str(RunIs) '_' TypeIs ', using first only']);
                                TypeRunIndices = TypeRunIndices(1);
                            end

                            if length(TypeRunIndices)==1 % if this scantype-run combination exists
                                % -----------------------------------------------
                                % 2e. Compile paths for copying
                                
                                TypeIs
                                RunIs
                                TypeRunIndices
                                ModalityIs
                                
                                % HERE WE COMPILE THE rawdata NIFTI PATH
                                % THEN WE COMPILE THE derivatives NIFTI PATH
                                % JSONs are always the full path with
                                % different extension, that we compile as
                                % well
                                % 
                                
                                    
                                    % -----------------------------------------------
                                    % 2f. Copy, check for overwriting, and
                                    % existance file
                                    % also
                                    % xASL_adm_CreateDir(fileparts(Path2Copy))

                                Path_Raw_Anat = fullfile(BIDS.subjects(iSubject).path, 'anat', ModalityFields(iScan).filename);
                                Path_xASL_Anat = fullfile(Dir_xASL_SubjectVisit, '

                                if ~xASL_exist(Path_Raw_Anat)
                                    warning(['Couldnt find ' Path_Raw_Anat, ' skipping']);
                                else
                                    xASL_Copy(


                        
                        
                        


end