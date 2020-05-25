function [Parms, x, Oldx] = xASL_adm_LoadParms(ParmsPath, x, bVerbose)
%xASL_adm_LoadParmsMat Loads parameters from .mat or .json sidecar
%
% FORMAT: [Parms, x, Oldx] = xASL_adm_LoadParms(ParmsPath[, x, bVerbose])
%
% INPUT:
%   ParmsPath   - path to sidecar (*.mat, *.json, OPTIONAL, DEFAULT=skip this)
%   x           - structure containing pipeline settings (OPTIONAL)
%   bVerbose    - boolean specifying whether verbose output on screen/log
%                 is desired (OPTIONAL, DEFAULT=true)
%
%   OUTPUT:
%   Parms       - Loaded parameters to be used in processing
%   x           - structure containing pipeline settings (OPTIONAL)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads the internal memory x struct, any
% legacy *_parms.mat sidecar, any *.json BIDS sidecar, to use scan-specific
% parameters for image processing/quantification. Also, per BIDS
% inheritance, any x.S.SetsID parameters (from participants.tsv) are loaded
% as well. This function performs the following steps:
%
% 1) Load .mat parameter file
% 2) Load JSON file
% 3) Deal with warnings
% 4) Find fields with scan-specific data in x.S.Sets, and use this if possible (per BIDS inheritance)
% 5) Sync Parms.* with x.(Q.)* (overwrite x/x.Q)
% 6) Fix M0 parameter if not set
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [~, x] = xASL_adm_LoadParms('/MyStudy/sub-001/ASL_1/ASL4D.json', x, bO);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% ------------------------------------------------------------------------
%% 0) Admin
Parms = struct; % default

if nargin<1 || isempty(ParmsPath)
    warning('ParmsPath was not specified');
    ParmsPath = '';
    Parms = struct;
end
if nargin<2 || isempty(x)
    warning('No x struct fields available, can be missing some information for quantification');
    x = struct; % dummy variable
	x.Q = [];
end
if nargin<3 || isempty(bVerbose)
    bVerbose = true;
end

[Fpath, Ffile, Fext] = fileparts(ParmsPath);

%% ------------------------------------------------------------------------
%% 1) Load .mat parameter file (if exists)
% if ~exist(ParmsPath, 'file') && isfield(x.Q,'ASL')
%     warning('No ParmsPath found, trying to use general parameters!');
%     
%     if ~isempty(regexp(Ffile,'ASL4D'))
%         if bVerbose; fprintf('Use generic ASL parameters\n'); end
%         ImType = 'ASL';
%     elseif ~isempty(regexp(Ffile,'M0'))
%         if bVerbose; fprintf('Use generic M0 parameters\n'); end
%         ImType = 'M0';
%     end
%     
%     FieldsGeneric = fields(x.Q.(ImType));
%     for iField=1:length(FieldsGeneric)
%         if ~isfield(x.Q,FieldsGeneric{iField})
%             x.Q.(FieldsGeneric{iField}) = x.Q.(ImType).(FieldsGeneric{iField});
%         end
%     end
    
if exist(ParmsPath, 'file') && strcmp(Fext,'.mat')
    Parms = load(ParmsPath,'-mat');
    if  isfield(Parms,'parms')
        Parms = Parms.parms;
    else
		warning(['Couldnt read parameter files from ' ParmsPath]);
    end
end

%% ------------------------------------------------------------------------
%% 2) Load JSON file (if exists)
% Define JSON path
% First try defining from input (parms.mat)
if ~isempty(ParmsPath)
    [Fpath, Ffile] = xASL_fileparts(ParmsPath);
	if ~isempty(strfind(Ffile(3:end),'_parms'))
		Ffile = Ffile(1:end-6);
	end
    JSONPath = fullfile(Fpath, [Ffile '.json']);
else
    if ~isempty(regexp(Ffile,'ASL4D'))
        JSONPath = fullfile(Fpath, 'ASL4D.json');
    elseif ~isempty(regexp(Ffile,'M0'))
        JSONPath = fullfile(Fpath, 'M0.json');
    elseif ~isempty(regexp(Ffile,'func'))
        JSONPath = xASL_adm_GetFileList(Fpath,'func.*bold\.json',  'FPList', [0 Inf]);
        if ~isempty(JSONPath)
            JSONPath = JSONPath{1};
        end
    elseif ~isempty(regexp(Ffile,'dwi'))
        JSONPath = xASL_adm_GetFileList(Fpath,'dwi.*dwi\.json',  'FPList', [0 Inf]);
        if ~isempty(JSONPath)
            JSONPath = JSONPath{1};
        end    
    else
        warning('Could not define JSON path');
    end
end

% Load JSON file
VendorRestore = false; % default
if exist(JSONPath,'file') % According to the BIDS inheritance principle, the JSON values overwrite the existing values
    % But this can have some exceptions, e.g. the specification of sequence
    % in x.Vendor
    
    if isfield(x,'Vendor') && (strcmp(x.Vendor,'GE_product') || strcmp(x.Vendor,'GE_WIP'))
        VendorBackup = x.Vendor;
        VendorRestore = true;
	end
   
	JSONParms = spm_jsonread(JSONPath);
	
	Parms = xASL_bids_parms2BIDS(Parms, JSONParms, 0, 1);
	
end
    

%% ------------------------------------------------------------------------
%% 3) Deal with warnings
if isempty(fields(Parms))
    warning('parms seem missing, something wrong with parmsfile?');
end

if VendorRestore
    x.Vendor = VendorBackup;
    Parms.Vendor = VendorBackup;
end


%% ------------------------------------------------------------------------
%% 4) Find fields with scan-specific data in x.S.Sets, and use this if possible (per BIDS inheritance)
% Note that x.S.Sets is filled with data from participants.tsv or e.g. qnt_T1a.mat in the analysis root folder

% Find current index
iSubject = find(strcmp(x.SUBJECTS, x.P.SubjectID));
iSession = find(strcmp(x.SESSIONS, x.P.SessionID));
iSubjSess = (iSubject-1)*x.nSessions + iSession;

if size(x.S.SetsID,1)~=iSubjSess
    warning('Inheritance x.S.SetsID data was not equal to numbers of Subjects/Sessions, skipping');
    return;
end

Sets2Find  = {'qnt_ATT' 'qnt_T1a' 'qnt_lab_eff'         'LabelingEfficiency' 'Hematocrit'};
FieldNames = {'ATT'     'BloodT1' 'LabelingEfficiency'  'LabelingEfficiency' 'Hematocrit'};

for iSet=1:length(Sets2Find)
    TempIndex = find(cellfun(@(x) strcmp(x, Sets2Find{iSet}), x.S.SetsName));
	if ~isempty(TempIndex)
		SetIndex(iSet) = TempIndex;
	else
		SetIndex(iSet) = NaN;
    end

    if ~isnan(SetIndex(iSet))
        % Use the data out SetsID
        Parms.(FieldNames{iSet}) = x.S.SetsID(iSubjSess, SetIndex(iSet));
        if bVerbose; fprintf('%s\n', ['Loaded ' FieldNames{iSet} ': ' xASL_num2str(Parms.(FieldNames{iSet}))]); end
        % But check if this is the true data content, or if this is an index (e.g. 1, 2, 3, 4)
        % If its not continuous (x.S.Sets_1_2Sample~=3), then ExploreASL believes that this is an ordinal data set (groups)
        if x.S.Sets1_2Sample(SetIndex(iSet))~=3
            % if data are saved as indices of x.S.SetsOptions, then use
            % the SetsOption field/ID/name
            if Parms.(FieldNames{iSet})<=length(x.S.SetsOptions{SetIndex(iSet)})
                Parms.(FieldNames{iSet}) = x.S.SetsOptions{SetIndex(iSet)}{Parms.(FieldNames{iSet})};
            end
        end
        if ischar(Parms.(FieldNames{iSet})) % convert string to float
            Parms.(FieldNames{iSet}) = str2num(Parms.(FieldNames{iSet}));
        end
    end
end


%% ------------------------------------------------------------------------
%% 5) Sync Parms.* with x.(Q.)* (overwrite x/x.Q)
% Input all fields from the Parms into the x structure, backup those that were already existing (inheritance principle)
% & backward compatibility
Oldx = struct;
ParmsFieldNames = fieldnames(Parms);
xFieldNames = fieldnames(x);

for iPar=1:length(ParmsFieldNames)
    % Check first if already exists, to backup
    for iS=1:length(xFieldNames)
        if strcmp(ParmsFieldNames{iPar},xFieldNames{iS})
            Oldx.(xFieldNames{iS}) = x.(xFieldNames{iS}); % backup
        end
    end
    x.(ParmsFieldNames{iPar}) = Parms.(ParmsFieldNames{iPar});
end

% Move quantification parameters to the Q (quantification) subfield, for
% backward compatibility
Qfields = {'BackGrSupprPulses' 'LabelingType' 'Initial_PLD' 'LabelingDuration' 'SliceReadoutTime' 'Lambda'...
    'T2art' 'BloodT1' 'TissueT1' 'nCompartments' 'NumberOfAverages'};
for iField=1:length(Qfields)
    if isfield(x,Qfields{iField})
        if isfield(x.Q,(Qfields{iField})) && ~min((x.Q.(Qfields{iField})==x.(Qfields{iField})))
            if bVerbose; warning(['Overwriting x.Q.' Qfields{iField} '=' xASL_num2str(x.Q.(Qfields{iField})) ', with x.' Qfields{iField} '=' xASL_num2str(x.(Qfields{iField}))]); end
        end
        
        x.Q.(Qfields{iField}) = x.(Qfields{iField});
        x = rmfield(x, Qfields{iField});
    end
end


%% ------------------------------------------------------------------------
%% 6) Fix M0 parameter if not set
if ~isfield(x,'M0')
    if xASL_exist(fullfile(Fpath, 'M0.nii'),'file') && (exist(fullfile(Fpath, 'M0.json'),'file') || exist(fullfile(Fpath, 'M0_parms.mat'),'file') )
        x.M0 = 'separate_scan';
        if bVerbose; fprintf('%s\n',['M0 parameter was missing, set to ' x.M0]); end
    elseif isfield(Parms,'BackGrSupprPulses') && Parms.BackGrSupprPulses==0
        x.M0 = 'UseControlAsM0';
        if bVerbose; fprintf('%s\n',['M0 parameter was missing, set to ' x.M0]); end
    else
        if bVerbose; fprintf('%s\n','M0 parameter was missing, OR didnt find M0 scan, AND BackGrSupprPulses wasnt set to 0...'); end
    end
end

if ~exist('Parms','var')
    Parms = struct; 
    if bVerbose; warning('parms seem missing, something wrong with parmsfile?'); end
end



end