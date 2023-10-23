function x = xASL_qc_CollectParameters(x, iSubject, ScanType, iSession)
%xASL_qc_CollectParameters Collect all parameters from structural & ASL, spread over the derivative folders
%
% FORMAT: x = xASL_qc_CollectParameters(x, iSubject, ScanType, CollectQCFunction [, iSession])
%
% INPUT:
%   x                   - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject            - index of current subject (REQUIRED)
%   ScanType            - string for ScanType, options = 'Structural' 'ASL' 'func' 'dwi' (REQUIRED)
%   iSession            - index of current session (OPTIONAL, default 1, but RECOMMENDED for ScanType == 'ASL')
%
% OUTPUT:
%   x                   - same as input
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function collects QC parameters for a module for a given subject and for ASL for a given subject+session
% 
% EXAMPLE: x = xASL_qc_CollectParameters(x, 10, 'func');
%          x = xASL_qc_CollectParameters(x, 10, 'ASL', 4);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


%% -----------------------------------------------------------------------------------------------
%% Admin
fclose all;

fprintf('Collecting QC parameters...\n');

if nargin<3 || isempty(ScanType)
    % Run collecting and saving for all
    error('Missing ScanType...');
end
if nargin<2 || isempty(iSubject)
    % Run collecting and saving for all
    error('Missing iSubject...');
end

% It is recommended to provide session number for ASL
if nargin<4 || isempty(iSession)
	iSession = 1;
	if strcmp(ScanType, 'ASL')
		warning('Missing iSession for ASL, setting to 1');
	end
end

if ~isfield(x,'Output')
    x.Output  = struct; 
end
if ~isfield(x.Output, ScanType)
    x.Output.(ScanType) = struct; 
end

%% -----------------------------------------------------------------------------------------------
%% Collect subject-specific (i.e. structural/anatomical) QC results

switch ScanType
    case 'Structural'
        x = xASL_qc_CollectQC_Structural(x, iSubject);
    case 'ASL'
		% Check for a session subfield and create when necessary
		if ~isfield(x.Output.ASL, x.SESSIONS{iSession})
			x.Output.ASL.(x.SESSIONS{iSession}) = struct;
		end
		x = xASL_qc_CollectQC_ASL(x, iSubject, iSession);
    case 'func'
		x = xASL_qc_CollectQC_func(x, iSubject);
    case 'dwi'
end
fprintf('\n');

if ~isfield(x,'DoWADQCDC')
    x.DoWADQCDC = false; % default
end

%% Remove empty fields
FN = fieldnames(x.Output);
TF = cellfun(@(c) isempty(fields(x.Output.(c))), FN);
x.Output = rmfield(x.Output, FN(TF));

%% Save QC data

% Backward compatibility (used to store multiple subjects in single x.mat)
FieldsAre = fields(x.Output);
for iField=1:length(FieldsAre)
    LengthFields = length(x.Output.(FieldsAre{iField}));
    if LengthFields>1 && iSubject<=LengthFields
        x.Output.(FieldsAre{iField}) = x.Output.(FieldsAre{iField})(iSubject);
    end
end

%% -----------------------------------------------------------------------------------------------
%% Query current software versions    
x = xASL_qc_CollectSoftwareVersions(x);

% Module-specific software versions:
switch ScanType
    case 'Structural'
        x.Output.(ScanType).Version_CAT12 = x.Output.SoftwareVersion.CAT12;
        x.Output.(ScanType).Version_LST = x.Output.SoftwareVersion.LST;
        
        % General software versions (put this in structural only, to avoid redundant output
        x.Output.(ScanType).Version_ExploreASL = x.Output.SoftwareVersion.ExploreASL;
        x.Output.(ScanType).Version_Matlab = x.Output.SoftwareVersion.Matlab;
        x.Output.(ScanType).Version_SPM12 = x.Output.SoftwareVersion.SPM12;        
    case {'dwi', 'func'}
        x.Output.(ScanType).Version_FSL = x.Output.SoftwareVersion.FSL;
    case {'ASL'}
        x.Output.ASL.(x.SESSIONS{iSession}).Version_FSL = x.Output.SoftwareVersion.FSL;
end
% now remove the SoftwareVersion field to avoid redundancy
x.Output = rmfield(x.Output,'SoftwareVersion');

%% -----------------------------------------------------------------------------------------------
%% Save QC output
QC_Path = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, ['QC_collection_' x.SUBJECTS{iSubject} '.json']);

% Load QC.json
if xASL_exist(QC_Path, 'file')
	oldOutput = xASL_io_ReadJson(QC_Path);

	% Delete the old QC.json
	xASL_delete(QC_Path);

	if strcmp(ScanType, 'ASL') && isfield(oldOutput, 'ASL')
		% Copy QC of other ASL sessions, but not the current one to the current QC
		listFields = fieldnames(oldOutput.ASL);
		for iField = 1:length(listFields)
			if ~isempty(regexp(listFields{iField}, 'ASL_\d+', 'once')) && ~strcmp(listFields{iField}, x.SESSIONS{iSession})
				x.Output.ASL.(listFields{iField}) = oldOutput.ASL.(listFields{iField});
			end
		end
	end
end

% Save current QC
xASL_io_WriteJson(QC_Path, x.Output);

% Generate WAD-QC Descriptor 
% Run once for each subject
xASL_qc_WADQC_GenerateDescriptor(x, iSubject, ScanType); % skipped when ~x.DoWADQCDC


end