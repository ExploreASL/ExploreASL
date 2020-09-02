function x = xASL_qc_CollectParameters(x, iSubject, ScanType)
%xASL_qc_CollectParameters Collect all parameters from structural & ASL, spread over the derivative folders
%
% FORMAT: x = xASL_qc_CollectParameters(x, iSubject, ScanType, CollectQCFunction)
%
% INPUT:
%   x                   - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject            - index of current subject (REQUIRED)
%   ScanType            - string for ScanType, options = 'Structural' 'ASL' 'func' 'dwi' (REQUIRED)
%
% OUTPUT:
%   x                   - same as input
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function collects QC parameters for a module
% 
% EXAMPLE: x = xASL_qc_CollectParameters_Structural(x, 10, 'ASL');
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


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

if ~isfield(x,'Output')
    x.Output  = struct; 
end
if ~isfield(x.Output,ScanType)
    x.Output.(ScanType) = struct; 
end

switch ScanType
    case 'Structural'
        QCCollectFunction = @xASL_qc_CollectQC_Structural;
    case 'ASL'
        QCCollectFunction = @xASL_qc_CollectQC_ASL;
    case 'func'
        QCCollectFunction = @xASL_qc_CollectQC_func;
    case 'dwi'
        QCCollectFunction = [];
end

if ~isfield(x,'DoWADQCDC')
    x.DoWADQCDC = false; % default
end



%% -----------------------------------------------------------------------------------------------
%% Collect subject-specific (i.e. structural/anatomical) QC results
if ~isempty(QCCollectFunction)
    x = QCCollectFunction(x, iSubject);
    fprintf('\n');
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
        if isdeployed
            x.Output.(ScanType).Version_Matlab = version;
            x.Output.(ScanType).Version_Matlab = x.Output.(ScanType).Version_Matlab(1:3);
        else
            x.Output.(ScanType).Version_Matlab = x.Output.SoftwareVersion.Matlab;
        end
        x.Output.(ScanType).Version_SPM12 = x.Output.SoftwareVersion.SPM12;        
    case {'ASL', 'dwi', 'func'}
        x.Output.(ScanType).Version_FSL = x.Output.SoftwareVersion.FSL;
end
% now remove the SoftwareVersion field to avoid redundancy
x.Output = rmfield(x.Output,'SoftwareVersion');

%% -----------------------------------------------------------------------------------------------
%% Save QC output
QC_Path = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, ['QC_collection_' x.SUBJECTS{iSubject} '.json']);
xASL_delete(QC_Path);
xASL_adm_SaveJSON(x.Output, QC_Path);

% Generate WAD-QC Descriptor 
% Run once for each subject
xASL_qc_WADQC_GenerateDescriptor(x, iSubject, ScanType); % skipped when ~x.DoWADQCDC


end