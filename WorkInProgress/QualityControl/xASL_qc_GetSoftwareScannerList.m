function [Sequence] = xASL_qc_GetSoftwareScannerList(RootIn, subjectRegexp, ScanType)
%xASL_qc_GetSoftwareScannerList Read JSON files for each NIfTI
% & compose a "sequence identifier" by concatenating
%[Manufacturer ?_? ManufacturersModelName DeviceSerialNumber ?_? SoftwareVersion]
% example: subjectRegexp = '^OAS\d*_\d*$';

% ScanType = ASL, T1 (default), FLAIR, M0

if nargin<2 || isempty(subjectRegexp)
    warning('Subject regular expression missing, defaulting to anything ".*"');
    subjectRegexp = '.*';
end
if nargin<3 || isempty(ScanType)
    ScanType = 'T1';
end
    

DirList = xASL_adm_GetFileList(RootIn, subjectRegexp, 'List', [0 Inf], true);

Sequence = '';

fprintf('Collecting ASL sequences:...');
for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    Sequence{end+1,1} = DirList{iDir};
    switch ScanType
        case 'T1'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^T1.*\.json$', 'FPList');
        case 'FLAIR'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^FLAIR.*\.json$', 'FPList');
		case 'T2'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^T2.*\.json$', 'FPList');
		case 'T1c'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^T1c.*\.json$', 'FPList');
        case 'ASL'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^ASL4D.*\.json$', 'FPListRec');
        case 'M0'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^M0.*\.json$', 'FPListRec');
        case 'dwi'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}, 'dwi'), '^dwi.*dwi\.json$', 'FPList');
        case 'func'
            SequenceList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}, 'func'), '^func.*bold\.json$', 'FPList');
        otherwise 
            warning('Undefined ScanType, quitting');
            return;
    end
    
    Sequence{end,2} = ScanType;
    Sequence{end,3} = 'n/a';
    if length(SequenceList)>0
        CurrentSequence = '';
        try
            json = xASL_io_ReadJson(SequenceList{1});
            FieldsAre = {'MagneticFieldStrength' 'Manufacturer' 'ManufacturersModelName' 'DeviceSerialNumber' 'SoftwareVersions'};
            for iField=1:length(FieldsAre)
                if isfield(json,FieldsAre{iField})
                    if iField==1
                        CurrentSequence = [xASL_num2str(json.(FieldsAre{iField})) 'T'];
                    else
                        CurrentSequence = [CurrentSequence '_' xASL_num2str(json.(FieldsAre{iField}))];
                    end
                end
            end
            if ~isempty(CurrentSequence)
                Sequence{end,3} = CurrentSequence;
            end
        catch ME
            warning(ME.message);
        end
    end
end

MatPath = fullfile(RootIn, 'Sequence.mat');
save(MatPath,'Sequence','-mat');


end
