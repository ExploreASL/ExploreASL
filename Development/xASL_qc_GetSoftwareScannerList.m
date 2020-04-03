function [Sequence] = xASL_qc_GetSoftwareScannerList(RootIn, subject_regexp)
%xASL_qc_GetSoftwareScannerList Read JSON files for each NIfTI
% & compose a "sequence identifier" by concatenating
%[Manufacturer ?_? ManufacturersModelName DeviceSerialNumber ?_? SoftwareVersion]
% example: subject_regexp = '^OAS\d*_\d*$';


DirList = xASL_adm_GetFileList(RootIn, subject_regexp, 'List', [0 Inf], true);

Sequence = '';

fprintf('Collecting ASL sequences:...');
for iDir=1:length(DirList)
    xASL_TrackProgress(iDir, length(DirList));
    Sequence{end+1,1} = DirList{iDir};
    ASLList = xASL_adm_GetFileList(fullfile(RootIn, DirList{iDir}), '^ASL_\d*$', 'List', [0 Inf], true);
    for iASL = 1:length(ASLList)
        Sequence{end,2} = ['ASL_' num2str(iASL)];
        Sequence{end,3} = 'n/a';
        ASLdir = fullfile(RootIn, DirList{iDir},['ASL_' num2str(iASL)]);
        PathJSON = fullfile(ASLdir, 'ASL4D.json');
        if exist(PathJSON, 'file')
            CurrentSequence = '';
            try
                json = spm_jsonread(PathJSON);
                FieldsAre = {'Manufacturer' 'ManufacturersModelName' 'DeviceSerialNumber' 'SoftwareVersions'};
                for iField=1:length(FieldsAre)
                    if isfield(json,FieldsAre{iField})
                        if iField==1
                            CurrentSequence = json.(FieldsAre{iField});
                        else
                            CurrentSequence = [CurrentSequence '_' json.(FieldsAre{iField})];
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
end

MatPath = fullfile(RootIn, 'Sequence.mat');
save(MatPath,'Sequence','-mat');


end
