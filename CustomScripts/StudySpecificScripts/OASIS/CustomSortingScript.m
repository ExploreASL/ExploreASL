%% Custom Sorting script

ExploreASL_Master('',0);

Odir = 'C:\BackupWork\ASL\OASIS\OASIS3';
Ddir = 'C:\BackupWork\ASL\OASIS\OASIS3_xASL';
xASL_adm_CreateDir(Ddir);
Slist1 = xASL_adm_GetFileList(Odir, '^OAS\d*', 'FPList',[0 Inf], true);

fprintf('Converting OASIS data to ExploreASL-compatible:   ');

for iS1=1:length(Slist1)
    xASL_TrackProgress(iS1, length(Slist1));
    MRlist = xASL_adm_GetFileList(Slist1{iS1}, '^OAS\d*_MR_.*', 'FPList',[0 Inf], true);
    for iMR=1:length(MRlist)
        if ~isempty(xASL_adm_GetFileList(MRlist{iMR}, '.*asl\.nii$', 'FPListRec',[0 Inf])) && ~isempty(xASL_adm_GetFileList(MRlist{iMR}, '.*T1w\.nii$', 'FPListRec',[0 Inf]))
            % if contains asl & T1w
            [~, SubjName] = fileparts(MRlist{iMR});
            [StartI, EndI] = regexp(SubjName, '^OAS\d*_');
            Name1 = SubjName(StartI:EndI-1);
            [StartI, EndI] = regexp(SubjName, '_d\d*');
            Name2 = SubjName(StartI+2:EndI);
            SubjName = [Name1 '_' Name2]; % name with longitudinal ID
            DestSubjDir = fullfile(Ddir, SubjName);
            xASL_adm_CreateDir(DestSubjDir);
            
            ScanList = xASL_adm_GetFileList(MRlist{iMR}, '^.*\d$', 'FPList',[0 Inf], true);
            for iScan=1:length(ScanList)
                NiftiDir = fullfile(ScanList{iScan}, 'NIFTI');
                BIDSDir = fullfile(ScanList{iScan}, 'BIDS');
                ScanType = {'asl' 'fieldmap' 'T1w', 'FLAIR'}; % DICOM name
                SubDirHas = [1 1 0 0];
                SubDirName = {'ASL_' 'FieldMap_' '' ''};
                FileName = {'ASL4D' 'FieldMap' 'T1' 'FLAIR'};
                for iType=1:length(ScanType)
                    NIIlist = xASL_adm_GetFileList(NiftiDir, ['.*' ScanType{iType} '\.nii$'], 'FPList',[0 Inf]);
                    if ~isempty(NIIlist)
                        for iNii=1:length(NIIlist)
                            clear DestDir DestFile
                            if SubDirHas(iType)
                                DestDir = fullfile(DestSubjDir, [SubDirName{iType} num2str(iNii)]);
                                DestFile = fullfile(DestDir, [FileName{iType} '.nii.gz']);
                            elseif ~SubDirHas(iType) && iNii==1
                                DestDir = DestSubjDir;
                                DestFile = fullfile(DestDir, [FileName{iType} '.nii.gz']);
                            elseif ~SubDirHas(iType) && iNii>1
                                DestDir = DestSubjDir;
                                DestFile = fullfile(DestDir, [FileName{iType} '_' num2str(iNii) '.nii.gz']);
                            else
                                warning(['Didnt know what to do for ' NIIlist{iNii}]);
                            end
                            
                            IndexIfExist=1;
                            while xASL_exist(DestFile, 'file')
                                IndexIfExist = IndexIfExist+1;
                                [~, Ffile] = xASL_fileparts(DestFile);
                                iStart = regexp(Ffile,'_\d*');
                                if ~isempty(iStart)
                                    Ffile = Ffile(1:iStart-1);
                                end
                                DestFile = fullfile(DestDir, [Ffile '_' num2str(IndexIfExist) Fext]);
                            end
                            
                            xASL_adm_CreateDir(DestDir);
                            xASL_Copy(NIIlist{iNii}, DestFile, true);
                            
                            % do the same for the JSON
                            JSONdir = fullfile(fileparts(fileparts(NIIlist{iNii})), 'BIDS');
                            JSON = xASL_adm_GetFileList(JSONdir, ['.*' ScanType{iType} '.*\.json$'], 'FPList', [0 Inf]);
                            if length(JSON)~=1
                                warning(['Too many or few jsons found in ' JSONdir]);
                            end
                            
                            [DestDir, Ffile] = xASL_fileparts(DestFile);
                            DestFile = fullfile(DestDir, [Ffile '.json']);

                            xASL_Copy(JSON{1}, DestFile, true);
                        end
                    end
                end
            end
        end
    end
end

end