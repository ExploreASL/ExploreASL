function EPAD_ReportMissingFiles(ROOT, DeleteNII)
%EPAD_ReportMissingFiles Creates report of missing files

if nargin<2 || isempty(DeleteNII)
    DeleteNII = false;
end    

AnalysisDir = fullfile(ROOT,'analysis');
RawDir = fullfile(ROOT,'raw');
SavePath = fullfile(AnalysisDir, 'ReportMissingRawNIfTI.csv');

%% Load the ScanType Labels Configuration
[~, ScanTypeConfig] = xASL_adm_csv2tsv(fullfile(RawDir, 'ScanType_LabelsConfig.csv'), false, false);
TypeList = ScanTypeConfig(2:end,4); % ExploreASL names of the scantypes
WhichRows = cellfun(@(x) ~isempty(x), TypeList); % which of the TypeLists is not empty
ScanTypeConfig = ScanTypeConfig(logical([1; WhichRows]),:);

TypeList = ScanTypeConfig(2:end,4);
SliceN = ScanTypeConfig(2:end,5:end); % HERE WE ASSUME THAT THE 5TH COLUMN CONTAINS THE FIRST SITE
SiteN = ScanTypeConfig(1,5:end);

% Now show for each ExploreASL tag in ScanTypeConfig (fieldnames here) what
% regular expression we use to search
Files2Check.anat_T2star = {'anat_T2star'};
Files2Check.FLAIR = {'FLAIR'};
Files2Check.anat_T2w = {'anat_T2w'};
Files2Check.swi_part_mag = {'swi.*(part_mag|part_phase)'};
Files2Check.func_RevPE = {'func.*RevPE'};
Files2Check.func_NormPE = {'func.*NormPE'};
Files2Check.func_bold = {'func.*bold'};
Files2Check.dwi_RevPE = {'dwi.*RevPE'}; % we create the NormPE ourselves
Files2Check.dwi = {'dwi(?!.*(ADC|RevPE))'};
Files2Check.ASL4D_RevPE = {'ASL4D.*RevPE'};
Files2Check.M0 = {'M0'};
Files2Check.ASL4D = {'ASL4D(?!.*RevPE)'};
Files2Check.T1 = {'T1'};

Fields2Check = fields(Files2Check);
nFields = length(Fields2Check);


%% Write the missing NIfTIs
fclose all;
xASL_delete(SavePath);
FID = fopen(SavePath,'wt'); % declare the file

% Print explanation
fprintf(FID, 'This table shows any missing scans, sorted per site and ScanType.\n');
fprintf(FID, 'Sites and ScanTypes are printed repeatedly for those that we expect to be there according to the ScanTypeConfig table\n');

fprintf(FID, 'Subject,Complete');

DeleteList = '';
NonDeletedList = '';

for iSite=1:length(SiteN)
    SiteCurrent = SiteN{iSite}(5:7);
    
    fprintf(FID, ['\n\n\n' SiteCurrent]);
    fprintf(['\nProcessing ' SiteCurrent ':   ']);
    SubjDirs = xASL_adm_GetFileList(AnalysisDir, ['^' SiteCurrent 'EPAD\d*$'], 'FPList', [0 Inf], true);
    for iT=1:nFields % for ScanTypes
    	xASL_TrackProgress(iT,nFields);
        % first check if this ScanType should exist
        if isempty(SliceN{iT,iSite})
            SliceN{iT,iSite} = 1;
        end
        
        String2Check = EPAD_RemoveExistingRegExp(SliceN{iT,iSite},'[');
        String2Check = EPAD_RemoveExistingRegExp(String2Check,']');
        String2Check = EPAD_RemoveExistingRegExp(String2Check,'_');
        
        if ~xASL_str2num(String2Check)==0
            % find scantype in Files2Check table
            fprintf(FID, ['\n' Fields2Check{iT}]);
            
            iField = find(cellfun(@(x) strcmp(x,TypeList{iT}), Fields2Check));
            RegExp2Check = Files2Check.(Fields2Check{iField});
            ExpectedN = length(RegExp2Check);
            
            for iS=1:length(SubjDirs) % for all subjects
                [~, SubjName] = fileparts(SubjDirs{iS});
                FoundN = 0;
                for iR=1:ExpectedN
                    if length(xASL_adm_GetFileList(SubjDirs(iS), ['.*' RegExp2Check{iR} '.*\.nii$'], 'FPListRec', [0 Inf], false))>0
                        FoundN = FoundN+1;
                    end
                end
                if FoundN<ExpectedN
                    fprintf(FID, [SubjName ',' xASL_num2str(FoundN/ExpectedN) '\n']);
                    DeleteList{end+1} = SubjName;
                end
            end
        end
    end
end
    
fclose(FID);

if DeleteNII
    warning('Deleting subjects that were incomplete');
    fprintf('Wait a minute:   ');
    for iD=1:length(DeleteList)
        xASL_TrackProgress(iD,length(DeleteList));
        try
            SubFolder = fullfile(AnalysisDir, DeleteList{iD});
            if exist(SubFolder,'dir')
                xASL_adm_DeleteFileList(SubFolder, '.*', true, [0 Inf]);
                xASL_delete(SubFolder);
            end
        catch
            NonDeletedList{end+1} = SubFolder;
        end
    end
end

if ~isempty(NonDeletedList)
    warning('Not all files could be deleted:');
    NonDeletedList
end

end