EPAD_CheckSameScannerSameSoftware

%% Admin

AnalysisDir = fullfile(ROOT,'analysis');
RawDir = fullfile(ROOT,'raw');
SavePath = fullfile(AnalysisDir, 'EPAD_CheckSameScannerSameSoftware.csv');

%% Store information in CSV file
fclose all;
xASL_delete(SavePath);
FID = fopen(SavePath,'wt'); % declare the file

JSONfields = {'MagneticFieldStrength' 'Manufacturer' 'ManufacturersModelName' 'DeviceSerialNumber' 'SoftwareVersions' 'CoilString'};
fprintf(FID, 'Subject');
for iJSON=1:length(JSONfields)
    fprintf(FID, [JSONfields{iJSON} ',']);
end
fprintf(FID, '\n');

%% Get list of scanners
Dlist = xASL_adm_GetFileList(AnalysisDir, '^\d{3}EPAD\d*$', 'List', [], true);
for iDir=1:length(Dlist)
    ScannerID{iDir} = Dlist{iDir}(1:3);
end
ScannerID = unique(ScannerID);

%% Loop over scanners
fprintf('Printing JSON fields to check same scanner same software:   ');
for iScanner=1:length(ScannerID)
    xASL_TrackProgress(iScanner, length(ScannerID));
    fprintf(FID, ['Scanner' num2str(ScannerID{iScanner}) '\n']);
    Dlist = xASL_adm_GetFileList(AnalysisDir, ['^' ScannerID{iScanner} 'EPAD\d*$'], 'FPList', [], true);

    for iDir=1:length(Dlist)
        % Print subjectID
        [~, Ffile] = fileparts(Dlist{iDir});
        fprintf(FID, [Ffile ',']);
        % Print JSON fields (if there are any)
        JSONlist = xASL_adm_GetFileList(Dlist{iDir}, '^(?!(QC|WAD)).*\.json$', 'FPListRec', [0 Inf]);
        if ~isempty(JSONlist)
            json = spm_jsonread(JSONlist{1});
            
            for iJSON=1:length(JSONfields)
                if isfield(json,JSONfields{iJSON})
                    fprintf(FID, xASL_num2str(json.(JSONfields{iJSON})));
                end
                fprintf(FID, ',');
            end
        end
        fprintf(FID, '\n'); % go to next line
    end
end
                    
fclose(FID);            