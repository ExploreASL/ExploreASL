function xASL_wrp_AssembleQC(x)
%xASL_wrp_AssembleQC Wrapper that 1) collects all QC information per site
% 2) creates database as CSV-file
% 3) Generates JPG-plots for each QC parameter

%% 1) Admin: load all QC JSONs

fprintf('Collecting and assembling QC parameters\n');

DirQC = fullfile(x.S.StatsDir, 'QC_stats');
xASL_adm_CreateDir(DirQC);

% a) make list of JSON files
FileList = xASL_adm_GetFileList(x.D.ROOT, '^QC_collection_.*\.json$', 'FPListRec');

TotalQC = struct;

fprintf('Loading data:   ');
% Load all data
for iFile=1:length(FileList)
    xASL_TrackProgress(iFile, length(FileList));
    
    try
        % First check JSON validity
        IsValid = xASL_qc_CheckValidityJSON(FileList{iFile});
        if IsValid
            [~, Ffile] = fileparts(FileList{iFile});
            SubjectName = ['QC' Ffile(15:end)];
            qc = spm_jsonread(FileList{iFile});
            
            Fields1 = fields(qc);
            for iField1=1:length(Fields1)
                Fields2 = fields(qc.(Fields1{iField1}));
                for iField2=1:length(Fields2)
                    TotalQC.(Fields1{iField1}).(Fields2{iField2}).(SubjectName) = qc.(Fields1{iField1}).(Fields2{iField2});
                end
            end
        end
    catch ME
        warning(['Something wrong for ' FileList{iFile}]);
        fprintf([ME.message '\n']);
    end
end
fprintf('\n');

%% 2) Create CSV file
PathCSV = fullfile(x.S.StatsDir, 'QC_Table.csv');
fclose all;
xASL_delete(PathCSV);
FID = fopen(PathCSV,'wt');

% write first row
fprintf(FID, 'Parameters,');
for iSubject=1:x.nSubjects
    fprintf(FID, [x.SUBJECTS{iSubject} ',']);
end
fprintf(FID, '\n');


Scantypes = fields(TotalQC);
for iScantype=1:length(Scantypes) % loop over scantypes
    Parameters = fields(TotalQC.(Scantypes{iScantype})); % loop over parameters
    for iParameter=1:length(Parameters)
        fprintf(FID, [Scantypes{iScantype} '_' Parameters{iParameter} ',']);
        
        for iSubject=1:x.nSubjects % loop over subjects
            if isfield(TotalQC.(Scantypes{iScantype}).(Parameters{iParameter}),['QC' x.SUBJECTS{iSubject}])
                fprintf(FID, [xASL_num2str( TotalQC.(Scantypes{iScantype}).(Parameters{iParameter}).(['QC' x.SUBJECTS{iSubject}]) ) ',']);
            else
                fprintf(FID, 'NaN,');
            end
        end
        fprintf(FID, '\n');
    end
end

fclose all;

%% 3) Print plots
% a) load it
DirQC = fullfile(x.S.StatsDir, 'QC_stats');
PathCSV = fullfile(x.S.StatsDir, 'QC_Table.csv');
[~, QC] = xASL_bids_csv2tsvReadWrite(PathCSV, false, false);

iSetSite = find(strcmp(x.S.SetsName,'Site'));
SiteIDs = unique(x.S.SetsID(:,iSetSite));
for iSiteID=1:length(SiteIDs) % loop over sites
    
    IndexSite = find(x.S.SetsID(:,iSetSite)==SiteIDs(iSiteID));
    SiteColumns = 1+IndexSite; % columns of Table having the subjects
    SiteName = x.S.SetsOptions{iSetSite}{iSiteID};
    fprintf(['Printing QC plots for Site ' SiteName ':   ']);
    
    for iParameter=2:size(QC,1) % loop over parameters
        xASL_TrackProgress(iParameter, size(QC,1));
        ParameterName = QC{iParameter, 1};
        
        if isnumeric(xASL_str2num(QC{iParameter, SiteColumns(1)})) && isempty(regexp(lower(ParameterName),'(matrix|voxelsize|id|version|readout|vendor|labelingtype)'))
            % here we skip the QC parameters that are not numerical,
            % e.g. the ID
            Values = cellfun(@(y) xASL_str2num(y), QC(iParameter, SiteColumns));

            % calculate stats
            MeanValue = xASL_stat_MeanNan(Values);
            SDValue = xASL_stat_StdNan(Values);

            close all;
            fig = figure('Visible','off');
            plot([1:length(Values)], Values,'k.');
            hold on
            plot([1,length(Values)], [MeanValue MeanValue],'g');
            plot([1,length(Values)], [MeanValue-SDValue MeanValue-SDValue],'--r');
            plot([1,length(Values)], [MeanValue+SDValue MeanValue+SDValue],'--r');
            plot([1,length(Values)], [MeanValue-1.96*SDValue MeanValue-1.96*SDValue],'-r');
            plot([1,length(Values)], [MeanValue+1.96*SDValue MeanValue+1.96*SDValue],'-r');
            xlabel('SubjectNumber (green=mean, red = SD 1 & 1.96)');
            ylabel(ParameterName);
            TitleName = [ParameterName '_Site_' SiteName];
            title(TitleName);
            SiteIDs = QC(1,SiteColumns);
            for iOutliers = 1:length(Values)
                if Values(iOutliers)> (MeanValue+(2*SDValue)) || Values(iOutliers)<(MeanValue-(2*SDValue))
                    text(iOutliers,Values(iOutliers), SiteIDs(iOutliers))
                end 
            end 

            PathPrint = fullfile(DirQC, [TitleName '.jpg']);
            saveas(fig, PathPrint, 'jpg');
        end
    end
    fprintf('\n');
end

