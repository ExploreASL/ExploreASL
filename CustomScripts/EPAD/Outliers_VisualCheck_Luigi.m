

function Outliers_VisualCheck(x, Parameter, Site)

% This Function can be used to select outliers  in a distribution of QC
% parameter. 


DirQC = fullfile(x.S.StatsDir, 'QC_stats');
PathCSV = fullfile(x.S.StatsDir, 'QC_Table.csv');
[~, QC] = xASL_bids_csvRead(PathCSV);
iSetSite = find(strcmp(x.S.SetsName,'Site'));
SiteIDs = unique(x.S.SetsID(:,iSetSite));
Analysisdir = '/data/RAD/share/EPAD500_new/analysis';
Populationdir = fullfile(Analysisdir, 'Population');
ParameterName = Parameter;
iParameter = QC(:,1)==string(Parameter);
iSiteID=find(SiteIDs == str2num(Site));
IndexSite = find(x.S.SetsID(:,iSetSite)==SiteIDs(iSiteID));
SiteColumns = 1+IndexSite; % columns of Table having the subjects
SiteName = x.S.SetsOptions{iSetSite}{iSiteID};
Site_QC = QC(:,SiteColumns);
LuigiDir = '/data/home/l.lorenzini/QC_EPAD';



Values = cellfun(@(y) xASL_str2num(y), QC(iParameter, SiteColumns));

if ~isempty(Values) 
% calculate stats
MeanValue = xASL_stat_MeanNan(Values);
SDValue = xASL_stat_StdNan(Values);

[MaxVal, idx_MaxVal] = maxk(Values, 2);
[~, idx_MinVal] = mink(Values, 2);
[~, idx_AvVal] = mink(abs(Values-MeanValue), 2);


if ~isnan(MaxVal(1))
Subjects = [{'High'}, Site_QC(1,idx_MaxVal(1)), Site_QC(1,idx_MaxVal(2));
    {'Average'}, Site_QC(1,idx_AvVal(1)),Site_QC(1,idx_AvVal(2));
    {'Low'}, Site_QC(1,idx_MinVal(1)), Site_QC(1,idx_MinVal(2))];

T1_files = [{'High'}, fullfile(Populationdir, strcat('rT1_', Site_QC(1,idx_MaxVal(1)), '.nii')), fullfile(Populationdir,strcat( 'rT1_', Site_QC(1,idx_MaxVal(2)), '.nii'));
    {'Average'}, fullfile(Populationdir,strcat('rT1_', Site_QC(1,idx_AvVal(1)), '.nii')),fullfile(Populationdir,strcat('rT1_', Site_QC(1,idx_AvVal(2)), '.nii'));
    {'Low'}, fullfile(Populationdir, strcat('rT1_',Site_QC(1,idx_MinVal(1)), '.nii')), fullfile(Populationdir, strcat('rT1_',Site_QC(1,idx_MinVal(2)), '.nii'))];

Flair_files = [{'High'}, fullfile(Populationdir, strcat( 'rFLAIR_' , Site_QC(1,idx_MaxVal(1)), '.nii')), fullfile(Populationdir, strcat( 'rFLAIR_', Site_QC(1,idx_MaxVal(2)), '.nii'));
    {'Average'}, fullfile(Populationdir,strcat( 'rFLAIR_',Site_QC(1,idx_AvVal(1)), '.nii')),fullfile(Populationdir,strcat( 'rFLAIR_',Site_QC(1,idx_AvVal(2)), '.nii'));
    {'Low'}, fullfile(Populationdir, strcat( 'rFLAIR_',Site_QC(1,idx_MinVal(1)), '.nii')), fullfile(Populationdir, strcat( 'rFLAIR_',Site_QC(1,idx_MinVal(2)), '.nii'))];

Func_files = [{'High'}, fullfile(Populationdir, strcat('SD_', Site_QC(1,idx_MaxVal(1)), 'func.nii')), fullfile(Populationdir,strcat('SD_',Site_QC(1,idx_MaxVal(2)), 'func.nii'));
    {'Average'}, fullfile(Populationdir,strcat('SD_',Site_QC(1,idx_AvVal(1)), 'func.nii')),fullfile(Populationdir,strcat('SD_',Site_QC(1,idx_AvVal(2)), 'func.nii'));
    {'Low'}, fullfile(Populationdir, strcat('SD_',Site_QC(1,idx_MinVal(1)), 'func.nii')), fullfile(Populationdir, strcat('SD_',Site_QC(1,idx_MinVal(2)),'func.nii'))];

%Func_SD_files = [{'High'}, fullfile(Populationdir, Site_QC(1,idx_MaxVal(1)), 'func',  'SD.nii'), fullfile(Populationdir,Site_QC(1,idx_MaxVal(2)), 'func','SD.nii');
 %   {'Average'}, fullfile(Populationdir,Site_QC(1,idx_AvVal(1)), 'func','SD.nii'),fullfile(Populationdir,Site_QC(1,idx_AvVal(2)), 'func', 'SD.nii');
  %  {'Low'}, fullfile(Populationdir, Site_QC(1,idx_MinVal(1)), 'func', 'SD.nii'), fullfile(Populationdir, Site_QC(1,idx_MinVal(2)),'func', 'SD.nii')];

DWI_files = [{'High'}, fullfile(Populationdir, strcat('rdwi_FA_', Site_QC(1,idx_MaxVal(1)), '.nii')), fullfile(Populationdir,strcat('rdwi_FA_',Site_QC(1,idx_MaxVal(2)), '.nii'));
    {'Average'}, fullfile(Populationdir,strcat('rdwi_FA_',Site_QC(1,idx_AvVal(1)), '.nii')),fullfile(Populationdir,strcat('rdwi_FA_',Site_QC(1,idx_AvVal(2)), '.nii'));
    {'Low'}, fullfile(Populationdir, strcat('rdwi_FA_',Site_QC(1,idx_MinVal(1)), '.nii')), fullfile(Populationdir, strcat('rdwi_FA_',Site_QC(1,idx_MinVal(2)),'.nii'))];

%ASL_files = [{'High'}, fullfile(Populationdir, Site_QC(1,idx_MaxVal(1)), 'ASL_1', 'CBF.nii'), fullfile(Populationdir,Site_QC(1,idx_MaxVal(2)), 'ASL_1', 'CBF.nii');
  %  {'Average'}, fullfile(Populationdir,Site_QC(1,idx_AvVal(1)), 'ASL_1', 'CBF.nii'),fullfile(Populationdir,Site_QC(1,idx_AvVal(2)), 'ASL_1', 'CBF.nii');
  %  {'Low'}, fullfile(Populationdir, Site_QC(1,idx_MinVal(1)), 'ASL_1', 'CBF.nii'), fullfile(Populationdir, Site_QC(1,idx_MinVal(2)),'ASL_1', 'CBF.nii')];


Parameterfolder = fullfile(LuigiDir, 'Outliers_inspection', Parameter, SiteName);
figurefolder = fullfile(LuigiDir, 'Outliers_inspection', 'Figures', Parameter, SiteName);
if ~exist(figurefolder)
mkdir(figurefolder)
end
x.S.TraSlices=[40 63];
x.S.SagSlices=[35 55];
x.S.CorSlices=[54 74];


for iValue = 1:3
    for iS = 2:3
        if startsWith(Parameter, 'Structural') && exist(T1_files{iValue,iS}, 'file')
            T1im = xASL_im_CreateVisualFig(x, T1_files{iValue,iS}); Flairim = xASL_im_CreateVisualFig(x, Flair_files{iValue,iS});
            T1im = fliplr(T1im) ; Flairim = fliplr(Flairim);
            imwrite(T1im, fullfile (figurefolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_T1', '.jpg')));
            imwrite(Flairim, fullfile (figurefolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_Flair', '.jpg')));
            %copyfile(T1_files{iValue,iS}, fullfile (Parameterfolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_T1', '.nii')))
            %copyfile(Flair_files{iValue,iS}, fullfile (Parameterfolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_Flair', '.nii')))
        elseif startsWith(Parameter, 'func') && exist(Func_files{iValue,iS}, 'file')
            boldIm = xASL_im_CreateVisualFig(x, Func_files{iValue,iS}) ; %SDIm = xASL_im_CreateVisualFig(x, Func_SD_files{iValue,iS});
            boldIm = fliplr(boldIm) ;% SDIm = fliplr(SDIm);
            imwrite(boldIm, fullfile (figurefolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_BOLD', '.jpg')))
            %imwrite(SDIm, fullfile (figurefolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_SD', '.jpg')))

            %copyfile(Func_files{iValue,iS}, fullfile (Parameterfolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_BOLD', '.nii')))
            %copyfile(Func_SD_files{iValue,iS}, fullfile (Parameterfolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_SD', '.nii')))
        elseif startsWith(Parameter, 'dwi') && exist(DWI_files{iValue,iS}, 'file')
            dwiIm = xASL_im_CreateVisualFig(x, DWI_files{iValue,iS}) ;
            dwiIm = fliplr(dwiIm);
            imwrite(dwiIm, fullfile (figurefolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_FA', '.jpg')));
            %copyfile(DWI_files{iValue,iS}, fullfile (Parameterfolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_FA', '.nii')))
%        elseif startsWith(Parameter, 'ASL') && exist(ASL_files{iValue,iS}, 'file')
           % aslIm = xASL_im_CreateVisualFig(x, ASL_files{iValue,iS});
            %aslIm = fliplr(aslIm);
           % imwrite(aslIm, fullfile (figurefolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_CBF', '.jpg')));
            %copyfile(ASL_files{iValue,iS}, fullfile (Parameterfolder, strcat(Subjects{iValue}, '_', Subjects{iValue,iS},'_CBF', '.nii')))
            
        end
    end
end

graph_file = [ParameterName '_Site_' SiteName '.jpg'];
graph_path = fullfile(DirQC, graph_file);
P = imread(graph_path);
image(P);    
copyfile(graph_path, Parameterfolder)
end
end
end
