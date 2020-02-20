Rdir='/Users/henk/ExploreASL/ASL/SydneyMS_Controls/raw';
ControlList = xASL_adm_GetFileList(Rdir,'^.*$','List',[0 Inf],1);

Adir = '/Users/henk/ExploreASL/ASL/SydneyMS_Controls/analysis';
AnalysisList = xASL_adm_GetFileList(Adir,'^.*$','List',[0 Inf],1);

for iA=1:length(AnalysisList)
    Cohort{iA,1} = AnalysisList{iA};
    if sum(strcmp(AnalysisList{iA},ControlList))==0
        Cohort{iA,2} = 'case';
    else
        Cohort{iA,2} = 'control';
    end
end

Cohort = Cohort(1:147,:);

save(fullfile(Adir,'Cohort.mat'),'Cohort');