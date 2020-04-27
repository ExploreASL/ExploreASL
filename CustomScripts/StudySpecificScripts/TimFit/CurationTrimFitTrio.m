% Curation TimFit Trio

ExploreASL_Master('',0);

Root = '/Users/henk/ExploreASL/ASL/TimFit';
RawDir = fullfile(Root, 'raw');

%% Rename Subjects dirs with different visits
SubjectList = xASL_adm_GetFileList(RawDir,'SUB\d*','List',[0 Inf],1);
for iSubject=1:length(SubjectList)
    for iVisit=1:2
        OriDir = fullfile(RawDir, SubjectList{iSubject}, ['V' num2str(iVisit)]);
        NewDir = fullfile(RawDir, [SubjectList{iSubject} '_' num2str(iVisit)]);
        
        if exist(OriDir,'dir') && ~exist(NewDir,'dir')
            xASL_Move(OriDir, NewDir);
        end
    end
    SubjectDir = fullfile(RawDir, SubjectList{iSubject});
    if exist(SubjectDir, 'dir')
        if isempty(xASL_adm_GetFileList(SubjectDir, '.*','FPListRec'))
            xASL_delete(SubjectDir);
        end
    end
end


ExploreASL_Import(ExploreASL_ImportConfig(Root), 0, 0, 0, 1, 0)