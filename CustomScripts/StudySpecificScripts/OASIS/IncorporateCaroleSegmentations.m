%% Incorporate Carole WMH segmentations

ExploreASL_Master('',0);

FLAIRdir = 'C:\BackupWork\ASL\OASIS\OASIS_WMH\FLAIR';
WMHdir = 'C:\BackupWork\ASL\OASIS\OASIS_WMH\Lesion';
xASLdir = 'C:\BackupWork\ASL\OASIS\OASIS3_xASL';

CleanupBeforeCompleteRerun(xASLdir);

FLAIRlist = xASL_adm_GetFileList(FLAIRdir,'^FLAIR_OAS.*\.nii$','List',[0 Inf]);
WMHlist = xASL_adm_GetFileList(WMHdir,'^Correct.*OAS.*\.nii$','List',[0 Inf]);
xASLlist = xASL_adm_GetFileList(xASLdir,'OAS\d*_\d*','List',[0 Inf], true);

fprintf('Moving Caroles WMH segmentations:   ');

for iWMH=1:length(WMHlist)
    xASL_TrackProgress(iWMH, length(WMHlist));
    [StartI, EndI] = regexp(WMHlist{iWMH}, 'OAS\d*_MR_d\d*_');
    SubjName = WMHlist{iWMH}(StartI:EndI-1);
    [StartI, EndI] = regexp(SubjName, '^OAS\d*_');
    Name1 = SubjName(StartI:EndI-1);
    [StartI, EndI] = regexp(SubjName, '_d\d*');
    Name2 = SubjName(StartI+2:EndI);
    SubjName2 = [Name1 '_' Name2]; % name with longitudinal ID
    
    % find the FLAIR with this index:
    IndexFLAIR = find(cellfun(@(x) ~isempty(regexp(x, SubjName)), FLAIRlist));
    % find the xASL folder
    IndexXASL = find(cellfun(@(x) ~isempty(regexp(x, SubjName2)), xASLlist));
    
    if ~isempty(IndexFLAIR) && ~isempty(IndexXASL)
        FLAIRpath = fullfile(FLAIRdir, FLAIRlist{IndexFLAIR(1)});
        WMHpath = fullfile(WMHdir, WMHlist{iWMH});
        xASLDir = fullfile(xASLdir, xASLlist{IndexXASL(1)});
        FLAIRnew = fullfile(xASLDir, 'FLAIR.nii.gz');
        WMHnew = fullfile(xASLDir, 'WMH_SEGM.nii.gz');
        if xASL_exist(FLAIRpath, 'file') && xASL_exist(WMHpath, 'file')
            xASL_delete(FLAIRnew);
            xASL_delete(WMHnew);
            xASL_Move(FLAIRpath, FLAIRnew);
            xASL_Move(WMHpath, WMHnew);
        elseif length(IndexFLAIR)>1
            FLAIRpath = fullfile(FLAIRdir, FLAIRlist{IndexFLAIR(end)});
            if xASL_exist(FLAIRpath, 'file') && xASL_exist(WMHpath, 'file')
                xASL_delete(FLAIRnew);
                xASL_delete(WMHnew);
                xASL_Move(FLAIRpath, FLAIRnew);
                xASL_Move(WMHpath, WMHnew);
            end
        end
    end
end