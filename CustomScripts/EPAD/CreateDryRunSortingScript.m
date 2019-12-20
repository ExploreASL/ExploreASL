function CreateDryRunSortingScript
%CreateDryRunSortingScript Summary of this function goes here
%   Detailed explanation goes here

if isunix
    ROOT = '/mnt/s4e_data/RAD/share/EPAD500_new/raw';
    DEST = '/mnt/s4e_data/RAD/share/EPAD500_new/DryRun/raw/'
    ExploreASLdir = '/mnt/s4e_data/RAD/share/ExploreASL/ExploreASL/'
else
    ROOT = 'c:\Backup\ASL\EPAD\raw2';
    DEST = 'c:\Backup\ASL\EPAD\Test';
    ExploreASLdir = 'c:\ExploreASL';
end

cd /
cd(ExploreASLdir);
ExploreASL_Master('',0); % load scripts
cd(fullfile('CustomScripts', 'EPAD'));
Dlist = xASL_adm_GetFileList(ROOT, '.*', 'FPListRec', [0 Inf], true);

fprintf('%s', 'Copying files:   ')
for iL=1:length(Dlist)
    xASL_TrackProgress(iL, length(Dlist));
    dcmList = xASL_adm_GetFileList(Dlist{iL},'^.*[^.(csv|tsv|gz)]$','FPList',[0 Inf]);
    if ~isempty(dcmList)
        % select half
        OriFile = dcmList{ceil(0.5*length(dcmList))};
        DestFile = fullfile(DEST, OriFile(length(ROOT)+2:end));
        xASL_Copy(OriFile, DestFile, true, false);
    end
end
        
fprintf('\n');
% exit

end