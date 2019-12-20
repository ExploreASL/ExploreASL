%% Create Site Code

RootDir     = 'C:\Backup\ASL\EPAD\analysis';
Dlist       = xASL_adm_GetFsList(RootDir,'^\d{3}-\d{5}',1);

SiteN               = {'010'         '011'         '012'        '020'          '022'          '030'          '031'       '040'            '050'          '060'};
SiteName            = {'010SiPrisma' '011SiPrisma' '012SiVerio' '020PhAchieva' '022PhIngenia' '030PhIngenia' '031SiTrio' '040PhIngenuity' '050PhIngenia' '060Skyra'};
clear Site
for iD=1:length(Dlist)
    Site{iD,1}      = Dlist{iD};
    
    IndexN = find(cellfun(@(x) strcmp(x,Dlist{iD}(1:3)), SiteN));
    
    Site{iD,2}      = SiteName{IndexN};
end

save( fullfile(RootDir,'Site.mat'),'Site');
