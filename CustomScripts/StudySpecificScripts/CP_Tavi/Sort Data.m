%% Sort Data Tavi

% 1) Rename _1 & _2

ROOT = 'C:\Backup\ASL\CP-Tavi\raw';
Dlist   = xASL_adm_GetFsList(ROOT,'^(AMC|HC)_\d{2}$',1);

for iD=1:length(Dlist)
    Dlist2  = xASL_adm_GetFsList(fullfile(ROOT,Dlist{iD}),'^(AMC|HC)_\d{2}_MR1$',1);
    for iD2=1:length(Dlist2)
        xASL_Rename(fullfile(ROOT,Dlist{iD}, Dlist2{iD2}),[Dlist2{iD2}(1:end-3) '1']);
    end
 
    Dlist2  = xASL_adm_GetFsList(fullfile(ROOT,Dlist{iD}),'^(AMC|HC)_\d{2}_MR2$',1);
    for iD2=1:length(Dlist2)
        xASL_Rename(fullfile(ROOT,Dlist{iD}, Dlist2{iD2}),[Dlist2{iD2}(1:end-3) '2']);
    end
end
