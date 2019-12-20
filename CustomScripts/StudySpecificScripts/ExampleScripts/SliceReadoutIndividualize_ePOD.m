ROOT            = 'D:\Backup\ASL_E\ePOD'; 

FPlist          = xASL_adm_GetFileList(ROOT, '^ASL4D_parms\.mat$','FPListRec')

ListAMC         = {'041_bl' '034_bl'};

for iList=1:length(FPlist)
    clear parms
    load(FPlist{iList})
    for iAMC=1:length(ListAMC)
        if  ~isempty(strfind(FPlist{iList},ListAMC{iAMC}))
            % This subject is from AMC
            parms.SliceReadoutTime  = 28.3;
        end
    end
    if ~isfield(parms,'SliceReadoutTime')
        parms.SliceReadoutTime  = 36.4;
    end
    
    save(FPlist{iList},'parms');
end
