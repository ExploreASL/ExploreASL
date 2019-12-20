%% APGEM_2

OriFolder   = 'C:\Backup\ASL\APGEM2\APGEM_2_folders_first';
OriNames    = xASL_adm_GetFsList(OriFolder,'^\d{3}B_A\d{5}$',1);
RawFolder   = 'C:\Backup\ASL\APGEM2\raw';

NewFolder   = 'C:\Backup\ASL\APGEM2\analysis';
LockFolder  = 'C:\Backup\ASL\APGEM2\analysis\lock\Struct';

for iO=1:length(OriNames)
    clear CurrName NewName OldDir NewDir FileList
    CurrName    = [OriNames{iO}(6:end) '_2'];
    NewName     = [OriNames{iO}(1:4) '_2'];
    
    % Directories
    OldDir      = fullfile(NewFolder, CurrName);
    if  isdir(OldDir)
        xASL_Rename(OldDir,NewName);
    end
    
    clear OldDir NewDir
    % LockStruct-dirs
    OldDir      = fullfile(LockFolder, CurrName);
    if  isdir(OldDir)
        xASL_Rename(OldDir,NewName);
    end    
    
    clear OldDir NewDir
    % raw-folder-dirs
    OldDir      = fullfile(RawFolder, CurrName);
    if  isdir(OldDir)
        xASL_Rename(OldDir,NewName);
    end        
    
    % Files
    FileList    = xASL_adm_GetFileList(NewFolder,['.*' CurrName '.*'],'FPListRec',[0 Inf]);
    
    for iL=1:length(FileList)
        clear Oldname Fpath Ffile Fext IndexN NewFile
        OldName     = FileList{iL};
        [Fpath Ffile Fext]  = fileparts(OldName);
        
        IndexN      = findstr(Ffile,CurrName);
        NewFile     = [Ffile(1:IndexN-1) NewName Ffile(IndexN+length(CurrName):end) Fext];
        xASL_Rename(OldName,NewFile);
    end
end

%% APGEM_2 into APGEM code
APGEM_ID{1,1}   = 'A80048_2';
APGEM_ID{1,2}   = 'A80109_2';

clear
NewFolder   = 'C:\Backup\ASL\APGEM\analysis';
LockFolder  = 'C:\Backup\ASL\APGEM\analysis\lock\ASL';
OriNames    = xASL_adm_GetFsList(NewFolder,'^A80048_2$',1);
load(fullfile('C:\Backup\ASL\APGEM','APGEM_ID.mat'));
RawFolder   = 'C:\Backup\ASL\APGEM\raw';
DARTELfolder= 'C:\Backup\ASL\APGEM\analysis';

for iS=1:length(OriNames)
    for iL=1:length(APGEM_ID)
        if  strcmp(OriNames{iS},APGEM_ID{iL,1})
            % Rename SubjectFolder
            clear OldDir NewName
            NewName     = APGEM_ID{iL,2};
            
            OldDir      = fullfile(NewFolder,OriNames{iS});
            if  exist(OldDir,'dir')
                xASL_Rename(OldDir,NewName);
            end

            % Rename LockFolder
            OldDir      = fullfile(LockFolder,OriNames{iS});
            if  exist(OldDir,'dir')
                xASL_Rename(OldDir,NewName);
            end        
            
            % Rename LockFolder
            OldDir      = fullfile(RawFolder,OriNames{iS});
            if  exist(OldDir,'dir')
                xASL_Rename(OldDir,NewName);
            end           
            
            % Rename files
            FileList    = xASL_adm_GetFileList(DARTELfolder,['.*' OriNames{iS} '.*'],'FPListRec',[0 Inf]);

            for iL=1:length(FileList)
                clear Oldname Fpath Ffile Fext IndexN NewFile
                OldName     = FileList{iL};
                [Fpath Ffile Fext]  = fileparts(OldName);

                IndexN      = findstr(Ffile,OriNames{iS});
                NewFile     = [Ffile(1:IndexN-1) NewName Ffile(IndexN+length(OriNames{iS}):end) Fext];
                xASL_Rename(OldName,NewFile);
            end            
            
        end
    end
end


%% APGEM_1 into APGEM code
clear
NewFolder   = 'C:\Backup\ASL\APGEM\analysis';
LockFolder  = 'C:\Backup\ASL\APGEM\analysis\lock\Struct';
OriNames    = xASL_adm_GetFsList(NewFolder,'^\d{3}B$',1);
load(fullfile('C:\Backup\ASL\APGEM','APGEM_ID.mat'));
RawFolder   = 'C:\Backup\ASL\APGEM\raw';
DARTELfolder= 'C:\Backup\ASL\APGEM\analysis';

for iS=1:length(OriNames)
    for iL=1:length(APGEM_ID)
        if  strcmp(OriNames{iS}(1:4),APGEM_ID{iL,1}(1:4))

            % Rename SubjectFolder
            clear OldDir NewName
            NewName     = APGEM_ID{iL,2};
            
            OldDir      = fullfile(NewFolder,OriNames{iS});
            if  exist(OldDir,'dir')
                xASL_Rename(OldDir,NewName);
            end

            % Rename LockFolder
            OldDir      = fullfile(LockFolder,OriNames{iS});
            if  exist(OldDir,'dir')            
                xASL_Rename(OldDir,NewName);            
            end
            
            % Rename RawFolder
            OldDir      = fullfile(RawFolder,[OriNames{iS} '_1']);
            if  exist(OldDir,'dir')            
                xASL_Rename(OldDir,NewName);            
            end              
            
            % Rename files
            FileList    = xASL_adm_GetFileList(DARTELfolder,['.*' OriNames{iS} '.*'],'FPListRec',[0 Inf]);

            for iL=1:length(FileList)
                clear Oldname Fpath Ffile Fext IndexN NewFile
                OldName     = FileList{iL};
                [Fpath Ffile Fext]  = fileparts(OldName);

                IndexN      = findstr(Ffile,OriNames{iS});
                NewFile     = [Ffile(1:IndexN-1) NewName Ffile(IndexN+length(OriNames{iS}):end) Fext];
                xASL_Rename(OldName,NewFile);
            end            
            
        end
    end
end


%% APGEM_rest raw into APGEM code
clear
NewFolder   = 'C:\Backup\ASL\APGEM\analysis';
LockFolder  = 'C:\Backup\ASL\APGEM\analysis\lock\Struct';
load(fullfile('C:\Backup\ASL\APGEM','APGEM_ID.mat'));
RawFolder   = 'C:\Backup\ASL\APGEM\raw';
OriNames    = xASL_adm_GetFsList(RawFolder,'^\d{3}B$',1);
DARTELfolder= 'C:\Backup\ASL\APGEM\analysis';

for iS=1:length(OriNames)
    for iL=1:length(APGEM_ID)
        if  strcmp(OriNames{iS}(1:4),APGEM_ID{iL,1}(1:4))
            
            
            % Rename SubjectFolder
            clear OldDir NewName
            NewName     = APGEM_ID{iL,2};   
            
            % Rename RawFolder
            OldDir      = fullfile(RawFolder,[OriNames{iS}]);
            xASL_Rename(OldDir,NewName);
            
        end
    end
end


%% Rename variable
load('C:\Backup\ASL\APGEM\analysis\Sex.mat');
load('C:\Backup\ASL\APGEM\analysis\age.mat');
load('C:\Backup\ASL\APGEM\analysis\cohort.mat');

clear iS iI iL

for iS=1:length(cohort)
    clear Indices

    for iI=1:length(APGEM_ID)
        if  strcmp(cohort{iS,1}(1:4),APGEM_ID{iI,1}(1:4))
            cohort{iS,1}   = APGEM_ID{iI,2};
        end
    end
end

save('C:\Backup\ASL\APGEM\analysis\Sex.mat','Sex');









