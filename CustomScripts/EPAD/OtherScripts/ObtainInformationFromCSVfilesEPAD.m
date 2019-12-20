%% Obtain CSV-data

ExploreASL_Master('',0);

rDir                = 'C:\Backup\ASL\EPAD\raw'; % root dir
SubjList            = xASL_adm_GetFsList(rDir,'^\d{3}-\d{5}$',1,[],[],[0 Inf]);
SiteN               = {'010'         '011'         '012'        '020'          '022'          '030'          '031'       '040'            '050'          '060'};
SiteName            = {'010SiPrisma' '011SiPrisma' '012SiVerio' '020PhAchieva' '022PhIngenia' '030PhIngenia' '031SiTrio' '040PhIngenuity' '050PhIngenia' '060Skyra'};


for iS=1:length(SubjList)
    xASL_TrackProgress(iS,length(SubjList));

    %%  Collect subject-wise metadata
    clear CSVfile rawData IndD DoB IndS SexInd

    CSVfile     = xASL_adm_GetFileList(fullfile(rDir,SubjList{iS}), '^.*\.csv$','FPListRec',[0 Inf]);
    [~, ~, rawData] = xlsread( CSVfile{1} );
    % Data is strangely stored, separated by '|', try to unwringle here
    rawData     = rawData{2,1};
    IndD        = strfind(rawData,'-07-01|');
    DoB         = rawData(IndD-4:IndD+5);
    IndS        = strfind(rawData,'|M|');
    if ~isempty(IndS)
        SexInd     = 'Male';
    else
        IndS        = strfind(rawData,'|F|');
        if ~isempty(IndS)
            SexInd     = 'Female';
        else
            error('Unknown sex');
        end
    end

    ScanDate            = rawData(IndS+3:IndS+12);
    AgeYears            = str2num(ScanDate(1:4))  - str2num(DoB(1:4));
    Months              = str2num(ScanDate(6:7))  - str2num(DoB(6:7));
    Days                = str2num(ScanDate(9:10)) - str2num(DoB(9:10));

    Age{iS,1}           = SubjList{iS};
    Sex{iS,1}           = SubjList{iS};
    Site{iS,1}          = SubjList{iS};

    Age{iS,2}           = AgeYears+Months/12+Days/365.25;
    Sex{iS,1}           = SubjList{iS};
    Sex{iS,2}           = SexInd;
    
    % Detect scanner (this is referred to as "site" in ExploreASL
    for iSite=1:length(SiteN)
        if strcmp(SiteN{iSite}, SubjList{iS}(1:3) )
            Site{iS,2}       = SiteName{iSite};
        end
    end
end

SaveList    = {'Age' 'Sex' 'Site'};
for iL=1:length(SaveList)
    SavePath    = fullfile(rDir,[SaveList{iL} '.mat']);
    save(SavePath,SaveList{iL});
end
