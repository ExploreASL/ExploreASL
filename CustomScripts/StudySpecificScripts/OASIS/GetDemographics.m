% Compute OASIS Stats
OASIS={}; % copy/paste data from excel CSV
ROOT = 'C:\BackupWork\ASL\OASIS';

MATpath = fullfile(ROOT, 'ASLoverview.mat');
save(MATpath, 'OASIS');

OASISnew = {};

% select subjects that have ASL data
for iS=2:size(OASIS,1) % data has headers
    if ~isempty(regexp(OASIS{iS,4},'asl'))
        OASISnew(end+1,1) = OASIS(iS,1);
        for iL=2:size(OASIS,2)
            OASISnew(end,iL) = OASIS(iS,iL);
        end
    end
end
        
% define TPs
SubjList = unique(OASISnew(:,2));
for iS=1:length(SubjList)
    SubjID = find(strcmp(OASISnew(:,2), SubjList{iS}));
    % check age order (later scans should have later age)
    for AgeN=2:length(SubjID)
        if floor(OASISnew{SubjID(AgeN),3})<floor(OASISnew{SubjID(AgeN-1),3})
            error([OASISnew{SubjID(AgeN),2} ' has wrong age order']);
        elseif OASISnew{SubjID(AgeN),3}==OASISnew{SubjID(AgeN-1),3} || OASISnew{SubjID(AgeN),3}+0.25==OASISnew{SubjID(AgeN-1),3}
            OASISnew{SubjID(AgeN),3} = OASISnew{SubjID(AgeN)-1,3}+0.5;
            % age is integer, therefore assume that follow-up was few
            % months at least
        end
    end
    
    for iI=1:length(SubjID)
        OASISnew{SubjID(iI),2} = [OASISnew{SubjID(iI),2} '_' num2str(iI)];
    end
end
   
OASIS=OASISnew;

Sex = {}; % copy/paste data from excel CSV
SexNew = Age(:,1);
for iS=1:length(SexNew)
    SexNew{iS,2} = 'n/a';
    Subj = SexNew{iS,1}(1:end-2);
    SubjID = find(strcmp(Sex(:,1), Subj));
    if ~isempty(SubjID)
        SexID = NaN;
        for iI=1:length(SubjID)
            if isempty(Sex{SubjID(iI),2}) || ~isfinite(Sex{SubjID(iI),2})
                SexID(iI+1) = NaN;
            else
                SexID(iI+1) = Sex{SubjID(iI),2};
            end
        end
        SexID = unique(SexID(isfinite(SexID)));
        if isempty(SexID)
            % skip this one
        elseif SexID(1)==1
            SexNew{iS,2} = 'M';
        elseif SexID(1)==2
            SexNew{iS,2} = 'F';
        else
            error('Unknown sex');
        end
    end
end
    
Sex = SexNew;

% stats age
for iS=1:length(Age)
    AgeN(iS) = Age{iS,2};
end

[X, N] = hist(AgeN);
figure(1);plot(N, X)
mean(AgeN)
std(AgeN)

% stats sex
for iS=1:length(Sex)
    if strcmp(Sex{iS,2}, 'M')
        SexN(iS) = 1;
    elseif strcmp(Sex{iS,2}, 'F')
        SexN(iS) = 2;
    else
        SexN(iS) = NaN;
    end
end

MPerc = sum(SexN==1)/numel(SexN)
FPerc = sum(SexN==2)/numel(SexN)
MPerc+(1-MPerc-FPerc)

% stats TP
for iS=1:length(Age)
    TP(iS,1) = str2num(Age{iS,1}(end));
end

numel(TP)

N5 = sum(TP==5);
N4 = sum(TP==4) - N5;
N3 = sum(TP==3) - N4;
N2 = sum(TP==2) - N3;
N1 = sum(TP==1) - N2;
FUnum = (N5*5+N4*4+N3*3+N2*2+N1*1)/numel(TP)+1;


%% Add new reconciliation key
ROOT = 'C:\BackupWork\ASL\OASIS';
MATpath = fullfile(ROOT, 'ASLoverview.mat');
AgePath = fullfile(ROOT, 'Age.mat');
SexPath = fullfile(ROOT, 'Sex.mat');
load(MATpath,'-mat');
load(SexPath,'-mat');

for iScan=1:size(OASIS,1)
    OASIS{iScan,5} = [OASIS{iScan,1}(1:8) '_' OASIS{iScan,1}(14:end)];
end
Age = OASIS(:,[5,3]);
save(AgePath,'Age');

% Same for sex
SexOld = Sex;
Sex = OASIS(:,[5]);

for iScan=1:size(OASIS,1)
    FirstPartName = Sex{iScan,1}(1:8);
    Indices = find(cellfun(@(y) ~isempty(regexp(y,FirstPartName)), SexOld(:,1)));
    Sex{iScan,2} = SexOld{Indices(1),2};
end

save(SexPath,'Sex');

