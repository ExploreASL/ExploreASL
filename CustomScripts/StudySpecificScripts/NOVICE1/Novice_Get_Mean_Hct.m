%% Sleep study get mean hematocrit for each session/cohort condition

Cohort              = x.S.SetsID;

for iH=1:size(hematocrit,1)
    hematocritOLD(iH,1)     = hematocrit{iH,2};
end
    
% with SPSS, mean hct's (significantly different)
% HC (cohort1) have higher Hct, so we underestimate CBF for HC relative to
% HIV.
hct_cohort1     = 0.392941176470588;
hct_cohort2     = 0.3700;

HCT(1:size(x.S.SetsID,1),1) = 0;
HCT(x.S.SetsID==1)          = hct_cohort1;
HCT(x.S.SetsID==2)          = hct_cohort2;

save( fullfile(x.D.ROOT,'Mean_Hct_Session_Cohort.mat') , 'HCT');

unique(HCT)
