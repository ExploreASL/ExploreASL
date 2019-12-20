clear cohort age sex
iS                              = 0;



%% AgeIV Cobra
iS                              = iS+1;

sex{iS}                         = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;1;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;0;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;1;0;1;0;0;1;0;0;1;0;1;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;1]';
% cohort{iS}  % for controls & diseases                    = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;-1;-1;-1;-1;-1;-1;-1;-1;-1;0;-1;-1;-1;0;-1;-1;0;-1;-1;-1;-1;-1;0;-1;0;-1;-1;0;-1;-1;-1;-1;-1;0;0;0;0;0;-1;-1;0;-1;-1;-1;-1;0;0;0;0;-1;-1;-1;0;0;-1;0;-1;-1;0;-1;-1;0;-1;-1;-1;0;-1;-1;-1;-1;0;-1;0;-1;0;0;0;-1;-1;-1;-1;-1;0;-1;-1;0;-1;-1]';
age{iS}                         = 57 + 7.*randn(size(sex{iS}));
% MultipleScansFactor(iS)         = 1.8708; % for follow-up

%% BrainAir
iS                              = iS+1;

cohort{iS}(1:12)                = 0; cohort{iS}(13:23)   = -1;
age{iS}                         = 43.5 + 7.*randn(size(cohort{iS}));
sex{iS}(1:length(cohort{iS}))   = 0; sex{iS}([5,10,15,20])  = 1;
MultipleScansFactor(iS)         = 1;

%% COCO
iS                              = iS+1;

age{iS}                         = [40.9397260300000;40.1452054800000;26.7452054800000;57.3643835600000;43.4739726000000;28.4054794500000;39.7424657500000;34.8219178100000;50.4849315100000;26.6904109600000;47.4136986300000;34.8739726000000;47.0438356200000;59.1780821900000;31.9123287700000;22.2383561600000;37.4821917800000;31.9424657500000;43.9287671200000;44.5342465800000;21.9178082200000;67.1506849300000;26.4931506800000;40.7287671200000;35.7698630100000;28.0054794500000;62.7671232900000;69.7917808200000;50.6712328800000;46.1479452100000;49.5287671200000;38.7232876700000;41.0602739700000;28.0794520500000;25.4876712300000;34.3726027400000;36.6821917800000;55.2904109600000;27.6438356200000;20.6575342500000;23.7424657500000;22.1287671200000;49.7726027400000;54.2082191800000;50.7013698600000;25.9232876700000;24.0794520500000;45.9315068500000]';
sex{iS}                         = [1,1,0,1,0,0,0,1,1,0,1,1,0,1,0,1,1,0,1,1,1,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,1,1,0,1,1,0,0,1,0,0,1,1];
cohort{iS}(1:length(age{iS}))   = -1;
MultipleScansFactor(iS)         = 1;

%% CRUISE patients
iS                              = iS+1;

cohort{iS}(1:40)                = 1;
age{iS}                         = 33.3 + 12.3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(40,13/(13+22) );
MultipleScansFactor(iS)         = 3;

%% CRUISE controls
iS                              = iS+1;

cohort{iS}(1:15)                = 0;
age{iS}                         = 33.3 + 14.3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(15,7/15 );
MultipleScansFactor(iS)         = 3;

%% FIND
iS                              = iS+1;

age{iS}                         = [10.2000000000000;13.5000000000000;19.1000000000000;10.8000000000000;11.6000000000000;14.8000000000000;17.5000000000000;17;13.7000000000000;10.2000000000000;16.7000000000000;17.8000000000000;14.4000000000000;13.3000000000000;10.6000000000000;16.6000000000000;17.2000000000000;14.7000000000000;11.2000000000000;11.6000000000000;14.9000000000000;19.4000000000000;13.9000000000000;10.3000000000000;18.2000000000000;18.8000000000000;15;16.7000000000000;16.8000000000000;15.2000000000000;12.8000000000000;17.3000000000000;10;11;13.8000000000000;14;20.1000000000000;20.5000000000000;9.80000000000000;11.8000000000000]';
sex{iS}                         = [1;1;1;1;2;1;1;1;2;2;2;1;1;1;2;2;2;1;2;1;2;1;1;1;1;1;1;2;1;1;1;2;2;2;2;2;2;1;2;1]';
cohort{iS}(1:length(age{iS}))   = 1;
MultipleScansFactor(iS)         = 1;

%% Hongerwinter
iS                              = iS+1;

age{iS}                         = [68.1500000000000;67.4100000000000;67.5800000000000;67.6500000000000;68.9900000000000;67.1200000000000;66.8900000000000;67.3300000000000;67.1300000000000;67.3300000000000;66.1000000000000;68.5100000000000;66.7100000000000;67.1000000000000;68.2300000000000;68.5500000000000;67.3900000000000;69.0200000000000;67.7200000000000;67.1800000000000;66.9900000000000;68.4700000000000;68.6600000000000;68.5200000000000;67.0900000000000;66.3900000000000;66.9700000000000;68.4800000000000;69.1300000000000;67.2800000000000;66.8500000000000;67.5100000000000;68.4000000000000;67.5500000000000;67.9900000000000;66.6400000000000;68.4600000000000;67.2600000000000;67.6700000000000;66.1800000000000;68.1500000000000;67.2100000000000;67.1600000000000;67.2600000000000;66.8600000000000;67.3700000000000;68.6600000000000;68.3800000000000;68.2100000000000;66.1500000000000;69.0200000000000;67.4500000000000;67.8800000000000;66.8100000000000;67.1800000000000;67.5100000000000;67.4200000000000;66;66.4700000000000;66.1600000000000;68.9200000000000;66.3200000000000;69.0300000000000;66.2900000000000;68.6200000000000;66.4200000000000;69.2900000000000;68.5400000000000;69.5900000000000;67.3400000000000;67.3500000000000;67.2200000000000;67.1600000000000;68.1600000000000;69.2200000000000;66.6400000000000;67.9100000000000;67.3800000000000;67.4200000000000;67.1600000000000;68.0400000000000;68.5300000000000;67.3400000000000;68.5000000000000;66.9400000000000;67.3400000000000;66.9400000000000;68.2500000000000;67.3000000000000;67.4300000000000;67.3100000000000;67.1800000000000;67.4300000000000;66.5100000000000;68.6400000000000;69.2200000000000;67.3900000000000;68.6700000000000;65.8900000000000;68.7300000000000;68.2900000000000;69.4700000000000;66.8600000000000;65.9000000000000;69.0200000000000;66.7500000000000;67.2700000000000;66.5200000000000;66.4100000000000;67.4300000000000;67.3100000000000;67.5300000000000;67.2300000000000;67.4500000000000;67.5200000000000;68.9800000000000;66.5600000000000;66.5300000000000;67.1200000000000;65.9600000000000;66.8100000000000;67.3200000000000;69.4000000000000;69.2900000000000;69.5000000000000;66.7600000000000;67.5800000000000;68.7600000000000;67.6500000000000;67.7200000000000;67.0300000000000;67.4900000000000;68.4800000000000;67.5400000000000;67.2500000000000]';
sex{iS}                         = [1;1;1;1;0;1;1;0;1;0;0;0;1;1;0;1;1;1;1;0;1;0;0;1;1;1;1;1;0;1;1;0;1;1;1;1;0;0;1;0;1;0;0;0;0;1;1;1;1;0;0;1;1;0;1;1;0;0;1;1;0;1;1;0;0;1;0;0;1;1;0;0;0;1;0;1;1;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;0;1;1;1;1;1;0;1;0;0;1;0;1;1;0;1;0;1;1;1;0;1;1;1;1;1;1;1;1;0;1;1;0;0;1;1;1;0;0]';
cohort{iS}                      = [0;-1;-1;-1;0;-1;0;-1;0;-1;0;0;0;0;0;0;-1;0;-1;-1;0;0;0;0;0;0;-1;0;0;-1;0;-1;0;-1;0;0;0;-1;-1;0;0;-1;0;0;0;0;0;0;0;0;0;-1;-1;0;0;-1;-1;0;0;0;0;0;0;0;0;0;0;0;0;-1;-1;0;0;0;0;0;0;-1;-1;-1;0;0;-1;0;0;-1;0;0;0;-1;-1;0;-1;0;0;0;-1;0;0;0;0;0;0;0;0;0;-1;0;0;-1;-1;-1;-1;-1;-1;0;0;0;0;0;0;-1;0;0;0;0;-1;0;-1;-1;0;-1;0;-1;0]';
MultipleScansFactor(iS)         = 1;

%% Insomnia
iS                              = iS+1;

age{iS}                         = [47;60;66;58;47;69;44;48;55;43;60;53;54;47;66;63;23;23;67;55;60;65;26;47;52;61;60;48;41;22;22;68;67;26;29;58;36;66;36;55;58;51;62;47.1214;45;45;61;56;43;59;51;30;35;43;42;44;27;70;61;65;59;62;62;67;63;47;70;27;39;68;54;42;67;30;23;27;62;30;44;59;61;69;53;37;25;43;59;67;42;28;36;58;55;27;62;35;34;66;68;27;60;34;59;67;27;41;56;53;20;53;32;34;40;32;37;23;63;56;20;39;25;54;30;48;29;39;29;21;62;63;29;24;27;32;51;55;24;53;31;25;37;31;46;30;28;44;62;36;49;57;54;54;56;64;60;63;56;68;26;63;40;40;45;57;62;68;29;52;53;33;59]';
age{iS}                         = age{iS}+(randn(size(age{iS}))./6);
sex{iS}                         = [1;1;0;1;1;1;1;0;1;0;1;1;1;1;0;1;1;0;1;1;0;1;1;1;1;1;1;1;1;1;0;1;1;0;1;1;1;1;0;1;1;1;0;1;1;1;0;1;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;0;1;0;1;1;1;1;0;0;0;0;0;0;1;0;0;1;1;1;1;1;1;0;0;0;1;1;1;0;0;1;0;1;1;1;1;1;1;1;0;1;0;1;1;1;1;1;0;0;1;1;1;0;0;0;0;0;1;0;1;0;1;0;1;0;0;1;0;1;0;1;1;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]';
cohort{iS}                      = [1;1;0;1;0;1;1;0;1;1;1;1;1;1;0;1;1;0;1;1;1;1;1;1;1;1;1;1;1;0;0;1;0;0;1;0;0;1;1;1;1;0;1;0;0;0;1;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;0;0;0;1;1;1;1;0;0;0;1;0;1;0;1;0;0;1;0;1;1;1;0;0;0;0;1;0;1;0;1;0;1;1;0;0;0;1;1;1;1;1;1;1;0;0;0;1;1;1;1;0;0;1;0;0;0;0;0;0;1;1;0;0;1;0;0;1;1;1;1;1;0;0;0;0;0;0;1;0;0;0;0;0;1;1;1;1;1;1;0;0;0;1;1;1;0;0;0;0;1;1;0;1;1;1]';
MultipleScansFactor(iS)         = 1;

%% Novice
iS                              = iS+1;

cohort{iS}(1:35)                = 0; cohort{iS}(36:71)  = 1;
age{iS}                         = 13 + 1.5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),(16+18)/(16+18+19+18) );
MultipleScansFactor(iS)         = 2;

%% Parelsnoer
iS                              = iS+1;

cohort{iS}(1:90)                = -1;
age{iS}                         = 75 + 10.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),48/(48+42) );
MultipleScansFactor(iS)         = 1;

%% PET-ASL
iS                              = iS+1;

cohort{iS}(1:15)                = 0;
age{iS}                         = 22 + 1.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),7/(7+8) );
MultipleScansFactor(iS)         = 12;

%% PreDiva
iS                              = iS+1;

cohort{iS}(1:187)               = -1;
age{iS}                         = 78 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),101/(101+86) );
MultipleScansFactor(iS)         = 1.6631;

%% Score
iS                              = iS+1;

cohort{iS}(1:28)                = 1;
age{iS}                         = 53 + 13.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.63 );
MultipleScansFactor(iS)         = 1.75;

%% Sickle cell Gevers
iS                              = iS+1;

cohort{iS}(1:15)                = 1; cohort{iS}(16:21)                = 0;
age{iS}                         = 16 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),(8+3)/(8+3+7+3) );
MultipleScansFactor(iS)         = 1.5;

%% Sickle cell Tweel
iS                              = iS+1;

cohort{iS}(1:24)                = 1; cohort{iS}(25:36)                = 0;
age{iS}                         = 14 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 1;

%% Sickle cell 7T
iS                              = iS+1;

cohort{iS}(1:10)                = 1;
age{iS}                         = 23 + 2.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 1;

%% Binnewijzend, VUmc
iS                              = iS+1;

cohort{iS}(1:150)               = -1;
age{iS}                         = 75 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 1;

%% T1DM reproducibility PET-ASL, VUmc
iS                              = iS+1;

cohort{iS}(1:11)                = 0; cohort{iS}(12:31)                = -1;
age{iS}                         = 35 + 11.3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0 );
MultipleScansFactor(iS)         = 1;

%% VESPA
iS                              = iS+1;

cohort{iS}(1:22)                = 0;
age{iS}                         = 23 + 2.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),13/(13+9) );
MultipleScansFactor(iS)         = 12;

%% Twins
iS                              = iS+1;

cohort{iS}(1:175)               = 0;
age{iS}                         = 77 + 10.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 1;

%% AMPH
iS                              = iS+1;

cohort{iS}(1:20)                = 0; cohort{iS}(21:40)                = -1;
age{iS}                         = 21 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0 );
MultipleScansFactor(iS)         = 1;

%% EPOD ADHD 1
iS                              = iS+1;

cohort{iS}(1:50)                = -1;
age{iS}                         = 11 + 1.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0 );
MultipleScansFactor(iS)         = 5;

%% EPOD ADHD 2
iS                              = iS+1;

cohort{iS}(1:49)                = 1;
age{iS}                         = 28 + 2.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0 );
MultipleScansFactor(iS)         = 5;

%% EPOD pharmo ADHD
iS                              = iS+1;

cohort{iS}(1:25)                = -1;
age{iS}                         = 27 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0 );
MultipleScansFactor(iS)         = 2;

%% EPOD pharmo MDD/ADHD
iS                              = iS+1;

cohort{iS}(1:50)                = -1;
age{iS}                         = 30 + 8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),1 );
MultipleScansFactor(iS)         = 1;

%% EPOD sert HC
iS                              = iS+1;

cohort{iS}(1:45)                = 0;
age{iS}                         = 22 + 6.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),1 );
MultipleScansFactor(iS)         = 1;

%% EPOD neuroshape HC
iS                              = iS+1;

cohort{iS}(1:52)                = 0;
age{iS}                         = 25 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),30/(30+22) );
MultipleScansFactor(iS)         = 2;

% %% Incorporate number of scans
% for iS=1:length(cohort)
%     cohort{iS}                  = repmat(cohort{iS},[1 round(MultipleScansFactor(iS))]);
%     age{iS}                     = repmat(age{iS},[1 round(MultipleScansFactor(iS))]);
%     sex{iS}                     = repmat(sex{iS},[1 round(MultipleScansFactor(iS))]);
% end    

%% Simulate CBF
close all
figure(1);

for iS=1:length(age)
    CBF{iS}         = 60-((age{iS}-40).*0.0075.*60); % age ~ CBF
    CBF{iS}         = CBF{iS} + sex{iS}.*CBF{iS}.*0.13; % add 13% perfusion for females
    CBF{iS}         = CBF{iS} + (cohort{iS}.*CBF{iS}.*0.2); % 20% CBF change for pathology
    CBF{iS}         = CBF{iS} + randn(size(CBF{iS})) .* 0.15 .* CBF{iS}; % add noise. randn on average adds 0, with SD = 1, which we multiply with 15% of bsCV CBF

    %    +randn(size(age)).*10+sex.*5;

    age1{iS}        = age{iS}(sex{iS}==1 & abs(cohort{iS})==1); % female pathology
    age2{iS}        = age{iS}(sex{iS}==0 & abs(cohort{iS})==1); %   male pathology
    age3{iS}        = age{iS}(sex{iS}==1 & abs(cohort{iS})==0); % female HC
    age4{iS}        = age{iS}(sex{iS}==0 & abs(cohort{iS})==0); %   male HC

    CBF1{iS}        = CBF{iS}(sex{iS}==1 & abs(cohort{iS})==1); % female pathology
    CBF2{iS}        = CBF{iS}(sex{iS}==0 & abs(cohort{iS})==1); %   male pathology
    CBF3{iS}        = CBF{iS}(sex{iS}==1 & abs(cohort{iS})==0); % female HC
    CBF4{iS}        = CBF{iS}(sex{iS}==0 & abs(cohort{iS})==0); %   male HC

    plot(age1{iS},CBF1{iS},'ro',age2{iS},CBF2{iS},'rx',age3{iS},CBF3{iS},'go',age4{iS},CBF4{iS},'gx');
    hold on
end

xlabel('Age (years)');
ylabel('CBF (mL/100g/min');
title('Overview ASL subjects Amsterdam Neuroscience (551 controls, 627 patients)');
axis([0 100 0 140])



%% BioCog
iS                              = iS+1;

cohort{iS}(1:339)               = 0;
age{iS}                         = 72 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),140/(140+199) );
MultipleScansFactor(iS)         = 1.7139;

%% 22q11
iS                              = iS+1;

cohort{iS}(1:27)                = -1;
age{iS}                         = 42 + 8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),11/(11+16) );
MultipleScansFactor(iS)         = 1;

%% 3CV
iS                              = iS+1;

cohort{iS}(1:11)                = 0;
age{iS}                         = 26 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 6;

%% ADHD HC Oslo
iS                              = iS+1;

cohort{iS}(1:6)                = 0;
age{iS}                         = 11 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 1;

%% ADHD patients Oslo
iS                              = iS+1;

cohort{iS}(1:50)                = 1;
age{iS}                         = 11 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 1;

%% ADNI go/2 AD
iS                              = iS+1;

cohort{iS}(1:32)                = -1;
age{iS}                         = 76 + 7.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),12/(12+21) );
MultipleScansFactor(iS)         = 1.5;

%% ADNI go/2 MCI
iS                              = iS+1;

cohort{iS}(1:152)                = -1;
age{iS}                         = 72 + 8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),68/(68+84) );
MultipleScansFactor(iS)         = 1.5;

%% ADNI go/2 NC
iS                              = iS+1;

cohort{iS}(1:66)                = 0;
age{iS}                         = 74 + 7.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),31/(31+35) );
MultipleScansFactor(iS)         = 1.5;

%% Antipsych
iS                              = iS+1;

cohort{iS}(1:53)                = 0;
age{iS}                         = 25 + 6.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),26/(26+27) );
MultipleScansFactor(iS)         = 2;

%% APGEM Czech controls
iS                              = iS+1;

cohort{iS}(1:56)                = 0;
age{iS}                         = 67 + 7.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),39/(39+17) );
MultipleScansFactor(iS)         = 1;

%% APGEM Czech AD
iS                              = iS+1;

cohort{iS}(1:14)                = -1;
age{iS}                         = 74 + 9.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),12/(12+2) );
MultipleScansFactor(iS)         = 1;

%% APGEM Czech MCI
iS                              = iS+1;

cohort{iS}(1:28)                = -1;
age{iS}                         = 70 + 7.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),17/(17+11) );
MultipleScansFactor(iS)         = 1;

%% APGEM Czech PD
iS                              = iS+1;

cohort{iS}(1:18)                = -1;
age{iS}                         = 61 + 8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),4/(4+14) );
MultipleScansFactor(iS)         = 1;

%% APGEM Czech PD/MCI
iS                              = iS+1;

cohort{iS}(1:25)                = -1;
age{iS}                         = 65 + 10.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),9/(9+16) );
MultipleScansFactor(iS)         = 1;

%% Bipolar
iS                              = iS+1;

cohort{iS}(1:45)                = 1;
age{iS}                         = 17 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),26/(26+19) );
MultipleScansFactor(iS)         = 1.64835164835165;

%% Bipolar HC
iS                              = iS+1;

cohort{iS}(1:46)                = 0;
age{iS}                         = 17 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),27/(27+19) );
MultipleScansFactor(iS)         = 1.64835164835165;

%% AAS Anabolic Steroid Users Oslo
iS                              = iS+1;

cohort{iS}(1:67)                = -1;
age{iS}                         = 33.3 + 8.8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),1 );
MultipleScansFactor(iS)         = 1;

%% AAS Anabolic Steroid non-Users Oslo
iS                              = iS+1;

cohort{iS}(1:46)                = 0;
age{iS}                         = 33.3 + 8.8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),1 );
MultipleScansFactor(iS)         = 1;

%% Dave UCL Young onset AD
iS                              = iS+1;

cohort{iS}(1:44)                = -1;
age{iS}                         = 62.3 + 4.9.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),19/(19+18));
MultipleScansFactor(iS)         = 1.25;

%% Dave UCL Young onset AD, HCs
iS                              = iS+1;

cohort{iS}(1:23)                = 0;
age{iS}                         = 60.7 + 6.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),13/(13+10) );
MultipleScansFactor(iS)         = 1.25;

%% Dave UCL Insight46, HCs
iS                              = iS+1;

cohort{iS}(1:450)               = 0;
age{iS}                         = 62 + 2.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5 );
MultipleScansFactor(iS)         = 1.5;

%% Hardy
iS                              = iS+1;

age{iS}                         = [58;63;87;72;64;73;52;79;77;77;51;76;81;62;88;62;64;60;68;86;62;87;61;52;76;68;82;88;68;69;60;77;82;68;78;68;79;72;52;59;67;81;66;85;63;80;56;70;75;61;61;59;65;85;65;69;65;67;77;68;76;63;61;62;58;68;65;62;69;78;81;82;77;82;58;65;77;60;59;53;74;73;78;74;73;64;72;77;76;55;63;73;75;82;72;71;62;54;80;81;75;85;70;77;77;68;82;63;77;62;66;73;52;71;74;75;77;87;74;72;71;78;57;58;79;77;66;71;70;82;74;65;65;81;72;73;77;68;77;74;80;80;86;78;95;80;71;81;67;72;77;91;65;79;64;66;77;74;76;89;84;63;80;67;52;75;65;71;67;69;71;65;63;65;73;69;65;62;71;73;67;64;86;82;78;63;68;84;83;53;74;74;71;65;85;65;85;75;66;74;68;85;84;52;78;78;50;83;55;84;79;71;78;75;69;73;63;86;75;68;80;62;65;69;84;73;62;84;71;81;73;84;73;81;65;67;76;78;81;57;75;88;83;77;84;77;66;60;78;78;65;78;73;57;83;61;68;66;81;71;77;76;78;83;67;69;83;69;78;74;70;71;74;78;74;72;74;72;62;74;61;75;72;78;79;66;75;72;75;84;72;75;76;88;76;77;87;76;58;79;81;76;77;66;79;83;79;73;74;78;68;61;62;63;74;77;65;85;54;73;62;92;60;81;72;65;70;85;80;73;74;72;70;67;72;76;83;87;77;84;68;76;74;70;83;69;64;69;65;69;67;73;67;75;66;78;82;79;67;81;77;73;70;69;73;75;70;69;79;73;67;78;70;63;72;74;71;73;69;76;75;70;61;69;73;74;71;83;73;82;74;74;69;69;73;73;66;71;74;71;79;72;74;78;85;62;72;72;85;63;85;67;65;65;79;70;71;77;78;76;70;77;80;69;70;71;62;71;60;75;82;75;75;81;77;74;55;74;83;75;75;83;62;76;74;71;82;68;77;59;82;74;77;77;81;78;82;65;78;74;83;77;73;73;76;71;70;74;63;77;80;77;73;75;76;70;79;76;76;68;78;62;73;78;91;80;75;70;75;80;77;77;79;65;78;73;80;76;75;70;75;77;77;82;70;77;81;71;68;69;76;75;68;80;84;77;80;78;74;84;70;79;80;60;89;63;67;80;69;80;87;76;73;90;85;73;78;74;75;68;73;79;73;80;66;75;76;79;74;78;70;75;76;82;87]';
age{iS}                         = age{iS}+(randn(size(age{iS}))./6);
sex{iS}                         = [0;1;1;0;0;1;0;0;0;0;1;1;0;1;1;1;1;0;1;1;1;1;1;0;0;1;0;1;1;0;0;1;1;0;0;1;1;0;0;0;1;1;1;1;0;0;1;0;1;1;0;1;1;0;0;1;1;1;0;0;1;1;1;1;1;0;1;1;1;1;0;1;0;1;1;1;0;0;0;0;1;1;1;1;1;0;0;1;0;0;0;0;0;1;1;1;0;0;1;0;0;1;1;0;1;1;1;1;0;0;0;0;0;1;0;1;0;0;1;0;0;1;1;0;1;0;0;0;0;0;1;1;0;0;1;1;1;0;0;1;0;0;0;0;0;1;0;0;0;1;1;1;0;1;0;0;0;0;1;0;1;0;1;0;0;0;0;0;0;0;0;0;0;0;1;0;0;1;0;0;0;0;1;0;1;0;0;1;0;0;1;1;1;1;1;0;0;1;1;0;1;1;0;1;0;0;1;1;0;1;1;1;0;0;1;1;1;1;1;1;1;1;0;1;1;0;1;1;0;1;1;0;0;1;1;1;1;0;1;0;0;0;1;1;1;0;1;1;0;1;1;1;1;1;0;1;0;1;0;0;0;1;0;1;1;0;0;1;1;1;0;1;0;0;1;1;1;0;1;1;1;0;1;0;0;0;0;1;1;0;1;1;0;1;0;1;1;1;0;0;1;1;1;0;1;0;1;1;1;0;1;0;1;1;1;1;0;1;0;1;1;1;0;0;0;1;0;1;1;1;0;1;1;1;0;1;1;1;1;1;0;1;0;0;0;1;1;0;1;1;1;0;1;0;1;0;1;1;1;0;0;1;1;1;1;1;1;1;0;1;0;0;1;1;0;0;0;0;1;0;1;1;0;1;0;1;0;0;1;0;1;0;1;1;0;0;1;1;0;1;0;0;1;0;1;0;1;0;0;0;0;0;1;0;1;1;1;0;1;0;1;1;0;1;1;1;1;1;0;1;0;1;1;1;1;1;1;1;1;1;0;1;0;1;0;0;1;0;1;1;0;1;0;1;0;1;1;1;0;0;1;1;0;1;1;1;1;1;1;1;1;0;0;0;1;1;1;0;0;0;0;1;0;1;1;0;1;0;1;0;1;1;1;0;0;0;0;1;1;1;1;0;0;0;1;1;0;1;1;0;1;1;1;1;0;1;0;0;1;1;1;0;1;0;0;1;1;0;0;1;0;1;0;0;1;1;0;1;1;0;1;1;1;0;0;0;1;0;1;1;1;1;0;0;0]';
cohort{iS}(1:length(age{iS}))   = -1; % all SVD cases
MultipleScansFactor(iS)         = 2;

%% Hardy population study (EDIS)
iS                              = iS+1;

cohort{iS}(1:861)               = 0;
age{iS}                         = 70.3 + 6.6.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.535);
MultipleScansFactor(iS)         = 1.5;

%% Epilepsia Czeck Republic
iS                              = iS+1;

cohort{iS}(1:31)                = 1;
age{iS}                         = 32.1 + 11.9.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),17/31);
MultipleScansFactor(iS)         = 1;

%% Freedivers
iS                              = iS+1;

cohort{iS}(1:15)                = 0;
age{iS}                         = 40 + 12.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0);
MultipleScansFactor(iS)         = 5;

%% GENFI carriers
iS                              = iS+1;

cohort{iS}(1:149)               = -1;
age{iS}                         = 46 + 12.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),108/(108+41));
MultipleScansFactor(iS)         = 1.6336;

%% GENFI non-carriers
iS                              = iS+1;

cohort{iS}(1:113)               = 0;
age{iS}                         = 50 + 14.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),72/(72+41));
MultipleScansFactor(iS)         = 1.6336;

%% KCL antipsych
iS                              = iS+1;

cohort{iS}(1:30)                = 0;
age{iS}                         = 50 + 14.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),72/(72+41));
MultipleScansFactor(iS)         = 2;

%% Mood Patricia
iS                              = iS+1;

cohort{iS}(1:14)                = 0;
age{iS}                         = 20 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),8/(8+6));
MultipleScansFactor(iS)         = 20;

%% Sleep study
iS                              = iS+1;

cohort{iS}(1:40)                = 0;
age{iS}                         = 22 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0);
MultipleScansFactor(iS)         = 3;

%% Sleep study
iS                              = iS+1;

cohort{iS}(1:40)                = 0;
age{iS}                         = 22 + 3.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0);
MultipleScansFactor(iS)         = 3;

%% Uppsala Dementia DF AD/FTD
iS                              = iS+1;

cohort{iS}(1:23)                = -1;
age{iS}                         = 75 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5);
MultipleScansFactor(iS)         = 1;

%% Uppsala Dementia DF HC
iS                              = iS+1;

cohort{iS}(1:40)                = -1;
age{iS}                         = 74 + 5.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5);
MultipleScansFactor(iS)         = 1;

%% Uppsala iNPH JV
iS                              = iS+1;

cohort{iS}(1:21)                = -1;
age{iS}                         = 74 + 8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5);
MultipleScansFactor(iS)         = 6;

%% Uppsala iNPH JV HC
iS                              = iS+1;

cohort{iS}(1:21)                = 0;
age{iS}                         = 74 + 8.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5);
MultipleScansFactor(iS)         = 1;

%% WM study Oslo
iS                              = iS+1;

cohort{iS}(1:8)                = 0;
age{iS}                         = 33 + 4.*randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0);
MultipleScansFactor(iS)         = 12;

%% Atle Dementia amyloid-
iS                              = iS+1;

cohort{iS}(1:51)                = 0;
age{iS}                         = 59.85 + 7.29 *randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),0.5);
MultipleScansFactor(iS)         = 1;

%% Atle Dementia amyloid+
iS                              = iS+1;

cohort{iS}(1:79)                = 0;
age{iS}                         = 63.77 + 6.87 *randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),19/(19+12));
MultipleScansFactor(iS)         = 1;

%% Brad Young HCs
iS                              = iS+1;

cohort{iS}(1:15)                = 0;
age{iS}                         = 25.7 + 4.6 *randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),7/15);
MultipleScansFactor(iS)         = 1;

%% Brad old HCs
iS                              = iS+1;

cohort{iS}(1:22)                = 0;
age{iS}                         = 69.1 + 7 *randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),17/22);
MultipleScansFactor(iS)         = 1;

%% Brad WMH HCs
iS                              = iS+1;

cohort{iS}(1:15)                = 0;
age{iS}                         = 72 + 8.2 *randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),8/(8+7));
MultipleScansFactor(iS)         = 2;

%% Brad WMH HCs
iS                              = iS+1;

cohort{iS}(1:30)                = 0;
age{iS}                         = 67.9 + 10.4 *randn(size(cohort{iS}));
sex{iS}                         = RandomizeSex(length(cohort{iS}),9/30);
MultipleScansFactor(iS)         = 2;

%% Simulate CBF
close all
figure(1);

for iS=1:length(age)
    CBF{iS}         = 60-((age{iS}-40).*0.0075.*60); % age ~ CBF
    CBF{iS}         = CBF{iS} + sex{iS}.*CBF{iS}.*0.13; % add 13% perfusion for females
    CBF{iS}         = CBF{iS} + (cohort{iS}.*CBF{iS}.*0.2); % 20% CBF change for pathology
    CBF{iS}         = CBF{iS} + randn(size(CBF{iS})) .* 0.15 .* CBF{iS}; % add noise. randn on average adds 0, with SD = 1, which we multiply with 15% of bsCV CBF

    %    +randn(size(age)).*10+sex.*5;

    age1{iS}        = age{iS}(sex{iS}==1 & abs(cohort{iS})==1); % female pathology
    age2{iS}        = age{iS}(sex{iS}==0 & abs(cohort{iS})==1); %   male pathology
    age3{iS}        = age{iS}(sex{iS}==1 & abs(cohort{iS})==0); % female HC
    age4{iS}        = age{iS}(sex{iS}==0 & abs(cohort{iS})==0); %   male HC

    CBF1{iS}        = CBF{iS}(sex{iS}==1 & abs(cohort{iS})==1); % female pathology
    CBF2{iS}        = CBF{iS}(sex{iS}==0 & abs(cohort{iS})==1); %   male pathology
    CBF3{iS}        = CBF{iS}(sex{iS}==1 & abs(cohort{iS})==0); % female HC
    CBF4{iS}        = CBF{iS}(sex{iS}==0 & abs(cohort{iS})==0); %   male HC

    plot(age1{iS},CBF1{iS},'ro',age2{iS},CBF2{iS},'rx',age3{iS},CBF3{iS},'go',age4{iS},CBF4{iS},'gx');
    hold on
end

xlabel('Age (years)');
ylabel('CBF (mL/100g/min)');
title('Overview ASL subjects');
axis([0 100 0 140])