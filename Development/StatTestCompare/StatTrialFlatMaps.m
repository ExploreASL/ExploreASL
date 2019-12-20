%% Administration
x.PrintDir    = fullfile( x.D.PopDir , 'StatsCompare');
xASL_adm_CreateDir(x.PrintDir);

%% Study 1 -> FLAT MAPS
% Create empty maps

FlatImage{1}      = zeros(121,145,121,15,'single');
FlatImage{2}      = zeros(121,145,121,15,'single');










%% Convert to single slice for time being
aslSlice{1}  = squeeze(FlatImage{1}(:,:,74,:));
aslSlice{2}  = squeeze(FlatImage{2}(:,:,74,:));

%% Create a BLOB signal
% Kernel = CreateBlob( FwHm, SliceDim,                                  NormValue, PeakValue, NoiseFactor, CutOffLevel, TranslVox )
% BLOB            = CreateBlob( 12, [size(aslSlice,1),size(aslSlice,2)], 1, xASL_stat_MeanNan(aslSlice(:))*4, 0, 0.01, [30 -10]);
BLOB            = CreateBlob( 16, [size(aslSlice{1},1),size(aslSlice{1},2)], 1, 100, 0, 0.01, [30 -10]); % 100=100%

% Plot2Dsurface(xASL_im_rotate(BLOB,90));

% Count number of voxels with signal
SigMask         = BLOB>0;

% Create random (normal) vector for variation of BLOB across scans, from -VarBLOBFactor : +VarBLOBFactor
VarBLOBFactor   = 0.33; % 25% variability between subjects
nScansTotal     = size(aslSlice{1},3) + size(aslSlice{2},3);
VarBLOB         = 1+(randn( nScansTotal )-0.5).*VarBLOBFactor;  

NoiseFactor     = 100;

for iScan=1:size(aslSlice{1},3)
    %% Add BLOB to images of group 1 only
    aslSliceBLOB{1}(:,:,iScan)      = aslSlice{1}(:,:,iScan) + BLOB.*VarBLOB(iScan);
    aslSliceBLOB{2}(:,:,iScan)      = aslSlice{2}(:,:,iScan) ;
    %% Add noise to image
    aslSliceBLOB{1}(:,:,iScan)      = aslSliceBLOB{1}(:,:,iScan) + ((randn( size(aslSlice{1},1), size(aslSlice{1},2) ) + 0.5).* NoiseFactor);
    aslSliceBLOB{2}(:,:,iScan)      = aslSliceBLOB{2}(:,:,iScan) + ((randn( size(aslSlice{1},1), size(aslSlice{1},2) ) + 0.5).* NoiseFactor);
    %% Print image
%     Plot2Dsurface(xASL_im_rotate(aslSliceBLOB{1}(:,:,iScan),90), fullfile(x.PrintDir, ['Flat1PlusBlob_' num2str(iScan) '.jpg']) );
%     Plot2Dsurface(xASL_im_rotate(aslSliceBLOB{2}(:,:,iScan),90), fullfile(x.PrintDir, ['Flat2PlusBlob_' num2str(iScan) '.jpg']) );
end

Plot2Dsurface(xASL_im_rotate(mean(aslSliceBLOB{1},3),90), fullfile(x.PrintDir, ['Flat1PlusBlob_Mean.jpg']) );
Plot2Dsurface(xASL_im_rotate(mean(aslSliceBLOB{2},3),90), fullfile(x.PrintDir, ['Flat2PlusBlob_Mean.jpg']) );

Plot2Dsurface( smooth( mean(aslSliceBLOB{1},3),8/2.2 ) , fullfile(x.PrintDir, ['Flat1PlusBlob_MeanSmoothed.jpg']) );
Plot2Dsurface( smooth( mean(aslSliceBLOB{2},3),8/2.2 ) , fullfile(x.PrintDir, ['Flat2PlusBlob_MeanSmoothed.jpg']) );


%% Ttest between groups
% Repeat test with different thresholds, and calculate sensitivity & specificity
% With blurring
%Group1          = dip_array( smooth( aslSliceBLOB{1}, 8/2.2 ) );
%Group2          = dip_array( smooth( aslSliceBLOB{2}, 8/2.2 ) );
Group1          = xASL_im_ndnanfilter(aslSliceBLOB{1},'gauss',[8 8 8],0);
Group2          = xASL_im_ndnanfilter(aslSliceBLOB{2},'gauss',[8 8 8],0);

% Without blurring
% Group1          = aslSliceBLOB{1};
% Group2          = aslSliceBLOB{2};

pMin            = 0; % minimal p-value (statistical threshold)
pMax            = 0.01;  % maximal p-value (statistical threshold)
ROCRes          = 100;   % resolution, or nSteps
StepSize        = (pMax - pMin) ./ ROCRes;

[H P]           = ttestExploreASL( Group1, Group2, 0.01, 'both', 3);
clear H StatCalc
StatCalc(:,1)   = [pMin:StepSize:pMax]';

% [X,Y]           = perfcurve(reshape(logical(BLOB),prod(size(P)),1) ,reshape(P,prod(size(P)),1),'numeric');

for iStat=1:size(StatCalc,1)

    H           = P<StatCalc(iStat,1);
    H(isnan(H)) = 0;

    GoldPos     = logical(BLOB);
    GoldNeg     = ~logical(BLOB);

    TruePos     = H .* GoldPos;
    FalsePos    = H .* GoldNeg;

    TPR         = sum(TruePos(:))  ./ sum(GoldPos(:));
    FPR         = sum(FalsePos(:)) ./ sum(GoldNeg(:));

    StatCalc(iStat,2)   = TPR;
    StatCalc(iStat,3)   = FPR;
    clear H TPR FPR TruePos FalsePos GoldPos GoldNeg
end

% Make last element [1,1], for area under curve (AUC) calculation
% CAVE: if last element will interpolate line from last point to [1,1], 
%       the above statistic threshold iteration was not adequate enough!
if  min(StatCalc(length(StatCalc),2:3))~=1
    StatCalc(length(StatCalc)+1,2:3)    = 1;
end

figure
plot(StatCalc(:,3),StatCalc(:,2)); % plot ROC
hold on
plot([0:0.01:1],[0:0.01:1],'k--'); % plot line of equality
xlabel('False positive rate (1-specificity)');
ylabel('True positive rate');
AUC     = trapz(StatCalc(:,3),StatCalc(:,2));
AUC     = round(AUC*100)/100; % round to 2 floating numbers
text(0.7,0.2,['AUC=' num2str(AUC)]);
