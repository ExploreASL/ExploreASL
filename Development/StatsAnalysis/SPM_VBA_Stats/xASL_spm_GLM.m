% function xASL_spm_GLM( x )
%xASL_spm_GLM xASL_spm_GLM but rewritten for one-way
% ANalysis Of VAriance (ANOVA) & t-tests
% x = from ExploreASL
%
% By HJMM Mutsaerts, ExploreASL 2016

% NB: with the use of an explicit (external) mask, SPM should not
% crash as much as it did using an implicit mask only
% (the infamous error: no inmask voxels)
% Therefore, some "try-catch" constructions may be removed

MainStatsDir    = cd;
LoadList        = {'xASL_session Session_3 vs. Session_2 within cohort 2' 'xASL_session Session_3 vs. Session_1 within cohort 2'};
for iX=1:length(LoadList)
    clearvars -except MainStatsDir LoadList iX
    load(fullfile(MainStatsDir,[LoadList{iX} '.mat']));
    fclose all;
    close all;
    Dir2Del     = fullfile(MainStatsDir,'SPMdir');
    Flist       = xASL_adm_GetFileList(Dir2Del,'.*','FPListRec',[0 Inf]);
    for iG=1:length(Flist)
        delete(Flist{iG});
    end




%% Administration
x.S.nSets               = length(x.S.DATASETS);
for ii=1:x.S.nSets
    x.S.nScans{ii}       = size(x.S.DATASETS{ii},1);
end


%% Options
% x.S.nSets=1, then do 1-sample t-test
% x.S.nSets=2, then do 2-sample t-test, either paired or unpaired
% x.S.nSets>2, do ANOVA, either paired or unpaired

% Permute different co-variates to account for them,
% & do regressions


%% Save dataset as NII files
x.S.SPMdir          = fullfile( x.S.StatsDir, 'SPMdir');
TemplateNii         = fullfile( x.D.MapsDir, 'rgrey.nii');
xASL_adm_CreateDir(x.S.SPMdir);
ReCreateNII         = 0;

for iSet=1:x.S.nSets
    clear ExpectedNIIcount FoundNIIcount
    ExpectedNIIcount    = size(x.S.DATASETS{iSet},1);
    FoundNIIcount       = length(xASL_adm_GetFileList(x.S.SPMdir, ['^Set' num2str(iSet) 'subject\d*\.nii$'],'FPList',[0 Inf]));
    if  ExpectedNIIcount~=FoundNIIcount
        ReCreateNII     = 1;
    end
end

if  ReCreateNII
    xASL_adm_DeleteFileList(x.S.SPMdir,'^.*\.nii$');
    fprintf('%s\n','Creating images...  ');
    for iSet=1:x.S.nSets
        for iSubject=1:size(x.S.DATASETS{iSet},1)
            xASL_TrackProgress((iSet-1).*size(x.S.DATASETS{iSet},1)+iSubject,x.S.nSets.*size(x.S.DATASETS{iSet},1));
            PathNew     = fullfile( x.S.SPMdir, ['Set' num2str(iSet) 'subject' num2str( iSubject, '%05.0f') '.nii']);

            if  length(size(x.S.DATASETS{iSet}))==2
                tIM         = xASL_im_Column2IM(x.S.DATASETS{iSet}(iSubject,:),x.S.VBAmask);
            else
                tIM         = squeeze(x.S.DATASETS{iSet}(iSubject,:,:,:));
            end

            xASL_io_SaveNifti( TemplateNii, PathNew, tIM, 16,0);
        end
    end
end
fprintf('\n')



%% Put breakpoint here for custom SPM  ----  CustomMixedEffectsSPMmodel.m










    %% 1. Start running SPM
    x.S.SPMmat          = fullfile(x.S.SPMdir, 'SPM.mat');
    x.S.SPMmat1         = fullfile(x.S.SPMdir, 'SPM1.mat');
    x.S.SPMmat2         = fullfile(x.S.SPMdir, 'SPM2.mat');
    x.S.SPMmat3         = fullfile(x.S.SPMdir, 'SPM3.mat');
    x.S.SPMmat4         = fullfile(x.S.SPMdir, 'SPM4.mat');

    spm('defaults','PET');
    % CopyFileList

    cd(x.S.SPMdir);

    %% Defaults
    if ~isfield(x.S,'PrintSPMOutput')
        x.S.PrintSPMOutput  = 1;
    end
    if ~isfield(x.S,'GlobalNormalization')
        x.S.GlobalNormalization  = 0;
    end
    if ~isfield(x.S,'MultiComparisonCorrType')
        x.S.MultiComparisonCorrType     = 'cluster';
    end
    if ~isfield(x.S,'clusterPthr')
        x.S.clusterPthr     = 0.001;
    end
    if ~isfield(x.S,'uncorrThresh')
        x.S.uncorrThresh     = 0.05;
    end
    if ~isfield(x.S,'RegressionCOVAR')
        x.S.RegressionCOVAR     = 0;
    end


    x.S.MaskPath        = fullfile(x.D.PopDir,'GMSliceMask.nii');




    %% -------------------------------------------------------------------------------------
    %% 1) GLMDesign

    IsTP1_TP2   = ~isempty(strfind(lower(LoadList{iX}),'session_1')) && ~isempty(strfind(lower(LoadList{iX}),'session_2'));
    IsTP2_TP3   = ~isempty(strfind(lower(LoadList{iX}),'session_2')) && ~isempty(strfind(lower(LoadList{iX}),'session_3'));
    IsTP1_TP3   = ~isempty(strfind(lower(LoadList{iX}),'session_1')) && ~isempty(strfind(lower(LoadList{iX}),'session_3'));

%     %% Fix Timing Diff Covariate
%     for iSet=1:length(x.S.SetsName) % First find the set
%         if  strcmp(x.S.SetsName{iSet},'AcquisitionTime')
%             for ii=1:length(x.S.CoVar)
%                 % 1) convert to continuous days
%                 TempSet             = x.S.CoVar{ii}(:,iSet); % Convert to hours
%                 TempSet             = xASL_adm_ConvertTime2Nr(TempSet./100);
%                 if  ii==length(x.S.CoVar) && (IsTP2_TP3 || IsTP1_TP3)
%                     TempSet         = TempSet+24; % Add 24 hours for the second dayConvert to day 2
%                 end
%
%                 % 2) Impute NaNs
%                 MeanTemp            = xASL_stat_MeanNan(TempSet);
%                 Indices             = find(isnan( TempSet ));
%                 TempSet(Indices)    = MeanTemp;
%
% %                 % 3) Demean per time point
% %                 TempSet             = TempSet-MeanTemp;
%
%                 x.S.CoVar{ii}(:,iSet)   = TempSet;
%             end
%         end
%     end

    clear CoVariate
    %% Mention here which covariates you want to use, these will be loaded in the GLMdesign
    % E.g. for sleep study

    IsWithinCohort          = ~isempty(strfind(LoadList{iX},'within cohort'));

    if  IsTP1_TP2 || IsWithinCohort
        Sets2Load           = ''; % ''AcquisitionTime'' 'DiffAcqTime'   'MeanMotion' % Motion removes too much
        CoVarInteract       = [1]; % 3= interact with t-test (time-point), 1=interact with prefix factor (==1)
        IfDemean            = [1];
    else
        Sets2Load           = {'cohort'}; % 'AcquisitionTime'
        CoVarInteract       = [3]; % 3= interact with t-test (time-point), 1=interact with prefix factor (==1)
        IfDemean            = [0];
    end


    NextN           = 1;
    for iSet=1:length(Sets2Load)
        for iS=1:length(x.S.SetsName)
            if  strcmp(x.S.SetsName{iS},Sets2Load{iSet})
                CoVariate(NextN).Name = x.S.SetsName{iS};
                SingleLength    = length(x.S.CoVar{1}); % length of covariates for 1 group (e.g. if 2 groups are compared)
                FullLength      = length(x.S.CoVar)*SingleLength; % length of all covariates
                % Fill data
                for iC=1:length(x.S.CoVar)
                    if  IfDemean
                        CoVariate(NextN).Data(iC:length(x.S.CoVar):FullLength,1)	 = x.S.CoVar{iC}(:,iS)-mean(x.S.CoVar{iC}(:,iS)); % assuming t-test
                    else
                        CoVariate(NextN).Data(iC:length(x.S.CoVar):FullLength,1)	 = x.S.CoVar{iC}(:,iS); % assuming t-test
                    end
                end
                CoVariate(NextN).iCFI  = CoVarInteract(NextN);
                figure(NextN);plot(1:SingleLength,CoVariate(NextN).Data(1:2:end),'b',1:SingleLength,CoVariate(NextN).Data(2:2:end),'r');
                title(CoVariate(NextN).Name);
                xlabel('Scans');
                NextN   = NextN+1;
            end
        end
    end

    if  exist('CoVariate','var')
        [x] = xASL_spm_GLMdesign(x, CoVariate); % possibility covariates 2nd entrance
    else
        [x] = xASL_spm_GLMdesign(x);
    end


    %% -------------------------------------------------------------------------------------
    %% 2 Estimate GLMmodel
    [x] = xASL_spm_GLMmodel(x);


     Conweights      = [-1 1 0 0]; % main effect
%     Conweights      = [0 0 -1 1]; % interaction effect












    %% From here, thresholds are being used
    %  We try two primary thresholds, p=0.01 (exploratory, previously FSL default)
    %  p=0.001 (true threshold, SPM default)
    %  Actual primary p-threshold should be in the middle, which would be when using TCFE or other permutation
    %  -based methods
    x.S.printTitleORI   = x.S.NAME;
    printTitleORI   = x.S.printTitleORI;
    ClusterPthr     = [0.001]; % 0.01
    for iCon=1:length( ClusterPthr)
        x.S.clusterPthr     = ClusterPthr(iCon);


        %% -------------------------------------------------------------------------------------
        %% 3    Contrast creation
        if  exist('CoVariate','var') && exist('Conweights','var')
            [x] = xASL_spm_GLMcontrast(x, CoVariate, Conweights); % possibility covariates 2nd entrance
        else
            [x] = xASL_spm_GLMcontrast(x);
        end


        %% -------------------------------------------------------------------------------------
        %% 4    Create significance masks
        x   = xASL_spm_GLMresults(x);





        %% ------------------------------------------------------------------------------
        %% Print overview clusters

        %% Load contrast map for analysis
        % Using SPM thresholds
    %     MINthresh       = 3.53; % p=0.001 with GM masking

        if  x.S.nSets>2
            ContrastPath            = fullfile( x.S.SPMdir, 'spmF_0001.nii');
            MaskPath                = fullfile( x.S.SPMdir, 'spmF_Masked.nii');
        else
            ContrastPath            = fullfile( x.S.SPMdir, 'spmT_0001.nii');
            MaskPath                = fullfile( x.S.SPMdir, 'spmT_0001_Masked.nii');
        end

        x.S.ContrastMap             = xASL_io_Nifti2Im( ContrastPath );
        x.S.MaskMap                 = xASL_io_Nifti2Im( MaskPath );

        diff_view_mean{1}           = x.S.ContrastMap;
        x.S.printTitleORI           = [printTitleORI 'P' num2str(x.S.clusterPthr)];

        if  x.S.GlobalNormalization
            x.S.printTitleORI      = [x.S.printTitleORI '_ScaleGlblMean'];
        end

        %% Print overview clusters (MaskMap)
        % Print overview cluster names according to Harvard-Oxford
        xASL_stat_PrintCluster_ROI_Stats(x.S.MaskMap,x,x.S.printTitleORI);
        % Take clusters as "atlas ROIs", & produce CBF
        xASL_wrp_GetClusterCBFstats(x,x.S.printTitleORI,'qCBF');


        %% Figure used in Sleep Paper
        x.S.CorSlices   = [];
        x.S.SagSlices   = [];
        x.S.TraSlices   = 34:4:34+15*4; % 34:4:34+17*4 was original, now removed 2 slices


        %% Take result & put it in nice figures
        [diff_view_mean H_ttestCONTRAST x]     = xASL_spm_GLMcreateStatsFigures(x);
        xASL_spm_GLMsaveMatlab_StatsFigures(diff_view_mean, H_ttestCONTRAST,x);


        %% Move all files to specific Dir
        %% Part to save contrasts for ConjunctionAnalysis
        NewContrastDir      = [x.S.MultiComparisonCorrType 'P' num2str(x.S.clusterPthr) 'P' num2str(x.S.uncorrThresh) '_' x.S.NAME];
        NewContrastPath     = ['spmT' NewContrastDir '.nii'];
        StatsDir2           = fullfile(x.S.StatsDir,NewContrastDir);
        xASL_adm_CreateDir(StatsDir2);

        if  xASL_exist(MaskPath,'file')
            xASL_Copy( MaskPath, fullfile(StatsDir2, NewContrastPath));
        end
        %% Move PDF & SPM.mat files
        PDFlist             = xASL_adm_GetFileList(x.S.SPMdir,'^spm.*\.pdf$');
        for iP=1:length(PDFlist)
            PDFpathNew      = fullfile(StatsDir2,[printTitleORI '_' num2str(iP) '.pdf']);
            xASL_Move(PDFlist{iP},PDFpathNew);
        end
        MatPathNew          = fullfile(StatsDir2,[printTitleORI '_SPM.mat']);
        if  exist(x.S.SPMmat,'file')
            xASL_Move(x.S.SPMmat,MatPathNew);
        end

        % Move residual files
        Flist               = xASL_adm_GetFileList(x.S.StatsDir,'.*(\.nii|\.pdf|\.eps|SPM\.mat|\.tiff|\.tsv)','FPList',[0 Inf]);
        for iL=1:length(Flist)
            [~, Ffile Fext] = xASL_fileparts(Flist{iL});
            NewPath         = fullfile(StatsDir2,[Ffile Fext]);
            if      length(StatsDir2)>180
                    error('Too long Path: make DirName shorter!');
            elseif length(NewPath)>200
                    NewPath     = [NewPath(1:184) NewPath(length(NewPath)-15:end)];
            end
            xASL_Move( Flist{iL}, NewPath );
        end
    end


end





%
% % end
