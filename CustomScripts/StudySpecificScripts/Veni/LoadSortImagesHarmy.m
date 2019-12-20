%% Load & sort images

%% Put T1w in standard space without normalization
T1IM    = zeros(121,145,121,x.nSubjects,'single');

for iS=1:x.nSubjects
    xASL_TrackProgress(iS,x.nSubjects);
    NativeFile          = fullfile(x.D.ROOT,x.SUBJECTS{iS},'T1.nii');
    LoadFileGZ          = [NativeFile '.gz'];

    StandardFile        = fullfile(x.D.PopDir,['rT1_' x.SUBJECTS{iS} '_ORI.nii']);




    if  (exist(NativeFile,'file') || exist(LoadFileGZ,'file')) && ~exist(StandardFile,'file')
        xASL_io_ReadNifti(NativeFile);
        xASL_spm_reslice( x.D.ResliceRef, NativeFile, [], [], x.Quality, 'native to standard space T1 ', StandardFile, 4 );
    end

    if  exist(StandardFile,'file')
        T1IM(:,:,:,iS)  = xASL_io_Nifti2Im(StandardFile);
    end
end

%% Sort them for GM Volume
for iS=1:length(x.S.SetsName) % get which set is the GM-ICV ratio
    if  strcmp(x.S.SetsName{iS},'GM_ICVRatio')
        iSet    = iS;
    end
end

ColumnSort          = x.S.SetsID(:,iSet);
ColumnSort(:,2)     = [1:1:length(ColumnSort)];
ColumnSort          = sortrows(ColumnSort,1);
Mask                = ~isnan(ColumnSort);
ColumnSort          = ColumnSort(Mask(:,1),:);

GMimMax             = T1IM(:,:,:,ColumnSort(end-10:end,2));
GMimMin             = T1IM(:,:,:,ColumnSort(1:10,2));

dip_image(GMimMax)

NiceIMhighest       = GMimMax(:,:,:,9); % 'HD023_1'
NiceIMlowest        = GMimMin(:,:,:,7); % 'HD048_1'
dip_image(NiceIMlowest)

% %% Put ASL in standard space without normalization
% ASLIM    = zeros(121,145,121,x.nSubjects,'single');
%
% for iS=1:x.nSubjects
%     xASL_TrackProgress(iS,x.nSubjects);
%     NativeFile          = fullfile(x.D.ROOT,x.SUBJECTS{iS},'ASL_1','PWI.nii');
%     LoadFileGZ          = [NativeFile '.gz'];
%
%     StandardFile        = fullfile(x.D.PopDir,['qCBF_untreated_' x.SUBJECTS{iS} '_ORI.nii']);
%
%
%     if  (exist(NativeFile,'file') || exist(LoadFileGZ,'file')) && ~exist(StandardFile,'file')
%         xASL_io_ReadNifti(NativeFile);
%         xASL_spm_reslice( x.D.ResliceRef, NativeFile, [], [], x.Quality, 'native to standard space T1 ', StandardFile, 4 );
%     end
%
%     if  exist(StandardFile,'file')
%         T1IM(:,:,:,iS)  = xASL_io_Nifti2Im(StandardFile);
%     end
% end

%% Load 10 ASL images with lowest & highest atrophy
ASLIM    = zeros(121,145,121,x.nSubjects,'single');

for iS=1:x.nSubjects
    xASL_TrackProgress(iS,x.nSubjects);
    StandardFile    = fullfile(x.D.PopDir,['qCBF_untreated_' x.SUBJECTS{iS} '_ASL_1.nii']);
    StandardFileGZ  = [StandardFile '.gz'];
    if  exist(StandardFile,'file') || exist(StandardFileGZ,'file')
        ASLIM(:,:,:,iS)  = xASL_io_Nifti2Im(StandardFile);
    end
end

ASLimMax             = ASLIM(:,:,:,ColumnSort(end-10:end,2));
ASLimMin             = ASLIM(:,:,:,ColumnSort(1:20,2));

ASLimMin             = xASL_stat_MeanNan(ASLIM(:,:,:,[1 3 4 5 6 7 9 10 12 13 17 18 19 20]),4);

dip_image(ASLimMin)

% Also load template

TemplateFile        = fullfile(x.D.TemplateDir,'Philips_2DEPI_noBsup_CBF_Template.nii')
TemplateIM          = xASL_io_Nifti2Im(TemplateFile);
dip_image([ASLimMin+TemplateIM])
VeniLowCBF_Example  = 'C:\Backup\ASL\Veni\Presentation\VeniLowCBF_Example.nii';
xASL_io_SaveNifti(TemplateFile,VeniLowCBF_Example,[ASLimMin+TemplateIM]);

%% Flow territories
FlowTerr            = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Atlases\VascularTerritories\CortVascTerritoriesTatu.nii');
LeftMask            = ones(121,145,121);
LeftMask(1:61,:,:)  = 0;
LeftMask            = logical(LeftMask);
LeftMask(~FlowTerr) = 0;
FlowTerr(LeftMask)  = FlowTerr(LeftMask)+3;
dip_image(FlowTerr(:,:,53))



%% Sort them for spatial CoV
for iS=1:length(x.S.SetsName) % get which set is the GM-ICV ratio
    if  strcmp(x.S.SetsName{iS},'CBF_spatial_CoV')
        iSet    = iS;
    end
end

ColumnSort          = x.S.SetsID(:,iSet);
ColumnSort(:,2)     = [1:1:length(ColumnSort)];
ColumnSort          = sortrows(ColumnSort,1);
Mask                = ~isnan(ColumnSort);
ColumnSort          = ColumnSort(Mask(:,1),:);
Mask                = ColumnSort(:,1)>0 & ColumnSort(:,1)<1;
ColumnSort          = ColumnSort(Mask,:);

ASLIMMax             = ASLIM(:,:,:,ColumnSort(end-10:end,2));
ASLIMMin             = ASLIM(:,:,:,ColumnSort(1:10,2));

dip_image(ASLIMMin)

%% VeniCBFExample2
VeniLowCBF_Example  = 'C:\Backup\ASL\Veni\Presentation\VeniLowCBF_Example.nii';
VeniCBFExample2     = 'C:\Backup\ASL\Veni\Presentation\VeniCBFExample2.nii';
VascSCfile          = 'C:\Backup\ASL\Veni\Presentation\VascMeanFile.nii';
VascMaxFile         = 'C:\ExploreASL\Maps\Templates\MaxVesselTemplate.nii';
IndFile             = 'C:\Backup\ASL\Harmy\analysis_Harmy\dartel\qCBF_HD001_1_ASL_1.nii';
IndFile             = xASL_io_Nifti2Im(IndFile).*(xASL_io_Nifti2Im(IndFile)>0);

CombiIM             = xASL_io_Nifti2Im(VeniLowCBF_Example)+0.5.*xASL_io_Nifti2Im(VascMaxFile)+0.25.*IndFile;
xASL_io_SaveNifti(VeniLowCBF_Example,VeniCBFExample2,CombiIM);
