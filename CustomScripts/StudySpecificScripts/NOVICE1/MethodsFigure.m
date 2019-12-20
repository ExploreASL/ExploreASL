%% Admin

x.MYPATH   = 'c:\ASL_pipeline_HJ';
AdditionalToolboxDir    = 'C:\ASL_pipeline_HJ_toolboxes'; % provide here ROOT directory of other toolboxes used by this pipeline, such as dip_image & SPM12
if ~isdeployed
    addpath(x.MYPATH);

    subfolders_to_add = { 'ANALYZE_module_scripts', 'ASL_module_scripts', fullfile('Development','dicomtools'), fullfile('Development','Filter_Scripts_JanCheck'), 'MASTER_scripts', 'spm_jobs','spmwrapperlib' };
    for ii=1:length(subfolders_to_add)
        addpath(fullfile(x.MYPATH,subfolders_to_add{ii}));
    end
end

addpath(fullfile(AdditionalToolboxDir,'DIP','common','dipimage'));

[x.SPMDIR, x.SPMVERSION] = xASL_adm_CheckSPM('FMRI',fullfile(AdditionalToolboxDir,'spm12') );
addpath( fullfile(AdditionalToolboxDir,'spm12','compat') );

if isempty(which('dip_initialise'))
    fprintf('%s\n','CAVE: Please install dip_image toolbox!!!');
else dip_initialise
end

%% AgeIV ASL paper

% For HIV
clear
ROOT            = 'D:\Backup\ASL_E\NoviceNew\MethodsFigure';

% this should be native space
Name{1}         = fullfile(ROOT, 'srT1_NOV001.nii');
Name{2}         = fullfile(ROOT, 'DARTEL_c1T1_NOV001.nii');
Name{3}         = fullfile(ROOT, 'DARTEL_c2T1_NOV001.nii');
Name{4}         = fullfile(ROOT, 'DARTEL_CBF_HctCorr_cohort_NOV001_ASL_1.nii');

DARTELdir                       = 'D:\Backup\ASL_E\NoviceNew\analysis\dartel';
HOSubMasks4DTotalFile           = fullfile( DARTELdir, 'HOSubMasks.dat');

HOSubMasks                      = memmapfile(HOSubMasks4DTotalFile,'Format',{'uint8' [121 145 121 62 28] 'data'});
Thalamus                        = HOSubMasks.Data.data(:,:,:,1,4)+HOSubMasks.Data.data(:,:,:,1,15);
Caudate                         = HOSubMasks.Data.data(:,:,:,1,5)+HOSubMasks.Data.data(:,:,:,1,16);
Putamen                         = HOSubMasks.Data.data(:,:,:,1,6)+HOSubMasks.Data.data(:,:,:,1,17);
Accumbens                       = HOSubMasks.Data.data(:,:,:,1,11)+HOSubMasks.Data.data(:,:,:,1,21);

% dip_image([Thalamus Caudate Putamen Accumbens])

for ii=1:4
    IM{ii}              = xASL_io_ReadNifti( Name{ii} );
end

IM{5}.dat   = Thalamus;
IM{6}.dat   = Caudate;
IM{7}.dat   = Putamen;
IM{8}.dat   = Accumbens;

clear Image
for ii=1:8
    Image{ii}(:,:,1)    = IM{ii}.dat(:,:,47);
    Image{ii}(:,:,2)    = IM{ii}.dat(:,:,53);
    Image{ii}(:,:,3)    = IM{ii}.dat(:,:,60);

    for iIm=1:3
        Temp(:,:,iIm)   = xASL_im_rotate(Image{ii}(:,:,iIm),90);
    end
    Image{ii}           = Temp;
    clear Temp
    VertIm{ii}                          = [Image{ii}(:,:,1);Image{ii}(:,:,2);Image{ii}(:,:,3)];
    VertIm{ii}                          = VertIm{ii}./max(VertIm{ii}(:));
    VertIm{ii}(VertIm{ii}<0)            = 0;
    VertIm{ii}(isnan(VertIm{ii}))       = 0;
end

jet256         = jet(256);
jet256(1,:)    = 0;

GMmask          = VertIm{2}>0.8*(max(VertIm{2}(:)));
WMmask          = VertIm{3}>0.8*(max(VertIm{3}(:)));
WMmask          = imerode(WMmask,strel('disk',5));

figure(1);imshow( [VertIm{1} VertIm{2} VertIm{3}],[],'InitialMagnification',200)
figure(2);imshow( [ GMmask  WMmask VertIm{5} VertIm{6} VertIm{7} VertIm{8}],[],'InitialMagnification',200)
figure(3);imshow( [VertIm{4} VertIm{4}.*GMmask VertIm{4}.*WMmask ]./mean(VertIm{4}(GMmask)).*65.3,[],'InitialMagnification',200,'Colormap',jet256)



    % 49 -> subcortex, anterior cingulate
    % 63 -> insula & posterior cingulate
    % 84 -> top of brain
