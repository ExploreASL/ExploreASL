x = ExploreASL_Master('',0);

%% Create non-deformed T1 in MNI
x.Quality=1;

if ismac
    Rdir = '/Volumes/nonshib-webdav/HolidayPics';
elseif ispc
    Rdir = 'C:\Users\kyrav\Desktop\SurfDrive\HolidayPics';
end

Rdir = fullfile(Rdir,'ExploreASL_ms','S2_PICTURE');


% PathIn = fullfile(Rdir, 'T1.nii');
% PathOut = fullfile(Rdir, 'rT1.nii');
% DeformationPath = fullfile(x.MyPath, 'Maps','Identity_Deformation_y_T1.nii');
% xASL_spm_deformations(x,PathIn,PathOut,[],[],[],DeformationPath);

%% Reslice other examples

% Path{1} = fullfile(Rdir, 'P01_S00631_enh_rDARTEL_pre_DEFmni_T1c.nii');
% Path{2} = fullfile(Rdir, 'P01_S00631_enh_rDARTELxasl_pre_DEFmni_T1c.nii');
%
% for ii=1:2
%     xASL_spm_reslice(x.D.ResliceRef, Path{ii});
% end


%% Create Figure
clear Path
Path{1} = 'rT1.nii';
Path{2} = 'rP01_S00631_enh_rDARTEL_pre_DEFmni_T1c.nii';
Path{3} = 'rP01_S00631_enh_rDARTELxasl_pre_DEFmni_T1c.nii';

for iPath=1:length(Path)
    IM(:,:,:,iPath) = xASL_io_Nifti2Im(fullfile(Rdir,Path{iPath}));
end

% x.S.TraSlices = [75:3:94];
x.S.TraSlices = [85]; % dummy
x.S.SagSlices = [];
x.S.CorSlices = [53:3:72];
x.S.CorSlices = x.S.CorSlices(4:5);
x.S.ConcatSliceDims = 1;
x.S.Square = 0;

clear Figure1
for ii=1:3
    Figure1{ii} = xASL_im_CreateVisualFig(x, {IM(:,:,:,ii)}, [], [1 0.75 0.75], [], {x.S.gray x.S.red x.S.blue}, 0);
end
Figure1 = [Figure1{1};Figure1{2};Figure1{3}];

%% Now create a sagittal registration slice
RootTemplate = fullfile(x.MyPath, 'External','SPMmodified','toolbox','cat12','templates_1.50mm');
Path{4} = 'Template_6_IXI555_MNI152.nii';
Template = xASL_io_Nifti2Im(fullfile(RootTemplate,Path{4}));
c1 = Template(:,:,:,1);
c2 = Template(:,:,:,2);

x.S.TraSlices = [53];
% x.S.CorSlices = [53]; % dummy
% x.S.SagSlices = []; % 63

clear Figure2
for ii=1:3
    Figure2{ii} = xASL_im_CreateVisualFig(x, {IM(:,:,:,ii) c2}, [], [1 0.75 0.75], [], {x.S.gray x.S.red x.S.blue }, 0);
end
Figure2 = [Figure2{1};Figure2{2};Figure2{3}];
figure(1);imshow([Figure1 Figure2],'InitialMagnification',250)
%
% figure(1);imshow([Figure1],'InitialMagnification',250)
% figure(1);imshow([Figure2],'InitialMagnification',250)
