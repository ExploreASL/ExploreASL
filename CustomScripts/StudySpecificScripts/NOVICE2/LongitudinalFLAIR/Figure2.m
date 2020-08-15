%% Figure 2
FLAIRpath = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis\Population\Templates\FLAIR_bs-mean.nii';
Path_WMH{1} = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis\Population\Templates\WMH_HC1.nii';
Path_WMH{2} = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis\Population\Templates\WMH_HC2.nii';
Path_WMH{3} = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis\Population\Templates\WMH_HIV1.nii';
Path_WMH{4} = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis\Population\Templates\WMH_HIV2.nii';

% FLAIRim = xASL_io_Nifti2Im(FLAIRpath);
MaskIM = xASL_io_Nifti2Im('C:\ExploreASL\Maps\SPM_xASLModified\rbrainmask.nii')>0.9;
% dip_image(FLAIRim)

for iPath=1:4
    WMHim{iPath} = xASL_io_Nifti2Im(Path_WMH{iPath});
    WMHim{iPath}(WMHim{iPath}<0) = 0;
    WMHim{iPath}(WMHim{iPath}>0.1) = 0.1;
end

x.S.TraSlices = x.S.slicesLarge(5:8);
x.S.Square = false;

for iPath=1:4
    ImOut{iPath} = xASL_vis_CreateVisualFig(x, {FLAIRpath WMHim{iPath}}, [], [1 0.4], [], {x.S.gray x.S.jet256}, 0, {MaskIM MaskIM}, 1, [], 5);
end

figure(1);imshow([ImOut{1};ImOut{2};ImOut{3};ImOut{4}])


%% Same as above, but with different FLAIR for each cohort
%% Figure 2
x = ExploreASL_Master('',0);
TemplateRoot = '/Users/henk/Google Drive/DataCopy/NoviceFigure2';
% TemplateRoot = '/home/hjmutsaerts/lood_storage/divi/Projects/novice/MRI/ExploreASLoutput/LongitudinalAnalysis/Population/Templates/';
% TemplateRoot = 'C:\BackupWork\ASL\Novice\LongitudinalAnalysis\Population\Templates';

Path_FLAIR{1} = fullfile(TemplateRoot, 'FLAIR_HC1.nii');
Path_FLAIR{2} = fullfile(TemplateRoot, 'FLAIR_HC2.nii');
Path_FLAIR{3} = fullfile(TemplateRoot, 'FLAIR_HIV1.nii');
Path_FLAIR{4} = fullfile(TemplateRoot, 'FLAIR_HIV2.nii');
Path_WMH{1} = fullfile(TemplateRoot, 'WMH_HC1.nii');
Path_WMH{2} = fullfile(TemplateRoot, 'WMH_HC2.nii');
Path_WMH{3} = fullfile(TemplateRoot, 'WMH_HIV1.nii');
Path_WMH{4} = fullfile(TemplateRoot, 'WMH_HIV2.nii');

% FLAIRim = xASL_io_Nifti2Im(FLAIRpath);
MaskIM = xASL_io_Nifti2Im(fullfile(x.MyPath,'External','SPMmodified','MapsAdded','rbrainmask.nii'))>0.9;
% dip_image(FLAIRim)

for iPath=1:4
    FLAIRim{iPath} = xASL_io_Nifti2Im(Path_FLAIR{iPath});
    WMHim{iPath} = xASL_io_Nifti2Im(Path_WMH{iPath});
    WMHim{iPath} = WMHim{iPath}>0;
    % let's dilate the regions a bit -> not required
%     ROIs = unique(WMHim{iPath});
%     ROIs = ROIs(2:end);
%     NewIm = zeros(size(WMHim{iPath}));
%     for iRoi=1:length(ROIs)
%         TempMask = WMHim{iPath}==ROIs(iRoi);
%         TempMask = xASL_im_DilateErodeFull(TempMask,'dilate',xASL_im_DilateErodeSphere(1));
%         NewIm = NewIm+(double(TempMask).*ROIs(iRoi));
%     end
%     WMHim{iPath} = NewIm;
end

x.S.TraSlices = x.S.slicesLarge(5:8);
x.S.Square = false;

for iPath=1:4
    if max(WMHim{iPath}(:))==0
%         ImOut{iPath} = xASL_vis_CreateVisualFig(x, {FLAIRim{iPath}}, [], [1], [], {x.S.gray}, 0, {MaskIM}, 0, [], 5);
        ImOut{iPath} = xASL_vis_CreateVisualFig(x, {FLAIRim{iPath}}, [], [1], [], {x.S.gray}, 0, [], 0, [], 5);
    else
%         ImOut{iPath} = xASL_vis_CreateVisualFig(x, {FLAIRim{iPath} WMHim{iPath}}, [], [1 1], [], {x.S.gray x.S.jet256}, 0, {MaskIM MaskIM}, 0, [], 5);
        ImOut{iPath} = xASL_vis_CreateVisualFig(x, {FLAIRim{iPath} WMHim{iPath}}, [], [1 1], [], {x.S.gray x.S.jet256}, 0, [], 0, [], 5);
    end
end

figure(1);imshow([ImOut{1};ImOut{2};ImOut{3};ImOut{4}],'InitialMagnification',250)
