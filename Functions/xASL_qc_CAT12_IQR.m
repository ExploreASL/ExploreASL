function [QA_Output] = xASL_qc_CAT12_IQR(InputImage, InputC1, InputC2, InputC3, bFLAIR)
%xASL_qc_CAT12_IQR Prepare and run CAT12s QC parameters (also for other
%images)
%
% InputImage = image that will be QCed
% InputC1, 2, 3 = posterior probability maps (segmentations) for GM, WM,
% CSF
% Input should be paths
%
% EXAMPLE:
% InputImage = '/Users/henk/ExploreASL/Obesitas_Nijmegen/analysis/sub-001/T1.nii';
% InputC1 = '/Users/henk/ExploreASL/Obesitas_Nijmegen/analysis/sub-001/p1T1.nii';
% InputC2 = '/Users/henk/ExploreASL/Obesitas_Nijmegen/analysis/sub-001/p2T1.nii';
% InputC3 = '/Users/henk/ExploreASL/Obesitas_Nijmegen/analysis/sub-001/p3T1.nii';


if nargin<5 || isempty(bFLAIR)
    bFLAIR = false;
end
    
if isempty(which('cat_vol_qa'))
    error('Please install/initialize CAT12 for this QC function');
end 

fprintf('Obtaining CAT12 IQR parameter\n');

[Fpath, Ffile, Fext] = xASL_fileparts(InputImage);
PathM = fullfile(Fpath, ['m' Ffile Fext]); % biasfield-corrected image
PathSegm{1} = InputC1;
PathSegm{2} = InputC2;
PathSegm{3} = InputC3;

% Some other stuff for CAT12
job=struct;
job.extopts = struct;
job.extopts.species = 'human';
job.opts = struct;

cat_warnings = struct('identifier',{},'message',{});

%% Load segmentations & get volumes
fprintf('Loading segmentations (& reslicing if needed)\n');
niiMain = xASL_io_ReadNifti(InputImage); % Obtain orientation matrix for main image

for iSegm=1:length(PathSegm)
    niiC{iSegm} = xASL_io_ReadNifti(PathSegm{iSegm});
    if ~isequal(niiC{iSegm}.mat,niiMain.mat)
        % if segmentations don't have an equal size resample them
        refPath = InputImage;
        srcPath = PathSegm{iSegm};
        [Fpath, Ffile, Fext] = xASL_fileparts(srcPath);
        PathSegm{iSegm} = fullfile(Fpath, ['r' Ffile Fext]);
        
        VoxelSizeSegment = prod([niiC{iSegm}.mat(1,1) niiC{iSegm}.mat(2,2) niiC{iSegm}.mat(3,3)]);
        VoxelSizeOriginal = abs(prod([niiMain.mat(1,1) niiMain.mat(2,2) niiMain.mat(3,3)]));
        
        if (VoxelSizeSegment*4)<VoxelSizeOriginal
            % make sure to correctly downsample, if the resolutions differ significantly
            resRef = niiMain.hdr.pixdim(2:4);
            resSrc = niiC{iSegm}.hdr.pixdim(2:4);
            xASL_im_PreSmooth(refPath, srcPath, PathSegm{iSegm}, resRef, resSrc, [], []);
            xASL_spm_reslice(refPath, PathSegm{iSegm}, [], [], 1, PathSegm{iSegm}, 1);
        else
            xASL_spm_reslice(refPath, srcPath, [], [], 1, PathSegm{iSegm}, 1);
        end
    end
    
    IM_C{iSegm} = xASL_io_Nifti2Im(PathSegm{iSegm});
    IM_C{iSegm}(IM_C{iSegm}<0) = 0;
    IM_C{iSegm}(IM_C{iSegm}>1) = 1;
	IM_C{iSegm}(isnan(IM_C{iSegm})) = 0;
    
	% Get voxel size and calculate the volume
	niiVol = xASL_io_ReadNifti(PathSegm{iSegm});
	resVol = sqrt(sum((niiVol.mat(1:3,1:3)).^2,1));
    VolumeSegm{iSegm} = prod(resVol).*xASL_stat_SumNan(IM_C{iSegm}(:))./1000;
end
vol_TIV = VolumeSegm{1}+VolumeSegm{2}+VolumeSegm{3};
vol_rel_CGW = [VolumeSegm{3} VolumeSegm{1} VolumeSegm{2} eps eps]./vol_TIV;

% Fill the volumes for CAT12
qa = struct;
qa.subjectmeasures = struct;
qa.subjectmeasures.vol_TIV = vol_TIV;
qa.subjectmeasures.vol_rel_CGW = vol_rel_CGW;

%% Get biasfield-corrected image
fprintf('Obtaining biasfield information\n');
xASL_spm_BiasfieldCorrection(InputImage, [], false, [], PathM);

% Load biasfield NIfTI for CAT12
res = struct;
res.image.fname = PathM;
res.image.private = xASL_io_ReadNifti(PathM);
res.image.dim = size(res.image.private.dat);
res.image.mat = res.image.private.mat;
res.image.descrip = res.image.private.descrip;
res.image.n = [1 1];
res.image.dt = double([res.image.private.hdr.datatype 0]);
res.image.pinfo = [1;0;res.image.private.hdr.vox_offset];
Ym = xASL_io_Nifti2Im(PathM);

%% create combined label/segmentation map
fprintf('Creating combined label/segmentation map\n');
Yp0 = (IM_C{3}.*255)/5 + (IM_C{1}.*255)/5*2 + (IM_C{2}.*255)/5*3;
WBmask = (IM_C{1}+IM_C{2}+IM_C{3})>0.5;

% low resolution Yp0b
[indx, indy, indz] = ind2sub(size(Yp0),find(WBmask));
indx = max((min(indx) - 1),1):min((max(indx) + 1),size(Yp0,1));
indy = max((min(indy) - 1),1):min((max(indy) + 1),size(Yp0,2));
indz = max((min(indz) - 1),1):min((max(indz) + 1),size(Yp0,3));
Yp0b = Yp0(indx,indy,indz);

Yp0 = zeros(size(Yp0),'single'); 
Yp0(indx,indy,indz) = single(Yp0b)/255*5; 
Yp0(Yp0>3.00001) = nan;

%% Run the CAT12 QC function
fprintf('Running CAT12 QC function\n');
if bFLAIR
    qa = cat_vol_qa('cat12',Yp0,InputImage,Ym,res,cat_warnings,job.extopts.species, ...
              struct('write_csv',0,'write_xml',1,'method','cat12','job',job,'qa',qa,'seqtype','flair'));
else
    qa = cat_vol_qa('cat12',Yp0,InputImage,Ym,res,cat_warnings,job.extopts.species, ...
              struct('write_csv',0,'write_xml',1,'method','cat12','job',job,'qa',qa));
end
      
QA_Output = qa.qualityratings;

%% Householding
xASL_delete(PathM);
if ~strcmp(PathSegm{1}, InputC1)
    xASL_delete(PathSegm{1});
end
if ~strcmp(PathSegm{2}, InputC2)
    xASL_delete(PathSegm{2});
end
if ~strcmp(PathSegm{3}, InputC3)
    xASL_delete(PathSegm{3});
end

end

