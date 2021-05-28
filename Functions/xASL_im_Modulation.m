function xASL_im_Modulation(x)
%xASL_im_Modulation Combines the transformations to create Jacobians, &
% multiplies the standard space segmentations with these to create volumetric
% images for volumetric analyses
%
% FORMAT:       xASL_im_Modulation(x)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Combines the transformations to create Jacobians, &
%               multiplies the standard space segmentations with these to create volumetric
%               images for volumetric analyses.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

fprintf('%s\n','Modulating PV images in standard space');

%% ------------------------------------------------------------------------------
%% 1) Create Jacobian determinants
% Assuming that the SPM flow field was calculated from the current position and only rigid-realign was applied

% Load the matrix of the T1w file
tNii = xASL_io_ReadNifti(x.P.Path_T1);
OriginalMat = tNii.mat0;
CurrentMat = tNii.mat;
DiffMat = CurrentMat/OriginalMat;

% Checks if there's no scaling between the mat0 and mat
for ii=1:3
	nr = max(norm(DiffMat(1:3,ii)),norm(DiffMat(ii,1:3)));
	if abs(nr-1) > 0.001
		warning('xASL_im_Modulation::Scaling during modulation');
	end
end

% Creates a blank image
%[fpath, fname] = xASL_fileparts(x.P.Path_T1);
%tNii.dat.fname = fullfile(fpath,[fname '_tmp.nii']);
%create(tNii);
%tNii.dat(:,:,:,:,:)  = single(ones(tNii.dat.dim));

% Run the SPM script for creating the determinants
xASL_adm_UnzipNifti(x.P.Path_y_T1);
matlabbatch{1}.spm.util.defs.comp{1}.def = {x.P.Path_y_T1};
matlabbatch{1}.spm.util.defs.out{1}.savejac.ofname = 'T1.nii';
matlabbatch{1}.spm.util.defs.out{1}.savejac.savedir.saveusr = {x.dir.SUBJECTDIR};
% For testing create also the transfile
%matlabbatch{2}.spm.util.defs.out{1}.pull.fnames = {tNii.dat.fname};
%matlabbatch{2}.spm.util.defs.out{1}.pull.savedir.saveusr = {x.dir.SUBJECTDIR};
%matlabbatch{2}.spm.util.defs.out{1}.pull.interp = 4;
%matlabbatch{2}.spm.util.defs.out{1}.pull.mask = 1;
%matlabbatch{2}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
%matlabbatch{2}.spm.util.defs.out{1}.pull.prefix = '';
spm_jobman('run', matlabbatch);


% Remove border artifacts
jIM = xASL_io_Nifti2Im(x.P.Path_j_T1);
jIM([1:5 end-5:end],:,:) = NaN; % assuming 121x145x121 this would be OK
jIM(:,[1:5 end-5:end],:) = NaN;
jIM(:,:,[1:5 end-5:end]) = NaN;
jIM(jIM<0)               = NaN;
while sum(sum(sum(isnan(jIM))))~=0 % the 6 window should be fine in 1 go, but just to be sure
    jIM = xASL_im_ndnanfilter(jIM, 'rect', [6 6 6], 2);
end
xASL_io_SaveNifti(x.P.Path_j_T1, x.P.Path_j_T1, jIM, [], 0); % this can be removed, just that you can see the result

%% ------------------------------------------------------------------------------
%% 2) Create modulated images by multiplying probability/PV images by Jacobian determinants
TissueType = {'GM' 'WM' 'CSF' 'WMH_SEGM'};
NativeList = {x.P.Path_c1T1 x.P.Path_c2T1 x.P.Path_c3T1 x.P.Path_WMH_SEGM};
StandardList = {x.P.Pop_Path_rc1T1 x.P.Pop_Path_rc2T1 x.P.Pop_Path_rc3T1 x.P.Pop_Path_rWMH_SEGM};
ModulateList = {x.P.Pop_Path_mrc1T1 x.P.Pop_Path_mrc2T1 x.P.Pop_Path_mrc3T1 x.P.Pop_Path_mrWMH_SEGM};
bDebug = true;

for iL=1:length(StandardList)
    if xASL_exist(StandardList{iL},'file')
        tIM = xASL_io_Nifti2Im(StandardList{iL}) .* jIM;

%         if bDebug
        if ~xASL_exist(NativeList{iL}, 'file')
            warning('Couldnt determine volumetric difference native<>standard spaces, native space image didnt exist');
        else
            tIMOR = xASL_io_Nifti2Im(NativeList{iL});
            tNii = xASL_io_ReadNifti(NativeList{iL});
            OriginalMat = tNii.mat0;
            volNat = round(xASL_stat_SumNan(tIMOR(:)) * norm(OriginalMat(1:3,1)) * norm(OriginalMat(1:3,2)) * norm(OriginalMat(1:3,3)));
            volStd = round(xASL_stat_SumNan(tIM(:)) * 1.5^3); % current default standard space resolution for ExploreASL
            Difference = xASL_round(100 * ((volStd-volNat) / volNat), 2);
            fprintf('%s\n', [TissueType{iL} ' volume in native/standard space: ' num2str(volNat) ', ' num2str(volStd) '. Difference ' num2str(Difference) '%']);
            % Do scaling
%             if abs(Difference)>5 % in percentages, correct when larger difference than 5%
%                 % This is not the case for pGM or pWM or pCSF, but is the case for e.g. pWMH_SEGM
%                 [N1, X1] = hist(nonzeros(tIM),100);
%                 [N2, X2] = hist(nonzeros(tIMOR),100);
%
%                 N1 = N1./sum(N1);
%                 N2 = N2./sum(N2);
%                 X1 = [0 X1]; N1 = [0 N1];
%                 X2 = [0 X2]; N2 = [0 N2];
%
%                 for iS=1:length(X1)-1
%                     ScaleF = X2(iS+1)/X1(iS+1);
%                     tROI = tIM>X1(iS) & tIM<X1(iS+1);
%                     tIM(tROI) = tIM(tROI).*ScaleF;
%                 end

%
%                 volStd = round(xASL_stat_SumNan(tIM(:)) * 1.5^3); % current default standard space resolution for ExploreASL
%                 Difference = xASL_round(100 * ((volStd-volNat) / volNat), 2);
%                 fprintf('%s\n', ['After scaling: ' num2str(volNat) ', ' num2str(volStd) '. Difference ' num2str(Difference) '%']);
%             end
             xASL_io_SaveNifti(StandardList{iL}, ModulateList{iL}, tIM, [], false);
        end
    end
end

if x.DELETETEMP
    xASL_delete(x.P.Path_j_T1);
end

end
