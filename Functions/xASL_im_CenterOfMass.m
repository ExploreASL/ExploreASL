function xASL_im_CenterOfMass(niiPath,OtherList)
%xASL_im_CenterOfMass Determines center of mass and applies it to the NIfTI header orientation matrix

%% Admin

nii = xASL_io_ReadNifti(niiPath);
IM = single(nii.dat(:,:,:,1)).^0.5; % restore contrast
% & use first image!
% nii.mat = nii.mat0; % restore original orientation
% Don't reload mat0, this can differ between derivatives & original files

%% ---------------------------------------------------------------------------------
%% Compute current center voxel position from orientation matrix

% orientation matrix is in mm, so translation stays same with voxel-size, rotation, zooming &
% shearing. It needs to be divided by voxel-size to translate mm to nVoxels

% the determinant should remain the same (this has some negatives) which is why
% the spm function spm_imatrix is used. This is sign-flipped in the translation from mm to nVoxels.
% X-Y-Z order is the same for nVoxels & mm orientation matrix: L->R P->A I->S (from low to high values)

% QC
if sum(isfinite(nii.mat(:)))<numel(nii.mat(:))
    nii.mat = nii.mat0; % if there is any error in the orientation matrix, revert to the original orientation
end


if ~exist('OtherList','var')
    OtherList{1} = niiPath;
else
    OtherList{end+1} = niiPath;
end

fprintf('%s','Automatic alignment with center of mass to XYZ:');

%% ---------------------------------------------------------------------------------
%% Compute 3D Center of Mass (in pixel coordinates)

Dims = [1 2 3 1 2];

for iD=1:3
	Hist1Dim = xASL_stat_MeanNan(xASL_stat_MeanNan(IM,Dims(iD+1)),Dims(iD+2)); % get distribution
	Hist1Dim = reshape(Hist1Dim,[max(size(Hist1Dim)) 1]);

	Hist1Dim = Hist1Dim./xASL_stat_SumNan(Hist1Dim); % sum of distribution should be 1
	Hist1Dim(:,2) = [1:1:length(Hist1Dim)]; % voxel indices

	CoM(iD) = round(xASL_stat_SumNan(Hist1Dim(:,1).*Hist1Dim(:,2))); % weighted-average voxel index (==center of mass)
end

%% Calculate Center of Mass (in world coordinates)
CoM2 = - (nii.mat(1:3,1:3) * CoM');

%% Add the difference of AnteriorCommissure and CenterOfMass of a skullstripped brain
CoM2 = CoM2 + [0.8462;-17.5297;15.3505];

%% Calculate shift (in world coordinate)
CoMshift = CoM2 - nii.mat(1:3,4);

for iO=1:length(OtherList)
    if xASL_exist(OtherList{iO},'file')
        nii = xASL_io_ReadNifti(OtherList{iO});

        % QC
        if sum(isfinite(nii.mat(:)))<numel(nii.mat(:))
            nii.mat = nii.mat0; % if there is any error in the orientation matrix, revert to the original orientation
        end

        nii.mat(1:3,4) = nii.mat(1:3,4) + CoMshift;
        create(nii);

        %% Apply also to the .mat (motion) files, if they exist
        [Fpath, Ffile] = xASL_fileparts(OtherList{iO});
        MatFile = fullfile(Fpath,[Ffile '.mat']);

        if exist(MatFile,'file')
            S = load(MatFile,'mat');
            mat = S.mat;
            for iM=1:size(mat,3)
				% Don't modify the first matrix if it is zero
				if (iM > 1) || (sum(sum(abs(mat(:,:,1))))>0)
					mat(1:3,4,iM) = mat(1:3,4,iM) + CoMshift;
				end
            end
			save(MatFile,'mat');
        end
    end
end


fprintf('%s\n',[num2str(CoM(1),3) ' ' num2str(CoM(2),3) ' ' num2str(CoM(3),3)]);



end
