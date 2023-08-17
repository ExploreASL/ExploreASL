function xASL_im_CenterOfMass(PathNIfTI, OtherList, AllowedDistance)
%xASL_im_CenterOfMass Determine and apply center of mass
%
% FORMAT: xASL_im_CenterOfMass(PathNIfTI, OtherList, AllowedDistance)
%
% INPUT:
%   PathNIfTI       - path to NIfTI file on which to estimate center of mass (REQUIRED)
%   OtherList       - Cell structure with paths to other NIfTI images that
%                     should remain in alignment (OPTIONAL)
%   AllowedDistance - scalar of distance (mm) above which realignment is
%                     applied (OPTIONAL, DEFAULT=50mm)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function estimates the center of mass of the image
%              matrix, and if this is too far off the current orientation
%              matrix center, the center will be reset.
%              This fixes any incorrect orientation outputted by the
%              scanner.
%              The realignment is only applied when any of the X/Y/Z
%              dimensions have a higher offset than AllowedDistance.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:      T1w: xASL_im_CenterOfMass('Path2Study/sub-001/T1.nii', {'Path2Study/sub-001/FLAIR.nii'}, 0);
%               ASL: xASL_im_CenterOfMass('Path2Study/sub-001/ASL_1/ASL4D.nii', {'Path2Study/sub-001/ASL_1/M0.nii'}, 50);
% __________________________________
% Copyright (C) 2015-2020 ExploreASL



%% Admin

if nargin<2 || isempty(OtherList) || isempty(OtherList{1})
    OtherList{1} = PathNIfTI;
else
    OtherList = xASL_adm_OtherListSPM(OtherList, false); 
    % here we want a NIfTI with multiple volumes mentioned only once in the
    % list, as below we deal with the "other/n>1" volumes
    for iList=1:length(OtherList)
        OtherList{iList} = OtherList{iList}(1:end-2);
    end
    if isempty(OtherList) || isempty(OtherList{1})
        OtherList{1} = PathNIfTI;
    else
        OtherList{end+1} = PathNIfTI;
    end
    % We only apply the transformation to the OtherList,
    % So the last of OtherList is PathNIfTI
end

if nargin<3 || isempty(AllowedDistance)
    % 50 mm is default minimal offset we need to apply the new center of mass
    % per reviewer comment to avoid unnecessary loosing initial alignment
    % between scans from the same subject/visit
    AllowedDistance = 50;
end

nii = xASL_io_ReadNifti(PathNIfTI);
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
    warning(['There was an error in the NIfTI orientation matrix, trying to repair: ' PathNIfTI]);
    nii.mat = nii.mat0; % if there is any error in the orientation matrix, revert to the original orientation
end



fprintf('%s','Automatic alignment with center of mass:');

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

%% Check if alignment to Center of Mass is needed
% if not, then skip
if norm(CoMshift)<AllowedDistance
    fprintf(' skipping, was too close\n');
    return;
end
    
%% Apply alignment

for iO=1:length(OtherList)
    if xASL_exist(OtherList{iO},'file')
        
        nii = xASL_io_ReadNifti(OtherList{iO});

        % QC
        if sum(isfinite(nii.mat(:)))<numel(nii.mat(:))
            warning(['There was an error in the NIfTI orientation matrix, trying to repair: ' OtherList{iO}]);
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


fprintf('%s\n',[' to XYZ:' num2str(CoM(1),3) ' ' num2str(CoM(2),3) ' ' num2str(CoM(3),3)]);



end
