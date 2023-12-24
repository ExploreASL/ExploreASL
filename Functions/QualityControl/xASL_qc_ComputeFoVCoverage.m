function [CoveragePerc] = xASL_qc_ComputeFoVCoverage(InputPath, x)
%xASL_qc_ComputeFoVCoverage Compute coverage of low resolution image
%
% FORMAT: [CoveragePerc] = xASL_qc_ComputeFoVCoverage(InputPath, x)
% 
% INPUT:
%   InputPath   - path to low resolution image (REQUIRED)
%   x           - struct containing pipeline environment parameters (REQUIRED)
%
% OUTPUT:
%   CoveragePerc - intersectin low resolution image with high resolution
%                  T1w image (%). 100% = full coverage. 0 = no coverage
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function computes the intersection/overlap between
%              brainmask on field-of-view (FoV) of low resolution image
%              (native space) & the same brainmask with expanded FoV.
%              It uses the pGM+pWM+pCSF as brainmask
%              This assumes that the structural reference image has full brain coverage,
%              and was properly segmented into GM, WM and CSF
%              Also, we assume that the InputPath contains a single 3D volume
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: CoveragePerc = xASL_qc_ComputeFoVCoverage('//MyDisk/MyStudy/sub-001/ASL_1/mean_control.nii', x);
% __________________________________
% Copyright 2015-2019 ExploreASL

    CoveragePerc = NaN; % default, or for if the script crashes

    if ~xASL_exist(InputPath, 'file')
        warning('Input image didnt exist, skipping');
        return;
    end

    if ~xASL_exist(x.P.Path_c1T1, 'file') || ~xASL_exist(x.P.Path_c2T1, 'file') || ~xASL_exist(x.P.Path_c3T1, 'file')
        warning('Not all structural segmentations existed, skipping QC coverage calculation');
        return;
    end


    [Fpath, Ffile] = xASL_fileparts(InputPath);
    InputPathMask = fullfile(Fpath, [Ffile 'Mask.nii']);
    InputPathMaskFull = fullfile(Fpath, [Ffile '_MaskFull.nii']);

    % Get resolution of low resolution image
    nii = xASL_io_ReadNifti(InputPath);
    ImRes = nii.hdr.pixdim(2:4);
    MatSize = size(nii.dat);
    
    PathBrainmask = fullfile(x.dir.SUBJECTDIR, 'Brainmask.nii');
    BrainMask = xASL_io_Nifti2Im(x.P.Path_c1T1)+xASL_io_Nifti2Im(x.P.Path_c2T1)+xASL_io_Nifti2Im(x.P.Path_c3T1)>0.5;
    % Save
    xASL_io_SaveNifti(x.P.Path_c1T1, PathBrainmask, logical(BrainMask), 8, 0);
    
    % Reslice: get partial brainmask in low resolution space
    xASL_spm_reslice(InputPath, PathBrainmask, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, InputPathMask, 0);    

    % Add voxels to this space
    xASL_im_Upsample(InputPathMask, InputPathMask, ImRes,[],MatSize);
    % Reslice: get full brainmask in low resolution space
    xASL_spm_reslice(InputPathMask, PathBrainmask, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.settings.Quality, InputPathMaskFull, 0);     
    
    
    % Compute DICE coefficient
    IM1 = xASL_io_Nifti2Im(InputPathMask)>0;
    IM2 = xASL_io_Nifti2Im(InputPathMaskFull)>0;

    Intersection = IM1 & IM2;
    TotalSum = sum(IM1(:));
    CoveragePerc = 100*sum(Intersection(:)) / TotalSum;

    % Admin
    xASL_delete(InputPathMaskFull);
    xASL_delete(InputPathMask);
    xASL_delete(PathBrainmask);

    fprintf('%s\n', ['Coverage ' xASL_num2str(CoveragePerc) '% was computed']);
end
