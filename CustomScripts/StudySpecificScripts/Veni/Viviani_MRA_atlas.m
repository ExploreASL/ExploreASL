%% Resample Viviani MRA atlas into xASL

FilePath     = 'C:\ExploreASL\Maps\Atlases\VesselDensityHR.nii';
FilePath2   = 'C:\ExploreASL\Maps\Atlases\VesselDensityHR2.nii';
nii          = xASL_io_Nifti2Im(FilePath);
nii          = nii./max(nii(:));

for ii=1:8
    nii      = xASL_im_ndnanfilter(nii,'gauss',[1.885 1.885 1.885]);
end

xASL_io_SaveNifti(FilePath,FilePath2,nii);

xASL_spm_reslice( x.D.ResliceRef, FilePath2, [], [], x.Quality, [], 'C:\ExploreASL\Maps\Atlases\new.nii', 4 )

FilePath    = 'C:\ExploreASL\Maps\Atlases\MRA_Viviani_xASL.nii';
nii         = xASL_io_Nifti2Im(FilePath);
nii         = nii./max(nii(:));
xASL_io_SaveNifti(FilePath,FilePath,nii);


FilePath    = 'C:\ExploreASL\Maps\Atlases\MRA_Viviani_xASL.nii';
FilePath2   = 'C:\ExploreASL\Maps\Atlases\MRA_Viviani_xASL2.nii';
nii         = xASL_io_Nifti2Im(FilePath);
nii         = exp(nii);
nii         = exp(nii);
nii(nii>exp(1)*2.5)     = exp(1)*2.5;
nii         = nii-min(nii(:));
nii         = nii./max(nii(:));
nii         = nii.^0.5;

for ii=1:2
    nii      = xASL_im_ndnanfilter(nii,'gauss',[1.885 1.885 1.885]);
end

nii         = nii./max(nii(:));

xASL_io_SaveNifti(FilePath,FilePath2,nii);
