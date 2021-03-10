%% Compress NIfTI masks for compilation
x = ExploreASL_Master('',0);

%% Convert to .nii.mat, keep .nii.gz for visualization (uint8 for masks)
Flist = {'C:\ExploreASL\Maps\WBmaskASL.nii';'C:\ExploreASL\Maps\WBmaskASLnarrow.nii';'C:\ExploreASL\Maps\Atlases\DeepWM.nii';'C:\ExploreASL\Maps\Atlases\HOcort_CONN.nii';'C:\ExploreASL\Maps\Atlases\HOsub_CONN.nii';'C:\ExploreASL\Maps\Atlases\Hammers.nii';'C:\ExploreASL\Maps\Atlases\LeftRight.nii';'C:\ExploreASL\Maps\Atlases\MNI_Structural.nii';'C:\ExploreASL\Maps\Atlases\Thalamus.nii';'C:\ExploreASL\Maps\Atlases\TotalGM.nii';'C:\ExploreASL\Maps\Atlases\TotalWM.nii';'C:\ExploreASL\Maps\Atlases\WholeBrain.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\ATTbasedFlowTerritories.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\BorderzoneMNI_Hartkamp2017JCBFM_3mmDilated.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\BorderzoneMNI_Hartkamp2017JCBFM_Original.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\CortVascTerritoriesTatu.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\LabelingTerritories.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\TatuICA_PCA.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\Watershed_Conventional.nii';'C:\ExploreASL\Maps\Atlases\VascularTerritories\vascRegionsAlternative.nii'};

for iFile=1:length(Flist)
    IM = uint8(xASL_io_Nifti2Im(Flist{iFile}));
    [Fpath, Ffile] = xASL_fileparts(Flist{iFile});
    MatPath = fullfile(Fpath, [Ffile '.nii.mat']);
    save(MatPath,'IM');
    gzip(Flist{iFile});
    delete(Flist{iFile});
end
    
%% Convert to .nii.mat, keep .nii.gz for visualization (single for probability maps)
Flist = {'C:\ExploreASL\Maps\Templates\Philips_2DEPI_Bsup_CBF.nii';'C:\ExploreASL\Maps\Templates\Philips_2DEPI_noBsup_Control.nii';'C:\ExploreASL\Maps\Templates\Susceptibility_pSignal_2D_EPI.nii';'C:\ExploreASL\Maps\Templates\Susceptibility_pSignal_3D_GRASE.nii'};
% first two are needed for xASL_qc_CollectQC_ASL ->
% xASL_qc_CompareTemplate, which uses these as templates

for iFile=1:length(Flist)
    IM = xASL_io_Nifti2Im(Flist{iFile});
    [Fpath, Ffile] = xASL_fileparts(Flist{iFile});
    MatPath = fullfile(Fpath, [Ffile '.nii.mat']);
    save(MatPath,'IM');
    gzip(Flist{iFile});
    delete(Flist{iFile});
end

%% Convert to .nii.gz, will be used by xASL_spm_Deformations, which will copy to data folder anyway
Flist = {'C:\ExploreASL\Maps\Templates\ATT_BiasField.nii';'C:\ExploreASL\Maps\Templates\MaxVesselTemplate.nii';'C:\ExploreASL\Maps\Templates\GE_3Dspiral_Product_CBF.nii';'C:\ExploreASL\Maps\Templates\GE_3Dspiral_QC_mask.nii';'C:\ExploreASL\Maps\Templates\Philips_2DEPI_Bsup_CBF.nii';'C:\ExploreASL\Maps\Templates\Philips_2DEPI_Bsup_QC_mask.nii';'C:\ExploreASL\Maps\Templates\Philips_2DEPI_noBsup_CBF.nii';'C:\ExploreASL\Maps\Templates\Philips_2DEPI_noBsup_Control.nii';'C:\ExploreASL\Maps\Templates\Siemens_2DEPI_PCASL_noBsup_CBF.nii';'C:\ExploreASL\Maps\Templates\Siemens_2DEPI_PCASL_noBsup_Control.nii';'C:\ExploreASL\Maps\Templates\Siemens_3DGRASE_PASL_CBF.nii';'C:\ExploreASL\Maps\Templates\Siemens_3DGRASE_PASL_QC_mask.nii';'C:\ExploreASL\Maps\Templates\Siemens_3DGRASE_PCASL_Control_BiasfieldCorr_MoodStudy.nii';'C:\ExploreASL\Maps\CentralWM_QC.nii';'C:\ExploreASL\Maps\Templates\brainmask_supratentorial.nii'};

for iFile=1:length(Flist)
    gzip(Flist{iFile});
    delete(Flist{iFile});
end