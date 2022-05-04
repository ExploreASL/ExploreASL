clear all; spm('defaults', 'FMRI')
spm_file = [pwd '\SPM.mat']; load(spm_file);
for contNum =1:size(SPM.xCon,2)
write_thresholded_owen(spm_file,contNum,'cluster',0.001)
spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'SPM-filtered',sprintf('%s%s%s','NEWcontrast_',int2str(contNum),'_FWE_clusterwise'))
write_thresholded_owen(spm_file,contNum,'fwe',[])
spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'SPM-filtered',sprintf('%s%s%s','NEWcontrast_',int2str(contNum),'_FWE_voxelwise'))
end

