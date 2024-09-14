% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

clear all; spm('defaults', 'FMRI')
spm_file = [pwd '\SPM.mat']; load(spm_file);
for contNum =1:size(SPM.xCon,2)
write_thresholded_owen(spm_file,contNum,'cluster',0.001)
spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'SPM-filtered',sprintf('%s%s%s','NEWcontrast_',int2str(contNum),'_FWE_clusterwise'))
write_thresholded_owen(spm_file,contNum,'fwe',[])
spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'SPM-filtered',sprintf('%s%s%s','NEWcontrast_',int2str(contNum),'_FWE_voxelwise'))
end

