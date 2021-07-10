%% Obtain results for Norman

% Set defaults
x.S.InputDataStr = 'qCBF_untreated';
x.S.bMasking = false;

AtlasPath{1} = fullfile(x.D.MapsSPMmodifiedDir,'TotalGM.nii');
AtlasPath{2} = fullfile(x.D.MapsSPMmodifiedDir,'DeepWM.nii');
AtlasPath{3} = fullfile(x.D.MapsSPMmodifiedDir,'MNI_Structural.nii');
AtlasPath{4} = fullfile(x.D.AtlasDir,'Hammers.nii');
AtlasPath{5} = fullfile(x.D.AtlasDir,'HOcort_CONN.nii');
AtlasPath{6} = fullfile(x.D.AtlasDir,'HOsub_CONN.nii');

x.Sequence = '3D_spiral'; % to fool the Population module for pooling sequences
x.Q.readoutDim = '3D';


% 1) For the reproducibility study
x = ExploreASL_Initialize('/Users/henk/ExploreASL/ASL/BioCog_Repro/analysis/DATA_PAR.m');

for iAtlas=1:length(AtlasPath)
    x.S.InputAtlasPath = AtlasPath{iAtlas};
    xASL_wrp_GetROIstatistics(x);
end

% 2) For the controls inn the main study
x = ExploreASL_Initialize('/Users/henk/ExploreASL/ASL/BioCog_Repro/analysis_Controls/DATA_PAR_All.m');

for iAtlas=1:length(AtlasPath)
    x.S.InputAtlasPath = AtlasPath{iAtlas};
    xASL_wrp_GetROIstatistics(x);
end