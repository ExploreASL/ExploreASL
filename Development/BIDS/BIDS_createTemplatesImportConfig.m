%% Create JSON files from ExploreASL_ImportConfig 

% Initialize ExploreASL
[x] = ExploreASL_Initialize;

% Define list
listOfConfigFlavors = { ...
                       'RACE5_Siemens',...
                       'RACE5_Philips',...
                       'TimFit',...
                       'Insight46',...
                       'Craniosynostosis',...
                       'SABRE',...
                       'SydneyMS_Controls',...
                       'SydneyMS_Controls',...
                       'MS_Sidney',...
                       'NOMARED',...
                       'Leiden',...
                       'WSU',...
                       'session_by_date',...
                       'Obesitas_Nijmegen',...
                       'RADAR_Sheffield',...
                       'RADAR_Glasgow',...
                       'RADAR_Bristol',...
                       'NoseWorthy',...
                       'InsomniaHC',...
                       'QSMSC',...
                       'QSMSC_MultiPhase',...
                       'HongerWinter2',...
                       'Dent',...
                       'Miami_GE',...
                       'Miami_Siemens_Rockland',...
                       'SPIR_artifact_Koen',...
                       'ExploreASLtest',...
                       'Trial_patient',...
                       'GENFI_Test',...
                       'RunDMC_inTENse',...
                       'MAGNIMS_Verona',...
                       'Frontier',...
                       'FRONTIER',...
                       'BrnoEpilepsy',...
                       'BoleStudien',...
                       'Sleep2',...
                       'HIFU_ASL_Vera',...
                       'Sleep_2018',...
                       'Sleep_2018_2',...
                       'CP_Tavi',...
                       'Iris4Jan_LongASL',...
                       'Maarten_Lequin',...
                       'Hardy',...
                       'Hardy2',...
                       'Chile',...
                       'APGEM',...
                       'Divers_Bonn',...
                       'GSP_perfusion_phantom',...
                       'Nijmegen_RunDMC_Serial_Imaging',...
                       'Nijmegen_RunDMC_Serial_Imaging2',...
                       'Michiel_Utrecht',...
                       'ACTL',...
                       'Atle_WIP_Siemens_3DGRASE',...
                       'MEDIRI',...
                       'Iris_unilateral_sclerosis',...
                       'Iris_unilateral_sclerosis2',...
                       'WRAP',...
                       'VICI',...
                       'BioCog_Repro',...
                       '22q11',...
                       'BioCog',...
                       'BioCogSiemens',...
                       'XNAT_BioCog',...
                       'DDI_IVS',...
                       'Hongerwinter',...
                       'CPC_NEURO_NPH',...
                       'Siemens_PASL_DCM',...
                       'GENFI_T1MatrixTrial',...
                       'CRUISE_pCASL_artefact',...
                       'GENFI_DF12',...
                       'GIFMI_STUMOOD',...
                       'FollowUp',...
                       'FInd',...
                       'CRUISE',...
                       '3D_Antipsychotics',...
                       'Trial_ePOD_Paul',...
                       'ASL_epod_MPH',...
                       'Score',...
                       'Score_2',...
                       'ePOD_AMPH',...
                       'Novice',...
                       'NOVICE',...
                       'NOVICE2',...
                       'Antipsychotics',...
                       'preDIVA_followUp1',...
                       'preDIVA_followUp2',...
                       'preDIVA_followUp3',...
                       'preDIVA_followUp4',...
                       'preDIVA_baselinePipeline',...
                       'preDIVA_baselinePipeline2',...
                       'GE_trial',...
                       'COBRA_ICL_BL',...
                       'COBRA_ICL_FU',...
                       'COBRA_AMC',...
                       'Parelsnoer_ASL_HJ',...
                       'INOX',...
                       'INOX2',...
                       'FIND_STUDIE',...
                       'WMH_AGEIV',...
                       'Sagittal_Sinus',...
                       'ADNI',...
                       'incoming',...
                        };

% Select CustomScripts directory
customScripts = uigetdir(x.MyPath,'Please select the CustomScripts directory...');

% Create JSON files
for iFlavor=1:length(listOfConfigFlavors)
    fprintf('%s\n', listOfConfigFlavors{iFlavor});
    imPar = ExploreASL_ImportConfig(listOfConfigFlavors{iFlavor});
    validFileName = [genvarname(listOfConfigFlavors{iFlavor}) '.json'];
    spm_jsonwrite(fullfile(customScripts, 'ConfigFiles', validFileName), imPar);
end

