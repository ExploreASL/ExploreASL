function xQC_Master (AnalysisDir, xASLDir, Regexp, bStructural, bFunctional, bDWI, bASL)


ExploreQC = pwd;

%handle REGexp

% First initialize exploreASL to obtain x, because some functions still
% need the x structure to work (deformations)
% Also add the paths to external tools(dti denoising)

cd (xASLDir);
x = ExploreASL_Master('', 0);
cd(ExploreQC)

% TEMPORARY : Add the path to scripts in development to test them
addpath(fullfile(ExploreQC, 'Development/mppca_denoise-master/mppca_denoise-master'))
addpath(fullfile(ExploreQC, 'Development/'))


Subjects = xASL_adm_GetFsList(AnalysisDir, '', 1);

%Preliminary settings, can be change by the user
SPMDir = fullfile(xASLDir, 'External/SPMmodified/');
AtlasDir = fullfile(xASLDir, 'External/SPMmodified/MapsAdded/');  %this folder should have a central WM region map for noise extraction







for iSub = 1:length(Subjects)
    
    Subject = Subjects{iSub};
    sSubject = ['s' Subject];
    
   
    
    SubjDir = fullfile(AnalysisDir, Subject);
    
    
    %% Structural
    if bStructural
        
        % Checking images   %Temporary just for T1
        fprintf('Checking Structural Images')
        T1Path = fullfile(SubjDir, 'T1.nii');
        c1T1Path = fullfile(SubjDir, 'c1T1.nii');
        c2T1Path = fullfile(SubjDir, 'c2T1.nii');
        c3T1Path = fullfile(SubjDir, 'c3T1.nii');
        yfile =  fullfile(SubjDir, 'y_T1.nii');
        
        if ~exist(T1Path) || ~exist(c1T1Path) || ~exist(c2T1Path) || ~exist(c3T1Path)
            warning(['One or more structural image is missing for subject ' Subject ' ...Skipping structural QC, please check your data'])
        else
            
            
            % Noise and Information Theory Domains
            NoisIT = xQC_Noise_IT(T1Path, c1T1Path, c2T1Path, SubjDir, AtlasDir, x, yfile);
            QC.(sSubject).Structural.Noise = NoisIT.Noise;
            QC.(sSubject).Structural.IT = NoisIT.IT;
            
            QC.(sSubject).Structural.Noise.CJV = xQC_CJV(T1Path, c1T1Path, c2T1Path);
            
            %Asymmetry Domain
            QC.(sSubject).Structural.Asymmetry.AI_perc = xQC_AsymmetryIndex(T1Path);
            
            % General Image Quality
            QA_Output = xASL_qc_CAT12_IQR(T1Path,  c1T1Path, c2T1Path, c3T1Path, []);
            QC.(sSubject).Structural.General.IQR = QA_Output.IQR;
            
            % Inhomogneity
             [BI_mean, BI_SD] = xQC_BiasIndex(SubjDir, SPMDir, 'T1');
             QC.(sSubject).Structural.Inhomogeneity.BI_mean = BI_mean;
             QC.(sSubject).Structural.Inhomogeneity.BI_SD = BI_SD;
                         
            % Stats
            QC.(sSubject).Structural.Stats = xQC_Stats(T1Path, 1,  c1T1Path, c2T1Path, c3T1Path);
            
            
            
           %missing for now 
            %ICVs = xASL_qc_CollectQC_Structural (x, iSubject)  % How to integrate this???????
            %xASL_qc_TanimotoCoeff(T1Im, Template) % still need to create templates

            
        end
    end
    
    
    
    %% Functional
    
end  


save(fullfile(AnalysisDir, 'QC.mat'), 'QC')

end 
    
    
