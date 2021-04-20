function xQC_Master

ExploreQC = mfilename('fullpath');
[fdir, ~, ~ ]= fileparts(ExploreQC);


%Administration : Read configuration file and extract parameters from there
configfile = fullfile(fdir,'ConfigFile.xlsx');
config = readtable(configfile);
domainconfig = readtable(configfile, 'Sheet', 'DomainsSheet', 'Range' , 'A:B');

% modules to run
bStructural = config.PE_properties{strcmp(config.ParameterExtractionModule, 'Structural')};
bFunctional = config.PE_properties{strcmp(config.ParameterExtractionModule, 'Functional')};
bDiffusion  = config.PE_properties{strcmp(config.ParameterExtractionModule, 'Diffusion')};
bASL        = config.PE_properties{strcmp(config.ParameterExtractionModule, 'ASL')};

% Folders and regular expressions
AnalysisDir = config.PE_properties{strcmp(config.ParameterExtractionModule, 'AnalysisDir')};
xASLDir     = config.PE_properties{strcmp(config.ParameterExtractionModule, 'xASLDir')};
SubjRegExp  = config.PE_properties{strcmp(config.ParameterExtractionModule, 'SubjRegExp')};
SiteRegExp  = config.PE_properties{strcmp(config.ParameterExtractionModule, 'SiteRegExp')};


%handle REGexp

% First initialize exploreASL to obtain x, because some functions still
% need the x structure to work (deformations)
% Also add the paths to external tools(dti denoising)

cd (xASLDir);
x = ExploreASL_Master('', 0);
cd(fdir);

% TEMPORARY : Add the path to scripts in development to test them
addpath(fullfile(fdir, 'xQC_Par_Extr_Tool', 'Development/mppca_denoise-master/mppca_denoise-master'))
addpath(fullfile(fdir, 'xQC_Par_Extr_Tool', 'Development/'))
addpath (fullfile(fdir, 'xQC_Par_Extr_Tool'))

%Get the list of subject with subjRegExpl
Subjects = xASL_adm_GetFsList(AnalysisDir, SubjRegExp, 1);

%Preliminary settings, can be change by the user
SPMDir = fullfile(xASLDir, 'External/SPMmodified/');
AtlasDir = fullfile(xASLDir, 'External/SPMmodified/MapsAdded/');  %this folder should have a central WM region map for noise extraction



%Load Back-up file if exist

backup = fullfile(AnalysisDir, 'QC_backup.mat');

if exist(backup)
    
    fprintf('Previous QC analysis found, starting from there. If you wish to rerun the full analysis please delete QC_backup.mat file \n')
    load(backup);
    % delete last subject to be sure
    allsub = fields(QC);
    QC = rmfield(QC, allsub(length(allsub)));
else
    QC={};
end





%% Structural

if bStructural
    fprintf('Running Structural QC  \n')
    
    
    % Iterate across subjects to run structural QC
    for iSub = 1:length(Subjects)
        
        Subject = Subjects{iSub};
        sSubject = ['s' Subject]; % needed to create the matlab struct
        
        % Check backup file
        if isstruct(QC) % there was a previous run, check wheter subject has been processed
            if isfield(QC, sSubject) % Some QC has been done on the subject --> check whether structural exist
                if isfield (QC.(sSubject), 'Structural')
                    fprintf(['Structural images of ' Subject ' have already been  QCed ... skippping  \n'])
                    continue
                end
            end
        end
        
        % Save back-up file every 50 subj
        if ~(mod(iSub,50))   % save every 50 subjects
            
            save(fullfile(AnalysisDir, 'QC_backup.mat'), 'QC')
        end
        
        SubjDir = fullfile(AnalysisDir, Subject);
        %if folder is not specified, folder is subject folder
        if isempty(config.PE_properties{strcmp(config.ParameterExtractionModule, 'StructFold')})
            StructFold =SubjDir;
        else
            StructFold = fullfile(SubjDir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'StructFold')});
        end
        
        % Checking images   
        fprintf('Checking Structural Images \n')
        
        T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'StructIm')});
        c1T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'GMim')});
        c2T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'WMim')});
        c3T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'CSFim')});
        yfile =  fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'Struct_Y_file')});
        
        if ~exist(T1Path) || ~exist(c1T1Path) || ~exist(c2T1Path) || ~exist(c3T1Path)
            fprintf(['One or more structural image is missing for subject ' Subject ' ...Skipping structural QC, please check your data'])
            
            QC = xQC_missing(QC, sSubject, 'Structural', configfile);
            continue
        end
        
        fprintf(['QCeing Structural images for ' Subject '\n'])
        
        % Noise Domain
        NoisIT = xQC_Noise_IT(T1Path, c1T1Path, c2T1Path, StructFold, AtlasDir, x, yfile);
        
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Structural_Noise'),2}
            QC.(sSubject).Structural.Noise = NoisIT.Noise;
            QC.(sSubject).Structural.Noise.CJV = xQC_CJV(T1Path, c1T1Path, c2T1Path);
        end
        
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Structural_Motion'),2} %motion in T1 is extracte with Information theory measures
            QC.(sSubject).Structural.Motion = NoisIT.IT;
        end
        
        
        %Asymmetry Domain
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Structural_Asymmetry'),2}
            QC.(sSubject).Structural.Asymmetry.AI_perc = xQC_AsymmetryIndex(T1Path);
        end
        
       
        
        % Inhomogeneity Domain
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Structural_Inhomogeneity'),2}
            [BI_mean, BI_SD] = xQC_BiasIndex(StructFold, SPMDir, 'T1');
            QC.(sSubject).Structural.Inhomogeneity.BI_mean = BI_mean;
            QC.(sSubject).Structural.Inhomogeneity.BI_SD = BI_SD;
            QA_Output = xASL_qc_CAT12_IQR(T1Path,  c1T1Path, c2T1Path, c3T1Path, []);
            QC.(sSubject).Structural.Inhomogeneity.IQR = QA_Output.IQR;
        end
        
        % Descriptives Domain
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Structural_Descriptives'),2}
            Descriptives_str = xQC_Descriptives(T1Path, 1,  c1T1Path, c2T1Path, c3T1Path);
            
            %Handle dimension of struct
            
            filedDescriptives = fields(Descriptives_str);
            for iSeg = 1:length(filedDescriptives)
                Seg = filedDescriptives{iSeg};
                
                if isstruct(Descriptives_str.(Seg))
                    Statistics = fields(Descriptives_str.(Seg));
                    for iStat = 1: length(Statistics)
                        S = Statistics{iStat};
                        
                        ThisStat = [Seg '_' S];
                        QC.(sSubject).Structural.Descriptives.(ThisStat) = Descriptives_str.(Seg).(S);
                    end
                else
                    QC.(sSubject).Structural.Descriptives.(Seg) = Descriptives_str.(Seg);
                    
                end
            end
            
        end
        
        
        
    end
end



%% Functional

if bFunctional
    
    % Iterate across subjects to run functional QC
    for iSub = 1:length(Subjects)
        
        Subject = Subjects{iSub};
        sSubject = ['s' Subject];
        
        % Check backup file
        if isstruct(QC) % there was a previous run --> check whether the subject was processed
            if isfield(QC, sSubject) % Some QC was already done on the subject --> check whether functional
                if isfield(QC.(sSubject), 'Functional')
                    fprintf(['Functional images of ' Subject ' have already been  QCed ... skippping  \n'])
                    continue
                end
            end
        end
        
        % Save back-up file every 50 subj
        if ~(mod(iSub,50))   % save every 50 subjects
            
            save(fullfile(AnalysisDir, 'QC_backup.mat'), 'QC')
        end
        
        % Setting Paths and Folders
        fprintf('Checking Functional Images \n')
        SubjDir = fullfile(AnalysisDir, Subject);
        
        % Handle functional and structural folders
        %Functional: if folder is not specified, folder is subject folder
        if isempty(config.PE_properties{strcmp(config.ParameterExtractionModule, 'FuncFold')})
            FuncFold =SubjDir;
        else
            FuncFold = fullfile(SubjDir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'FuncFold')});
        end
        %Structural if folder is not specified, folder is subject folder
        if isempty(config.PE_properties{strcmp(config.ParameterExtractionModule, 'StructFold')})
            StructFold =SubjDir;
        else
            StructFold = fullfile(SubjDir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'StructFold')});
        end
        
        
        
        
        
        boldnii = fullfile(FuncFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'FuncIm')});
        GhostTemplate = fullfile(AtlasDir, 'GhostSignalRatio.nii');
        c1T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'GMim')});
        c2T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'WMim')});
        c3T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'CSFim')});
        yfile =  fullfile(FuncFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'Func_Y_file')});
        % Check if file exist
        if ~exist(boldnii, 'file') || ~exist(c1T1Path , 'file') || ~exist(c2T1Path, 'file') || ~exist(c3T1Path, 'file') || ~exist(yfile, 'file')
            fprintf(['One or more necessary images is missing for subject ' Subject ' ...Skipping Functional QC, please check your data'])
            
            QC = xQC_missing(QC, sSubject, 'Functional', configfile);
            
            continue
        end
        
        fprintf(['QCeing Functional Images for Subject' Subject '\n'])
        
        nii4D     = xASL_io_Nifti2Im(boldnii);
        
        if size(nii4D,4)<6
            fprintf(['Functional Image of subject ' Subject ' has less than 6 volumes ...skipping'])
            continue
        end
        
        %Motion Domain
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Functional_Motion'),2}
            NoisIT = xQC_Noise_IT(boldnii, c1T1Path, c2T1Path, StructFold, AtlasDir, x, yfile);
            
            QC.(sSubject).Functional.Motion = NoisIT.IT;
            
            QC.(sSubject).Functional.Motion.GhostToSignal = xQC_GhostToSignal(boldnii, GhostTemplate);
       
            QC.(sSubject).Functional.Motion.GlobalCorrelation = xQC_Global_correlation(boldnii, c1T1Path, StructFold);
        
            if isempty(config.PE_properties{strcmp(config.ParameterExtractionModule, 'PathToMotion')})
                PathToMotion =FuncFold;
            else
                PathToMotion = fullfile(config.PE_properties{strcmp(config.ParameterExtractionModule, 'PathToMotion')});
            end
            
            Motion = xQC_Motion_Functional(PathToMotion, Subject);
            
            % Take the second value output of xQC_Motion_Functional
            QC.(sSubject).Functional.Motion.MotionMean_mm = Motion.MotionMean_mm(2);
            QC.(sSubject).Functional.Motion.MotionSD_mm = Motion.MotionSD_mm(2);
            QC.(sSubject).Functional.Motion.MotionMax_mm = Motion.MotionMax_mm(2);
            QC.(sSubject).Functional.Motion.MotionExcl_Perc = Motion.MotionExcl_Perc;
        
        end
        
        
        % Noise Domain
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Functional_Noise'),2}
            %NoisIT = XQC_Noise_IT(boldnii, c1T1Path, c2T1Path, SubjDir, AtlasDir, x, yfile);
            %QC.(sSubject).Functional.Temporal.tSNR_GM_Ratio = NoisIT.Noise.SNR_GM_Ratio;
            %QC.(sSubject).Functional.Temporal.CNR_GM_WM_ratio = NoisIT.Noise.CNR_GM_WM_Ratio;
            %For the moment we use xASL tSNRs which are more, discuss this
            %and decide whether xQC_Noise_IT is okay for functional
            %Reslice GM and WM to Nifti
            xASL_spm_reslice(boldnii,c1T1Path , [], [], [], fullfile(StructFold, 'tmp_GM_mask.nii'), 0)
            xASL_spm_reslice(boldnii,c2T1Path , [], [], [], fullfile(StructFold, 'tmp_WM_mask.nii'), 0)
            
            tSNR = xASL_qc_temporalSNR(boldnii,{fullfile(StructFold, 'tmp_GM_mask.nii') fullfile(StructFold, 'tmp_WM_mask.nii')});
            
            xASL_delete(fullfile(StructFold, 'tmp_GM_mask.nii'))
            xASL_delete(fullfile(StructFold, 'tmp_WM_mask.nii'))
            
            QC.(sSubject).Functional.Noise = tSNR;
           
        end
        
                
        % SUMMARY STATISTIC
        % Descriptives  Descriptives = QC.(sSubject).Structural.Descriptives
        if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Functional_Descriptives'),2}
            Descriptives_fun = xQC_Descriptives(boldnii, 0,  c1T1Path, c2T1Path, c3T1Path);
            
            
            %Handle dimension of Descriptives
            
            filedDescriptives = fields(Descriptives_fun);
            for iSeg = 1:length(filedDescriptives)
                Seg = filedDescriptives{iSeg};
                
                if isstruct(Descriptives_fun.(Seg))
                    Statistics = fields(Descriptives_fun.(Seg));
                    for iStat = 1: length(Statistics)
                        S = Statistics{iStat};
                        
                        ThisStat = [Seg '_' S];
                        QC.(sSubject).Functional.Descriptives.(ThisStat) = Descriptives_fun.(Seg).(S);
                    end
                else
                    QC.(sSubject).Functional.Descriptives.(Seg) = Descriptives_fun.(Seg);
                    
                end
            end
            
            
        end
        
        
        
        
    end
    
    
    
    
    
    %% DWI QC
    if bDiffusion
        
        for iSub = 1:length(Subjects)
            
            Subject = Subjects{iSub};
            sSubject = ['s' Subject];
            
            % Check backup file
            if isstruct(QC) % there was a previous run --> check whether the subject was processed
                if isfield(QC, sSubject) % Some QC was already done on the subject --> check whether diffusion
                    if isfield(QC.(sSubject), 'Diffusion')
                        fprintf(['Diffusion images of ' Subject ' have already been  QCed ... skippping  \n'])
                        continue
                    end
                end
            end
            
            % Save back-up file every 50 subj
            if ~(mod(iSub,50))   % save every 50 subjects
                
                save(fullfile(AnalysisDir, 'QC_backup.mat'), 'QC')
            end
            
            fprintf('Checking Diffusion Images \n')
            % Setting Paths and Folders
            SubjDir = fullfile(AnalysisDir, Subject);
            
            % Handle Structural and Diffusion Folder
            %if folder is not specified, folder is subject folder
            if isempty(config.PE_properties{strcmp(config.ParameterExtractionModule, 'DiffFold')})
                DWIdir =SubjDir;
            else
                DWIdir = fullfile(SubjDir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'DiffFold')});
            end
            if isempty(config.PE_properties{strcmp(config.ParameterExtractionModule, 'StructFold')})
                StructFold =SubjDir;
            else
                StructFold = fullfile(SubjDir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'StructFold')});
            end
            
            
            c1T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'GMim')});
            c2T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'WMim')});
            c3T1Path = fullfile(StructFold, config.PE_properties{strcmp(config.ParameterExtractionModule, 'CSFim')});
            FApath = fullfile(DWIdir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'FAIm')});
            ADCpath = fullfile(DWIdir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'ADCIm')});
            TopUp_Path = fullfile(DWIdir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'dwi_TopUp_im')});
            
            %Check for the existence of these files that are necessary for the
            %analysis
            
            if ~exist(FApath, 'file') ||~exist(ADCpath, 'file') || ~exist(DWIdir, 'dir') || ~exist(c1T1Path, 'file') || ~exist(c2T1Path, 'file') || ~exist(c3T1Path, 'file') || ~exist(TopUp_Path, 'file')
                
                fprintf(['One or more necessary files not found for subject ' Subject 'skipping' ])
                QC = xQC_missing(QC, sSubject, 'Diffusion', configfile);
                
                continue
            end
            
            fprintf(['QCeing Diffusion Images for subject ' Subject '\n'])
            % Motion Domain
            if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Diffusion_Motion'),2}
                [EddyQC] = xQC_dwi_EddyMotion(DWIdir, Subject, c1T1Path, c2T1Path, c3T1Path, TopUp_Path);
                
                QC.(sSubject).Diffusion.Motion = EddyQC.motion;
                QC.(sSubject).Diffusion.Motion.Outlier_Perc = EddyQC.Outlier_Perc;
                QC.(sSubject).Diffusion.Motion.Induced_Distortion = EddyQC.Induced_Distortion;
                QC.(sSubject).Diffusion.Motion.Susc_Induced_Distortion = EddyQC.Susc_Induced_Distortion;
                
            end
            
            % Noise Domain
         
            
           
            if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Diffusion_Noise'),2}
                [Noise_mean , Noise_SD ] = xQC_dwi_Noise(StructFold, 'dwi_run-1_dwi.nii', 0 , 0) ; % don't save the noise image for now
                
                QC.(sSubject).Diffusion.Noise.Mean_Noise_dwi = Noise_mean;
                
                QC.(sSubject).Diffusion.Noise.SD_Noise_dwi = Noise_SD;
                if ~isempty(config.PE_properties{strcmp(config.ParameterExtractionModule, 'SSE_im')})
                    SSEpath = fullfile(DWIdir, config.PE_properties{strcmp(config.ParameterExtractionModule, 'SSE_im')});
                    SSE_WM_Mean = xQC_SSE_WM_mean(SSEpath, c2T1Path);            
                    QC.(sSubject).Diffusion.Noise.SSE_WM_Mean = SSE_WM_Mean;
                    
                end            
            
            
            end
            
            % DERIVATIVES % Descriptives %
            if domainconfig{strcmp(domainconfig.ScanType_Domain, 'Diffusion_Descriptives'),2}
                QC.(sSubject).Diffusion.Descriptives.FA_Outliers_mL = xASL_qc_FA_Outliers(FApath);
                
                % FA Descriptives
                FA_Descriptives = xQC_Descriptives(FApath, 0, c1T1Path, c2T1Path, c3T1Path);
                
                filedDescriptives = fields(FA_Descriptives);
                for iSeg = 1:length(filedDescriptives)
                    Seg = filedDescriptives{iSeg};
                    
                    if isstruct(FA_Descriptives.(Seg))
                        Statistics = fields(FA_Descriptives.(Seg));
                        for iStat = 1: length(Statistics)
                            S = Statistics{iStat};
                            
                            ThisStat = ['FA_' Seg '_' S];
                            QC.(sSubject).Diffusion.Descriptives.(ThisStat) = FA_Descriptives.(Seg).(S);
                        end
                    else
                        QC.(sSubject).Diffusion.Descriptives.(Seg) = FA_Descriptives.(Seg);
                        
                    end
                end
                
                % ADC Descriptives
                ADC_Descriptives = xQC_Descriptives(ADCpath, 0, c1T1Path, c2T1Path, c3T1Path);
                
                filedDescriptives = fields(ADC_Descriptives);
                for iSeg = 1:length(filedDescriptives)
                    Seg = filedDescriptives{iSeg};
                    
                    if isstruct(ADC_Descriptives.(Seg))
                        Statistics = fields(ADC_Descriptives.(Seg));
                        for iStat = 1: length(Statistics)
                            S = Statistics{iStat};
                            
                            ThisStat = ['ADC_' Seg '_' S];
                            QC.(sSubject).Diffusion.Descriptives.(ThisStat) = ADC_Descriptives.(Seg).(S);
                        end
                    else
                        QC.(sSubject).Diffusion.Descriptives.(Seg) = ADC_Descriptives.(Seg);
                        
                    end
                end
                
                
                
            end
            
            
            
            
            
            
        end
        
        
    end
    
    
    
    
    
    
    %% Save File
    
    save(fullfile(AnalysisDir, 'QC.mat'), 'QC')
    save(fullfile(AnalysisDir, 'QC_backup.mat'), 'QC')
    
    %% Write CSV for R visualization Tool
    [list , QC_Parameters] = xQC_create_parameter_Table(configfile);
    QC_Parameters = ['Subject' 'Site' QC_Parameters];
    Subjects = fields(QC);
   
    
    QC_T = array2table(zeros(length(Subjects),(length(QC_Parameters))), 'VariableNames', QC_Parameters);
    QC_T.Subject = Subjects;
    

    for iSubject = 1:length(Subjects)
        Subj = QC_T.Subject{iSubject};
        subjori = Subj(2:end); % Take the original name of the subject, we added an s for the struct
        QC_T.Site(iSubject)= string(regexp(subjori, SiteRegExp, 'Match'));  % Extract Site ith regular expression from configuration
       
        for iPar = 1:height(list)
            if list.Visualize(iPar) %just if include
                st= list.Scantype{iPar};
                dom= list.Domain{iPar};
                param = list.Parameter{iPar};
                fullparam= list.ParameterName{iPar};
                if isfield(QC.(Subj), st)
                    if isfield(QC.(Subj).(st), dom)
                        if isfield(QC.(Subj).(st).(dom), param)
                            QC_T.(fullparam)(iSubject) = QC.(Subj).(st).(dom).(param);
                        end
                    end
                end
            end
        end
        QC_T.Subject(iSubject) = {subjori};
    end
    
       
       
        
        
        
        
%     for iSubject = 1:length(Subjects)
%         
%         Subj = QC_T.Subject{iSubject};
%         
%         Modalities=fields(QC.(Subj));
%         
%         for iModality = 1:length(Modalities)
%             
%             Mod = Modalities{iModality};
%             
%             if isfield(QC.(Subj), Mod)
%                 
%                 Domains = fields(QC.(Subj).(Mod));
%                 for iDomain = 1:length(Domains)
%                     Dom = Domains{iDomain};
%                     
%                     if isfield(QC.(Subj).(Mod), Dom)
%                         
%                         Parameters = fields(QC.(Subj).(Mod).(Dom));
%                         
%                         for iParameter = 1: length(Parameters)
%                             
%                             Par = Parameters{iParameter};
%                             
%                             if isfield(QC.(Subj).(Mod).(Dom), Par)
%                                 
%                                 ThisPar =[Mod '_' Dom '_' Par];
%                                 
%                                 QC_T.(ThisPar)(iSubject) = QC.(Subj).(Mod).(Dom).(Par);
%                                 
%                                 
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
    writetable(QC_T,fullfile(AnalysisDir, 'QC.csv'));
    writetable(QC_T,fullfile(fdir, 'xQC_Visualization_Tool', 'dataframes','QC.csv' )) % save it for the visualization
    
    
    
end


