function xQC_Master (AnalysisDir, xASLDir, RegExp, bStructural, bFunctional, bDWI, bASL)


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
        sSubject = ['s' Subject];
        
        % Check backup file
        if isstruct(QC) % there was a previous run --> check whther the subject was processed
            if isfield(QC, sSubject) %subject was already QCed
                
                fprintf(['Structural images of ' Subject ' have been already QCed ... skippping  \n'])
                continue
                
            end
        end
        
        % Save back-up file every 50 subj
        if ~(mod(iSub,50))   % save every 50 subjects
            
            save(fullfile(AnalysisDir, 'QC_backup.mat'), 'QC')
        end
        
        SubjDir = fullfile(AnalysisDir, Subject);
        
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
            
            fprintf(['QCeing Structural images for ' Subject])
            
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
            
            % Stats  Stats = QC.(sSubject).Structural.Stats
            Stats = xQC_Stats(T1Path, 1,  c1T1Path, c2T1Path, c3T1Path);
            
            %Handle dimension of struct
            
            filedstats = fields(Stats);
            for iSeg = 1:length(filedstats)
                Seg = filedstats{iSeg};
                
                if isstruct(Stats.(Seg))
                    Statistics = fields(Stats.(Seg));
                    for iStat = 1: length(Statistics)
                        S = Statistics{iStat};
                        
                        ThisStat = [Seg '_' S];
                        QC.(sSubject).Structural.Stats.(ThisStat) = Stats.(Seg).(S);
                    end
                else
                    QC.(sSubject).Structural.Stats.(Seg) = Stats.(Seg);
                    
                end
            end
            
            
            
            %missing for now
            %ICVs = xASL_qc_CollectQC_Structural (x, iSubject)  % How to integrate this???????
            %xASL_qc_TanimotoCoeff(T1Im, Template) % still need to create templates
            
            
        end
    end
    
    
    
    %% Functional
    
end


save(fullfile(AnalysisDir, 'QC.mat'), 'QC')


%% Write CSV for R visualization Tool
[QC_Parameter]= xQC_create_parameter_list(QC);
QC_Parameter = ['Subject' QC_Parameter];
Subjects = fields(QC);


QC_T = array2table(zeros(length(Subjects),(length(QC_Parameter))), 'VariableNames', QC_Parameter);
QC_T.Subject = Subjects;





for iSubject = 1:length(Subjects)
    
    Subj = QC_T.Subject{iSubject};
    
    Modalities=fields(QC.(Subj));
    
    for iModality = 1:length(Modalities)
        
        Mod = Modalities{iModality};
        
        if isfield(QC.(Subj), Mod)
            
            Domains = fields(QC.(Subj).(Mod));
            for iDomain = 1:length(Domains)
                Dom = Domains{iDomain};
                
                if isfield(QC.(Subj).(Mod), Dom)
                    
                    Parameters = fields(QC.(Subj).(Mod).(Dom));
                    
                    for iParameter = 1: length(Parameters)
                        
                        Par = Parameters{iParameter};
                        
                        if isfield(QC.(Subj).(Mod).(Dom), Par)
                            
                            ThisPar =[Mod '_' Dom '_' Par];
                            
                            QC_T.(ThisPar)(iSubject) = QC.(Subj).(Mod).(Dom).(Par);
                            
                            
                        end
                    end
                end
            end
        end
    end
end

writetable(QC_T,fullfile(AnalysisDir, 'QC.csv'));
end


