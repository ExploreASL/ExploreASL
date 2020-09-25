function [x] = xASL_spm_GLMcontrast(x, CoVariate, Conweights)
%xASL_spm_GLMcontrast %xASL_spm_GLMcontrast ExploreASL wrapper for SPM GLM - contrast part



%% ------------------------------------------------------------------------------------------------------
%% Run contrast creation
% Remove previous contrast
xASL_adm_DeleteFileList(x.S.SPMdir,'^spm(T|F)_.*\.(nii|nii\.gz)$');
xASL_adm_DeleteFileList(x.S.SPMdir,'^con_.*\.(nii|nii\.gz)$');
% Reload mat
load(x.S.SPMmat2);
save(x.S.SPMmat,'SPM');


%% ------------------------------------------------------------------------------------------------------
%% Get cluster extent (crashes if no inmask voxels)
if      strcmpi(x.S.MultiComparisonCorrType,'FWE') || strcmpi(x.S.MultiComparisonCorrType,'uncorrected')
        x.S.ClusterExtent   = 0;
elseif  strcmpi(x.S.MultiComparisonCorrType,'cluster')
        load(x.S.SPMmat);
        [k,Pc] = CorrClusTh(SPM,x.S.clusterPthr,x.S.uncorrThresh,1:50000);
        x.S.ClusterExtent   = k;
end



%% ------------------------------------------------------------------------------------------------------
%% Get type of multiple comparison correction
if      strcmpi(x.S.MultiComparisonCorrType,'FWE')
        x.S.ThreshType  = 'FWE';
        x.S.PrintTitleStats = ['p=' num2str(x.S.uncorrThresh) ', Bonferroni FWE'];

elseif  strcmpi(x.S.MultiComparisonCorrType,'cluster')
        x.S.ThreshType      = 'none';
        x.S.PrintTitleStats = ['primThr p=' num2str(x.S.clusterPthr) ', cluster FWE p=' num2str(x.S.uncorrThresh) ', clustersize ' num2str(k) ' voxels'];

elseif  strcmpi(x.S.MultiComparisonCorrType,'uncorrected')
        x.S.ThreshType      = 'none';
        x.S.PrintTitleStats = ['p=' num2str(x.S.uncorrThresh) ', unc.'];
end

if ~exist('Conweights','var') % allow inputting contrast

    if      x.S.RegressionCOVAR % regression
            % Then the maps shouldn't be predictors, only the
            % covariates
            Conweights                                            = zeros(1,x.S.nSets);

    else    % otherwise maps are predictors

        if      x.S.nSets==1 % one-sample t-test
                Conweights    = [1];

                % this creates t-test contrast [1 -1] if x.S.nSets==2
                % & ANOVA contrasts when x.S.nSets>2
        else
                Conweights    = diff(eye(x.S.nSets));
                % Conweights    = [-1 1 -1 1];
        end
    end



%     %% ------------------------------------------------------------------------------------------------------
%     %% Add covariates
%     if  exist('CoVariate','var')
%         for iCoV=1:length(CoVariate)
%             Conweights(:,end+1)    = 0;
%         end
%     end
end


% %% Fill/complete conweights for all contrasts
% load(x.S.SPMmat);
% nB                      = size(SPM.xX.X,2);
% Conweights(end+1:nB)    = 0;





%% ------------------------------------------------------------------------------------------------------
%% SPM part
matlabbatch{1}.spm.stats.con.spmmat                                 = { x.S.SPMmat };
matlabbatch{1}.spm.stats.con.delete                                 = 0;


if  x.S.nSets>2 % f-contrast, such as for ANOVA
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights            = Conweights;
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name               = 'Test';
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep            = 'none';
else % t-contrast
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights            = Conweights;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name               = 'Test';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep            = 'none';
end

spm_jobman('run',matlabbatch);

% Save SPMmat
load(x.S.SPMmat);
save(x.S.SPMmat3,'SPM');




end
