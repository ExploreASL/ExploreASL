function [x] = xASL_spm_GLMdesign(x, CoVariate)
%xASL_spm_GLMdesign ExploreASL wrapper for SPM GLM - statistical design part


%% -----------------------------------------------------------------
%% Admin

    fprintf('%s\n','Running GLM design')

    if ~isfield(x.S,'SPMdir') && isfield(x.S,'StatsDir')
        x.S.SPMdir  = fullfile(x.S.StatsDir,'SPMdir');
    end

    % Delete all previous runs first
    xASL_adm_DeleteFileList(x.S.SPMdir,'^.*\.mat$');
    xASL_adm_DeleteFileList(x.S.SPMdir,'^*_\d*\.(nii|nii\.gz)$');
    xASL_adm_DeleteFileList(x.S.SPMdir,'^*_\d*_.*\.(nii|nii\.gz)$');

    xASL_adm_DeleteFileList(x.S.SPMdir,'^beta.*\.(nii|nii\.gz)$');
    xASL_adm_DeleteFileList(x.S.SPMdir,'^mask\.(nii|nii\.gz)$');
    xASL_adm_DeleteFileList(x.S.SPMdir,'^ResMS\.(nii|nii\.gz)$');
    xASL_adm_DeleteFileList(x.S.SPMdir,'^RPV\.(nii|nii\.gz)$');

%     if ~isfield(x.S,'VBAmask') This VBAmask should be there, is also needed to decompress
%       the Columns into Images
%         ExplicitMask    = fullfile(x.D.PopDir,'VBA_mask_final.nii');
%
%         if  xASL_exist(ExplicitMask)
%             x.S.VBAmask     = xASL_io_ReadNifti(ExplicitMask);
%         end
%     end


    GlobalNormalization     = 0;
    if  isfield(x.S,'GlobalNormalization')
        if  x.S.GlobalNormalization  == 1
            GlobalNormalization     = 1;
        end
    end
    if  GlobalNormalization==1
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_mean        = 1; % calculate mean
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm       = 2; % normalize each image (proportional)
    else
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit        = 1; % don't calculate mean
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm       = 1; % no global normalization
    end

%         matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm       = 3; % covarying for mean
    % scaling by mean also scales variance, so this is good when variance is proportional of mean
    % mean as covariate "models mean" as covariate, hence doesn't have to be proportional

    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no= 1; % no grand mean scaling
    matlabbatch{1}.spm.stats.factorial_design.dir                   = {x.S.SPMdir}; % output dir for design matrix
    matlabbatch{1}.spm.stats.factorial_design.multi_cov             = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none    = 1; % threshold masking
    matlabbatch{1}.spm.stats.factorial_design.masking.im            = 0; % implicit masking doesn't always work, can crash
    xASL_io_ReadNifti(x.S.MaskPath); % unzipping if needed
    matlabbatch{1}.spm.stats.factorial_design.masking.em            = { x.S.MaskPath }; % explicit masking always works




    %% -----------------------------------------------------------------
    %% Load scans

    for iSet=1:x.S.nSets
        scans{iSet} = sort(xASL_adm_GetFileList(x.S.SPMdir, ['^Set' num2str(iSet) 'subject\d*\.nii\.gz$']));
        scans{iSet} = cellfun(@(y) [y '_1'],scans{iSet},'UniformOutput',false);
    end
    if x.S.nSets>1 && x.S.ONE_TWO_SAMPLE_TEST==1 % Check whether the compared datasets have equal size (only required for paired sample tests)
        for iSet1=1:x.S.nSets
            for OtherSets2=iSet1+1:x.S.nSets
                if  length(scans{iSet1})~=length(scans{OtherSets2})
                    error('Datasets do not have equal number of scans!');
                end
            end
        end
    end

    GeneralSPMSettings.variance   = 0; % unequal variance -> ASSUMING EQUAL VARIANCE!!!!!!!
    GeneralSPMSettings.gmsca      = 0; % no additional grand mean scaling
    GeneralSPMSettings.ancova     = 0; % allows different subjects to have different relationships between local & global measures (not for 2nd level)





    %% -----------------------------------------------------------------
    %% Design options

    if      x.S.nSets==1 % 1-sample t-test
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans      =    scans{1};
            x.S.printTitleORI          = [x.S.OutputID ' 1-t-test (n=' num2str(x.S.nScans{1}) ')'];

    elseif  x.S.nSets==2 && x.S.ONE_TWO_SAMPLE_TEST==1 % paired 2-sample t-test (dependent samples)

            matlabbatch{1}.spm.stats.factorial_design.des.pt            = GeneralSPMSettings;
            for iSubject=1:length(scans{1})
                matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(iSubject).scans = {scans{1}{iSubject};scans{2}{iSubject}};
            end

            x.S.printTitleORI          = [x.S.OutputID ' 1-t-test ' x.S.SetsName{x.S.iCurrentSet} ' ' x.S.SetsOptions{x.S.iCurrentSet}{x.S.iSet{2}} ' (n=' num2str(x.S.nScans{2}) ') - ' x.S.SetsOptions{x.S.iCurrentSet}{x.S.iSet{1}} ' (n=' num2str(x.S.nScans{1}) ')'];

    elseif  x.S.nSets==2 && x.S.ONE_TWO_SAMPLE_TEST==2 % unpaired 2-sample t-test (independent samples)

            matlabbatch{1}.spm.stats.factorial_design.des.t2            = GeneralSPMSettings;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1     = scans{1};
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2     = scans{2};
            matlabbatch{1}.spm.stats.factorial_design.des.t2.dept       = 0; % independence

            x.S.printTitleORI          = [x.S.OutputID ' 2-t-test ' x.S.SetsName{x.S.iCurrentSet} ' ' x.S.SetsOptions{x.S.iCurrentSet}{x.S.iSet{2}} ' (n=' num2str(x.S.nScans{2}) ') - ' x.S.SetsOptions{x.S.iCurrentSet}{x.S.iSet{1}} ' (n=' num2str(x.S.nScans{1}) ')'];

    elseif  x.S.nSets>2 && x.S.ONE_TWO_SAMPLE_TEST==1 % within-subject ANOVA (dependent samples)

            matlabbatch{1}.spm.stats.factorial_design.des.anovaw        = GeneralSPMSettings;
            for iScan=1:length(scans)     % for each set
                for iSubject=1:length(scans{1})
                    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(iSubject).scans{iScan,1}  = scans{iScan}{iSubject};
                    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(iSubject).conds           = [1:1:length(scans)];
                end
            end

            matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept       = 1; % dependence
            x.S.printTitleORI          = [x.S.OutputID ' 1w ws-ANOVA ' x.S.SetsName{x.S.iCurrentSet}];


    elseif  x.S.nSets>2 && x.S.ONE_TWO_SAMPLE_TEST==2 % ANOVA (independent samples)

            matlabbatch{1}.spm.stats.factorial_design.des.anova             = GeneralSPMSettings;

            for iSet=1:x.S.nSets
                matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(iSet).scans = scans{iSet};
            end

            matlabbatch{1}.spm.stats.factorial_design.des.anova.dept        = 0; % independence
            x.S.printTitleORI          = [x.S.OutputID ' 1w ANOVA ' x.S.SetsName{x.S.iCurrentSet}];
    end





    %% -----------------------------------------------------------------
    %% Put in co-variate(s), this should be permuted in level above this function
    if ~exist('CoVariate','var') % no co-variate
        matlabbatch{1}.spm.stats.factorial_design.cov                          = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    else
        for iC=1:length(CoVariate)
            matlabbatch{1}.spm.stats.factorial_design.cov(iC).c                 = double(CoVariate(iC).Data);
            matlabbatch{1}.spm.stats.factorial_design.cov(iC).cname             = CoVariate(iC).Name;
            matlabbatch{1}.spm.stats.factorial_design.cov(iC).iCFI              = CoVariate(iC).iCFI; % 3; % interactions (1=with factor 1 (==none), 2 with factor 2 (which is 1), etc)
            matlabbatch{1}.spm.stats.factorial_design.cov(iC).iCC               = 1; % mean-centering (overall mean, default) 5 = no mean-centering
        end
    end

    spm_jobman('run',matlabbatch);
    % Save SPM1
    load(x.S.SPMmat);
    save(x.S.SPMmat1,'SPM');




    %% -----------------------------------------------------------------
    %% Model review
    if  x.S.PrintSPMOutput
        try
            clear matlabbatch
            matlabbatch{1}.spm.stats.review.spmmat              = { x.S.SPMmat };
            matlabbatch{1}.spm.stats.review.display.matrix      = 1;
            matlabbatch{1}.spm.stats.review.print               = 'pdf';
            spm_jobman('run',matlabbatch);
            close all
        catch
        end
    else
        % Skip review
    end


end
