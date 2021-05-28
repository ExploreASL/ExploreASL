function [result, x] = xASL_module_DARTEL(x)
%xASL_module_DARTEL ExploreASL module for non-linear registration
%
% FORMAT: [result, x] = xASL_module_DARTEL(x)
%
% INPUT:
%   x  - x structure containing all input parameters (REQUIRED)
%   x.dir.SUBJECTDIR  -  anatomical directory, containing the derivatives of anatomical images (REQUIRED)
%
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters (REQUIRED)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This module performs a non-linear registration between subjects, using first time points only (saving the rest for longitudinal registration).
% While this is now included in CAT12, using existing templates, it can still be useful to recreate population-specific templates that don't fit the
% existing templates (e.g. for studies with children).
%
% SPM12 DARTEL                             : DARTEL providing full non-linear freedom, could potentially overfit, but more robust to lesions
%                                            (when properly lesion-masked) or extracranial growth.
% This module will run only for the first T1w/time point of each subject, if there are multiple time points
% Submodules:
% 1) DARTEL Estimate = run DARTEL to create new templates
% 2) Rerun DARTEL for newly found subjects (e.g. if they are <5%, otherwise we should simply rerun submodule 1).
%    This runs DARTEL for existing templates, using the templates created in submodule 1
% 3) Repetition of resampling of images to standard space & visual QC (as performed as submodules 9 & 12 in the structural module)
%
% EXAMPLE: [result, x] = xASL_module_DARTEL(x)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
% Ashburner J. A fast diffeomorphic image registration algorithm. Neuroimage. 2007;38(1):95
% __________________________________
% Copyright ? 2015-2019 ExploreASL

% PM: Replace this by Geodesic Shooting



%% -------------------------------------------------------------------------------------------------------
%% Administration
x = xASL_init_InitializeMutex(x, 'DARTEL'); % starts mutex locking process to ensure that everything will run only once
x.IMAGES_CREATE     = {['rc1' x.P.STRUCT] ['rc2' x.P.STRUCT]};
x.ESTIMATE          = true; % 1 = estimate DARTEL flowfields, if not already processed. 0 = skip estimating
x.WARP              = true; % 1 = apply DARTEL, creating warped maps; 0 = skip warping
x.forceAllSubjects  = true;
x.DARTEL_TEMPLATE   = 'Template';

bFirstRun = 0;

if ~x.mutex.HasState('999_ready')
    %% -------------------------------------------------------------------------------------------------------
    %% 1    DARTEL estimate
    if ~x.mutex.HasState('001_create_templates') && x.ESTIMATE
		bFirstRun = 1;

        [warp, TotalnEstimate] = GetWarpList(x);
        warp.settings.template = x.DARTEL_TEMPLATE;

        fprintf(['\n' xASL_num2str(length(x.IMAGES_CREATE) ) ' image types will be warp estimated by DARTEL,']);
        fprintf('%s\n',[xASL_num2str(TotalnEstimate) ' files in total']);

        fprintf('\n\n');

        if  x.Quality==1
            slam = [ 8 4 2 1 0.5 0.25];
            % with CAT12 segmentation, we can achieve less blurry end-results
            % and need a less blurry start
        else
            slam = [16 8 4 2   1  0.5];
        end

        % Run the DARTEL steps 1 by 1, or skip when partly ran previously
        for iIT=1:length(slam)
            if ~x.mutex.HasState(['001_create_templates_' num2str(iIT)])

                warp                                 = GeneralDARTELsettings(x,warp,iIT);
                warp.settings.param(1).slam          = slam(iIT);
                matlabbatch{1}.spm.tools.dartel.warp = warp;

                spm_jobman('run',matlabbatch);

                % Copy created templates, to avoid overwriting by iteration
                tempTemplate1 = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_0.nii']);
                tempTemplate2 = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_' num2str(iIT-1) '_temp.nii']);
				if exist(fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_0.nii.gz']),'file')
					delete(fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_0.nii.gz']));
				end
                xASL_Copy(tempTemplate1, tempTemplate2, true);
                tempTemplate1 = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_1.nii']);
				if exist(fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_1.nii.gz']),'file')
					delete(fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_1.nii.gz']));
				end
                tempTemplate2 = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_' num2str(iIT  ) '_temp.nii']);
                xASL_Copy(tempTemplate1, tempTemplate2, true);

                x.mutex.AddState(['001_create_templates_' num2str(iIT)]);
            end
        end

        for iIT=0:length(slam)
            tempTemplate1 = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_' num2str(iIT) '_temp.nii']);
            tempTemplate2 = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_' num2str(iIT) '.nii']);
            xASL_Move(tempTemplate1, tempTemplate2, true); % force overwriting!
        end

            x.mutex.AddState('001_create_templates');
            fprintf('%s\n','001_create_templates was performed');
    else
            fprintf('%s\n','001_create_templates           has already been performed, skipping...');
    end

end



%% -------------------------------------------------------------------------------------------------------
%% 2 ReRun DARTEL if additional T1 scans are found
% This should now be re-pogrammed, on the basis of lock files
FList_FlowField = xASL_adm_GetFileList(x.D.PopDir, ['^u_rc1T1_.*_Template\.nii$']);

% Create list
ToWarpList = '';

for iS=1:x.nSubjects

    % Here, check longitudinal registration.
    % If we run longitudinal registration, only use first volumes
    x.P.SubjectID = x.SUBJECTS{iS};
    [~, ~, IsSubject, SubjectID_FirstVolume] = xASL_init_LongitudinalRegistration(x);

    % To check whether or not we will run longitudinal registration
	
    if  IsSubject==1 % only run DARTEL for first volumes
        ToWarpList{end+1}    = SubjectID_FirstVolume;
    end
end


if  ~isempty(FList_FlowField) && ~bFirstRun
    % if there new subjects at all but not the first run
    if  0.05*length(ToWarpList)<length(FList_FlowField)
        % if there are more than 5% new subjects
        error('More than 5% new subjects, have to rerun total T1 DARTEL module!');

    else    % create missing flowfields using existing T1 DARTEL population templates


        %% First resurrect previous DARTEL templates if they are missing

        for ii=1:7
            if xASL_exist( fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_' num2str(ii-1) '.nii']), 'file')
                ExistTemplate(ii) = true;
            else ExistTemplate(ii) = f;
            end
        end

        if  sum(ExistTemplate)<7

            StartSmooth     = 1;
            EndSmooth       = 8;
            StepsN          = 4;

            OUTPUTim        = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_0.nii']);
            OUTPUTim2       = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_5.nii']);
            xASL_Move( OUTPUTim ,OUTPUTim2,1);
            OUTPUTim        = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_1.nii']);
            OUTPUTim2       = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_6.nii']);
            xASL_Move( OUTPUTim ,OUTPUTim2,1);

            for ii=1:5

                fwhm(ii,:)      = repmat(EndSmooth- (ii-1)*((EndSmooth-StartSmooth)/StepsN),[1 3]);

                clear matlabbatch
                INPUTim{1,1} = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_5.nii,1']);
                INPUTim{2,1} = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_5.nii,2']);

                matlabbatch{1}.spm.spatial.smooth.data      = INPUTim;
                matlabbatch{1}.spm.spatial.smooth.fwhm      = fwhm(ii,:);
                matlabbatch{1}.spm.spatial.smooth.dtype     = 0;
                matlabbatch{1}.spm.spatial.smooth.im        = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix    = 's';

                spm_jobman('run',matlabbatch);

                OUTPUTim = fullfile(x.D.PopDir,['s' x.DARTEL_TEMPLATE '_5.nii']);
                OUTPUTim2 = fullfile(x.D.PopDir,[x.DARTEL_TEMPLATE '_' num2str(ii-1) '.nii']);

                xASL_Move(OUTPUTim , OUTPUTim2, true);
            end
        end

        %% Create list of subjects with missing FlowFields

        MisSubList  = '';

        for iSub=1:length(ToWarpList)
            Fname = fullfile( x.D.PopDir, ['u_rc1' x.P.STRUCT '_' ToWarpList{iSub} '_' x.DARTEL_TEMPLATE '.nii']);
            if ~xASL_exist(Fname, 'file')
                MisSubList{end+1,1} = ToWarpList{iSub};
            end
        end

        warp = GetWarpList(x, MisSubList);

        % Run the DARTEL steps 1 by 1, or skip when partly ran previously
        for iIT=1:6 %  depends on settings!
            warp    = GeneralDARTELsettings(x,warp,iIT);

            clear matlabbatch

            % For existing templates
            warp.settings.param(1).template = {fullfile( x.D.PopDir, [x.DARTEL_TEMPLATE '_' num2str(iIT) '.nii'])} ;
            xASL_io_ReadNifti(warp.settings.param(1).template);
            % template_0 is first average/template, template_1 is first created average/template
            % likewise, template_5 is last average/template BEFORE rerunning DARTEL outer iteration, template_6 is final created average/template
            matlabbatch{1}.spm.tools.dartel.warp1 = warp;

            spm_jobman('run',matlabbatch);
        end

        % Rename newly created files
        for iSub=1:length(MisSubList)
            tempTemplate1 = fullfile(x.D.PopDir,['u_rc1' x.P.STRUCT '_' MisSubList{iSub} '.nii']);
            tempTemplate2 = fullfile(x.D.PopDir,['u_rc1' x.P.STRUCT '_' MisSubList{iSub} '_T1_template.nii']);
            xASL_Move(tempTemplate1, tempTemplate2, true);
        end

        % If there are old and new files for the pre-existing ones, keep
        % new ones
        for iS=1:x.nSubjects
            tempTemplate1 = fullfile(x.D.PopDir,['u_rc1' x.P.STRUCT '_' x.SUBJECTS{iS} '.nii']);
            tempTemplate2 = fullfile(x.D.PopDir,['u_rc1' x.P.STRUCT '_' x.SUBJECTS{iS} '_T1_template.nii']);
            if xASL_exist(tempTemplate1,'file')
                xASL_Move(tempTemplate1, tempTemplate2, true);
            end
        end

    end
end

%% -------------------------------------------------------------------------------------------------------
%% Normalize the DARTEL templates to MNI (those used by CAT12)
% CAT12 uses both DARTEL & Geodesic shooting, both templates are registered
% in the same space.
% Since we ran the DARTEL here, we will register to the DARTEL template

PathDARTEL = fullfile(x.D.PopDir, 'Template_6.nii');
PathDARTEL_snMat = fullfile(x.D.PopDir, 'Template_6_sn.mat');
[~,catVer] = cat_version();
if str2double(catVer) > 1500
	catTempDir = 'templates_volumes';
else
	catTempDir = 'templates_1.50mm';
end
PathMNI = fullfile(x.D.SPMDIR,'toolbox','cat12',catTempDir,'Template_6_IXI555_MNI152.nii');

xASL_delete(PathDARTEL_snMat); % make sure that this is always repeated for new DARTEL flow fields

clear matlabbatch
matlabbatch{1}.spm.tools.oldnorm.est.subj.source = {[PathDARTEL ',1']};
matlabbatch{1}.spm.tools.oldnorm.est.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.template = {[PathMNI ',1']};
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smoref = 8;
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.reg = 1;

spm_jobman('run',matlabbatch);





%% -------------------------------------------------------------------------------------------------------
%% Combine flow fields of DARTEL & CAT12, & rerun resampling & visual QC
%  PM: create status files as well to keep track here what was already done
fprintf('%s\n', 'Concatenating SPM & DARTEL flowfields:  0%');
for iS=1:x.nSubjects
    xASL_TrackProgress(iS, x.nSubjects);
    y_file  = fullfile(x.D.ROOT, x.SUBJECTS{iS}, ['y_' x.P.STRUCT '.nii']);
    u_file  = fullfile(x.D.PopDir, ['u_rc1' x.P.STRUCT '_' x.SUBJECTS{iS} '_' x.DARTEL_TEMPLATE '.nii']);
    y_y_file= fullfile(x.D.ROOT, x.SUBJECTS{iS},['y_y_' x.P.STRUCT '.nii']);

    if xASL_exist(y_file, 'file')

		xASL_adm_UnzipNifti(y_file, 1);
        xASL_im_FillNaNs(y_file, 3, x.Quality, [], x);

        clear matlabbatch
        matlabbatch{1}.spm.util.defs.comp{1}.def = {y_file};
        
        if xASL_exist(u_file, 'file')
            % first unzip & fill NaNs, just to be sure
            xASL_adm_UnzipNifti(u_file, 1);
            xASL_im_FillNaNs(u_file, 2);

            matlabbatch{1}.spm.util.defs.comp{end+1}.dartel.flowfield = {u_file};
            matlabbatch{1}.spm.util.defs.comp{end}.dartel.times = [1 0];
            matlabbatch{1}.spm.util.defs.comp{end}.dartel.K = 6;
            matlabbatch{1}.spm.util.defs.comp{end}.dartel.template = {''};
        end
        if exist(PathDARTEL_snMat,'file')
            matlabbatch{1}.spm.util.defs.comp{end+1}.sn2def.matname = {PathDARTEL_snMat};
            matlabbatch{1}.spm.util.defs.comp{end}.sn2def.vox = [NaN NaN NaN];
            matlabbatch{1}.spm.util.defs.comp{end}.sn2def.bb = [NaN NaN NaN; NaN NaN NaN];
        end
        
        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'y_T1.nii';
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {fullfile(x.D.ROOT, x.SUBJECTS{iS})};

        spm_jobman('run',matlabbatch);

        xASL_Move(y_y_file, y_file, true);
        xASL_im_FillNaNs(y_file, 3);
        xASL_delete(u_file);

        x.P.SubjectID = x.SUBJECTS{iS};
        x.dir.SUBJECTDIR = fullfile(x.D.ROOT, x.SUBJECTS{iS});
        x = xASL_init_FileSystem(x);

        xASL_wrp_Resample2StandardSpace(x);
        xASL_wrp_VisualQC_Structural(x);
    end
end






%% -------------------------------------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
result = true;



end




%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------
function [warp, TotalnEstimate] = GetWarpList(x, SubjectList)
%GETWARPLIST Create input list from predefined filetypes & subject-list
% Lists per subject for both DARTEL_create & DARTEL_warp, to make sure all file types are grouped correctly to same subject
% Throws error when a subject in subject_list cannot be found in DARTEL path for all required filetypes (DARTEL_create & DARTEL_warp)

    if nargin<2
        SubjectList = x.dataset.TimePointSubjects{1};
    end

    TotalnEstimate = 0;
    NextN = 1;

    for iSubject = 1:length(SubjectList)
        CellFillNr = 1;

        for ii=1:length(x.IMAGES_CREATE)
            TempPath = fullfile(x.D.PopDir,[x.IMAGES_CREATE{ii} '_' SubjectList{iSubject} '.nii']);

            if ~xASL_exist(TempPath,'file')
                error([TempPath ' was expected but did not exist!']);
            else
                xASL_io_ReadNifti(TempPath);
                warp.images{CellFillNr}(NextN,1) = {[TempPath ',1']};
                fprintf('%s\n',['Fileset ' num2str(ii) ' file ' char(TempPath) ' will be estimated']);
                TotalnEstimate = TotalnEstimate + 1;
                CellFillNr = CellFillNr + 1;
            end
        end
        NextN = NextN+1;
    end

end




%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------
function [warp] = GeneralDARTELsettings(x,warp,iIT)
%GeneralDARTELsettings Common DARTEL estimate settings (with or without existing templates)

    warp.settings.optim.lmreg = 0.01;
    warp.settings.optim.cyc = 3;

    if  x.Quality==1
        rparam1 = [2   1    0.5   0.25   0.25  0.125];
        rparam2 = [1   0.5  0.25  0.125  0.125 0.0625];
        K       = [0   2    4     6      8     16];
    else
        rparam1 = [4  2   1  0.5 0.25   0.25];
        rparam2 = [2  1 0.5 0.25 0.125 0.125];
        K       = [0  0   2    4    6      8];
    end

    fprintf('\n\n\n');
    fprintf(['\nDARTEL is run, creating template ' x.DARTEL_TEMPLATE '\n\n']);

    if isfield(warp.settings,'param')
        warp.settings = rmfield(warp.settings,'param');
    end

    warp.settings.rform = 0;
    warp.settings.optim.its = 3;

    warp.settings.param(1).its = 3;
    warp.settings.param(1).rparam = [rparam1(iIT) rparam2(iIT) 1e-06];
    warp.settings.param(1).K = K(iIT);

end
