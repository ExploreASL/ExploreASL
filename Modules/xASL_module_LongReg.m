function [result, x] = xASL_module_LongReg(x)
%xASL_module_LongReg ExploreASL module for Longitudinal Registration
%
% FORMAT: [result, x] = xASL_module_LongReg(x)
%
% INPUT:
%   x  - x structure containing all input parameters (REQUIRED)
%   x.dir.SUBJECTDIR  -  anatomical directory, containing the derivatives of anatomical images (REQUIRED)
%   x.WhichLongReg  - (OPTIONAL, DEFAULT = 'LongReg'), when this is set to 'DARTEL', DARTEL is performed instead of SPM12 Longitudinal Registration
%
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters (REQUIRED)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This module performs a non-linear registration within-subject, between time points.
% This is very useful for longitudinal studies, where care needs to be taken that longitudinal analyses are not
% biased by registration errors within subjects. The idea is that registration will always make a (small) error,
% and also interpolation will induce some error. Different time points from the same subject should only
% differ in local or global volumetric changes, whereas differences between subjects can be geometrically much larger, requiring
% larger transformations and interpolations. By splitting the between-subject non-linear registration (e.g. DARTEL/Geodesic Shooting within
% the structural module/CAT12) and the within-subject non-linear registration performed in this module, the statistical bias of the non-linear
% registration can be controlled for longitudinal analyses.
%
% Two algorithms can be used for longitudinal registration:
% SPM12 Longitudinal Registration (default): modified version of DARTEL, restricting the non-linear registration between time points
%                                            to growing & shrinking only. Works very well for atrophy changes only, but fails in the presence of lesions
%                                            or when there are extracranial changes (e.g. growing children). Attempt to mitigate this (by brainmasking before
%                                            running Longitudinal Registration) did not fully work
% SPM12 DARTEL                             : DARTEL providing full non-linear freedom, could potentially overfit, but more robust to lesions
%                                            (when properly lesion-masked) or extracranial growth.
% This module will run only for the first T1w/time point of each subject, if there are multiple time points
% Submodules:
% 1) Longitudinal registration = after repetition of co-registrations, T1w images are resampled to 1.5x1.5x1.5 mm and registration is started
% 2) Visual QC of the longitudinal registration, in which also any existing transformations (e.g. from T1w to MNI) are concatenated/combined with the newly
%    obtained transformations
% 3) Repetition of resampling of images to standard space & visual QC (as performed as submodules 9 & 12 in the structural module)
% EXAMPLE: [result, x] = xASL_module_LongReg(x)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
% Ashburner J, Ridgway GR. Symmetric diffeomorphic modeling of longitudinal structural MRI. Front Neurosci. 2012. p. 197.
% Ashburner J. A fast diffeomorphic image registration algorithm. Neuroimage. 2007;38(1):95
% __________________________________
% Copyright ? 2015-2019 ExploreASL



%% --------------------------------------------------------------
%% Administration
[x] = xASL_init_SubStructs(x);
x = xASL_init_InitializeMutex(x, 'LongReg' ); % starts mutex locking process to ensure that everything will run only once

% 0.9 change working directory to make sure that unspecified output will go there...
oldFolder = cd(x.dir.SUBJECTDIR);

% Check T1 presence

structList = xASL_adm_GetFileList(x.dir.SUBJECTDIR, ['^' x.P.STRUCT '\.nii$'],'List');

if isempty(structList)
    error('AslPipeline:missingFile', 'No file found that matches "%s" in "%s"', x.P.STRUCT, x.dir.SUBJECTDIR);
elseif length(structList)>2
    error('AslPipeline:incorrectNrOfFiles', 'Expected exactly one nifti file "%s", but found %d', x.P.STRUCT, length(structList));
end

[~, x.P.SubjectID] = fileparts(x.dir.SUBJECTDIR);
[~, ~, ~, ~, VolumeList, VolumeN] = xASL_init_LongitudinalRegistration(x);

VolumeList = logical(VolumeList(x.dataset.nSessions:x.dataset.nSessions:end,1));
CurrentSub = x.SUBJECTS(VolumeList);

iSubject = find(strcmp(x.SUBJECTS,x.P.SubjectID));
iSession = 1; % append sessions to accommodate sessions in x.S.SetsID
iSubjSess = ((iSubject-1)*x.dataset.nSessions)+iSession;

if ~isfield(x,'WhichLongReg') || isempty(x.WhichLongReg) || ~strcmpi(x.WhichLongReg,'DARTEL')
    x.WhichLongReg = 'LongReg'; % default
end

% To check whether or not we will run longitudinal registration
% length(VolumeN)==1 -> no Longitudinal Registration
% length(VolumeN)>1  -> do Longitudinal Registration
% VolumeN provides the numbers of the TimePoints (this can be [1 2 3] for
% three subsequent TimePoints, but could also be [2 3] if the first
% TimePoint was excluded from ExploreASL




%% --------------------------------------------------------------
%% 1    Longitudinal registration
% cave: check time intervals
% This part does not change the T1.nii, neither does CAT12 or SPM12
% segmentation, so longitudinal registration results can be kept if
% segmentation needs to be run again

% Backwards compatibility check
% Either the LongReg transformations should exist, or they are already combined into a y_T1.nii.
% Otherwise we should rerun this module

if  strcmp(x.P.SubjectID,CurrentSub{1}) && length(VolumeN)>1 % only perform if this is the first time point of multiple volumes
    if  x.mutex.HasState('001_LongitudinalRegistration')
        fprintf('%s\n','001_LongitudinalRegistration           has already been performed, skipping...');
    else
        fprintf('%s\n','Performing SPM12 longitudinal registration');


        %% ---------------------------------------------------------------------------------
        %% Check first if all T1w volumes were aligned to MNI ACPC
        RegLockFileCount    = 0;

        for iS=1:length(CurrentSub) % count aligned T1w volumes
            LockFile = '';
            LockFile = fullfile(x.D.ROOT, 'lock', 'xASL_module_Structural',CurrentSub{iS}, 'xASL_module_Structural', '010_LinearReg_T1w2MNI.status');

            if ~isempty(LockFile) && xASL_exist(LockFile, 'file')
                RegLockFileCount = RegLockFileCount+1;
            end
        end

        if  RegLockFileCount~=length(CurrentSub)
            warning('Not all T1w volumes were aligned to MNI ACPC, skipping longitudinal registration');
            return;
        else

            fprintf('%s\n','Running longitudinal registration');

			% Unzip & remove redundant files
            xASL_io_ReadNifti(fullfile(x.dir.SUBJECTDIR, [x.P.STRUCT '.nii']));



            %% ---------------------------------------------------------------------------------
            %% Now register all other volumes to the first
            %  Here we also reslice the first TimePoint, but below we only register the later TimePoints
            % We



            for iV=1:length(CurrentSub)
                % Define paths
                Path_T1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'T1.nii');
                Path_rT1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'rT1.nii');
%                 Path_c1T1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'c1T1.nii');
%                 Path_rc1T1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'rc1T1.nii');
%                 Path_c2T1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'c2T1.nii');
%                 Path_rc2T1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'rc2T1.nii');
%                 Path_c3T1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'c3T1.nii');
%                 Path_rc3T1{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'rc3T1.nii');
%                 Path_Bmask{iV} = fullfile(x.D.ROOT, CurrentSub{iV},'rBrainMask.nii');

                % Reslice
                xASL_spm_reslice(x.D.ResliceRef, Path_T1{iV}, [], [], x.settings.Quality);
%                 xASL_spm_reslice(x.D.ResliceRef, Path_c1T1{iV}, [], [], x.settings.Quality, [], 1);
%                 xASL_spm_reslice(x.D.ResliceRef, Path_c2T1{iV}, [], [], x.settings.Quality, [], 1);
%                 if xASL_exist(Path_c3T1{iV},'file')
%                     xASL_spm_reslice(x.D.ResliceRef, Path_c3T1{iV}, [], [], x.settings.Quality, [], 1);
%                 end

%                 % Create brainmask
%                 BmaskIM = xASL_io_Nifti2Im(Path_rc1T1{iV}) + xASL_io_Nifti2Im(Path_rc2T1{iV});
%                 if exist(Path_rc3T1{iV},'file')
%                     BmaskIM = BmaskIM + xASL_io_Nifti2Im(Path_rc3T1{iV});
%                     BmaskIM = BmaskIM>0.25;
%                 else
%                     BmaskIM = BmaskIM>0.25;
%                     BmaskIM = xASL_im_DilateErodeFull(BmaskIM, 'dilate', xASL_im_DilateErodeSphere(1));
%                 end


                % Create common brainmask

                xASL_im_SkullStrip(Path_rT1{iV}, fullfile(x.D.MapsSPMmodifiedDir, 'rbrainmask.nii'), x, Path_rT1{iV}); % skullstrip
                tIM = xASL_io_Nifti2Im(Path_rT1{iV});
                tIM(tIM<=0) = NaN;
                xASL_io_SaveNifti(Path_rT1{iV}, Path_rT1{iV}, tIM, [], false);

                FirstStruct = fullfile(x.D.ROOT, CurrentSub{1},['r'  x.P.STRUCT '.nii']);

                OtherList{1,1} = fullfile(x.D.ROOT,CurrentSub{iV}, [x.P.STRUCT '.nii']);
                OtherList{end+1,1} = fullfile(x.D.ROOT,CurrentSub{iV}, ['c1' x.P.STRUCT '.nii']);
                OtherList{end+1,1} = fullfile(x.D.ROOT,CurrentSub{iV}, ['c2' x.P.STRUCT '.nii']);
                OtherList{end+1,1} = fullfile(x.D.ROOT,CurrentSub{iV}, [x.P.FLAIR '.nii']);
                OtherList{end+1,1} = fullfile(x.D.ROOT,CurrentSub{iV}, [x.P.WMH_SEGM '.nii']);
				OtherList{end+1,1} = fullfile(x.D.ROOT,CurrentSub{iV}, [x.P.T1c '.nii']);
				OtherList{end+1,1} = fullfile(x.D.ROOT,CurrentSub{iV}, [x.P.T2 '.nii']);

                for iSess = 1:x.dataset.nSessions
                    SessionDir = fullfile(x.D.ROOT,CurrentSub{iV},x.SESSIONS{iSess});
                    OtherList{end+1,1} = fullfile(SessionDir, [x.P.ASL4D '.nii']);
                    OtherList{end+1,1} = fullfile(SessionDir, [x.P.M0 '.nii']);
                end

                if  iV~=1
                    xASL_spm_coreg(FirstStruct, Path_rT1{iV}, OtherList, x);
                end
            end



            %% ---------------------------------------------------------------------------------
            %% Default settings longitudinal registration
            LongRegBatch{1}.spm.tools.longit{1}.series.noise     = NaN;
            LongRegBatch{1}.spm.tools.longit{1}.series.wparam    = [0 0 100 25 100];

            LongRegBatch{1}.spm.tools.longit{1}.series.write_avg = false; % save average image
            LongRegBatch{1}.spm.tools.longit{1}.series.write_jac = false; % save Jacobians -> this is repeated separately for the combined flowfield in xASL_wrp_Resample
            LongRegBatch{1}.spm.tools.longit{1}.series.write_div = false; % save divergence rate (contraction=atrophy/expansion)
            LongRegBatch{1}.spm.tools.longit{1}.series.write_def = true; % save deformation fields

            % Specify biasfield regularization (if GE, no regularization, otherwise
            % keep default)
            LongRegBatch{1}.spm.tools.longit{1}.series.bparam    = 1000000;
            if isfield(x,'Q')
                if isfield(x.Q,'Vendor')
                    if ~isempty(regexpi(x.Q.Vendor,'GE'))
                        LongRegBatch{1}.spm.tools.longit{1}.series.bparam    = 0;
                    end
                end
            end

            % Get age, if exists, otherwise use 4.5 time change as default.
            % This is an optimum for modeling deformations & computation time spent

            if  isfield(x, 'iSetAge')
                IndicesN    = [iSubjSess:x.dataset.nSessions:iSubjSess+x.dataset.nSessions*(length(CurrentSub)-1)];
                AgeN        = x.S.SetsID( IndicesN ,x.iSetAge);
            else
                AgeN        = [50:4.5:(50+4.5*(length(CurrentSub)-1))];
            end

            %% ---------------------------------------------------------------------------------
            %% Default settings DARTEL
            DARTELbatch{1}.spm.tools.dartel.warp.settings.template = 'Template';
            DARTELbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
            DARTELbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;

            %% ---------------------------------------------------------------------------------
            %% Firstly, clone 1st y_T1.nii, reslice & define all volumes:
            for iV=1:length(CurrentSub)
                Path_T1 = fullfile(x.D.ROOT, CurrentSub{iV},'T1.nii');
                Path_c1T1 = fullfile(x.D.ROOT, CurrentSub{iV},'c1T1.nii');
                Path_c2T1 = fullfile(x.D.ROOT, CurrentSub{iV},'c2T1.nii');
                Path_c3T1 = fullfile(x.D.ROOT, CurrentSub{iV},'c3T1.nii');
                Path_rT1 = fullfile(x.D.ROOT, CurrentSub{iV},'rT1.nii');
                Path_rc1T1 = fullfile(x.D.ROOT, CurrentSub{iV},'rc1T1.nii');
                Path_rc2T1 = fullfile(x.D.ROOT, CurrentSub{iV},'rc2T1.nii');
                Path_rc3T1 = fullfile(x.D.ROOT, CurrentSub{iV},'rc3T1.nii');

                xASL_delete(Path_rT1);

                if strcmpi(x.WhichLongReg,'DARTEL')
                    xASL_spm_reslice(x.D.ResliceRef, Path_c1T1, [], [], x.settings.Quality, [], 1);
                    xASL_spm_reslice(x.D.ResliceRef, Path_c2T1, [], [], x.settings.Quality, [], 1);
                    if exist(Path_c3T1,'file')
                        xASL_spm_reslice(x.D.ResliceRef, Path_c3T1, [], [], x.settings.Quality, [], 1);
                    end

                    % Remove NaNs for DARTEL
                    Files2 = {Path_rc1T1 Path_rc2T1 Path_rc3T1};
                    for iL=1:length(Files2)
                        tIM = xASL_io_Nifti2Im(Files2{iL});
                        tIM(isnan(tIM)) = 0;
                        xASL_io_SaveNifti(Files2{iL}, Files2{iL}, tIM, [], false);
                    end
                else
                    xASL_spm_reslice(x.D.ResliceRef, Path_T1, [], [], x.settings.Quality);
                    xASL_im_SkullStrip(Path_rT1, fullfile(x.D.MapsSPMmodifiedDir, 'rbrainmask.nii'), x, Path_rT1); % skullstrip
                    tIM = xASL_io_Nifti2Im(Path_rT1);
                    tIM(tIM<=0) = NaN;
                    xASL_io_SaveNifti(Path_rT1, Path_rT1, tIM, [], false);
                end
                % we have to repeat the reslicing here, because of the previous rigid-body registration

                LongRegBatch{1}.spm.tools.longit{1}.series.vols{iV,1} = Path_rT1;
                LongRegBatch{1}.spm.tools.longit{1}.series.times(iV,1) = AgeN(iV);
                DARTELbatch{1}.spm.tools.dartel.warp.images{1}{iV,1} = [Path_rc1T1 ',1'];
                DARTELbatch{1}.spm.tools.dartel.warp.images{2}{iV,1} = [Path_rc2T1 ',1'];
                if exist(Path_rc3T1,'file')
                    DARTELbatch{1}.spm.tools.dartel.warp.images{3}{iV,1} = [Path_rc3T1 ',1'];
                end
            end

            %% ---------------------------------------------------------------------------------
            % Backup deformation field y_T1.nii for the first time point
            TempRegFile = fullfile( x.D.ROOT, CurrentSub{1}, ['Temp_y_' x.P.STRUCT '.nii']);
            DeformationFile = fullfile( x.D.ROOT, CurrentSub{1}, ['y_' x.P.STRUCT '.nii']);

            if xASL_exist(DeformationFile, 'file')
                xASL_Move(DeformationFile, TempRegFile, true);
            end

            if strcmpi(x.WhichLongReg, 'DARTEL')
                spm_jobman('run',DARTELbatch);
            else % default
                spm_jobman('run',LongRegBatch);
            end

            %% ---------------------------------------------------------------------------------
            % Rename deformation file
            for iS=1:length(CurrentSub)
                if strcmpi(x.WhichLongReg,'DARTEL')
                    DeformationFile = fullfile(x.D.ROOT, CurrentSub{iS}, 'u_rc1T1_Template.nii');
                else
                    DeformationFile = fullfile(x.D.ROOT, CurrentSub{iS}, 'y_rT1.nii');
                end
                LongRegFile = fullfile(x.D.ROOT, CurrentSub{iS}, 'LongReg_y_T1.nii');
                xASL_Move(DeformationFile, LongRegFile, true);

                % return backed up deformation field
                TempRegFile = fullfile( x.D.ROOT, CurrentSub{iS}, ['Temp_y_' x.P.STRUCT '.nii']);
                TempRegFile2 = fullfile( x.D.ROOT, CurrentSub{iS}, ['y_' x.P.STRUCT '.nii']);
                if  xASL_exist(TempRegFile,'file')
                    xASL_Move(TempRegFile, TempRegFile2, true);
                end

                TempFiles = {'rT1.nii' 'rc1T1.nii' 'rc2T1.nii' 'rc3T1.nii' ...
                    'Template_0.nii' 'Template_1.nii' 'Template_2.nii' 'Template_3.nii' ...
                    'Template_4.nii' 'Template_5.nii' 'Template_6.nii' 'u_rc1T1_Template.nii'};
                for iT=1:length(TempFiles)
                    TempPath = fullfile(x.D.ROOT, CurrentSub{iS}, TempFiles{iT});
                    xASL_delete(TempPath);
                end
            end


            %% ---------------------------------------------------------------------------------
            %% 2 Visual check % this is combined with 001, as they cannot run independently
            fprintf('%s\n','Print visual QC checks');
            % First resample images

            % NB: non-linear between-subject registration (e.g. DARTEL) should be completed before this,
            % as the between-subject & within-subject transformations are combined here

            % Here we compare the effect of LongReg, by comparing
            % 1) T1w_2 without LongReg -> T1w_1
            % 2) T1w_2 with    LongReg -> T1w_1
            % In both cases, we do not apply the affine/DARTEL warp of CAT12 to
            % MNI, as to compare the images in the "linear registered native space"

            x.P.SubjectID   = CurrentSub{1};
            LongRegFile1    = fullfile(x.D.ROOT, CurrentSub{1}, ['LongReg_y_' x.P.STRUCT '.nii']);
            y_File1         = fullfile(x.D.ROOT, CurrentSub{1}, ['y_' x.P.STRUCT '.nii']);
            INPUTname       = fullfile(x.D.ROOT, CurrentSub{1}, [x.P.STRUCT '.nii']);
            Path_rT1     = fullfile(x.D.ROOT, CurrentSub{1}, ['r' x.P.STRUCT '.nii']);
            xASL_spm_deformations(x, INPUTname, Path_rT1, 1); % TimePoint 1 without LongReg

            for iS=2:length(CurrentSub)

                % Without LongReg
                y_File2         = fullfile( x.D.ROOT, CurrentSub{iS}, ['y_' x.P.STRUCT '.nii']);
                xASL_Copy(y_File1, y_File2, true);

                INPUTname       = fullfile( x.D.ROOT, CurrentSub{iS}, [    x.P.STRUCT '.nii']);
                OUTPUTname      = fullfile( x.D.ROOT, CurrentSub{iS}, ['r' x.P.STRUCT '.nii']);
                x.P.SubjectID   = CurrentSub{iS};
                % Create rT1.nii with the same y_T1.nii as the first TimePoint, i.e. without LongReg
                xASL_spm_deformations(x, INPUTname, OUTPUTname, 1);
                xASL_delete(y_File2);


                %% ---------------------------------------------------------------------------------
                %% Combine transformations

                % names other time point
                LongRegFile2    = fullfile( x.D.ROOT, CurrentSub{iS}, ['LongReg_y_' x.P.STRUCT '.nii']);

                if  xASL_exist(LongRegFile1,'file') && xASL_exist(LongRegFile2,'file') && xASL_exist(y_File1,'file') % if all files are there

                    fprintf('%s\n',['Combining longitudinal registration transformations for ' CurrentSub{iS}]);

                    % combine the longitdinal transformations & the y_T1.nii from the first TimePoint
                    % & then remove the original longitudinal transformations

                    if strcmpi(x.WhichLongReg,'DARTEL')
                        matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {LongRegFile2};
                        matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
                        matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
                        matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {''};

                        matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {LongRegFile1};% go from "average time point" to first time point
                        matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [0 1];
                        matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
                        matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {''};
                    else
                        matlabbatch{1}.spm.util.defs.comp{1}.def                    = {LongRegFile2}; % go from other time point to "average time point"
                        matlabbatch{1}.spm.util.defs.comp{2}.inv.comp{1}.def        = {LongRegFile1}; % go from "average time point" to first time point
                        matlabbatch{1}.spm.util.defs.comp{2}.inv.space              = {x.D.ResliceRef};
                    end

                    matlabbatch{1}.spm.util.defs.comp{3}.def                    = {y_File1}; % go from first time point to common space

                    matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname          = [x.P.STRUCT '.nii']; % y_ is prepended by default
                    matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {fullfile(x.D.ROOT,CurrentSub{iS})};

                    spm_jobman('run',matlabbatch);

                    xASL_delete(LongRegFile2);

                    xASL_im_FillNaNs(y_File2, 3, x.settings.Quality, [], x);
                end

                OUTPUT2name      = fullfile( x.D.ROOT, CurrentSub{iS}, ['r2' x.P.STRUCT '.nii']);
                xASL_spm_deformations(x, INPUTname, OUTPUT2name); % with LongReg

                % Perform QC check for all volumes, compare with first volume
                xASL_vis_CreateVisualLongReg(x, CurrentSub);
                xASL_delete(OUTPUTname);
                xASL_delete(OUTPUT2name);
            end
            xASL_delete(LongRegFile1);
            xASL_delete(Path_rT1);

            x.mutex.AddState('001_LongitudinalRegistration');
            x.mutex.DelState('002_Resample');

        end
    end


    %% ---------------------------------------------------------------------------------
    %% 3) Rerun Structural submodules 9 resampling & 12 visual QC
    fprintf('%s\n', 'Now rerunning structural submodules for resampling & visual QC');

    if ~x.mutex.HasState('002_Resample')
        fprintf('%s\n','Resample images of other volumes, corrected by longitudinal registration');
        %% 3) Replace standard space images
        % Resample standard spaces niftis for other time points
        % Here we only resample the other time points, whose transformations have changed with LongReg

        for iS=2:length(CurrentSub)
            x.P.SubjectID = CurrentSub{iS};
            x.dir.SUBJECTDIR = fullfile(x.D.ROOT, CurrentSub{iS});
            x = xASL_init_FileSystem(x);

            xASL_wrp_Resample2StandardSpace(x);
            xASL_wrp_VisualQC(x);
        end

        x.mutex.AddState('002_Resample');
    else
        fprintf('%s\n','002_Resample              has already been performed, skipping...');
    end

end


%% ---------------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
cd(oldFolder);
result = true;

end





%% ======================================================================================================================
%% ======================================================================================================================




function xASL_vis_CreateVisualLongReg( x, CurrentSub)
% xASL_vis_CreateVisualLongReg Creates for each Other TimePoint (TP):
% First row 1) First TP 2) Other TP
% Second row First TP with normalized difference image between TPs without (1) and with (2) Longitudinal Registration

    fprintf('%s\n','Creating comparison image for images with and without longitudinal registration');

    xASL_adm_CreateDir(x.D.LongRegCheckDir);

    x.S.CorSlices = [];
    x.S.SagSlices = [];
    x.S.TraSlices = x.S.slicesLarge;

    % Create brainmask with TimePoint 1 pGM & pWM
    pGMpath  = fullfile(x.D.PopDir,['rc1T1_' CurrentSub{1} '.nii']);
    pWMpath  = fullfile(x.D.PopDir,['rc2T1_' CurrentSub{1} '.nii']);
    Bmask    = (xASL_io_Nifti2Im(pGMpath) + xASL_io_Nifti2Im(pWMpath))>0.1;

    % This is before LongReg for the first TimePoint,
    % which we will compare all other TimePoints against,
    % for both before {1} & after {2} LongReg
    T1wOri = fullfile( x.D.ROOT, CurrentSub{1}, ['r' x.P.STRUCT '.nii']);
    imOri = RobustScaling(xASL_io_Nifti2Im(T1wOri));

    % Load files
    for iC=2:length(CurrentSub)

        % before LongReg = {1}, after LongReg = {2}
        T1wFile{1} = fullfile( x.D.ROOT, CurrentSub{iC}, ['r' x.P.STRUCT '.nii']);
        T1wFile{2} = fullfile( x.D.ROOT, CurrentSub{iC}, ['r2' x.P.STRUCT '.nii']);

        for ii=1:2
            im{ii} = RobustScaling(xASL_io_Nifti2Im(T1wFile{ii})); % scale
            DiffIM{ii} = abs(imOri - im{ii}); % subtract original from other TimePoint
            MeanIM{ii} = (imOri + im{ii})./2; % mean original & other TimePoint
            AI{ii} = DiffIM{ii} ./ MeanIM{ii} .* 100; % normalize difference image into asymmetry image
            AI{ii} = AI{ii} .*single(Bmask); % mask asymmetry image
        end

        % Make sure both AIs will get the same scale (as they are scaled separately in xASL_vis_CreateVisualFig below)
        SortValues = max(AI{1},AI{2});
        SortValues = sort(SortValues(isfinite(SortValues) & Bmask));
        ThreshMax = SortValues(round(0.99*length(SortValues)));
        ThreshMin = median(SortValues);
        for ii=1:2
            % Clip both with the same max value
            AI{ii}(AI{ii}>ThreshMax) = ThreshMax;
            % Remove all below the median asymmetries
            AI{ii}(AI{ii}<ThreshMin) = 0;
            % Square asymmetries, this helps to remove the noise and show only important differences
            AI{ii} = AI{ii}.^2;
        end




        % Plot first TimePoint & other TimePoint
        OutIM1 = xASL_vis_CreateVisualFig(x, fullfile(x.D.ROOT, CurrentSub{1}, ['r' x.P.STRUCT '.nii']));
        OutIM2 = xASL_vis_CreateVisualFig(x, fullfile(x.D.ROOT, CurrentSub{iC}, ['r2' x.P.STRUCT '.nii']));
        % Plot first TimePoint + AI image
        OutIM3 = xASL_vis_CreateVisualFig(x, {im{1}, AI{1}}, [], [], [], {x.S.gray x.S.jet256}, false);
        OutIM4 = xASL_vis_CreateVisualFig(x, {im{1}, AI{2}}, [], [], [], {x.S.gray x.S.jet256}, false);

        xASL_vis_Imwrite([OutIM1,OutIM2;OutIM3,OutIM4], fullfile( x.D.LongRegCheckDir,[CurrentSub{1} '_vs_' CurrentSub{iC} '.jpg']));

    end
end





function [IMout]    = RobustScaling( IMin)

    SortedIM    = sort(IMin(isfinite(IMin)));
    Value95     = SortedIM(round(1*length(SortedIM))); % 0.95*(length....
    IMout       = IMin./Value95.*500;

end
