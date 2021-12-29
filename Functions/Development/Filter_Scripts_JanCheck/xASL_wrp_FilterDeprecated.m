function xASL_wrp_FilterDeprecated(InputFile, OutputFile, x, MaskingOption )
%xASL_wrp_Filter





    %% ------------------------------------------------------------------------------------------
    %% Administration
    path    = fileparts(OutputFile);

    if ~exist('MaskingOption','var')
        MaskingOption   = 0; % default = no masking
    end

    [Fpath Ffile Fext]              = xASL_fileparts(InputFile);
    PreFix                          = 'TempFilter_';
    TempFile                        = fullfile(Fpath,[PreFix Ffile Fext]);

    % Create PWI & mean control images in native space
    xASL_io_PairwiseSubtraction( InputFile,OutputFile,0,0);

    if  strcmp(InputFile,x.P.Path_ASL4D) || strcmp(InputFile,x.P.Path_despiked_ASL4D)
        xASL_Copy(OutputFile,x.P.Path_PWI,1);
    else
        xASL_Copy(OutputFile,[OutputFile(1:end-4) '_noClipping.nii'],1);
    end

    xASL_im_NoiseRemoval(OutputFile);

    tIM                         = xASL_io_Nifti2Im(OutputFile);
    SortedInt                   = sort(tIM(tIM~=0 & ~isnan(tIM) ));
    MinIn                       = min(SortedInt);
    MaxIn                       = max(SortedInt);

    INTsearch                   = SortedInt(round(0.6*length(SortedInt))); % lot of voxels are noise, more than half
    tIM(tIM<INTsearch)          = 0;

    tIM                         = xASL_im_ClipExtremes( tIM,0.95,0,0);
    tIM                         = tIM - min(tIM(:));

    SortedInt                   = sort(tIM(isfinite(tIM) ));
    MinOut                      = min(SortedInt);
    MaxOut                      = max(SortedInt);

    xASL_io_SaveNifti( OutputFile, OutputFile, tIM, [], 0 );

	fprintf('%s\n',[OutputFile ' clipped from [' num2str(MinIn,3) ':' num2str(MaxIn,3) '] to [' num2str(MinOut,3) ':' num2str(MaxOut,3) ']']);

    tNII                        = xASL_io_ReadNifti(InputFile);

    %% ------------------------------------------------------------------------------------------
    %% If we need to run filter first, to remove artifact, if we have sufficient number of averages
    if  x.settings.BILAT_FILTER>0 && size(tNII.dat(:,:,:,:),4)>9
        % First reslice time-series


        VoxelSize                   = tNII.hdr.pixdim(2:4);

        for iD=1:size(tNII.dat,4)
            matlabbatch{1}.spm.spatial.realign.write.data{iD,1}     = [InputFile ',' num2str(iD)];
        end

        matlabbatch{1}.spm.spatial.realign.write.roptions.which     = [2 0];
        matlabbatch{1}.spm.spatial.realign.write.roptions.interp    = 4;
        matlabbatch{1}.spm.spatial.realign.write.roptions.wrap      = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.write.roptions.mask      = 1;
        matlabbatch{1}.spm.spatial.realign.write.roptions.prefix    = PreFix;
        spm_jobman('run', matlabbatch);

        volIM                       = xASL_io_Nifti2Im(TempFile); % load resliced time-series

        xASL_im_RegisterSkullMask(Fpath,x, OutputFile);

        mask                        = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(OutputFile));

        maskFile                    = fullfile(Fpath, 'wmask_ICV.nii');
        mask                        = xASL_io_Nifti2Im(maskFile);
        mask(isnan(mean(volIM,4)))  = 0; % remove outside FoV voxels

        if  sum(mask(:))==0 % for backward compatibility & robustness.
            % we need to go to an ASL-based mask, making the ASL processing
            % independent of a structural prior. however, this should first
            % be tested

            mask = xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(OutputFile));
        end
        %%
		%%
		%% HERE WE NEED TO RUN THE FILTER


        ovol = xASL_im_BilateralFilter(volIM, mask, VoxelSize, x); % run filter
		%%
		%%

        if  size(ovol,4)==size(tNII.dat,4) % not yet subtracted
            [control_im label_im OrderContLabl]     = xASL_quant_GetControlLabelOrder( ovol );
            ovol                                    = control_im-label_im;
        end

        xASL_io_SaveNifti( OutputFile, OutputFile, mean(ovol,4),[],0); % store in file

        xASL_delete(maskFile);
        xASL_delete(TempFile);

        % Repeat the noise removal
        NoiseRemoval(OutputFile);

    end






end




%% ------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------
function xASL_im_RegisterSkullMask(Fpath,x, OutputFile)

    ICVName     = fullfile( x.D.SPMDIR,'tpm','mask_ICV.nii'); % this is the mask it has been tested with
    GMName      = fullfile( x.D.MapsDir, 'rgrey.nii'   );

    xASL_Copy( ICVName, fullfile(Fpath, 'mask_ICV.nii') ,1);
    xASL_Copy(  GMName, fullfile(Fpath, 'rgrey.nii'   ) ,1);

    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source       = { fullfile(Fpath, 'rgrey.nii,1'   )};
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc        = '';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample     = { fullfile(Fpath, 'mask_ICV.nii,1')};
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = { [OutputFile ',1']};
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight   = '';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc   = 8;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref   = 8;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype  = 'mni';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff   = 25;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits     = 16;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg      = 1;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb       = [NaN NaN NaN
                                                                   NaN NaN NaN];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox      = [NaN NaN NaN];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp   = 0;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap     = [0 0 0];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix   = 'w';

    spm_jobman('run',matlabbatch);

    xASL_delete( fullfile(Fpath, 'mask_ICV.nii') );
    xASL_delete( fullfile(Fpath, 'rgrey.nii') );
    xASL_delete( fullfile(Fpath, 'rgrey_sn.mat') );
end





%% ------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------
