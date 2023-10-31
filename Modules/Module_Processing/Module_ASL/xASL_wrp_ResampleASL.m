function xASL_wrp_ResampleASL(x)
%xASL_wrp_ResampleASL Submodule of ExploreASL ASL Module, that reslices native
%space images to standard space
%
% FORMAT: xASL_wrp_ResampleASL(x)
%
% INPUT:
%   x  - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT:
% OUTPUT FILES: NIfTIs in native & standard space
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule resamples native space NIfTIs to standard space, using the deformation fields computed in the structural module
%              after smoothing these transformation fields to the ASL resolution.
%              The applied interpolation is a combination of all transformations (e.g. motion correction, registration to
%              T1w, and transformation of T1w to standard space. This submodule performs the following steps:
%
%               0. Administration   
%               1. Warp TopUp QC files
%               2. Create slice gradient image for quantification reference, in case of 2D ASL
%               3. Resample ASL time series to MNI space
%               4. Resample to native space (applying any motion correction or registration)
%               5. Bilateral filter (currently disabled)
%               6. Create mean control image, if available, in native & standard space
%               7. Clone mean control image to be used as pseudo-M0 (if x.Q.M0==UseControlAsM0)
%               8. Pair-wise subtraction & saving PWI & PWI4D in both spaces
%               9. Save PWI NIfTI & time-series-related maps (SD, SNR)
%               10. Delete temporary files
%               11. Report spatial CoV as QC
%
% EXAMPLE: xASL_wrp_ResampleASL(x);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL



%% ------------------------------------------------------------------------------------------
% 0. Administration

% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D, 'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end
tempnii = xASL_io_Nifti2Im(x.P.Path_despiked_ASL4D);
nVolumes = size(tempnii, 4);

xASL_im_CreateASLDeformationField(x); % make sure we have the deformation field in ASL resolution

%% ------------------------------------------------------------------------------------------
% 1. Warp TopUp QC files
PathB0 = fullfile(x.dir.SESSIONDIR ,'B0.nii');
PathUnwarped = fullfile(x.dir.SESSIONDIR ,'Unwarped.nii');
PathPopB0 = fullfile(x.D.PopDir, ['rASL_B0_' x.P.SubjectID '_' x.P.SessionID '.nii']);
PathPopUnwarped = fullfile(x.D.PopDir, ['rASL_Unwarped_' x.P.SubjectID '_' x.P.SessionID '.nii']);

InputPaths = {PathB0, PathUnwarped};
OutputPaths = {PathPopB0, PathPopUnwarped};
xASL_spm_deformations(x, InputPaths, OutputPaths, [], [], [], x.P.Path_y_ASL);

%% ------------------------------------------------------------------------------------------
% 2. Create slice gradient image for quantification reference, in case of 2D ASL
xASL_im_CreateSliceGradient(x);


%% ------------------------------------------------------------------------------------------
% 3. Resample ASL time series to MNI space (currently 1.5 mm^3)
fprintf('%s\n','Convert ASL time series to single precision');
% Convertion to single precision for precious data, to avoid
% digitization artifacts in the spatial processing
xASL_io_SaveNifti(x.P.Path_despiked_ASL4D, x.P.Path_despiked_ASL4D, tempnii, 32, 0);

% Resample to standard space
if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % BACKWARDS COMPATIBILITY, AND NOW NEEDED ALSO WHEN DCT APPLIED ON TOP OF AFFINE
    AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
else
    AffineTransfPath = [];
end

xASL_spm_deformations(x, x.P.Path_despiked_ASL4D, x.P.Path_rtemp_despiked_ASL4D, 4, [], AffineTransfPath, x.P.Path_y_ASL);


%% ------------------------------------------------------------------------------------------
% 4. Resample to native space (applying any motion correction or registration)
if nVolumes>1
    for iS=1:nVolumes
        matlabbatch{1}.spm.spatial.realign.write.data{iS,1} = [x.P.Path_despiked_ASL4D ',' num2str(iS)];
    end

    matlabbatch{1}.spm.spatial.realign.write.roptions.which     = [2 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp    = 4;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap      = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask      = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix    = 'r';

    spm_jobman('run',matlabbatch); % this applies the motion correction in native space
    
    [Fpath, Ffile] = xASL_fileparts(x.P.Path_despiked_ASL4D);
    x.P.Path_rdespiked_ASL4D = fullfile(Fpath, ['r' Ffile '.nii']);
else
    x.P.Path_rdespiked_ASL4D = x.P.Path_despiked_ASL4D;
end


%% ------------------------------------------------------------------------------------------
% 5. Bilateral filter (currently disabled)
    if x.settings.BILAT_FILTER>0
        warning('The bilateral filter is currently disabled, as it needs more testing');
    end
%     %% ------------------------------------------------------------------------------------------
%     %% Bilateral filter to remove smooth artifacts
%     if ~isdeployed % skip this part for compilation, to avoid DIP image issues
%         if  x.settings.BILAT_FILTER>0 && nVolumes>9
%             volIM = xASL_io_ReadNifti(x.P.Path_rtemp_despiked_ASL4D); % load resliced time-series
%             VoxelSize = double(volIM.hdr.pixdim(2:4));
%             volIM = single(volIM.dat(:,:,:,:,:,:,:,:));
% 
%             mask = x.S.masks.skull; % get standard space mask
%             mask(isnan(mean(volIM,4))) = 0; % remove outside FoV voxels
%             ovol = xASL_im_BilateralFilter(volIM, mask, VoxelSize, x); % run filter
% 
%             xASL_io_SaveNifti( x.P.Path_rtemp_despiked_ASL4D, x.P.Path_rtemp_despiked_ASL4D, ovol,32,0 ); % store in file
%         end
%     end



%% ------------------------------------------------------------------------------------------
% 6. Create mean control image, if available, in native & standard space
if  nVolumes>1    
	% Obtain the mean control image
	imMeanControl = xASL_quant_GetMeanControl(x, xASL_io_Nifti2Im(x.P.Path_rdespiked_ASL4D));
	xASL_io_SaveNifti(x.P.Path_rdespiked_ASL4D, x.P.Path_mean_control, imMeanControl, [], 0);
	
    % Transform mean control to standard space
    if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') 
        % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
        AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
    else
        AffineTransfPath = [];
    end

    xASL_spm_deformations(x, {x.P.Path_mean_control}, {x.P.Pop_Path_mean_control}, [], [], AffineTransfPath, x.P.Path_y_ASL);
    
    % Visual check of M0-pGM registration for masking
    xASL_vis_CreateVisualFig(x, {x.P.Pop_Path_mean_control x.P.Pop_Path_rc1T1}, x.D.M0regASLdir,[0.5 0.2]);
else
    fprintf('%s\n', 'Creation mean control image was skipped, because there was only 1 volume');
end


%% ------------------------------------------------------------------------------------------
% 7. Clone mean control image to be used as pseudo-M0 (if
% x.Q.M0==UseControlAsM0) -> %%%%% check if the 4d value doesnt break this

% If x.Q.M0 is set as UseControlAsM0, this mean control NIfTI will be
% cloned to an M0 NIfTI (and processed in the M0 submodule)
if strcmpi(x.Q.M0, 'UseControlAsM0')
    if nVolumes==1
        warning('Cant clone mean control NIfTI as M0.nii, timeseries missing');
        % we assume that a single ASL image is already subtracted, so does not
        % contain control (== raw image) information (e.g. GE 3D spiral)            
    elseif ~xASL_exist(x.P.Path_mean_control)
        error('Cant clone mean control NIfTI as M0.nii, NIfTI missing. Something went wrong earlier in this function');
    else
        % Here we clone the mean control image as M0 image, which has perfect registration with PWI
        % Clone the M0, overwrite in case of previous processing            
        fprintf('%s\n','Cloning mean control image for use as M0 (overwriting if needed)');
        fprintf('%s\n',[x.P.Path_mean_control '->' x.P.Path_M0]);

        % First backup any pre-existing M0, if it is not already performed in previous run
        if xASL_exist(x.P.Path_M0, 'file')
            fprintf('%s\n', ['Pre-existing M0 NIfTI detected: ' x.P.Path_M0]);

            if ~xASL_exist(x.P.Path_M0_backup, 'file')
                fprintf('%s\n', ['Backing up ' x.P.Path_M0 ' to ' x.P.Path_M0_backup]);
                xASL_Move(x.P.Path_M0, x.P.Path_M0_backup);
            else
                fprintf('%s\n', ['Not backing this up, backup already existed:' x.P.Path_M0_backup]);
            end
        end

        % Do the same for the JSON sidecar & legacy parms.mat
        if exist(x.P.Path_M0_parms_mat, 'file') && ~exist(x.P.Path_M0_backup_parms_mat, 'file')
            xASL_Move(x.P.Path_M0_json, x.P.Path_M0_backup_json);
        end
        if exist(x.P.Path_M0_parms_mat, 'file') && ~exist(x.P.Path_M0_backup_parms_mat, 'file')
            xASL_Move(x.P.Path_M0_parms_mat, x.P.Path_M0_backup_parms_mat);
        end

        xASL_Copy(x.P.Path_mean_control, x.P.Path_M0, true);
        % Do the same for the JSON sidecar & legacy parms.mat
        if exist(x.P.Path_ASL4D_json, 'file')
            xASL_Copy(x.P.Path_ASL4D_json, x.P.Path_M0_json, true);
        end
        if exist(x.P.Path_ASL4D_parms_mat, 'file')
            xASL_Copy(x.P.Path_ASL4D_parms_mat, x.P.Path_M0_parms_mat, true);
        end
    end
end


%% ------------------------------------------------------------------------------------------
% 8. Pair-wise control-label subtraction
% for both native & standard space

% Nomenclature:
% ASL4D = original timeseries
% PWI4D = subtracted timeseries -> for quantification purposes
% PWI3D = averaged per TE—PLD—LabDur combination -> for QC purposes
% PWI   = single perfusion-weighted volume -> for registration purposes

% 

% Load ASL time series (after being pre-processed)
StringSpaceIs = {'native' 'standard'};
PathASL4D = {x.P.Path_rdespiked_ASL4D x.P.Path_rtemp_despiked_ASL4D};
PathPWI = {x.P.Path_PWI x.P.Pop_Path_PWI};
PathPWI3D = {x.P.Path_PWI3D x.P.Pop_Path_PWI3D};
PathPWI4D = {x.P.Path_PWI4D x.P.Pop_Path_PWI4D};

for iSpace=1:2
    fprintf('%s\n', ['Saving in ' StringSpaceIs{iSpace} ' space:']);
    
    ASL_im = xASL_io_Nifti2Im(PathASL4D{iSpace}); % Load time-series nifti
    dim4 = size(ASL_im, 4);

    
    % =====================================================================
    % A) Check subtraction
    if numel(size(ASL_im))>4
        warning('In BIDS ASL NIfTIs should have [X Y Z n/PLD/TE/etc] as 4 dimensions');
        error(['NIfTI has more than 4 dimensions: ' PathASL4D{iSpace}]);
        
    elseif dim4>1 && round(dim4/2)~=dim4/2
        warning('Odd (i.e., not even) number of control-label pairs, skipping');
        return;
        
        
    % =====================================================================
    % B) Time-Encoded subtraction
    elseif x.modules.asl.bTimeEncoded
        % Decoding of TimeEncoded data - it outputs a decoded image and it also saves as a NII
        ASL_im = xASL_quant_HadamardDecoding(PathASL4D{iSpace}, x.Q);
		
		% Calculate Hadamard block size (number of unique PLDs * multiTEs) = number of volumes per block
        vectorSizePLD = length(unique(x.Q.Initial_PLD));
        blockSize = x.Q.NumberEchoTimes * (vectorSizePLD - (vectorSizePLD/x.Q.TimeEncodedMatrixSize)); 
		
        if blockSize ~= size(ASL_im,4) % In case there is more than 1 block, we will average across repetitions
			if mod(size(ASL_im,4), blockSize)
				% Number of total volumes cannot be divided by blockSize
				error(['Total number of volumes ' xASL_num2str(size(ASL_im,4)) ' cannot be composed of blocks of size ' xASL_num2str(blockSize)]);
			else
				PWI4D = zeros(size(ASL_im,1), size(ASL_im,2), size(ASL_im,3), blockSize); % preallocate PWI
				for iBlock = 1:blockSize
					PWI4D(:,:,:,iBlock) = xASL_stat_MeanNan(ASL_im(:,:,:,iBlock:blockSize:end), 4); % Averaged PWI4D across repetitions
				end
			end
        else % if Hadamard Block size == unique PLDs * number of multi-TEs -> a single repetitions -> we don't want to average
            PWI4D = ASL_im;
        end
        
        % Create single PWI for further steps in ASL module
        PWI3D = xASL_stat_MeanNan(PWI4D(:,:,:,1:x.Q.NumberEchoTimes:end),4); % Average across PLDs from each first TE
        
        
    % =====================================================================
    % C) Single- and multi-PLD subtraction
    else
        % Firstly, equalize PLD for single- and multi-PLD
        nPLD = length(x.Q.Initial_PLD);
        factorPLD = nVolumes/nPLD;
        initialPLD = repmat(x.Q.Initial_PLD, [factorPLD 1]);
        % this vector now has the length of the number of volumes
        % either all values are identical (in the case of single-PLD)
        % or a combination of multi-PLD
        
        % Secondly, equalize LabDur for single- and multi-PLD
        nLabDur = length(x.Q.LabelingDuration);
        factorLabDur = nVolumes/nLabDur;
        LabDuration = repmat(x.Q.LabelingDuration, [factorLabDur 1]);


        %% 1) Create PWI4D
        if dim4>1 && ~x.modules.asl.bContainsDeltaM
            % Paired subtraction
            [ControlIm, LabelIm] = xASL_quant_GetControlLabelOrder(ASL_im);
            PWI4D = ControlIm - LabelIm;
            
            % Skip every other value in x.Q.Initial_PLD as it was stored for both control and label images 
            % we need the PLD vector now for the pairwise subtractions only            
            initialPLD_PWI4D = initialPLD(1:2:end);
            % Do the same for LabDur
            LabDuration_PWI4D = LabDuration(1:2:end);
        else % the same but then without subtraction
            PWI4D = ASL_im;
            
            initialPLD_PWI4D = initialPLD;
            LabDuration_PWI4D = LabDuration;
        end
       
        %% 2) Create PWI3D

		% After averaging across PLDs, we'll obtain these unique PLDs+LD combinations
		% indexAverage_PLD_LabDur lists for each original position to where it should be averaged
		[~, ~, iUnique_PLD_LabDur] = unique(unique([initialPLD_PWI4D, LabDuration_PWI4D], 'stable', 'rows'), 'stable', 'rows');

        % MultiPLD-multiLabDur PWI3D after averaging
        for iPLD_LabDur = 1:max(iUnique_PLD_LabDur)
            PWI3D(:, :, :, iPLD_LabDur) = xASL_stat_MeanNan(PWI4D(:, :, :, iUnique_PLD_LabDur == iPLD_LabDur), 4); % Averaged PWI4D 
        end

        %% 3) Create PWI
        % We create a dummy CBF image for registration purposes
        % The earliest echo, and latest PLD-labeling duration combination
        % are the best for this, having most SNR and CBF-weighting,
        % respectively

        PWI = PWI3D(:, :, :, iPLD_LabDur==max(iUnique_PLD_LabDur));

    end
    
    
    % =====================================================================
    % D) Save subtracted to disk

    % Save PWI (single volume)
    fprintf('%s\n', PathPWI{iSpace});
    xASL_io_SaveNifti(PathASL4D{iSpace}, PathPWI{iSpace}, PWI, 32, false);

    % Save PWI3D (averaged volume for each PLD-labdur combination)
    fprintf('%s\n', PathPWI{iSpace});
    xASL_io_SaveNifti(PathASL4D{iSpace}, PathPWI3D{iSpace}, PWI3D, 32, false);
    
    % Save PWI4D (subtracted only, not yet averaged)
    fprintf('%s\n', PathPWI4D{iSpace});
    xASL_io_SaveNifti(PathASL4D{iSpace}, PathPWI4D{iSpace}, PWI4D, 32, false);    
end


%% ------------------------------------------------------------------------------------------
% 9. Save PWI NIfTI & time-series-related maps (SD, SNR)
MaskIM = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir, 'rgrey.nii'));
MaskIM = MaskIM>(0.7*max(MaskIM(:)));

if size(PWI4D, 4)>3
    %% PM: here we need to use the positions for the earliest echo, and the latest PLD-labdur
    
    % Save SD & SNR maps
    SD = xASL_stat_StdNan(PWI4D, 0, 4); 
    SNR = PWI4D./SD;
    SNR(SNR<0) = 0; % clip @ zero    
    
    xASL_io_SaveNifti(x.P.Path_PWI, x.P.Pop_Path_SD, SD ,[], 0);
    xASL_io_SaveNifti(x.P.Path_PWI, x.P.Pop_Path_SNR, SNR, [], 0);
    fprintf('%s\n','Standard space SD & SNR saved, & part1 & part2 for reproducibility');

    
    % Calculate two halves, for reproducibility (wsCV) CBF & spatial CoV
    HalfVol                         = round(size(PWI4D,4)/2);
    Part1                           = xASL_stat_MeanNan(PWI4D(:,:,:,1:HalfVol),4);
    Part2                           = xASL_stat_MeanNan(PWI4D(:,:,:,HalfVol+1:end),4);

    % Determine mask
    if isequal(size(MaskIM),size(Part1),size(Part2))
        MaskIM = MaskIM & isfinite(Part1) & isfinite(Part2);
    else
        error('Dimension mismatch between mask and part 1 & part 2...');
    end
    
    % Calculate mean, standard deviation and covariance
    mean1                           = mean(Part1(MaskIM));
    std1                            = std(Part1(MaskIM));
    CoV1                            = std1/mean1;
    mean2                           = mean(Part2(MaskIM));
    std2                            = std(Part2(MaskIM));
    CoV2                            = std2/mean2;

    wsCV_mean                       = 100*abs(mean1-mean2)/(0.5*(mean1+mean2));
    wsCV_CoV                        = 100*abs(CoV1-CoV2)/(0.5*(CoV1+CoV2));

    fprintf('%s\n',['pGM>0.7: wsCV mean CBF = ' num2str(wsCV_mean,3) '%, wsCV spatial CoV = ' num2str(wsCV_CoV,3) '%'])

    xASL_io_SaveNifti(x.P.Path_PWI,fullfile(x.D.PopDir,['PWI_part1_'  x.P.SubjectID '_' x.P.SessionID '.nii']), Part1, 32);
    xASL_io_SaveNifti(x.P.Path_PWI,fullfile(x.D.PopDir,['PWI_part2_'  x.P.SubjectID '_' x.P.SessionID '.nii']), Part2, 32);
    fprintf('%s\n','Also saved part1 & part2 for reproducibility');
    
else
    fprintf('%s\n',['Standard space SD, SNR & reproducibility part1 & part2 maps were skipped, had only ' num2str(size(PWI4D,4)) ' volumes(s)']);
end


%% ------------------------------------------------------------------------------------------
% 10. Delete temporary files
if x.settings.DELETETEMP
    if ~strcmp(x.P.Path_ASL4D, x.P.Path_rdespiked_ASL4D)
        % in case of single volumes, these can be set to the same NIfTI
        xASL_adm_DeleteFilePair(x.P.Path_rdespiked_ASL4D, 'mat');
    end
    xASL_delete(x.P.Path_rtemp_despiked_ASL4D);
    xASL_adm_DeleteFilePair(x.P.Path_temp_despiked_ASL4D, 'mat');
else
    xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D, x.P.Path_rtemp_despiked_ASL4D, ASL_im, 32);
end


%% ------------------------------------------------------------------------------------------
% 11. Report spatial CoV as QC
%% PM: here we also need PWI4D for the first TE and last PLD/labdur

CoV0 = xASL_stat_ComputeSpatialCoV(PWI4D, MaskIM)*100;
fprintf('%s\n',['Standard space spatial CoV pGM>0.7 = ' num2str(CoV0,3) '%']);


end