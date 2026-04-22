function xASL_wrp_RealignASL(x, bASL)
%xASL_wrp_RealignASL Submodule of ExploreASL ASL Module, that realigns
%volumes
%
% FORMAT: xASL_wrp_RealignASL(x[, bASL])
%
% INPUT:
%   x                               - structure containing fields with all information required to run this submodule (REQUIRED)
%   bASL                            - boolean that the input is a ASL-based imaging (true) or not (false, e.g. fMRI, DTI etc). (OPTIONAL, DEFAULT = true) 
%   x.modules.asl.SpikeRemovalAbsoluteThreshold   - absolute threshold for removing
%                                                  motion spike volumes, in mm (OPTIONAL, DEFAULT = 0 = disabled; e.g., 0.05)
%                                     
%
% OUTPUT: n/a (registration changes the NIfTI orientation header only,
%              with the exception of the affine transformation, which is
%              saved separately as x.P.Path_mean_PWI_Clipped_sn_mat
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule estimates motion by spm_realign, which uses a
% rigid-body registration (3 translations, 3 rotations). It runs ENABLE to
% reject outliers and provides a visualization. ENABLE, QC and visualizations
% are based on the Net Displacement Vector (NDV) (in mm):
% according to Pythagorean/Euclydian RMS
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1211&L=fsl&P=R34458&1=fsl&9=A&J=on&d=No+Match%3BMatch%3BMatches&z=4
% view this link for image of rotation roll, pitch and yaw https://www.google.nl/search?q=rotation+pitch+yaw+roll&espv=2&tbm=isch&imgil=LW3Nn1K-L6Oc7M%253A%253B-aSyykkRityJoM%253Bhttp%25253A%25252F%25252Fwww.grc.nasa.gov%25252FWWW%25252Fk-12%25252Fairplane%25252Frotations.html&source=iu&usg=__MlLQ5VuyRbm6kZP0vBJlPxmfbkw%3D&sa=X&ei=TWfjU4WcK4bqyQPqu4Fo&ved=0CD8Q9QEwBQ&biw=1680&bih=946#facrc=_&imgdii=_&imgrc=LW3Nn1K-L6Oc7M%253A%3B-aSyykkRityJoM%3Bhttp%253A%252F%252Fwww.grc.nasa.gov%252FWWW%252Fk-12%252Fairplane%252FImages%252Frotations.gif%3Bhttp%253A%252F%252Fwww.grc.nasa.gov%252FWWW%252Fk-12%252Fairplane%252Frotations.html%3B709%3B533
% 
% This submodule performs the following steps:
%
% 1. Estimate motion
%    Several options are distinguished for ASL and non-ASL image with the following options
%    A. Standard ASL with control-label pairs (opposed to non-ASL images such as fMRI/DTI), uses the zig-zag approach in which the average control vs label intensity differences are disregarded (not erroneously seen as motion)
%    B. If multiple TEs are acquired, only the shortest TE (with the highest SNR) is used for motion estimation and the same motion is applied to longer TEs in the same block
%    C. Motion outlier detection is skipped for multiple PLDs and/or multiple TEs
% 2. Calculate and plot position and motion parameters
% 3. Threshold-free spike definition (based on ENABLE, but with t-stats rather than the threshold p<0.05)
% 4. Remove spike frames from nifti
%
% EXAMPLE: xASL_wrp_RealignASL(x);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________

%% ----------------------------------------------------------------------------------------
%% Administration

if nargin<2 || isempty(bASL) || ~isnumeric(bASL) || bASL<0 || bASL>1
    bASL = 1; % this tells script that timeseries contain subtractive/pair-wise data
end

if  bASL
    InputPath = x.P.Path_ASL4D;
else
    InputPath = x.P.Path_func_bold;
end

[Fpath, Ffile, Fext] = xASL_fileparts(InputPath);
rpfile = fullfile( Fpath, ['rp_' Ffile '.txt']);
rInputPath = fullfile( Fpath, ['r' Ffile Fext]);
InputPathJson = fullfile( Fpath, [Ffile '.json']);
rInputPathJson = fullfile( Fpath, ['r' Ffile '.json']);
matFile = fullfile(Fpath, [Ffile '.mat']);


%% Set defaults
exclusion = NaN;
PercExcl = NaN;
MinimumtValue = NaN;


%% Read basic image information
tempnii = xASL_io_ReadNifti(InputPath);
nFrames = double(tempnii.hdr.dim(5)); % Total number of frames
nFramesPerTE=nFrames/numel(unique(x.Q.EchoTime)); % Number of frames per unique TE. Note that for most cases, where we only have a single TE, nFramesPerTE==nFrames

minVoxelSize = double(min(tempnii.hdr.pixdim(2:4)));

%% Define motion correction options
% bMoCoPossible boolean states if it is possible to perform motion correction with the given data 
% x.asl.module.motionCorrection states if the motion correction is wanted by the user
if isfield(x.Q,'LookLocker') && x.Q.LookLocker
	bMoCoPossible = false;
	fprintf('%s\n',['Skipping motion correction for ' x.P.SubjectID '_' x.P.SessionID ' as Look-Locker correction is not implemented.']);
elseif x.modules.asl.bContainsSubtracted
	% Motion correction is disabled, potentially insufficient contrast
	bMoCoPossible = false;
	fprintf('%s\n',['Skipping motion correction for ' x.P.SubjectID '_' x.P.SessionID ' because it only has DeltaM volumes, which may contain insufficient contrast']);
elseif nFramesPerTE > 1
	% We only do motion correction when the number of frames is higher than one. For multi-TE, we count number of frames per TE.
	bMoCoPossible = true;
else
	bMoCoPossible = false;
	fprintf('%s\n',['Skipping motion correction for ' x.P.SubjectID '_' x.P.SessionID ' because it had only ' num2str(nFramesPerTE) ' 3D frames.']);
end

% Set flags for multi-TE and multi-PLD datasets

if x.Q.nUniqueInitial_PLD>1
	% Here, we consider only simple multiPLD for averaging across PLDs
    bMultiPLD = 1;

	% Note that TimeEncoded is a special case and needs to be treated separately. Therefore, we make sure this flag is initialized
	if ~isfield(x.modules.asl, 'bTimeEncoded')
		x.modules.asl.bTimeEncoded = false;
	end
else
    bMultiPLD = 0;
end

if x.Q.nUniqueEchoTime>1
    bMultiTE = 1;
else
    bMultiTE = 0;
end

%% File management
% Here we define the files created in this wrapper, and delete them if they
% pre-exist
pathSave_NDV = fullfile(x.D.MotionDir, ['motion_correction_NDV_' x.P.SubjectID '_' x.P.SessionID '.mat']);
jpgfile_Motion = fullfile(x.D.MotionDir, ['rp_' x.P.SubjectID '_' x.P.SessionID '_motion.jpg']);
jpgfile_ThresholdFree = fullfile( x.D.MotionDir,['rp_' x.P.SubjectID '_' x.P.SessionID '_threshold_free_spike_detection.jpg']);
jpgfile_MotionSorted = fullfile( x.D.MotionDir,['rp_' x.P.SubjectID '_' x.P.SessionID '_PWI_motion_sorted.jpg']);

xASL_delete(pathSave_NDV);
xASL_delete(jpgfile_Motion);
xASL_delete(jpgfile_ThresholdFree);
xASL_delete(jpgfile_MotionSorted);

if ~bMoCoPossible
    return; % no sense to run this function without motion correction
end

    
%% Manage spike removal options: threshold-free (bENABLE) & fixed threshold (bSpikeRemoval)
% bSpikeRemoval = boolean for removing spikes with a fixed threshold
% (OPTIONAL, DEFAULTED to false)
% SpikeRemovalAbsoluteThreshold = the fixed threshold (in mm; OPTIONAL, defaulted to 0)
%
% bENABLE = boolean for running ENABLE (OPTIONAL, defaulted to true)
% SpikeRemovalThreshold = ENABLE's relative threshold (t-stats, OPTIONAL, DEFAULT =
% 0.01)

% SpikeRemovalThreshold is an optional field, by default we use ENABLE
% So, here we check if it exists, but otherwise we default to disabling it

if bMultiPLD || bMultiTE
    % outlier exclusion is temporarily disabled for multiPLD (including Time-encoded) or multi-TE
    % as we are still developing this feature
    fprintf('%s\n', 'multi-PLD or multi-TE detected, disabling outlier exclusion, not yet implemented');
    bSpikeRemoval = false;
    bENABLE = false;

elseif ~isfield(x.modules.asl,'SpikeRemovalThreshold') && ~isfield(x.modules.asl,'SpikeRemovalAbsoluteThreshold')

    fprintf('%s\n','x.modules.asl.SpikeRemovalThreshold was not defined yet, default setting = 0.01 used');
    x.modules.asl.SpikeRemovalThreshold = 0.01; % default threshold, decreased this from 0.05 to 0.01,
    % since we want to remove Spikes, perhaps except very small spikes
    x.modules.asl.SpikeRemovalAbsoluteThreshold = 0; % disable by default
    bSpikeRemoval = false;
    bENABLE = true;

elseif ~isfield(x.modules.asl,'SpikeRemovalAbsoluteThreshold')
    x.modules.asl.SpikeRemovalAbsoluteThreshold = 0; % disable by default
    bSpikeRemoval = false;
    bENABLE = true;

elseif ~isnumeric(x.modules.asl.SpikeRemovalAbsoluteThreshold) || x.modules.asl.SpikeRemovalAbsoluteThreshold<0 || x.modules.asl.SpikeRemovalAbsoluteThreshold>1
    % If the field exists (as evidenced by the previous elseif) we check if
    % it has a correct value
    warning('Invalid x.modules.asl.SpikeRemovalAbsoluteThreshold');
    fprintf('%s\n', ['set to ' xASL_num2str(x.modules.asl.SpikeRemovalAbsoluteThreshold)]);
    fprintf('%s\n', 'Should be numeric and set between 0:1');
    bSpikeRemoval = true;
    bENABLE = false;    

else % Here it is clear that SpikeRemovalThreshold exists, and has a correct value
    % If this is the case, we disable ENABLE
	bSpikeRemoval = true;
    bENABLE = false;
end

if nFramesPerTE < 10 % Only execute ENABLE if we have at least 5 control-label pairs (==10 volumes)
	bENABLE = false;
end


%% Manage zig-zag in motion correction
if bMultiTE || x.modules.asl.bTimeEncoded
    % ZigZag are temporarily disabled for multiTE and TimeEncoded
    % as we are still developing this feature
	% Note that for standard multi-PLD, Zig-zag can be applied because standard multi-PLD still has controls and labels
    fprintf('%s\n', 'multi-PLD or multi-TE detected, disabling zig-zag motion estimation, not yet implemented');    
    bZigZag = false;

elseif bASL && nFramesPerTE > 2
    % we use zig-zag motion regression for ASL
    % Minimum number of frames for ZigZag is > 2 (1 control-label pair)
	bZigZag = true;
    
else
    % for non-ASL data (e.g., fMRI) we disable zig-zag, but we still keep ENABLE
	bZigZag = false;
end


if ~usejava('jvm')
    fprintf('%s\n', 'No JavaVM detected, skipping plotting motion & exclusion matrix');
end


%% Manage quality settings
switch x.settings.Quality
	case 1 % normal quality
		flags.quality = 1;
		flags.sep = minVoxelSize;
		flags.rtm = 1; % realign to mean
	case 0 % low quality for fast try-out
		flags.quality = 0.01;
		flags.rtm = 0; % disable realign to mean
		flags.sep = minVoxelSize*2;
end

flags.interp = 1;
flags.graphics = 0;


%% ----------------------------------------------------------------------------------------
%% 1. Estimate motion
fprintf('\nSPM motion estimation:\n');
fprintf('ExploreASL estimates motion and aligns based on the first TEs (in the case of multiTE)\n');
fprintf('ExploreASL estimates motion and aligns within PLD only (in the case of multiPLD (excluding Hadamard))\n');


% Issue warning if empty image
if max(max(max(max(tempnii.dat(:)))))==0 || numel(unique(tempnii.dat(:)))==1
	warning('Invalid input image, skipping');
	return;
end


% If previous realign parameters exist, delete them
xASL_delete(rpfile);

% Run motion correction for corresponding case
% Note that this is the adapted spm_realign, including zig-zag
% regression to account for ASL's potential control-label difference in
% average head position

V = spm_vol(InputPath); % load the path
Y = spm_read_vols(V); % Read the image
rp_all = zeros(nFrames, 6); % rp_all is what ends up in the rp*.txt sidecar, rp=realign parameters
mat_all = zeros(4, 4, size(Y,4)); % mat_all is what ends up in the ASL*.mat sidecar, containing the orientation matrices for each frame/volume


if bMultiTE
    fprintf('Multi-TE data detected: aligning based on shortest TEs only\n');
    % Handles Multi-TE dataset regardless of PLD
	% Registers all frames with the shortest TEs and then applies the same transformation to all the longer TEs
	% It assumes that the volume is sorted in the order of acquisition with blocks of increasing TEs
    
	% Indexes that flags the first/shortest TEs frames
    idx_minTE = find((x.Q.EchoTime == min(x.Q.EchoTime)));
	idx_maxTE =  find((x.Q.EchoTime == max(x.Q.EchoTime)));
    
	if length(idx_minTE) ~= length(idx_maxTE)
		% TEs should all have the same number of blocks (e.g., control-label repetitions and/or PLDs)
		error('Number of shortest TEs and longest TEs do not match, check if there are missing frames/volumes');
	end

	% Realigns only the flagged frames
    spm_realign(V(idx_minTE), flags, bZigZag); 
	% This affects the ASL4D.mat file, rp_ASL4D.txt file and MAT within the ASL4D.nii volume
    
    rp_temp = load(rpfile); % Load the rp file written by spm_realign

	V = spm_vol(InputPath); % Read the updated volumes
    % Loop through TE groups
    for idxTE = 1:length(idx_minTE)
        % Repeats the motion estimates extracted from the rp file for the rest of the TEs and writes to the corresponding rows in rp_all
        rp_all(idx_minTE(idxTE):idx_maxTE(idxTE), :) = repmat(rp_temp(idxTE,:), x.Q.nUniqueEchoTime, 1); 
    
		for idxAllTE = idx_minTE(idxTE):idx_maxTE(idxTE)
			mat_all(:, :, idxAllTE) = V(idx_minTE(idxTE)).mat;
		end
        
    end
    
elseif bMultiPLD && ~x.modules.asl.bTimeEncoded
    fprintf('MultiPLD detected, aligning within PLDs only\n');
    % Handles only Multi-PLD datasets - aligns only between the same PLDs
	% Motion correction across all PLDs is not really necessary as that can be done with the simple motion correction
	% Note that we handle normal mutli-PLD (not TimeEncoded), so there are still control and label images and we can thus do ZigZag

    for pld = x.Q.uniqueInitial_PLD(:)'
        idxSinglePLD = find(x.Q.Initial_PLD == pld); % Finds the same PLDs
        spm_realign(V(idxSinglePLD), flags, bZigZag);
        
        rp_temp = load(rpfile); % Load the rp file written by spm_realign
		V = spm_vol(InputPath); % Read the updated volumes

		% Loop through TE groups
		rp_all(idxSinglePLD, :) = rp_temp;

		for idxPLD = idxSinglePLD(:)'
			mat_all(:, :, idxPLD) = V(idxPLD).mat;
		end
    end
else
    fprintf('Standard SPM motion estimation\n');
    % Handles simple datasets
    spm_realign(V, flags, bZigZag);
end
 

% For these special cases, we need to save the updated transformation matrices
if bMultiTE || (bMultiPLD && ~x.modules.asl.bTimeEncoded)
	% Save the updated matrix TXT
	writematrix(rp_all, rpfile, 'delimiter', '\t'); 

	% Also generate the MAT file - so remove the MAT file first. The values are already stored in mat_all
	xASL_delete(matFile);

	% Save the updated volume, one by one
	for iVolume = 1:size(Y,4)
		Vt = V(iVolume);                
		Vt.n = [iVolume 1];
		Vt.mat = mat_all(:,:,iVolume);
		spm_write_vol(Vt, Y(:,:,:,iVolume)); % Save one 3D volume
	end
	mat = mat_all;
	save(matFile, 'mat'); % Save also the MAT file - NIfTI header and MAT-file contain the same information, but that's what normally happens after spm_realign
end


%% ----------------------------------------------------------------------------------------
%% 2. Calculate position and motion parameters
fprintf('%s\n','Calculate & plot position & motion parameters');

% Summarize real-world realign parameters into net displacement vector (NDV)
rp = load(rpfile, '-ascii'); % load the 3 translation and 3 rotation values
MeanRadius = 50; % typical distance center head to cerebral cortex (Power et al., NeuroImage 2012)
% PM: assess this from logical ASL EPI mask? This does influence the weighting of rotations compared to translations

if max(rp(:))==0
	warning('Something wrong with motion parameters, skipping');
	return;
end

% Calculate the mean displacements
FD{1}=rp; % position (absolute displacement)
FD{2} = diff(rp); % motion (relative displacement)
	
[NDV, median_NDV, mean_NDV, max_NDV, SD_NDV, MAD_NDV] = xASL_wrp_RealignASL_compute_NDV(FD, MeanRadius);

if bMultiTE
	% For MultiTE, recalculate the means and SD from the first echo only. But keep the vectors of all displacements as they were
	idx_minTE = find((x.Q.EchoTime == min(x.Q.EchoTime)));
    FD_firstTE{1} = rp(idx_minTE, :); % position (absolute displacement) for first TEs only
    FD_firstTE{2} = diff(rp(idx_minTE, :)); % motion (relative displacement) for first TEs only

	[~, median_NDV, mean_NDV, max_NDV, SD_NDV, MAD_NDV] = xASL_wrp_RealignASL_compute_NDV(FD_firstTE, MeanRadius);
end

%% ----------------------------------------------------------------------------------------
%% 3. Threshold-free spike definition (based on ENABLE, but with t-stats rather than the threshold p<0.05)

if bENABLE || bSpikeRemoval
    % Resample ASL image (apply motion estimation)
    xASL_adm_DeleteFilePair(rInputPath, 'json');

    matlabbatch{1}.spm.spatial.realign.write.data = {InputPath};
    matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
    
    xASL_Copy(InputPathJson, rInputPathJson);
    
    % Create a mask from the mean PWI
    [PWI, ~, PWI4D] = xASL_im_ASLSubtractionAveraging(x, rInputPath);
    xASL_io_SaveNifti(rInputPath, x.P.Path_mean_PWI_Clipped, PWI, 32, false);
    
    MaskIm = xASL_im_ClipExtremes(x.P.Path_mean_PWI_Clipped, 0.95, 0.7);
    MaskIm = MaskIm>min(MaskIm(:));
    xASL_delete(x.P.Path_mean_PWI_Clipped);
end


if ~bENABLE && ~bSpikeRemoval
    fprintf('%s\n', 'Skipping ENABLE');
elseif bENABLE
    % Sort motion of control-label pairs
    MotionTime = NDV{2}; % motion
    MotionTime = MotionTime(1:2:end-1)+MotionTime(2:2:end); % additive motion for each control-label pair
    MotionTime(:,2) = 1:length(MotionTime);
    MotionTimeSort = sortrows(MotionTime,1);
    
    
    fprintf('Running ENABLE:   ');
    tValue(1,1) = 0;
    SortIM = PWI4D(:,:,:,MotionTimeSort(:,2)); % Sort ASL-pairs by motion
    for iVolume = 1:size(SortIM,4)
        TempIm = SortIM(:,:,:,iVolume);
        SortMask(:,iVolume) = TempIm(MaskIm); % create data columns
    end
    
    for iVolume=2:size(SortIM,4)
        xASL_TrackProgress(iVolume,size(SortMask, 2));
        TempTimeSeries = SortMask(:, 1:iVolume);
        [~, ~, ~, stats] = xASL_stat_ttest(TempTimeSeries, 0, 0.05, 'both', 2);
        tValue(iVolume,1) = xASL_stat_MedianNan(stats.tstat(:));
    end
    fprintf('\n');
    
    INDEXn = round(0.5*length(tValue));
    OptimumV = max(tValue(INDEXn:end));
    mintValue = max(find(tValue(INDEXn:end)==OptimumV)+INDEXn-1);
    % max() is added here in coincidental case where there are 2 identical t-values, max() errs on the conservative side
    MinimumtValue = max(tValue(INDEXn:end));
    
    if tValue(end)>(1-x.modules.asl.SpikeRemovalThreshold)*MinimumtValue
        % only exclude frames if the optimal t-value
        % is more than x% higher than including all frames (default
        % x.modules.asl.SpikeRemovalThreshold= 0.01;
        mintValue = length(tValue);
    end
    
    mintValuePlot = zeros(1,length(tValue));
    mintValuePlot(mintValue+1:end)=min(tValue);
end

xASL_adm_DeleteFilePair(rInputPath, 'json'); % delete temporary image


%% ----------------------------------------------------------------------------------------
%% 4. Set volumes to exclude
if bSpikeRemoval
    fprintf('%s\n', ['Running spike removal with threshold ' xASL_num2str(x.modules.asl.SpikeRemovalAbsoluteThreshold) ' mm']);
    exclusion = (NDV{2}>x.modules.asl.SpikeRemovalAbsoluteThreshold)';
    % Merge control and label
    exclusionPairs = (exclusion(1:2:end-1) + exclusion(2:2:end))>0;
    % Recreate controls
    exclusion(1:2:end-1) = exclusionPairs;
    % Recreate labels
    exclusion(2:2:end) = exclusionPairs;

elseif bENABLE
    % Detect frames for exclusion
    if bASL % if ASL
        exclusion = zeros(1, length(MotionTimeSort)*2);
    else  % if no ASL (e.g. fMRI)
        exclusion = zeros(1, length(MotionTimeSort));
    end
    
    if  mintValue<length(tValue)
        for iFrame=mintValue+1:length(MotionTimeSort)
            % Exclude pair
            ExcludePair = MotionTimeSort(iFrame, 2);
            if  bASL % if ASL
                ExcludeFrames = [ExcludePair*2-1 ExcludePair*2];
            else % if no ASL (e.g. fMRI)
                ExcludeFrames = ExcludePair;
            end
            exclusion(ExcludeFrames) = 1;
        end
    end
end

%% ----------------------------------------------------------------------------------------
%% 5. Plot motion 

% Assigning titles to the plots depending on motion correction approach and plotting the Position and Motion plots
pTitle = ['Position plot of ' x.P.SubjectID '-' x.P.SessionID ' relative to first frame'];
mTitle = ['Motion plot of ' x.P.SubjectID '-' x.P.SessionID];

xASL_wrp_RealignASL_plot_motion(NDV, mean_NDV, bENABLE, bSpikeRemoval, exclusion, jpgfile_Motion, pTitle, mTitle);


%% ----------------------------------------------------------------------------------------
%% 6. Save ENABLE sorting plot
if bENABLE
    tValue(1:3) = tValue(4); % for nicer plotting

    if usejava('jvm') % only if JVM loaded
        fig = figure('Visible','off');
        plot([1:length(tValue)],tValue,'b',[1:length(tValue)],mintValuePlot * max(tValue),'r');
        xlabel('control-label pairs sorted by motion');
        ylabel('mean voxel-wise 1-sample t-test p-value');
        PercExcl    = round((sum(exclusion)/length(exclusion)*100)*10)/10;
        title(['Threshold free motion spike exclusion (red, ' num2str(PercExcl) '%) for ' x.P.SubjectID '_' x.P.SessionID]);
        
        fprintf('Saving motion plot to %s\n', jpgfile_ThresholdFree);
        
        xASL_adm_CreateDir(fileparts(jpgfile_ThresholdFree));
        saveas(fig, jpgfile_ThresholdFree, 'jpg');
        close all;
        clear fig;
    end
end

%% ----------------------------------------------------------------------------------------
%% 7. Save QC images before and after volume-spikes exclusion

if bSpikeRemoval || bENABLE
	Slice2Show = floor(size(PWI4D,3)*0.67); % e.g. slice 11/17
end

if bSpikeRemoval
    % display full timeseries without despiking
    ExampleIm_Full = xASL_im_rotate(xASL_stat_MeanNan(PWI4D(:,:,Slice2Show,:), 4), 90);
    % display timeseries with removing spike volumes
    ExampleIm_noSpikes = xASL_im_rotate(xASL_stat_MeanNan(PWI4D(:,:,Slice2Show,~exclusionPairs), 4), 90);
    TotalCheck = [ExampleIm_Full ExampleIm_noSpikes];

elseif bENABLE
    % Save 7 images, 3 before & 3 after exclusion
    IndexIs = [1 round(mintValue/3)  round(mintValue/2) mintValue];
    diffIndex = (length(tValue)-mintValue)/3;
    IndexIs(5:7) = [mintValue+diffIndex mintValue+2*diffIndex length(tValue)];
    IndexIs = round(IndexIs);
    
    % pre-allocation for more efficient memory usage
    ExampleIM = zeros(size(SortIM,1), size(SortIM,2), length(IndexIs));
    ExampleIM = single(ExampleIM);
    for iVolume=1:length(IndexIs)
        ExampleIM(:,:,iVolume) = xASL_stat_MeanNan(SortIM(:,:,Slice2Show,1:IndexIs(iVolume)), 4);
    end
    TotalCheck = xASL_vis_TileImages(xASL_im_rotate(ExampleIM,90), 4);
end
   
if bSpikeRemoval || bENABLE
    % Find intensities
    SortValues = sort(TotalCheck(isfinite(TotalCheck)));
    MinValue = SortValues(max(1,round(0*length(SortValues))));
    MaxValue = SortValues(round(0.975*length(SortValues)));
    
    TotalCheck(TotalCheck<MinValue) = MinValue;
    TotalCheck(TotalCheck>MaxValue) = MaxValue;
    
    fprintf('Saving motion plot to %s\n', jpgfile_MotionSorted);
    xASL_vis_Imwrite(TotalCheck, jpgfile_MotionSorted);
end


%% ----------------------------------------------------------------------------------------
%% 8. Remove spike volumes from NIfTI

if bENABLE || bSpikeRemoval

    if sum(exclusion)==0
        fprintf('No spike volumes detected for removal from NIfTI\n');
    elseif sum(exclusion)<0
        warning('Illegal exclusion matrix');
    else % only if spikes have been detected
        fprintf('Remove spike volumes from NIfTI\n');

        % Load nifti
        TempIm = xASL_io_Nifti2Im(InputPath);
        
        % Remove spikes
        % skip this for fMRI, which is more complicated
        % due to tissue T1 effects (incomplete saturation, so temporal relation between volumes)
        NewIm = TempIm(:,:,:, ~exclusion);
        
        % Do same for *.mat motion sidecars
        LoadParms = load(x.P.Path_ASL4D_mat, '-mat');
        mat = LoadParms.mat;
        mat = mat(:,:, ~exclusion);
        save(x.P.Path_despiked_ASL4D_mat, 'mat');

        % Also change parameter vectors
        jsonFields.Q.EchoTime = x.Q.EchoTime(~exclusion);
        jsonFields.Q.Initial_PLD = x.Q.Initial_PLD(~exclusion);
		if isfield(x.Q, 'LabelingDuration') && ~isempty(x.Q.LabelingDuration)
			jsonFields.Q.LabelingDuration = x.Q.LabelingDuration(~exclusion);
		end

        xASL_io_SaveNifti(x.P.Path_ASL4D, x.P.Path_despiked_ASL4D, NewIm, 32, 0, [], 1, jsonFields, true, [true true, false]);
    end
    
else
    exclusion = 0;
    PercExcl = 0;
    MinimumtValue = 0;
end


%% ----------------------------------------------------------------------------------------
%% 9. Save motion statistics before excluding motion spikes
% Save results for later summarization in analysis module
xASL_adm_CreateDir(x.D.MotionDir);
save(pathSave_NDV, 'NDV','median_NDV','mean_NDV','max_NDV','SD_NDV','MAD_NDV','exclusion','PercExcl','MinimumtValue');

%% ----------------------------------------------------------------------------------------
%% 10. Save motion statistics after excluding motion spikes
    
    if sum(exclusion)>0
        for ii=1:2
            NDV_SpikesRemoved{ii} = NDV{ii}(~exclusion);
            
            median_NDV_SpikesRemoved{ii} = median(NDV_SpikesRemoved{ii});
            mean_NDV_SpikesRemoved{ii} = mean(NDV_SpikesRemoved{ii});
            max_NDV_SpikesRemoved{ii} = max(NDV_SpikesRemoved{ii});
            SD_NDV_SpikesRemoved{ii} = std(NDV_SpikesRemoved{ii});
            MAD_NDV_SpikesRemoved{ii} = xASL_stat_MadNan(NDV_SpikesRemoved{ii},0); % median absolute deviation from median
        end
        
        save(pathSave_NDV, 'NDV','median_NDV','mean_NDV','max_NDV','SD_NDV','MAD_NDV','exclusion','PercExcl','MinimumtValue', ...
        'NDV_SpikesRemoved','median_NDV_SpikesRemoved','mean_NDV_SpikesRemoved','max_NDV_SpikesRemoved','SD_NDV_SpikesRemoved','MAD_NDV_SpikesRemoved');
    end

    
end

function [NDV, median_NDV, mean_NDV, max_NDV, SD_NDV, MAD_NDV] = xASL_wrp_RealignASL_compute_NDV(FD, MeanRadius)
    for ii = 1:2 % 1 = absolute displacement 2 = relative displacement==motion
	    tx{ii} = FD{ii}(:,1); ty{ii} = FD{ii}(:,2); tz{ii}  = FD{ii}(:,3); % translations
	    rx{ii} = FD{ii}(:,4); ry{ii} = FD{ii}(:,5); rz{ii}  = FD{ii}(:,6); % rotations (pitch, roll, yaw)
	    
	    PartTranslation{ii} = tx{ii}.^2 + ty{ii}.^2 + tz{ii}.^2;
	    PartRotation{ii} = 0.2*MeanRadius^2* ((cos(rx{ii})-1).^2 + (sin(rx{ii})).^2 + (cos(ry{ii})-1).^2 + (sin(ry{ii})).^2 + (cos(rz{ii})-1).^2 + (sin(rz{ii})).^2);
	    try
		    NDV{ii} = sqrt(PartTranslation{ii} + PartRotation{ii});
	    catch
		    
	    end
	    
	    if ii==2
		    NDV{2} = [0; NDV{2}]; % add leading zero difference
	    end
	    
	    % Descriptives
	    median_NDV{ii} = median(NDV{ii});
	    mean_NDV{ii} = mean(NDV{ii});
	    max_NDV{ii} = max(NDV{ii});
	    SD_NDV{ii} = std(NDV{ii});
	    MAD_NDV{ii} = xASL_stat_MadNan(NDV{ii},0); % median absolute deviation from median
    end
end

function fig = xASL_wrp_RealignASL_plot_motion(NDV, mean_NDV, bENABLE, bSpikeRemoval, exclusion, outFile, pTitle, mTitle)
    if usejava('jvm') % only if JVM loaded
        fig = figure('Visible','off');
        for FD_idx = 1:2 % 1 = absolute displacement 2 = relative displacement==motion
		    subplot(3,1,FD_idx); % plot position (subplot 1) & motion (subplot 2)
		    plot(NDV{FD_idx},'Color',[0.4,0.4,0.4]); % lines between frames
		    hold on
		    plot(NDV{FD_idx},'o','MarkerSize',5); % circles for frames
		    hold on
		    
		    plot(repmat(mean_NDV{FD_idx},length(NDV{FD_idx}),1),'Color',[0,0,1]); % mean NDV in blue
		    hold on

		    if FD_idx==1
			    title(pTitle);
			    ylabel('NDV (mm)');
			    
		    elseif FD_idx==2
			    title(mTitle);
			    ylabel('NDV/frame (mm//frame)');
		    end
		    
			xlabel('Frame #');
			axis([1 length(NDV{FD_idx}) 0 max(NDV{FD_idx})*1.05]); % axis fixing
        end
        
        if (bENABLE || bSpikeRemoval)
            subplot(3,1,3);
			hold on
            plot(exclusion,'r');
            ylabel('Exclusion matrix');
            axis([1 length(NDV{FD_idx}) 0 max(NDV{FD_idx})*1.05]);
        end
        
        fprintf('Saving motion plot to %s\n', outFile);
        
        xASL_adm_CreateDir(fileparts(outFile));
        saveas(fig, outFile, 'jpg');
        close (fig);
    end
end
