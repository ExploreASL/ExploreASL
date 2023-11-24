function [PWI, PWI3D, PWI4D, x, Control, Control3D, Control4D] = xASL_im_ASLSubtractionAveraging(x, path_ASL4D, PWI4D, PWI3D, Control4D, Control3D)
%xASL_im_ASLSubtractionAveraging Central function for ASL subtraction and averaging, for both PWI & controls
%
% FORMAT: [PWI, PWI3D, PWI4D, x] = xASL_im_ASLSubtractionAveraging(x, ASL4D [, PWI4D, PWI3D])
%
% INPUT:
%   x       - structure containing fields with all information required to run this function (REQUIRED)
%             with the following information:
%   x.modules.asl.bTimeEncoded - if we perform time decoding instead of
%             pair-wise subtraction (boolean, REQUIRED)
%             NB: see xASL_quant_HadamardDecoding for fields required by
%             time encoding
%   x.Q     - field containing quantification parameters (REQUIRED)
%   x.Q.EchoTime - numerical vector. If this is a scalar, it will be extended to the number of image volumes (REQUIRED)
%   x.Q.Initial_PLD - numerical vector. If this is a scalar, it will be extended to the number of image volumes (REQUIRED)
%   x.Q.LabelingDuration - numerical vector. If this is a scalar, it will be extended to the number of image volumes (REQUIRED)
% 
%   path_ASL4D - Path to the encoded ASL4D image we want to decode (STRING OR 4D MATRIX, REQUIRED)
%             Alternatively, this can contain the image matrix
%             (xASL_io_Nifti2Im below allows both path and image matrix inputs).
%             Best if this is a path, so any motion correction can be applied.
%   PWI4D   - likewise, see explanation below. If this is provided, it will override the PWI4D calculation below (STRING OR 4D MATRIX, OPTIONAL,
%               DEFAULT = generated below from ASL4D)
%   PWI3D   - likewise, see explanation below. If this is provided, it will override the PWI4D calculation below (STRING OR 4D MATRIX, OPTIONAL,
%               DEFAULT = generated below from PWI4D)
%
% OUTPUT:
%   PWI     - 3D image matrix, see explanation below
%   PWI3D   - 4D image matrix, see explanation below
%   PWI4D   - 4D image matrix, see explanation below
%   x       - updated structure, with the following quantification parameters:
%             fitting to the nomenclature explained below
%   Control - same as PWI but then non-subtracted controls
%   Control3D - same as PWI3D but then non-subtracted controls
%   Control4D - same as PWI4D but then non-subtracted controls
%
%   x.Q.EchoTime_PWI4D
%   x.Q.InitialPLD_PWI4D
%   x.Q.LabelingDuration_PWI4D
%
%   x.Q.EchoTime_PWI3D
%   x.Q.InitialPLD_PWI3D
%   x.Q.LabelingDuration_PWI3D
%
%   x.Q.EchoTime_PWI
%   x.Q.InitialPLD_PWI
%   x.Q.LabelingDuration_PWI
%
% DESCRIPTION: This function performs subtraction and averaging all across the ExploreASL ASL modules, such as xASL_wrp_Realign,
%              and xASL_wrp_RegisterASL to compute temporary PWI images, xASL_wrp_ResampleASL for creating QC images and preparing
%              quantification, and xASL_wrp_Quantify if PWI4D is not used by the quantification algorithm (such as BASIL).
%              It takes any PWI and control phase as input and will output any phase that it computed.
%              The phases are: ASL4D->PWI4D    ->PWI3D    ->PWI
%                         and: ASL4D->Control4D->Control3D->Control
%
% This function performs the following steps:
% 1. Admin & checks
% 2. Apply motion correction if this is needed
%              
% Nomenclature:
% ASL4D = original timeseries
% PWI4D = subtracted timeseries -> for quantification purposes (also applies to control4D)
% PWI3D = averaged per TE—PLD—LabDur combination -> for QC purposes  (also applies to control3D)
% PWI   = single perfusion-weighted volume -> for registration purposes  (also applies to control)
% 
% % we use the same nomenclature for control: control4D, control3D, control, as these are the same control-label pairs but without subtraction
%
% The quantification submodule xASL_wrp_Quantify only uses PWI4D from here, and recalculates PWI3D and PWI for quantification, potentially improving
% quantification. So the user can, after QC, potentially remove some volumes from PWI4D, but should be aware that the information in the JSON sidecar
% (TE, PLD, LabDur) should be adapted accordingly.
%
% EXAMPLE used by xASL_wrp_RegisterASL: PWI = xASL_im_ASLSubtractionAveraging(x, ASL_im);
% EXAMPLE used by xASL_wrp_ResampleASL: [PWI, PWI3D, PWI4D, x] = xASL_im_ASLSubtractionAveraging(x, ASL4D, PWI4D, PWI3D);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


%% ========================================================================================
%% 1. Admin & checks
% Input image volumes should be correct
% Creation of the PWI3D or PWI4D image volumes can be skipped if they are provided as input

if nargin<1 || isempty(x)
    error('x input missing');
end

if nargin<2 || isempty(path_ASL4D)
    error('path_ASL4D input missing');
elseif numel(path_ASL4D)~=1 && ~isnumeric(path_ASL4D) && numel(path_ASL4D)<1000
    % this is a path, not an image
else
    error('path_ASL4D cannot be an image matrix, needs to be a path');
end

% By default, we create 4D stuff, unless they are provided as input
bCreatePWI4D = true;
bCreateControl4D = true;

% The 3D stuff is not created by default, because we need PWI4D OR Control4D for this
bCreatePWI3D = false;
bCreateControl3D = false;

% The PWI/Control are not created by default, because we need PWI3D OR Control3D for this
bCreatePWI = false;
bCreateControl = false;

% We skip the computation of PWI4D if it is provided
if nargin>=3 && ~isempty(PWI4D)
    if ~isnumeric(PWI4D)
        error('Illegal PWI4D, it should be numerical');
    elseif ndims(PWI4D)>4 || ndims(PWI4D)<3 % In theory, PWI4D can be a 3D matrix when it contains a single deltaM
        error('PWI4D has an incorrect number of dimensions');
	else
        bCreatePWI4D = false; % All is correct and we don't have to compute PWI4D
        bCreatePWI3D = true; % if PWI4D is provided, we compute PWI3D
    end
else
    PWI4D = []; % dummy output
end

% We skip the computation of Control4D if it is provided
if nargin>=5 && ~isempty(Control4D) % same as above
    if ~isnumeric(Control4D)
        error('Illegal Control4D, it should be numerical');
    elseif ndims(Control4D)>4 || ndims(Control4D)<3
        error('Control4D has an incorrect number of dimensions');
	else
        bCreateControl4D = false; % All is correct and we don't have to compute Control4D
        bCreateControl3D = true; % if Control4D is provided, we compute Control3D
    end
else
    Control4D = []; % dummy output  
end

% If PWI4D is provided but Control4D is not, we only want to compute PWI3D and PWI
if bCreatePWI3D && ~bCreateControl3D
    bCreateControl4D = false; % If only PWI4D is provided, we only want to compute PWI3D and PWI
end
% Same vice versa
if bCreateControl3D && ~bCreatePWI3D
    bCreatePWI4D = false; % If only Control4D is provided, we only want to compute Control3D and Control
end


% We skip the computation of PWI3D if it is provided
if nargin>=4 && ~isempty(PWI3D)
    if ~isnumeric(PWI3D)
        error('Illegal PWI3D, it should be numerical');
    elseif ndims(PWI3D)>4 || ndims(PWI3D)<3
        error('PWI3D has an incorrect number of dimensions');
	else
        bCreatePWI3D = false;
        bCreatePWI = true; % If PWI3D is provided, we compute PWI
    end
else
    PWI3D = []; % dummy output
end

% We skip the computation of Control3D if it is provided
if nargin>=6 && ~isempty(Control3D)
    if ~isnumeric(Control3D)
        error('Illegal Control3D, it should be numerical');
    elseif ndims(Control3D)>4 || ndims(Control3D)<3
        error('Control3D has an incorrect number of dimensions');
	else
        bCreateControl3D = false;
        bCreateControl = true; % If Control3D is provided, we compute PWI
    end
else
    Control3D = []; % dummy output
end

PWI = []; % dummy output
Control = []; % dummy output



%% ========================================================================================
%% 2. PWI4D
if bCreatePWI4D || bCreateControl4D

    % Load ASL4D.nii
    ASL4Dnii = xASL_io_ReadNifti(path_ASL4D);
    sizeASL4D = size(ASL4Dnii.dat);
    nVolumes = sizeASL4D(4);
    nDims = length(sizeASL4D);
    
    if nDims>4
        fprintf('\n\nIn BIDS ASL NIfTIs should have [X Y Z n/PLD/TE/LD/etc] as 4 dimensions\n\n');
        error('ASL4D.nii contains more than 4 dimensions');
        
    elseif nVolumes>1 && mod(nVolumes, 2) ~= 0 && ~x.modules.asl.bContainsDeltaM
        warning('Odd (i.e., not even) number of control-label pairs, skipping');
        % mod(nVolumes, 2) ~= 0 here means that nVolumes is not divisible by 2
        % unsubtracted control-label pairs should result in an even number of volumes
        return;
    end


    %% 2A. Apply motion correction if needed
    % At the beginning of the ASL module we estimate motion, and optionally remove outliers/spikes
    % This can already be applied, but in some cases this is not yet applied and by just loading ASL4D.nii
    % this is then forgotten. So here, we check if a .mat (e.g., ASL4D.mat) sidecar exists, and then run
    % resampling to apply motion correction. This creates a temporary r-prefixed file (e.g., rASL4D.nii)
    % which is deleted right after we reload its NIfTI as image volumes.
    
    [Fpath, Ffile] = xASL_fileparts(path_ASL4D);
    sideCarMat = fullfile(Fpath, [Ffile '.mat']);

    if exist(sideCarMat, 'file')
        % we check if the mat-sidecar exists, with the motion estimation parameters
        % If this file exists, we assume that the original ASL4D input file has not been
        % interpolated yet (and the motion estimation has not been applied yet)

        % Apply motion correction & resample
        if nVolumes>1
            for iS=1:nVolumes
                matlabbatch{1}.spm.spatial.realign.write.data{iS,1} = [path_ASL4D ',' num2str(iS)];
            end
        
            matlabbatch{1}.spm.spatial.realign.write.roptions.which     = [2 0];
            matlabbatch{1}.spm.spatial.realign.write.roptions.interp    = 4;
            matlabbatch{1}.spm.spatial.realign.write.roptions.wrap      = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.write.roptions.mask      = 1;
            matlabbatch{1}.spm.spatial.realign.write.roptions.prefix    = 'r';
        
            spm_jobman('run',matlabbatch); % this applies the motion correction in native space
            
            path_rASL4D = fullfile(Fpath, ['r' Ffile '.nii']);
            
            % RELOAD ASL4D.nii
            % Note that we only reload ASL4D.nii here if we have motion corrected timeseries.
            % If path_ASL4D was not a path, or if sideCarMat didn't exist (motion estimation was not performed)
            % or if nVolumes==1 (we did not have multiple ASL volumes) there is no reason to resample and we keep using
            % the already above loaded ASL4D.nii
            ASL4D = xASL_io_Nifti2Im(path_rASL4D);
            xASL_delete(path_rASL4D);
        else
            ASL4D = xASL_io_Nifti2Im(path_ASL4D);
        end
    end


    %% 2B. Apply time decoding if needed
    if x.modules.asl.bTimeEncoded
        %% #997 Check that the motion correction of Hadamard is correctly performed here,
        % should we use path_rASL4D from here onwards, which defaults to path_ASL4D?

        % Decoding of TimeEncoded data - it outputs the decoded image as well as the parameters
        [PWI4D, Control4D, x.Q] = xASL_quant_HadamardDecoding(ASL4D, x.Q);
        %% #997 JAN here we have to load ASL4D which is resampled for motion correction in 2A
	    nVolumes = size(PWI4D, 4);

        bCreatePWI3D = true; % because we can do this now
        bCreateControl3D = true; % because we can do this now        


    %% 2C. Non-time encoded subtraction
    else
        if nVolumes>1 && ~x.modules.asl.bContainsDeltaM
            % Paired subtraction
            [Control4D, Label4D] = xASL_quant_GetControlLabelOrder(ASL4D);
            PWI4D = Control4D - Label4D;
            
            % Skip every other value in the vectors as they were stored for both control and label images 
            x.Q.EchoTime_PWI4D = x.Q.EchoTime_PWI4D(1:2:end);
            x.Q.InitialPLD_PWI4D = x.Q.InitialPLD_PWI4D(1:2:end);
            x.Q.LabelingDuration_PWI4D = x.Q.LabelingDuration_PWI4D(1:2:end);
            % In this case, the vectors are equal for PWI4D and Control4D
            x.Q.EchoTime_Control4D = x.Q.EchoTime_PWI4D;
            x.Q.InitialPLD_Control4D = x.Q.InitialPLD_PWI4D;
            x.Q.LabelingDuration_Control4D = x.Q.LabelingDuration_PWI4D;

            bCreatePWI3D = true; % because we can do this now
            bCreateControl3D = true; % because we can do this now
        elseif ~x.modules.asl.bContainsDeltaM
            warning('Single volume detected that was not subtracted, cannot create PWI images');
            PWI4D = [];
            Control4D = ASL4D;
            bCreateControl3D = true; % because we can do this now
            % but we cannot create PWI4D, PWI3D, PWI

        else % the same but then without subtraction
            % either if we have only a single volume
            % or if we have multiple volumes that are already subtracted (bContainsDeltaM)
            % in both cases we can create PWI4D, PWI3D, PWI, but not Control4D, etc
            PWI4D = ASL4D;
            Control4D = [];

            % Now the vectors stay the same as for ASL4D (there is no subtraction)
            x.Q.EchoTime_PWI4D = x.Q.EchoTime;
            x.Q.InitialPLD_PWI4D = x.Q.Initial_PLD;
            x.Q.LabelingDuration_PWI4D = x.Q.LabelingDuration;
            x.Q.EchoTime_Control4D = x.Q.EchoTime;
            x.Q.InitialPLD_Control4D = x.Q.InitialPLD;
            x.Q.LabelingDuration_Control4D = x.Q.LabelingDuration;

            bCreatePWI3D = true; % because we can do this now
        end
    end
end

%% ========================================================================================
%% 5. Create PWI3D AND/OR Control3D
if bCreatePWI3D
    
    % After averaging across PLDs, we'll obtain these unique PLDs+LD combinations
    % indexAverage_PLD_LabDur lists for each original position to where it should be averaged
    [~, ~, iUnique_TE_PLD_LabDur_PWI4D] = unique([x.Q.EchoTime_PWI4D(:), x.Q.InitialPLD_PWI4D(:), x.Q.LabelingDuration_PWI4D(:)], 'stable', 'rows');
    
    % MultiPLD-multiLabDur PWI3D after averaging
    for iTE_PLD_LabDur = 1:max(iUnique_TE_PLD_LabDur_PWI4D)
        indicesAre = iUnique_TE_PLD_LabDur_PWI4D == iTE_PLD_LabDur;
        PWI3D(:, :, :, iTE_PLD_LabDur) = xASL_stat_MeanNan(PWI4D(:, :, :, indicesAre), 4); % Averaged PWI4D

        x.Q.EchoTime_PWI3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.EchoTime_PWI4D(indicesAre));
        x.Q.InitialPLD_PWI3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.InitialPLD_PWI4D(indicesAre));
        x.Q.LabelingDuration_PWI3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.LabelingDuration_PWI4D(indicesAre));
    end

    bCreatePWI = true; % because we can do this now
end

if bCreateControl3D
    % Now we do the same for Control4D
    [~, ~, iUnique_TE_PLD_LabDur_Control4D] = unique([x.Q.EchoTime_Control4D(:), x.Q.InitialPLD_Control4D(:), x.Q.LabelingDuration_Control4D(:)], 'stable', 'rows');
    
    % MultiPLD-multiLabDur PWI3D after averaging
    for iTE_PLD_LabDur = 1:max(iUnique_TE_PLD_LabDur_Control4D)
        indicesAre = iUnique_TE_PLD_LabDur_Control4D == iTE_PLD_LabDur;
        PWI3D(:, :, :, iTE_PLD_LabDur) = xASL_stat_MeanNan(Control4D(:, :, :, indicesAre), 4); % Averaged Control4D
        Control3D(:, :, :, iTE_PLD_LabDur) = xASL_stat_MeanNan(Control4D(:, :, :, indicesAre), 4); % Averaged control4D (same indices as PWI4D)

        x.Q.EchoTime_Control3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.EchoTime_Control4D(indicesAre));
        x.Q.InitialPLD_Control3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.InitialPLD_Control4D(indicesAre));
        x.Q.LabelingDuration_Control3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.LabelingDuration_Control4D(indicesAre));
    end

    bCreateControl = true; % because we can do this now
end


%% 6. Create PWI AND/OR Control
% We create a dummy CBF image for registration purposes
% The earliest echo, the latest PLD and the longest labeling duration
% are the best for this, having most SNR, CBF-weighting, and SNR,
% respectively. 
% But from a single volume, we have little signal.
% So we take a weighted average accordingly

%% This would be taking a single volume only
% % The earliest echo has the most SNR for perfusion-weighting
% iMinTE = x.Q.EchoTime_PWI3D == min(x.Q.EchoTime_PWI3D(:));
% % The latest PLD has the most perfusion-weighting
% iMaxPLD = x.Q.InitialPLD_PWI3D == max(x.Q.InitialPLD_PWI3D(:));
% % The longest labeling duration has the most SNR and most perfusion-weighting
% iMaxLabelingDuration = x.Q.LabelingDuration_PWI3D == max(x.Q.LabelingDuration_PWI3D(:));
% % We take the index that has all these        
% i3D = iMinTE & iMaxPLD & iMaxLabelingDuration;
% if isempty(i3D) || sum(i3D)~=1
%     error('Illegal index for PWI image'); % THIS ERROR SHOULD NOT BE GIVEN TO THE USER
% end
% PWI = PWI3D(:, :, :, i3D);
% x.Q.EchoTime_PWI = min(x.Q.EchoTime_PWI3D(:));
% x.Q.initialPLD_PWI = max(x.Q.InitialPLD_PWI3D(:));
% x.Q.LabelingDuration_PWI = max(x.Q.LabelingDuration_PWI3D(:));

%% PM: do a weighted average here?
%
%% Get mean control does it this way, choosing the PLD closest to 2000 ms
    % 
	% % Get unique PLDs
	% idealPLD = unique(Initial_PLD);
    % 
	% % Find the index of the one closest to 2000 ms
	% [~, iPLD] = min(abs(idealPLD-2000));
    % 
	% % Pick up the ideal PLD as the one closest to 2000 ms
	% idealPLD = idealPLD(iPLD(1));
    % 
	% if (isfield(x.Q,'LookLocker') && x.Q.LookLocker) || x.modules.asl.bContainsDeltaM
	% 	% For Look-Locker, get the middle one
	% 	idealPLD = idealPLD(round(numel(idealPLD)/2));
	% else
	% 	% For normal multi-PLD, get the latest PLD
	% 	idealPLD = idealPLD(end);
	% end
    % 
	% % Find all dynamics with that PLD
	% imMeanControl = imMeanControl(:,:,:,Initial_PLD==idealPLD);

if bCreatePWI
    %% Current quick&dirty fix was the average of all volumes with the lowest echo time
    iMinTE = x.Q.EchoTime_PWI3D == min(x.Q.EchoTime_PWI3D(:));
    PWI = xASL_stat_MeanNan(PWI3D(:, :, :, iMinTE), 4);
    x.Q.EchoTime_PWI = min(x.Q.EchoTime_PWI3D(:));
    x.Q.initialPLD_PWI = mean(x.Q.InitialPLD_PWI3D(:)); %% Jan is this correct?
    x.Q.LabelingDuration_PWI = mean(x.Q.LabelingDuration_PWI3D(:)); %% Jan is this correct?
end

if bCreateControl
    % The same for Control
    iMinTE = x.Q.EchoTime_Control3D == min(x.Q.EchoTime_Control3D(:));
    Control = xASL_stat_MeanNan(Control3D(:, :, :, iMinTE), 4);
    x.Q.EchoTime_Control = min(x.Q.EchoTime_Control3D(:));
    x.Q.initialPLD_Control = mean(x.Q.InitialPLD_Control3D(:)); %% Jan is this correct?
    x.Q.LabelingDuration_Control = mean(x.Q.LabelingDuration_Control3D(:)); %% Jan is this correct?
end

end