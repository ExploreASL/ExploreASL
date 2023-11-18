function [PWI, PWI3D, PWI4D, x, Control, Control3D, Control4D] = xASL_im_ASLSubtractionAveraging(x, path_ASL4D, PWI4D, PWI3D)
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
bCreatePWI3D = true;
bCreatePWI4D = true;


% By default, we calculate PWI3D, the calculation of PWI3D is only skipped if a relevant non-empty image is provided
if nargin>=4 && ~isempty(PWI3D) % If PWI3D parameter is not provided or is missing, then we skip to its calculation
    if ~isnumeric(PWI3D) % If PWI3D is provided, it has to be numeric
        error('Illegal PWI3D, it should be numerical');
    elseif ndims(PWI3D)>4 || ndims(PWI3D)<3 % And have the correct dimensions
        error('PWI3D has an incorrect number of dimensions');
	else
        bCreatePWI3D = false; % All is correct and we can use the provided PWI3D
    end
end

% Similar checks are done for PWI4D
if nargin>=3 && ~isempty(PWI4D)
    if ~isnumeric(PWI4D)
        error('Illegal PWI4D, it should be numerical');
    elseif ndims(PWI4D)>4 || ndims(PWI4D)<3 % In theory, PWI4D can be a 3D matrix when it contains a single deltaM
        error('PWI4D has an incorrect number of dimensions');
	else
        bCreatePWI4D = false; % All is correct and we can use the provided PWI4D
    end
end


%% ========================================================================================
%% 2. PWI4D
if bCreatePWI4D

    if nargin<2 || isempty(path_ASL4D)
        error('path_ASL4D input missing');
    else
        % LOAD ASL4D.nii
        ASL4D = xASL_io_Nifti2Im(path_ASL4D); % Load time-series nifti from path if input was a path
        % (input can also be an image matrix)
        
        % here we perform several dimensional checks before we can continue

        nVolumes = size(ASL4D, 4);
        
        if ndims(ASL4D)>4
            fprintf('\n\nIn BIDS ASL NIfTIs should have [X Y Z n/PLD/TE/LD/etc] as 4 dimensions\n\n');
            error('ASL4D.nii contains more than 4 dimensions');
            
        elseif nVolumes>1 && mod(nVolumes, 2) ~= 0
            warning('Odd (i.e., not even) number of control-label pairs, skipping');
            % mod(nVolumes, 2) ~= 0 here means that nVolumes is not divisible by 2
            return;
        end
    end
    if isempty(x)
        error('x input missing');
    end



    %% 2A. Apply motion correction if needed
    % At the beginning of the ASL module we estimate motion, and optionally remove outliers/spikes
    % This can already be applied, but in some cases this is not yet applied and by just loading ASL4D.nii
    % this is then forgotten. So here, we check if a .mat (e.g., ASL4D.mat) sidecar exists, and then run
    % resampling to apply motion correction. This creates a temporary r-prefixed file (e.g., rASL4D.nii)
    % which is deleted right after we reload its NIfTI as image volumes.
    
    if numel(path_ASL4D)~=1 && ~isnumeric(path_ASL4D) && numel(path_ASL4D)<1000
        % here we detect a path
        [Fpath, Ffile] = xASL_fileparts(path_ASL4D);
        sideCarMat = fullfile(Fpath, [Ffile '.mat']);
    
        if exist(sideCarMat, 'file')
            % we check if the mat-sidecar exists, with the motion estimation parameters
            % If this file exists, we assume that the original ASL4D input file has not been
            % interpolated yet (and the motion estimation has not been applied yet)
    
            % Apple motion correction & resample
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
            end
        end
    end


    %% 2B. Apply time decoding if needed
    if x.modules.asl.bTimeEncoded
        %% #997 Check that the motion correction of Hadamard is correctly performed here,
        % should we use path_rASL4D from here onwards, which defaults to path_ASL4D?

        % Decoding of TimeEncoded data - it outputs the decoded image as well as the parameters
        [PWI4D, Control4D] = xASL_quant_HadamardDecoding(path_ASL4D, x.Q);
	    nVolumes = size(PWI4D, 4);


    %% 2C. Non-time encoded subtraction
    else
        if nVolumes>1 && ~x.modules.asl.bContainsDeltaM
            % Paired subtraction
            [Control4D, Label4D] = xASL_quant_GetControlLabelOrder(ASL4D);
            PWI4D = Control4D - Label4D;
            
            % Skip every other value in the vectors as they were stored for both control and label images 
            x.Q.EchoTime_PWI4D = x.Q.EchoTime(1:2:end);
            x.Q.InitialPLD_PWI4D = x.Q.Initial_PLD(1:2:end);
            x.Q.LabelingDuration_PWI4D = x.Q.LabelingDuration(1:2:end);
        else % the same but then without subtraction
            PWI4D = ASL4D;
            Control4D = NaN;
        end
    end


%% ========================================================================================
%% 5. Create PWI3D
if bCreatePWI3D
    
    % After averaging across PLDs, we'll obtain these unique PLDs+LD combinations
    % indexAverage_PLD_LabDur lists for each original position to where it should be averaged
    [~, ~, iUnique_TE_PLD_LabDur] = unique([x.Q.EchoTime_PWI4D(:), x.Q.InitialPLD_PWI4D(:), x.Q.LabelingDuration_PWI4D(:)], 'stable', 'rows');
    
    % MultiPLD-multiLabDur PWI3D after averaging
    for iTE_PLD_LabDur = 1:max(iUnique_TE_PLD_LabDur)
        indicesAre = iUnique_TE_PLD_LabDur == iTE_PLD_LabDur;
        PWI3D(:, :, :, iTE_PLD_LabDur) = xASL_stat_MeanNan(PWI4D(:, :, :, indicesAre), 4); % Averaged PWI4D
        Control3D(:, :, :, iTE_PLD_LabDur) = xASL_stat_MeanNan(Control4D(:, :, :, indicesAre), 4); % Averaged control4D (same indices are PWI4D)

        x.Q.EchoTime_PWI3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.EchoTime_PWI4D(indicesAre));
        x.Q.InitialPLD_PWI3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.InitialPLD_PWI4D(indicesAre));
        x.Q.LabelingDuration_PWI3D(iTE_PLD_LabDur) = xASL_stat_MeanNan(x.Q.LabelingDuration_PWI4D(indicesAre));
    end
end


%% 6. Create PWI
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
%     error('Illegal index for PWI image');
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



%% Current quick&dirty fix was the average of all volumes with the lowest echo time
iMinTE = x.Q.EchoTime_PWI3D == min(x.Q.EchoTime_PWI3D(:));

PWI = xASL_stat_MeanNan(PWI3D(:, :, :, iMinTE), 4);
Control = xASL_stat_MeanNan(Control3D(:, :, :, iMinTE), 4);

x.Q.EchoTime_PWI = min(x.Q.EchoTime_PWI3D(:));
x.Q.initialPLD_PWI = mean(x.Q.InitialPLD_PWI3D(:)); %% Jan is this correct?
x.Q.LabelingDuration_PWI = mean(x.Q.LabelingDuration_PWI3D(:)); %% Jan is this correct?


end