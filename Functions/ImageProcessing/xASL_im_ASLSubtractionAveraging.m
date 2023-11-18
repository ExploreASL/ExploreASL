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

if bCreatePWI4D

    if nargin<2 || isempty(path_ASL4D)
        error('path_ASL4D input missing');
    else
        ASL4D = xASL_io_Nifti2Im(path_ASL4D); % Load time-series nifti from path if input was a path
        % (input can also be an image matrix)
        nVolumes = size(ASL4D, 4);
        
        
        % =====================================================================
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



    %% 2. Apply motion correction if this is needed
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
    
            %% ------------------------------------------------------------------------------------------
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
                % Reload motion corrected timeseries:
                ASL4D = xASL_io_Nifti2Im(path_rASL4D);
                xASL_delete(path_rASL4D);
            end
        end
    end
end
    
% =====================================================================
% -> THIS WILL BE MOVED TO xASL_quant_HadamardDecoding
%% 3. Time-Encoded subtraction
% xASL_quant_HadamardDecoding should not only decode the image volumes but also the quantification parameters
% that accompany the different volumes, i.e. EchoTime, PLD, and LabelingDuration.
% When this is done appropriately, PWI derivative volumes are easily created based on these parameters

if x.modules.asl.bTimeEncoded
    
    if bCreatePWI4D
        %% 1) Create PWI4D
        % Decoding of TimeEncoded data - it outputs a decoded image and it also saves as a NII
        [PWI4D, Control4D] = xASL_quant_HadamardDecoding(path_ASL4D, x.Q);
	    nVolumes = size(PWI4D, 4);
    end

    if bCreatePWI3D

    %% 2) Create PWI3D
	% Calculate Hadamard block size (number of unique PLDs.*labeling durations.* echo times) = number of volumes per repetition
    

    %% 1) x.Q.EchoTime vectors (identical to below)
    % First we deal with a too long vector (e.g., cutting off control volumes for Hadamard)
    nTE = length(x.Q.EchoTime);
    if nTE>nVolumes
        x.Q.EchoTime_PWI4D = x.Q.EchoTime(1:nVolumes); % can we do this based on the ASL4Dcontext.tsv ?
    else
        x.Q.EchoTime_PWI4D = x.Q.EchoTime;
    end
    nTE = length(x.Q.EchoTime_PWI4D); % Obtain the number of echo times

    
    nUniqueTE = length(uniqueTE); % Obtain the number of echo times
    % instead of nUniqueTE which we don't use anymore


    % % 
    % % % We do this here now for each sequence, but we could also do this
    % % % at the start of the ASL module
    % % 
    % % iUniqueTE = 1:nTEs;
    % % iUniqueTE = repmat(iUniqueTE, [1 nVolumes/nTEs]);


    
    %% 2) PLD vectors
    unique_InitialPLD = unique(x.Q.Initial_PLD);

    % with time encoded, we always skip the latest PLD
    % because that is a dummy PLD (it is in reality the control
    % scan)
    unique_InitialPLD = unique_InitialPLD(1:end-1);
    nUnique_InitialPLD = length(unique_InitialPLD);

    for iPLD=1:nUnique_InitialPLD
        startIndex = (iPLD-1).*nUniqueTE+1;
        endIndex = iPLD.*nUniqueTE;
        iUnique_InitialPLD(startIndex:endIndex) = iPLD;
    end
    iUnique_InitialPLD = repmat(iUnique_InitialPLD, [1 nVolumes/length(iUnique_InitialPLD)]);
    x.Q.InitialPLD_PWI4D = iUnique_InitialPLD;

    % %% 2) PLD vectors (identical to below)
    %% ->> THIS is what it should be
    % % First we deal with a too long vector (e.g., cutting off control volumes for Hadamard)
    % nPLD = length(x.Q.Initial_PLD);
    % if nPLD>nVolumes
    %     x.Q.InitialPLD_PWI4D = x.Q.Initial_PLD(end-nVolumes+1:end); % can we do this based on the ASL4Dcontext.tsv ?
    % else
    %     x.Q.InitialPLD_PWI4D = x.Q.Initial_PLD;
    % end
    % nPLD = length(x.Q.InitialPLD_PWI4D);
    % 


    

    %% 3) Labeling duration (LD) vectors (identical to below)
    % First we deal with a too long vector (e.g., cutting off control volumes for Hadamard)
    nLabDur = length(x.Q.LabelingDuration);
    if nLabDur>nVolumes
        x.Q.LabelingDuration_PWI4D = x.Q.LabelingDuration(1:nVolumes); % can we do this based on the ASL4Dcontext.tsv ?
    else
        x.Q.LabelingDuration_PWI4D = x.Q.LabelingDuration;
    end
    nLabDur = length(x.Q.LabelingDuration_PWI4D);


    end


    %% TIME DECODING WARNING -> MOVE TO xASL_quant_HadamardDecoding (check the original code)
    %%  BECAUSE TIME DECODING WILL ALSO ADAPT THE QUANTIFICATION PARAMETERS
    %%  ORDER
    % %% This part needs to be synced with the below part still
    % %% That we give equal warnings for all ASL flavor, if the number of TE/PLD/labdurs 
    % %% don't fit within the number of volumes
    % %% unless we already give this in the import, then we can remove it here
    % 
    % % nVolumesPerRepetition = nTEs .* (nUnique_InitialPLD - (nUnique_InitialPLD./x.Q.TimeEncodedMatrixSize));
    % nVolumesPerRepetition = nTEs .* nUnique_InitialPLD;
    % 
    % % nVolumesPerRepetition is also called blockSize by some
	% nRepetitions = nVolumes ./ nVolumesPerRepetition;
    % 
    % % First check if the number of volumes fits with the number of
    % % expected volumes per repetition (for the amount of x.Q.EchoTimes, PLDs, and labeling durations)
    % if nRepetitions~=round(nRepetitions)
    %     error(['nVolumes ' xASL_num2str(size(PWI4D,4)) ' cannot be composed of ' xASL_num2str(nVolumesPerRepetition) ' volumes per repetition']);
    % end



    
    
% =====================================================================
% C) Single- and multi-PLD subtraction
elseif bCreatePWI4D
    %% 4. Create PWI4D
    if nVolumes>1 && ~x.modules.asl.bContainsDeltaM
        % Paired subtraction
        [Control4D, LabelIm] = xASL_quant_GetControlLabelOrder(ASL4D);
        PWI4D = Control4D - LabelIm;
        
        % Skip every other value in the vectors as they were stored for both control and label images 
        x.Q.EchoTime_PWI4D = x.Q.EchoTime(1:2:end);
        x.Q.InitialPLD_PWI4D = x.Q.Initial_PLD(1:2:end);
        x.Q.LabelingDuration_PWI4D = x.Q.LabelingDuration(1:2:end);
    else % the same but then without subtraction
        PWI4D = ASL4D;
        Control4D = NaN;
    end

end
   
if bCreatePWI3D
    %% 5. Create PWI3D
    
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