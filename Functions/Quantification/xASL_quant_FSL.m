function [CBF_nocalib, ATT_map, ABV_map, Tex_map, resultFSL] = xASL_quant_FSL(path_PWI4D, x)
%xASL_quant_FSL Perform quantification using FSL BASIL/FABBER
%
% FORMAT: [CBF_nocalib, ATT_map, ABV_map, Tex_map, resultFSL] = xASL_quant_FSL(path_PWI4D, x)
% 
% INPUT:
%   path_PWI4D      - path to PWI4D (OPTIONAL, defaults to x.P.Path_PWI4D)
%                     e.g., alternatives could be another space, or a concatenated PWI4D
%   x               - struct containing pipeline environment parameters (REQUIRED)
%
% OUTPUT:
% CBF_nocalib       - Quantified CBF image
%                     (if there is no FSL/BASIL installed, we return the original PWI)
% ATT_map           - ATT map (if possible to calculate with multi-PLD, otherwise empty)
% ABV_map           - arterial blood volume map (if possible to calculate with multi-PLD, otherwise empty)
% Tex_map           - Time of exchange map of transport across BBB (if possible to calculate with multi-TE, otherwise empty)
% resultFSL         - describes if the execution was successful
%                     (0 = successful, NaN = no FSL/BASIL found, 1 or other = something failed)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script performs quantification of the PWI using the FSL Basil/Fabber pipeline. Final calibration to
%              physiological units is performed by dividing the quantified PWI by the M0 image/value.
%              Fabber is used instead of Basil for multiTE data.
%
%              This function performs the following steps:
%
% 1. Define paths
% 2. Delete previous BASIL/Fabber output
% 3. Write the PWI as Nifti file for BASIL/Fabber to read as input
% 4. Create option_file that contains options which are passed to the FSL command
% 5. Run BASIL and retrieve CBF output
% 6. Scaling to physiological units
% 7. Householding
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: CBF_nocalib = xASL_quant_FSL(PWI, x);
%
% __________________________________
% Copyright 2015-2024 ExploreASL 

    %% 0. Admin
    fprintf('%s\n','Quantification CBF using FSL BASIL/FABBER:');   

    if nargin<1 || isempty(path_PWI4D)
        path_PWI4D = x.P.Path_PWI4D;
    end

    % Define defaults
	Tex_map = [];
	ATT_map = [];
	ABV_map = [];

	if ~isfield(x.modules.asl, 'bCleanUpBASIL') || isempty(x.modules.asl.bCleanUpBASIL)
		x.modules.asl.bCleanUpBASIL = true;
	end

	% Define if BASIL or FABBER is used - multiTE needs FABBER 
	if isfield(x.modules.asl, 'bQuantifyMultiTE') && x.modules.asl.bQuantifyMultiTE
		bUseFabber = 1;
        FSLfunctionName = 'fabber_asl';
	else
		bUseFabber = 0;
        FSLfunctionName = 'basil';
	end
    
    %% 1. Define temporary paths for FSL
	% Create FSL output directory
    dirFSLOutput = 'FSL_Output';
	pathFSLOutput = fullfile(x.dir.SESSIONDIR, dirFSLOutput);
	xASL_adm_CreateDir(pathFSLOutput)

    % For input, output, and options
    pathFSLInput = fullfile(pathFSLOutput, 'PWI4D_FSLInput.nii');
    pathFSLOptions = fullfile(pathFSLOutput, 'FSL_ModelOptions.txt');

    %% 2. Delete previous output
    xASL_adm_DeleteFileList(x.dir.SESSIONDIR, ['(?i)^' dirFSLOutput '.*$'], 1, [0 Inf]);
    FolderList = xASL_adm_GetFileList(x.dir.SESSIONDIR, ['(?i)^' dirFSLOutput '.*$'], 'FPList', [0 Inf], 1);
    for iFolder=1:numel(FolderList)
        xASL_delete(FolderList{iFolder}, 1);
    end
    fprintf('%s\n', 'Note that any file not found warnings can be ignored, this pertains to the use of symbolic links by BASIL/FABBER');
    
    % Remove residual BASIL-related files
    %xASL_delete(pathFSLOptions); % These files are now deleted as part of the FSLOutput directory
    %xASL_delete(pathFSLInput);
	xASL_delete(pathFSLOutput, 1);
    
    %% 3. Write the PWI4D as Nifti file for BASIL/FABBER to read as input
    [PWI4D, PWI4D_json] = xASL_io_Nifti2Im(path_PWI4D, [], [], true);
	PWI4D_size = size(PWI4D);
    
    % First, we extrapolate values to fill NaNs with a small kernel only, inside the brainmask
    % We don't have a brainmask here yet, so now we just run this small kernel once
    voxelSize = PWI4D_nii.hdr.pixdim(2:4);
    kernelSize = round([8 8 8]./voxelSize);
	for i4D = 1:size(PWI4D,4)
        PWI4D(:,:,:,i4D) = xASL_im_ndnanfilter(PWI4D(:,:,:,i4D), 'gauss', double(kernelSize), 2);
	end

	% Then, we extrapolate all outside the brain mask to ensure that there are no NaNs left
	PWI4D = xASL_im_FillNaNs(PWI4D, 1, 1, voxelSize);

    xASL_io_SaveNifti(path_PWI4D, pathFSLInput, PWI4D);

    %% 4. Create option_file that contains options which are passed to the FSL command
    % FSLOptions is a character array containing CLI args for the BASIL/FABBER command
	FSLOptions = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, PWI4D_size, PWI4D_json, pathFSLInput, pathFSLOutput);

    %% 5. Run BASIL and retrieve CBF output
    [~, resultFSL] = xASL_fsl_RunFSL([FSLfunctionName ' ' FSLOptions], x);
    
    % Check if FSL failed
    if isnan(resultFSL)
        error([FSLfunctionName ' was not found, exiting...']);
    elseif resultFSL~=0
		error(['Something went wrong running ' FSLfunctionName '...']);
    end
    
    fprintf('%s\n', 'The following warning (if mentioned above) can be ignored:');
    fprintf('%s\n', '/.../fsl/bin/basil: line 124: imcp: command not found');

    % CBF/nocalib, mean fit (->> is this what "ftiss" means?)
    pathBasilCBF = xASL_adm_GetFileList(pathFSLOutput, '^mean_ftiss\.nii$', 'FPListRec');
    
	if isempty(pathBasilCBF)
        error([FSLfunctionName ' failed']);
	end
   
    pathBasilCBF = pathBasilCBF{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
    CBF_nocalib = xASL_io_Nifti2Im(pathBasilCBF);
    
    % ATT
	pathBasilATT = xASL_adm_GetFileList(pathFSLOutput, '^mean_delttiss\.nii$', 'FPListRec');
	if ~isempty(pathBasilATT)
		ATT_map = xASL_io_Nifti2Im(pathBasilATT{end}); % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
	end
    
    % ABV
    pathBasilABV = xASL_adm_GetFileList(pathFSLOutput, '^mean_fblood\.nii$', 'FPListRec');
	if ~isempty(pathBasilABV)
		ABV_map = xASL_io_Nifti2Im(pathBasilABV{end}); % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
	end
    
    % Tex
    pathFabberTex = xASL_adm_GetFileList(pathFSLOutput, '^mean_T_exch\.nii$', 'FPListRec');
	if ~isempty(pathFabberTex)
		Tex_map = xASL_io_Nifti2Im(pathFabberTex{end}); % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
	end
	
    %% 6. Scaling to physiological units
    % Note different to xASL_quant_ASL since Fabber has T1 in seconds
    % and does not take into account labeling efficiency
    
    CBF_nocalib = CBF_nocalib .* 6000 .* x.Q.Lambda ./ x.Q.LabelingEfficiency;
    % (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
    % & old Siemens sequence also didn't need the 1 gr->100 gr conversion

	ABV_map = ABV_map ./ x.Q.LabelingEfficiency;
    
    %% 7. Householding
	% Basils Output is in the subfolder '/FSL_Output' which contains multiple steps if there are multiple iterations, and always contains
    % a symbolic link (symlink) to the foldername of the latest iteration/step ('stepX_latest').
	
	if x.modules.asl.bCleanUpBASIL
		%xASL_delete(pathFSLInput); % These files are now part of FSLOutput DIR, so are deleted with the rest of the dir
		%xASL_delete(pathFSLOptions);
		xASL_delete(pathFSLOutput, 1);
	end
    
end

function [FSLOptions] = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, sizePWI4D, jsonPWI4D, pathFSLInput, pathFSLOutput)
%xASL_sub_FSLOptions generates the options and saves them in a file and returns some commandline options as well
%
% FORMAT: [FSLOptions] = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, sizePWI4D, jsonPWI4D, pathFSLInput, pathFSLOutput)
% 
% INPUT:
%   pathFSLOptions  - filepath to the options file (REQUIRED)
%   x               - struct containing pipeline environment parameters (REQUIRED)
%   bUseFabber      - Use FABBER, alternative BASIL (REQUIRED)
%   sizePWI4D       - size of the PWI4D matrix (REQUIRED)
%   jsonPWI4D       - JSON in Legacy of the PWI4D containing LD, PLD, JSON (REQUIRED)
%   pathFSLInput    - Path to the data input file (REQUIRED)
%   pathFSLOutput   - Path to the output directory (REQUIRED)
%
% OUTPUT:
% FSLOptions      - command-line options
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Options-file is saved and commandline options returned in a single string
%
% 0. Admin
% 1. Create the options file
% 2. Basic model and tissue parameters
% 3. Basic acquisition parameters
% 4. Model fiting parameters
% 5. Extra BASIL fitting options
% 6. Save and close the options file
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [FSLOptions] = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, sizePWI4D, jsonPWI4D, pathFSLInput, pathFSLOutput)
%
% __________________________________
% Copyright 2015-2024 ExploreASL 

%% 0. Admin
if nargin<7 
	error('Require 7 input parameters.');
end

% Set BASIL dataPar options and their defaults
if ~isfield(x.modules.asl,'bMaskingBASIL') || isempty(x.modules.asl.bMaskingBASIL)
		fprintf('BASIL: Setting default option bMasking = true\n');
		x.modules.asl.bMaskingBASIL = true;
end

% Set basic parameters newly as they might differ in case of a merged sequence
%bQuantifyMultiPLD = x.modules.asl.bQuantifyMultiPLD;
if length(unique(jsonPWI4D.Q.Initial_PLD))>1 && length(unique(jsonPWI4D.Q.LabelingDuration))>1
	bQuantifyMultiPLD = true;
else
	bQuantifyMultiPLD = false;
end

if ~bUseFabber
	% On Low quality settings, turn off all extra processing options
	if isfield(x, 'settings') && isfield(x.settings, 'Quality') && ~x.settings.Quality
		x.modules.asl.bSpatialBASIL = false;
		x.modules.asl.bInferT1BASIL = false;
		x.modules.asl.bInferArtBASIL = false;
		x.modules.asl.ExchBASIL = 'simple';
		x.modules.asl.DispBASIL = 'none';
	end

	% Setting defaults for BASIL specific options
	if ~isfield(x.modules.asl,'bSpatialBASIL') || isempty(x.modules.asl.bSpatialBASIL)
		fprintf('BASIL: Setting default option bSpatial = false\n');
		x.modules.asl.bSpatialBASIL = false;
	end

	if ~isfield(x.modules.asl,'bInferT1BASIL') || isempty(x.modules.asl.bInferT1BASIL)
		fprintf('BASIL: Setting default option bInferT1 = false\n');
		x.modules.asl.bInferT1BASIL = false;
	end

	if ~isfield(x.modules.asl,'bInferArtBASIL') || isempty(x.modules.asl.bInferArtBASIL)
		fprintf('BASIL: Setting default option bInferArt = true\n');
		x.modules.asl.bInferArtBASIL = true;
	end

	if ~isfield(x.modules.asl,'ExchBASIL') || isempty(x.modules.asl.ExchBASIL)
		fprintf('BASIL: Setting default option Exch = simple\n');
		x.modules.asl.ExchBASIL = 'simple';
	end

	if ~isfield(x.modules.asl,'DispBASIL') || isempty(x.modules.asl.DispBASIL)
		fprintf('BASIL: Setting default option Disp = none\n');
		x.modules.asl.DispBASIL = 'none';
	end

	if ~isfield(x.modules.asl, 'ATTSDBASIL') || isempty(x.modules.asl.ATTSDBASIL)
		x.modules.asl.ATTSDBASIL = 1.0;
	end
end

%% 1. Create the options file
% FSLOptions is a character array containing CLI args for the Basil command
% Path to the options file
if bUseFabber
	FSLOptions = ['-@ ' xASL_adm_UnixPath(pathFSLOptions, ispc)];
else
	FSLOptions = ['--optfile ' xASL_adm_UnixPath(pathFSLOptions, ispc)];
end

FIDoptionFile = fopen(pathFSLOptions, 'w+');
if bUseFabber
	fprintf(FIDoptionFile, '# FABBER options written by ExploreASL\n');
else
	fprintf(FIDoptionFile, '# BASIL options written by ExploreASL\n');
end

% Define basic paths
if bUseFabber
	fprintf(FIDoptionFile, '--output=%s\n', xASL_adm_UnixPath(pathFSLOutput, ispc));
    fprintf(FIDoptionFile, '--data=%s\n', xASL_adm_UnixPath(pathFSLInput, ispc));

else
	% Path to input and output
	FSLOptions = [FSLOptions ' -o ' xASL_adm_UnixPath(pathFSLOutput, ispc)];
	FSLOptions = [FSLOptions ' -i ' xASL_adm_UnixPath(pathFSLInput, ispc)];
end

% Define masking
if x.modules.asl.bMaskingBASIL
	% Check for uninitialized Mask variable or file
	if ~isfield(x.P, 'Path_BrainMaskProcessing')
		warning('BASIL masking set to TRUE, but the mask variable x.P.Path_BrainMaskProcessing is not initialized.');
	elseif ~xASL_exist(x.P.Path_BrainMaskProcessing, 'file')
		warning('BASIL masking set to TRUE, but the mask is missing: %s\n', x.P.Path_BrainMaskProcessing);
	else
		% Add the mask to the options file
		if bUseFabber
			fprintf(FIDoptionFile, '--mask=%s\n', xASL_adm_UnixPath(x.P.Path_BrainMaskProcessing, ispc));
		else
			FSLOptions = [FSLOptions ' -m ' xASL_adm_UnixPath(x.P.Path_BrainMaskProcessing, ispc)];
		end
	end
end

%% 2. Basic model and tissue parameters
% Basic model options
if bUseFabber
    fprintf(FIDoptionFile, '--method=vb\n');
	fprintf(FIDoptionFile, '--model=asl_multite\n');
	fprintf(FIDoptionFile, '--infertexch\n'); % Fit Tex
	fprintf(FIDoptionFile, '--inferitt\n');   % Fit ATT
end

% Basic fitting and output options
if bUseFabber
	fprintf(FIDoptionFile, '--save-var\n');
	fprintf(FIDoptionFile, '--save-residuals\n');
	fprintf(FIDoptionFile, '--allow-bad-voxels\n');
	fprintf(FIDoptionFile, '--save-model-fit\n');
	fprintf(FIDoptionFile, '--noise=white\n');
end

% Basic tissue parameters
fprintf(FIDoptionFile, '--t1b=%f\n', x.Q.BloodT1/1000);
fprintf(FIDoptionFile, '--t1=%f\n', x.Q.TissueT1/1000);

if bUseFabber
	% T2-times needed for multi-TE quantification
	fprintf(FIDoptionFile, '--t2b=%f\n', x.Q.T2art/1000);
	fprintf(FIDoptionFile, '--t2=%f\n', x.Q.T2/1000);
end

%% 3. Basic acquisition parameters
switch lower(x.Q.LabelingType)
	% PASL quantification
	case 'pasl'
		% PASL model is assumed by default and does not need to be specified in the config file
		fprintf('BASIL: PASL model\n');

		% For PASL, there can be only a single LabelingDuration, so unique PLD+LabDur combinations are uniquely based on PLDs
		TIs = jsonPWI4D.Q.Initial_PLD'/1000;

		% Print all the TIs
		if bQuantifyMultiPLD	
			for iTI = 1:length(TIs)
				fprintf(FIDoptionFile, '--ti%d=%.2f\n', iTI, TIs(iTI));
			end
		else
			fprintf(FIDoptionFile, '--ti=%.2f\n', TIs);
		end

		% Either print bolus duration or unspecify it
		if isfield(jsonPWI4D.Q, 'LabelingDuration') && jsonPWI4D.Q.LabelingDuration
			if length(unique(jsonPWI4D.Q.LabelingDuration))>1
				warning('PASL multi-PLD currently supports only a single Labeling Duration');
			end
			fprintf(FIDoptionFile, '--tau=%.2f\n', x.Q.LabelingDuration(1)/1000);
		else
			% Bolus duration not know. If multi-TI, then try to infer it
			if length(TIs) > 1
				fprintf(FIDoptionFile, '--infertau\n');
				fprintf('BASIL: Infer bolus duration component\n')
			end
		end
vvvvvvvvvvvvvvvvvvv
	% CASL and PCASL quantification
	case {'casl','pcasl'}
		% Prepare unique PLDs+LabDur combinations
		
		% First create a labeling duration vector of the same length
		if length(x.Q.LabelingDuration)>1
			LabDurs = x.Q.LabelingDuration/1000;
		else
			LabDurs = ones(size(x.Q.Initial_PLD))*x.Q.LabelingDuration/1000;
		end

		PLDs = x.Q.Initial_PLD/1000;

        %% #1543 TEMPORARY SOLUTION, MERGING GOES OUT & to xASL_im_MergePWI4D

		% For Time-encoded, we skip the first volume per block
		if x.modules.asl.bTimeEncoded
			% In case we are merging sessions, we have to read all JSONs again from files
			if bMergingSessions
				PLDs = [];
				LabDurs = [];
				nTE = [];
				TEs = [];
				for iSession = 1: numel(x.modules.asl.sessionsToMerge)
					sessionsToMergeJSONdir = fullfile(x.dir.SUBJECTDIR, x.modules.asl.sessionsToMerge{iSession}, 'ASL4D.json');

					sessionsToMergeJSON = xASL_io_ReadJson(sessionsToMergeJSONdir);
					NewPLDs = sessionsToMergeJSON.PostLabelingDelay;
					NewLabDurs = sessionsToMergeJSON.LabelingDuration;

					NewEchoTimes = unique(sessionsToMergeJSON.EchoTime);
					if length(NewEchoTimes) < length(NewPLDs)
						NewEchoTimes = repmat(NewEchoTimes, [length(NewPLDs)/length(NewEchoTimes), 1]);
					end

					% For Time-encoded, we skip the first volume per block -> the same as above
					if x.modules.asl.bTimeEncoded
						NewTimeEncodedMatrixSize = sessionsToMergeJSON.TimeEncodedMatrixSize;
						[NewPLDs, ~] = unique(NewPLDs, 'stable');

						numberBlocks = numel(NewPLDs)/NewTimeEncodedMatrixSize;
						index = (ones(numberBlocks,1)*(2:NewTimeEncodedMatrixSize) + (0:(numberBlocks-1))' * NewTimeEncodedMatrixSize * ones(1,NewTimeEncodedMatrixSize-1))';
						NewPLDs = NewPLDs(index(:));
						NewEchoTimes = NewEchoTimes((numberBlocks*length(unique(NewEchoTimes))+1):end);
					else
						% For normal multi-timepoint, we look for unique PLD+LabDur combinations
						[~, indexNew, ~] = unique([NewPLDs(:), NewLabDurs(:)], 'stable', 'rows');

						NewPLDs = NewPLDs(indexNew);
					end
					if length(NewLabDurs)==1
						NewLabDurs = ones(size(NewPLDs))*NewLabDurs;
					end
					NewNTE = ones(size(NewPLDs))*length(unique(NewEchoTimes));

					PLDs = [PLDs; NewPLDs];
					LabDurs = [LabDurs; NewLabDurs];
					nTE = [nTE; NewNTE];
					TEs = [TEs; NewEchoTimes];
				end
			else
				[PLDs, index] = unique(PLDs, 'stable');
				LabDurs = LabDurs(index)';

				numberBlocks = numel(PLDs)/x.Q.TimeEncodedMatrixSize;
				index = (ones(numberBlocks,1)*(2:x.Q.TimeEncodedMatrixSize) + (0:(numberBlocks-1))' * x.Q.TimeEncodedMatrixSize * ones(1,x.Q.TimeEncodedMatrixSize-1))';
				PLDs = PLDs(index(:));
				LabDurs = LabDurs(index(:));
			end
		else
			% For normal multi-timepoint, we look for unique PLD+LabDur combinations
			[~, indexNew, ~] = unique([PLDs(:), LabDurs(:)], 'stable', 'rows');

			PLDs = PLDs(indexNew);
			LabDurs = LabDurs(indexNew);
		end

		if bUseFabber
			%Echo Times 
			nPLD = length(PLDs); % Number of PLDs in the PLD vector
			nVolume = size(PWI,4); % Number of volumes in PWI

			nTE = length(unique(x.Q.EchoTime)); % Calculate the number of Echo Times
			TEs = round(x.Q.EchoTime'/1000,3); % Convert Echo Times to seconds and keep 4 decimal digits
			
			if (nPLD*nTE) ~= nVolume
				error('The number of volumes %d does not match number of PLDs %d * number of TEs %d', nVolume, nPLD, nTE);
			end
			nTE = ones(size(PLDs))*nTE;

			% Printing the values in the FSL option file (PLD=ti, LD=tau)
			for iPLD = 1:length(PLDs)
				fprintf(FIDoptionFile, '--ti%d=%.2f\n', iPLD, PLDs(iPLD) + LabDurs(iPLD));
			end
			for iNTE = 1:length(nTE)
				fprintf(FIDoptionFile, '--nte%d=%d\n', iNTE, nTE(iNTE)); % --nte1=8 --nte2=8 --nte3=8 (if nTE=8)
			end

			if length(nTE) == 1 && nTE == 1
				% For a single-TE, we have to repeat it for each volume
				for iTE = 1:nVolume %We need a TE for each volume
					fprintf(FIDoptionFile, '--te%d=%.3f\n', iTE, TEs(1));
				end
			else
				% For multi-TE, we print all of them
				for iTE = 1:nVolume %We need a TE for each volume
					fprintf(FIDoptionFile, '--te%d=%.3f\n', iTE, TEs(iTE));
				end
			end

			% Right now, we assume that we have averaged over PLDs
			%fprintf(FIDoptionFile, '--repeats=%i\n', size(PWI, 4)/PLDAmount);
			%fprintf(FIDoptionFile, '--repeats=1\n');
		else
			% Specify that we run the PCASL/CASL model
			fprintf(FIDoptionFile, '--casl\n');
			fprintf('BASIL: (P)CASL model\n');

			% For BASIL, PLDs are specified
			if bQuantifyMultiPLD
				for iPLD = 1:length(PLDs)
					fprintf(FIDoptionFile, '--pld%d=%.2f\n', iPLD, PLDs(iPLD));
				end
			else
				fprintf(FIDoptionFile, '--pld=%.2f\n', PLDs);
			end
		end

		% Print labeling durations
		if bQuantifyMultiPLD
			for iLabDurs = 1:length(LabDurs)
				fprintf(FIDoptionFile, '--tau%d=%.2f\n', iLabDurs, LabDurs(iLabDurs));
			end
		else
			fprintf(FIDoptionFile, '--tau=%.2f\n', LabDurs);
		end
end
^^^^^^^^^^^^^^^^^^^^^^^^^^
if ~bUseFabber
	% Right now, we assume that we have averaged over PLDs
	%fprintf(FIDoptionFile, '--repeats=%i\n', size(PWI, 4)/PLDAmount);
	fprintf(FIDoptionFile, '--repeats=1\n');

	% Slice-timing
	fprintf(FIDoptionFile, '--slicedt=%f\n', x.Q.BasilSliceReadoutTime/1000);

	if isfield(x.Q,'LookLocker') && x.Q.LookLocker
		if isfield(x.Q,'FlipAngle')
			if length(unique(x.Q.FlipAngle))>1
				warning('Look-Locker quantification with multiple flip angles, e.g. QUASAR, is not implemented yet');
			end
			fprintf(option_file, '--FA=%f\n', x.Q.FlipAngle(1));
			fprintf('BASIL: Flip angle for Look-Locker readout: %f\n', x.Q.FlipAngle(1));
		else
			warning('BASIL: Unknown flip angle for Look-Locker\n');
		end
	end
end

%% 4. Model fiting parameters
if ~bUseFabber
	switch lower(x.Q.LabelingType)
		case 'pasl'
			% Default initial ATT for PASL is 0.7
			fprintf(FIDoptionFile, '--bat=0.7\n');
		case {'pcasl','casl'}
			% Default initial ATT for PCASL is 1.3
			fprintf(FIDoptionFile, '--bat=1.3\n');
	end

	if bQuantifyMultiPLD
		% Multi-PLD or Time Encoded data allows to fit arrival times
		fprintf(FIDoptionFile, '--batsd=%f\n', x.modules.asl.ATTSDBASIL);
	end
end

%% 5. Extra BASIL fitting options
if ~bUseFabber
	if x.modules.asl.bSpatialBASIL
		fprintf('BASIL: Use automated spatial smoothing\n');
		FSLOptions = [FSLOptions ' --spatial'];
	end

	if x.modules.asl.bInferT1BASIL
		if bQuantifyMultiPLD
			fprintf('BASIL: Infer variable T1 values\n');
			FSLOptions = [FSLOptions ' --infert1'];
		end
	end

	if x.modules.asl.bInferArtBASIL
		if bQuantifyMultiPLD
			fprintf('BASIL: Infer arterial BV and arrival time\n');
			FSLOptions = [FSLOptions ' --inferart'];
		end
	end

	switch (x.modules.asl.ExchBASIL)
		case 'simple'
			fprintf('BASIL Exchange model: Simple single compartment with T1 of blood, per white paper\n');
			FSLOptions = [FSLOptions ' --exch=simple'];
		case 'mix'
			fprintf('BASIL Exchange model: Well-mixed\n');
			FSLOptions = [FSLOptions ' --exch=mix'];
		case '2cpt'
			fprintf('BASIL Exchange model: A two compartment exchange model following Parkes & Tofts\n');
			FSLOptions = [FSLOptions ' --exch=2cpt'];
		case 'spa'
			fprintf('BASIL Exchange model: A single pass approximation from St. Lawrence\n');
			FSLOptions = [FSLOptions ' --exch=spa'];
		otherwise
			warning(['BASIL Exchange model: ' x.modules.asl.ExchBASIL ' not recognized.'])
	end

	if bQuantifyMultiPLD
		switch (x.modules.asl.DispBASIL)
			case 'none'
				fprintf('BASIL Dispersion model: none\n');
				FSLOptions = [FSLOptions ' --disp=none'];
			case 'gamma'
				fprintf('BASIL Dispersion model: Gamma\n');
				FSLOptions = [FSLOptions ' --disp=gamma'];
			case 'gauss'
				fprintf('BASIL Dispersion model: Temporal Gaussian dispersion kernel\n');
				FSLOptions = [FSLOptions ' --disp=gauss'];
			case 'sgauss'
				fprintf('BASIL Dispersion model: Spatial Gaussian dispersion kernel\n');
				FSLOptions = [FSLOptions ' --disp=sgauss'];
			otherwise
				warning(['BASIL Dispersion model: ' x.modules.asl.DispBASIL ' not recognized.'])
		end
	else
		fprintf('BASIL Dispersion model: none\n');
		FSLOptions = [FSLOptions ' --disp=none'];
	end


	% 	%% Aquisition options we might be able to use in the future
	%   fprintf(option_file, '--sliceband=%i\n', sliceband);
	%   fprintf('BASIL: Multi-band setup with number of slices per band: %i\n', slicedband);
	%
	% 	fprintf(option_file, '--t1im=%s\n', t1im)
	%   fprintf('BASIL: Using supplied T1 (tissue) image in BASIL: %s\n', $t1im)
	%
end

%% 6. Close options file
fclose(FIDoptionFile);

end
