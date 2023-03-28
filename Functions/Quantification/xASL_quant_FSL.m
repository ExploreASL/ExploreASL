function [CBF_nocalib, ATT_map, ABV_map, Tex_map, resultFSL] = xASL_quant_FSL(PWI, x)
%xASL_quant_FSL Perform quantification using FSL BASIL/FABBER
%
% FORMAT: [CBF_nocalib, ATT_map, ABV_map, Tex_map, resultFSL] = xASL_quant_FSL(PWI, x)
% 
% INPUT:
%   PWI             - image matrix of perfusion-weighted image (REQUIRED)
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
% 3. Write the PWI as Nifti file for Basil/Fabber to read as input
% 4. Create option_file that contains options which are passed to the FSL command
% 5. Run Basil and retrieve CBF output
% 6. Scaling to physiological units
% 7. Householding
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: CBF_nocalib = xASL_quant_FSL(PWI, x);
%
% __________________________________
% Copyright 2015-2023 ExploreASL 
    

    %% Admin
    fprintf('%s\n','Quantification CBF using FSL BASIL/FABBER:');   

	Tex_map = [];
	ATT_map = [];
	ABV_map = [];
    
    %% 1. Define paths
    % For input, output, and options
    pathFSLInput = fullfile(x.dir.SESSIONDIR, 'PWI4D_FSLInput.nii');
    pathFSLOptions = fullfile(x.dir.SESSIONDIR, 'FSL_ModelOptions.txt');
    dirFSLOutput = 'FSL_Output';
	pathFSLOutput = fullfile(x.dir.SESSIONDIR, dirFSLOutput);
	
	% Define if BASIL or FABBER is used - multiTE needs FABBER 
	if isfield(x.modules.asl, 'bMultiTE') && x.modules.asl.bMultiTE == 1
		bUseFabber = 1;
	else
		bUseFabber = 0;
	end
        
    %% 2. Delete previous output
    xASL_adm_DeleteFileList(x.dir.SESSIONDIR, ['(?i)^' dirFSLOutput '.*$'], 1, [0 Inf]);
    FolderList = xASL_adm_GetFileList(x.dir.SESSIONDIR, ['(?i)^' dirFSLOutput '.*$'], 'FPList', [0 Inf], 1);
    for iFolder=1:numel(FolderList)
        xASL_delete(FolderList{iFolder}, 1);
    end
    fprintf('%s\n', 'Note that any file not found warnings can be ignored, this pertains to the use of symbolic links by BASIL/FABBER');
    
    % Remove residual BASIL-related files
    xASL_delete(pathFSLOptions);
    xASL_delete(pathFSLInput);
	xASL_delete(pathFSLOutput);
    %xASL_adm_DeleteFileList(x.dir.SESSIONDIR, '(?i)^.*basil.*$', 1, [0 Inf]);
    
    %% 3. Write the PWI as Nifti file for Basil/Fabber to read as input
    % FIXME would be good to have a brain mask at this point -> PM: if this would be a brainmask as well, we can skip creating a dummy input image here
    PWI(isnan(PWI)) = 0;
    
    if ~x.modules.asl.bMultiPLD
        % SinglePLD
        xASL_io_SaveNifti(x.P.Path_PWI, pathFSLInput, PWI, [], 0); % use PWI path
    else
        % MultiPLD
        xASL_io_SaveNifti(x.P.Path_PWI4D, pathFSLInput, PWI, [], 0); % use PWI4D path
    end

    %% 4. Create option_file that contains options which are passed to the FSL command
    % Define basic paths

	% Path to the options file
	cmdlineOptions = [' --optfile ' xASL_adm_UnixPath(pathFSLOptions, ispc)];

	% Path to input and output
	if bUseFabber
		cmdlineOptions = [cmdlineOptions ' --output=' xASL_adm_UnixPath(pathFSLOutput, ispc)];
		cmdlineOptions = [cmdlineOptions ' --data=', xASL_adm_UnixPath(pathFSLInput, ispc)];
	else
		cmdlineOptions = [cmdlineOptions ' -o ' xASL_adm_UnixPath(pathFSLOutput, ispc)];
		cmdlineOptions = [cmdlineOptions ' -i ' xASL_adm_UnixPath(pathFSLInput, ispc)];
	end
	
	% basil_options is a character array containing CLI args for the Basil/Fabber command
	BasilOptions = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, PWI, pathFSLInput, dirFSLOutput);
        
    %% 5. Run Basil and retrieve CBF output
    if bUseFabber
        [~, resultFSL] = xASL_fsl_RunFSL(['fabber_asl' cmdlineOptions ' ' BasilOptions], x);
    else
        [~, resultFSL] = xASL_fsl_RunFSL(['basil' cmdlineOptions ' ' BasilOptions], x);
    end
    
    % Check if FSL failed
    if isnan(resultFSL) 
		if bUseFabber
			error('FSL FABBER was not found, exiting...');
		else
			error('FSL BASIL was not found, exiting...');
		end
    elseif resultFSL~=0
		if bUseFabber
			error('Something went wrong running FSL FABBER...');
		else
			error('Something went wrong running FSL BASIL...');
		end
    end
    
    fprintf('%s\n', 'The following warning (if mentioned above) can be ignored:');
    fprintf('%s\n', '/.../fsl/bin/basil: line 124: imcp: command not found');
    
    pathBasilCBF = xASL_adm_GetFileList(pathFSLOutput, '^mean_ftiss\.nii$', 'FPListRec');
    
	if isempty(pathBasilCBF)
		if bUseFabber
			error('FSL FABBER failed');
		else
			error('FSL BASIL failed');
		end
	end
   
    pathBasilCBF = pathBasilCBF{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
       
    CBF_nocalib = xASL_io_Nifti2Im(pathBasilCBF);
        
	pathBasilATT = xASL_adm_GetFileList(pathFSLOutput, '^mean_delttiss\.nii$', 'FPListRec');
	if ~isempty(pathBasilATT)
		pathBasilATT = pathBasilATT{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
		ATT_map = xASL_io_Nifti2Im(pathBasilATT);
	end
    
    pathBasilABV = xASL_adm_GetFileList(pathFSLOutput, '^mean_fblood\.nii$', 'FPListRec');
	if ~isempty(pathBasilABV)
		pathBasilABV = pathBasilABV{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
		ABV_map = xASL_io_Nifti2Im(pathBasilABV);
	end
    
    pathFabberTex = xASL_adm_GetFileList(pathFSLOutput, '^mean_T_exch\.nii$', 'FPListRec');
	if ~isempty(pathFabberTex)
		pathFabberTex = pathFabberTex{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
		Tex_map = xASL_io_Nifti2Im(pathFabberTex);
	end
	
    %% 6. Scaling to physiological units
    % Note different to xASL_quant_ASL since Fabber has T1 in seconds
    % and does not take into account labeling efficiency
    
    CBF_nocalib = CBF_nocalib .* 6000 .* x.Q.Lambda ./ x.Q.LabelingEfficiency;
    % (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
    % & old Siemens sequence also didn't need the 1 gr->100 gr conversion

	if numel(ABV_map) > 1
		ABV_map = ABV_map .* 6000 .* x.Q.Lambda ./ x.Q.LabelingEfficiency;
	end
    
    %% 7. Householding
	% Basils Output is in the subfolder '/FSL_Output' which contains multiple steps if there are multiple iterations, and always contains
    % a symbolic link (symlink) to the foldername of the latest iteration/step ('stepX_latest').
    xASL_delete(pathFSLInput);
	xASL_delete(pathFSLOptions);
	xASL_delete(pathFSLOutput);
	%xASL_adm_DeleteFileList(x.dir.SESSIONDIR, '(?i)^.*basil.*$', 1, [0 Inf]);
    
end




function [BasilOptions] = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, PWI)
if bUseFabber
	%% FABBER
	% Save a Fabber options file for FSL command
	% 1. Create an option file
	% 2. Basic tissue parameters
	% 3. Basic acquisition parameters
	% 4. Model fiting parameters
	% 5. Save Fabber options file

	%% 1. Create option_file that contains options which are passed to Fabber
	% basil_options is a character array containing CLI args for the Basil command
	BasilOptions = '';

	FIDoptionFile = fopen(pathFSLOptions, 'w+');
	fprintf(FIDoptionFile, '# Fabber options written by ExploreASL\n');
	fprintf(FIDoptionFile, '--method=vb\n');

	fprintf(FIDoptionFile, '--model=asl_multite\n');
	fprintf(FIDoptionFile, '--infertexch\n');
	fprintf(FIDoptionFile, '--save-var\n');
	fprintf(FIDoptionFile, '--save-residuals\n');
	fprintf(FIDoptionFile, '--allow-bad-voxels\n');
	fprintf(FIDoptionFile, '--save-model-fit\n');
	fprintf(FIDoptionFile, '--noise=white\n');

	%% 2. Basic tissue parameters
	fprintf(FIDoptionFile, '--t1b=%f\n', x.Q.BloodT1/1000);
	fprintf(FIDoptionFile, '--t1=%f\n', x.Q.TissueT1/1000);

	fprintf(FIDoptionFile, '--t2b=%f\n', x.Q.T2art/1000);
	fprintf(FIDoptionFile, '--t2=%f\n', x.Q.T2/1000);

	%% 3. Basic acquisition parameters

	% Prepare unique PLDs+LabDur combinations

	% First create a labeling duration vector of the same length
	if length(x.Q.LabelingDuration)>1
		LabDurs = x.Q.LabelingDuration/1000;
	else
		LabDurs = ones(size(x.Q.Initial_PLD))*x.Q.LabelingDuration/1000;
	end

	PLDs = x.Q.Initial_PLD/1000;

	% For Time-encoded, we skip the first volume per block
	if x.modules.asl.bTimeEncoded
		[PLDs, index] = unique(PLDs, 'stable');
		LabDurs = LabDurs(index)';

		numberBlocks = numel(PLDs)/x.Q.TimeEncodedMatrixSize;
		index = (ones(numberBlocks,1)*(2:x.Q.TimeEncodedMatrixSize) + (0:(numberBlocks-1))' * x.Q.TimeEncodedMatrixSize * ones(1,x.Q.TimeEncodedMatrixSize-1))';
		PLDs = PLDs(index(:)');
		LabDurs = LabDurs(index(:)');
	else
		% For normal multi-timepoint, we look for unique PLD+LabDur combinations
		[~, indexNew, ~] = unique([PLDs(:), LabDurs(:)], 'stable', 'rows');

		PLDs = PLDs(indexNew);
		LabDurs = LabDurs(indexNew);
	end

	%Echo Times
	nTE = length(unique(x.EchoTime));
	NVolumes = size(PWI,4);
	TEs = round(x.EchoTime(1:NVolumes)'/1000,4); % To keep 4 decimal digits

	% Plotting the values into the doc (PLD=ti, LD=tau)
	for iPLD = 1:length(PLDs)
		fprintf(FIDoptionFile, '--ti%d=%.2f\n', iPLD, PLDs(iPLD));
		fprintf(FIDoptionFile, '--nte%d=%d\n', iPLD, nTE); % --nte1=8 --nte2=8 --nte3=8 (if nTE=8)
	end
	for iLabDurs = 1:length(LabDurs)
		fprintf(FIDoptionFile, '--tau%d=%.2f\n', iLabDurs, LabDurs(iLabDurs));
	end

	for iTE = 1:length(TEs) %We need a TE for each volume
		fprintf(FIDoptionFile, '--te%d=%.2f\n', iTE, TEs(iTE));
	end

	% Right now, we assume that we have averaged over PLDs
	%fprintf(FIDoptionFile, '--repeats=%i\n', size(PWI, 4)/PLDAmount);
	%fprintf(FIDoptionFile, '--repeats=1\n');

	% 4. Model fiting parameters (ATT map)
	fprintf(FIDoptionFile, '--inferitt');

	%% 5. Save Fabber options file
	fclose(FIDoptionFile);

else
	%% BASIL
	% Save a Basil options file and store CLI options for Basil
	% 1. Create an option file
	% 2. Basic tissue parameters
	% 3. Basic acquisition parameters
	% 4. Model fiting parameters
	% 5. Extra features on demand
	% 6. Save BASIL options file

	%% 1. Create option_file that contains options which are passed to Fabber
	% basil_options is a character array containing CLI args for the Basil command

	FIDoptionFile = fopen(pathFSLOptions, 'w+');
	BasilOptions = '';

	fprintf(FIDoptionFile, '# Basil options written by ExploreASL\n');

	%% 2. Basic tissue parameters
	fprintf(FIDoptionFile, '--t1b=%f\n', x.Q.BloodT1/1000);
	fprintf(FIDoptionFile, '--t1=%f\n', x.Q.TissueT1/1000);

	%% 3. Basic acquisition parameters

	% Labelling type - PASL or pCASL
	switch lower(x.Q.LabelingType)
		case 'pasl'
			% PASL model is assumed by default and does not need to be specified in the config file
			fprintf('BASIL: PASL model\n');

			% For PASL, there can be only a single LabelingDuration, so unique PLD+LabDur combinations are uniquely based on PLDs
			TIs = (unique(x.Q.Initial_PLD, 'stable'))'/1000;

			% Print all the TIs
			if x.modules.asl.bMultiPLD
				% For Time-encoded, we skip the first volume
				if x.modules.asl.bTimeEncoded
					numberBlocks = numel(TIs)/x.Q.TimeEncodedMatrixSize;
					ind = (ones(numberBlocks,1)*(2:x.Q.TimeEncodedMatrixSize) + (0:(numberBlocks-1))' * x.Q.TimeEncodedMatrixSize * ones(1,x.Q.TimeEncodedMatrixSize-1))';
					TIs = TIs(ind(:)');
				end
				for iTI = 1:length(TIs)
					fprintf(FIDoptionFile, '--ti%d=%.2f\n', iTI, TIs(iTI));
				end

			else
				fprintf(FIDoptionFile, '--ti=%.2f\n', TIs);
			end

			% Either print bolus duration or unspecify it
			if isfield(x.Q,'LabelingDuration') && x.Q.LabelingDuration
				if length(unique(x.Q.LabelingDuration))>1
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

		case {'casl','pcasl'}
			% Specify that we run the PCASL/CASL model
			fprintf(FIDoptionFile, '--casl\n');
			fprintf('BASIL: (P)CASL model\n');

			% Prepare unique PLDs+LabDur combinations

			% First create a labeling duration vector of the same length
			if length(x.Q.LabelingDuration)>1
				LabDurs = x.Q.LabelingDuration/1000;
			else
				LabDurs = ones(size(x.Q.Initial_PLD))*x.Q.LabelingDuration/1000;
			end

			PLDs = x.Q.Initial_PLD/1000;

			% For Time-encoded, we skip the first volume per block
			if x.modules.asl.bTimeEncoded
				[PLDs, index] = unique(PLDs, 'stable');
				LabDurs = LabDurs(index)';

				numberBlocks = numel(PLDs)/x.Q.TimeEncodedMatrixSize;
				index = (ones(numberBlocks,1)*(2:x.Q.TimeEncodedMatrixSize) + (0:(numberBlocks-1))' * x.Q.TimeEncodedMatrixSize * ones(1,x.Q.TimeEncodedMatrixSize-1))';
				PLDs = PLDs(index(:)');
				LabDurs = LabDurs(index(:)');
			else
				% For normal multi-timepoint, we look for unique PLD+LabDur combinations
				[~, indexNew, ~] = unique([PLDs(:), LabDurs(:)], 'stable', 'rows');

				PLDs = PLDs(indexNew);
				LabDurs = LabDurs(indexNew);
			end

			if x.modules.asl.bMultiPLD
				for iPLD = 1:length(PLDs)
					fprintf(FIDoptionFile, '--pld%d=%.2f\n', iPLD, PLDs(iPLD));
				end
				for iLabDurs = 1:length(LabDurs)
					fprintf(FIDoptionFile, '--tau%d=%.2f\n', iLabDurs, LabDurs(iLabDurs));
				end
			else
				fprintf(FIDoptionFile, '--tau=%.2f\n', LabDurs);
				fprintf(FIDoptionFile, '--pld=%.2f\n', PLDs);
			end
	end

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

	%% 4. Model fiting parameters

	% This helps avoid failure on the structural-space image (we do not perform non-native space quantification yet)
	% fprintf(FIDoptionFile, '--allow-bad-voxels\n');

	switch lower(x.Q.LabelingType)
		case 'pasl'
			% Default initial ATT for PASL is 0.7
			fprintf(FIDoptionFile, '--bat=0.7\n');
		case {'pcasl','casl'}
			% Default initial ATT for PCASL is 1.3
			fprintf(FIDoptionFile, '--bat=1.3\n');
	end

	if x.modules.asl.bMultiPLD
		% Multi-PLD or Time Encoded data allows to fit arrival times

		% Set the variance of ATT estimation
		if ~isfield(x.Q, 'BasilATTSD')
			x.Q.BasilATTSD = 1.0;
		end
		fprintf(FIDoptionFile, '--batsd=%f\n', x.Q.BasilATTSD);
	end


	% 	%% FIXME Aquisition options we might be able to use in the future
	%     %fprintf(option_file, '--sliceband=%i\n', sliceband);
	%     %fprintf('BASIL: Multi-band setup with number of slices per band: %i\n', slicedband);
	%
	% 	%fprintf(option_file, '--t1im=%s\n', t1im)
	%     %fprintf('BASIL: Using supplied T1 (tissue) image in BASIL: %s\n', $t1im)

	%     %% Noise specification
	%     % For small numbers of time points we need an informative noise prior.
	%     % The user can specify an assumed SNR for this, or give prior estimated noise standard deviation below.
	%     if ~isfield(x.Q,'BasilSNR') || ~x.Q.BasilSNR
	%         x.Q.BasilSNR = 10;
	%     end
	%
	%     if size(PWI, 4) < 5
	%         x.Q.BasilNoisePrior = 1;
	%         fprintf('BASIL: Small number of volumes (%i < 5): informative noise prior will be used\n', size(PWI, 4));
	%     end
	%     if isfield(x.Q,'BasilNoisePrior') && x.Q.BasilNoisePrior
	%         % Use an informative noise prior
	%         if ~isfield(x.Q,'BasilNoiseSD') || ~x.Q.BasilNoiseSD
	%             fprintf('BASIL: Using SNR of %f to set noise std dev\n', x.Q.BasilSNR);
	%             % Estimate signal magntiude FIXME brain mask assume half of voxels
	%             mag_max = max(PWI, [], 4);
	%             brain_mag = 2*xASL_stat_MeanNan(mag_max(:));
	%             fprintf('BASIL: Mean maximum signal across brain: %f\n', brain_mag);
	%             % This will correspond to whole brain CBF (roughly) - about 0.5 of GM
	%             x.Q.BasilNoiseSD = sqrt(brain_mag * 2 / x.Q.BasilSNR);
	%         end
	%         fprintf('BASIL: Using a prior noise std.dev. of: %f\n', x.Q.BasilNoiseSD);
	%         fprintf(FIDoptionFile, '--prior-noise-stddev=%f\n', x.Q.BasilNoiseSD);
	%     end


	%% 5. Extra features on demand
	if ~isfield(x,'BasilSpatial') || isfield(x,'BasilSpatial') && isequal(x.Q.BasilSpatial,1)
		fprintf('BASIL: Instructing BASIL to use automated spatial smoothing by default\n');
		BasilOptions = [BasilOptions ' --spatial'];
	end

	if isfield(x.Q,'BasilInferT1') && x.Q.BasilInferT1
		fprintf('BASIL: Instructing BASIL to infer variable T1 values\n');
		BasilOptions = [BasilOptions ' --infert1'];
	end

	if isfield(x.Q,'BasilInferATT') && x.Q.BasilInferATT
		fprintf('BASIL: Infer arterial component');
		fprintf('BASIL: Variable arterial component arrival time');
		BasilOptions = [BasilOptions ' --inferart'];
	end

	if ~isfield(x.Q,'BasilExch')
		fprintf('BASIL: Using exchange model: simple single compartment with T1 of blood, per white paper as default%s\n', x.Q.BasilExch);
		BasilOptions = [BasilOptions ' --exch=simple'];
	elseif isfield(x.Q,'BasilExch') && strcmp(x.Q.BasilSpatial,'simple')
		fprintf('BASIL: Using exchange model: simple single compartment with T1 of blood, per white paper%s\n', x.Q.BasilExch);
		BasilOptions = [BasilOptions ' --exch=simple'];
	elseif isfield(x.Q,'BasilExch') && strcmp(x.Q.BasilSpatial,'mix')
		fprintf('BASIL: Using exchange model: well-mixed%s\n', x.Q.BasilExch);
		BasilOptions = [BasilOptions ' --exch=mix'];
	elseif isfield(x.Q,'BasilExch') && strcmp(x.Q.BasilSpatial,'mix')
		fprintf('BASIL: Using exchange model: a two compartment exchange model following Parkes & Tofts%s\n', x.Q.BasilExch);
		BasilOptions = [BasilOptions ' --exch=2cpt'];
	elseif isfield(x.Q,'BasilExch') && strcmp(x.Q.BasilSpatial,'mix')
		fprintf('BASIL: Using exchange model: isngle pass approximiation from St. Lawrence%s\n', x.Q.BasilExch);
		BasilOptions = [BasilOptions ' --exch=spa'];
	end

	if ~isfield(x.Q,'BasilDisp')
		fprintf('BASIL: Using no dispersion model as default %s\n', x.Q.BasilDisp);
		BasilOptions = [BasilOptions ' --disp=none'];
	elseif isfield(x.Q,'BasilDisp') && strcmp(x.Q.BaselDisp,'none')
		fprintf('BASIL: Using no dispersion model %s\n', x.Q.BasilDisp);
		BasilOptions = [BasilOptions ' --disp=none'];
	elseif isfield(x.Q,'BasilDisp') && strcmp(x.Q.BaselDisp,'gamma')
		fprintf('BASIL: Using dispersion model: Gamma %s\n', x.Q.BasilDisp);
		BasilOptions = [BasilOptions ' --disp=gamma'];
	elseif isfield(x.Q,'BasilDisp') && strcmp(x.Q.BaselDisp,'gauss')
		fprintf('BASIL: Using dispersion model: Temporal Gaussian dispersion kernel%s\n', x.Q.BasilDisp);
		BasilOptions = [BasilOptions ' --disp=gauss'];
	elseif isfield(x.Q,'BasilDisp') && strcmp(x.Q.BaselDisp,'gamma')
		fprintf('BASIL: Using dispersion model: Spatial Gaussian dispersion kernel%s\n', x.Q.BasilDisp);
		BasilOptions = [BasilOptions ' --disp=sgauss'];
	end

	if isfield(x.Q,'BasilDebug') && x.Q.BasilDebug
		BasilOptions = [BasilOptions ' --devel'];
	end


	%% 6. Save Basil options file
	fclose(FIDoptionFile);

end
end
