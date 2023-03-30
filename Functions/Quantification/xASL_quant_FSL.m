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
	
	% FSLOptions is a character array containing CLI args for the BASIL/FABBER command
	FSLOptions = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, PWI);
        
    %% 5. Run Basil and retrieve CBF output
    if bUseFabber
        [~, resultFSL] = xASL_fsl_RunFSL(['fabber_asl' cmdlineOptions ' ' FSLOptions], x);
    else
        [~, resultFSL] = xASL_fsl_RunFSL(['basil' cmdlineOptions ' ' FSLOptions], x);
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
		ABV_map = ABV_map ./ x.Q.LabelingEfficiency;
	end
    
    %% 7. Householding
	% Basils Output is in the subfolder '/FSL_Output' which contains multiple steps if there are multiple iterations, and always contains
    % a symbolic link (symlink) to the foldername of the latest iteration/step ('stepX_latest').
    xASL_delete(pathFSLInput);
	xASL_delete(pathFSLOptions);
	xASL_delete(pathFSLOutput);
	%xASL_adm_DeleteFileList(x.dir.SESSIONDIR, '(?i)^.*basil.*$', 1, [0 Inf]);
    
end

function [FSLOptions] = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, PWI)
%xASL_sub_FSLOptions generates the options and saves them in a file and returns some commandline options as well
%
% FORMAT: [FSLOptions] = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, PWI)
% 
% INPUT:
%   pathFSLOptions  - filepath to the options file (REQUIRED)
%   x               - struct containing pipeline environment parameters (REQUIRED)
%   bUseFabber      - Use FABBER, alternative BASIL (REQUIRED)
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
% 5. Extra BASIL fitting options
% 6. Save and close the options file
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [FSLOptions] = xASL_sub_FSLOptions(pathFSLOptions, x, bUseFabber, PWI)
%
% __________________________________
% Copyright 2015-2023 ExploreASL 

%% 0. Admin
if ~bUseFabber
	% Set BASIL dataPar options and their defaults
	if ~isfield(x.Q,'BASIL')
		x.Q.BASIL = [];
	end

	if ~isfield(x.Q.BASIL,'bSpatial') || isempty(x.Q.BASIL.bSpatial)
		fprintf('BASIL: Setting default option bSpatial = true\n');
		x.Q.BASIL.bSpatial = true;
	end

	if ~isfield(x.Q.BASIL,'bInferT1') || isempty(x.Q.BASIL.bInferT1)
		fprintf('BASIL: Setting default option bInferT1 = false\n');
		x.Q.BASIL.bInferT1 = false;
	end

	if ~isfield(x.Q.BASIL,'bInferArt') || isempty(x.Q.BASIL.bInferArt)
		fprintf('BASIL: Setting default option bInferArt = true\n');
		x.Q.BASIL.bInferArt = true;
	end

	if ~isfield(x.Q.BASIL,'Exch') || isempty(x.Q.BASIL.Exch)
		fprintf('BASIL: Setting default option Exch = simple\n');
		x.Q.BASIL.Exch = 'simple';
	end

	if ~isfield(x.Q.BASIL,'Disp') || isempty(x.Q.BASIL.Disp)
		fprintf('BASIL: Setting default option Disp = none\n');
		x.Q.BASIL.Disp = 'none';
	end

	if ~isfield(x.Q.BASIL, 'ATTSD')
		x.Q.BASIL.ATTSD = 1.0;
	end

end

%% 1. Create the options file
% FSLOptions is a character array containing CLI args for the Basil command
FSLOptions = '';

FIDoptionFile = fopen(pathFSLOptions, 'w+');
if bUseFabber
	fprintf(FIDoptionFile, '# FABBER options written by ExploreASL\n');
else
	fprintf(FIDoptionFile, '# BASIL options written by ExploreASL\n');
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
		if bUseFabber
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
		else

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

end

if ~bUserFabber
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

	if x.modules.asl.bMultiPLD
		% Multi-PLD or Time Encoded data allows to fit arrival times
		fprintf(FIDoptionFile, '--batsd=%f\n', x.Q.BASIL.ATTSD);
	end
end


%% 5. Extra BASIL fitting options
if ~bUseFabber
	if x.Q.BASIL.bSpatial
		fprintf('BASIL: Use automated spatial smoothing\n');
		FSLOptions = [FSLOptions ' --spatial'];
	end

	if x.Q.BASIL.bInferT1
		fprintf('BASIL: Infer variable T1 values\n');
		FSLOptions = [FSLOptions ' --infert1'];
	end

	if x.Q.BASIL.bInferArt
		fprintf('BASIL: Infer arterial BV and arrival time\n');
		FSLOptions = [FSLOptions ' --inferart'];
	end

	switch (x.Q.BASIL.Exch)
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
			warning(['BASIL Exchange model: ' x.Q.BASIL.Exch ' not recognized.'])
	end

	switch (x.Q.BASIL.Disp)
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
			warning(['BASIL Dispersion model: ' x.Q.BASIL.Disp ' not recognized.'])
	end
	
	% 	%% Aquisition options we might be able to use in the future
	%   fprintf(option_file, '--sliceband=%i\n', sliceband);
	%   fprintf('BASIL: Multi-band setup with number of slices per band: %i\n', slicedband);
	%
	% 	fprintf(option_file, '--t1im=%s\n', t1im)
	%   fprintf('BASIL: Using supplied T1 (tissue) image in BASIL: %s\n', $t1im)
	%
	%   % For small numbers of time points we need an informative noise prior.
	%   % The user can specify an assumed SNR for this, or give prior estimated noise standard deviation below.
	%   if ~isfield(x.Q.BASIL,'SNR') || isempty(x.Q.BASIL.SNR)
	%       x.Q.BASIL.SNR = 0; %Otherwise, a reasonable default is 10
	%   end
	%
	%   if x.Q.BASIL.SNR
	%       fprintf('BASIL: Using SNR of %f to set noise std dev\n', x.Q.BasilSNR);
	%       % Estimate signal magntiude FIXME brain mask assume half of voxels
	%       mag_max = max(PWI, [], 4);
	%       brain_mag = 2*xASL_stat_MeanNan(mag_max(:));
	%       fprintf('BASIL: Mean maximum signal across brain: %f\n', brain_mag);
	%       % This will correspond to whole brain CBF (roughly) - about 0.5 of GM
	%       imNoiseSD = sqrt(brain_mag * 2 / x.Q.BasilSNR);
	%       fprintf('BASIL: Using a prior noise std.dev. of: %f\n', imNoiseSD);
	%       fprintf(FIDoptionFile, '--prior-noise-stddev=%f\n', imNoiseSD);
	%     end
end

%% 6. Close options file
fclose(FIDoptionFile);

end
