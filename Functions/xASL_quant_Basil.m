function [CBF_nocalib, resultFSL] = xASL_quant_Basil(PWI, x)
%xASL_quant_Basil Perform quantification using FSL BASIL

% FORMAT: [CBF_nocalib] = xASL_quant_Basil(PWI, x)
% 
% INPUT:
%   PWI             - image matrix of perfusion-weighted image (REQUIRED)
%   x               - struct containing pipeline environment parameters (REQUIRED)
%
% OUTPUT:
% CBF_nocalib       - Quantified CBF image
%                     (if there is no FSL/BASIL installed, we return the original PWI)
% resultFSL         - describes if the execution was successful
%                     (0 = successful, NaN = no FSL/BASIL found, 1 or other = something failed)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script performs quantification of the PWI using the FSL Basil pipeline. Final calibration to
%              physiological units is performed by dividing the quantified PWI by the M0 image/value.
%              This function performs the following steps:
%
% 1. Define paths
% 2. Delete previous BASIL output
% 3. Write the PWI as Nifti file for Basil to read as input
% 4. Create option_file that contains options which are passed to Fabber
% 5. Run Basil and retrieve CBF output
% 6. Scaling to physiological units
% 7. Householding
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: CBF_nocalib = xASL_quant_Basil(PWI, x);
%
% __________________________________
% Copyright 2015-2021 ExploreASL 
    

    %% Admin
    fprintf('%s\n','Quantification CBF using FSL Basil:');    
    
    %% 1. Define paths
    pathBasilInput = fullfile(x.dir.SESSIONDIR, 'PWI4D_BasilInput.nii');
    pathBasilOptions = fullfile(x.dir.SESSIONDIR, 'Basil_ModelOptions.txt');
    dirBasilOutput = fullfile(x.dir.SESSIONDIR, 'BasilOutput');
    
    
    %% 2. Delete previous BASIL output
    xASL_adm_DeleteFileList(x.dir.SESSIONDIR, '(?i)^basilOutput.*$', 1, [0 Inf]);
    FolderList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(?i)^basilOutput.*$', 'FPList', [0 Inf], 1);
    for iFolder=1:numel(FolderList)
        xASL_delete(FolderList{iFolder}, 1);
    end
    fprintf('%s\n', 'Note that any file not found warnings can be ignored, this pertains to the use of symbolic links by BASIL');
    
    % Remove residual BASIL-related files
    xASL_delete(pathBasilOptions);
    xASL_delete(pathBasilInput);
    xASL_adm_DeleteFileList(x.dir.SESSIONDIR, '(?i)^.*basil.*$', 1, [0 Inf]);
    
    %% 3. Write the PWI as Nifti file for Basil to read as input
    % FIXME would be good to have a brain mask at this point -> PM: if this would be a brainmask as well, we can skip creating a dummy input image here
    PWI(isnan(PWI)) = 0;
    
    if size(PWI,4) == 1 
        % singlePLD
        xASL_io_SaveNifti(x.P.Path_PWI, pathBasilInput, PWI, [], 0); % use PWI path
    else
        % multiPLD
        xASL_io_SaveNifti(x.P.Path_PWI4D, pathBasilInput, PWI, [], 0); % use PWI4D path
    end


    %% 4. Create option_file that contains options which are passed to Fabber
    % basil_options is a character array containing CLI args for the Basil command
    if size(PWI,4) < 2
        bMultiPLD = false; 
    else
        bMultiPLD = true;
    end
    
    BasilOptions = xASL_quant_Basil_Options(pathBasilOptions, x, PWI, bMultiPLD);
    
    %% 5. Run Basil and retrieve CBF output
    % args.bAutomaticallyDetectFSL=1;
    [~, resultFSL] = xASL_fsl_RunFSL(['basil -i ' xASL_adm_UnixPath(pathBasilInput) ' -@ ' xASL_adm_UnixPath(pathBasilOptions) ' -o ' xASL_adm_UnixPath(dirBasilOutput) ' ' BasilOptions], x);
    
    % Check if FSL failed
    if isnan(resultFSL)
        error('FSL BASIL was not found, exiting...');
    elseif resultFSL~=0
        error('Something went wrong running FSL BASIL...');       
    end
    
    fprintf('%s\n', 'The following warning (if mentioned above) can be ignored:');
    fprintf('%s\n', '/.../fsl/bin/basil: line 124: imcp: command not found');
    
    pathBasilMean = xASL_adm_GetFileList(dirBasilOutput, '^mean_ftiss\.nii$', 'FPListRec');
    pathBasilMean = pathBasilMean{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
       
    CBF_nocalib = xASL_io_Nifti2Im(pathBasilMean);
        
    %% 6. Scaling to physiological units
    % Note different to xASL_quant_SinglePLD since Fabber has T1 in seconds
    % and does not take into account labeling efficiency
    
    CBF_nocalib = CBF_nocalib .* 6000 .* x.Q.Lambda ./ x.Q.LabelingEfficiency;
    % (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
    % & old Siemens sequence also didn't need the 1 gr->100 gr conversion

    
    %% 7. Householding
    xASL_delete(pathBasilInput);

    %% Output of BASIL
    % Basils Output is in the subfolder '/BasilOutput' which contains
    % multiple steps if there are multiple iterations, and always contains
    % a symbolic link (symlink) to the foldername of the latest
    % iteration/step ('stepX_latest').
    %
    % When we know what we want with this output, we can remove this
    % symbolic link or some folders/output, using xASL_delete.
    % 
    % delete, xASL_delete & xASL_adm_DeleteFileList are able to delete symbolic
    % links as well (in Unix-systems)
    
end




function [BasilOptions] = xASL_quant_Basil_Options(pathBasilOptions, x, PWI, bMultiPLD)
% Save a Basil options file and store CLI options for Basil

    %% Create option_file that contains options which are passed to Fabber
    % basil_options is a character array containing CLI args for the Basil command
    
    FIDoptionFile = fopen(pathBasilOptions, 'w+');
    BasilOptions = '';

    %% Basic acquisition and tissue parameters
     % create CSV style text for TIs input

     % print parameters
     fprintf(FIDoptionFile, '# Basil options written by ExploreASL\n');
     fprintf(FIDoptionFile, '--iaf=diff\n'); % as input is PWI
     
     TIs = (x.Q.LabelingDuration + x.Q.Initial_PLD)'/1000;
     PLDAmount = length(TIs);
     for TIsingle = 1:PLDAmount
         fprintf(FIDoptionFile, ['--ti%d=%.2f\n'], TIsingle, TIs(TIsingle));
     end
     fprintf(FIDoptionFile, '--repeats=%i\n', size(PWI, 4)/PLDAmount);
     fprintf(FIDoptionFile, '--t1b=%f\n', x.Q.BloodT1/1000);
     fprintf(FIDoptionFile, '--tau=%f\n', x.Q.LabelingDuration/1000); % FIXME tau could be list
     fprintf(FIDoptionFile, '--slicedt=%f\n', x.Q.BasilSliceReadoutTime/1000);
     fprintf(FIDoptionFile, '--save-model-fit\n');


    %% FIXME Aquisition options we might be able to use in the future
    %fprintf(option_file, '--FA=%f\n', fa);
    %fprintf('Basil: Flip angle for look-locker readout: %f\n', fa);
    %fprintf(option_file, '--sliceband=%i\n', sliceband);
    %fprintf('Basil: Multi-band setup with number of slices per band: %i\n', slicedband);

    % This helps avoid failure on the structural-space image
    fprintf(FIDoptionFile, '--allow-bad-voxels\n');

    %% FIXME is a user-specified T1 map possible in ExploreASL?
    %fprintf(option_file, '--t1im=%s\n', t1im)
    %fprintf('Basil: Using supplied T1 (tissue) image in BASIL: %s\n', $t1im)

    % Labelling type - PASL or pCASL
    switch lower(x.Q.LabelingType)
        case 'pasl'
            fprintf('Basil: PASL model\n');
        case 'casl'
            fprintf(FIDoptionFile, '--casl\n');
            fprintf('Basil: CASL/PcASL model\n');
    end


    %% Model option - 1 or 2 compartment
    % The 1-compartment model is what we call 'White paper mode' in oxford-asl.
    % This means zero ATT (all bolus delivered by imaging time) and blood T1 only
    % The 2-compartment model is the 'standard' buxton model which takes into account
    % the ATT and tissue T1 value.

    switch x.Q.nCompartments
        case 1
            fprintf(FIDoptionFile, '--bat=0\n');
            fprintf(FIDoptionFile, '--t1=%f\n', x.Q.BloodT1/1000);
            fprintf('Basil: Single-compartment (white paper mode)\n');
        case 2
            fprintf(FIDoptionFile, '--bat=%f\n', x.Q.ATT/1000);
            fprintf(FIDoptionFile, '--t1=%f\n', x.Q.TissueT1/1000);
            fprintf('Basil: 2-compartment - ATT=%fs\n', x.Q.ATT/1000);
    end

    
    %%  Bolus duration\labeling duration typically fixed but can be inferred
    if isfield(x.Q,'BasilInferTau') && x.Q.BasilInferTau
         fprintf(FIDoptionFile, '--infertau\n');
         fprintf('Basil: Infer bolus duration component\n')
    else
        fprintf('Basil: Fixed bolus duration component\n')
    end

    
    %% ATT and arterial component inference only possible with multi-PLD
    if ~bMultiPLD
        fprintf('Basil: Single-delay data - cannot infer ATT or arterial component\n');
        x.Q.BasilInferATT = 0;
        x.Q.BasilInferArt = 0;
    end


    %% Infer arterial transit time
    if isfield(x.Q,'BasilInferATT') && x.Q.BasilInferATT
        if ~isfield(x.Q, 'BasilATTSD')
            x.Q.BasilATTSD = 1.0;
        end
        fprintf(FIDoptionFile, '--inferbat\n');
        fprintf(FIDoptionFile, '--batsd=%f\n', x.Q.BasilATTSD);
	    fprintf('Basil: Setting std dev of the (tissue) BAT prior std.dev. to %f\n', x.Q.BasilATTSD);
    else
        BasilOptions = [BasilOptions ' --fixbat'];
        fprintf('Basil: Fixed arterial arrival time\n');
    end

    
    %% Infer arterial component
    if isfield(x.Q,'BasilInferATT') && x.Q.BasilInferATT
        fprintf(FIDoptionFile, '--inferart\n');
        fprintf('Basil: Infer arterial component');
        fprintf('Basil: Variable arterial component arrival time');
    end


    %% Noise specification
    % For small numbers of time points we need an informative noise prior. 
    % The user can specify an assumed SNR for this, or give noise standard deviation directly.
    if ~isfield(x.Q,'BasilSNR') || ~x.Q.BasilSNR
        x.Q.BasilSNR = 10;
    end

    if size(PWI, 4) < 5
        x.Q.BasilNoisePrior = 1;
        fprintf('Basil: Small number of volumes (%i < 5): informative noise prior will be used\n', size(PWI, 4));
    end

    if isfield(x.Q,'BasilNoisePrior') && x.Q.BasilNoisePrior
        % Use an informative noise prior
        if ~isfield(x.Q,'BasilNoiseSD') || ~x.Q.BasilNoiseSD
            fprintf('Basil: Using SNR of %f to set noise std dev\n', x.Q.BasilSNR);
            % Estimate signal magntiude FIXME brain mask assume half of voxels
            mag_max = max(PWI, [], 4);
            brain_mag = 2*xASL_stat_MeanNan(mag_max(:));
            fprintf('Basil: Mean maximum signal across brain: %f\n', brain_mag);
            % This will correspond to whole brain CBF (roughly) - about 0.5 of GM
            x.Q.BasilNoiseSD = sqrt(brain_mag * 2 / x.Q.BasilSNR);
        end
        fprintf('Basil: Using a prior noise std.dev. of: %f\n', x.Q.BasilNoiseSD);
        fprintf(FIDoptionFile, '--prior-noise-stddev=%f\n', x.Q.BasilNoiseSD);
    end

    
    %% Miscellaneous features
    if isfield(x,'BasilSpatial') && x.Q.BasilSpatial
        fprintf('Basil: Instructing BASIL to use automated spatial smoothing\n');
        BasilOptions = [BasilOptions ' --spatial'];
    end

    if isfield(x.Q,'BasilInferT1') && x.Q.BasilInferT1
        fprintf(FIDoptionFile, '--infert1\n');
        fprintf('Basil: Instructing BASIL to infer variable T1 values\n');
    end

    if isfield(x.Q,'BasilExch')
        fprintf('Basil: Using exchange model: %s\n', x.Q.BasilExch);
        fprintf(FIDoptionFile, '--exch=%s\n', x.Q.BasilExch);
    end

    if isfield(x.Q,'BasilDisp')
        fprintf('Basil: Using dispersion model: %s\n', x.Q.BasilDisp);
        fprintf(FIDoptionFile, '--disp=%s\n', x.Q.BasilDisp);
    end

    if isfield(x.Q,'BasilDebug') && x.Q.BasilDebug
        BasilOptions = [BasilOptions ' --devel'];
    end
    
    
    %% Save Basil options file
    fclose(FIDoptionFile);
    
end


