function [CBF_nocalib, ATT_map, resultFSL] = xASL_quant_Fabber(PWI, x)
%xASL_quant_Basil Perform quantification using FSL FABBER
%
% FORMAT: [CBF_nocalib, ATT_map, resultFSL] = xASL_quant_Fabber(PWI, x)
% 
% INPUT:
%   PWI             - image matrix of perfusion-weighted image (REQUIRED)
%   x               - struct containing pipeline environment parameters (REQUIRED)
%
% OUTPUT:
% CBF_nocalib       - Quantified CBF image
%                     (if there is no FSL/BASIL installed, we return the original PWI)
% ATT_Map           - ATT map (if possible to calculate with multi-PLD)
% resultFSL         - describes if the execution was successful
%                     (0 = successful, NaN = no FSL/BASIL found, 1 or other = something failed)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script performs quantification of the PWI using the FSL Basil pipeline. Final calibration to
%              physiological units is performed by dividing the quantified PWI by the M0 image/value.
%              This function performs the following steps:
%
% 1. Define paths
% 2. Delete previous Fabber output
% 3. Write the PWI as Nifti file for Basil to read as input
% 4. Create option_file that contains options which are passed to Fabber
% 5. Run Basil and retrieve CBF output
% 6. Scaling to physiological units
% 7. Householding
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: CBF_nocalib = xASL_quant_Fabber(PWI, x);
%
% __________________________________
% Copyright 2015-2021 ExploreASL 
    

    %% Admin
    fprintf('%s\n','Quantification CBF using FSL FABBER:');    
    
    %% 1. Define paths
    pathFabberInput = fullfile(x.dir.SESSIONDIR, 'PWI4D_BasilInput.nii'); % Fabber uses the same input as Basil
    pathFabberOptions = fullfile(x.dir.SESSIONDIR, 'Fabber_ModelOptions.txt');
    dirFabberOutput = fullfile(x.dir.SESSIONDIR, 'Fabber_outputs');
        
    %% 2. Delete previous BASIL output
    xASL_adm_DeleteFileList(x.dir.SESSIONDIR, '(?i)^FabberOutput.*$', 1, [0 Inf]);
    FolderList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '(?i)^FabberOutput.*$', 'FPList', [0 Inf], 1);
    for iFolder=1:numel(FolderList)
        xASL_delete(FolderList{iFolder}, 1);
    end
    fprintf('%s\n', 'Note that any file not found warnings can be ignored, this pertains to the use of symbolic links by BASIL');
    
    % Remove residual BASIL-related files
    xASL_delete(pathFabberOptions);
    xASL_delete(pathFabberInput);
    xASL_adm_DeleteFileList(x.dir.SESSIONDIR, '(?i)^.*Fabber.*$', 1, [0 Inf]);
    
    %% 3. Write the PWI as Nifti file for Basil to read as input
    % FIXME would be good to have a brain mask at this point -> PM: if this would be a brainmask as well, we can skip creating a dummy input image here
    PWI(isnan(PWI)) = 0;
    
    if ~x.modules.asl.bMultiPLD
        % SinglePLD
        xASL_io_SaveNifti(x.P.Path_PWI, pathFabberInput, PWI, [], 0); % use PWI path
    else
        % MultiPLD
        xASL_io_SaveNifti(x.P.Path_PWI4D, pathFabberInput, PWI, [], 0); % use PWI4D path
    end

    %% 4. Create option_file that contains options which are passed to Fabber
    % basil_options is a character array containing CLI args for the Basil command
    
    BasilOptions = xASL_quant_Fabber_Options(pathFabberOptions, x, PWI, pathFabberInput);
    
    %% 5. Run Fabber and retrieve CBF output
    % args.bAutomaticallyDetectFSL=1;
    [~, resultFSL] = xASL_fsl_RunFSL(['fabber_asl -@ ' xASL_adm_UnixPath(pathFabberOptions)], x);
    %FSLCommand = 'fabber_asl --output=Fabber_outputs --method=vb --data=PWI4D_BasilInput.nii.gz --model=asl_multite --infertexch --save-var --save-residuals --save-model-fit --noise=white --tau1=1.0 --tau2=1.0 --tau3=1.0 --ti1=1.6 --ti2=2.6 --ti3=3.6 --nte1=8 --nte2=8 --nte3=8  --te1=0.01312 --te2=0.03936 --te3=0.0656 --te4=0.09184 --te5=0.11808 --te6=0.14432 --te7=0.17056 --te8=0.1968 --te9=0.01312 --te10=0.03936 --te11=0.0656 --te12=0.09184 --te13=0.11808 --te14=0.14432 --te15=0.17056 --te16=0.1968 --te17=0.01312 --te18=0.03936 --te19=0.0656 --te20=0.09184 --te21=0.11808 --te22=0.14432 --te23=0.17056 --te24=0.1968 --inferitt';

    
    % Check if FSL failed
    if isnan(resultFSL)
        error('FSL Fabber was not found, exiting...');
    elseif resultFSL~=0
        error('Something went wrong running FSL Fabber...');       
    end
    
    fprintf('%s\n', 'The following warning (if mentioned above) can be ignored:');
    fprintf('%s\n', '/.../fsl/bin/basil: line 124: imcp: command not found');
    
    pathFabberCBF = xASL_adm_GetFileList(dirFabberOutput, '^mean_ftiss\.nii$', 'FPListRec');
	if isempty(pathFabberCBF)
		error('FSL BASIL failed');
	end
    pathFabberCBF = pathFabberCBF{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
       
    CBF_nocalib = xASL_io_Nifti2Im(pathFabberCBF);
        
	pathFabberATT = xASL_adm_GetFileList(dirFabberOutput, '^mean_delttiss\.nii$', 'FPListRec');
	if ~isempty(pathFabberATT)
		pathFabberATT = pathFabberATT{end}; % we assume the latest iteration (alphabetically) is optimal. also converting cell to char array
		ATT_map = xASL_io_Nifti2Im(pathFabberATT);
	else
		ATT_map = [];
	end
	
    %% 6. Scaling to physiological units
    % Note different to xASL_quant_SinglePLD since Fabber has T1 in seconds
    % and does not take into account labeling efficiency
    
    CBF_nocalib = CBF_nocalib .* 6000 .* x.Q.Lambda ./ x.Q.LabelingEfficiency;
    % (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
    % & old Siemens sequence also didn't need the 1 gr->100 gr conversion
    
    %% 7. Householding
    xASL_delete(pathFabberInput);

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




function [FIDoptionFile] = xASL_quant_Fabber_Options(pathFabberOptions, x, PWI, pathFabberInput)
% Save a Basil options file and store CLI options for Basil
% 1. Create an option file
% 2. Basic tissue parameters
% 3. Basic acquisition parameters
% 4. Model fiting parameters
% 5. Extra features on demand
% 6. Save BASIL options file

%% 1. Create option_file that contains options which are passed to Fabber
% basil_options is a character array containing CLI args for the Basil command

FIDoptionFile = fopen(pathFabberOptions, 'w+');

[~,PWIfileName, ext] = fileparts(pathFabberInput);
PWIfile = [PWIfileName ext];

fprintf(FIDoptionFile, '# Fabber options written by ExploreASL\n');
%fprintf(FIDoptionFile, 'fabber_asl\n');
fprintf(FIDoptionFile, '--output=Fabber_outputs\n');
fprintf(FIDoptionFile, '--method=vb\n');
fprintf(FIDoptionFile, '--data=%s\n', PWIfile);
fprintf(FIDoptionFile, '--model=asl_multite\n');
fprintf(FIDoptionFile, '--infertexch\n');
fprintf(FIDoptionFile, '--save-var\n');
fprintf(FIDoptionFile, '--save-residuals\n');
fprintf(FIDoptionFile, '--save-model-fit\n');
fprintf(FIDoptionFile, '--noise=white\n');

%% 2. Basic tissue parameters
fprintf(FIDoptionFile, '--t1b=%f\n', x.Q.BloodT1/1000);
fprintf(FIDoptionFile, '--t1=%f\n', x.Q.TissueT1/1000);

%% 3. Basic acquisition parameters

%PCASL


% Print all the PLDs and LabDurs

[PLDs, ind] = unique(x.Q.Initial_PLD);
PLDs = PLDs'/1000;

% For Time-encoded, we skip the first volume
if x.modules.asl.bTimeEncoded
    PLDs = PLDs(2:end);
end

if length(x.Q.LabelingDuration)>1
    LDs = (x.Q.LabelingDuration(ind))'/1000;
else
    LDs = ones(size(PLDs))*x.Q.LabelingDuration/1000;
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
for iLD = 1:length(LDs)
    fprintf(FIDoptionFile, '--tau%d=%.2f\n', iLD, LDs(iLD));
end

for iTE = 1:length(TEs) %We need a TE for each volume
    fprintf(FIDoptionFile, '--te%d=%.2f\n', iTE, TEs(iTE));
end

% Right now, we assume that we have averaged over PLDs
%fprintf(FIDoptionFile, '--repeats=%i\n', size(PWI, 4)/PLDAmount);
%fprintf(FIDoptionFile, '--repeats=1\n');
	
% 4. Model fiting parameters (ATT map)
fprintf(FIDoptionFile, '--inferitt');
	 


%% 6. Save Basil options file
fclose(FIDoptionFile);
    
end


