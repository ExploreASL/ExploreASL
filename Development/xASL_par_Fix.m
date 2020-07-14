function xASL_par_Fix(pathASL4D,pathM0)
%xASL_par_Fix Script which tries to figure out the missing fields of the
% JSON files and adds a few defaults.
%
% FORMAT:       xASL_par_Fix(pathASL4D,pathM0);
%
% INPUT:        pathASL4D - Path to ASL4D JSON file
%               pathM0 - Path to M0 JSON file
%
% OUTPUT:       n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Warning: This function tries to add missing values which it
%               does not now. We use simplified scenarios to make the
%               ExploreASL pipeline run nevertheless.
%
% EXAMPLE:      xASL_par_Fix('.../ASL4D.json','.../M0.json');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Parameters
fprintf('Check parameters...\n');

% Parameters (from TestDataSet)
parametersASL4D = [ "Modality","MagneticFieldStrength","ImagingFrequency","Manufacturer","ManufacturersModelName","SoftwareVersions",...
                    "MRAcquisitionType","ScanningSequence","SequenceVariant","ScanOptions","ImageType","SeriesNumber","AcquisitionTime",...
                    "AcquisitionNumber","SliceThickness","SpacingBetweenSlices","SAR","EchoTime","RepetitionTime","FlipAngle",...
                    "CoilString","PartialFourier","PercentPhaseFOV","EchoTrainLength","PhaseEncodingSteps","AcquisitionMatrixPE",...
                    "ReconMatrixPE","PixelBandwidth","PhaseEncodingAxis","ImageOrientationPatientDICOM","InPlanePhaseEncodingDirectionDICOM",...
                    "ConversionSoftware","ConversionSoftwareVersion","M0","readout_dim","Vendor","bReproTesting",...
                    "Q",...
                    "Q.BackGrSupprPulses","Q.LabelingType","Q.Initial_PLD","Q.LabelingDuration","Q.SliceReadoutTime"]';

% Parameters (from TestDataSet)
parametersM0 = [    "Modality","MagneticFieldStrength","ImagingFrequency","Manufacturer","ManufacturersModelName","SoftwareVersions",...
                    "MRAcquisitionType","ScanningSequence","SequenceVariant","ScanOptions","AcquisitionTime",...
                    "SliceThickness","SpacingBetweenSlices","SAR","EchoTime","RepetitionTime","FlipAngle","CoilString",...
                    "PartialFourier","PercentPhaseFOV","EchoTrainLength","PhaseEncodingSteps","AcquisitionMatrixPE",...
                    "ReconMatrixPE","PixelBandwidth","PhaseEncodingAxis","ImageOrientationPatientDICOM","InPlanePhaseEncodingDirectionDICOM",...
                    "ConversionSoftware","ConversionSoftwareVersion","M0","readout_dim","Vendor","bReproTesting",...
                    "Q",...
                    "Q.BackGrSupprPulses","Q.LabelingType","Q.Initial_PLD","Q.LabelingDuration","Q.SliceReadoutTime"]';
                
%% Fix ASL JSON file
try
    if xASL_exist(pathASL4D,'file')
        % Read JSON file
        valASL = jsondecode(fileread(pathASL4D));
        % If there is no Q then add it directly
        if ~isfield(valASL,'Q')
            fprintf('Added Q field...\n');
            valASL.Q = struct;
        end
        % If there is no M0 then add "separate_scan" option as a default
        if ~isfield(valASL,'M0')
            fprintf('Added M0 = separate_scan...\n');
            valASL.M0 = "separate_scan";
        end
        % If there is no Vendor but a Manufacturer set the Manufacturer as a default
        if ~isfield(valASL,'Vendor') && isfield(valASL,'Manufacturer')
            fprintf('Added Vendor = Manufacturer...\n');
            valASL.Vendor = valASL.Manufacturer;
        end
        % If there is no readout_dim set 2D as a default
        if ~isfield(valASL,'readout_dim')
            fprintf('Added readout_dim = 2D...\n');
            valASL.readout_dim = "2D";
        end
        % Iterate through all other parameters
        for i=1:length(parametersASL4D)
            % Q fields
            if contains(parametersASL4D(i,1),'Q.') && isfield(valASL,'Q')
                if ~isfield(valASL.Q,strrep(parametersASL4D(i,1),'Q.',''))
                    fprintf('Parameter %s is missing...\n', parametersASL4D(i,1)); 
                end
            else % Other fields
                if ~isfield(valASL,parametersASL4D(i,1))
                    fprintf('Parameter %s is missing...\n', parametersASL4D(i,1)); 
                end
            end        
        end
    end

    % Save modified JSON file
    txtASL = jsonencode(valASL);
    fID = fopen(pathASL4D,'w');
    fwrite(fID, txtASL, 'char');
    fclose(fID);
catch
   fprintf('Something went wrong trying to fix the ASL4D JSON file...\n'); 
end

%% Fix M0 JSON file
try
    if xASL_exist(pathM0,'file')
        % Read JSON file
        valM0 = jsondecode(fileread(pathM0));
        % If there is no Q then add it directly
        if ~isfield(valM0,'Q')
            fprintf('Added Q field...\n');
            valM0.Q = struct;
        end
        % If there is no M0 then add "separate_scan" option as a default
        if ~isfield(valM0,'M0')
            fprintf('Added M0 = separate_scan...\n');
            valM0.M0 = "separate_scan";
        end
        % If there is no Vendor but a Manufacturer set the Manufacturer as a default
        if ~isfield(valM0,'Vendor') && isfield(valM0,'Manufacturer')
            fprintf('Added Vendor = Manufacturer...\n');
            valM0.Vendor = valM0.Manufacturer;
        end
        % If there is no readout_dim set 2D as a default
        if ~isfield(valM0,'readout_dim')
            fprintf('Added readout_dim = 2D...\n');
            valM0.readout_dim = "2D";
        end
        % Iterate through all other parameters
        for i=1:length(parametersM0)
            % Q fields
            if contains(parametersM0(i,1),'Q.') && isfield(valASL,'Q')
                if ~isfield(valM0.Q,strrep(parametersM0(i,1),'Q.',''))
                    fprintf('Parameter %s is missing...\n', parametersM0(i,1)); 
                end
            else % Other fields
                if ~isfield(valM0,parametersM0(i,1))
                    fprintf('Parameter %s is missing...\n', parametersM0(i,1)); 
                end
            end
        end
    end

    % Save modified JSON file
    txtM0 = jsonencode(valM0);
    fID = fopen(pathM0,'w');
    fwrite(fID, txtM0, 'char');
    fclose(fID);
catch
   fprintf('Something went wrong trying to fix the M0 JSON file...\n'); 
end



