function xASL_par_Fix(pathASL4D,pathM0)
%xASL_par_Fix Script which tries to handle some missing parameters.
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
%               Right now only the Q field is added and the LabelingType
%               PCASL is set.
%
% EXAMPLE:      xASL_par_Fix('.../ASL4D.json','.../M0.json');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Parameters
fprintf('Check parameters...\n');

% % Parameters (from TestDataSet)
% parametersASL4D = [ "Modality","MagneticFieldStrength","ImagingFrequency","Manufacturer","ManufacturersModelName","SoftwareVersions",...
%                     "MRAcquisitionType","ScanningSequence","SequenceVariant","ScanOptions","ImageType","SeriesNumber","AcquisitionTime",...
%                     "AcquisitionNumber","SliceThickness","SpacingBetweenSlices","SAR","EchoTime","RepetitionTime","FlipAngle",...
%                     "CoilString","PartialFourier","PercentPhaseFOV","EchoTrainLength","PhaseEncodingSteps","AcquisitionMatrixPE",...
%                     "ReconMatrixPE","PixelBandwidth","PhaseEncodingAxis","ImageOrientationPatientDICOM","InPlanePhaseEncodingDirectionDICOM",...
%                     "ConversionSoftware","ConversionSoftwareVersion","M0","readout_dim","Vendor","bReproTesting",...
%                     "Q",...
%                     "Q.BackGrSupprPulses","Q.LabelingType","Q.Initial_PLD","Q.LabelingDuration","Q.SliceReadoutTime"]';
% 
% % Parameters (from TestDataSet)
% parametersM0 = [    "Modality","MagneticFieldStrength","ImagingFrequency","Manufacturer","ManufacturersModelName","SoftwareVersions",...
%                     "MRAcquisitionType","ScanningSequence","SequenceVariant","ScanOptions","AcquisitionTime",...
%                     "SliceThickness","SpacingBetweenSlices","SAR","EchoTime","RepetitionTime","FlipAngle","CoilString",...
%                     "PartialFourier","PercentPhaseFOV","EchoTrainLength","PhaseEncodingSteps","AcquisitionMatrixPE",...
%                     "ReconMatrixPE","PixelBandwidth","PhaseEncodingAxis","ImageOrientationPatientDICOM","InPlanePhaseEncodingDirectionDICOM",...
%                     "ConversionSoftware","ConversionSoftwareVersion","M0","readout_dim","Vendor","bReproTesting",...
%                     "Q",...
%                     "Q.BackGrSupprPulses","Q.LabelingType","Q.Initial_PLD","Q.LabelingDuration","Q.SliceReadoutTime"]';
                
%% Fix ASL JSON file
try
    if xASL_exist(pathASL4D,'file')
        % Read JSON file
        val = jsondecode(fileread(pathASL4D));

        % Add missing fields if they are not defined
        if ~isfield(val,'Vendor'),              val.Vendor = val.Manufacturer;  end
        if ~isfield(val,'Q'),                   val.Q = struct;                 end
        if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";   end
        
        % Get Manufacturer
        if strcmp(val.Manufacturer,'Siemens')
            if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
            if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 5;        end % = BackgroundSuppressionNumberPulses?
            if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 1800;      end
            if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 2000;           end % = PostLabelingDelay?
            if ~isfield(val.Q,'SliceReadoutTime'),  val.Q.SliceReadoutTime = 40;        end % = LabelingDistance?
        elseif strcmp(val.Vendor,'Philips')
            if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
            if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 2;        end % = BackgroundSuppressionNumberPulses?
            if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 1800;      end
            if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 2000;           end % = PostLabelingDelay?
            if ~isfield(val.Q,'SliceReadoutTime'),  val.Q.SliceReadoutTime = 40;        end % = LabelingDistance?
        elseif strcmp(val.Vendor,'GE')
            if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
            if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 5;        end % = BackgroundSuppressionNumberPulses?
            if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 2025;      end
            if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 1450;           end % = PostLabelingDelay?
            if ~isfield(val.Q,'SliceReadoutTime'),  val.Q.SliceReadoutTime = 40;        end % = LabelingDistance?
        end
    end

    % Save modified JSON file
    txtASL = jsonencode(val);
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
        val = jsondecode(fileread(pathM0));

        % Add missing fields if they are not defined
        if ~isfield(val,'Vendor'),              val.Vendor = val.Manufacturer;  end
        if ~isfield(val,'Q'),                   val.Q = struct;                 end
        if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";   end
        
        % Get Manufacturer
        if strcmp(val.Manufacturer,'Siemens')
            if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
            if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 5;        end % = BackgroundSuppressionNumberPulses?
            if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 1800;      end
            if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 2000;           end % = PostLabelingDelay?
            if ~isfield(val.Q,'SliceReadoutTime'),  val.Q.SliceReadoutTime = 40;        end % = LabelingDistance?
        elseif strcmp(val.Vendor,'Philips')
            if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
            if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 2;        end % = BackgroundSuppressionNumberPulses?
            if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 1800;      end
            if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 2000;           end % = PostLabelingDelay?
            if ~isfield(val.Q,'SliceReadoutTime'),  val.Q.SliceReadoutTime = 40;        end % = LabelingDistance?
        elseif strcmp(val.Vendor,'GE')
            if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
            if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 5;        end % = BackgroundSuppressionNumberPulses?
            if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 2025;      end
            if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 1450;           end % = PostLabelingDelay?
            if ~isfield(val.Q,'SliceReadoutTime'),  val.Q.SliceReadoutTime = 40;        end % = LabelingDistance?
        end
    end

    % Save modified JSON file
    txtM0 = jsonencode(val);
    fID = fopen(pathM0,'w');
    fwrite(fID, txtM0, 'char');
    fclose(fID);
catch
   fprintf('Something went wrong trying to fix the M0 JSON file...\n'); 
end



