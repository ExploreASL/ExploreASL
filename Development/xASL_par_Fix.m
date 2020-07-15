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
for i=1:2
    if i==1, path=pathASL4D; end
    if i==2, path=pathM0;    end
    
    try
        if xASL_exist(path,'file')
            % Read JSON file
            val = jsondecode(fileread(path));

            % Add missing fields if they are not defined
            if ~isfield(val,'Vendor'),              val.Vendor = val.Manufacturer;  end
            if ~isfield(val,'Q'),                   val.Q = struct;                 end

            % Get Manufacturer
            if strcmp(val.Manufacturer,'Siemens')
                if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
                if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 5;        end
                if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 1800;      end
                if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 2000;           end
            elseif strcmp(val.Vendor,'Philips')
                if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
                if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 2;        end
                if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 1800;      end
                if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 2000;           end
            elseif strcmp(val.Vendor,'GE')
                if ~isfield(val.Q,'LabelingType'),      val.Q.LabelingType = "PCASL";       end
                if ~isfield(val.Q,'BackGrSupprPulses'), val.Q.BackGrSupprPulses = 5;        end
                if ~isfield(val.Q,'LabelingDuration'),  val.Q.LabelingDuration = 2025;      end
                if ~isfield(val.Q,'Initial_PLD'),       val.Q.Initial_PLD = 1450;           end
            end

            % SliceReadoutTime only necessary for 2D datasets
            if isfield(val,'MRAcquisitionType')
                if strcmp(val.MRAcquisitionType,"2D")
                    if ~isfield(val.Q,'SliceReadoutTime'),  val.Q.SliceReadoutTime = 40;    end
                end
            end
        end

        % Save modified JSON file
        txt = jsonencode(val);
        fID = fopen(path,'w');
        fwrite(fID, txt, 'char');
        fclose(fID);
    catch
       fprintf('Something went wrong trying to fix the JSON file...\n'); 
    end
end






