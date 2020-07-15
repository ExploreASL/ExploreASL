function xASL_par_Fix(DataParFile,pathASL4D,pathM0)
%xASL_par_Fix Script which tries to handle some missing parameters.
%
% FORMAT:       xASL_par_Fix(pathASL4D,pathM0);
%
% INPUT:        DataParFile - Path to DataParFile
%               pathASL4D - Path to ASL4D JSON file
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

% Checkout the ExploreASL/DataParFile.m for detailed explanations
                
%% Fix ASL and M= JSON files
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

%% Fix DataParFile

% Open DataParFile
data = jsondecode(fileread(DataParFile));

switch val.Manufacturer
    case 'Siemens'
        data.x.M0 = 'separate_scan';
        data.x.Sequence = '3D_GRASE';
    case 'Philips'
        data.x.M0 = 'separate_scan';
        data.x.Sequence = '2D_EPI';
    case 'GE'
        data.x.M0 = 'UseControlAsM0';
        data.x.Sequence = '3D_spiral';
end

% Save modified JSON file
txt = jsonencode(data);
fID = fopen(DataParFile,'w');
fwrite(fID, txt, 'char');
fclose(fID);





