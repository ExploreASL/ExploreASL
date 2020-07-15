function xASL_par_Fix(DataParFile,pathASL4D)
%xASL_par_Fix Script which tries to handle some missing parameters.
%
% FORMAT:       xASL_par_Fix(pathASL4D,pathM0);
%
% INPUT:        DataParFile - Path to DataParFile
%               pathASL4D - Path to ASL4D JSON file
%
% OUTPUT:       n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Warning: This function tries to add missing values which it
%               does not now. We use simplified scenarios to make the
%               ExploreASL pipeline run nevertheless.
%               Right now only the Q field is added and the LabelingType
%               PCASL is set.
%
% EXAMPLE:      xASL_par_Fix('.../DataParFile.json','.../ASL4D.json');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Parameters
fprintf('Check parameters...\n');

% Checkout the ExploreASL/DataParFile.m for detailed explanations
                
%% Get vendor
try
    % Read JSON file
    if xASL_exist(pathASL4D,'file')
        val = jsondecode(fileread(pathASL4D));
        % Get vendor
        if ~isfield(val,'Vendor')
            Vendor = val.Manufacturer;
        end
        % SliceReadoutTime only necessary for 2D datasets -> Get AcquisitionType
        if isfield(val,'MRAcquisitionType')
            Acquisition = val.MRAcquisitionType;
        end
    end
catch
    fprintf('Something went wrong...\n');
end


%% Fix DataParFile

% Open DataParFile
data = jsondecode(fileread(DataParFile));

% Add Q field
if ~isfield(data.x,'Q'), data.x.Q = struct; end

switch Vendor
    case 'Siemens'
        % Define temporarily static parameters
        data.x.M0 = 'separate_scan';
        data.x.Sequence = '3D_GRASE';
        if ~isfield(data.x.Q,'LabelingType'),      data.x.Q.LabelingType = "PCASL";       end
        if ~isfield(data.x.Q,'BackGrSupprPulses'), data.x.Q.BackGrSupprPulses = 5;        end
        if ~isfield(data.x.Q,'LabelingDuration'),  data.x.Q.LabelingDuration = 1800;      end
        if ~isfield(data.x.Q,'Initial_PLD'),       data.x.Q.Initial_PLD = 2000;           end
    case 'Philips'
        % Define temporarily static parameters
        data.x.M0 = 'separate_scan';
        data.x.Sequence = '2D_EPI';
        if ~isfield(data.x.Q,'LabelingType'),      data.x.Q.LabelingType = "PCASL";       end
        if ~isfield(data.x.Q,'BackGrSupprPulses'), data.x.Q.BackGrSupprPulses = 2;        end
        if ~isfield(data.x.Q,'LabelingDuration'),  data.x.Q.LabelingDuration = 1800;      end
        if ~isfield(data.x.Q,'Initial_PLD'),       data.x.Q.Initial_PLD = 2000;           end
    case 'GE'
        % Define temporarily static parameters
        data.x.M0 = 'UseControlAsM0';
        data.x.Sequence = '3D_spiral';
        if ~isfield(data.x.Q,'LabelingType'),      data.x.Q.LabelingType = "PCASL";       end
        if ~isfield(data.x.Q,'BackGrSupprPulses'), data.x.Q.BackGrSupprPulses = 5;        end
        if ~isfield(data.x.Q,'LabelingDuration'),  data.x.Q.LabelingDuration = 2025;      end
        if ~isfield(data.x.Q,'Initial_PLD'),       data.x.Q.Initial_PLD = 1450;           end
end

if strcmp(Acquisition,"2D")
    if ~isfield(data.x.Q,'SliceReadoutTime'),  data.x.Q.SliceReadoutTime = 40; end
end

% Save modified JSON file
txt = jsonencode(data);
fID = fopen(DataParFile,'w');
fwrite(fID, txt, 'char');
fclose(fID);





