function [xasl,parameters,parameterList,phoenixProtocol] = xASL_bids_GetPhoenixProtocol(pathData,bUseDCMTK)
%xASL_bids_GetPhoenixProtocol Wrapper function for xASL_bids_PhoenixProtocolReader and xASL_bids_PhoenixProtocolAnalyzer
%
% FORMAT: [xasl,parameters,parameterList,phoenixProtocol] = xASL_bids_GetPhoenixProtocol(pathData,bUseDCMTK)
%
% INPUT:
%        pathData           - Path to the DICOM dataset [CHAR ARRAY] (REQUIRED)
%                             Please insert a single path as a char array (not a string). The file extension should be
%                             either ".dcm" or ".IMA". The script is going to check this.
%        bUseDCMTK          - Use DCMTK toolbox to get the phoenix protocol (OPTIONAL, DEFAULT: true)
%                             If you select "true", then the DICOM tag will be read using the DCMTK toolbox.
%                             Alternatively, you can select "false" to read the DICOM tag using python.
%                             We highly recommend using the default option.
%
% OUTPUT:
%        xasl               - list of xasl parameters
%        parameters         - list of phoenix parameters you want to know
%        parameterList      - List of parameters from the reduced phoenix protocol. This parameter list is created to 
%                             make it easier to analyze and interpret the individual phoenix protocol parameters.
%                             You can alternatively check out the actual protocol as a cell array or char array.
%        phoenixProtocol    - Cell array containing the Siemens phoenix protocol [CELL ARRAY]
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Function that reads raw DICOM data (".dcm" or ".IMA") and extracts the phoenix protocol parameters.
%                   Only works for Siemens DICOM data with phoenix protocol (tag = [0x29,0x1020]).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          % Define the path to your DICOM file and use DCMTK
%                   pathData = 'example.dcm';
%                   bUseDCMTK = true;
%                   % Run the workflow
%                   [xasl,parameters,parameterList,phoenixProtocol] = xASL_bids_GetPhoenixProtocol(pathData,bUseDCMTK);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL

    %% Input Check
    
    % Check number of input parameters
    if nargin < 1
        error('Missing input parameters...');
    end

    % Check bUseDCMTK
    if nargin < 2
        bUseDCMTK = true;
    end

    % Check if the input path is a character array or a string
    if ~(ischar(pathData))
        error('Input path is neither a char array not a string...');
    end
    
    % Get file parts
    [~,~,fileExtension] = fileparts(pathData);
    
    % Check extension
    if ~strcmp(fileExtension,'.dcm') && ~strcmp(fileExtension,'.IMA')
        fprintf('It is possible that the input path is not a DICOM file...\n');
    end
    
    %% Run phoenix functions
    rawPhoenixProtocol = parsePhoenixProtocol(pathData,bUseDCMTK);
    [parameterList,phoenixProtocol] = xASL_bids_PhoenixProtocolReader(rawPhoenixProtocol);
    [xasl,parameters] = xASL_bids_PhoenixProtocolAnalyzer(parameterList);

end

% Function to parse the phoenix protocol
function rawPhoenixProtocol = parsePhoenixProtocol(pathData,bUseDCMTK)
    rawPhoenixProtocol = ''; % Fallback
    if bUseDCMTK
        headerDCMTK = xASL_io_DcmtkRead(pathData);
        % Check if protocol exists
        if isfield(headerDCMTK,'PhoenixProtocol') && ~isempty(headerDCMTK.PhoenixProtocol)
            phoenixProtocol = headerDCMTK.PhoenixProtocol;
            rawPhoenixProtocol = phoenixProtocol;                   % Raw phoenix protocol as a char array
        end
    else
        py.importlib.import_module('pydicom');
        ds = py.pydicom.dcmread(pathData,false,true);
        phoenixProtocol = char(ds{2691104}.value); % 2691104 = (0x29,0x1020)
        rawPhoenixProtocol = phoenixProtocol;                       % Raw phoenix protocol as a char array
    end
end



