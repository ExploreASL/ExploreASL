function phoenixProtocol = xASL_bids_PhoenixProtocolReader(pathDataset)
%xASL_bids_PhoenixProtocolReader Function that reads raw DICOM data and extracts the phoenix protocol parameters.
%
% FORMAT: phoenixProtocol = xASL_bids_PhoenixProtocolReader(pathDataset);
%
% INPUT:
%        pathDataset        - path to DICOM dataset (REQUIRED)
%
% OUTPUT:
%        phoenixProtocol    - structure containing the phoenix protocol parameters
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Function that reads raw DICOM data and extracts the phoenix protocol parameters.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          pathDataset = '...\test-data';
%                   phoenixProtocol = xASL_bids_PhoenixProtocolReader(pathDataset);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL


    %% Input Check

    % Check if both root folders are valid char arrays or strings
    if ~(ischar(pathDataset) || isstring(pathDataset))
        error('Input path is neither a char array not a string...');
    end

    % Default value for bPrintReport
    if nargin < 1
        error('Missing input parameters...');
    end


    %% Defaults

    % Initialize phoenix protocol
    phoenixProtocol = struct;


    %% Read DICOM file


    %% Extract parameters

    % function hdr = spm_dicom_headers(P, essentials)


end






