function [s, FieldNames] = xASL_imp_AppendParmsParameters(parms)
%xASL_imp_AppendParmsParameters Append Parms Parameters.
%
% FORMAT: [s, FieldNames] = xASL_imp_AppendParmsParameters(parms)
% 
% INPUT:
%   parms      - Parameters (REQUIRED, STRUCT)
%
% OUTPUT:
%   s          - String
%   FieldNames - List of strings
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Append Parms Parameters.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [s, FieldNames] = xASL_imp_AppendParmsParameters(parms);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


    %% Append Parms Parameters

    % This function outputs s=fields of _parms.mat
    s = [];

    FieldNames = {'RepetitionTimePreparation' 'RepetitionTime' 'EchoTime' 'NumberOfAverages' 'RescaleSlope' 'RescaleSlopeOriginal'...
        'MRScaleSlope' 'RescaleIntercept' 'AcquisitionTime' 'AcquisitionMatrix' 'TotalReadoutTime'...
        'EffectiveEchoSpacing'};

    if ~isempty(parms)
        for iField=1:length(FieldNames)
            if isfield(parms,FieldNames{iField})
                s = [s ',' xASL_num2str(parms.(FieldNames{iField}))];
            else
                s = [s ',n/a'];
            end
        end
    end
    
end



