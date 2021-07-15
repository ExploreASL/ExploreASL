function x = xASL_io_CheckDeprecatedFieldsX(x, bVerbose)
% xASL_io_CheckDeprecatedFieldsX Check deprecated fields of x and fix them based on a conversion table
%
% FORMAT:   x = xASL_io_CheckDeprecatedFieldsX(x)
%
% INPUT:
%   x               - x structure
%   bVerbose        - Print extra information
%
% OUTPUT:
%   x               - x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check deprecated fields of x and fix them based on a conversion table.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:      n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL

    if nargin<2 || isempty(bVerbose)
        bVerbose = false;
    end

    % Get conversion table
    conversionTable = xASL_adm_GetDeprecatedFields();
    
    nOutdatedParameter = 0;
    bOldFieldsDetected = 0;
    for iField = 1:size(conversionTable,1)
        if isfield(x,conversionTable{iField,1})
            bOldFieldsDetected = 1;
            if ~isfield(x,conversionTable{iField,2}) || ...
                    ~isfield(x.(conversionTable{iField,2}),conversionTable{iField,3}) || ...
                    ~isfield(x.(conversionTable{iField,2}),conversionTable{iField,4})
                % Debugging
                if bVerbose
                    fprintf('%s %s %s\n',conversionTable{iField,2},conversionTable{iField,3},conversionTable{iField,4});
                end
                % Deprecated field detected
                nOutdatedParameter = nOutdatedParameter+1;
                detectedFields{nOutdatedParameter,1} = conversionTable{iField,1};
                % Basic field does not exist
                if ~isfield(x,conversionTable{iField,2})
                    x.(conversionTable{iField,2}) = struct;
                end
                % Actual deprecated field
                if ~isfield(x.(conversionTable{iField,2}),conversionTable{iField,3})
                    % Sub-Structure (x.FIELD -> x.STRUCT.FIELD)
                    x.(conversionTable{iField,2}).(conversionTable{iField,3}) = x.(conversionTable{iField,1});
                elseif ~isfield(x.(conversionTable{iField,2}).(conversionTable{iField,3}),conversionTable{iField,4}) ...
                        && isstruct(x.(conversionTable{iField,2}).(conversionTable{iField,3}))
                    % Sub-Sub-Structure (x.FIELD -> x.STRUCT.STRUCT.FIELD)
                    x.(conversionTable{iField,2}).(conversionTable{iField,3}).(conversionTable{iField,4}) = x.(conversionTable{iField,1});
                else
                    if bVerbose
                        fprintf('%s.%s.%s already exists...\n',conversionTable{iField,2},conversionTable{iField,3},conversionTable{iField,4});
                    end
                end
            end
            % Remove old field
            x = rmfield(x,conversionTable{iField,1});
        end
    end
    
    if bOldFieldsDetected && bVerbose
        % Print warning and individual fields
        fprintf('Detected deprecated fields in dataPar.json ...\n');
        for iField = 1:size(detectedFields,1)
            fprintf('Deprecated field: %s\n', detectedFields{iField,1});
        end
    end

end


