function x = xASL_io_CheckDeprecatedFieldsX(x, bVerbose)
% xASL_io_CheckDeprecatedFieldsX Check deprecated fields of x and fix them based on a conversion table
%
% FORMAT:   x = xASL_io_CheckDeprecatedFieldsX(x[, bVerbose])
%
% INPUT:
%   x               - x structure (REQUIRED)
%   bVerbose        - Print extra information (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%   x               - x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Check deprecated fields of x and fix them based on a conversion table.
%               This table is used within:
%
%                - xASL_bids_parms2BIDS
%                - xASL_io_ReadDataPar
%                - xASL_adm_LoadParms
%                - xASL_adm_LoadX
%
%               It is not only used to convert deprecated x structure
%               fields to fields within up-to-date substructures of x, but
%               also to rename fields and to move them back and forwards
%               for the comparison with BIDS parameters within
%               xASL_bids_parms2BIDS e.g., which is why it is important to
%               make sure that if a row within the table is used to move &
%               rename, that there is also another row where the new
%               fieldname is moved to the same substructure.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:      n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Add default substructs
    [x] = xASL_init_SubStructs(x);

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
					renamedFields{nOutdatedParameter,1} = [conversionTable{iField,2} '.' conversionTable{iField,3}];
                elseif ~isfield(x.(conversionTable{iField,2}).(conversionTable{iField,3}),conversionTable{iField,4}) ...
                        && isstruct(x.(conversionTable{iField,2}).(conversionTable{iField,3}))
                    % Sub-Sub-Structure (x.FIELD -> x.STRUCT.STRUCT.FIELD)
                    x.(conversionTable{iField,2}).(conversionTable{iField,3}).(conversionTable{iField,4}) = x.(conversionTable{iField,1});
					renamedFields{nOutdatedParameter,1} = [conversionTable{iField,2} '.' conversionTable{iField,3} '.' conversionTable{iField,4}];
				else
					if bVerbose
                        fprintf('%s.%s.%s already exists...\n',conversionTable{iField,2},conversionTable{iField,3},conversionTable{iField,4});
					end
					renamedFields{nOutdatedParameter,1} = 'duplicit field';
                end
            end
            % Remove old field
            x = rmfield(x,conversionTable{iField,1});
        end
    end
    
    if bOldFieldsDetected && bVerbose
        % Print warning and individual fields
        fprintf('Detected and corrected deprecated fields in dataPar.json. Please modify the file accordingly:\n');
        for iField = 1:size(detectedFields,1)
            fprintf('%s --> %s\n', detectedFields{iField,1}, renamedFields{iField,1});
        end
    end

end
