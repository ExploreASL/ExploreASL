function [x] = xASL_qc_AddLoggingInfo(x, loggingEntry)
%xASL_qc_AddLoggingInfo Logging of errors and warnings within the x structure
%
% FORMAT: [x] = xASL_qc_AddLoggingInfo(x, loggingEntry)
%
% INPUT:
%   x             - Struct containing pipeline environment parameters (STRUCT, REQUIRED)
%   loggingEntry  - Matlab exception (MException, REQUIRED)
%                   (should contain the message and the stack field)
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Logging of errors and warnings within the x structure.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    %% Try logging of errors and warnings
    try
        % Logging field (for warnings and errors)
        if ~isfield(x, 'logging')
            x.logging = struct;
        end
        % Get number of current errors & warnings
        numOfLoggingEntries = size(x.logging,2);
        % Determine current logging entry number
        if numOfLoggingEntries==1
            if ~isfield(x.logging, 'message')
                numOfLoggingEntries = 0;
            end
        end
        entryNumber = numOfLoggingEntries+1;
        % Add logging entry
        x.logging(entryNumber).message = loggingEntry.message;
        x.logging(entryNumber).stack = loggingEntry.stack;
        
    catch
        fprintf('xASL_qc_AddLoggingInfo failed...\n');
    end
    
end
