function [x] = xASL_adm_CleanUpX(x)
%xASL_adm_CleanUpX Clean-Up before processing pipeline
%
% FORMAT: x = xASL_adm_CleanUpX(x)
%
% INPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Clean-Up before processing pipeline.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_adm_CleanUpX(x);
%
% __________________________________
% Copyright 2015-2023 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________




    %% Clean-up of the x structure

    % We want to reload the derivatives data correctly, which is why we delete the following 
    % fields before we run the processing pipeline, we need them for xASL_init_Iteration though.
    if isfield(x,'D')
        x = rmfield(x,'D');
        x.D = struct;
    end
    if isfield(x,'SUBJECT')
        x = rmfield(x,'SUBJECT');
    end
    if isfield(x,'SUBJECTS')
        x = rmfield(x,'SUBJECTS');
    end
    if isfield(x,'SESSION')
        x = rmfield(x,'SESSION');
    end
    if isfield(x,'SESSIONS')
        x = rmfield(x,'SESSIONS');
    end
    


end


