function x = xASL_init_Process(x)
%xASL_init_Process Initialization before ExploreASL_Process
%
% FORMAT: x = xASL_init_Process(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Initialization before ExploreASL_Process.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% Initialization

    % Initialize x struct
    x = xASL_init_SubStructs(x);

    
    %% Determine subject/session/run structure from rawdata
    x = xASL_imp_DetermineStructureFromRawdata(x);
    
    % Create logging directory if it does not exist already
    xASL_adm_CreateDir(fullfile(x.dir.xASLDerivatives,'log'));

    % Force x.D.ROOT to be derivatives
    x.D.ROOT = x.dir.xASLDerivatives;
    x.ROOT = x.dir.xASLDerivatives;

end



function [x] = xASL_imp_DetermineStructureFromRawdata(x)
    %xASL_imp_DetermineStructureFromRawdata Determine structure from rawdata
    %
    % FORMAT: [x] = xASL_imp_DetermineStructureFromRawdata(x)
    %
    % INPUT:
    %   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
    %
    % OUTPUT:
    %   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION:    Determine structure from rawdata.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % EXAMPLE:        n/a
    % __________________________________
    % Copyright 2015-2022 ExploreASL
    
    
        %% Load rawdata subjects -> x.SUBJECTS
        if xASL_exist(x.dir.RawData,'dir')
            
            % SUBJECTS
            x.SUBJECTS = xASL_adm_GetFileList(x.dir.RawData,[],false,[],true);
    

        %% Check if data can be loaded
        %  This is non-specific for BIDS2Legacy
            if isempty(x.SUBJECTS)
                warning('Unable to find subjects...');
                x.opts.bLoadData = false;
                x.opts.bLoadableData = false;
            else
                % We can probably load the data
                x.opts.bLoadableData = true;
            end
            
        else
            % Maybe a user does not have a BIDS sourcedata or rawdata directory
            % and only runs the ExploreASL workflow on derivatives data.
            % Maybe BIDS2Legacy is turned on, but it actually shouldn't be.
            fprintf(2,'There is no rawdata directory...\n');
            x.opts.bSkipBIDS2Legacy = true;
            % We need to try to load the data anyway right now, otherwise the
            % workflow will skip the processing for datasets without rawdata
            % completely.
            x.opts.bLoadableData = true;
        end
    



    end