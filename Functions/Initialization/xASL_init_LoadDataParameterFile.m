function [x] = xASL_init_LoadDataParameterFile(x)
%xASL_init_LoadDataParameterFile Load data parameter file
%
% FORMAT: [x] = xASL_init_LoadDataParameterFile(x, DataParPath)
%
% INPUT:
%   x                 - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Load data parameter file.
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% Get file extension and back up the x structure
    [pathstr, ~, Dext] = fileparts(x.dir.dataPar);
    xBackup = x;

    % Reading the .json file
    if strcmp(Dext,'.json')
        x = xASL_io_ReadDataPar(x.dir.dataPar, false);
    elseif strcmp(Dext,'.m')
        error('No .m file backwards compatibility starting v2.0.0...');
    elseif strcmp(Dext,'.mat')
        error('No .mat file backwards compatibility starting v2.0.0...');
    end
    
    % Put x fields back from backup
	x = xASL_adm_MergeStructs(x,xBackup);
    
    % Directories
    if ~isfield(x,'D')
        x.D = struct;
    end
    
    % ROOT
    if ~isfield(x.D,'ROOT')
        if isfield(x,'ROOT')
            x.D.ROOT = x.ROOT;
        else
            % Default
            x.D.ROOT = pathstr;
        end
    end
    
    % Check if x.D.ROOT does not exist
    if ~exist(x.D.ROOT, 'dir')
        % Check if x.D.ROOT was defined as a relative path
        if exist(fullfile(pathstr,x.D.ROOT), 'dir')
            x.D.ROOT = fullfile(pathstr,x.D.ROOT);
            x.ROOT = x.D.ROOT;
        else
            warning([x.D.ROOT ' didnt exist as folder, trying path of DataPar file']);
            x.D.ROOT = pathstr;
            x.ROOT = pathstr;
        end
    end

end

