function [x] = xASL_init_LoadDataParameterFile(x, DataParPath, SelectParFile)
%xASL_init_LoadDataParameterFile Load data parameter file
%
% FORMAT: [x] = xASL_init_LoadDataParameterFile(x, DataParPath, SelectParFile)
%
% INPUT:
%   x             - ExploreASL x structure (STRUCT, REQUIRED)
%   DataParPath   - Path to the data parameter file (CHAR ARRAY, PATH, REQUIRED)
%   SelectParFile - Variable which tells the import workflow if we have to ask the user for the study root directory a second time (BOOLEAN, REQUIRED)
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
%
% Copyright 2015-2021 ExploreASL


    if SelectParFile
        fprintf('ExploreASL requires a data parameter file, which can be either a .m or .json file\n');
        DataParPath = input('Please provide the full path of the DataParameter file, including ": ');
    end

    [pathstr, ~, Dext] = fileparts(DataParPath);
    xBackup = x;

    if strcmp(Dext,'.mat') % Now reading a mat-file, soon to be .json table complying with BIDS
        %% Mat file can only contain single variable, that should contain the x/parameters
        TempVar = load(DataParPath);
        FieldN = fields(TempVar);
        x = TempVar.(FieldN{1});
    elseif strcmp(Dext,'.json')
        % JSON import
        try
            x = xASL_io_ReadDataPar(DataParPath);
        catch ME1
            % if this fails, try to recreate the json file from an .m file,
            % if it exists
            [Fpath, Ffile] = fileparts(DataParPath);
            DataParPath = fullfile(Fpath, [Ffile '.m']);
            if exist(DataParPath, 'file')
                warning('Invalid DataPar JSON file, trying to repair from detected .m file');
                fprintf('A common issue is needing escaping e.g. "\\d" instead of "\d"\n');
                try
                    x = xASL_io_ReadDataPar(DataParPath);
                catch ME2
                    fprintf('%s\n', ME1.message);
                    fprintf('%s\n', ME2.message);
                    fprintf('A common issue is needing escaping e.g. "\\d" instead of "\d"\n');
                    error('Something went wrong loading the DataParFile');
                end
            else
                warning('Invalid DataPar JSON file');
                fprintf('A common issue is needing escaping e.g. "\\d" instead of "\d"\n');
                fprintf('%s\n', ME1.message);
                error('Something went wrong loading the DataParFile');
            end
        end
    elseif strcmp(Dext,'.m')
        try
            %% Backward compatibility
            x = xASL_io_ReadDataPar(DataParPath); % convert .m to .json and read it
        catch ME1
            try
                % Bypass eval error stuff with long names, spaces etc
                %% BACKWARD COMPATIBILITY FOR NOW, REPLACE WITH THE xASL_init_ConvertM2JSON ABOVE
                TempDataParPath = fullfile(pathstr,'TempDataPar.m');
                copyfile(DataParPath, TempDataParPath,'f' ); % overwrite ,if exist

                pathstr = fileparts(TempDataParPath);
                BackupCD = pwd;
                cd(pathstr);
                x = TempDataPar;
                cd(BackupCD);
                if exist(TempDataParPath,'file')
                    delete(TempDataParPath);
                end
            catch ME2
                fprintf('%s',ME1.message);
                fprintf('%s',ME2.message);
                error('Something went wrong loading the DataParFile');
            end
        end
    end
    % Put x fields back from backup
    FieldsAre = fields(xBackup);
    for iField=1:length(FieldsAre)
        if ~isfield(x,(FieldsAre{iField}))
            x.(FieldsAre{iField}) = xBackup.(FieldsAre{iField});
        end
    end
    
    if ~isfield(x,'D')
        x.D = struct;
    end
    
    if ~isfield(x.D,'ROOT')
        if isfield(x,'ROOT')
            x.D.ROOT = x.ROOT;
        else
            x.D.ROOT = pathstr; % default
        end
    end
    
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

