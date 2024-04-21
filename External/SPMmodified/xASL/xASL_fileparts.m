function [Fpath, Ffile, Fext, SuffixSPM] = xASL_fileparts(InputPath)
%xASL_fileparts Wrapper around fileparts for special extensions
% FORMAT: [Fpath, Ffile, Fext] = xASL_fileparts(InputPath)
%
% INPUT:
%   InputPath - input file name with path, ending either .ext(.gz) or .ext(.gz),n
% OUTPUT:
%   Fpath     - path or empty string
%   Ffile     - file name
%   Fext      - file extension or empty string
%   SuffixSPM - volume extension for NIfTI files (per spm usage)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Returns the path, file name, and file extension for InputPath using the fileparts.m function.
%              If a file ending at nii.gz is given, then the whole nii.gz is returned as the extension.
%              Does not verify the existence of the file, or existence of .nii or .nii.gz.
%              This function uses the following steps:
%              1. Catch SPM suffixes
%              2. Manage .gz (double) extensions  
%              3. Manage folders differently than files
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES: ['/path' 'file' '.nii']         = xASL_fileparts('/path/file.nii');
%          ['/path' 'file' '.nii.gz']      = xASL_fileparts('/path/file.nii.gz');
%          ['/path' 'file' '.nii.gz' ',2'] = xASL_fileparts('/path/file.nii.gz,2');
%          ['c:\path' 'file' '.nii']       = xASL_fileparts('c:\path\file.nii') (on Windows)
% __________________________________
% Copyright 2015-2024 ExploreASL

    %% Admin
    % First catch cell
    if iscell(InputPath)
        if length(InputPath)>1
            warning('InputPath contained multiple cells, using the first only');
        end
        InputPath = InputPath{1};
    end

    %% 1. Catch SPM suffixes
    [Fpath, Ffile, Fext] = fileparts(InputPath);

    [Ind1, Ind2] = regexp(Fext, ',\d+');
    if ~isempty(Ind1)
        if Ind2~=length(Fext)
            error('Something wrong with this file extension');
        else
            SuffixSPM = Fext(Ind1:Ind2);
            InputPath = fullfile(Fpath, [Ffile Fext(1:Ind1-1)]);
        end
    else
        SuffixSPM = '';
    end

    %% 2. Manage .gz (double) extensions
    [Fpath, Ffile, Fext] = fileparts(InputPath);
    [~, Ffile2, Fext2] = fileparts(Ffile);
    
    if strcmpi(Fext, '.gz') && ~isempty(Fext2)
        Ffile = Ffile2;
        Fext = [Fext2 Fext];
    end
    
    %% 3. Manage folders differently than files
    if exist(InputPath, 'dir')
        % when running this for a folder, there is no extension. So any 
        % extension should be ignored and put back to the Ffile
        Ffile = [Ffile Fext];
        Fext = '';
    end

end