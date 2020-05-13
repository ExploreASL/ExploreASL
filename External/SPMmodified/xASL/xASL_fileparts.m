function [Fpath, Ffile, Fext] = xASL_fileparts(InputPath)
% Wrapper around the fileparts.m that treats nii.gz as an extension
% (this also works with other extensions)
% FORMAT: [Fpath, Ffile, Fext] = xASL_fileparts(InputPath)
%
% INPUT:
%   InputPath - input file name with path, ending either .ext or .ext.gz
% OUTPUT:
%   Fpath     - path or empty string
%   Ffile     - file name
%   Fext      - file extension or empty string
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Returns the path, file name, and file extension for InputPath using the fileparts.m function.
%              If a file ending at nii.gz is given, then the whole nii.gz is returned as the extension.
%              Does not verify the existence of the file, or existence of .nii or .nii.gz
% EXAMPLE: xASL_fileparts('/path/file.nii');
%          xASL_fileparts('/path/file.nii.gz','file');
%          xASL_fileparts('c:\path\file.nii');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (R) 2015-2019 ExploreASL
%
% 2015-01-01 HJ

    % First catch cell
    if iscell(InputPath)
        if length(InputPath)>1
            warning('InputPath contained multiple cells, using the first only');
        end
        InputPath = InputPath{1};
    end
        
    % Now catch SPM suffixes
    [Ind1, Ind2] = regexp(InputPath,',\d');
    if ~isempty(Ind1)
        SuffixSPM = InputPath(Ind1:Ind2);
        InputPath = InputPath(1:Ind1-1);
    else
        SuffixSPM = '';
    end

    % Now do our core business
    [Fpath, Ffile, Fext] = fileparts(InputPath);
    [~, Ffile2, Fext2]   = fileparts(Ffile);
    
    if  strcmpi(Fext,'.gz') && ~isempty(Fext2)
        Ffile           = Ffile2;
        Fext            = [Fext2, Fext];
    end
    
    % Put the SuffixSPM back
    Fext = [Fext SuffixSPM];

    if exist(InputPath, 'dir')
        % when running this for a folder, there is no extension, and any
        % dot in the name should be ignored and put in the Ffile
        Ffile = [Ffile Fext];
        Fext = '';
    end

end

