function [x] = xASL_qc_CollectSoftwareVersions(x)
%xASL_qc_CollectSoftwareVersions Collect currently used software versions
%
% FORMAT: [x] = xASL_qc_CollectSoftwareVersions(x)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT:
%   x       - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects software versions for Matlab, SPM, CAT, LST & ExploreASL
%              If FSL is installed, it will obtain its version as well.
%              These are stored in x.Output.Software.
%
% EXAMPLE: x = xASL_qc_CollectSoftwareVersions(x);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


    %% Admin
    if nargin<1
        error('This function needs the ExploreASL x-struct as input');
    end

    %% Get Matlab version
    Software.Matlab = version;
    [startIndex, endIndex] = regexp(Software.Matlab, '^\d*\.\d*\.');
	if isempty(startIndex) || isempty(endIndex) || endIndex(1) < 2
		warning('Something wrong with the Matlab version');
	else
		Software.Matlab = Software.Matlab(startIndex(1):endIndex(1)-1);
	end
    % Extract matlab version until just before the second dot

    
    %% Get SPM version
    [~, rev_spm] = spm('Ver');
    Software.SPM12 = rev_spm;

    
    %% Get CAT version
    [~, rev_cat] = cat_version;
    Software.CAT12 = rev_cat;

    
    %% Get LST version
    TxtPath = fullfile(x.D.SPMDIR,'toolbox','LST','lst-version.txt');
    FID = fopen(TxtPath);
    LSTversion = textscan(FID,'%s');
    Software.LST = LSTversion{1}{1};

    
    %% Get FSL version
    Software.FSL = 'Not used/found'; % default

    if ispc % Manage PC WSL
        wslString = 'wsl '; % windows subsystem for linux
    else
        wslString = '';
    end
    
    % Print FSL directory, if used
    if isfield(x,'FSLdir') && ~min(isnan(x.FSLdir)) && exist(x.FSLdir, 'dir')
        [Result1, Result2] = system([wslString ' cat ' xASL_adm_UnixPath(fullfile(x.FSLdir,'etc','fslversion'))]);
        if Result1==0 && ischar(Result2)
            Software.FSL = Result2;
        end
    end

    
    %% Get ExploreASL version
    % PM: make similar function as spm('Ver')/cat_version

    VersionPath = xASL_adm_GetFileList(x.opts.MyPath, '^VERSION.*$', 'FPList', [0 Inf]);
    if isempty(VersionPath)
        warning('Could not obtain ExploreASL version, version file missing');
    else
        [~, Fname, Fext] = fileparts(VersionPath{1});
        Software.ExploreASL = [Fname(9:end) Fext];
    end
    

    %% Get ExploreASL commit for optimal provenance
    [~, xASL_gitCommit] = xASL_system(['cd ' x.opts.MyPath '; git rev-parse HEAD']);
    Software.ExploreASL_git = strtrim(xASL_gitCommit);

    %% Add software field to x output
    try
        x.Output.SoftwareVersion(1) = Software;
    catch
        x.Output = rmfield(x.Output,'SoftwareVersion');
        x.Output.SoftwareVersion(1) = Software;
    end
    
    
end