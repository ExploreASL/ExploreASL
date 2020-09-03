function [x, IsLoaded] = xASL_adm_LoadX(x, Path_xASL, bOverwrite)
%xASL_adm_LoadX Load x.mat file that keeps track of QC output
% FORMAT: [x[, IsLoaded]] = xASL_adm_LoadX(x[, Path_xASL, bOverwrite])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   Path_xASL   - path to the x.mat that contains the QC output (OPTIONAL, DEFAULT = x.SUBJECTDIR/x.mat)
%   bOverwrite  - true to overwrite the current x structure with the
%                 x.Output & x.Output_im from the x.mat (OPTIONAL, DEFAULT=false)
%
% OUTPUT:
%   x           - as input
%   IsLoaded  - if x.mat was succesfully loaded & succesfully added/replaced x.Output/x.Output_im
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads x.Output & x.Output_im struct fields
%              from the x.mat on the hard drive & adds them to the current x struct
%              located in memory. If it didnt exist in the x.mat, it will
%              set IsLoaded to false, which can be catched externally & a warning issued if managed so
%              in the calling function. If it didnt exist in the memory x
%              struct, or bOverwrite was requested, the contents of x.mat
%              will be loaded to the memory x struct
%
% EXAMPLE: [x, IsLoaded] = xASL_adm_LoadX(x, fullfile(x.D.ROOT,'x.mat'), true);
%
% __________________________________
% Copyright (C) 2015-2020 ExploreASL



%% -------------------------------------
%  Admin
IsLoaded = false; % default

if nargin<2 || isempty(Path_xASL)
	if ~isfield(x,'SUBJECTDIR')
		error('Path_xASL not provided and x.SUBJECTDIR not defined.');
	end
	Path_xASL = fullfile(x.SUBJECTDIR,'x.mat');
end

if nargin<3 || isempty(bOverwrite)
    bOverwrite = false;
end
FieldNames = {'Output', 'Output_im'};

%% -------------------------------------
%  Load OldX
if ~exist(Path_xASL, 'file')
    fprintf('%s\n',['Couldnt load ' Path_xASL]);
    return;
else
    OldX = load(Path_xASL, '-mat');
    IsLoaded = true;
end

%% -------------------------------------
%  Try adding
for iField=1:length(FieldNames)
    if isfield(OldX.x,FieldNames{iField}) % if x.mat contained the field
        if ~isfield(x, FieldNames{iField}) || bOverwrite % if x not contained the field OR overwrite
            x.(FieldNames{iField}) = OldX.x.(FieldNames{iField});
        end
    else
        IsLoaded = false;
    end
end


end
