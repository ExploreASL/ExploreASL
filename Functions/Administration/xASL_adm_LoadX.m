function [x, IsLoaded] = xASL_adm_LoadX(x, Path_xASL, bOverwrite)
%xASL_adm_LoadX Load x.mat file that keeps track of QC output
% FORMAT: [x[, IsLoaded]] = xASL_adm_LoadX(x[, Path_xASL, bOverwrite])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   Path_xASL   - path to the x.mat that contains the QC output (OPTIONAL, DEFAULT = x.dir.SUBJECTDIR/x.mat)
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
% 1. Admin
% 2. Load X-struct from disc
% 3. Look for and update deprecated fields
% 4. Add fields from disc to the current x-struct
%
% EXAMPLE: [x, IsLoaded] = xASL_adm_LoadX(x, fullfile(x.dir.xASLDerivatives,'x.mat'), true);
%
% __________________________________
% Copyright (C) 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

%% 1. Admin
IsLoaded = false; % default
% Check for deprecated field SUBJECTDIR
if isfield(x,'SUBJECTDIR')
	warning('Deprecated field. Please use x.dir.SUBJECTDIR instead of x.SUBJECTDIR');
	if ~isfield(x,'dir') || ~isfield(x.dir,'SUBJECTDIR')
		x.dir.SUBJECTDIR = x.SUBJECTDIR;
	end
	x = rmfield(x,'SUBJECTDIR');
end
if nargin<2 || isempty(Path_xASL)
	if ~isfield(x.dir,'SUBJECTDIR')
		error('Path_xASL not provided and x.dir.SUBJECTDIR not defined.');
	end
	Path_xASL = fullfile(x.dir.SUBJECTDIR,'x.mat');
end
if nargin<3 || isempty(bOverwrite)
    bOverwrite = false;
end
FieldNames = {'Output', 'Output_im'};
%% 2. Load X-struct from disc
if ~exist(Path_xASL, 'file')
    fprintf('%s\n',['Couldnt load ' Path_xASL]);
    return;
else
    OldX = load(Path_xASL, '-mat');
    IsLoaded = true;
end
%% 3. Look for and update deprecated fields
OldX.x = xASL_io_CheckDeprecatedFieldsX(OldX.x,true);
%% 4. Add fields from disc to the current x-struct
for iField=1:length(FieldNames)
    % Check if x.mat contained the field
    if isfield(OldX.x,FieldNames{iField})
        % Check if x does not contain the field OR overwrite
        if ~isfield(x, FieldNames{iField}) || bOverwrite
            x.(FieldNames{iField}) = OldX.x.(FieldNames{iField});
        end
    else
        IsLoaded = false;
    end
end
end
