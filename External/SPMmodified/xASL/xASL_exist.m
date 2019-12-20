function OutputArg = xASL_exist(PathIn, Type)
% ExploreASL wrapper of Matlab exist function. For NIFTII file, checks also the GZ version
%
% FORMAT: xASL_exist(PathIn[,Type])
%
% INPUT:
%   PathIn 	- path to be checked if it exists (REQUIRED)
%   Type    - optional type of if .nii or .nii.gz should be checked. 
%             If input does not have a .nii or .nii.gz extension, this argument is ignored (OPTIONAL, DEFAULT = check both)
% OUTPUT:
%   OutputArg - true or false depending on if exists
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check if the given path exists, wrapper around the Matlab
%              exist function, to allow checking for either .nii or .nii.gz
%              Otherwise, exist is used normally.
% EXAMPLE: xASL_exist('/path/file.nii');
%          xASL_exist('/path/file.nii.gz','file');
%          xASL_exist('image','var');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

% It can have one or two input parameters
if nargin>2
	error('xASL_exist:Two input parameters is maximum.');
end

if nargin<2 
	Type = [];
end

% It only accepts strings on input
if ~ischar(PathIn)
    OutputArg = false; %%%% >>>>> difficult to debug
%     warning('xASL_exist: Cannot check PathIn');
    return;
end

% Checks for both NII and NII.GZ if a file is given
if strcmp(PathIn(end-3:end),'.nii')
	PathInGZ = [PathIn '.gz'];
	
	OutputArg = exist(PathIn,'file') || exist(PathInGZ,'file');
elseif  strcmp(PathIn(end-2:end),'.gz')
	PathInNII = PathIn(1:end-3);
	OutputArg = exist(PathIn,'file') || exist(PathInNII,'file');
elseif isempty(Type)
	% In case another input is given from NII and NII.GZ, just check for the given type
	OutputArg = exist(PathIn);
else
	OutputArg = exist(PathIn,Type);
end
        
end

