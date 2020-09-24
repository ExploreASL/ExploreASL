function OutputArg = xASL_exist(PathIn, Type)
%xASL_exist ExploreASL wrapper of Matlab exist function, checking for
%.nii & .nii.gz at same time
%
% FORMAT: xASL_exist(PathIn[,Type])
%
% INPUT:
%   PathIn 	- path to be checked if it exists (REQUIRED)
%   Type    - type of checking (e.g. 'var', 'file', 'dir', as with Matlab's built in exist function (OPTIONAL, DEFAULT = [])
% OUTPUT:
%   OutputArg - true or false depending on if exists
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Check if the given path exists, wrapper around the Matlab
%              exist function, to allow checking for either .nii or .nii.gz
%              Otherwise, exist is used normally.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_exist('/path/file.nii');
%          xASL_exist('/path/file.nii.gz','file');
%          xASL_exist('image','var');
% __________________________________
% Copyright (C) 2015-2020 ExploreASL

% It can have one or two input parameters
if nargin>2
	error('Two input parameters is maximum');
end

if nargin<2 || isempty(Type)
	Type = [];
elseif ~strcmp(Type, 'file')
    OutputArg = exist(PathIn, Type);
    return; % bypass rest if not a file
end

% Function only accepts strings on input
if ~ischar(PathIn)
    OutputArg = false; %%%% >>>>> difficult to debug
%     warning('xASL_exist: Cannot check PathIn');
    return;
end

% Check for both NII and NII.GZ if a file is given
if strcmp(PathIn(end-3:end),'.nii')
	PathInGZ = [PathIn '.gz'];
	
	OutputArg = exist(PathIn, 'file') || exist(PathInGZ, 'file');
elseif  strcmp(PathIn(end-2:end),'.gz')
	PathInNII = PathIn(1:end-3);
	OutputArg = exist(PathIn, 'file') || exist(PathInNII, 'file');
elseif ~isempty(Type)
	OutputArg = exist(PathIn, Type);
else
    OutputArg = exist(PathIn);
end
        
end

