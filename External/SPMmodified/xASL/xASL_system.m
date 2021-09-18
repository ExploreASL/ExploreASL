function [result1, result2] = xASL_system(Command, bVerbose)
%xASL_system Properly run a system call from Matlab
%
% FORMAT: 
%   [result1, result2] = xASL_system(Command[, bVerbose])
%
% INPUT:
%   Command     - string containing command we want to call to the system (REQUIRED)
%   bVerbose    - boolean specifying if we want to print system output on the
%                 screen (OPTIONAL, DEFAULT = true)
%
% OUTPUT:
%   Result1     - 0 for correct run, nonzero for error
%   Result2     - screen output (useful for debugging an error)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: 
%
% This function allows running a system call from Matlab in an optimized fashion.
% E.g., it will use user-specific CLI initializations, which are in ~/.bashrc 
% or ~/.zshrc (depending on the CLI used, Linux by default uses
% bash, macOS by default uses zsh).
%
% It runs the following steps:
% 1. Initialize the user-specific startup lines
% 2. Run the command
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:
% [Result1, Result2] = xASL_system('git fetch;git pull');
% xASL_system('echo hello');
%
% __________________________________
% Copyright 2015-2021 ExploreASL



if nargin<2 || isempty(bVerbose)
    bVerbose = true;
end

if bVerbose
    VerbosityString = '-echo';
else
    VerbosityString = '';
end


%% 1. Initialize the user-specific startup lines
if ismac
    Command = ['source ~/.zshrc;' Command];
elseif isunix
    Command = ['source ~/.bashrc;' Command];
else % e.g., Windows
    % then we don't change the command line
end


%% 2. Run the command
[result1, result2] = system(Command, VerbosityString);



end