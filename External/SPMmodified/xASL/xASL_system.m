function [result1, result2] = xASL_system(Command, bVerbose)
%xASL_system Properly run a system call from Matlab

% E.g., some user-specific initializations are in ~/.bashrc 
% or ~/.zshrc (depending on the CLI used, Linux by default uses
% bash, macOS by default uses zsh).

% INPUT     Command = system command to run
% bVerbose  boolean specifying if we want to print system output on the
%           screen (OPTIONAL, DEFAULT = true)

% OUTPUT    Result1 = 0 for correct run, nonzero for error
%           Result2 = screen output (useful for debugging an error)


if nargin<2 || isempty(bVerbose)
    bVerbose = true;
end

if bVerbose
    VerbosityString = '-echo';
else
    VerbosityString = '';
end


%% Add the startup loading
if ismac
    Command = ['source ~/.zshrc;' Command];
elseif isunix
    Command = ['source ~/.bashrc;' Command];
else % e.g., Windows
    % then we don't change the command line
end


%% Run the command
[result1, result2] = system(Command, VerbosityString);






end

