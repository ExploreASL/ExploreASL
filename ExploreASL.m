function [x] = ExploreASL(varargin)
%ExploreASL Alias to ExploreASL_Master
%
% FORMAT: [x] = ExploreASL([DataParPath, ImportArray, ProcessArray, SkipPause, iWorker, nWorkers])
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function calls the masterscript ExploreASL_Master.
% Please type help ExploreASL_Master for more explanation.
% The function ExploreASL_Master will be renamed into ExploreASL with
% version 2.0.
% __________________________________
% Copyright 2015-2021 ExploreASL

    x = ExploreASL_Master(varargin{:});

    if x.ProcessData==2 && nargout==0
        warning('Data loading requested but no output structure defined');
        fprintf('%s\n', 'Next time, try adding "x = " to the command to load data into the x structure');
        fprintf('%s\n', 'Hint: you can now run "x = ans" to load the data into the x structure');
        fprintf('%s\n', 'Note that this only works if you did not run another command and generated a new ans variable in between');
    end

end
