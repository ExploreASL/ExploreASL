function [x] = ExploreASL(varargin)
%ExploreASL Alias to ExploreASL_Master
%
% FORMAT: [x] = ExploreASL([StudyRoot, ImportModules, ProcessModules, bPause, iWorker, nWorkers])
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function calls the masterscript ExploreASL_Master.
% Please type help ExploreASL_Master for more explanation.
% The function ExploreASL_Master will be renamed into ExploreASL with
% version 2.0.0
% __________________________________
% Copyright 2015-2021 ExploreASL

    x = ExploreASL_Master(varargin{:});


end
