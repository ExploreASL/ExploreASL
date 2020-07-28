function [x] = ExploreASL(var1, var2, var3, var4, var5, var6)
%ExploreASL Alias to ExploreASL_Master
%
% FORMAT: [x] = ExploreASL([DataParPath, ProcessData, SkipPause, iWorker, nWorkers, iModules])
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function calls the masterscript ExploreASL_Master.
% Please type help ExploreASL_Master for more explanation.
% The function ExploreASL_Master will be renamed into ExploreASL with
% version 2.0.
% __________________________________
% Copyright 2015-2020 ExploreASL

if nargin<6
    var6 = [];
end
if nargin<5
    var5 = [];
end
if nargin<4
    var4 = [];
end
if nargin<3
    var3 = [];
end
if nargin<2
    var2 = [];
end
if nargin<1
    var1 = [];
end

x = ExploreASL_Master(var1, var2, var3, var4, var5, var6);


end
