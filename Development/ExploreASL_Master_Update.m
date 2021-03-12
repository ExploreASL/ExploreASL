function [x] = ExploreASL_Master_Update(varargin)

    % Test: [x] = ExploreASL_Master_Update('MY_PATH','[1 1 1]','1','1','1','[1 2 3]');
    % Test: [x] = ExploreASL_Master_Update('MY_PATH',[1 1 1],1,1,1,[1 2 3]);

    % Define input parser
    p = inputParsing(varargin{:});
    
    % Convert parsed input
    parameters = convertParsedInput(p.Results);
    
    % Store parsed input
    x = storeParsedInput(parameters);
    
    % Print chosen settings
    printSettings(x);

end

%% Define input parser
function p = inputParsing(varargin)

    % Initialize input parser
    p = inputParser;
    
    % Define valid input variables
    validDataParPath = @(variable) ischar(variable) || isempty(variable);
    validProcessArray = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validSkipPause = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validiWorker = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validnWorkers = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validiModules = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    
    % Define defaults
    defaultDataParPath = [];
    defaultProcessArray = [0 0 0];
    defaultSkipPause = 1;
    defaultiWorker = 1;
    defaultnWorkers = 1;
    defaultiModules = [1 2 3];
    
    % Add definitions to the input parser
    addOptional(p, 'DataParPath', defaultDataParPath, validDataParPath);
    addOptional(p, 'ProcessArray', defaultProcessArray, validProcessArray);
    addOptional(p, 'SkipPause', defaultSkipPause, validSkipPause);
    addOptional(p, 'iWorker', defaultiWorker, validiWorker);
    addOptional(p, 'nWorkers', defaultnWorkers, validnWorkers);
    addOptional(p, 'iModules', defaultiModules, validiModules);
    
    % Parse input
    parse(p,varargin{:});

end

%% Convert parsed input
function parameters = convertParsedInput(parameters)

    % Check if inputs are chars
    if ischar(parameters.ProcessArray),     parameters.ProcessArray = str2num(parameters.ProcessArray);     end
    if ischar(parameters.SkipPause),        parameters.SkipPause = str2num(parameters.SkipPause);           end
    if ischar(parameters.iWorker),          parameters.iWorker = str2num(parameters.iWorker);               end
    if ischar(parameters.nWorkers),         parameters.nWorkers = str2num(parameters.nWorkers);             end
    if ischar(parameters.iModules),         parameters.iModules = str2num(parameters.iModules);             end

end

%% Store parsed input
function x = storeParsedInput(parameters)

    % Store input
    x.DataParPath = parameters.DataParPath;
    x.RunDCM2BIDS = parameters.ProcessArray(1);
    x.RunBIDS2RAW = parameters.ProcessArray(2);
    x.RunEXPLOREASL = parameters.ProcessArray(3);
    x.SkipPause = parameters.SkipPause;
    x.iWorker = parameters.iWorker;
    x.nWorkers = parameters.nWorkers;
    x.iModules = parameters.iModules;
    
end

%% Print chosen settings
function printSettings(x)

    fprintf('========== ExploreASL Settings ========== \n');
    fprintf('DataParPath\t\t%s\n', x.DataParPath);
    fprintf('RunDCM2BIDS\t\t%d\n', x.RunDCM2BIDS);
    fprintf('RunBIDS2RAW\t\t%d\n', x.RunBIDS2RAW);
    fprintf('RunEXPLOREASL\t%d\n', x.RunEXPLOREASL);
    fprintf('SkipPause\t\t%d\n', x.SkipPause);
    fprintf('iWorker\t\t\t%d\n', x.iWorker);
    fprintf('nWorkers\t\t%d\n', x.nWorkers);
    fprintf('iModules\t\t%d %d %d\n', x.iModules(1), x.iModules(2), x.iModules(3));
    fprintf('========================================= \n');

end




