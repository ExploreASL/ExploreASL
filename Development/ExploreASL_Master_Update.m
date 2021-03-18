function [x] = ExploreASL_Master_Update(varargin)

    % Test: [x] = ExploreASL_Master_Update('MY_PATH','[1 1 1]','[1 1 1]','1','1','1');
    % Test: [x] = ExploreASL_Master_Update('MY_PATH',[1 1 1],[1 1 1],1,1,1);
    % Test: [x] = ExploreASL_Master_Update('MY_PATH','[1 1 1]','[1 1 1]','1','1','1');
    % Test: [x] = ExploreASL_Master_Update('MY_PATH','[1 1 1]','[1 1 1]','1','1','1');

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
    validImportArray = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validProcessArray = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validSkipPause = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validiWorker = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    validnWorkers = @(variable) ischar(variable) || isempty(variable) || isnumeric(variable);
    
    % Define defaults
    defaultDataParPath = [];
    defaultImportArray = [0 0 0];
    defaultProcessArray = [0 0 0];
    defaultSkipPause = 1;
    defaultiWorker = 1;
    defaultnWorkers = 1;
    
    % Add definitions to the input parser
    addOptional(p, 'DataParPath', defaultDataParPath, validDataParPath);
    addOptional(p, 'ImportArray', defaultImportArray, validImportArray);
    addOptional(p, 'ProcessArray', defaultProcessArray, validProcessArray);
    addOptional(p, 'SkipPause', defaultSkipPause, validSkipPause);
    addOptional(p, 'iWorker', defaultiWorker, validiWorker);
    addOptional(p, 'nWorkers', defaultnWorkers, validnWorkers);
    
    % Parse input
    parse(p,varargin{:});

end

%% Convert parsed input
function parameters = convertParsedInput(parameters)

    % Check if inputs are empty or chars
    if isempty(parameters.DataParPath),     parameters.DataParPath = '';                                    end
    if ischar(parameters.ImportArray),      parameters.ImportArray = str2num(parameters.ImportArray);       end
    if ischar(parameters.ProcessArray),     parameters.ProcessArray = str2num(parameters.ProcessArray);     end
    if ischar(parameters.SkipPause),        parameters.SkipPause = str2num(parameters.SkipPause);           end
    if ischar(parameters.iWorker),          parameters.iWorker = str2num(parameters.iWorker);               end
    if ischar(parameters.nWorkers),         parameters.nWorkers = str2num(parameters.nWorkers);             end
    
    % Check length of arrays (single digit input)
    if length(parameters.ImportArray)<3
        parameters.ImportArray = [parameters.ImportArray(1),...
                                  parameters.ImportArray(1),...
                                  parameters.ImportArray(1)];
    end
    if length(parameters.ProcessArray)<3
        parameters.ProcessArray = [parameters.ProcessArray(1),...
                                   parameters.ProcessArray(1),...
                                   parameters.ProcessArray(1)];
    end
    
    % Different default for deployed mode
    if isdeployed
        parameters.SkipPause = 0;
    end


end

%% Store parsed input
function x = storeParsedInput(parameters)

    % Store input
    x.DataParPath = parameters.DataParPath;
    x.ImportArray = parameters.ImportArray;
    x.ProcessArray = parameters.ProcessArray;
    x.SkipPause = parameters.SkipPause;
    x.iWorker = parameters.iWorker;
    x.nWorkers = parameters.nWorkers;
    
end

%% Print chosen settings
function printSettings(x)

    fprintf('======================== ExploreASL Settings =========================\n');
    fprintf('DataParPath\t\t\t\t%s\n', x.DataParPath);
    fprintf('Run DCM2NII\t\t\t\t%d\n', x.ImportArray(1));
    fprintf('Run NII2BIDS\t\t\t%d\n', x.ImportArray(2));
    fprintf('Run BIDS2LEGACY\t\t\t%d\n', x.ImportArray(3));
    fprintf('Run Structural Module\t%d\n', x.ProcessArray(1));
    fprintf('Run ASL Module\t\t\t%d\n', x.ProcessArray(2));
    fprintf('Run Population Module\t%d\n', x.ProcessArray(3));
    fprintf('====================================================================== \n');

end




