% Generate value ranges for QC keys for ASPIRE PDF report

rootDir = '/Users/hjmutsaerts/ExploreASL/ASL/EPAD';

x = ExploreASL(rootDir, 0, 0); % load the data

QCstructural = struct;
QCasl = struct;

for iSubject=1:x.dataset.nSubjects
    xASL_TrackProgress(iSubject, x.dataset.nSubjects);
    pathQC = fullfile(x.ROOT, x.SUBJECTS{iSubject}, 'x.mat');

    if ~exist(pathQC, 'file')
        fprintf('%s\n', ['Missing x.mat for subject: ' x.SUBJECTS{iSubject}]);
    else
        mat = load(pathQC);

        % Load structural fields
        fieldsAre = fieldnames(mat.x.Output.Structural);

        for iField=1:length(fieldsAre)
            fieldValue = mat.x.Output.Structural.(fieldsAre{iField});
            if isnumeric(fieldValue)
                if ~isfield(QCstructural, fieldsAre{iField}) % first create struct with empty matrix
                    QCstructural.(fieldsAre{iField}) = [];
                end
                QCstructural.(fieldsAre{iField})(iSubject, 1) = fieldValue; % add value for current subject
            end
        end
    
        % Load ASL fields, assuming single session only
        fieldsAre = fieldnames(mat.x.Output.ASL.ASL_1);

        for iField=1:length(fieldsAre)
            fieldValue = mat.x.Output.ASL.ASL_1.(fieldsAre{iField});
            if isnumeric(fieldValue)
                if ~isfield(QCasl, fieldsAre{iField}) % first create struct with empty matrix
                    QCasl.(fieldsAre{iField}) = [];
                end
                QCasl.(fieldsAre{iField})(iSubject, 1) = fieldValue; % add value for current subject
            end
        end
    end
end

json = struct;


%% Determine structural QC stats
% Initialization
QCstructuralMean = struct;
QCstructuralMedian = struct;
QCstructuralSD = struct;
QCstructuralMAD = struct;

QCstructuralMin = struct;
QCstructuralMax = struct;

fieldsAre = fieldnames(QCstructural);

for iField=1:length(fieldsAre)

    valuesAre = QCstructural.(fieldsAre{iField});

    nonfiniteValues = sum(~isfinite(valuesAre));
    zeroesValues = sum(valuesAre==0);
    if nonfiniteValues~=0
        fprintf('%s\n', [num2str(nonfiniteValues) ' non-finite values (NaN, Inf, etc) detected for ' fieldsAre{iField}]);
    elseif zeroesValues~=0
        fprintf('%s\n', [num2str(zeroesValues) ' zeroes detected for ' fieldsAre{iField}]);
    end

    QCstructuralMean.(fieldsAre{iField}) = xASL_stat_MeanNan(valuesAre(valuesAre~=0));
    QCstructuralMedian.(fieldsAre{iField}) = xASL_stat_MedianNan(valuesAre(valuesAre~=0));
    QCstructuralSD.(fieldsAre{iField}) = xASL_stat_StdNan(valuesAre(valuesAre~=0));
    QCstructuralMAD.(fieldsAre{iField}) = xASL_stat_MadNan(valuesAre(valuesAre~=0));

    % Determine robust min & max
    sortValues = sort(valuesAre);
    indexMin = round(0.05*length(sortValues));
    indexMax = round(0.95*length(sortValues));
    
    QCstructuralMin.(fieldsAre{iField}) = sortValues(indexMin);
    QCstructuralMax.(fieldsAre{iField}) = sortValues(indexMax);

    % QCstructuralMin.(fieldsAre{iField}) = QCstructuralMean.(fieldsAre{iField})-1.96.*QCstructuralSD.(fieldsAre{iField});
    % QCstructuralMax.(fieldsAre{iField}) = QCstructuralMean.(fieldsAre{iField})+1.96.*QCstructuralSD.(fieldsAre{iField});
end

% Do some specific fixes
QCstructuralMin.FLAIR_WMH_vol_mL = 0; % Healthy = no lesions
QCstructuralMin.FLAIR_WMH_n = 0;
QCstructuralMin.T1w_LR_flip_YesNo = 0;  % we want no flip
QCstructuralMax.T1w_LR_flip_YesNo = 0;

% Store in JSON
json.Structural.Mean = QCstructuralMean;
json.Structural.Median = QCstructuralMedian;
json.Structural.SD = QCstructuralSD;
json.Structural.MAD = QCstructuralMAD;
json.Structural.Min = QCstructuralMin;
json.Structural.Max = QCstructuralMax;

%% Determine ASL QC stats
% Initialization
QCaslMean = struct;
QCaslMedian = struct;
QCaslSD = struct;
QCaslMAD = struct;

QCaslMin = struct;
QCaslMax = struct;

fieldsAre = fieldnames(QCasl);
skipFields = {'BackgroundSuppressionNumberPulses', 'NumberOfAverages', 'uniqueEchoTime', 'nUniqueEchoTime', ...
    'uniqueInitial_PLD', 'nUniqueInitial_PLD', 'uniqueLabelingDuration', 'nUniqueLabelingDuration'};

for iField=1:length(fieldsAre)
    if sum(contains(skipFields, fieldsAre{iField}))==0 % skip skipFields

        valuesAre = QCasl.(fieldsAre{iField});
    
        nonfiniteValues = sum(~isfinite(valuesAre));
        zeroesValues = sum(valuesAre==0);
        if nonfiniteValues~=0
            fprintf('%s\n', [num2str(nonfiniteValues) ' non-finite values (NaN, Inf, etc) detected for ' fieldsAre{iField}]);
        elseif zeroesValues~=0
            fprintf('%s\n', [num2str(zeroesValues) ' zeroes detected for ' fieldsAre{iField}]);
        end
    
        keepZeroes = {'MotionExcl_Perc', 'LR_flip_YesNo'};

        if sum(contains(keepZeroes, fieldsAre{iField}))==0 % skip keepZeroes
            valuesAre = valuesAre(valuesAre~=0); % remove zeroes
        end

        QCaslMean.(fieldsAre{iField}) = xASL_stat_MeanNan(valuesAre);
        QCaslMedian.(fieldsAre{iField}) = xASL_stat_MedianNan(valuesAre);
        QCaslSD.(fieldsAre{iField}) = xASL_stat_StdNan(valuesAre);
        QCaslMAD.(fieldsAre{iField}) = xASL_stat_MadNan(valuesAre);
    
    % Determine robust min & max
    sortValues = sort(valuesAre);
    indexMin = round(0.05*length(sortValues));
    indexMax = round(0.95*length(sortValues));
    
    QCaslMin.(fieldsAre{iField}) = sortValues(indexMin);
    QCaslMax.(fieldsAre{iField}) = sortValues(indexMax);

        % QCaslMin.(fieldsAre{iField}) = QCaslMean.(fieldsAre{iField})-1.96.*QCaslSD.(fieldsAre{iField});
        % QCaslMax.(fieldsAre{iField}) = QCaslMean.(fieldsAre{iField})+1.96.*QCaslSD.(fieldsAre{iField});
    end
end

% Do some specific fixes
QCaslMin.RMSE_Perc = 0; % zero difference with template is also fine
QCaslMin.nRMSE_Perc = 0; 
QCaslMin.MotionExcl_Perc = 0; % zero exclusion is also fine
QCaslMin.MotionMean_mm = 0; % zero motion is best
QCaslMin.LR_flip_YesNo = 0; % we want no flip
QCaslMax.LR_flip_YesNo = 0;
QCaslMax.TC_CBF2template = 1; % Tanimoto coefficient ranges from 0-1, 1 is optimal
QCaslMax.TC_M02template = 1;

% Store in JSON
json.ASL.Mean = QCaslMean;
json.ASL.Median = QCaslMedian;
json.ASL.SD = QCaslSD;
json.ASL.MAD = QCaslMAD;
json.ASL.Min = QCaslMin;
json.ASL.Max = QCaslMax;

%% Store reference values in the same way as we store individual subject values
pathMatMin = fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'xMinimalReferenceValues.mat');
pathMatMax = fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'xMaximalReferenceValues.mat');

oldX = x;

clear x
x.Output.Structural = QCstructuralMin;
x.Output.ASL = QCaslMin;
save(pathMatMin, 'x');

clear x
x.Output.Structural = QCstructuralMax;
x.Output.ASL = QCaslMax;
save(pathMatMax, 'x');

x = oldX;


%% Also dump all in a JSON file
pathJson = fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'QC_collectionReference.json');
xASL_io_WriteJson(pathJson, json);