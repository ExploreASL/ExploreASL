% __________________________________
% SPDX-License-Identifier: Apache-2.0
% __________________________________

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
end


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
    end
end


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


%% Add to /Functions/QualityControl/qc_glossary.tsv
x = ExploreASL;
pathTSV = fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'qc_glossary.tsv');
QCglossary = xASL_tsvRead(pathTSV);

pathMatMin = fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'xMinimalReferenceValues.mat');
pathMatMax = fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'xMaximalReferenceValues.mat');

MatMin = load(pathMatMin, '-mat');
MatMax = load(pathMatMax, '-mat');

structsAre = {'Structural' 'ASL'};

for iStruct=1:length(structsAre)

    fieldsMin = fields(MatMin.x.Output.(structsAre{iStruct}) );
    fieldsMax = fields(MatMax.x.Output.(structsAre{iStruct}) );
    
    if ~isequal(fieldsMin, fieldsMax)
        warning('fieldsMin & fieldsMax are unequal');
    end

    for iField=1:length(fieldsMin)
        keyIs = fieldsMin{iField};
        valueMin = MatMin.x.Output.(structsAre{iStruct}).(fieldsMin{iField});
        valueMax = MatMax.x.Output.(structsAre{iStruct}).(fieldsMin{iField});
    
        % Do some specific fixes for fields where min should be 0
        zeroFields = {'T1w_LR_flip_YesNo', 'FLAIR_WMH_vol_mL' 'FLAIR_WMH_n' 'MotionExcl_Perc' 'MotionMax_mm' 'MotionMean_mm' 'MotionSD_mm'};
        % explanation: no left-right flip, no WMH lesions, no WMH lesions zero motion exclusion is best
        if sum(strcmp(zeroFields, keyIs))>0
            valueMin = 0;
        end

        % Do some specific fixes for fields where max should be 1
        onesFields = {'TC_ASL2T1w_Perc' 'Mean_SSIM_Perc' 'T1w_IQR_Perc' 'TC_CBF2template' 'TC_M02template' 'tSNR_Slope_Corr' 'ASL_tSNR_Slope_Corr'};
        % explanation: no left-right flip, maximum overlap, perfect similarity, perfect IQR score
        if sum(strcmp(onesFields, keyIs))>0
            if valueMax<1
                valueMax = 1;
            else
                valueMax = 100;
            end
        end        

        if strcmp(keyIs, 'T1w_LR_flip_YesNo')
            valueMax = 0; % we don't want a left-right flip
        end

        % Parameters where higher is better, can have a higher maximum value
        maxList = {'_vol_' '_tSNR_' '_SNR_' 'CBF_WM_' 'CNR_' 'FBER'};

        bContainsMaxList = sum(cellfun(@(y) contains(keyIs, y), maxList))>0;
        if bContainsMaxList
            valueMax = 2 .* valueMax;
        end

        % Parameters where lower is better, can have a lower minimum value
        minList = {'AI_' '_SD_' 'tSD_' 'motion' 'SpatialCoV' 'RMSE' 'RigidBody' 'CBF_GM_WM_Ratio'};
        bContainsMinList = sum(cellfun(@(y) contains(keyIs, y), minList))>0;
        if bContainsMinList
            valueMin = 0.5 .* valueMin;
        end

        % Get the key-index from QCglossary
        indexGlossary = find(strcmp(QCglossary(:,1), keyIs));
        if numel(indexGlossary)<1
            warning(['Key missing: ' keyIs]);
        elseif numel(indexGlossary)>1
            warning(['Multiple keys present with the same name in QCglossary: ' keyIs]);
        else
    
            % Convert the numerical values to strings with the correct amount of floating points/depending on the order of magnitude
            valueMin = xASL_adm_formatWithRoundedMagnitude(valueMin, 3);
            valueMax = xASL_adm_formatWithRoundedMagnitude(valueMax, 3);
            
            rangeString = [valueMin '-' valueMax];
            QCglossary{indexGlossary, 5} = rangeString;
        end
    end
end      

% Store in new QC glossary
xASL_tsvWrite(QCglossary, pathTSV, 1);