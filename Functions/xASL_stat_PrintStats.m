function [x] = xASL_stat_PrintStats(x, bFollowSubjectSessions)
%xASL_stat_PrintStats Print overview of (ROI-) data from all subjects/sessions in
%TSV file
%
% FORMAT: [x] = xASL_stat_PrintStats(x)
% 
% INPUT:
%   x                   - struct containing statistical pipeline environment parameters (REQUIRED)
%   x.S.SaveFile        - name of the TSV output file
%   x.S.SetsName        - names of the covariates/sets
%   x.S.NamesROI        - names of the ROIs
%   x.S.SubjectSessionID       - names of the Subject/Sessions in the order of x.S.SetsID/x.S.DAT (as specified by bFollowSubjectSessions)
%   x.S.SetsID          - values of the sets
%   x.S.SetsOptions     - options for sets (if the sets are ordinal)
%   bFollowSubjectSessions - boolean specifying which source determines the
%                            subjects/sessions that we process:
%                            TRUE: follow x.SUBJECTS & x.SESSIONS
%                            FALSE: follow data in x.S.DAT
%                            (OPTIONAL, DEFAULT: FALSE)
%                         
% OUTPUT:
%   x                   - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function prints an overview of statistics from
%              data that were acquired per ROI, in a TSV file. It starts by
%              printing covariates (called "Sets"). Rows will be
%              subjects/sessions, columns will be the sets and
%              ROI-statistics.
%              Any missing data will be skipped (setting them to NaN should
%              have happened in a previous function).
%
%              This function performs the following steps:
%
%              1. First remove previous TSV-file, if already existed
%                 printing to a TSV file can be tricky if it is opened by
%                 Excel. Make sure to close previous versions first,
%                 otherwise this part will crash.
%              2. Print overview of sets to TSV
%                 as explained above. Uses subfunction
%                 xASL_stat_CreateLegend to put legends. Aim is to create a
%                 single TSV file that has a proper overview of the data,
%                 & is self-explanatory to those reading/using it.
%              3. Define number of ASL sessions, force to 1 in case of TT or volume metrics
%              4. Print the overview
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_stat_PrintStats(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% Admin
if nargin<2 || isempty(bFollowSubjectSessions)
    bFollowSubjectSessions = false;
end


%% -----------------------------------------------------------------------------------------------
%% 1) First remove previous TSV-file, if already existed
fclose all; % safety procedure
if isfield(x.S,'summary_fid')
    x.S = rmfield(x.S,'summary_fid'); 
end
    
try
    xASL_delete(x.S.SaveFile);
catch ME
    warning(['Couldnt delete ' x.S.SaveFile ', if it was opened, please close this file first']);
    fprintf('%s\n',['Message: ' ME.message]);
    return;
end




%% -----------------------------------------------------------------------------------------------
%% 2) Print overview of sets to TSV
try
    xASL_adm_CreateDir(fileparts(x.S.SaveFile));
catch ME
    warning(['Couldnt delete ' x.S.SaveFile ', if it was opened, please close this file first']);
    fprintf('%s\n',['Message: ' ME.message]);
    return;    
end

% Build cell array 'SUBJECT' & x.S.SetsName{iSet} & x.S.NamesROI{ii}
[~, thisFileName] = fileparts(x.S.SaveFile);
thisFile = matlab.lang.makeValidName(thisFileName);
[x, statCell] = xASL_stat_PrintStats_GetStatCellArray(x);
x.modules.population.(thisFile) = statCell;

%% -----------------------------------------------------------------------------------------------
%% 3) Define number of sessions, force to 1 in case of TT or volume metrics
SingleSessions = {'volume' 'TT' 'PV_pGM' 'PV_pWM' 'pGM' 'pWM' 'mrc1T1' 'mrc2T1'};

HasSingleSessionOnly = sum(cellfun(@(y) ~isempty(regexpi(y,['^' x.S.output_ID])), SingleSessions));

if HasSingleSessionOnly
    nSessions = 1;
else
    nSessions = xASL_adm_GetPopulationSessions(x); % obtain number of Sessions by determining amount of input files present in the Population folder
end


%% -----------------------------------------------------------------------------------------------
%% 4) Print overview

if bFollowSubjectSessions
    
    for iSubject=1:x.nSubjects
        for iSession=1:nSessions
            iSubjectSession = (iSubject-1)* nSessions +iSession;

            % check first if this SubjectSession has data, otherwise skip &
            % issue a warning
            if length(x.S.SubjectSessionID)<iSubjectSession || size(x.S.DAT, 1)<iSubjectSession || size(x.S.SetsID,1)<iSubjectSession
                warning(['Missing data, skipping printing data for ' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession}]);
            else
                % print subject name
                x.modules.population.(thisFile){iSubjectSession+2,1} = x.S.SubjectSessionID{iSubjectSession, 1}; % I'm not 100% sure this is correct, how do I test this?

                %% Print the covariates and data
                iSubjectSession_SetsID = iSubjectSession;
                iSubjectSession_DAT = iSubjectSession;
                bPrintSessions = true;
                x.modules.population.(thisFile) = xASL_stat_PrintStats_FillStatCellArray(x,x.modules.population.(thisFile),iSubjectSession,iSubjectSession_SetsID, iSubjectSession_DAT, bPrintSessions);

            end
        end
    end

else
    
    % Get Subject regular expression
    SubjectExpression = x.dataset.subjectRegexp;
    if strcmp(SubjectExpression(end), '$') % remove this for allowing an ASL suffix
        SubjectExpression = SubjectExpression(1:end-1);
    end

    printedSessionN = 0;
    
    % Initialize empty table
    numElements = size(x.modules.population.(thisFile),2);
    numSubjectsSessions = length(x.S.SubjectSessionID);
    for iSubjSess=1:numSubjectsSessions
        x.modules.population.(thisFile)(2+iSubjSess,:) = repmat({nan(1, 1)}, 1, numElements);
    end
    
    for iSubjSess=1:numSubjectsSessions
        % x.S.SubjectSessionID == subject/session IDs created in xASL_stat_GetROIStatistics

        % Get subject ID
        [startSubjectIndex, endSubjectIndex] = regexp(x.S.SubjectSessionID{iSubjSess}, SubjectExpression);
        
        if isempty(startSubjectIndex) || isempty(endSubjectIndex)
            warning(['Could not find subject for ' x.S.SubjectSessionID{iSubjSess}]);
		else
			% Look also for the session ID including the ASL_\d substring
			[startSessionIndex, endSessionIndex] = regexp(x.S.SubjectSessionID{iSubjSess}, 'ASL_\d+');
            
            if isempty(startSessionIndex) || isempty(endSessionIndex)
				% If the sessionID was not found, then report a warning and use the entire expression for subjectID
				SubjectID = x.S.SubjectSessionID{iSubjSess}(startSubjectIndex:endSubjectIndex);
                warning(['Could not find session for ' x.S.SubjectSessionID{iSubjSess}]);
			else
				% If session ID was identified, we have to double-check that session ID is not part of subject ID and exclude if necessary
				SubjectID = x.S.SubjectSessionID{iSubjSess}(startSubjectIndex:min(endSubjectIndex, startSessionIndex-2));
				
                % Fix visit suffix
                [StartIndexTemp, EndIndexTemp] = find(regexp(SubjectID, '_\d+'));
                if isempty(StartIndexTemp) && isempty(EndIndexTemp)
                    [StartIndexTemp, EndIndexTemp] = regexp(x.S.SubjectSessionID{iSubjSess}(endSubjectIndex+1:startSessionIndex-2), '_\d+');
                    if ~isempty(StartIndexTemp)
                        SubjectID = [SubjectID x.S.SubjectSessionID{iSubjSess}(endSubjectIndex+StartIndexTemp:endSubjectIndex+EndIndexTemp)];
                    end
                end
                
				% Get session ID
                SessionID = x.S.SubjectSessionID{iSubjSess}(startSessionIndex:endSessionIndex);
                
                % Write fields to dataset
                x.dataset.currentSubjectID = SubjectID;
                x.dataset.currentSessionID = SessionID;
                
                % Fill stat cell array
                x.modules.population.(thisFile) = xASL_stat_PrintStats_AddSubjectSessionStatCellArray(x,x.modules.population.(thisFile),iSubjSess);
                

                %% print values for other covariates
                if isfield(x.S,'SetsID')
                    % Ensure to match subject/session
                    SubjectIndex = find(strcmp(x.SUBJECTS, SubjectID));
                    SessionIndex = find(strcmp(x.SESSIONS, SessionID));

                    if isempty(SubjectIndex)
                        warning(['Could not find subject ' SubjectID ', skipping']);
                    else
                        SessionColumn = find(strcmpi(x.S.SetsName, 'session'));
                        if isempty(SessionColumn) || length(SessionColumn)>1
                            warning('Could not find session data');
                        else
                            SessionN = iSubjSess;
                            if isempty(SessionN) || ~isnumeric(SessionN)
                                warning(['Something wrong with session ' SessionID]);
                            elseif SessionIndex>x.dataset.nSessions
                                if ~max(printedSessionN==SessionN)
                                    warning('Could not find values for other covariates');
                                end
                                printedSessionN = [printedSessionN SessionN];
                            else
                                
                                %% Print the covariates and data
                                iSubjectSession_SetsID = x.dataset.nSessions*(SubjectIndex-1)+SessionN;
                                iSubjectSession_DAT = iSubjSess;
                                bPrintSessions = false;
                                
                                % Write it to the cell array instead
                                x.modules.population.(thisFile) = xASL_stat_PrintStats_FillStatCellArray(x,x.modules.population.(thisFile),iSubjSess,iSubjectSession_SetsID, iSubjectSession_DAT, bPrintSessions);
                            end
                        end
                    end
                end
            end
        end
    end
end

% Write table to cell array
fprintf('Saving %s...\n', thisFile);
xASL_tsvWrite(x.modules.population.(thisFile),x.S.SaveFile,1);
x.S = rmfield(x.S,'SaveFile');


end


%% Legend
function [Legend] = xASL_stat_CreateLegend(x)
%xASL_stat_CreateLegend Create a row with legends to be printed in the TSV file
%containing ASL analysis stats

    if ~isfield(x,'S')
        x.S = struct;
    end

    Legend{1} = 'StudyID';
    
    if isfield(x.S,'SetsName')
        for iSet=1:length(x.S.SetsName)
            switch x.S.SetsName{iSet}
                case {'LongitudinalTimePoint', 'SubjectNList', 'Site'}
                    Legend{iSet+1} = 'integer';
                case 'AcquisitionTime'
                    Legend{iSet+1} = 'hhmmss';
                case {'GM_vol', 'GM_L', 'WM_vol', 'WM_L', 'CSF_vol', 'CSF_L'}
                    Legend{iSet+1} = 'Liter'; % 1 Liter = 10 cm*10 cm * 10 cm
                case {'WMH_vol'}
                    Legend{iSet+1} = 'mL';
                case 'GM_ICVRatio'
                    Legend{iSet+1} = 'ratio GM/ICV';
                case 'GMWM_ICVRatio'    
                    Legend{iSet+1} = 'ratio (GM+WM)/ICV';
                case {'CBF_spatial_CoV', 'PseudoCBF_spatial_CoV'}
                    Legend{iSet+1} = '(SD/mean)';
                case 'CBF_spatial_CoV_norm'
                    Legend{iSet+1} = 'ratio (spatial CoV/pseudoCBF spatial CoV)';
                case 'Age'
                    Legend{iSet+1} = 'years';
                case 'WMH_count'
                    Legend{iSet+1} = 'n lesions (integer)';
                case 'WMH_vol'
                    Legend{iSet+1} = 'mL'; % 1 mL = 10 mm * 10 mm * 10 mm (so 1000 1 mm^3 voxels)
                case 'MeanMotion'
                    Legend{iSet+1} = 'mm';
                otherwise
                    Legend{iSet+1} = '...';
            end
        end
        
        for iSet=1:length(x.S.SetsName)
            if strcmp(Legend{iSet+1},'...') && ~isempty(strfind(x.S.SetsName{iSet},'count'))
                Legend{iSet+1} = 'integer';
            end
        end
    end
    
    if isfield(x.S,'NamesROI') && isfield(x.S,'SaveFile')
        for iM=1:length(x.S.NamesROI)
            if ~isempty(findstr(x.S.SaveFile,'spatialCoV'))
                Legend{end+1} = '% SD/mean';
            elseif ~isempty(findstr(x.S.SaveFile,'CBF'))
                Legend{end+1} = 'mL/100g/min';
            else
                Legend{end+1} = '...';
            end
        end
    else
        warning('Could not find ROI names or filename to save, tsv header may be invalid');
    end

end

%% Function to write statistics to cell array so that we can use xASL_tsvWrite later on
function [x, statCell] = xASL_stat_PrintStats_GetStatCellArray(x)

    x.S.Legend = xASL_stat_CreateLegend(x);

    % Columns
    iCell = 1;
    iCell2 = 1;
    
    % First row
    statCell{iCell,1} = 'SUBJECT';
    iCell = iCell+1;
    
    % Sets
    if isfield(x.S,'SetsName')
        for iSet=1:length(x.S.SetsName)
            statCell{1,iCell} = x.S.SetsName{iSet};
            iCell = iCell+1;
        end
    end

    % ROIs
    if isfield(x.S,'NamesROI')
        for ii=1:length(x.S.NamesROI)
            statCell{1,iCell} = x.S.NamesROI{ii};
            iCell = iCell+1;
        end
    end
    
    % Second row
    
    % Legend
    if isfield(x.S,'Legend')
        for ii=1:length(x.S.Legend)
            statCell{2,iCell2} = x.S.Legend{ii};
            iCell2 = iCell2+1;
        end
    end

end

%% Fill the stat cell array with all subjects and sessions xASL_tsvWrite later on
function statCell = xASL_stat_PrintStats_AddSubjectSessionStatCellArray(x,statCell,rowNum)

    % Skip the first two rows (Labels & Legend)
    rowNum = rowNum+2;
    
    % Add subject and session of current row
    statCell{rowNum,1} = x.dataset.currentSubjectID;
    statCell{rowNum,2} = x.dataset.currentSessionID;

end

%% Fill the stat cell array with all subjects and sessions xASL_tsvWrite later on
function statCell = xASL_stat_PrintStats_FillStatCellArray(x,statCell, rowNum, iSubjectSession_SetsID, iSubjectSession_DAT, bPrintSessions)

    % Skip the first two rows and columns (Labels & Legend, Subjects & Sessions)
    rowNum = rowNum+2;
    iCell = 3;
    
    % 1. print values for other covariates
    % here we always have nSubjects*nSessions values so use "SubjectSession"
    
    if bPrintSessions
        printMatrix = x.S.SetsID;
        optionsMatrix = x.S.SetsOptions;
        sampleMatrix = x.S.Sets1_2Sample;
    else
        SessionColumn = find(strcmpi(x.S.SetsName, 'session'));
        vector2Print = [1:size(x.S.SetsName,2)];
        vector2Print = vector2Print(vector2Print~=SessionColumn);
        
        printMatrix = x.S.SetsID(:, vector2Print);
        optionsMatrix = x.S.SetsOptions(:, vector2Print);
        sampleMatrix = x.S.Sets1_2Sample(:, vector2Print);
    end
    
    for iPrint=1:size(printMatrix,2)
        String2Print = printMatrix(iSubjectSession_DAT, iPrint); % Previously we used iSubjectSession_SetsID

        if length(optionsMatrix{iPrint})>1 % we need options
            if length(optionsMatrix{iPrint}) >= length(unique(printMatrix(:,iPrint)))-2 % allow for zeros & NaNs
                if isnumeric(String2Print) && (int16(String2Print) == String2Print) && String2Print>0 && x.S.Sets1_2Sample(iPrint)~=3
                    String2Print = optionsMatrix{iPrint}{String2Print};
                end
            end
        end
        statCell{rowNum,iCell} = xASL_num2str(String2Print);
        iCell = iCell+1;
    end

    % 2. Print data in x.S.DAT
    % This part is different for volume or TT, since there will be only 1 value per subject (this will be done by the above in which nSessions is set to 1
    for iPrint=1:size(x.S.DAT,2) % print actual data
        statCell{rowNum,iCell} = xASL_num2str(x.S.DAT(iSubjectSession_DAT, iPrint));
        iCell = iCell+1;
    end


end
