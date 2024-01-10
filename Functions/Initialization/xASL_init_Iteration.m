function [bAborted, xOut] = xASL_init_Iteration(x, moduleName, dryRun, stopAfterErrors, SelectedSubjects)
%xASL_init_Iteration Iterates through all combinations of the parameter sets defined in x and calls the xASL module for each.
%
% FORMAT: [bAborted, xOut] = xASL_init_Iteration(x, moduleName[, dryRun, stopAfterErrors])
%
% INPUT:
%   x                - ExploreASL x-struct (STRUCT, REQUIRED)
%   moduleName       - Should be a string xASL_module_xxxxx (STRING, REQUIRED)
%   dryRun           - Dry run - does not execute the module (BOOLEAN, OPTIONAL, DEFAULT = 0)
%                      This argument can also be defined as a field of x.
%                      This settings overrides what is in the x-struct.
%   stopAfterErrors  - Number of allowed errors before job iteration is stopped. (OPTIONAL, DEFAULT = INF)
%                      This argument can also be defined as a field of x.
%                      This settings overrides what is in the x-struct.
%   SelectedSubjects - If this field exist, it will replace x.SUBJECTS in this script, to allow 
%                      running certain subjects only when debugging. (OPTIONAL, DEFAULT = x.SUBJECTS)
% OUTPUT:
%   bAborted         - Report if the run was aborted (BOOLEAN)
%   xOut             - x-struct on the output (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:       Parses the settings and runs the DatabaseLoop sub-function.
%
% EXAMPLE:           [~,x] = xASL_init_Iteration(x,'xASL_module_ASL');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% General settings
    dbSettings.x                     = x;
    dbSettings.x.settings.RERUN      = false;
    dbSettings.x.settings.MUTEXID    = moduleName;

    % Lock dir
    dbSettings.x.dir.LockDir = fullfile(x.dir.xASLDerivatives, 'lock', moduleName);
    
	% Set the dryRun field
	if nargin < 3 || isempty(dryRun)
		if ~isfield(x.settings,'dryRun') || isempty(x.settings.dryRun)
			dbSettings.settings.dryRun = 0;
		else
			dbSettings.settings.dryRun = x.settings.dryRun;
		end
	else
		dbSettings.settings.dryRun = dryRun;
	end
	
    % Set the stopAfterErrors field
    if nargin < 4 || isempty(stopAfterErrors)
        if ~isfield(x.settings,'stopAfterErrors') || isempty(x.settings.stopAfterErrors)
            dbSettings.settings.stopAfterErrors = Inf;
        else
            dbSettings.settings.stopAfterErrors = x.settings.stopAfterErrors;
        end
    else
        dbSettings.settings.stopAfterErrors = stopAfterErrors;
    end
    
    % Check if x.SUBJECTS and x.SESSIONS exist (they should at least be
    % initialized, could also be a dummy session)
    % Note that "SelectedSubjects" here is a variable only used when
    % debugging
    if nargin<5 || isempty(SelectedSubjects)
        if ~isfield(x,'SUBJECTS')
            error('x.SUBJECTS field missing: loading subjects went wrong, verify that ExploreASL was correctly initialized');
        elseif isempty(x.SUBJECTS)
            error('x.SUBJECTS was empty: loading subjects went wrong, verify that ExploreASL was correctly initialized');
        end

        try
            SelectedSubjects = x.SUBJECTS;
        catch ME
            error('Something wrong with x.SUBJECTS: loading subjects went wrong, verify that ExploreASL was correctly initialized: \n%s',ME.message);
        end
    end

    
    %% Print module name
    dbSettings.jobfn = str2func(moduleName);
    fprintf('\nRunning %s...\n', moduleName);
    
    %% Specific settings
    
    % Get module name string
    if  strcmp(moduleName(1:12),'xASL_module_')
        ModName = moduleName(13:end);
    elseif strcmp(moduleName(1:7),'module_')
        % backward compatibility
        ModName = moduleName(8:end);
    else
        ModName = moduleName;
    end
    
    % SUBJECT, SUBJECTDIR & LockDir
    if ~isempty(regexpi(ModName,'(Import|BIDS2Legacy)', 'once'))
        dbSettings.x.dir.SUBJECTDIR = fullfile('<ROOT>', '<SUBJECT>');
        dbSettings.x.dir.LockDir = fullfile(dbSettings.x.dir.LockDir, '<SUBJECT>');
    elseif ~isempty(regexpi(ModName,'Population', 'once'))
        % For population module, we do not specify SUBJECTDIR
        % For population module, we do not add subject name to the lock-folder
    else
        dbSettings.x.dir.SUBJECTDIR = fullfile(x.dir.xASLDerivatives, '<SUBJECT>');
        dbSettings.x.dir.LockDir = fullfile(dbSettings.x.dir.LockDir, '<SUBJECT>');
    end
    
    if isempty(regexpi(ModName,'Population', 'once'))
        % we don't want the population module to run for multiple sets (subjects or sessions)
        dbSettings.sets.SUBJECT = SelectedSubjects; % x.SUBJECTS
    end

    
    % Checking SESSIONS
    if isempty(regexpi(ModName,'(BIDS2Legacy|Population)', 'once'))
        % We don't use this for the BIDS2Legacy or population modules
        if ~isfield(x,'SESSIONS')
            error('x.SESSIONS field missing: loading sessions went wrong, verify that ExploreASL was correctly initialized');
        elseif isempty(x.SESSIONS)
            error('x.SESSIONS was empty: loading sessions went wrong, verify that ExploreASL was correctly initialized');
        end
    end    

    
    % SESSION, MUTEXID & SESSIONDIR
    if ~isempty(regexpi(ModName,'(ASL|func|dwi)', 'once'))
        dbSettings.sets.SESSION = x.SESSIONS;
        dbSettings.x.settings.MUTEXID = [dbSettings.x.settings.MUTEXID '_<SESSION>'];
        dbSettings.x.dir.SESSIONDIR = fullfile(x.dir.xASLDerivatives, '<SUBJECT>', '<SESSION>');
    elseif ~isempty(regexpi(ModName,'(Import)', 'once'))
        dbSettings.sets.SESSION = x.SESSIONS;
        dbSettings.x.settings.MUTEXID = [dbSettings.x.settings.MUTEXID]; % Currently we do not support session-wise import
        dbSettings.x.dir.SESSIONDIR = fullfile('<ROOT>', '<SUBJECT>', '<SESSION>');
    end
    
    %% Diary file
    if ~isempty(regexpi(ModName, '(DARTEL|Population|Analyze)', 'once'))
        dbSettings.diaryFile = fullfile(x.dir.xASLDerivatives, 'log', [moduleName '.log']);
    elseif ~isempty(regexpi(ModName,'Import', 'once'))
        dbSettings.diaryFile = fullfile(x.dir.xASLDerivatives ,'log', [moduleName '_sub-<SUBJECT>.log']);
    elseif ~isempty(regexpi(ModName,'(BIDS2Legacy|Structural|LongReg)', 'once'))
        dbSettings.diaryFile = fullfile(x.dir.xASLDerivatives, 'log', [moduleName '_<SUBJECT>.log']);
    elseif ~isempty(regexpi(ModName,'(ASL|func|dwi)', 'once'))        
        dbSettings.diaryFile = fullfile(x.dir.xASLDerivatives, 'log', [moduleName '_<SUBJECT>_<SESSION>.log']);
    else
        warning(['Unknown module name: ' moduleName]);
        % default
        dbSettings.diaryFile = fullfile(x.dir.xASLDerivatives, 'log', [moduleName '_<SUBJECT>_<SESSION>.log']);
    end
    
    %% Actually run the iteration
    [bAborted, xOut] = runIteration(dbSettings);
end


%% Run Iteration
function [bAborted, x] = runIteration(db)
% Iterates through all combinations of the parameter sets defined in db and calls a job function for each set.
%
% FORMAT: bAborted = runIteration(db)
%
% INPUT:
%   db        - a structure defining paramters sets and x-struct. 
%               This can also be a cell array containing multiple definitions.
%               format:
%                   db.sets.SETNAME = { ... }
%                   db.x.SYMBOLNAME = VALUE
%                   db.job          optional; Name of spm job template file or user defined for non-spm batches
%                   db.jobfn        optional; A user defined function that will be called on each dataset iteration
%                   db.diaryFile    optional; logfile
%                   db.settings.stopAfterErrors - optional: number of allowed errors before job iteration is stopped
% OUTPUT:
%   bAborted        - Report if the run was aborted
%   xOut            - x-struct on the output

    % Admin
    bAborted = false;

    % Run multi-level batch
    if numel(db)>1
        for iCell=1:numel(db)
            if iscell(db)
                dbi = db{iCell};
            else
                dbi = db(iCell);
            end
            if iCell==1
                % remember definitions from first run; might need them for following jobs
                db1 = dbi;
            else
                % clone x def's from first 'job' if none explicitly specified
                if isempty(dbi.x)
                    dbi.x = db1.x;
                end
                % clone set def's from first 'job' if none explicitly specified
                if isempty(dbi.sets)
                    dbi.sets = db1.sets;
                end
            end
            if runIteration(dbi)
                bAborted = true;
                break;
            end
        end
        
        return
    end

    % Make sure minimum required input is defined 
    if ~isfield(db,'sets') % add a dummy set to be able to run the script using parfor without sets
        db.sets.DUMMY = { '' };
    elseif numel(fieldnames(db.sets))<1
        db.sets.DUMMY = { '' };
    end
    
    % Parse arguments
    job = []; % not required for non-spm stuff
    jobfn = [];
	diaryFile = [];
    
	% Arguments can also be specified as structure field (for multi-level batches)
	if isfield(db, 'job') && ~isempty(db.job)
		job = db.job;
	end
	
	if isfield(db, 'jobfn') && ~isempty(db.jobfn)
		jobfn = db.jobfn;
	end
	
	if isfield(db, 'diaryFile') && ~isempty(db.diaryFile)
		diaryFile = db.diaryFile;
	end
	
	stopAfterErrors = db.settings.stopAfterErrors;
   
    % Which parameters sets should be iterated?
    setNames = fieldnames(db.sets);
    setCount = length(setNames);
    
    % Create an indexing and count vector
    I = zeros(1,setCount); % NextIter will fill with one's on first call
    N = zeros(1,setCount);
    for iSet=1:setCount
        setName = setNames{iSet};
        N(iSet) = length(db.sets.(setName));
    end
    
    %% Then loop through all permutations
    bAborted = false;
    iIter = 0; % just count the number of iterations
    CountErrors = 0;
    while ~bAborted 
        % do not use for iIter=1:prod(N) anymore because we now support dynamic sets with a fluctuating number of values
        % Iterate to the next permutation by increment one of the indices (possibly resetting others)
        % Remember the previous indices (Iprev) because dynamic sets have to be re-eavaluated when the 'parent' set changes
        Iprev=I;
        I = NextIter(I,N);
        if (sum(I)==0), break; end
        iIter = iIter+1;
        
        % Mark start time
        tocID = tic;
     
        % copy x-struct from user defined db
		x = db.x;
        
        % Backward compatibility
        Fields = {'ROOT' 'DARTELDIR'};

        if ~isfield(x,'D')
            x.D = struct;
        end
        
        for iL=1:length(Fields)
            if ~isfield(x.D,Fields{iL}) && isfield(x,Fields{iL})
                x.D.(Fields{iL})    = x.(Fields{iL});
            end
        end
        
        for iSet=1:setCount
            setName = setNames{iSet};
            set = db.sets.(setName);
            classid = class(db.sets.(setName));
            switch classid
                case 'cell'  % just a regular cell array with values to use one-by-one
                    if  isfield(x,setName)
                        x = rmfield(x,setName);
                    end
                    x.(setName) = set{I(iSet)}; % copy the indexed value from the set as regular (new) symbol
                    
                case 'function_handle'  % advanced usage: a user defined function will be called to get the set values
                    % A user defined function will be called when the 'parent' set iterates to the next item in its set.
                    % The dynamic set contents will be stored in db.x so it can be used at the next iteration
                    % without the need to call the user defined function again. So, it is important to realise that
                    % the ordering of the sets is crucial.
                    whichSetsChanged = find(I-Iprev); 
                    % only refresh the dynamic set when the parent set changes or when the dynamic set has no parent and is empty
                    if (iSet>1 && ~isempty(find(whichSetsChanged==iSet-1,1))) || (iSet==1 && isempty(db.x.(setName)))
                        S = set(x,db); % just call the specified function and pass the current symbol table as arg.
                        if ~iscell(S), S={S}; end % just in case: make sure the dynamic set is a cell array
                        db.x.(setName) = S;
                        I(iSet) = 1; % restart dynamic set
                        N(iSet) = length(db.x.(setName)); % set count
                        fprintf('\n[\b==============================================================================================]\b\n');
                        fprintf('new dynamic set %s with %d values: \n',setName,N(iSet)); disp(S);
                        if N(iSet)==0, N(iSet)=1; end % must allow one iteration when set is empty
                    end
                    % copy the indexed value from the dynamic set
                    if I(iSet)<=length(db.x.(setName))
                        x.(setName) = db.x.(setName){I(iSet)};
                    else
                        x.(setName) = []; % should not happen though! (Empty sets will be forced with N=1 above)
                    end
                otherwise
                    error('xASL_Iteration:illegalSymbolType', 'Dynamic parameter %s contains unsupported value type: %s!',setName,classid);
            end
        end
        
        % First check if this iteration has been fully processed, then we will skip the logging for this iteration
        if  exist(xASL_adm_ReplaceSymbols(fullfile(x.dir.LockDir,x.settings.MUTEXID,'999_ready.status'),x),'file')
            AlreadyProcessed = 1;
        else
            AlreadyProcessed = 0;
        end
        
        
        %% Optional diary logfile of command window
        if ~isempty(diaryFile)
            % NB. Don't overwrite diaryFile itself because it might be expanded
            % with different values in the next iteration.
            diaryFileEx = xASL_adm_ReplaceSymbols(diaryFile,x);
            diaryFileEx = fullfile(diaryFileEx);
                
            % Make sure folder exists, but limit to 2 missing levels
            status = xASL_adm_CreateDir(fileparts(diaryFileEx), 1);
            
            % Write diary file
            if status>=0
                diary(diaryFileEx);
            else
                error('xASL_Iteration:cannotCreateDirectory', 'Cannot create folder for diary file: %s!',diaryFileEx);
            end
        end
        

        %% Replace symbols
        if ischar(job)
            job_ex = xASL_adm_ReplaceSymbols(job,x);
                % NB. Don't overwrite job itself because it might be expanded
                % with different values in the next iteration.
            fprintf('Job file: %s\n',job_ex);
        else
            job_ex = job;
		end
		% current set parameters
		symbol_names = fieldnames(x);
		symbol_count = length(symbol_names);
		for iSymbol=1:symbol_count
			symbol_name = symbol_names{iSymbol};
			symbol_value_org = x.(symbol_name);
			if ischar(symbol_value_org)
				x.(symbol_name) = xASL_adm_ReplaceSymbols(symbol_value_org, x);
			end
        end
        % x.dir
        symbol_names_xdir = fieldnames(x.dir);
        symbol_count_xdir = length(symbol_names_xdir);
        for iSymbol=1:symbol_count_xdir
            symbol_name = symbol_names_xdir{iSymbol};
			symbol_value_org = x.dir.(symbol_name);
            if ischar(symbol_value_org)
				x.dir.(symbol_name) = xASL_adm_ReplaceSymbols(symbol_value_org, x);
            end
            % Make sure all paths have a correct format
            x.dir.(symbol_name) = fullfile(x.dir.(symbol_name));
        end
        % x.settings.MUTEXID
        symbol_names_xsettings = fieldnames(x.settings);
        symbol_count_xsettings = length(symbol_names_xsettings);
        for iSymbol=1:symbol_count_xsettings
            symbol_name = symbol_names_xsettings{iSymbol};
			symbol_value_org = x.settings.(symbol_name);
            if ischar(symbol_value_org)
				x.settings.(symbol_name) = xASL_adm_ReplaceSymbols(symbol_value_org, x);
            end
        end
        
        %% Some feedback about this iteration (after opening diary log)
        if ~AlreadyProcessed
            % Print the ExploreASL separator line
            fprintf('\n[\b==============================================================================================]\b\n');
            printModule = ['Module: ' x.settings.MUTEXID ', ' datestr(now)];
            
            [StartIndex, EndIndex] = regexp(diaryFileEx, '(?i)\/ASL_\d+\/'); %to find the name of the session inside diaryFileEx: ASL_with any digit after
             if ~isempty(StartIndex) %writes the session only for ASL module
                session = diaryFileEx(StartIndex(end)+1:EndIndex(end)-1); %isolate "ASL_1 or ASL_2 etc
                fprintf('%s\n',['Subject: ' x.SUBJECT ', Session: ' session ', ' printModule]);
             elseif ~isempty(regexp(diaryFileEx, 'xASL_module_Population.log'))
                 fprintf('%s\n', printModule); % Just print without a subject or session notation
             elseif isfield(x,'SUBJECT') 
                     fprintf('%s\n',['Subject: ' x.SUBJECT ', ' printModule]);
             else
                  %some datasets don't have any field that specifies the current subject
             end
            fprintf('\n');
        end        


        %% Start the job with all expanded x
        pwd_org = pwd; % remember the current working dir, so it can be restored if job fails
        try
            if  db.settings.dryRun
                fprintf('\nDRYRUN: skipping %s\n',func2str(jobfn));
            else 
				[result, x] = jobfn(x); % [result, x] = jobfn(x, job_ex);
                if ~result
                    error('xASL_Iteration:jobfnFailed', 'job function returned failed');
                end
            end
        catch ME
            
            % Obtain current subject ID from subject dir
            if isfield(x,'SUBJECT') && isfield(x, 'SUBJECTS')
                iSubject = find(strcmp(x.SUBJECTS, x.SUBJECT));
                subjectString = ['for subject ' xASL_num2str(iSubject) ': ' x.SUBJECT];
            else
                subjectString = '';
            end
                
            % Obtain current module
            ModuleString = func2str(jobfn);
            String2Remove = 'xASL_module';
            nString2Remove = length(String2Remove);
            IndexString2Remove = regexpi(ModuleString, String2Remove);
            
            if ~isempty(IndexString2Remove) && length(ModuleString)>nString2Remove
                ModuleString = ModuleString(nString2Remove+2:end);
            end
            
            fprintf('\n%s\n',['ERROR: ' ModuleString ' module terminated ' subjectString]);

            % Catch error in logging field
            [x] = xASL_qc_AddLoggingInfo(x, ME);
            
            ME.getReport()
            CountErrors     = CountErrors+1;
            if stopAfterErrors>0
                stopAfterErrors = stopAfterErrors -1;
                fprintf('\nCONT: but continue with next iteration!   \n');
            else
                % assume that there is a serious error in the template or
                % batch if it occurs at the very first run
                bAborted        = true; % bailing out
                % this should stop the inner DB loop in this script
                
            end
        end
        cd(pwd_org); % always go back to where we came from
        
        %% Print CPU time
        if ~AlreadyProcessed
            fprintf('\nJob-iteration %i stopped at %s and took %u seconds\n',iIter,datestr(now),ceil(toc(tocID)));
            fprintf('[\b==============================================================================================]\b\n');
        end
        
        %% Close logfile
        if ~isempty(diaryFile)
            diary off
        end
        
        if AlreadyProcessed
            fprintf('   ');
            xASL_TrackProgress(I,N);
            fprintf('\n');
        end
    end % next DB iteration
    fprintf('%s completed 100%s\n',x.settings.MUTEXID,'%');
   
end    

% Helper function to increment indices of parameter sets
function I = NextIter(I,N)
% This function increments the indices in vector I and by
% changing the 'least' significant indices first (=fastest running).
% Vector N should hold the maximum index values for all levels in I.
% The easiest way to start is to initialize N with the required values
% and then use I=ones(size(N)) as initial indexing vector.

% First check if all indices are zero.
if sum(I)==0
	I = I + 1; % set all to 1
	return
end

% First flip ordering of index vectors because we consider the set that was defined first (lowest index)
% as most significant (=slowest running index). The resulting vector I will be flipped back at the end.
I=fliplr(I);
N=fliplr(N);

% Which indices are about to overflow?
C = I==N;
if size(C,2)==1
	C = C';
end

% Convert this 0-1 pattern to a decimal number and add one
d = bin2dec(char(char('0')+fliplr(C))); % bin2dec accepts strings, but remember that the rightmost (least significant) char get highest index
d = d + 1; % increment, might overflow!
D = fliplr(dec2bin(d)-char('0')); % ascii pattern back to 0-1 integers
if length(D)>length(I)
	% all indices overflowed
	I = zeros(size(I));
else
	% check length of D; leading zeros will be missing
	if length(D)<length(I)
		D(length(I)) = 0; % pad zeros
	end
	I(D<C)=1; % set 'overflowed' indices back to 1
	ii=find(D>C); % this one must increment
	I(ii) = I(ii)+1;
end

% Reverse the flipping that was applied above
I=fliplr(I);
end % function NextIter
