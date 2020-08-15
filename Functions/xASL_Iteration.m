function [bAborted, xOut] = xASL_Iteration(x, moduleName, dryRun, stopAfterErrors, SelectedSubjects)
% Iterates through all combinations of the parameter sets defined in x and calls the xASL module for each.
%
% FORMAT: [bAborted, xOut] = xASL_Iteration(x, moduleName[, dryRun, stopAfterErrors])
%
% INPUT:
%   x               - ExploreASL x-struct
%   moduleName      - Should be a string xASL_module_xxxxx
%   dryRun          - Dry run - does not execute the module (DEFAULT 0). 
%                     This argument can also be defined as a field of x. This settings overrides what is in the x-struct.
%   stopAfterErrors - Number of allowed errors before job iteration is stopped (DEFAULT INF). 
%                     This argument can also be defined as a field of x. This settings overrides what is in the x-struct.
%   SelectedSubjects - if this field exist, it will replace x.SUBJECTS in
%                      this script, to allow running certain subjects only
%                      when debugging (OPTIONAL, DEFAULT=x.SUBJECTS)
% OUTPUT:
%   bAborted        - Report if the run was aborted
%   xOut            - x-struct on the output
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Parses the settings and runs the DatabaseLoop sub-function.
% EXAMPLE:  [~,x] = xASL_Iteration(x,'xASL_module_ASL');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

    % General settings
    dbSettings.x                     = x;
    dbSettings.x.RERUN               = false;
    dbSettings.x.MUTEXID             = moduleName;
    dbSettings.x.LockDir             = ['<ROOT>/lock/' moduleName];
    
	% Set the dryRun field
	if nargin < 3 || isempty(dryRun)
		if ~isfield(x,'dryRun') || isempty(x.dryRun)
			dbSettings.dryRun = 0;
		else
			dbSettings.dryRun = x.dryRun;
		end
	else
		dbSettings.dryRun = dryRun;
	end
	
	% Set the stopAfterErrors field
	if nargin < 4 || isempty(stopAfterErrors)
		if ~isfield(x,'stopAfterErrors') || isempty(x.stopAfterErrors)
			dbSettings.stopAfterErrors = Inf;
		else
			dbSettings.stopAfterErrors = x.stopAfterErrors;
		end
	else
		dbSettings.stopAfterErrors = stopAfterErrors;
    end
	    
    if nargin<5 || isempty(SelectedSubjects)
        try
            SelectedSubjects = x.SUBJECTS;
        catch ME
            warning('Loading subjects didnt work, verify that ExploreASL was correctly initialized with the correct DataParameter file');
            fprintf('%s\n', ME.message);
        end
    end
    
    dbSettings.jobfn                 = str2func(moduleName);
    fprintf(['Running ' moduleName ' ...   ']);

    if  strcmp(moduleName(1:12),'xASL_module_')
        ModName     = moduleName(13:end);
    elseif strcmp(moduleName(1:7),'module_') % backward compatibility
        ModName     = moduleName(8:end);
    else
        ModName     = moduleName;
    end

    % Specific settings
    if  strcmp(ModName,'DARTEL') || strcmp(ModName,'LongReg')
        dbSettings.x.LockDir         = [dbSettings.x.LockDir '_' x.P.STRUCT ];
	end
    
	if ~isempty(regexp(ModName,'(Struct|ASL|func|LongReg|dwi)'))
		dbSettings.sets.SUBJECT      = SelectedSubjects; % x.SUBJECTS
		dbSettings.x.SUBJECTDIR      = '<ROOT>/<SUBJECT>';
		dbSettings.x.LockDir         = [dbSettings.x.LockDir '/<SUBJECT>'];
	end
	
	if ~isempty(regexp(ModName,'(ASL|func|dwi)'))
		dbSettings.sets.SESSION      = x.SESSIONS;
		dbSettings.x.MUTEXID         = [dbSettings.x.MUTEXID '_<SESSION>'];
		dbSettings.x.SESSIONDIR      = '<ROOT>/<SUBJECT>/<SESSION>';
	end
	
	if ~isempty(regexp(ModName, '(DARTEL|Population|Analyze)'))
		dbSettings.diaryFile               = ['<PopDir>/' moduleName '.log'];
	elseif ~isempty(regexp(ModName,'(ASL|func|dwi)'))
		dbSettings.diaryFile               = ['<SESSIONDIR>/' moduleName '.log'];
	else
		dbSettings.diaryFile               = ['<SUBJECTDIR>/' moduleName '.log'];
	end
	
	[bAborted, xOut] = runIteration(dbSettings);
end

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
%                   db.stopAfterErrors - optional: number of allowed errors before job iteration is stopped
% OUTPUT:
%   bAborted        - Report if the run was aborted
%   xOut            - x-struct on the output

    % Admin
    bAborted = false;

    % Check if parallel computing toolbox is available
    poolsize = 0;
    if license('test','Distrib_Computing_Toolbox') && license('checkout','Distrib_Computing_Toolbox')
        if verLessThan('matlab','8.2') 
            % < R2013b
            if exist('matlabpool','file')
                poolsize = matlabpool('size');
            end
        else
            % >= R2013b
            if  exist('gcp','file') % bugfix, some even have the Distrib_Computing_Toolbox, but don't have the Parallel_Computing_Toolbox
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if ~isempty(poolobj)
                    poolsize = poolobj.NumWorkers;
                end
            end
        end
    end
    
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
    if     ~isfield(db,'sets') % add a dummy set to be able to run the script using parfor without sets
            db.sets.DUMMY = { '' };
    elseif  numel(fieldnames(db.sets))<1
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
	
	stopAfterErrors = db.stopAfterErrors;
   
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

    % Check if a job has to be converted for parallel execution...
    if poolsize>0 && setCount>0
        % find largest set and use it to submit parallel
        [maxSetSize, maxIndex] = max(N);
        if maxSetSize>1
            setName = setNames{maxIndex};
            setValues = db.sets.(setName);
            parfor iPar=1:maxSetSize
                % copy set value as symbol and remove the set before launching a separate job
                dbi = db;
                dbi.sets = rmfield(dbi.sets,setName);
                if isfield(dbi.x,setName)
                    error('xASL_Iteration:symbolExists', 'Dynamic parameter %s already exists as symbol!',setName);
                end
                dbi.parindex = iPar; % add a variable to indicate parallel execution
                dbi.parcount = maxSetSize; % add a variable to indicate parallel execution
                fprintf('submitting job #%d: %s\n', iPar, setValues{iPar});
                dbi.x.(setName) = setValues{iPar}; % copy the indexed value from the set as regular (new) symbol
                % and finally submit...
                runIteration(dbi); % ignore result
            end
            return
        end
    end
    
    % Then loop through all permutations
    bAborted = false;
    iIter = 0; % just count the number of iterations
    CountErrors     = 0;
    while ~bAborted    % do not use for iIter=1:prod(N) anymore because we now support dynamic sets with a fluctuating number of values
        % Iterate to the next permutation by increment one of the indices (possibly resetting others)
        Iprev=I; % remember the previous indices because dynamic sets have to be re-eavaluated when the 'parent' set changes
        I = NextIter(I,N);
        if (sum(I)==0), break; end
        iIter = iIter+1;
        
        % Mark start time
        tocID = tic;
     
        % copy x-struct from user defined db
		x     = db.x;
        
        % Backward compatibility
        Flds    = {'ROOT' 'DARTELDIR'};

        if ~isfield(x,'D')
            x.D = struct;
        end
        
        for iL=1:length(Flds)
            if ~isfield(x.D,Flds{iL}) && isfield(x,Flds{iL})
                x.D.(Flds{iL})    = x.(Flds{iL});
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
                        fprintf('\n%s\n',repmat('+',1,72)); % just draw a separator line
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
        if  exist(xASL_adm_ReplaceSymbols(fullfile(x.LockDir,x.MUTEXID,'999_ready.status'),x),'file')
            AlreadyProcessed    = 1;
        else
            AlreadyProcessed    = 0;
        end
        
        
        % optional diary logfile of command window

        if ~isempty(diaryFile)
            diaryFileEx = xASL_adm_ReplaceSymbols(diaryFile,x);
                % NB. Don't overwrite diaryFile itself because it might be expanded
                % with different values in the next iteration.
            status = xASL_adm_CreateDir(fileparts(diaryFileEx),2); % make sure folder exists, but limit to 2 missing levels
            if status>=0
%                 fprintf('Opening diary log %s\n',diaryFileEx);
                diary(diaryFileEx);
            else
                error('xASL_Iteration:cannotCreateDirectory', 'Cannot create folder for diary file: %s!',diaryFileEx);
            end
        end
        
        % Some feedback about this iteration (after opening diary log)
        if ~AlreadyProcessed
            fprintf('\n%s\n',repmat('+',1,72)); % just draw a separator line
            fprintf('%s\n',['Starting ' x.MUTEXID ' at ' datestr(now)]);
            fprintf('\n');
        end
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
        
        % Start the job with all expanded x
        pwd_org = pwd; % remember the current working dir, so it can be restored if job fails
        try
            if  db.dryRun
                fprintf('\nDRYRUN: skipping %s\n',func2str(jobfn));
            else 
				[result, x] = jobfn(x); % [result, x] = jobfn(x, job_ex);
                if ~result
                    error('xASL_Iteration:jobfnFailed', 'job function returned failed');
                end
            end
        catch ME
            fprintf('\nERROR: Job iteration terminated!\n');
            ME.getReport()
            CountErrors     = CountErrors+1;
            if stopAfterErrors>0
                stopAfterErrors = stopAfterErrors -1;
                fprintf('\nCONT: but continue with next iteration!\n');
            else
                % assume that there is a serious error in the template or
                % batch if it occurs at the very first run
                bAborted        = true; % bailing out
                % this should stop the inner DB loop in this script
                
            end
        end
        cd(pwd_org); % always go back to where we came from
        
        % Print CPU time
        if ~AlreadyProcessed
            fprintf('\nJob-iteration %i stopped at %s and took %u seconds\n',iIter,datestr(now),ceil(toc(tocID)));
            fprintf('%s\n',repmat('+',1,72)); % just draw a separator line
        end
        
        % Close logfile
        if ~isempty(diaryFile)
            diary off
        end
        xASL_TrackProgress(I,N);
    end % next DB iteration
    fprintf('\n%s\n',[x.MUTEXID ' completed']);
%     if  CountErrors>0 
%         bAborted = true; % don't start other modules, this shouldn"t involve the inner DB loop in this script, 
%         % but should avoid running other modules
%     end    
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
