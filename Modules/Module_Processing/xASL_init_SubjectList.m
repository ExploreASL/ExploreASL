function [x] = xASL_init_SubjectList(x)
%xASL_init_SubjectList Generate subject & visit list and adding to x struct
%
% FORMAT: [x] = xASL_init_SubjectList(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x          - ExploreASL x structure
%                         
% 1. Define folder to read subjects from
% 2. Advanced parameter to enforce a subject list, without querying folder names
% 3. Load the BIDS structure of all subjects & translate to ExploreASL legacy
% 4. Create legacy subject list from folders
% 5. Then load subjects
% 6. Manage exclusions
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_SubjectList(x);
% __________________________________
% Copyright (c) 2015-2024 ExploreASL


%% 1. Define folder to read subjects from
if isfield(x.opts, 'subjectFolder') && ~exist(x.opts.subjectFolder, 'dir')
    warning(['subjectFolder' x.opts.subjectFolder ' is no valid directory']);
    fprintf('Using default folder instead\n');
elseif isfield(x.opts, 'subjectFolder')
    subjectFolder = x.opts.subjectFolder;
elseif x.opts.bReadRawdata
    subjectFolder = x.dir.RawData;
else
    subjectFolder = x.dir.xASLDerivatives;
end

x.ROOT = subjectFolder;

% ------------------------------------------------------------------------------------------------
%% 2. Advanced parameter to enforce a subject list, without querying folder names
x.dataset.TotalSubjects = cell(0);

if isfield(x.dataset,'ForceInclusionList')
    % This is an option if you want to select subjects yourself,
    % instead of using all the subjects that comply with the regular expression
    x.dataset.TotalSubjects = x.dataset.ForceInclusionList(:);
    warning('Using custom list of subjects, on your own risk');

    if ~isfield(x.dataset, 'subjectRegexp')
        x.dataset.subjectRegexp = '^sub-.*$';
        warning('Ensure that the ForceInclusionList is in accordance with the BIDS subjectRegexp');
        fprintf('%s\n', 'Otherwise, this will go wrong further down the road');
    end

    if isempty(x.dataset.TotalSubjects)
        error('No subjects found, check x.dataset.ForceInclusionList in dataPar.json');
    end

elseif x.opts.bReadRawdata

    %% 3. Load the BIDS structure of all subjects & translate to ExploreASL legacy
    x.modules.bids2legacy.BIDS = bids.layout(subjectFolder);
    
    % Now we translate the BIDS naming convention to ExploreASL's legacy name convention
    for iSubj=1:length(x.modules.bids2legacy.BIDS.subjects)
        subjectName = x.modules.bids2legacy.BIDS.subjects(iSubj).name;
        visitName = x.modules.bids2legacy.BIDS.subjects(iSubj).session;
        if isempty(visitName)
            visitName = '1'; % we default to visit _1
        elseif ~strcmp(visitName(1:4), 'ses-')
            warning(['Illegal BIDS session definition for ' subjectName '_' visitName]);
        else
            visitName = visitName(5:end); % if we find a BIDS session, we add it as _n visit suffix in ExploreASL legacy format
        end
        
        x.dataset.TotalSubjects{iSubj} = [subjectName '_' visitName];
        x.dataset.subjectRegexp = '^sub-.*$'; % defaulting to BIDS subject regexp
    end

    if isempty(x.dataset.TotalSubjects)
        error('No subjects found, check your rawdata folder for BIDS compatibility');
    end

else
    %% 4. Create legacy subject list from folders

   % Manage subjectRegexp
   if ~isfield(x.dataset, 'subjectRegexp')
       warning('Please define the subjectRegexp of your dataset, defaulting to BIDS definition');
       x.dataset.subjectRegexp = '^sub-.*$';
   else    
        % Do escapes and fixes
        % First escape double escapes
        x.dataset.subjectRegexp = strrep(x.dataset.subjectRegexp,'\\','\');
        % add the visit-postfix as option
        x.dataset.subjectRegexp = strrep(x.dataset.subjectRegexp,'$','');
        x.dataset.subjectRegexp = [x.dataset.subjectRegexp '(|_\d+)$'];
   end

    % Then load subjects
    x.dataset.TotalSubjects = sort(xASL_adm_GetFileList(subjectFolder, x.dataset.subjectRegexp, 'List', [0 Inf], true)); % find dirs

    if isempty(x.dataset.TotalSubjects)
        fprintf('%s\n', 'No subjects found, subjectRegexp in dataPar.json should match with subjectFolder');
        fprintf(2,'This was %s ...\n\n', x.opts.subjectFolder);
        error('No subjects defined...');
    end

    % And check that this doesn't accidentally contain pipeline directories
    ListPipelineDir = logical([]);
    for iSubject=1:x.dataset.nTotalSubjects
        % This is not a subject but a pipeline folder
        ListPipelineDir(iSubject) = ~isempty(regexpi(x.dataset.TotalSubjects{iSubject},'^(log|dartel|lock|Population)$'));
    end

    % Remove pipeline dirs from list
    x.dataset.TotalSubjects = x.dataset.TotalSubjects(~ListPipelineDir);

    if sum(ListPipelineDir)>0
        warning('Somehow your subjectRegexp was too lax, try making your subjectRegexp as specific as possible');
    end

end

x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);
x.SUBJECTS = x.dataset.TotalSubjects;
x.dataset.nSubjects = length(x.SUBJECTS);

% ------------------------------------------------------------------------------------------------
%% 6. Manage exclusions

if isfield(x.dataset, 'exclusion')
    warning('Exclusion list may not work with longitudinal registration');
end

if ~isfield(x.dataset,'exclusion')
    x.dataset.exclusion = {''};
    nExclusion = 0;
elseif isempty(x.dataset.exclusion)
    x.dataset.exclusion = {''};    
    nExclusion = 0;
else
    nExclusion = length(x.dataset.exclusion);
end

if ~iscell(x.dataset.exclusion)
    x.dataset.exclusion = {x.dataset.exclusion};
end

if  nExclusion==1 && isempty(x.dataset.exclusion{1})
    nExclusion = 0;
end

x.dataset.ExcludedSubjects = cell(0);
x.SUBJECTS = cell(0);
for iSubject=1:x.dataset.nTotalSubjects
    excludeThisSubject = false;
    for j=1:nExclusion % Check if subject should be excluded
        if strcmp(x.dataset.TotalSubjects{iSubject},x.dataset.exclusion{j})
            excludeThisSubject = true;
            x.dataset.ExcludedSubjects{end+1} = x.dataset.TotalSubjects{iSubject};
        end
    end
    
    if ~excludeThisSubject
         x.SUBJECTS{end+1}  = x.dataset.TotalSubjects{iSubject};
    end
end

if isempty(x.SUBJECTS)
    fprintf('%s\n','No subjects found (after exclusion)');
end

x.dataset.nSubjects = length(x.SUBJECTS);
x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);
x.dataset.nExcluded = x.dataset.nTotalSubjects - x.dataset.nSubjects;


end