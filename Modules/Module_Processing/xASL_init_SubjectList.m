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
    subjectFolder = x.dir.Rawdata;
else
    subjectFolder = x.dir.xASLDerivatives;
end


% ------------------------------------------------------------------------------------------------
%% 2. Advanced parameter to enforce a subject list, without querying folder names
if isfield(x.dataset,'ForceInclusionList')
    % This is an option if you want to select subjects yourself,
    % instead of using all the subjects that comply with the regular expression
    x.dataset.TotalSubjects = x.dataset.ForceInclusionList(:);
    warning('Using custom list of subjects, on your own risk');

elseif x.opts.bReadRawdata

    %% 3. Load the BIDS structure of all subjects & translate to ExploreASL legacy
    x.modules.bids2legacy.BIDS = bids.layout(subjectFolder);
    
    % Now we translate the BIDS naming convention to ExploreASL's legacy name convention
    for iSubj=1:length(x.modules.bids2legacy.BIDS.subjects)
        subjectName = x.modules.bids2legacy.BIDS.subjects(iSubj).name;
        visitName = x.modules.bids2legacy.BIDS.subjects(iSubj).session;
        if ~strcmp(visitName(1:4), 'ses-')
            warning(['Illegal BIDS session definition for ' subjectName '_' visitName]);
        else
            visitName = visitName(5:end);
        end
        
        x.dataset.TotalSubjects{iSubj} = [subjectName '_' visitName];
    end

else
    %% 4. Create legacy subject list from folders

   % Manage subjectRegexp
   if ~isfield(x.dataset, 'subjectRegexp')
       warning('Please define the subjectRegexp of your dataset, defaulting to BIDS definition');
       subjectRegexp = '^sub-.*$';
   else
       subjectRegexp = x.dataset.subjectRegexp;
   end
    
    % Do escapes and fixes
    % First escape double escapes
    subjectRegexp = strrep(subjectRegexp,'\\','\');
    % add the visit-postfix as option
    subjectRegexp = strrep(subjectRegexp,'$','');
    subjectRegexp = [subjectRegexp '(|_\d+)$'];

    %% 5. Then load subjects
    x.dataset.TotalSubjects = sort(xASL_adm_GetFileList(subjectFolder, subjectRegexp, 'List', [0 Inf], true)); % find dirs

end

x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);


% ------------------------------------------------------------------------------------------------
%% 6. Manage exclusions

if ~isfield(x,'SUBJECTS')
    x.SUBJECTS = '';
end
x.dataset.nSubjects = length(x.SUBJECTS);

if isempty(x.SUBJECTS)
    fprintf(2,'No subjects found...\n');
    fprintf(2,'Please check the sub_regexp in your Data Parameter File.\n');
    fprintf(2,'This should match with the subject folders inside the provided subjectFolder\n');
    fprintf(2,'This was %s ...\n\n', x.opts.subjectFolder);
    error('No subjects defined, x.SUBJECTS was empty...');
elseif isfield(x.dataset, 'exclusion')
    warning('Exclusion list may not work with longitudinal registration');
end

if     ~isfield(x.dataset,'exclusion')
        x.dataset.exclusion = {''};
        nExclusion = 0;
elseif  isempty(x.dataset.exclusion)
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

x.dataset.ExcludedSubjects = '';
x.SUBJECTS = '';
for iSubject=1:x.dataset.nTotalSubjects
    excl=0;
    for j=1:nExclusion % Check if subject should be excluded
        if  strcmp(x.dataset.TotalSubjects{iSubject},x.dataset.exclusion{j})
            excl=1;
            x.dataset.ExcludedSubjects{end+1} = x.dataset.TotalSubjects{iSubject};
        end
    end
    
    if ~isempty(regexp(x.dataset.TotalSubjects{iSubject},'^(log|dartel|lock|Population)$'))
         % This is not a subject but a pipeline folder
         ListNoPipelineDir(iSubject) = 0;
    elseif ~excl
         % include subject if not to be excluded, and if it doesn't have
         % same name as pipeline directories     
         x.SUBJECTS{end+1}  = x.dataset.TotalSubjects{iSubject};
         x.dataset.TotalInclusionList( (iSubject-1)*x.dataset.nSessions+1:iSubject*x.dataset.nSessions,1) = 1;
         ListNoPipelineDir(iSubject) = 1;
    else
        x.dataset.TotalInclusionList( (iSubject-1)*x.dataset.nSessions+1:iSubject*x.dataset.nSessions,1) = 0;
        ListNoPipelineDir(iSubject) = 1;
    end
end

% Default definition of ListNoPipelineDir
if ~exist('ListNoPipelineDir','var')
    ListNoPipelineDir = [];
end

% Remove pipeline dirs from list
if sum(ListNoPipelineDir)>0 % if there are any subjects
    x.dataset.TotalSubjects = x.dataset.TotalSubjects(logical(ListNoPipelineDir));
end

if ~isfield(x,'SUBJECTS')
    fprintf('%s\n','No subjects found');
    x.SUBJECTS = [];
elseif isempty(x.SUBJECTS)
    fprintf('%s\n','No subjects found');
end

x.dataset.nSubjects = length(x.SUBJECTS);
x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);
x.dataset.nExcluded = x.dataset.nTotalSubjects - x.dataset.nSubjects;


end