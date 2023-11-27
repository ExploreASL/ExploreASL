function [x] = xASL_init_Parallelization(x)
%xASL_init_Parallelization If running parallel, select cases for this worker
%
% FORMAT: [x] = xASL_init_Parallelization(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x          - ExploreASL x structure
%                         
%               Parallelization is allowed here by calling ExploreASL different times,
%               where it divides the subjects/images for processing across the nWorkers,
%               using iWorker as the reference for the part that the current ExploreASL
%               call will process. This requires having a Matlab license that can be
%               started multiple times on a server, or alternatively running the
%               ExploreASL compilation, and doesn't require the Matlab parallel toolbox.
%
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_Parallelization(x);
% __________________________________
% Copyright (c) 2015-2024 ExploreASL

% ------------------------------------------------------------------------------------------------
%% 2) Parallelization: If running parallel, select cases for this worker
if x.opts.nWorkers>1
    nSubjPerWorker = x.dataset.nTotalSubjects/x.opts.nWorkers; % ceil to make sure all subjects are processed

    % e.g., if nWorkers=3 & nTotalSubjects=10
    % iWorker 1 does [1 2 3]
    % iWorker 2 does [4 5 6 7]
    % iWorker 3 does [8 9 10]

    iStartSubject = round((x.opts.iWorker-1)*nSubjPerWorker+1);
    iEndSubject = min( round(x.opts.iWorker*nSubjPerWorker), x.dataset.nTotalSubjects);

    if iStartSubject>x.dataset.nTotalSubjects
        warning('Closing down this worker, had too many workers');
        exit;
    end
    
    % Adapt SUBJECTS
    x.dataset.TotalSubjects = x.dataset.TotalSubjects(iStartSubject:iEndSubject);
    x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);

    
    fprintf(['I am worker ' num2str(x.opts.iWorker) '/' num2str(x.opts.nWorkers) '\n']);
    fprintf(['I will process subjects ' num2str(iStartSubject) '-' num2str(iEndSubject) '\n']);
end


end