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
%               IMPORTANT: xASL_init_Iteration processes x.SUBJECTS (not
%               x.dataset.TotalSubjects). Both lists must be sliced for this worker,
%               otherwise every worker still iterates the full subject list and
%               collides on mutex locks.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_Parallelization(x);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


% ------------------------------------------------------------------------------------------------
%% Parallelization: If running parallel, select cases for this worker
if x.opts.nWorkers>1
    nSubjPerWorker = x.dataset.nTotalSubjects/x.opts.nWorkers; % ceil to make sure all subjects are processed

    % e.g., if nWorkers=3 & nTotalSubjects=10
    % iWorker 1 does [1 2 3]
    % iWorker 2 does [4 5 6 7]
    % iWorker 3 does [8 9 10]

    iStartSubject = round((x.opts.iWorker-1)*nSubjPerWorker+1);
    iEndSubject = min( round(x.opts.iWorker*nSubjPerWorker), x.dataset.nTotalSubjects);

    % Empty slice (iStart>iEnd) happens when nWorkers > nTotalSubjects and
    % rounding assigns this worker no subjects. Treat like "too many workers".
    if iStartSubject>x.dataset.nTotalSubjects || iStartSubject>iEndSubject
        warning('Closing down this worker, had too many workers');
        exit;
    end
    
    % Adapt both subject lists for this worker
    x.dataset.TotalSubjects = x.dataset.TotalSubjects(iStartSubject:iEndSubject);
    x.dataset.nTotalSubjects = length(x.dataset.TotalSubjects);

    % Keep x.SUBJECTS in sync. Iteration uses x.SUBJECTS, not TotalSubjects.
    if isfield(x, 'SUBJECTS') && ~isempty(x.SUBJECTS)
        x.SUBJECTS = x.SUBJECTS(ismember(x.SUBJECTS, x.dataset.TotalSubjects));
    else
        x.SUBJECTS = x.dataset.TotalSubjects;
    end
    x.dataset.nSubjects = length(x.SUBJECTS);
    
    fprintf(['I am worker ' num2str(x.opts.iWorker) '/' num2str(x.opts.nWorkers) '\n']);
    fprintf(['I will process subjects ' num2str(iStartSubject) '-' num2str(iEndSubject) '\n']);
end


end
