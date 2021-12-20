function xASL_init_PrintUserFeedback(x, outputArguments, currentState)
%xASL_init_PrintUserFeedback Print user feedback
%
% FORMAT: xASL_init_PrintUserFeedback(x)
% 
% INPUT:
%   x               - ExploreASL x structure (STRUCT, REQUIRED)
%   outputArguments - nargout of ExploreASL_Master (INTEGER, REQUIRED)
%   currentState    - State before/after processing (before = 0, after = 1, INTEGER, REQUIRED)
%
% OUTPUT:        n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Print user feedback.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        xASL_init_PrintUserFeedback(x);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Print user feedback
    if currentState==0
        % If this is true, we break here
        if x.opts.bProcessData && ~isdeployed && x.opts.bPause
            fprintf('%s\n','Press any key to start processing & analyzing');
            fprintf('Please ensure you have a read-only copy of your original data as they may be overwritten\n');
            fprintf('%s\n','Or press CTRL/command-C to cancel...  ');
            pause;
        end
        % Tell user to add output argument, if data loading was requested
        % without processing
        if ~x.opts.bProcessData && x.opts.bLoadData && outputArguments==0
            fprintf('Data loading requested but no output structure defined...\n');
            fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
        end
    else
        % Only print final feedback if data import or processing was performed
        if x.opts.bImportData || x.opts.bProcessData
            fprintf('Many thanks for using <a href="https://github.com/ExploreASL" rel="nofollow">ExploreASL</a>, ');
            fprintf('please don''t forget to cite <a href="https://pubmed.ncbi.nlm.nih.gov/32526385/" rel="nofollow">https://pubmed.ncbi.nlm.nih.gov/32526385/</a>.\n');
            fprintf('Note that ExploreASL is a collaborative effort.\n');
            fprintf('Therefore, please don''t hesitate to contribute by feedback, adding code snippets, or clinical experience!\n');
        end
    end


end


