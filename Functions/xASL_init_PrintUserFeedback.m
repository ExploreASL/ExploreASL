function xASL_init_PrintUserFeedback(x)
%xASL_init_PrintUserFeedback Print user feedback
%
% FORMAT: xASL_init_PrintUserFeedback(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
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
    if ~x.opts.bProcessData || x.opts.bOnlyLoad
        if x.opts.bOnlyLoad && nargout==0
            warning('Data loading requested but no output structure defined');
            fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
        end
        return; % skip processing
    elseif ~isdeployed && x.opts.bPause % if this is true, we skip the break here
        fprintf('%s\n','Press any key to start processing & analyzing');
        fprintf('Please ensure you have a read-only copy of your original data as they may be overwritten\n');
        fprintf('%s\n','Or press CTRL/command-C to cancel...  ');
        pause;
    end


end


