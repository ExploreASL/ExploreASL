function xASL_init_printSettings(x)
%xASL_init_printSettings Print chosen settings
%
% FORMAT: xASL_init_printSettings(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Print chosen settings.
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright 2015-2021 ExploreASL
    
    %% Printing
    fprintf('==================================== ExploreASL Settings =====================================\n');
	% Dataset root
    if length(x.opts.DatasetRoot)>70
        fprintf('Dataset Root        ...%s\n', x.opts.DatasetRoot(end-70:end));
    else
        fprintf('Dataset Root        %s\n', x.opts.DatasetRoot); 
	end
	
	% Import modules
	textPrint = 'Import Modules      ';
	if x.opts.ImportModules(1)==1
		textPrint = [textPrint 'DCM2NII '];
	end
	if x.opts.ImportModules(2)==1
		textPrint = [textPrint 'NII2BIDS '];
	end
	if x.opts.ImportModules(3)==1
		textPrint = [textPrint 'DEFACE '];
	end
	if x.opts.ImportModules(4)==1
		textPrint = [textPrint 'BIDS2LEGACY'];
	end
	fprintf([textPrint '\n']);
	
    % Process modules
	textPrint = 'Process Modules     ';
	if x.opts.ProcessModules(1)==1
		textPrint = [textPrint 'STRUCTURAL '];
	end
	if x.opts.ProcessModules(2)==1
		textPrint = [textPrint 'ASL '];
	end
	if x.opts.ProcessModules(3)==1
		textPrint = [textPrint 'POPULATION '];
	end
	fprintf([textPrint '\n']);
	
	% Pause before processing
	if x.opts.bPause==1
        fprintf('bPause              %s\n', 'True');
    else
        fprintf('bPause              %s\n', 'False');
	end
	
	% Worker numbers
    fprintf('iWorker             %d\n', x.opts.iWorker);
    fprintf('nWorkers            %d\n', x.opts.nWorkers);
    fprintf('==============================================================================================\n');

end


