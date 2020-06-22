function xASL_TrackProgress(iCurrent, iMax)
% Performs an in-screen counting in percentages
%
% FORMAT: xASL_TrackProgress(iCurrent[, iMax])
%
% INPUT:
%   iCurrent  - The part of the work done. Either between 0 and 100, or smaller than iMax when iMax is given.
%   iMax      - Optional - total number of iterations.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Counts the percentage of the work done and display on the screen. Either iCurrent of iMax are given. Or only the percentages are given.
% EXAMPLE: xASL_TrackProgress(23,80);
%          xASL_TrackProgress(0.5);
%          xASL_TrackProgress(80,80);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% NOTE FOR DEBUGGING:
% When this function is breaked manually/crashes,
% the diary can be off and diary file lost, so the diary file needs to be
% set manually on again. To fix this, first get the diary file:
% PathDiary = get(0, 'DiaryFile');
% And reset the diary back to this after the break/crash:
% diary(PathDiary);
% __________________________________
% Copyright 2015-2020 ExploreASL

    HasDiary = strcmp(get(0, 'Diary'), 'on');
    
    if HasDiary
        % we don't want to print this to the diary, if there is any
        diary off;
    end

    if nargin < 2
		% If iMax is not given, then we assume that iCurrent is between 0 and 1
		% But we set it to the maximum of 1 and iCurrent to avoid having more than 100%
		iMax = max(100,iCurrent);
	end
	
	% Avoid exceeding 100%
	if iCurrent > iMax
		iCurrent = iMax;
	end
	
	% Avoid negative percentages
	if iCurrent < 0
		iCurrent = 0;
	end
	
	% Avoid division by zero and negative percentages
	if iMax <= 0
		iMax = 100;
	end
	
	if isunix && ~usejava('desktop')
		fprintf('\r%02d',floor(100*iCurrent/iMax));
	else
		fprintf('\b\b\b%02d',floor(100*iCurrent/iMax));
	end
	fprintf('%s','%');

    if HasDiary
        diary on; % set the diary back on again
    end
    
end
