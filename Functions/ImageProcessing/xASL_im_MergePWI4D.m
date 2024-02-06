function [path_PWI4D, path_PWI4D_Pop] = xASL_im_MergePWI4D(x)
% Merge several PWIs based on the list in x.modules.asl.sessionsToMerge
PWI = []; % this becomes the total concatenated image matrix

fprintf('Concatenating PWI4D for sessions/runs (nVolumes):');

for iSession = 1:numel(x.modules.asl.sessionsToMerge)
	% Here we get the path of the first session, by replacing the foldername ASL_X
    pathCurrentPWI4D = replace(x.P.Path_PWI4D, x.modules.asl.sessionsToMerge{end}, x.modules.asl.sessionsToMerge{iSession});
    
	if xASL_exist(pathCurrentPWI4D, 'file')
		imPWI4Dcurrent = xASL_io_Nifti2Im(pathCurrentPWI4D);
        
        dim4(iSession) = size(imPWI4Dcurrent, 4);

        fprintf([' ' x.modules.asl.sessionsToMerge{iSession} ' (' xASL_num2str(dim4(iSession)) ')']);

		if iSession == 1
			% Take the first volume as it is
			PWI = imPWI4Dcurrent;
		else
			% For the following ones, check dimensions
			if isequal(size(imPWI4Dcurrent, 1:3), size(PWI, 1:3))
				% Here we concatenate (cat) over the 4rd dimension (4), PWI (total concatenated image matrix) with the new current PWI
                PWI = cat(4, PWI, imPWI4Dcurrent);
			else
				error(['Cannot concatenate sessions for subject ' x.SUBJECT ' as session '  x.modules.asl.sessionsToMerge{iSession} ' has a different matrix size from the previous sessions.']);
			end
		end
	else
		% If one of the sessions are missing, then we issue an error
		error(['Cannot concatenate sessions for subject ' x.SUBJECT ': session '  x.modules.asl.sessionsToMerge{iSession} ' is missing.']);
	end
end

fprintf('\n%s\n', ['Concatenated PWI4Ds into a new PWI4D with name XXXXX #1543 and ' xASL_num2str(sum(dim4)) ' volumes']);

end