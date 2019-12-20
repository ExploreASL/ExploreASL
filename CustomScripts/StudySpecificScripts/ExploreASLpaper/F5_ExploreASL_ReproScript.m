mainPath = '/pet/projekte/asl/data/ExploreASLrepro';

studyName = {'EPAD','Novice','Sleep'};

versionName = 'NEW';'OLD';

matlabName = 'Lin18b';

% Repeat twice
for rS = 1:2
	% All studies
	for iS = 1:length(studyName)
		locDir = fullfile(mainPath,[studyName{iS} '_' matlabName '_' versionName '_' num2str(rS)],'DataPar.m')
		ExploreASL_Master(locDir,1,1);
	end
end
