mainPath = '/pet/projekte/asl/data/ExploreASLreproAll';

studyName = {'EPAD','Novice','Sleep'};

versionName = {'NEW','OLD'};

matlabName = {'Lin18b','Win18b','Win16b','Win15a','WinLin16b','Mac14b','Mac16b'};
pathXasl = '/home/janpetr/code/ExploreASL';
% Mat1, try1, ver1, Mat2, try2, ver2
compareVec = [1,1,1,1,2,1; % compare twice the same on linux NEW 18b
	          1,1,1,2,1,1; % Compare linux and windows 18b on NEW
			  1,1,2,2,1,2; % Compare linux and windows 18b on OLD
			  1,1,2,1,2,2;% compare twice the same on linux OLD 18b
		      2,1,2,3,1,2; % Compare windows 18bvs 16b on OLD
		      2,1,1,3,1,1; % Compare windows 18bvs 16b on NEW
		      2,1,2,4,1,2; % Compare windows 18bvs 15a on OLD
		      2,1,1,4,1,1; % Compare windows 18bvs 15a on NEW
			  5,1,1,2,1,1; % Compare WinLin16b with Win18b NEW
		      5,1,2,2,1,2; % Compare WinLin16b with Win18b OLD
			  5,1,1,1,1,1; % Compare WinLin16b with Lin18b NEW
			  5,1,2,1,1,2; % Compare WinLin16b with Lin18b OLD
			  7,1,1,1,1,1; % Compare Mac16b with Lin18b NEW
			  6,1,1,1,1,1; % Compare Mac14b with Lin18b NEW
		      7,1,2,1,1,2; % Compare Mac16b with Lin18b OLD
			  6,1,2,1,1,2; % Compare Mac14b with Lin18b OLD
		      7,1,1,2,1,1; % Compare Mac16b with Win18b NEW
		      7,1,2,2,1,2; % Compare Mac16b with Win18b OLD
			  4,1,1,1,1,1; % Compare Win15a with Lin18b NEW
		      4,1,2,1,1,2; % Compare Win15a with Lin18b OLD
			  7,2,1,7,3,1]; % Compare Mac16b repeat NEW

locDir = pwd;

% All studies
for iS = 1:length(studyName)
	for iC = size(compareVec,1)
		nameRef = [studyName{iS} '_' matlabName{compareVec(iC,1)} '_' versionName{compareVec(iC,3)} '_' num2str(compareVec(iC,2))];
		nameSrc = [studyName{iS} '_' matlabName{compareVec(iC,4)} '_' versionName{compareVec(iC,6)} '_' num2str(compareVec(iC,5))];

		% Rename forward
		system(['mv ' fullfile(mainPath,nameRef)  ' ' fullfile(pathXasl,'TestDataSet_Processed')]);
		system(['mv ' fullfile(mainPath,nameSrc) ' ' fullfile(mainPath,'TestDataSet_Result')]);

		% Run the comparison script
		cd(pathXasl);
		ExploreASL_Master(fullfile(mainPath,'TestDataSet_Result','DataPar.m'),1,1);
		cd(locDir);

		% Copy out the result
		system(['cp ' fullfile(mainPath,'TestDataSet_Result','RMS_Reproducibility.mat') ' ' fullfile(mainPath,['RMS_Repro_' nameRef '_' nameSrc '.mat'])]);

		% Rename back
		system(['mv ' fullfile(pathXasl,'TestDataSet_Processed') ' ' fullfile(mainPath,nameRef)]);
		system(['mv ' fullfile(mainPath,'TestDataSet_Result') ' ' fullfile(mainPath,nameSrc)]);
	end
end
