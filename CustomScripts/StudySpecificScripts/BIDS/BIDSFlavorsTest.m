%% This functions makes a copy of the source and reference flavors database and runs the conversion and verifies
% the results of the conversion with respect to the reference
% Set the paths
srcPath = '/pet/projekte/asl/data/BIDS/';

flavorsPath  = fullfile(srcPath,'Flavors');
conversionPath = fullfile(srcPath,'TmpConversion');
referencePath = fullfile(srcPath,'TmpReference');
%% Copy the source data apart
if ~exist(conversionPath,'dir')
	mkdir(conversionPath);
end

fList = xASL_adm_GetFileList(flavorsPath,[],'List',[],1);
for iList = 1:length(fList)
	if ~exist(fullfile(conversionPath,fList{iList}),'dir')
		mkdir(fullfile(conversionPath,fList{iList}));
	end
	if ~exist(fullfile(conversionPath,fList{iList},'sourcedata'),'dir')
		mkdir(fullfile(conversionPath,fList{iList},'sourcedata'));
	end
	system(['cp -r ' fullfile(flavorsPath,fList{iList},'sourcedata','*') ' ' fullfile(conversionPath,fList{iList},'sourcedata')]);
	system(['cp -r ' fullfile(flavorsPath,fList{iList},'imPar.json') ' ' fullfile(conversionPath,fList{iList},'imPar.json')]);
	system(['cp -r ' fullfile(flavorsPath,fList{iList},'studyPar.json') ' ' fullfile(conversionPath,fList{iList},'studyPar.json')]);
end

%% Create a copy of the reference data
if ~exist(referencePath,'dir')
	mkdir(referencePath);
end

fList = xASL_adm_GetFileList(flavorsPath,[],'List',[],1);
for iList = 1:length(fList)
	if ~exist(fullfile(referencePath,fList{iList}),'dir')
		mkdir(fullfile(referencePath,fList{iList}));
	end
	if ~exist(fullfile(referencePath,fList{iList},'rawdata'),'dir')
		mkdir(fullfile(referencePath,fList{iList},'rawdata'));
	end
	system(['cp -r ' fullfile(flavorsPath,fList{iList},'rawdata','*') ' ' fullfile(referencePath,fList{iList},'rawdata')]);
end


%% Run the conversion
xASL_bids_TestBidsConversion(conversionPath,referencePath,1,0);

%% Run the comparison
xASL_bids_TestBidsConversion(conversionPath,referencePath,0,1);
