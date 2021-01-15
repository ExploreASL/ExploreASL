%% Take output of ASLCheck from several folders, export to a new folder, then merge coronal and axial view,
%  randomly order and anonymize and create a table for evaluation

% First start with folder names

basicPath = '/pet/projekte/asl/data/Craniosynostosis';
methodPaths = {'analysisASLonlyAffineFinal'
	           'analysisASLonlyDCTFinal'
			   'analysisASLonlyRigidFinal'
			   'analysisSegm1y'};
		   
viewNames = {'Cor', 'Tra'};
methodRegexpNormal = {'Reg_pWMPV_pWM_.*.jpg$'
	                  'Reg_pWMPV_pWM_.*.jpg$'
					  'Reg_pWMPV_pWM_.*.jpg$'
					  'Temp_pWM_.*rc2T1_ASL_res.jpg$'};
methodRegexpContour = {'Reg_pWM_.*pWMC.jpg$'
	                   'Reg_pWM_.*pWMC.jpg$'
					   'Reg_pWM_.*pWMC.jpg$'
					   'Temp_pWM_.*pWMC.jpg$'};

outputPath = fullfile(basicPath,'QCcomparison');

%% Create all the directories
mkdir(outputPath);

% Path for the merged files
for iPath = 1:length(methodPaths)
	mkdir(fullfile(outputPath,'merged',methodPaths{iPath}));
end

% Path for the anonymized files
mkdir(fullfile(outputPath,'anonymized','normal'));
mkdir(fullfile(outputPath,'anonymized','contour'));

%% Copy all the files and merge coronal+axial

for iMethod = 1:length(methodPaths)
	% List all the relevant files (normal and contour)
	fListNormal  = xASL_adm_GetFileList(fullfile(basicPath,methodPaths{iMethod},'Population','ASLCheck'),...
	                     ['^' viewNames{1} '_' methodRegexpNormal{iMethod}], 'List', [], 0);
	fListContour = xASL_adm_GetFileList(fullfile(basicPath,methodPaths{iMethod},'Population','ASLCheck'),...
	                     ['^' viewNames{1} '_' methodRegexpContour{iMethod}], 'List', [], 0);
					 
	% From all the 'coronal' files, create also the axial
	for iFile = 1:length(fListNormal)
		fListNormalAxial{iFile} = [viewNames{2} fListNormal{iFile}((length(viewNames{1})+1):end)];
		fListContourAxial{iFile} = [viewNames{2} fListContour{iFile}((length(viewNames{1})+1):end)];
	end
	
	% Load all the files, merge them and save merged to the correct directory
	for iFile = 1:length(fListNormal)
		imA = imread(fullfile(basicPath,methodPaths{iMethod},'Population','ASLCheck',fListNormal{iFile}));
		imB = imread(fullfile(basicPath,methodPaths{iMethod},'Population','ASLCheck',fListNormalAxial{iFile}));
		imC = zeros(max(size(imA,1),size(imB,1)),size(imA,2)+size(imB,2),3);
		imC(1:size(imA,1),1:size(imA,2),1:3) = imA;
		imC(1:size(imB,1),size(imA,2)+1:end,1:3) = imB;
		imwrite(imC/256,fullfile(outputPath,'merged',methodPaths{iMethod},fListNormal{iFile}((length(viewNames{1})+2):end)));
		
		imA = imread(fullfile(basicPath,methodPaths{iMethod},'Population','ASLCheck',fListContour{iFile}));
		imB = imread(fullfile(basicPath,methodPaths{iMethod},'Population','ASLCheck',fListContourAxial{iFile}));
		imC = zeros(max(size(imA,1),size(imB,1)),size(imA,2)+size(imB,2),3);
		imC(1:size(imA,1),1:size(imA,2),1:3) = imA;
		imC(1:size(imB,1),size(imA,2)+1:end,1:3) = imB;
		imwrite(imC/256,fullfile(outputPath,'merged',methodPaths{iMethod},fListContour{iFile}((length(viewNames{1})+2):end)));
	end
end

%% Make reordering, rename to output and create pseudonymization tables

% Go over all methods
for iMethod = 1:length(methodPaths)
	% List all the relevant files (normal and contour)
	fListNormal  = xASL_adm_GetFileList(fullfile(outputPath,'merged',methodPaths{iMethod}),...
	                     ['^' methodRegexpNormal{iMethod}], 'List', [], 0);
	fListContour = xASL_adm_GetFileList(fullfile(outputPath,'merged',methodPaths{iMethod}),...
	                     ['^' methodRegexpContour{iMethod}], 'List', [], 0);
    
	% Make a list of all files, including a list of directories and methods in the same order
	if iMethod == 1
		fListNormalTotal = fListNormal;
		fListContourTotal = fListContour;
		fMethodsTotal = repmat({methodPaths{iMethod}},length(fListNormal),1);
		fDirsTotal = repmat({fullfile(outputPath,'merged',methodPaths{iMethod})},length(fListNormal),1);
	else
		fListNormalTotal = [fListNormalTotal; fListNormal];
		fListContourTotal = [fListContourTotal; fListContour];
		fMethodsTotal = [fMethodsTotal; repmat({methodPaths{iMethod}},length(fListNormal),1)];
		fDirsTotal = [fDirsTotal; repmat({fullfile(outputPath,'merged',methodPaths{iMethod})},length(fListNormal),1)];
	end
end

% Random number permutation
listPerm = randperm(length(fListNormalTotal))';

% Save a table - file name, methods name, number
resultsTable = fMethodsTotal;
resultsTable(:,2) = fListNormalTotal;
resultsTable(:,3) = fListContourTotal;
resultsTable(:,4) = cellstr(num2str(listPerm,'%.4d'));
xASL_tsvWrite(resultsTable,fullfile(outputPath,'pseudonymTable.tsv'), 1);

% Copy the file name to the order for both contour and normal
for iFile = 1:length(fListNormalTotal)
	xASL_Copy(fullfile(fDirsTotal{iFile},fListNormalTotal{iFile}),fullfile(outputPath,'anonymized','normal',[resultsTable{iFile,4} '.jpg']),1,1);
	xASL_Copy(fullfile(fDirsTotal{iFile},fListContourTotal{iFile}),fullfile(outputPath,'anonymized','contour',[resultsTable{iFile,4} '.jpg']),1,1);
end
