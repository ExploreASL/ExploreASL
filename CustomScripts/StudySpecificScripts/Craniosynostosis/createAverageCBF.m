%% Set the paths
% Main path
rawDir = '/pet/projekte/asl/data/Craniosynostosis';

% Focus dataset
referenceName = 'analysisASLonlyDCTExtra';

% Directory for the results
resultsDir = 'templates';

groupName = {'Controle','Cranio'};
% Exclude datasets
excludeNames = {'Controle-0010','Cranio-0014','Cranio-0012',};
%% Create and save the average images

% Load the GM and WM segmentations
imGM = xASL_io_Nifti2Im(fullfile(rawDir,referenceName,'Population','rc1T1_Controle-0001.nii.gz'));
imWM = xASL_io_Nifti2Im(fullfile(rawDir,referenceName,'Population','rc2T1_Controle-0001.nii.gz'));

imTemplate = [];
for ii = 1:length(groupName)
	% List controls and cranios
	
	% Load all the subjects in the directory
	patientNameList = xASL_adm_GetFileList(fullfile(rawDir,referenceName,'Population'), ['^qCBF_' groupName{ii} '-\d{4}.*$'], 'FPList', [], 0);
	cnt(ii) = 0;
	for iF = 1:length(patientNameList)
		
		% Check if the subject shouldn't be excluded
		bCont = 1;
		for iN = 1:length(excludeNames)
			if ~isempty(strfind(patientNameList{iF},excludeNames{iN}))
				bCont = 0;
			end
		end
		
		% The subject should be included in the atlas
		if bCont
			% Load the subject
			imSubject = xASL_io_Nifti2Im(patientNameList(iF));
			
			% Normalize to 50
			meanSubject = xASL_stat_MeanNan(imSubject((imGM(:)+imWM(:) > 0.5)));
			imSubject = imSubject/meanSubject*50;
			
			% Increase count
			cnt(ii) = cnt(ii) + 1;
			
			% Initialize for the first subjects
			if cnt(ii) == 1
				if ii == 1
					imTemplate = imSubject;
				else
					imTemplate(:,:,:,ii) = imSubject;
				end
			else
				% Add to the template
				imTemplate(:,:,:,ii) = imTemplate(:,:,:,ii) + imSubject;
			end
		end
	end
end

% Create an average of all control+patients
imTemplate(:,:,:,length(groupName)+1) = (sum(imTemplate(:,:,:,1:length(groupName)),4))/sum(cnt);

% Create the directory for the results
resultsDirFull = fullfile(rawDir,resultsDir);
if ~exist(resultsDirFull,'dir')
	mkdir(resultsDirFull);
end

% Divide by number of files
for ii =  1:length(groupName)
	imTemplate(:,:,:,ii) = imTemplate(:,:,:,ii)/cnt(ii);
	xASL_io_SaveNifti(patientNameList{ii},fullfile(resultsDirFull,[referenceName '-' groupName{ii} '.nii']),imTemplate(:,:,:,ii));
end

xASL_io_SaveNifti(patientNameList{ii},fullfile(resultsDirFull,[referenceName '-all.nii']),imTemplate(:,:,:,length(groupName)+1));

