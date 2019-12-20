rootPath = '/pet/projekte/asl/data/BrnoEpilepsy';

% Make a list of all directories
allSubjects = xASL_adm_GetFsList(fullfile(rootPath,'analysis'),'^\d{4}A$',true);

for i = 1:length(allSubjects)
	% If there's ASL and no M0 file:
	pathASL = fullfile(rootPath,'analysis',allSubjects{i},'ASL_1');
	if xASL_exist(fullfile(pathASL,'ASL4D.nii')) && ~xASL_exist(fullfile(pathASL,'M0.nii'))
		% Load the ASL
		imASL = xASL_io_Nifti2Im(fullfile(pathASL,'ASL4D.nii'));
	
		% Take the 2 last and save as M0
		xASL_io_SaveNifti(fullfile(pathASL,'ASL4D.nii'),fullfile(pathASL,'M0.nii'),imASL(:,:,:,(end-1):end),32,0,[]);
		
		% Throw away the last four and save as ASL4D
		xASL_io_SaveNifti(fullfile(pathASL,'ASL4D.nii'),fullfile(pathASL,'ASL4D.nii'),imASL(:,:,:,1:(end-4)),32,0,[]);
	end
	
	if ~xASL_exist(fullfile(rootPath,'analysis',allSubjects{i},'ASL_1','ASL4D_parms.mat'))
		% Go to raw and locate the .mat file corresponding to ASL
		pathASL = xASL_adm_GetFsList(fullfile(rootPath,'raw',allSubjects{i}),'.*PCASL.*',true);
		pathMAT = xASL_adm_GetFsList(fullfile(rootPath,'raw',allSubjects{i},pathASL{1}),['.*' allSubjects{i} '.*header\.mat$'],false);
		% Load it
		mat = load(fullfile(rootPath,'raw',allSubjects{i},pathASL{1},pathMAT{1}));

		% Extract relevant parameters and rename them
		parms.EchoTime = mat.dicom_header.EchoTime;
		parms.RepetitionTime = mat.dicom_header.RepetitionTime;
		parms.NumberOfTemporalPositions = mat.dicom_header.NumberOfTemporalPositions-4;
		
		% Save the mat to the analysis directory
		save(fullfile(rootPath,'analysis',allSubjects{i},'ASL_1','ASL4D_parms.mat'),'parms');
	end
	
	if ~xASL_exist(fullfile(rootPath,'analysis',allSubjects{i},'ASL_1','M0_parms.mat'))
		% Go to raw and locate the .mat file corresponding to ASL
		pathASL = xASL_adm_GetFsList(fullfile(rootPath,'raw',allSubjects{i}),'.*PCASL.*',true);
		pathMAT = xASL_adm_GetFsList(fullfile(rootPath,'raw',allSubjects{i},pathASL{1}),['.*' allSubjects{i} '.*header\.mat$'],false);
		% Load it
		mat = load(fullfile(rootPath,'raw',allSubjects{i},pathASL{1},pathMAT{1}));

		% Extract relevant parameters and rename them
		parms.EchoTime = mat.dicom_header.EchoTime;
		parms.RepetitionTime = 8000;%mat.dicom_header.RepetitionTime;
		parms.NumberOfTemporalPositions = 2;
		
		% Save the mat to the analysis directory
		save(fullfile(rootPath,'analysis',allSubjects{i},'ASL_1','M0_parms.mat'),'parms');
	end
	
end
