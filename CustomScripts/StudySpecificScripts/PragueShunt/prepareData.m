rootPath = '/pet/projekte/asl/data/PragueShunt';

% Make a list of all directories
subjectList = {'_100' '_120' '_140'};
TElist = [18.9, 21.2, 21.4];

for i = 1:length(subjectList)
	% If there's ASL and no M0 file:
	pathASL = xASL_adm_GetFsList(fullfile(rootPath,'rawnew'),['.*' subjectList{i} '\d{2}.nii'],false);
	
	for jj = 1:length(pathASL)
		imASL = xASL_io_Nifti2Im(fullfile(rootPath,'rawnew',pathASL{jj}));
		if jj == 1
			imASL4D = imASL;
		else
			imASL4D(:,:,:,jj) = imASL;
		end
	end
	
	% Save the ASL
	xASL_io_SaveNifti(fullfile(rootPath,'rawnew',pathASL{i}),fullfile(rootPath,'analysis','001',['ASL_' num2str(i)],'ASL4D.nii'),imASL4D,32,0,[]);

	parms.EchoTime = TElist(i);
	parms.RepetitionTime = 5600;
	parms.NumberOfTemporalPositions = length(subjectList);
	save(fullfile(rootPath,'analysis','001',['ASL_' num2str(i)],'ASL4D_parms.mat'),'parms');
end


