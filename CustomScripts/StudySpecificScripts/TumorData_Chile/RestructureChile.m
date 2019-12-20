%InputDir = '/pet/projekte/asl/data/ChileImported/analysis';
%OutputDir = '/pet/projekte/asl/data/ChileExtra/analysis';
%bOrderOnly = 1;

%InputDir = '/pet/projekte/asl/data/ChileNewBeforeFix/analysis';
%OutputDir = '/pet/projekte/asl/data/ChileNew/analysis';
%bOrderOnly = 0;

InputDir = '/pet/projekte/asl/data/ChileNewBeforeFix2/analysis';
OutputDir = '/pet/projekte/asl/data/ChileNew2/analysis';
bOrderOnly = 1;

nSessions = 2;

bOrderOnly = 1;

% Make a list of all directories
allSubjects = xASL_adm_GetFsList(InputDir,'^.+$',true);

% Create the dir for test-files
if ~exist(fullfile(OutputDir,'test'),'dir')
	mkdir(fullfile(OutputDir,'test'));
end

for i = 1:length(allSubjects)

		% Create the subject and ASL dir
		if ~exist(fullfile(OutputDir,allSubjects{i}),'dir')
			mkdir(fullfile(OutputDir,allSubjects{i}));
		end

		for iS = 1:nSessions
			if ~exist(fullfile(OutputDir,allSubjects{i},[ 'ASL_' num2str(iS)]),'dir')
				mkdir(fullfile(OutputDir,allSubjects{i},[ 'ASL_' num2str(iS)]));
			end
		end

		% Copy the T1 files
		if exist(fullfile(InputDir, allSubjects{i}, 'T1.nii'),'file')
			system(['cp ' fullfile(InputDir, allSubjects{i},'T1.nii') ' ' fullfile(OutputDir, allSubjects{i}, 'T1.nii')]);
		else
			disp(['Missing T1 file: ' fullfile(InputDir, allSubjects{i}, 'T1.nii')]);
		end

		% Copy the Lesion files
		if exist(fullfile(InputDir, allSubjects{i}, 'Lesion_T1_1.nii.gz'),'file')
			system(['cp ' fullfile(InputDir, allSubjects{i}, 'Lesion_T1_1.nii.gz') ' ' fullfile(OutputDir, allSubjects{i}, 'Lesion_T1_1.nii.gz')]);
		else
			if strcmp(allSubjects{i}(1:3),'Sub')
				disp(['Missing lesion file: ' fullfile(InputDir, allSubjects{i}, 'Lesion_T1_1.nii.gz')]);
			end
		end

		% Copy the ASL params file
		for iS = 1:nSessions
			if exist(fullfile(InputDir,allSubjects{i}, ['ASL_' num2str(iS)],'ASL4D_parms.mat'),'file')
				system(['cp ' fullfile(InputDir, allSubjects{i}, ['ASL_' num2str(iS)],'ASL4D_parms.mat') ' ' fullfile(OutputDir, allSubjects{i},['ASL_' num2str(iS)],'ASL4D_parms.mat')]);
			else
				disp(['Missing ASL file: ' fullfile(InputDir, allSubjects{i},['ASL_' num2str(iS)],'ASL4D_parms.mat')]);
			end


			% Copy the ASL file
			if exist(fullfile(InputDir,allSubjects{i}, ['ASL_' num2str(iS)],'ASL4D.nii'),'file')
				% Load the image
				tIM     = xASL_io_Nifti2Im(fullfile(InputDir,allSubjects{i},['ASL_' num2str(iS)],'ASL4D.nii'));
				cIM = zeros(size(tIM));

				if (size(tIM,4)==1)
					% Already subtracted images are just copied
					cIM = tIM;
					aslIM = cIM;
				else
					if bOrderOnly
						timeHalf = size(tIM,4)/2;
						cIM(:,:,:,1:2:end) = tIM(:,:,:,1:timeHalf);
						cIM(:,:,:,2:2:end) = tIM(:,:,:,(timeHalf+1):end);
						aslIM = mean(cIM(:,:,:,1:2:end)-cIM(:,:,:,2:2:end),4);
					else
						if mod(size(tIM,3),2)
							sliceHalf = floor(size(tIM,3)/2);
							timeHalf = size(tIM,4)/2;
							cIM(:,:,1:2:end,1:2:end) = tIM(:,:,1:(sliceHalf+1),1:timeHalf);
							cIM(:,:,2:2:end,2:2:end) = tIM(:,:,(sliceHalf+2):end,1:timeHalf);
							
							cIM(:,:,2:2:end,1:2:end) = tIM(:,:,1:sliceHalf,(timeHalf+1):end);
							cIM(:,:,1:2:end,2:2:end) = tIM(:,:,(sliceHalf+1):end,(timeHalf+1):end);
							
							aslIM = mean(cIM(:,:,:,1:2:end)-cIM(:,:,:,2:2:end),4);
						else
							sliceHalf = size(tIM,3)/2;
							timeHalf = size(tIM,4)/2;
							cIM(:,:,1:2:end,1:2:end) = tIM(:,:,1:sliceHalf,1:timeHalf);
							cIM(:,:,1:2:end,2:2:end) = tIM(:,:,(sliceHalf+1):end,1:timeHalf);
							
							cIM(:,:,2:2:end,1:2:end) = tIM(:,:,1:sliceHalf,(timeHalf+1):end);
							cIM(:,:,2:2:end,2:2:end) = tIM(:,:,(sliceHalf+1):end,(timeHalf+1):end);
							
							aslIM = mean(cIM(:,:,:,1:2:end)-cIM(:,:,:,2:2:end),4);
						end
					end
				end
				xASL_io_SaveNifti(fullfile(InputDir,allSubjects{i},['ASL_' num2str(iS)],'ASL4D.nii'), fullfile(OutputDir,allSubjects{i},['ASL_' num2str(iS)],'ASL4D.nii'), cIM);%, BITS, bGZIP )
				xASL_io_SaveNifti(fullfile(InputDir,allSubjects{i},['ASL_' num2str(iS)],'ASL4D.nii'), fullfile(OutputDir,'test',[allSubjects{i} 'ASL4D_' num2str(iS) '.nii']),aslIM);%, BITS, bGZIP )
			else
				disp(['Missing ASL file: ' fullfile(InputDir, allSubjects{i},['ASL_' num2str(iS)],'ASL4D.nii')]);
			end
		end

end


%size(tIM)
%clear cIM
%cIM(:,:,1:2:15,[1:2:60-1])     = tIM(:,:,1:8,1:30);
%cIM(:,:,2:2:16,[1:2:60-1])     = tIM(:,:,1:8,31:60);

cIM(:,:,1:2:15,[2:2:60-0])     = tIM(:,:,9:16,1:30);
cIM(:,:,2:2:16,[2:2:60-0])     = tIM(:,:,9:16,31:60);
