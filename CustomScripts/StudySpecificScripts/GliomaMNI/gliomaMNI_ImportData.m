% Read the data as given by Martin and crunch them to xASL structure
inPath  = '/pet/projekte/asl/data/gliomaMNI/origData';
outPath = '/pet/projekte/asl/data/gliomaMNI/analysis';

if ~exist(outPath,'dir')
	mkdir(outPath);
end

% Run it for GBM and NON-GBM
for ii = 1:2
	if ii == 1
		inPathList = fullfile(inPath,'GBM');
		inPathExt = 'enh';
	else
		inPathList = fullfile(inPath,'NONGBM');
		inPathExt = 'nonenh';
	end
	
	fileList = dir(inPathList);
	
	% Read all the patients
	jj = 0;
	for kk = 1:length(fileList)
		if (fileList(kk).isdir) && (fileList(kk).name(1) ~= '.')
			jj = jj+1;
			% Save the name for re-creation
			listNames{ii,jj} = fileList(kk).name;
			
			% Create the pre and post directories
			if ~exist(fullfile(outPath,[listNames{ii,jj} '_01']),'dir')
				mkdir(fullfile(outPath,[listNames{ii,jj} '_01']));
			end
			
			if ~exist(fullfile(outPath,[listNames{ii,jj} '_02']),'dir')
				mkdir(fullfile(outPath,[listNames{ii,jj} '_02']));
			end
			
			% Copy the T1 files
			system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_pre_T1c.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_01'],'T1.nii.gz')]);
			system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_post_T1c.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_02'],'T1.nii.gz')]);
			
			% Copy the Chollet landmarks
			system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_pre_chollet.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_01'],'ROI_T1_1.nii.gz')]);
			system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_post_chollet.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_02'],'ROI_T1_1.nii.gz')]);
			
			% Read all these files so that they would unzip
			xASL_adm_UnzipNifti(fullfile(outPath,[listNames{ii,jj} '_01'],'T1.nii.gz'));
			xASL_adm_UnzipNifti(fullfile(outPath,[listNames{ii,jj} '_02'],'T1.nii.gz'));
			xASL_adm_UnzipNifti(fullfile(outPath,[listNames{ii,jj} '_01'],'ROI_T1_1.nii.gz'));
			xASL_adm_UnzipNifti(fullfile(outPath,[listNames{ii,jj} '_02'],'ROI_T1_1.nii.gz'));
			
			% Check and copy TMR and CAV
			if exist(fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_pre_TMR.nii.gz']),'file')
				system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_pre_TMR.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_01'],'Lesion_T1_1.nii.gz')]);
				indx = 2;
				listTMRCAV{ii,jj,1,1} = 1;
			else
				indx = 1;
				listTMRCAV{ii,jj,1,1} = 0;
			end
			
			if exist(fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_pre_CAV.nii.gz']),'file')
				system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_pre_CAV.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_01'],['Lesion_T1_' num2str(indx) '.nii.gz'])]);
				listTMRCAV{ii,jj,1,2} = indx;
			else
				listTMRCAV{ii,jj,1,2} = 0;
			end
			
			if exist(fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_post_TMR.nii.gz']),'file')
				system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_post_TMR.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_02'],'Lesion_T1_1.nii.gz')]);
				indx = 2;
				listTMRCAV{ii,jj,2,1} = 1;
			else
				indx = 1;
				listTMRCAV{ii,jj,2,1} = 0;
			end
			
			if exist(fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_post_CAV.nii.gz']),'file')
				system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_post_CAV.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_02'],['Lesion_T1_' num2str(indx) '.nii.gz'])]);
				listTMRCAV{ii,jj,2,2} = indx;
			else
				listTMRCAV{ii,jj,2,2} = 0;
			end
			
			% Session
			for li = 1:2
				% Reslices it to the T1 space
				matlabbatch = [];
				matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],'T1.nii')};
				matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],['ROI_T1_1.nii'])};
				matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
				matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
				matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
				matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
				matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
				spm_jobman('run',matlabbatch);
						
				% Delete the original and replace it with the resliced file
				system(['rm ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],'ROI_T1_1.nii')]);
				system(['mv ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],'wROI_T1_1.nii') ' ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],'Chol.nii')]);
				%system(['mv ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],'wROI_T1_1.nii') ' ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],'ROI_T1_1.nii')]);
				
				% Lesion nr
				for lj = 1:2
					% If exists
					if listTMRCAV{ii,jj,li,lj}
						% Read the file so it gets unzipped
						
						xASL_adm_UnzipNifti(fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],['Lesion_T1_' num2str(listTMRCAV{ii,jj,li,lj}) '.nii']));
						
						% Reslices it to the T1 space
						matlabbatch = [];
						matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],'T1.nii')};
						matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],['Lesion_T1_' num2str(listTMRCAV{ii,jj,li,lj}) '.nii'])};
						matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
						matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
						matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
						matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
						matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
						spm_jobman('run',matlabbatch);
						
						% Delete the original and replace it with the resliced file
						system(['rm ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],['Lesion_T1_' num2str(listTMRCAV{ii,jj,li,lj}) '.nii'])]);
						system(['mv ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],['wLesion_T1_' num2str(listTMRCAV{ii,jj,li,lj}) '.nii']) ' ' fullfile(outPath,[listNames{ii,jj} '_0' num2str(li)],['Lesion_T1_' num2str(listTMRCAV{ii,jj,li,lj}) '.nii'])]);
					end
				end
			end
		end
	end
end

save(fullfile(outPath,'list.mat'),'listNames','listTMRCAV');