%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select the output type
regType = 'SPMAFF';
%regType = 'SPM5';
%regType = 'DARTEL';
%regType = 'DARTELxasl';
%regType = 'GS';
%regType = 'DARTELmixed';
%regType = 'GSmixed';
%regType = 'DARTEL2';
%regType = 'DARTEL4';

%regTypeVec = {'SPMAFF','SPM5','DARTEL','DARTELmixed','DARTELxasl','DARTEL2','DARTEL4','GS','GS2','GS4','GSmixed'};
regTypeVec = {'GS','GS2','GS4','GSmixed'};%%;
bSaveGM = 0;

for ir = 1:length(regTypeVec)
	regType = regTypeVec{ir};
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Configure everything according to the output type
	pathRoot = '/pet/projekte/asl/data/gliomaMNI/';
	segStr = '_DEFmni';

	switch (regType)
		case 'SPMAFF'
			inPath = fullfile(pathRoot,'analysisDARTEL');
			coregPath = fullfile(pathRoot,'template-coreg/DARTEL/');
			segStr = '_LINmni';
		case 'SPM5'
			inPath = fullfile(pathRoot,'analysisDARTEL');
			coregPath = fullfile(pathRoot,'template-coreg/DARTEL/');
			deformName = 'y_T1_SPM.nii';
		case 'DARTEL'
			inPath = fullfile(pathRoot,'analysisDARTEL');
			coregPath = fullfile(pathRoot,'template-coreg/DARTEL/');
			deformName = 'y_T1.nii';
		case 'DARTEL2'
			inPath = fullfile(pathRoot,'analysisDARTEL');
			coregPath = fullfile(pathRoot,'template-coreg/DARTEL/');
			deformName = 'y_T1_2.nii';
		case 'DARTEL4'
			inPath = fullfile(pathRoot,'analysisDARTEL');
			coregPath = fullfile(pathRoot,'template-coreg/DARTEL/');
			deformName = 'y_T1_4.nii';
		case 'GS'
			inPath = fullfile(pathRoot,'analysisGS');
			coregPath = fullfile(pathRoot,'template-coreg/shooting/');
			deformName = 'y_T1_orig.nii';
		case 'GS2'
			inPath = fullfile(pathRoot,'analysisGS');
			coregPath = fullfile(pathRoot,'template-coreg/shooting/');
			deformName = 'y_T1_2.nii';
		case 'GS4'
			inPath = fullfile(pathRoot,'analysisGS');
			coregPath = fullfile(pathRoot,'template-coreg/shooting/');
			deformName = 'y_T1_4.nii';
		case 'DARTELmixed'
			inPath = fullfile(pathRoot,'analysisDARTEL');
			coregPath = fullfile(pathRoot,'template-coreg/DARTEL/');
			deformName = 'y_T1_mixed.nii';
		case 'DARTELxasl'
			inPath = fullfile(pathRoot,'analysisDARTEL');
			coregPath = fullfile(pathRoot,'template-coreg/DARTEL/');
			deformName = 'y_T1.nii';
		case 'GSmixed'
			inPath = fullfile(pathRoot,'analysisGS');
			coregPath = fullfile(pathRoot,'template-coreg/shooting/');
			deformName = 'y_T1_mixed.nii';
	end

	outPath  = fullfile(pathRoot,['resData' regType]);

	if ~exist(outPath,'dir')
		mkdir(outPath);
	end

	load(fullfile(inPath,'list.mat'));

	for ii = 1:2
		% Run it for GBM and NON-GBM
		if ii == 1
			outPathList = fullfile(outPath,'GBM');
			outGBMExt = 'enh';
		else
			outPathList = fullfile(outPath,'NONGBM');
			outGBMExt = 'nonenh';
		end

		% Read all the patients
		for jj = 1:length(listNames(ii,:))

			% Session
			for kk = 1:2
				%if ( (kk==2) && (ii==1) && (jj==7))

					if kk == 1
						outSessExt = 'pre';
					else
						outSessExt = 'post';
					end

					% Save the name for re-creation
					patName = listNames{ii,jj};

					inPathPat = fullfile(inPath,[listNames{ii,jj} '_0' num2str(kk)]);
					outPathPat = fullfile(outPathList,[listNames{ii,jj}]);

					% Create the pre and post directories
					if ~exist(outPathPat,'dir')
						mkdir(outPathPat);
					end

					if strcmp(regType,'SPMAFF')
						matAffine = load(fullfile(inPath,'Population','TissueVolume',['cat_T1_', listNames{ii,jj} '_0' num2str(kk), '.mat']));
					end

					% Transform the Lesions, GM and T1
					matlabbatch = [];
					switch (regType)
						case 'SPMAFF'
							matlabbatch{1}.spm.util.defs.comp{1}.id.space = {[coregPath 'Template_6_IXI555_MNI152.nii']};
						case {'SPM5','DARTEL','DARTELxasl','DARTEL2','DARTEL4','GS','GS2','GS4','DARTELmixed','GSmixed'}
							xASL_adm_UnzipNifti(fullfile(inPathPat,deformName));
							matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(inPathPat,deformName)};
					end
					switch (regType)
						case {'SPMAFF', 'SPM5','DARTEL','DARTEL2','DARTEL4','DARTELxasl','DARTELmixed'}
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {[coregPath 'u_wmni_icbm152_gm_tal_nlin_sym_09a.nii']};
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [0 1];
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {[coregPath 'Template_6_IXI555_MNI152.nii']};
						case {'GS','GS2','GS4','GSmixed'}
							matlabbatch{1}.spm.util.defs.comp{2}.comp{1}.inv.comp{1}.def = {[coregPath 'y_wmni_icbm152_gm_tal_nlin_sym_09a_Template.nii']};
							matlabbatch{1}.spm.util.defs.comp{2}.comp{1}.inv.space = {[coregPath 'wTemplate_T1_IXI555_MNI152_GS_sym.nii']};
					end
					matlabbatch{1}.spm.util.defs.comp{3}.id.space = {[coregPath 'mni_icbm152_gm_tal_nlin_sym_09a.nii']};
					matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {
						fullfile(inPathPat,'pT1.nii')
						fullfile(inPathPat,'pc1T1.nii')
						};
					matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
					matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
					matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
					matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
					matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';

					% The index also correctly gives the lesion number. So  listTMRCAV{ii,jj,kk,2} == 1 if CAV exists and TMR not, and == 2 if both CAV and TMR exist
					if (listTMRCAV{ii,jj,kk,1} + listTMRCAV{ii,jj,kk,2}) > 0
						matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{end+1} = fullfile(inPathPat,'pLesion_T1_1.nii');
					end

					if (listTMRCAV{ii,jj,kk,1} + listTMRCAV{ii,jj,kk,2}) > 1
						matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{end+1} = fullfile(inPathPat,'pLesion_T1_2.nii');
					end

					% We need to make a copy with a 'p' prefix
					ptLen = length(inPathPat);
					for ll = 1:length(matlabbatch{1}.spm.util.defs.out{1}.pull.fnames)
						xASL_adm_UnzipNifti(fullfile(inPathPat,matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{ll}((ptLen+3):end)));
						switch(regType)
							case 'SPMAFF'
								tmpVol = spm_vol(fullfile(inPathPat,matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{ll}((ptLen+3):end)));
								tmpIm  = spm_read_vols(tmpVol);
								tmpVol.fname = matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{ll};
								tmpVol.mat = matAffine.S.parameter.spm.Affine*tmpVol.mat;
								tmpVol.private.mat = tmpVol.mat;
								tmpVol.private.mat0 = tmpVol.mat;
								spm_write_vol(tmpVol,tmpIm);

							case {'SPM5','DARTEL','DARTEL2','DARTEL4','GS','GS2','GS4','DARTELxasl','DARTELmixed','GSmixed'}
								system(['cp ' fullfile(inPathPat,matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{ll}((ptLen+3):end)) ' ' matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{ll}]);
						end
					end

					spm_jobman('run',matlabbatch);

					% And then delete this copy
					for ll = 1:length(matlabbatch{1}.spm.util.defs.out{1}.pull.fnames)
						system(['rm ' matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{ll}]);
					end

					% Copy the T1 files results
					system(['mv ' fullfile(inPathPat,'wpT1.nii') ' ' fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_T1c.nii'])]);

					Fname = fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_GM.nii']);
					system(['mv ' fullfile(inPathPat,'wpc1T1.nii') ' ' Fname]);
					% If keeping GM, then zip it and delete the origin NII, otherwise, just delete the original without prior zipping
					if bSaveGM
						gzip(Fname);
					end
					delete(Fname);

					% Check and copy TMR and CAV
					if ~exist(fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_PWIseg']),'dir')
						mkdir(fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_PWIseg']));
					end

					if listTMRCAV{ii,jj,kk,1}
						Fname = fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_PWIseg'],[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_TMR.nii']);
						system(['mv ' fullfile(inPathPat,'wpLesion_T1_1.nii') ' ' Fname]);
						gzip(Fname);
						delete(Fname);
					end

					if listTMRCAV{ii,jj,kk,2}
						Fname = fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_PWIseg'],[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_CAV.nii']);
						system(['mv ' fullfile(inPathPat,['wpLesion_T1_' num2str(listTMRCAV{ii,jj,kk,2}) '.nii']) ' ' Fname]);
						gzip(Fname);
						delete(Fname);
					end

					% Transforms the Chollet landmarks - deprecated direct transformation
					% 			matlabbatch = [];
					% 			matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(inPathPat,deformName)};
					% 			switch (regType)
					% 				case {'SPMAFF', 'SPM5','DARTEL'}
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {[coregPath 'u_wmni_icbm152_gm_tal_nlin_sym_09a.nii']};
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [0 1];
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {[coregPath 'Template_6_IXI555_MNI152.nii']};
					% 				case 'GS'
					% 					matlabbatch{1}.spm.util.defs.comp{2}.comp{1}.inv.comp{1}.def = {[coregPath 'y_wmni_icbm152_gm_tal_nlin_sym_09a_Template.nii']};
					% 					matlabbatch{1}.spm.util.defs.comp{2}.comp{1}.inv.space = {[coregPath 'wTemplate_T1_IXI555_MNI152_GS_sym.nii']};
					% 			end
					% 			matlabbatch{1}.spm.util.defs.comp{3}.id.space = {[coregPath 'mni_icbm152_gm_tal_nlin_sym_09a.nii']};
					% 			matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(inPathPat,'Chol.nii')};
					% 			matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
					% 			matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
					% 			matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
					% 			matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
					% 			matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
					% 			spm_jobman('run',matlabbatch);
					%
					% 			system(['mv ' fullfile(inPathPat,'wChol.nii') ' ' fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_PWIseg'],[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt '_chollet.nii'])]);


					% First calculates the joint transformation of everything - deprecated
					% 			matlabbatch = [];
					% 			matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(inPathPat,deformName)};
					% 			switch (regType)
					% 				case {'SPMAFF', 'SPM5','DARTEL'}
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {[coregPath 'u_wmni_icbm152_gm_tal_nlin_sym_09a.nii']};
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [0 1];
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
					% 					matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {[coregPath 'Template_6_IXI555_MNI152.nii']};
					% 				case 'GS'
					% 					matlabbatch{1}.spm.util.defs.comp{2}.comp{1}.inv.comp{1}.def = {[coregPath 'y_wmni_icbm152_gm_tal_nlin_sym_09a_Template.nii']};
					% 					matlabbatch{1}.spm.util.defs.comp{2}.comp{1}.inv.space = {[coregPath 'wTemplate_T1_IXI555_MNI152_GS_sym.nii']};
					% 			end
					% 			matlabbatch{1}.spm.util.defs.comp{3}.id.space = {[coregPath 'mni_icbm152_gm_tal_nlin_sym_09a.nii']};
					% 			matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'joint.nii';
					% 			matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {inPathPat};
					% 			spm_jobman('run',matlabbatch);

					% Transform Chollet landmarks - the coordinates version
					% Then inverts the joint transformation so that we know for each entry point where it is transformed
					matlabbatch = [];
					matlabbatch{1}.spm.util.defs.comp{1}.id.space = {[coregPath 'wmni_icbm152_gm_tal_nlin_sym_09a.nii']};
					switch (regType)
						case {'SPMAFF', 'SPM5','DARTEL','DARTELxasl','DARTEL2','DARTELxasl','DARTEL4','DARTELmixed'}
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {[coregPath 'u_wmni_icbm152_gm_tal_nlin_sym_09a.nii']};
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
							matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {[coregPath 'Template_6_IXI555_MNI152.nii']};
						case {'GS','GS2','GS4','GSmixed'}
							matlabbatch{1}.spm.util.defs.comp{2}.def = {[coregPath 'y_wmni_icbm152_gm_tal_nlin_sym_09a_Template.nii']};
					end
					switch (regType)
						case 'SPMAFF'
							tmpVol = spm_vol(fullfile(inPathPat,'T1.nii'));
							tmpIm  = spm_read_vols(tmpVol);
							tmpVol.fname = fullfile(inPathPat,'pT1.nii');
							tmpVol.mat = matAffine.S.parameter.spm.Affine*tmpVol.mat;
							tmpVol.private.mat = tmpVol.mat;
							tmpVol.private.mat0 = tmpVol.mat;
							spm_write_vol(tmpVol,tmpIm);
							matlabbatch{1}.spm.util.defs.comp{3}.id.space = {fullfile(inPathPat,'pT1.nii')};
						case {'SPM5','DARTEL','DARTEL2','DARTELxasl','DARTEL4','GS','GS2','GS4','DARTELmixed','GSmixed'}
							matlabbatch{1}.spm.util.defs.comp{3}.comp{1}.inv.comp{1}.def = {fullfile(inPathPat,deformName)};
							matlabbatch{1}.spm.util.defs.comp{3}.comp{1}.inv.space = {fullfile(inPathPat,'T1.nii')};
					end
					matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'inversejoint.nii';
					matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {inPathPat};
					spm_jobman('run',matlabbatch);

					% And then delete this copy
					if strcmp(regType,'SPMAFF')
						system(['rm ' fullfile(inPathPat,'pT1.nii')]);
					end

					% Loads the image with landmarks
					imChollet = xASL_io_Nifti2Im(fullfile(inPathPat,'Chol.nii'));

					% Loads the transformation file
					imInverseJoint = xASL_io_Nifti2Im(fullfile(inPathPat,'y_inversejoint.nii'));
					delete(fullfile(inPathPat,'y_inversejoint.nii'));
					% Loads vol of the MNI output file
					volCholletMNI = spm_vol(fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_T1c.nii']));

					% Goes through all landmarks
					nChollet = max(imChollet(:));
					imCholletMNI = zeros(volCholletMNI.dim);
					for iL = 1 : nChollet
						% Finds all voxels of the landmark and calculate the center/mean coordinates
						X = find(imChollet(:) == iL);
						if ~isempty(X)
							[X,Y,Z] = ind2sub(size(imChollet),X);

							mX = round(mean(X));
							mY = round(mean(Y));
							mZ = round(mean(Z));

							% Transforms the coordinates to MNI world coordinates
							%imChollet(mX,mY,mZ)
							xMNI = squeeze(imInverseJoint(mX,mY,mZ,1,1:3));

							% Transforms the coordinates to MNI voxel coordinates
							cMNI = ((volCholletMNI.mat)^-1)*[xMNI; 1];

							% NIFTI indexes the voxel from 0 to dim-1 (center of the voxels)
							% So we need to add 1, as Matlab indexes from 1
							cMNI = round(cMNI)+1;

							% Put the voxel to the images + creates a 3x3x3 cube
							imCholletMNI((-1:1)+cMNI(1),(-1:1)+cMNI(2),(-1:1)+cMNI(3)) = iL;
						end
					end
					% Save the transformed Chollet marks image
					xASL_io_SaveNifti(fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_T1c.nii']),fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_PWIseg'],[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_POINTS.nii']),imCholletMNI,16,1);

					% Copy the Chollet landmarks
					%system(['mv ' fullfile(inPathPat,'wROI_T1_1.nii') ' ' fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_PWIseg'],[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt '_chollet.nii'])]);
					%system(['cp ' fullfile(inPathList,listNames{ii,jj},[listNames{ii,jj} '_' inPathExt '_PWIseg'],[listNames{ii,jj} '_' inPathExt '_post_chollet.nii.gz']) ' ' fullfile(outPath,[listNames{ii,jj} '_02'],'ROI_T1_1.nii.gz')]);

					Fname = fullfile(outPathPat,[listNames{ii,jj} '_' outGBMExt '_r' regType '_' outSessExt segStr '_T1c.nii']);
					gzip(Fname);
					delete(Fname);
				end
			end
%		end
	end
end
