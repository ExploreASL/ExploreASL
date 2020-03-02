%% Initialize ExploreASL
% First import data with ImportStep
% Then run the normal xASL analysis
% Then run the preprocess step to do the additional resampling and co-registration - something to be later integrated directly in xASL
% Now run the analysis part
%% Initialize the paths
rawDir    = '/pet/projekte/asl/data/FRONTIER';
modDir = 'ASL_1';
%modDir = 'PET_1';
gmth = 0.7;

%% Do all the comparisons

% Load all the ROIs, Maps, and CBFs
patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P\d{2}$', 'List', [], 1);

tmpIm = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{1},modDir,'CBF.nii'));

% Prepare the empty matrices
imGM = zeros([size(tmpIm) length(patientNameList)]);
imWM = zeros(size(imGM));

imPET = zeros([size(tmpIm) 2 length(patientNameList)]);
imASL = zeros(size(imPET));
imDSC = zeros(size(imPET));

imFLAIR = zeros([size(tmpIm) 6 length(patientNameList)]);
imT1 = zeros(size(imFLAIR));

vecASL = ones(length(patientNameList),1);
vecPET = zeros(length(patientNameList),1);
vecDSC = zeros(length(patientNameList),1);
vecFLAIR = zeros(length(patientNameList),1);
vecT1  = zeros(length(patientNameList),1);


for iL = 1:length(patientNameList)
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_CBF.nii'))
		imTmp = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_GM.nii'));imGM(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),iL) = imTmp;
		imWM(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_WM.nii'));

		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_PET.nii'))
			vecPET(iL) = 1;
			imPET(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),1,iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_PET.nii'));
			imPET(:,:,:,2,iL) = imPET(:,:,:,1,iL);
		end

		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_Lesion_FLAIR.nii'))
			vecFLAIR(iL) = 1;
			imFLAIR(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),:,iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_Lesion_FLAIR.nii'));
		end

		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_Lesion_T1.nii'))
			vecT1(iL) = 1;
			imT1(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),:,iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_Lesion_T1.nii'));
		end

		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_DSC.nii'))
			vecDSC(iL) = 1;
			imDSC(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),1,iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_DSC.nii'));
			imDSC(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),2,iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_DSC_Deform.nii'));
		end

		imASL(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),1,iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_CBF.nii'));
		imASL(1:size(imTmp,1),1:size(imTmp,2),1:size(imTmp,3),2,iL) = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_CBF_Deform.nii'));
	end
end

% Take FLAIR ROI for the patients 3,5,6
vecT1([3,5,6]) = vecFLAIR([3,5,6]);
imT1(:,:,:,:,[3,5,6]) = imFLAIR(:,:,:,:,[3,5,6]);

resVec = zeros(3,3,4,8);
resMean = zeros(3,3,4,8,2);
resMeanGM = zeros(3,3,4,8,2);
resMax = zeros(3,3,4,8,2);
resHist = zeros(3,3,4,8,60,60);
% Compare PET, ASL, DSC
for iMod = 1:3
	switch (iMod)
		case 1
			imRef = imPET;
			imSrc = imASL;
		case 2
			imRef = imPET;
			imSrc = imDSC;
		case 3
			imRef = imDSC;
			imSrc = imASL;
	end

	% Compare CBF, CBFnorm, CBFdefnorm
	for iVal = 1:3

		% Compare whole brain normal, contralateral normal, T1-ROI, FLAIR-ROI
		for iRoi = 1:4
			% For all patients
			for iL = 1:length(patientNameList)
				% See of the modality exists
				switch(iMod)
					case 1
						resVec(iMod,iVal,iRoi,iL) = vecPET(iL);
					case 2
						resVec(iMod,iVal,iRoi,iL) = vecPET(iL)*vecDSC(iL);
					case 3
						resVec(iMod,iVal,iRoi,iL) = vecDSC(iL);
				end
				% See if the ROI is present
				if iRoi == 3
					resVec(iMod,iVal,iRoi,iL) = resVec(iMod,iVal,iRoi,iL)*vecT1(iL);
				end

				% Select normalization and deformation
				switch(iVal)
					case 1
						imLocRef = imRef(:,:,:,1,iL);
						imLocSrc = imSrc(:,:,:,1,iL);
					case 2
						% Obtain the GM>70% in contralateral ROI
						imLocMask = (imFLAIR(:,:,:,6,iL).*(imGM(:,:,:,iL)>gmth).*(imRef(:,:,:,1,iL)>0).*(imSrc(:,:,:,1,iL)>0))>0;
						% Normalize by the healthy cortex value
						imLocRef = imRef(:,:,:,1,iL);
						imLocRef = imLocRef./(mean(imLocRef(imLocMask)));
						imLocSrc = imSrc(:,:,:,1,iL);
						imLocSrc = imLocSrc./(mean(imLocSrc(imLocMask)));
					case 3
						% Normalized+deform
						% Obtain the GM>70% in contralateral ROI
						imLocMask = (imFLAIR(:,:,:,6,iL).*(imGM(:,:,:,iL)>gmth).*(imRef(:,:,:,2,iL)>0).*(imSrc(:,:,:,2,iL)>0))>0;
						% Normalize by the healthy cortex value
						imLocRef = imRef(:,:,:,2,iL);
						imLocRef = imLocRef./(mean(imLocRef(imLocMask)));
						imLocSrc = imSrc(:,:,:,2,iL);
						imLocSrc = imLocSrc./(mean(imLocSrc(imLocMask)));
				end

				% Select the correct ROI
				switch(iRoi)
					case 1
						imROI = (imFLAIR(:,:,:,3,iL) + imFLAIR(:,:,:,6,iL) + imFLAIR(:,:,:,5,iL))>0;
					case 2
						imROI = imFLAIR(:,:,:,6,iL)>0;
					case 3
						imROI = imT1(:,:,:,1,iL)>0;
					case 4
						imROI = imFLAIR(:,:,:,1,iL)>0;
				end
				% Calculate the values, or do the histograms

				imROI = imROI.*(imLocRef>0).*(imLocSrc>0);
				if iVal == 1
					imROI = imROI.*(imLocRef<350).*(imLocSrc<350);
				else
					imROI = imROI.*(imLocRef<3.8).*(imLocSrc<3.8);
				end
				imROI = imROI>0;

				% Normal mean
				resMean(iMod,iVal,iRoi,iL,1) = mean(imLocRef(imROI));
				resMean(iMod,iVal,iRoi,iL,2) = mean(imLocSrc(imROI));

				% Mean over GM
				resMeanGM(iMod,iVal,iRoi,iL,1) = mean(imLocRef((imROI.*(imGM(:,:,:,iL)>gmth))>0));
				resMeanGM(iMod,iVal,iRoi,iL,2) = mean(imLocSrc((imROI.*(imGM(:,:,:,iL)>gmth))>0));
				
				% ASL vs PET
				if iMod == 1 && (iVal == 1 || iVal == 2) && (iRoi == 2 || iRoi == 3)
					listValRef{iMod,iVal,iRoi,iL} = imLocRef(imROI);
					listValSrc{iMod,iVal,iRoi,iL} = imLocSrc(imROI);
				end
				
				% DSC vs PET
				if iMod == 2 && iVal == 2 && (iRoi == 2 || iRoi == 3)
					listValRef{iMod,iVal,iRoi,iL} = imLocRef(imROI);
					listValSrc{iMod,iVal,iRoi,iL} = imLocSrc(imROI);
				end

				% 95% percentile
				imTmp1 = imLocRef(imROI);
				if ~isempty(imTmp1)
					imTmp = sort(imTmp1);
					resMax(iMod,iVal,iRoi,iL,1) = imTmp(floor(0.95*length(imTmp)));
					imTmp2 = imLocSrc(imROI);
					imTmp = sort(imTmp2);
					resMax(iMod,iVal,iRoi,iL,2) = imTmp(floor(0.95*length(imTmp)));
					if iVal == 1
						resHist(iMod,iVal,iRoi,iL,:,:) = xASL_im_JointHist(imTmp1,imTmp2,ones(size(imTmp1)),0,150,0,150,60);
					else
						resHist(iMod,iVal,iRoi,iL,:,:) = xASL_im_JointHist(imTmp1,imTmp2,ones(size(imTmp1)),0,3,0,3,60);
					end
				end


			end
		end

	end

end

%% Create graphs
% Compare PET, ASL, DSC
nMod = 2;
for iMod = 1:nMod
	switch (iMod)
		case 1
			strMod = 'PET-ASL';
		case 2
			strMod = 'PET-DSC';
		case 3
			strMod = 'DSC-ASL';
	end
	
	% Compare whole brain normal, contralateral normal, T1-ROI, FLAIR-ROI
	for iRoi = 1:4
		switch(iRoi)
			case 1
				strRoi = 'Whole brain';
			case 2
				strRoi = 'Contralateral';
			case 3
				strRoi = 'T1-ROI';
			case 4
				strRoi = 'FLAIR-ROI';
		end
		% Scatter plots of the comparisons for mean, meanGM, and max
		
		figure(1);subplot(nMod,4,4*(iMod-1)+iRoi);ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		%groupMean = mean(resMean(iMod,1,iRoi,ind,1))/mean(resMean(iMod,2,iRoi,ind,1));
		%plot(squeeze(resMean(iMod,1,iRoi,ind,1))/groupMean,squeeze(resMean(iMod,1,iRoi,ind,2))/groupMean,'r+');hold on
		plot(squeeze(resMean(iMod,2,iRoi,ind,1)),squeeze(resMean(iMod,2,iRoi,ind,2)),'go');hold on
		plot(squeeze(resMean(iMod,3,iRoi,ind,1)),squeeze(resMean(iMod,3,iRoi,ind,2)),'bx');
		plot([0.4,3.2],[0.4,3.2],'k--');
		title(['mean CBF ' strMod ' ' strRoi]);
		
		figure(2);subplot(nMod,4,4*(iMod-1)+iRoi);ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		plot(squeeze(resMeanGM(iMod,1,iRoi,ind,1)),squeeze(resMeanGM(iMod,1,iRoi,ind,2)),'r+');hold on
		%plot(squeeze(resMeanGM(iMod,2,iRoi,ind,1)),squeeze(resMeanGM(iMod,2,iRoi,ind,2)),'go');
		%plot(squeeze(resMeanGM(iMod,3,iRoi,ind,1)),squeeze(resMeanGM(iMod,3,iRoi,ind,2)),'bx');
		plot([30,80],[30,80],'k--');
		title(['mean GM-CBF ' strMod ' ' strRoi]);
		
		if ((iMod == 1)||(iMod == 2)) && (iRoi == 2 || iRoi == 3)
			% iVal 1 non-normalized
			% iVal 2 normalized
			% iMod 2 PET vs DSC
			if iMod == 1
				iValN = 1:2;
			else
				iValN = 2;
			end
			
			for iVal = iValN
				ind = find(squeeze(resVec(iMod,iVal,iRoi,:)));
				if (iVal == 1) || (iVal == 2 && iRoi == 3)
					[b,CI,pval,stats] = xASL_stat_MultipleLinReg(squeeze(resMeanGM(iMod,iVal,iRoi,ind,1)),squeeze(resMeanGM(iMod,iVal,iRoi,ind,2)),true);
					[bnoInt,CInoInt,pvalnoInt,statsnoInt] = xASL_stat_MultipleLinReg(squeeze(resMeanGM(iMod,iVal,iRoi,ind,1)),squeeze(resMeanGM(iMod,iVal,iRoi,ind,2)),false);
					% Mean relative difference
					mrd = mean(2*abs(squeeze(resMeanGM(iMod,iVal,iRoi,ind,1))-squeeze(resMeanGM(iMod,iVal,iRoi,ind,2)))./(squeeze(resMeanGM(iMod,iVal,iRoi,ind,1))+squeeze(resMeanGM(iMod,iVal,iRoi,ind,2))))*100;
				end
				if iRoi == 3
					[bmax,CImax,pvalmax,statsmax] = xASL_stat_MultipleLinReg(squeeze(resMax(iMod,iVal,iRoi,ind,1)),squeeze(resMax(iMod,iVal,iRoi,ind,2)),true);
					mrdmax = mean(2*abs(squeeze(resMax(iMod,iVal,iRoi,ind,1))-squeeze(resMax(iMod,iVal,iRoi,ind,2)))./(squeeze(resMax(iMod,iVal,iRoi,ind,1))+squeeze(resMax(iMod,iVal,iRoi,ind,2))))*100;
				end
				xx = [];
				yy = [];
				for iL = 1:size(listValRef,4)
					xx = [xx;squeeze(listValRef{iMod,iVal,iRoi,iL})];
					yy = [yy;squeeze(listValSrc{iMod,iVal,iRoi,iL})];
				end
				[bAll,CIAll,pvalAll,statsAll] = xASL_stat_MultipleLinReg(xx,yy,true);
				[bAllnoint,CIAllnoint,pvalAllnoint,statsAllnoint] = xASL_stat_MultipleLinReg(xx,yy,false);
				if iVal == 1
					indxx = (xx>3).*(yy>3);
				else
					indxx = (xx>0.06).*(yy>0.06);
				end
				xxnonlow = xx(indxx>0);
				yynonlow = yy(indxx>0);
				
				mrdvoxel = mean(2*abs(xxnonlow-yynonlow)./(xxnonlow+yynonlow))*100;
				
				str = 'Comparison in';
				if iRoi == 2
					str = [str ' contralateral hemisphere'];
					figure(15);
				else
					str = [str ' tumor region'];
					figure(16);
				end
				
				if iMod == 1
					if iVal == 1
						str = [str ' for PET vs ASL:\n'];
						subplot(1,3,1);
					else
						str = [str ' for normalized PET vs ASL:\n'];
						subplot(1,3,2);
					end
				else
					str = [str ' for normalized PET vs DSC:\n'];
					subplot(1,3,3);
				end
				fprintf(str);
				if (iVal == 1) || (iVal == 2 && iRoi == 3)
					fprintf('Mean GM CBF. %.5f*X+%.5f, CI %.5f -- %.5f, p=%.5f. Adjusted R2 %.3f\n',b(1),b(2),CI(1,1),CI(1,2),pval(1),stats.rSQadj);
					fprintf('Mean GM CBF. %.5f*X, CI %.5f -- %.5f, p=%.5f. Adjusted R2 %.3f\n',bnoInt(1),CInoInt(1,1),CInoInt(1,2),pvalnoInt(1),statsnoInt.rSQadj);
					fprintf('Mean relative difference in mean GM CBF: %.2f%%\n',mrd);
				end
				if iRoi == 3
					fprintf('Max CBF. %.5f*X+%.5f, CI %.5f -- %.5f, p=%.5f. Adjusted R2 %.3f\n',bmax(1),bmax(2),CImax(1,1),CImax(1,2),pvalmax(1),statsmax.rSQadj);
					fprintf('Mean relative difference in max CBF: %.2f%%\n',mrdmax);
				end
				fprintf('Voxel-wise CBF. %.5f*X+%.5f, CI %.5f -- %.5f, p=%.5f. Adjusted R2 %.3f\n',bAll(1),bAll(2),CIAll(1,1),CIAll(1,2),pvalAll(1),statsAll.rSQadj);
				fprintf('Voxel-wise CBF. %.5f*X, CI %.5f -- %.5f, p=%.5f. Adjusted R2 %.3f\n',bAllnoint(1),CIAllnoint(1,1),CIAllnoint(1,2),pvalAllnoint(1),statsAllnoint.rSQadj);
				fprintf('Mean relative difference in voxel-wise CBF: %.2f%%\n\n',mrdvoxel);
				%[b,bint,r,rint,stats] = regress(yy,[xx,ones(length(yy),1)]);
				
				
				hold on;
				underRat = 0.005; % undersampling factor for the scatter plot
				for iL = 1:size(listValRef,4)
					% Plot all scatter plots for different subjects in a different color
					xx = squeeze(listValRef{iMod,iVal,iRoi,iL});
					yy = squeeze(listValSrc{iMod,iVal,iRoi,iL});
					indXX = randi(length(xx),ceil(length(xx)*underRat),1);
					xx = xx(indXX);
					yy = yy(indXX);
					clrInd = 'rgbcmykr';
					plot(xx,yy,[clrInd(iL) 'o']);
				end
				hold off
				
			end
		end
		
		figure(3);subplot(2,4,4*(iMod-1)+iRoi);
		ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		if iRoi == 3 && iMod == 1
			indFit = find(squeeze(resVec(iMod,1,iRoi,:)).*[1; 1; 1; 1; 1; 0; 1; 1]);
			indOut = 6;
		else
			indOut = [];
			indFit = ind;
		end
		%groupMean = mean(resMean(iMod,1,iRoi,ind,1))/mean(resMean(iMod,2,iRoi,ind,1));
		%plot(squeeze(resMax(iMod,1,iRoi,ind,1))/groupMean,squeeze(resMax(iMod,1,iRoi,ind,2))/groupMean,'r+');hold on
		plot(squeeze(resMax(iMod,2,iRoi,ind,1)),squeeze(resMax(iMod,2,iRoi,ind,2)),'rx');hold on
		if ~isempty(indOut)
			plot(squeeze(resMax(iMod,2,iRoi,indOut,1)),squeeze(resMax(iMod,2,iRoi,indOut,2)),'rx');
		end
		X = [ones(length(indFit),1),squeeze(resMax(iMod,2,iRoi,indFit,1))];
		Y = squeeze(resMax(iMod,2,iRoi,indFit,2));
		sol = pinv(X)*Y;
		%imPseudoCBF = imPseudoCBF*sol(2)+sol(1);
		%plot([0.4,3.2],[0.4*sol(2)+sol(1),3.2*sol(2)+sol(1)],'r-');
		%plot(squeeze(resMax(iMod,3,iRoi,ind,1)),squeeze(resMax(iMod,3,iRoi,ind,2)),'bx');
		plot([0.4,3.2],[0.4,3.2],'k--');
		title(['max CBF ' strMod ' ' strRoi]);
		
		figure(4);sp=subplot(nMod,4,4*(iMod-1)+iRoi);ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		imagesc((squeeze(sum(resHist(iMod,1,iRoi,ind,:,:),4))).^0.4);hold on
		set(sp,'Layer','top','XTickLabel',{'25','50','75','100','125','150'});
		set(sp,'Layer','top','YTickLabel',{'25','50','75','100','125','150'});
		plot([1,60],[1,60],'r-');
		axis(sp,'xy');
		title(['hist CBF ' strMod ' ' strRoi]);
		
		figure(5);sp=subplot(nMod,4,4*(iMod-1)+iRoi);ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		imagesc((squeeze(sum(resHist(iMod,2,iRoi,ind,:,:),4))).^0.4);hold on
		set(sp,'Layer','top','XTickLabel',{'0.5','1','1.5','2','2.5','3'});
		set(sp,'Layer','top','YTickLabel',{'0.5','1','1.5','2','2.5','3'});
		plot([1,60],[1,60],'r-');
		axis(sp,'xy');
		title(['hist CBF ' strMod ' ' strRoi]);
		
		figure(6);sp=subplot(nMod,4,4*(iMod-1)+iRoi);ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		imagesc((squeeze(sum(resHist(iMod,3,iRoi,ind,:,:),4))).^0.4);hold on
		set(sp,'Layer','top','XTickLabel',{'0.5','1','1.5','2','2.5','3'});
		set(sp,'Layer','top','YTickLabel',{'0.5','1','1.5','2','2.5','3'});
		plot([1,60],[1,60],'r-');
		axis(sp,'xy');
		title(['hist CBF ' strMod ' ' strRoi]);
		
		if (iMod == 1) && (iRoi == 2)
			for iL = 1:8
				if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_PET.nii'))
					figure(7);sp=subplot(2,4,iL);
					imagesc((squeeze((resHist(iMod,1,iRoi,iL,:,:)))).^0.4);hold on
					set(sp,'Layer','top','XTickLabel',{'25','50','75','100','125','150'});
					set(sp,'Layer','top','YTickLabel',{'25','50','75','100','125','150'});
					plot([1,60],[1,60],'r-');
					axis(sp,'xy');
					title(['hist CBF ' strMod ' ' strRoi ' P' num2str(iL)]);
					
					figure(8);sp=subplot(2,4,iL);
					imagesc((squeeze((resHist(iMod,2,iRoi,iL,:,:)))).^0.4);hold on
					set(sp,'Layer','top','XTickLabel',{'0.5','1','1.5','2','2.5','3'});
					set(sp,'Layer','top','YTickLabel',{'0.5','1','1.5','2','2.5','3'});
					plot([1,60],[1,60],'r-');
					axis(sp,'xy');
					title(['hist CBF ' strMod ' ' strRoi ' P' num2str(iL)]);
				end
			end
		end
		
		if (iMod == 1) && (iRoi == 3)
			for iL = 1:8
				if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_PET.nii'))
					figure(11);sp=subplot(2,4,iL);
					imagesc((squeeze((resHist(iMod,1,iRoi,iL,:,:)))).^0.4);hold on
					set(sp,'Layer','top','XTickLabel',{'25','50','75','100','125','150'});
					set(sp,'Layer','top','YTickLabel',{'25','50','75','100','125','150'});
					plot([1,60],[1,60],'r-');
					axis(sp,'xy');
					title(['hist CBF ' strMod ' ' strRoi ' P' num2str(iL)]);
					
					figure(12);sp=subplot(2,4,iL);
					imagesc((squeeze((resHist(iMod,2,iRoi,iL,:,:)))).^0.4);hold on
					%set(sp,'Layer','top','XTickLabel',{'0.5','1','1.5','2','2.5','3'});
					set(sp,'Layer','top','XTickLabel',{'1','2','3'});
					set(sp,'Layer','top','YTickLabel',{'0.5','1','1.5','2','2.5','3'});
					plot([1,60],[1,60],'r-');
					axis(sp,'xy');
					title(['hist CBF ' strMod ' ' strRoi ' P' num2str(iL)]);
				end
			end
		end
		
		if (iMod == 2) && (iRoi == 3)
			for iL = 1:8
				if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},modDir,'Final_DSC.nii'))
					figure(13);sp=subplot(2,4,iL);
					imagesc((squeeze((resHist(iMod,1,iRoi,iL,:,:)))).^0.4);hold on
					set(sp,'Layer','top','XTickLabel',{'25','50','75','100','125','150'});
					set(sp,'Layer','top','YTickLabel',{'25','50','75','100','125','150'});
					plot([1,60],[1,60],'r-');
					axis(sp,'xy');
					title(['hist CBF ' strMod ' ' strRoi ' P' num2str(iL)]);
					
					figure(14);sp=subplot(2,4,iL);
					imagesc((squeeze((resHist(iMod,2,iRoi,iL,:,:)))).^0.4);hold on
					%set(sp,'Layer','top','XTickLabel',{'0.5','1','1.5','2','2.5','3'});
					set(sp,'Layer','top','XTickLabel',{'1','2','3'});
					set(sp,'Layer','top','YTickLabel',{'0.5','1','1.5','2','2.5','3'});
					plot([1,60],[1,60],'r-');
					axis(sp,'xy');
					title(['hist CBF ' strMod ' ' strRoi ' P' num2str(iL)]);
				end
			end
		end
		
		
		figure(9);subplot(nMod,4,4*(iMod-1)+iRoi);ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		plot(squeeze(resMean(iMod,1,iRoi,ind,1)),squeeze(resMean(iMod,1,iRoi,ind,2)),'r+');hold on
		plot([1,80],[1,80],'k--');
		title(['mean CBF ' strMod ' ' strRoi]);
		figure(10);subplot(nMod,4,4*(iMod-1)+iRoi);ind = find(squeeze(resVec(iMod,1,iRoi,:)));
		plot(squeeze(resMax(iMod,1,iRoi,ind,1)),squeeze(resMax(iMod,1,iRoi,ind,2)),'r+');hold on
		plot([1,200],[1,200],'k--');
		title(['max CBF ' strMod ' ' strRoi]);
	end
end
