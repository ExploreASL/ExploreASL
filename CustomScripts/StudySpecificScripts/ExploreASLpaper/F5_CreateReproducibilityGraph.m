%% CreateReproducibilityGraph
doStd = 0;
if doStd
	% standard settings
	ylims = [0,0,0,0];
	sizeLine = 0.5;
	sizeMarker = 6;
	replaceName{1,1} = 'None';
	plotSize = [2,2];
	fixedStepName = {};
else
	% paper settings
	ylims = [2.6,2.6,0.08,2];%[2,0.08,2.6,2];
	sizeLine = 3;
	sizeMarker = 10;
	replaceName = {'010_LinearReg_T1w2MNI','2.1_1';
		           '020_LinearReg_FLAIR2T1w','2.1_2';
				   '040_LST_Segment_FLAIR_WMH','2.1_3';
				   '060_Segment_T1w','2.2';
				   '080_Resample2StandardSpace','2.3';
				   '020_realign_ASL','3.1';
				   '025_register_ASL','3.3';
				   '040_realign_reslice_M0','3.4';
				   '050_quantification','3.5'};
    fixedStepName = {{'2.1_1','2.1_2','2.1_3','2.2','2.3'},{'3.1','3.3','3.4','3.5'}};
    plotSize = [1,4];
end
% MATLAB
% 4 on one line, denser
% Adjust size

% INKSCAPE
% reformat
% Headlines to columns/rows
% Legend out


%% ExploreASL_Master('',0);
pathRep = '/pet/projekte/asl/data/ExploreASLpaper/';
pathPre = 'RMS_Repro_';

comb(1,:) = {'EPAD_Lin18b_NEW_1_EPAD_Lin18b_NEW_2'};
comb(2,:) = {'EPAD_Lin18b_NEW_1_EPAD_Win18b_NEW_1'};
comb(3,:) = {'EPAD_Lin18b_OLD_1_EPAD_Lin18b_OLD_2'};
comb(4,:) = {'EPAD_Lin18b_OLD_1_EPAD_Win18b_OLD_1'};
comb(5,:) = {'EPAD_Mac14b_NEW_1_EPAD_Lin18b_NEW_1'};
comb(6,:) = {'EPAD_Mac16b_NEW_1_EPAD_Lin18b_NEW_1'};
comb(7,:) = {'EPAD_Mac16b_OLD_1_EPAD_Lin18b_OLD_1'};
comb(8,:) = {'EPAD_Win18b_NEW_1_EPAD_Win15a_NEW_1'};
comb(9,:) = {'EPAD_Win18b_NEW_1_EPAD_Win16b_NEW_1'};
comb(10,:) = {'EPAD_Win18b_OLD_1_EPAD_Win15a_OLD_1'};
comb(11,:) = {'EPAD_Win18b_OLD_1_EPAD_Win16b_OLD_1'};
comb(12,:) = {'EPAD_WinLin16b_NEW_1_EPAD_Lin18b_NEW_1'};
comb(13,:) = {'EPAD_WinLin16b_NEW_1_EPAD_Win18b_NEW_1'};
comb(14,:) = {'EPAD_WinLin16b_OLD_1_EPAD_Lin18b_OLD_1'};
comb(15,:) = {'EPAD_WinLin16b_OLD_1_EPAD_Win18b_OLD_1'};
comb(16,:) = {'Novice_Lin18b_NEW_1_Novice_Lin18b_NEW_2'};
comb(17,:) = {'Novice_Lin18b_NEW_1_Novice_Win18b_NEW_1'};
comb(18,:) = {'Novice_Lin18b_OLD_1_Novice_Lin18b_OLD_2'};
comb(19,:) = {'Novice_Lin18b_OLD_1_Novice_Win18b_OLD_1'};
comb(20,:) = {'Novice_Mac14b_OLD_1_Novice_Lin18b_OLD_1'};
comb(21,:) = {'Novice_Mac16b_NEW_1_Novice_Lin18b_NEW_1'};
comb(22,:) = {'Novice_Mac16b_OLD_1_Novice_Lin18b_OLD_1'};
comb(23,:) = {'Novice_Win18b_NEW_1_Novice_Win15a_NEW_1'};
comb(24,:) = {'Novice_Win18b_NEW_1_Novice_Win16b_NEW_1'};
comb(25,:) = {'Novice_Win18b_OLD_1_Novice_Win15a_OLD_1'};
comb(26,:) = {'Novice_Win18b_OLD_1_Novice_Win16b_OLD_1'};
comb(27,:) = {'Novice_WinLin16b_NEW_1_Novice_Lin18b_NEW_1'};
comb(28,:) = {'Novice_WinLin16b_NEW_1_Novice_Win18b_NEW_1'};
comb(29,:) = {'Novice_WinLin16b_OLD_1_Novice_Lin18b_OLD_1'};
comb(30,:) = {'Novice_WinLin16b_OLD_1_Novice_Win18b_OLD_1'};
comb(31,:) = {'Sleep_Lin18b_NEW_1_Sleep_Lin18b_NEW_2'};
comb(32,:) = {'Sleep_Lin18b_NEW_1_Sleep_Win18b_NEW_1'};
comb(33,:) = {'Sleep_Lin18b_OLD_1_Sleep_Lin18b_OLD_2'};
comb(34,:) = {'Sleep_Lin18b_OLD_1_Sleep_Win18b_OLD_1'};
comb(35,:) = {'Sleep_Mac14b_NEW_1_Sleep_Lin18b_NEW_1'};
comb(36,:) = {'Sleep_Mac14b_OLD_1_Sleep_Lin18b_OLD_1'};
comb(37,:) = {'Sleep_Mac16b_NEW_1_Sleep_Lin18b_NEW_1'};
comb(38,:) = {'Sleep_Mac16b_OLD_1_Sleep_Lin18b_OLD_1'};
comb(39,:) = {'Sleep_Win18b_NEW_1_Sleep_Win15a_NEW_1'};
comb(40,:) = {'Sleep_Win18b_NEW_1_Sleep_Win16b_NEW_1'};
comb(41,:) = {'Sleep_Win18b_OLD_1_Sleep_Win15a_OLD_1'};
comb(42,:) = {'Sleep_Win18b_OLD_1_Sleep_Win16b_OLD_1'};
comb(43,:) = {'Sleep_WinLin16b_NEW_1_Sleep_Lin18b_NEW_1'};
comb(44,:) = {'Sleep_WinLin16b_NEW_1_Sleep_Win18b_NEW_1'};
comb(45,:) = {'Sleep_WinLin16b_OLD_1_Sleep_Lin18b_OLD_1'};
comb(46,:) = {'Sleep_WinLin16b_OLD_1_Sleep_Win18b_OLD_1'};
comb(47,:)= {'EPAD_Mac16b_NEW_1_EPAD_Win18b_NEW_1'};
comb(48,:)= {'EPAD_Mac16b_OLD_1_EPAD_Win18b_OLD_1'};
comb(49,:)= {'Sleep_Mac16b_NEW_1_Sleep_Win18b_NEW_1'};
comb(50,:)= {'Sleep_Mac16b_OLD_1_Sleep_Win18b_OLD_1'};
comb(51,:)= {'Novice_Mac16b_NEW_1_Novice_Win18b_NEW_1'};
comb(52,:)= {'Novice_Mac16b_OLD_1_Novice_Win18b_OLD_1'};
comb(53,:)= {'EPAD_Win15a_NEW_1_EPAD_Lin18b_NEW_1'};
comb(54,:)= {'EPAD_Win15a_OLD_1_EPAD_Lin18b_OLD_1'};
comb(55,:)= {'Sleep_Win15a_NEW_1_Sleep_Lin18b_NEW_1'};
comb(56,:)= {'Sleep_Win15a_OLD_1_Sleep_Lin18b_OLD_1'};
comb(57,:)= {'Novice_Win15a_NEW_1_Novice_Lin18b_NEW_1'};
comb(58,:)= {'Novice_Win15a_OLD_1_Novice_Lin18b_OLD_1'};
comb(59,:)= {'EPAD_Mac16b_NEW_2_EPAD_Mac16b_NEW_3'};
comb(60,:)= {'Sleep_Mac16b_NEW_2_Sleep_Mac16b_NEW_3'};
comb(61,:)= {'Novice_Mac16b_NEW_2_Novice_Mac16b_NEW_3'};


for iF = 1:size(comb,1)%all
	%59:61
	%47:58%table2
	%[6,7,8,10,21,22,23,25,37,38,39,41]%table1
	%[2,4,17,19,32,34]%graphs
	%1:size(comb,1)%all
	clear Mat; clear RMS;clear RemoveStrings;clear ModInd;
	clear RMSplot;clear FieldNames;
	close all
	Mat         = load(fullfile(pathRep,[pathPre,comb{iF,1} '.mat']));
	RMS         = Mat.RMS;
	
	%% First remove all pipeline steps that are not significant for reproducibility (i.e. the QC steps)
	StepsToSkip = zeros(length(RMS),1);
	RemoveStrings = {'070_CleanUpWMH_SEGM','090_GetVolumetrics','100_VisualQC','999_ready','060_VisualQC','030_reslice_ASL','035_PreparePV'};
	for ii = 1:length(RemoveStrings)
		idxRemove = find(strcmp(RMS(:,3),RemoveStrings{ii}));
		if ~isempty(idxRemove)
			StepsToSkip(idxRemove) = 1;
		end
	end
	
	%StepsToSkip(17:end)     = 1; % remove the Population module
	%StepsToSkip([6 8 9])   = 1; % remove the Structural module Get_WMH_volume, TissueVolume & Visualize steps
	%StepsToSkip(16)         = 1; % remove the ASL module QC step
	
	RMS                     = RMS(~StepsToSkip,:);
	
	% Go through the certain measures (in setZeroString), and set them to 0 if for all steps prior to setZeroStep
	setZeroString = {'WMH_SEGM','c1T1','c2T1','PV_pGM','qCBF','rc1T1'};
	setZeroStep   = [40          60     60     35       50     80];
	% Cycle through all steps
	for iR = 1:size(RMS,1)
		stepNum = str2num(RMS{iR,3}(1:3));
		% Cycle through all removed String
		for iS = 1:length(setZeroString)
			% Too small step, start removing
			if stepNum<setZeroStep(iS)
				% Cycle through the 3 subfields of RMS
				for iA = 4:6
					% Cycle through the list
					for iU = 1:size(RMS{iR,iA},1)
						% Found the matching string - start zeroing...
						if ~isempty(strfind(RMS{iR,iA}{iU,1},setZeroString{iS}))
							RMS{iR,iA}{iU,2} = 0;
							RMS{iR,iA}{iU,3} = 0;
						end
					end
				end
			end
		end
	end
	
	%% Define the files we want to track
	% cells from each module
	% clear iMod
	ModInd{1}    = find(cellfun(@(x) ~isempty(regexp(x,'module_Structural')), RMS(:,2))); % Module Indices
	ModInd{2}    = find(cellfun(@(x) ~isempty(regexp(x,'module_ASL')),        RMS(:,2))); % Module Indices
	
	% Define files
	%Files2Track{1}      = {'T1'  'c1T1'  'FLAIR' 'WMH_SEGM' 'rc1T1'};% this can be automated later
	Files2Track{1,1}      = {'rc1T1'  'rT1' 'rFLAIR' 'rWMH_SEGM'};
	Files2Track{1,2}      = {'c1T1'  'T1' 'FLAIR' 'WMH_SEGM'};
	Files2Track{2,1}      = {'rM0' 'qCBF' 'PV_pGM'};
	Files2Track{2,2}      = {'M0'  'CBF'  'PV_pGM'};
	
	%Files2Track{1}      = {'T1'  'c1T1'  'c2T1' 'y_T1'};  % this can be automated later
	%Files2Track{2}      = {'ASL4D' 'M0' 'mean_control' 'PWI' 'y_ASL' 'CBF'};
	% Files2Track{3}      = {'rFLAIR' 'rT1' 'rWMH_SEGM' 'rc1T1' 'rc2T1'};
	% Files2Track{4}      = {'noSmooth_M0' 'mean_control' 'SD' 'PWI' 'PV_pGM' 'PV_pWM' 'qCBF_untreated'};
	
	ModuleName          = {'Structural' 'ASL'};
	
	close;
	figure('units','normalized','outerposition',[0 0 1 1]);
	
	for iMod=1:2 % loop across modules (Structural or ASL)
		% >>>>>>>>>>>>>>>>>> First build RMSplot structure with filenames (keys) & values
		for iT = 1:2
			for iN=1:length(Files2Track{iMod,iT})  % loop across files
				NiiFileName    = Files2Track{iMod,iT}{iN};
				for iC=1:length(ModInd{iMod}) % loop over pipeline steps
					% Find name of file in the cell
					
					Index2          = ModInd{iMod}(iC);
					
					CellContent     = RMS{Index2,3+iMod}(:,1); % find FileNames for each pipeline step
					FileIndex       = find(cellfun(@(x) ~isempty(regexp(x,['^(\d{3}|)' NiiFileName '$'])),CellContent));
					if  isempty(FileIndex)
						% Does not exist in the correct module - let's check the Population module
						CellContent     = RMS{Index2,6}(:,1); % find FileNames for each pipeline step
						FileIndex       = find(cellfun(@(x) ~isempty(regexp(x,['^' NiiFileName ])),CellContent));
						if isempty(FileIndex)
							% Not in the population module
							RMSplot{iMod}.(NiiFileName)(iC,1:2)   = NaN;
						else
							if ~isempty(strfind(RMS{Index2,6}{FileIndex(1),2},'Reference '))
								RMSplot{iMod}.(NiiFileName)(iC,1:2)   = 0;
							else
								RMSplot{iMod}.(NiiFileName)(iC,1)     = RMS{Index2,6}{FileIndex(1),2};
								RMSplot{iMod}.(NiiFileName)(iC,2)     = RMS{Index2,6}{FileIndex(1),3};
							end
						end
						
					elseif ~isempty(strfind(RMS{Index2,3+iMod}{FileIndex,2},'Reference '))
						RMSplot{iMod}.(NiiFileName)(iC,1:2)   = 0;
					else
						RMSplot{iMod}.(NiiFileName)(iC,1)     = RMS{Index2,3+iMod}{FileIndex,2};
						RMSplot{iMod}.(NiiFileName)(iC,2)     = RMS{Index2,3+iMod}{FileIndex,3};
					end
				end
			end
		end
		% >>>>>>>>>>>>>>>>>> Now start plotting
		WhichRMS  = {'Voxel-wise' 'transformation'};
		FieldNames{iMod}  = fieldnames(RMSplot{iMod});
		
		StepName                = RMS(ModInd{iMod},3)';
		for iS=1:length(StepName)
			StepName{iS}        = char(StepName{iS});
			for iR = 1:size(replaceName,1)
				if strcmp(StepName{iS},replaceName{iR,1})
					StepName{iS} = replaceName{iR,2};
				end
			end
		end
		
		for iW=1:2
			if doStd
				plotNr = (iMod-1)*2+iW;
			else
				plotNr = (iW-1)*2 + iMod;
			end
			plot1 = subplot(plotSize(1),plotSize(2),plotNr);
			
			
			%          b     blue          .     point              -     solid
			%          g     green         o     circle             :     dotted
			%          r     red           x     x-mark             -.    dashdot
			%          c     cyan          +     plus               --    dashed
			%          m     magenta       *     star             (none)  no line
			%          y     yellow        s     square
			%          k     black         d     diamond
			%          w     white
			
			Colors          = {'b*-'  'g*-'   'r*-' 'c*-'  'm*-'     'k*-'};
			PrintColors     = {'blue' 'green' 'red' 'cyan' 'magenta' 'black'};
			xl = 1; 
			for iL=1:length(Files2Track{iMod,iW})%FieldNames{iMod,iW})
				dispName = Files2Track{iMod,iW}{iL};%FieldNames{iMod,iW}{iL};
				% Replace the first underscore
				indUnd = strfind(dispName,'_');
				if ~isempty(indUnd)
					dispName = [dispName(1:(indUnd(1)-1)), '\_', dispName((indUnd(1)+1):end)];
				end
				
				xl = max(xl,size(RMSplot{iMod}.(Files2Track{iMod,iW}{iL}),1));
				% Display the transformations y_ only for the transformation step
				doPlot = 1;
				if (iW == 1) 
					if strcmp(Files2Track{iMod,iW}{iL}(1:2),'y_')
						doPlot = 0;
					end
				else
					if strcmp(dispName,'c1T1') || strcmp(dispName,'WMH\_SEGM') || strcmp(dispName,'PV\_pGM')
						doPlot = 0;
					end
				end
					
				if doPlot
					if doStd
						plotLoc = plot(RMSplot{iMod}.(Files2Track{iMod,iW}{iL})(:,iW),Colors{iL},'MarkerSize',sizeMarker,'LineWidth',sizeLine);
					else
						pltTmp = zeros(size(fixedStepName{iMod},2),1);
						xl = max(xl,size(pltTmp,1));
						for iS = 1:size(pltTmp,1)
							% For each output FixedStepName - find if the recorded variant exists
							idx = find(strcmp(fixedStepName{iMod}{iS},StepName));
							if ~isempty(idx)
								% Found it
								pltTmp(iS,1) = RMSplot{iMod}.(Files2Track{iMod,iW}{iL})(idx(1),iW);
							else
								% Otherwise
								if iS == 1
									% Set to zero for the first
									pltTmp(iS,1) = 0;
								else
									% Set to previous for the next
									pltTmp(iS,1) = pltTmp(iS-1,1);
								end
							end
						end
						plotLoc = plot(pltTmp,Colors{iL},'MarkerSize',sizeMarker,'LineWidth',sizeLine);
					end
					set(plotLoc(1),'DisplayName',dispName);
				end
				hold on
			end
			
			
			ax                      = gca;
			ax.XTickLabelRotation   = 90;
			if doStd
				ax.XTick                = 1:length(StepName);
				ax.XTickLabel           = StepName;
			else
				ax.XTick                = 1:length(fixedStepName{iMod});
				ax.XTickLabel           = fixedStepName{iMod};
				ax.XTick = [];
			end
			ax.TickLabelInterpreter = 'none';
			
			if doStd
				xlabel('ExploreASL pipeline steps');
			end
			if  iW==1
				if doStd || (iMod == 1)
					ylabel('RMS difference (%)');
				end
			else
				if doStd || (iMod == 1)
					ylabel('Mean displacement (mm)');
				end
			end
			set(plot1,'FontSize',16);
			if doStd
				legend(plot1);
			else
				if length(ylims) >= plotNr
					ylim(plot1,[0 ylims(plotNr)]);
				end
				box(plot1,'off');
			end
			
			xlim(plot1,[1 xl]);
			if doStd
				title([ModuleName{iMod} ' module, ' WhichRMS{iW}]);
			end
		end % for iN=1:length(Files2Track{iMod})
	end % for iMod=1:2
	
	%fprintf('%s\n','<<<<<<<<<<<<<<<Legend>>>>>>>>>>>>>>>>>');
	%for iMod=1:2
	%	fprintf('%s\n',[ModuleName{iMod} '  module:']);
	%	for iL=1:length(FieldNames{iMod})
	%		fprintf('%s\n',[PrintColors{iL} ' = ' FieldNames{iMod,iW}{iL}]);
	%	end
	%end
	if ~doStd
		set(gcf, 'PaperUnits', 'inches');
		set(gcf, 'PaperSize' , [8.2677 4 ]);
		set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.33]);
	end
	saveas(gcf,fullfile(pathRep,[pathPre,comb{iF,1} '.png']));
	saveas(gcf,fullfile(pathRep,[pathPre,comb{iF,1} '.epsc']));
	xASL_Move(fullfile(pathRep,[pathPre,comb{iF,1} '.epsc']),fullfile(pathRep,[pathPre,comb{iF,1} '.eps']),1);
end	




%% Check what CBF difference this accounts for
% JP - Don't know where these numbers come from. Now I calculate the change as a mean difference divided by the mean value of both 
% = mean percentage change - percentages are calculated voxel-wise and added.

% For y_T1 - 1/3 of the volumes are cut and mean mm displacement between both transformations is calculated and given in mm
% For transformation matrices - 5 points on the middle plane are transformed (center + 1/3 in all directions) and their mean displacement between
% the two transformations is given in mm

% In MNI, 5.9 stands for 3.5 mL/100g/min difference -> 3.5/56 = 6.3 % difference

% In Native space, 1.9 with registration 0.8 stands for 8.2 mL/100g/min difference -> 8.2/56 = 15%

