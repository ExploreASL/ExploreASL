% Compare different combinations of processed data and store the RMS-matrices apart

% Specify the input directory and directory pre-fix
pathRep = '/home/janpetr/code/ExploreASL';

pathPre = 'TestDataSet_';

% Run CompareDataSets on these combinations
comb(1,:) = {'NEW_Linux18b','NEW_Linux18b_2'};
comb(2,:) = {'NEW_Linux18b','NEW_Win18b'};
comb(3,:) = {'NEW_Linux18b_2','NEW_Win18b_2'};
comb(4,:) = {'NEW_Win18b_2','NEW_Win18b'};
comb(5,:) = {'NEW_Win18b','NEW_Win11a'};
comb(6,:) = {'NEW_Win18b','NEW_Win14a'};
comb(7,:) = {'NEW_Linux18b','NEW_Win11a'};
comb(8,:) = {'NEW_Linux18b','NEW_Win14a'};
comb(9,:) = {'OLD_Linux18b','OLD_Linux18b_2'};
comb(10,:) = {'OLD_Linux18b','OLD_Win18b'};
comb(11,:) = {'OLD_Linux18b_2','OLD_Win18b_2'};
comb(12,:) = {'OLD_Win18b_2','OLD_Win18b'};
comb(13,:) = {'OLD_Win18b','OLD_Win11a'};
comb(14,:) = {'OLD_Win18b_2','OLD_Win14a'};
comb(15,:) = {'OLD_Linux18b','OLD_Win11a'};
comb(16,:) = {'OLD_Linux18b','OLD_Win14a'};
comb(17,:) = {'NEW_Linux18b_3','NEW_Win18b_3'};
comb(18,:) = {'OLD_Linux18b_3','OLD_Win18b_3'};
comb(19,:) = {'NEW_Linux18b_4','NEW_Win18b_4'};
comb(20,:) = {'OLD_Linux18b_4','OLD_Win18b_4'};

locDir = pwd;

for iC = [2,3,10,11,17:20]%1:size(comb,1)
	% Rename the first to _processed
	system(['mv ' fullfile(pathRep,[pathPre comb{iC,1}]) ' ' fullfile(pathRep,[pathPre 'Processed'])]);

	% Run the comparison
	cd('/home/janpetr/code/ExploreASL');
	ExploreASL_Master(fullfile(pathRep,[pathPre comb{iC,2}],'DATA_PAR.m'),1,1);
	cd(locDir);
	% Rename the directories back
	system(['mv ' fullfile(pathRep,[pathPre 'Processed']) ' ' fullfile(pathRep,[pathPre comb{iC,1}])]);
	
	% Copy out the RMS_Reproducibility.mat and rename it accordingly
	
	system(['cp ' fullfile(pathRep,[pathPre comb{iC,2}],'RMS_Reproducibility.mat') ' ' fullfile(pathRep,['RMS_Repro_' comb{iC,1} '_' comb{iC,2} '.mat'])]);
end