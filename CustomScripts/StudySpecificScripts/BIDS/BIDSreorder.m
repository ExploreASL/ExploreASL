srcPath = '/pet/projekte/asl/data/BIDS/BIDS';
newPath = '/pet/projekte/asl/data/BIDS/newBIDS';

if ~exist(newPath, 'dir')
	mkdir(newPath);
end

fList = xASL_adm_GetFileList(srcPath, [], 'List', [], 1);

for iList = 1:length(fList)
	if ~exist(fullfile(newPath,fList{iList}), 'dir')
		mkdir(fullfile(newPath,fList{iList}));
	end
	system(['cp -r ' fullfile(srcPath,fList{iList}, 'rawdata', '*') ' ' fullfile(newPath,fList{iList})]);
end