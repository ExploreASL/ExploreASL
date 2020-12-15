function [DiceCoeff] = xASL_stat_PairwiseDice(GroupA, GroupB)
%PairwiseDice Obtain pairwise Dice coefficient for all permutations between groupA & groupB
%
% FORMAT: [DiceCoeff] = xASL_stat_PairwiseDice(GroupA, GroupB)
%
% INPUT:
%   GroupA      - cell structure with paths or images (REQUIRED)
%   GroupB      - cell structure with paths or images (REQUIRED)
%
% OUTPUT:
%   DiceCoeff   - numerical list of Dice coefficients
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function obtains for two lists of images Dice
% coefficients, for all possible permutations of both lists, by the
% following steps:
% 1. Admin (check cell, image exist etc)
% 2. Obtain matrix of pair-wise permutations
% 3. Obtain DICE scores
%
% PM: Allow entering one group only
% PM: could extend with xASL_qc_TanimotoCoeff
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: DiceCoeff = xASL_stat_PairwiseDice({'/MyPath/Image1.nii';'/MyPath/Image2.nii'}, {'/MyPath/Image3.nii';'/MyPath/Image4.nii'});
% 
% __________________________________
% Copyright 2015-2020 ExploreASL

%% 0. For testing with CICERO data, remove later
%
% e.g. A = healthy maskVascular
% e.g. B = sten-ASYM maskVascular
%
% SubjectsA = x.SUBJECTS(x.S.SetsID([1:2:end],4)==1);
% SubjectsB = x.SUBJECTS(x.S.SetsID([1:2:end],4)==4);
% 
% for iSubject=1:length(SubjectsA)
%     GroupA{iSubject,1} = fullfile(x.D.PopDir, ['MaskVascular_' SubjectsA{iSubject} '_ASL_1.nii']);
% end
% for iSubject=1:length(SubjectsB)
%     GroupB{iSubject,1} = fullfile(x.D.PopDir, ['MaskVascular_' SubjectsB{iSubject} '_ASL_1.nii']);
% end


%% 1. Admin (check cell, image exist etc)
if ~iscell(GroupA) && isnumeric(GroupA) && ndims(GroupA)==4
    fprintf('Group A matrix detected, converting to cell structure\n');
    TempMatrix = GroupA;
    GroupA = '';
    for iCell=1:size(TempMatrix, 4)
        GroupA{iCell} = TempMatrix(:,:,:,iCell);
    end
end

if ~iscell(GroupB) && isnumeric(GroupB) && ndims(GroupB)==4
    fprintf('Group B matrix detected, converting to cell structure\n');
    TempMatrix = GroupB;
    GroupB = '';
    for iCell=1:size(TempMatrix, 4)
        GroupB{iCell} = TempMatrix(:,:,:,iCell);
    end
end


if ~iscell(GroupA) || ~iscell(GroupB)
    error('Input images/paths should be a cell structure');
end

InitialLengths = [length(GroupA), length(GroupB)];

ExcludeGroupA = zeros(InitialLengths(1), 1);
ExcludeGroupB = zeros(InitialLengths(2), 1);

for iSubject=1:length(GroupA)
    if ~isnumeric(GroupA{iSubject}) && ~xASL_exist(GroupA{iSubject}, 'file')
        warning([GroupA{iSubject} ' did not exist, excluding']);
        ExcludeGroupA(iSubject, 1) = 1;
    end
end
for iSubject=1:length(GroupB)
    if ~isnumeric(GroupA{iSubject}) && ~xASL_exist(GroupB{iSubject}, 'file')
        warning([GroupB{iSubject} ' did not exist, excluding']);
        ExcludeGroupB(iSubject, 1) = 1;
    end
end

fprintf('%s\n', ['Group A: found ' xASL_num2str(InitialLengths(1)-sum(ExcludeGroupA)) ' out of ' xASL_num2str(InitialLengths(1)) ' scans']);
fprintf('%s\n', ['Group B: found ' xASL_num2str(InitialLengths(2)-sum(ExcludeGroupB)) ' out of ' xASL_num2str(InitialLengths(2)) ' scans']);
fprintf('%s\n', ['In total: ' xASL_num2str(sum(ExcludeGroupA)+sum(ExcludeGroupB)) ' missing scans, excluded from calculations']);

% Exclude scans
PathsGroup{1} = GroupA(~ExcludeGroupA);
PathsGroup{2} = GroupB(~ExcludeGroupB);
GroupA = NaN;
GroupB = NaN;

%% 2. Obtain matrix of pair-wise permutations
IndicesA = (1:length(PathsGroup{1}))';
IndicesB = (1:length(PathsGroup{2}))';

PermutationList = xASL_stat_UniquePairwisePermutations(IndicesA, IndicesB);

%% 3. Obtain DICE scores
fprintf('%s\n', 'Obtaining DICE scores:   ');
for iPerm=1:length(PermutationList)
    xASL_TrackProgress(iPerm, length(PermutationList));
    for iMask=1:size(PermutationList,2)
        % Load mask
        PathMask{iMask} = PathsGroup{iMask}{PermutationList(iPerm,iMask)};
        MaskImage{iMask} = xASL_io_Nifti2Im(PathMask{iMask});
        % Check if it is a binary mask
        NonLogical{iMask} = sum(MaskImage{iMask}(:)<0) + sum(MaskImage{iMask}(:)>1);
        if NonLogical{iMask}~=0
            warning([PathMask{iMask} was not a binary mask']);
        end
        % Convert to a binary mask
        MaskImage{iMask} = logical(MaskImage{iMask});
    end

    DiceCoeff(iPerm,1) = xASL_im_ComputeDice(MaskImage{1}, MaskImage{2});
end
fprintf('\n');
DiceMean = xASL_stat_MeanNan(DiceCoeff);
DiceSD = xASL_stat_StdNan(DiceCoeff);
fprintf('%s\n', ['Dice coeff was (mean +/- SD) : ' xASL_num2str(DiceMean) ' +/- ' xASL_num2str(DiceSD)]);


end