function [PermutationList] = xASL_stat_UniquePairwisePermutations(GroupA, GroupB)
%xASL_stat_UniquePairwisePermutations Create list of pairwise permutations
%
% FORMAT: [PermutationList] = xASL_stat_UniquePairwisePermutations(GroupA, GroupB)
%
% INPUT:
%   GroupA      - numerical single column of indices (REQUIRED)
%   GroupB      - numerical single column of indices for two sample comparisons (OPTIONAL, default is one sample)
%
% OUTPUT:
%   PermutationList   - numerical matrix of two columns with indices, each row
%   containing a pairwise comparison of indices
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function lists for one or two samples of indices
% all possible permutations of indices, 
% performing the following steps:
%
% 1. One sample permutations
% 2. Two sample permutations    
% 3. Print conclusion
%
% PM: Allow entering one group only
% PM: could extend with xASL_qc_TanimotoCoeff
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE SINGLE SAMPLE: PermutationList = xASL_stat_UniquePairwisePermutations([1 2 3]);
% EXAMPLE TWO SAMPLES: PermutationList = xASL_stat_UniquePairwisePermutations([1 2 3], [4 5 6]);
% __________________________________
% Copyright 2015-2020 ExploreASL



%% 0. Admin
if nargin>2
    error('Too many input arguments');
elseif nargin>1 && isempty(GroupB)
    error('Second group should not be empty');
elseif nargin>1 && ~isnumeric(GroupB)
    error('Second group should be numerical');
end

if nargin<1 || isempty(GroupA)
    error('First input group missing');
elseif ~isnumeric(GroupA)
    error('First input group should be numerical');
end

%% 1. One sample permutations
if nargin<2
    fprintf('Obtain list of unique pair-wise permutations for single sample (1 group)\n');
    
    if length(GroupA)>33
        fprintf('This can take a while (this code needs optimization)\n');
    end

    % Create all partial permutations
    % A=[1:1:12]
    r=2;

     c = nchoosek(GroupA,r)';
     ncr = size(c,2);
     p = perms([1:r]);
     pr = size(p,1);
     p = reshape(p',1,[]);
     a = zeros(ncr*pr,r);
     for k = 1:ncr
       a((k-1)*pr+1:k*pr,:) = reshape(c(p,k),r,[])';
     end

    % Select the first one
    NoDoubles     = a(1,:);

    for iD=2:size(a,1)
        FOUND   = 0;
        for iF=1:size(NoDoubles,1)
            if  min(fliplr(a(iD,:))==NoDoubles(iF,:))
                FOUND   = 1;
            end
        end
        if ~FOUND % Add to NoDoubles
            NoDoubles(end+1,:)  = a(iD,:);
        end
    end

    PermutationList = fliplr(NoDoubles);
    
else
%% 2. Two sample permutations    
    fprintf('Obtain list of unique pair-wise permutations for two independent samples (2 groups)\n');
    
    IndicesA = (1:length(GroupA))';
    IndicesB = (1:length(GroupB))';

    ComparisonList = zeros(length(IndicesA)*length(IndicesB), 2);

    for iIndex=1:length(IndicesA)
        StartIndex = (iIndex-1)*length(IndicesB)+1;
        EndIndex = iIndex*length(IndicesB);
        ComparisonList(StartIndex:EndIndex, 1) = iIndex;
        ComparisonList(StartIndex:EndIndex, 2) = IndicesB;
    end

    PermutationList(:,1) = GroupA(ComparisonList(:,1))';
    PermutationList(:,2) = GroupB(ComparisonList(:,2))';
end

%% 3. Print conclusion
if nargin<2
    fprintf('%s\n', ['Created ' xASL_num2str(length(PermutationList)) ' permutations/combinations for ' xASL_num2str(length(GroupA)) ' indices']);
else
    fprintf('%s\n', ['Created ' xASL_num2str(length(PermutationList)) ' permutations/combinations for ' xASL_num2str(length(IndicesA)) ' (groupA) & ' xASL_num2str(length(IndicesB)) ' (groupB) indices']);
end


end