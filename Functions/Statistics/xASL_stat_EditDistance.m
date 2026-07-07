function d = xASL_stat_EditDistance(s1, s2, varargin)
%xASL_stat_EditDistance Damerau-Levenshtein edit distance between two strings
%
% FORMAT: d = xASL_stat_EditDistance(s1, s2)
%         d = xASL_stat_EditDistance(s1, s2, DelCost, InsCost, ReplCost)
%
% INPUT:
%   s1, s2         - horizontal char arrays (REQUIRED)
%   DelCost        - cost of deletion (OPTIONAL, DEFAULT=1)
%   InsCost        - cost of insertion (OPTIONAL, DEFAULT=1)
%   ReplCost       - cost of replacement (OPTIONAL, DEFAULT=1)
%
% OUTPUT:
%   d              - edit distance (minimum cost to convert s1 to s2),
%                    allowing deletion, insertion, substitution, and
%                    transposition of two adjacent characters
%
% DESCRIPTION: Implementation of Damerau-Levenshtein distance
%               using the standard recurrence relation.
%               Useful for fuzzy-matching typos in BIDS sidecar string
%               fields (e.g., PulseSequenceType, Manufacturer).
%
% EXAMPLE: xASL_stat_EditDistance('sprial','spiral') returns 1
%          xASL_stat_EditDistance('Seimen','Siemens') returns 2
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% Parse optional costs (default 1 for each operation)
%   Called as: EditDistance(s1, s2) or EditDistance(s1, s2, del, ins, sub)
if nargin < 3
    costDel = 1; costIns = 1; costSub = 1;
elseif nargin == 5
    [costDel, costIns, costSub] = deal(varargin{1:3});
else
    error('Expected 2 or 5 inputs: (s1, s2) or (s1, s2, DelCost, InsCost, ReplCost).');
end

%% Validate inputs
if ~(ischar(s1) && (isrow(s1) || isempty(s1))) || ~(ischar(s2) && (isrow(s2) || isempty(s2)))
    error('Inputs must be horizontal character vectors.');
end

n1 = length(s1);
n2 = length(s2);

%% Allocate distance matrix
% D(r, c) holds the edit distance between s1(1:r-1) and s1(1:c-1).
% Row 1 / column 1 correspond to the empty-prefix case.
D = zeros(n1 + 1, n2 + 1);

%% Base cases: transforming to/from the empty string
for r = 1:n1
    D(r+1, 1) = D(r, 1) + costDel;
end
for c = 1:n2
    D(1, c+1) = D(1, c) + costIns;
end

%% Fill the matrix row by row (Damerau-Levenshtein recurrence)
for r = 1:n1
    for c = 1:n2

        % Substitution cost: 0 when characters already match
        if s1(r) == s2(c)
            subCost = 0;
        else
            subCost = costSub;
        end

        % Three standard edit operations
        fromAbove  = D(r,   c+1) + costDel;   % delete s1(r)
        fromLeft   = D(r+1, c)   + costIns;   % insert s2(c)
        fromDiag   = D(r,   c)   + subCost;   % match or substitute

        D(r+1, c+1) = min([fromAbove, fromLeft, fromDiag]);

        % Transposition of two adjacent characters (Damerau extension):
        % s1(r) matches s2(c-1) AND s1(r-1) matches s2(c)
        if r >= 2 && c >= 2 && s1(r) == s2(c-1) && s1(r-1) == s2(c)
            D(r+1, c+1) = min(D(r+1, c+1), D(r-1, c-1) + costSub);
        end
    end
end

d = D(n1+1, n2+1);

end