function INDEX = xASL_adm_FindStrIndex(ARRAY, STRING)
% xASL_adm_FindStrIndex Similar to find, but then for a cell array filled with strings
% Only takes 4 dimensions
%
% FORMAT:       INDEX = xASL_adm_FindStrIndex(ARRAY, STRING)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Similar to find, but then for a cell array filled with strings.
%               Only takes 4 dimensions.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

nDim    = length(size(ARRAY));

if nDim==2 % correct that, if nDim should be 1, it also gives 2
    if size(ARRAY,2)==1
        nDim=1;
    end
end

next    = 1;

for i1=1:size(ARRAY,1)
    for i2=1:size(ARRAY,2)
        for i3=1:size(ARRAY,3)
            for i4=1:size(ARRAY,4)
                for i5=1:size(ARRAY,5)
                    for i6=1:size(ARRAY,6)
                        for i7=1:size(ARRAY,7)
                
                            if  strcmp(ARRAY{i1,i2,i3,i4,i5,i6,i7}, STRING)
                                INDEX(next,:)   = [i1 i2 i3 i4 i5 i6 i7];
                                next            = next+1;
                            end

                        end
                    end
                end
            end
        end
    end
end

% remove redundant dimensions

if  ~exist('INDEX','var')
    INDEX   = 0; % allocation, if no matches
else
    INDEX   = INDEX(:,1:nDim);
end




end

