function [ NewList ] = xASL_adm_CompareLists( list1, list2 )
%xASL_adm_CompareLists Compare 2 single dimension lists

    % Checks
    if  length(size(list1))>2 || length(size(list2))>2
        error('One or more lists have too many dimensions!');
    elseif size(list1,2)>1 || size(list2,2)>1
        error('One or more lists have too many columns, or should be transposed!');
    end

    % Get longest list
    if      length(list1)>length(list2)
            PrimList    = list1;
            SecList     = list2;
    else
            PrimList    = list2;
            SecList     = list1;
    end

    for iList   = 1:length(PrimList)
        if  ischar(PrimList{iList})      % skip if it doesn't contain characters, e.g. NaN
            NewList{iList,1}        = PrimList{iList};
            NewList{iList,2}        = 0;
            iList2  = 1;

            while   iList2<=length(SecList) && ~NewList{iList,2} % search for identical query
                if  strcmp( PrimList{iList} , SecList{iList2} )
                    NewList{iList,2}    = 1; % note found
                else
                    iList2 = iList2+1; % check next query
                end
            end
        end
    end

end

