function [ OtherList ] = xASL_adm_Remove_1_SPM( OtherList )
%Append_1_SPM % Remove ,1 at end of OtherLists, if exists.
% These are appended in CoregInit, OldNormalizeWrapper etc,
% since this should allow 4rd dim (e.g. as in ASL4D)


    for iL=1:size(OtherList,1)
        FoundStr    = strfind( OtherList{iL}, ',1');
        if ~isempty(FoundStr)
            FoundStr    = FoundStr(length(FoundStr));
            if  length(OtherList{iL})==FoundStr+1
                OtherList{iL}   = OtherList{iL}(1:end-2);
            end
        end
    end


end

