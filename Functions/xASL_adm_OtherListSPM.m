function [OtherListSPM, OtherListOut] = xASL_adm_OtherListSPM(OtherList, bList4D)
%xASL_adm_OtherListSPM bPadComma1 is to add the ,1 to the end of the pathstring, which SPM uses
% to assign the first image of a 4D image array (OPTIONAL, DEFAULT = true)
% bList4D: boolean, true for listing multiple 4D volumes separately in the
% list (OPTIONAL, DEFAULT=true).
%
% FORMAT:       [OtherListSPM, OtherListOut] = xASL_adm_OtherListSPM(OtherList, bList4D)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  bPadComma1 is to add the ,1 to the end of the pathstring, which SPM uses
%               to assign the first image of a 4D image array (OPTIONAL, DEFAULT = true)
%               bList4D: boolean, true for listing multiple 4D volumes separately in the
%               list (OPTIONAL, DEFAULT=true).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

    if nargin<2 || isempty(bList4D)
        bList4D = true;
    end

    OtherListSPM = '';
    OtherListOut = '';

    if ~iscell(OtherList)
        OtherList = {OtherList};
    end
    
    for iList=1:numel(OtherList)
        if ~isempty(OtherList{iList})
            if ~xASL_exist(OtherList{iList},'file')
%                 fprintf('%s\n',[OtherList{iL} ' not existing, not coregistered']);
            else
                % xASL_spm_admin
                OtherList(iList) = xASL_spm_admin(OtherList{iList}, true);
                OtherList{iList} = OtherList{iList}(1:end-2);
                tIM = xASL_io_ReadNifti(OtherList{iList});
                if length(size(tIM.dat))>4
                    warning([OtherList{iList} ' had more than 4 dimensions']);
                else % add to list
                    OtherListOut{end+1,1} = [OtherList{iList}];
                    
                    if bList4D
                        for iS=1:size(tIM.dat,4)
                            OtherListSPM{end+1,1} = [OtherList{iList} ',' num2str(iS)];
                        end
                    else
                        OtherListSPM{end+1,1} = [OtherList{iList} ',1'];
                    end
                    
                    fprintf('%s\n',[OtherList{iList} ' is co-registered']);
                end
            end
        end
    end

    if isempty(OtherListSPM)
        OtherListSPM = {''};
    end
    if isempty(OtherListOut)
        OtherListOut = {''};
    end



end
