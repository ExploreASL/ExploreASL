function [OtherListSPM, OtherListOut] = xASL_adm_OtherListSPM(OtherList)
%xASL_adm_OtherListSPM Takes care of the others list for registration
%functions
% bPadComma1 is to add the ,1 to the end of the pathstring, which SPM uses
% to assign the first image of a 4D image array (OPTIONAL, DEFAULT = true)

    OtherListSPM = '';
    OtherListOut = '';

    if ~iscell(OtherList)
        OtherList = {OtherList};
    end
    
    for iL=1:numel(OtherList)
        if ~isempty(OtherList{iL})
            if ~xASL_exist(OtherList{iL},'file')
%                 fprintf('%s\n',[OtherList{iL} ' not existing, not coregistered']);
            else
                % xASL_spm_admin
                OtherList(iL) = xASL_spm_admin(OtherList{iL}, true);
                OtherList{iL} = OtherList{iL}(1:end-2);
                tIM = xASL_io_ReadNifti(OtherList{iL});
                if length(size(tIM.dat))>4
                    warning([OtherList{iL} ' had more than 4 dimensions']);
                else % add to list
                    OtherListOut{end+1,1} = [OtherList{iL}];
                    for iS=1:size(tIM.dat,4)
                        OtherListSPM{end+1,1} = [OtherList{iL} ',' num2str(iS)];
                    end
                    fprintf('%s\n',[OtherList{iL} ' is co-registered']);
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
