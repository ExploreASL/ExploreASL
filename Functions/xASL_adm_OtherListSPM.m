function [xASL_adm_OtherListSPM, OtherListOut] = xASL_adm_OtherListSPM(OtherList)
%xASL_adm_OtherListSPM Takes care of the others list for rigid-body & affine registration

    xASL_adm_OtherListSPM = '';
    OtherListOut = '';

    for iL=1:numel(OtherList)
        if ~isempty(OtherList{iL})
            if ~exist(OtherList{iL},'file')
%                 fprintf('%s\n',[OtherList{iL} ' not existing, not coregistered']);
            else
                % xASL_spm_admin
                OtherList(iL) = xASL_spm_admin(OtherList{iL});
                OtherList{iL} = OtherList{iL}(1:end-2);
                tIM = xASL_io_ReadNifti(OtherList{iL});
                if  length(size(tIM.dat))>4
                    warning([OtherList{iL} ' had more than 4 dimensions']);
                else % add to list
                    OtherListOut{end+1,1} = [OtherList{iL}];
                    for iS=1:size(tIM.dat,4)
                        xASL_adm_OtherListSPM{end+1,1} = [OtherList{iL} ',' num2str(iS)];
                    end
                    fprintf('%s\n',[OtherList{iL} ' is co-registered']);
                end
            end
        end
    end

    if isempty(xASL_adm_OtherListSPM)
        xASL_adm_OtherListSPM = {''};
    end
    if isempty(OtherListOut)
        OtherListOut = {''};
    end



end
