function xASL_adm_DeleteManyTempFiles(x)
%xASL_adm_DeleteManyTempFiles This function removes as many files as possible
%
% FORMAT:       xASL_adm_DeleteManyTempFiles(x)
% 
% INPUT:        x          - x structure
%
% OUTPUT:       n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function removes as many files as possible.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_adm_DeleteManyTempFiles(x);
% __________________________________
% Copyright (c) 2015-2021 ExploreASL


    %% Delete temporary files
    Files2Del   = {'ATT_BiasField.nii' 'Mask_Template.nii' 'Mean_CBF_Template.nii' 'PseudoCBF.nii' 'RawTemplate.nii' 'VascularArtifact_Template.nii' 'mean_PWI_Clipped.nii' 'SliceGradient_extrapolated.nii' 'FoV.nii'};

    % Check if root directory is defined
    if isfield(x,'D') && isfield(x.D,'ROOT')
        fprintf('Deleting temporary files:    ');
        for iP=1:length(Files2Del)
            % Show progress
            xASL_TrackProgress(iP,length(Files2Del));
            [~, File2Search, Fext] = xASL_fileparts(Files2Del{iP});
            % Path which files to delete
            PathList = xASL_adm_GetFileList(x.D.ROOT,['^' File2Search Fext '$'],'FPListRec',[0 Inf]);
            % Iterate over file list
            for iL=1:length(PathList)
                xASL_delete(PathList{iL});
            end
        end
        fprintf('\n');
    else
        fprintf('Missing root directory, deleting temporary files failed...\n');
    end

end
