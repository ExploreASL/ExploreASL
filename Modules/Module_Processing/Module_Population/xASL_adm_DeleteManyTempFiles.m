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
% Copyright 2015-2023 ExploreASL

Files2Del = {'ATT_BiasField.nii' 'Mask_Template.nii' 'Mean_CBF_Template.nii' 'PseudoCBF.nii' 'RawTemplate.nii' 'VascularArtifact_Template.nii' 'mean_PWI_Clipped.nii' 'SliceGradient_extrapolated.nii' 'FoV.nii'};

if isfield(x,'D') && isfield(x.dir,'xASLDerivatives')
    fprintf('Deleting temporary files in native space subject folders:    ');
    for iP=1:length(Files2Del)
        xASL_TrackProgress(iP,length(Files2Del));
        [~, File2Search, Fext] = xASL_fileparts(Files2Del{iP});
        
        PathList = xASL_adm_GetFileList(x.dir.xASLDerivatives,['^' File2Search Fext '$'],'FPListRec',[0 Inf]);
        for iL=1:length(PathList)
            xASL_delete(PathList{iL});
        end
    end
    fprintf('\n');
else
    fprintf('Missing root directory, deleting temporary files failed...\n');
end

%% Now we do the same for the standard space files in the population folder
Files2Del_Population = {'rT1_ORI_' 'rT1_' 'rFLAIR_' '(m|)rc\dT1_' 'noSmooth_M0_' 'mean_control_' 'PWI_' 'SliceGradient_' 'SNR'};

fprintf('Deleting temporary files in population folder:    ');
for iPath=1:length(Files2Del_Population)
    xASL_TrackProgress(iP,length(Files2Del_Population));
    
    xASL_adm_DeleteFileList(x.D.PopDir,['^' Files2Del_Population{iPath} '.*$'], 1,[0 Inf]);
end
fprintf('\n');


end