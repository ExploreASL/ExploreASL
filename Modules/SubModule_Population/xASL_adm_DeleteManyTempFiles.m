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
% Copyright 2015-2020 ExploreASL

Files2Del   = {'ATT_BiasField.nii' 'Mask_Template.nii' 'Mean_CBF_Template.nii' 'PseudoCBF.nii' 'RawTemplate.nii' 'VascularArtifact_Template.nii' 'mean_PWI_Clipped.nii' 'SliceGradient_extrapolated.nii' 'FoV.nii'};


for iP=1:length(Files2Del)
    xASL_TrackProgress(iP,length(Files2Del));
    [~, Ffile, Fext]    = xASL_fileparts(Files2Del{iP});
    
    if  strcmp(Fext,'.nii')
        File2Search     = [Ffile '\.(nii|nii\.gz)'];
    else
        File2Search     = [Ffile Fext];
    end
        
    PathList    = xASL_adm_GetFileList(x.D.ROOT,['^' File2Search '$'],'FPListRec',[0 Inf]);
    for iL=1:length(PathList)
        xASL_delete(PathList{iL});
    end

end


end
