function xASL_spm_deface(PathIn, bReplace)
%xASL_spm_deface Defaces an anatomical image to preserve anonimity

% FORMAT: xASL_spm_deface(PathIn, bReplace)
% 
% INPUT:
%   PathIn   - path to anatomical NIfTI file that will be defaced. .nii & .nii.gz are both allowed (REQUIRED)
%   bReplace - boolean to specify whether the original anatomical NIfTI will be
%              overwritten (true) or whether the defaced result will be
%              stored under 'anon_FileName' (false) (OPTIONAL,
%              DEFAULT=false)
% OUTPUT n/a
% OUTPUT file:
%   Defaced NIfTI (but always .nii, not .nii.gz format)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function removes the face from an anatomical NIfTI
%              image, e.g. T1w or FLAIR, for disidentification/privacy purposes.
%              When this script is run after the ExploreASL structural
%              module, it does a pretty good job even for 2D images.
%              However, note that this can always fail, strip part of
%              the brain, or change the output of pipelines. So best not
%              to compare results from defaced and non-defaced images.
%              Also, note that defacing makes it difficult to ensure that
%              the FLAIR and T1w are from the same subject.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_deface('PathToMyStudy/T1.nii.gz');
% __________________________________
% Copyright 2015-2019 ExploreASL

    % ===========================================================
    %% Admin
    if nargin<2 || isempty(bReplace)
        bReplace = 0;
    end
    
    % get .nii (without .gz) name
    [Fpath, Ffile, Fext] = fileparts(PathIn);
    PathIn = fullfile(Fpath, [Ffile '.nii']);
    
    % Get the path to SPM defaced output file
    PathAnon = fullfile( Fpath, ['anon_' Ffile Fext]);
    xASL_delete(PathAnon);
    
    fprintf('%s\n','Defacing image for privacy reasons...');
    
    % ===========================================================
    %% Unzipping if needed
    xASL_io_ReadNifti(PathIn);
    
    % ===========================================================
    %% Run spm_deface
    matlabbatch{1}.spm.util.deface.images = {PathIn};
    spm_jobman('run',matlabbatch); close all;

    % ===========================================================
    %% Replace the original file with the defaced one
    if bReplace
        xASL_Move(PathAnon , PathIn, 1);
    end

end
