function xASL_spm_deface( InPath )
%xASL_spm_deface Defaces the anatomical images, to preserve anonimity
% However, this makes it more difficult to check whether FLAIR & T1w images
% belong to the same subject

	fprintf('%s\n','Defacing image for privacy reasons...  ');

    [Fpath Ffile Fext]      = fileparts( InPath );
    anonName                = fullfile( Fpath, ['anon_' Ffile Fext]);
    xASL_delete(anonName);

    xASL_io_ReadNifti(InPath); % make sure to unzip
    matlabbatch{1}.spm.util.deface.images   = {InPath};
    spm_jobman('run',matlabbatch); close all;

    % Replace the original file with the defaced one
%     xASL_Move( anonName, InPath, 1);

end
