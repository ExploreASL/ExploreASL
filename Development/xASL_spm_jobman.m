function xASL_spm_jobman(x, matlabbatch)
%xASL_spm_jobman Wrapper around spm_jobman, allowing running it deployed as
%well

if isdeployed
    % First save matlabbatch in m-file
    PathJob = fullfile(pwd, 'Matlabbatch.mat');
    save(PathJob,'matlabbatch');
    system([x.SPMpath ' batch ' PathJob]);
    xASL_delete(PathJob);
else
    spm('defaults','PET');
    spm_jobman('run',matlabbatch);
    close all
end

end