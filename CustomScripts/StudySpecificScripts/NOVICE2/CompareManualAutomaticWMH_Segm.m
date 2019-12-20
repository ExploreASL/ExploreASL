%% Get WMH volume from NOVICE

load('C:\BackupWork\ASL\Novice\ComparisonSubjects\ComparisonSubjects.mat','-mat');

RDir='C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_NOVICE';

clear List
for iSubject=1:length(ComparingSubject)
    xASL_TrackProgress(iSubject,length(ComparingSubject));
    SubjectID = [ComparingSubject{iSubject,1} '_1'];
    WMHfile = fullfile(RDir,SubjectID,'WMH_SEGM.nii');
    List{iSubject,1} = SubjectID;

    if xASL_exist(WMHfile,'file')
        nii = xASL_io_ReadNifti(WMHfile);
        im = xASL_io_Nifti2Im(WMHfile);
        dim = nii.hdr.pixdim(2:4);
        vol = xASL_stat_SumNan(im(:));

        vol = prod(dim).*vol;
        List{iSubject,2} = vol;
    else
        List{iSubject,2} = 'n/a';
    end
end
