function xASL_wrp_LST_lga( x )
% input (1) T1 & (2) rFLAIR
% output (1) ModulateResample FLAIR (2) lession probability map ples_lga_[k]_rm[FLAIR].nii
% (3) mat file

clear matlabbatch

matlabbatch{1}.spm.tools.LST.lga.data_T1            = { x.P.Path_T1 };
matlabbatch{1}.spm.tools.LST.lga.data_F2            = { x.P.Path_rFLAIR };
InitialParm                                         = 0.3; % default
matlabbatch{1}.spm.tools.LST.lga.opts_lga.initial   = InitialParm; 
matlabbatch{1}.spm.tools.LST.lga.opts_lga.mrf       = 1; % default
matlabbatch{1}.spm.tools.LST.lga.html_report        = 0; % no HTML report generation, takes too long
matlabbatch{1}.spm.tools.LST.lga.xasl_quality       = x.Quality;

if  x.Quality
    matlabbatch{1}.spm.tools.LST.lga.opts_lga.maxiter   = 100; % default=50
else
    matlabbatch{1}.spm.tools.LST.lga.opts_lga.maxiter   = 3;
end

 fprintf('\n');    
 spm_jobman('run',matlabbatch); close all  

            
 % save output BIDS files 
 % (1) File: rmFLAIR 
 [pthf2, namf2, extf2] = fileparts(x.P.Path_rFLAIR); %namef2 =  sub-001_ses-ASL01_desc-resample_FLAIR
 Fileoutput = ['rm' namf2 extf2];
 cd(pthf2);
 if xASL_exist( Fileoutput,'file') %replace 'resample' with 'modulate' 
     C = strsplit(namf2,'_'); %get subID, sesID and suffix for BIDS file naming
     subID = strsplit(C{1},'-');
     sesID = strsplit(C{2},'-');
     suffix = C{4};
     new_mFLAIR = fullfile(pthf2, BIDS_FileName(subID{2}, sesID{2},'','','','ResampleModulate', suffix, 'nii'));
     old_mFLAIR = fullfile(pthf2, Fileoutput);
     copyfile(old_mFLAIR, new_mFLAIR);
    
     %delete msub-001 output of ps_LST_lpa.m
     xASL_delete(old_mFLAIR);
 end
 %(2) lession probability map ples_lga_[k]_rm[FLAIR].nii 
 File_ples = ['ples_lga_' k '_rm' namf2 extf2];
 if xASL_exist(File_ples)
     C = strsplit(namf2,'_'); %get subID, sesID and suffix for BIDS file naming
     subID = strsplit(C{1},'-');
     sesID = strsplit(C{2},'-');
     suffix = 'dseg';
     % label = 'L' Lession (common tissue classes)  atlas = lpa
     new_File_ples = fullfile(pthf2, BIDS_FileName(subID{2}, sesID{2},'ACPC','L','lga','', suffix, 'nii'));
     old_File_ples = fullfile(pthf2, File_ples);
     copyfile(old_File_ples, new_File_ples);
     copyfile( new_File_ples, x.P.rWMH_path);
     xASL_delete(old_File_ples);
 end
 % (3) mat-File: LST_lga_m[r][FLAIR].mat
 % contains all necessary
 % components that are needed for the longitudinal pipeline or for lesion filling.
 extf ='.mat';
 File_mat = ['LST_lga_m' namf2 extf];
 if xASL_exist(File_mat)
     C = strsplit(namf2,'_'); %get subID, sesID and suffix for BIDS file naming
     subID = strsplit(C{1},'-');
     sesID = strsplit(C{2},'-');
     suffix = 'dseg';
     % label = 'WMH' white matter hyperintensities  atlas = lga
     new_File_mat = fullfile(pthf2, BIDS_FileName(subID{2}, sesID{2},'ACPC','WMH','lga','', suffix, 'mat'));
     old_File_mat = fullfile(pthf2, File_mat);
     copyfile(old_File_mat, new_File_mat);
     xASL_delete(old_File_mat);
 end

end