function new_File_ples = xASL_wrp_LST_lpa( x )
new_File_ples ='';
%xASL_wrp_LST_lpa ExploreASL wrapper for LST lpa function in BIDS format 
% input is FLAIR image
% it outputs 2 files (1) coregistered mFLAIR (2) lesson probability map (ples_lpa_m[r][FLAIR].nii)
% (3) mat file

 clear matlabbatch
 fprintf('\n%s\n','----------------------------------------');      
 fprintf('%s\n','WMH segmentation performed using LST LPA')
 % Lesion Prediction Algorithm (LPA)
 matlabbatch{1}.spm.tools.LST.lpa.data_F2        = { x.P.Path_rFLAIR }; % FLAIR resampled to T1 native space
 matlabbatch{1}.spm.tools.LST.lpa.data_coreg     = {''};
 matlabbatch{1}.spm.tools.LST.lpa.html_report    = 0; % no HTML report output (takes a long time)
 matlabbatch{1}.spm.tools.LST.lpa.xasl_quality   = x.Quality;
 fprintf('\n');    
 spm_jobman('run',matlabbatch); close all  
 

 % save output BIDS files 
 % (1) File: mFLAIR 
 [pthf2, namf2, extf2] = fileparts(x.P.Path_rFLAIR); %namef2 =  sub-001_ses-ASL01_desc-resample_FLAIR
 Fileoutput = ['m' namf2 extf2];
 cd(pthf2);
 if xASL_exist( Fileoutput,'file') %replace 'resample' with 'modulate' 
     C = strsplit(namf2,'_'); %get subID, sesID and suffix for BIDS file naming
     subID = strsplit(C{1},'-');
     sesID = strsplit(C{2},'-');
     suffix = C{4};
     new_mFLAIR = fullfile(pthf2, BIDS_FileName(subID{2}, sesID{2},'','','','modulate', suffix, 'nii'));
     old_mFLAIR = fullfile(pthf2, Fileoutput);
     copyfile(old_mFLAIR, new_mFLAIR);
     
    
     %delete msub-001 output of ps_LST_lpa.m
     xASL_delete(old_mFLAIR);
 end
 % (2) File: lession prediction map (ples_lpa_mFLAIR) =rWMHPath
 File_ples = ['ples_lpa_m' namf2 extf2];
 if xASL_exist(File_ples)
     C = strsplit(namf2,'_'); %get subID, sesID and suffix for BIDS file naming
     subID = strsplit(C{1},'-');
     sesID = strsplit(C{2},'-');
     suffix = 'dseg';
     % label = 'WMH' white matter hyperintensities,  atlas = lpa
     new_File_ples = fullfile(pthf2, BIDS_FileName(subID{2}, sesID{2},'ACPC','WMH','lpa','', suffix, 'nii'));
     old_File_ples = fullfile(pthf2, File_ples);
     copyfile(old_File_ples, new_File_ples);
     %copyfile( new_File_ples, x.P.rWMH_path);
     xASL_delete(old_File_ples);
 end
 % (3) mat-File: LST_lpa_m[r][FLAIR].mat
 % contains all necessary
 % components that are needed for the longitudinal pipeline or for lesion filling.
     extf ='.mat';
     File_mat = ['LST_lpa_m' namf2 extf];
     C = strsplit(namf2,'_'); %get subID, sesID and suffix for BIDS file naming
     subID = strsplit(C{1},'-');
     sesID = strsplit(C{2},'-');
     suffix = 'dseg';
     % label = 'WMH' white matter hyperintensities  atlas = lpa
     new_File_mat = fullfile(pthf2, BIDS_FileName(subID{2}, sesID{2},'ACPC','WMH','lpa','', suffix, 'mat'));
     old_File_mat = fullfile(pthf2, File_mat);
     copyfile(old_File_mat, new_File_mat);
     xASL_delete(old_File_mat);


