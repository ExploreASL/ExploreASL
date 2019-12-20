function xASL_wrp_LST_lesfill( x )

%xASL_wrp_LST_lesfill ExploreASL wrapper for LST Lesion filling function in BIDS format 
% input is two images: 1- image to be filled in native space(T1w) 2-
% probability lesion map obtained by obtained either by the LGA, LPA or
%longitudinal pipeline)
% it outputs 1 file: This procedure saves the filled image in native space 
%([image]_filled_lga_[k]_rm[FLAIR].nii for LGA or [image]_filled_lpa_m[r][FLAIR].nii
       clear matlabbatch
            matlabbatch{1}.spm.tools.LST.filling.data           = {x.P.Path_T1};
            matlabbatch{1}.spm.tools.LST.filling.data_plm       = {rWMHPath};
            matlabbatch{1}.spm.tools.LST.filling.html_report    = 0; % saves time
            spm_jobman('run',matlabbatch); close all
            
     %rename output T1 filled according to BIDS specification       
     %file name: sub-001_ses-ASL01_space-ACPC_label-WMH_desc-filled_T1w.nii
     [pthf2, namf2, extf2] = fileparts(x.P.Path_T1); %namef2 =  sub-001_ses-ASL01_T1
     Fileoutput = x.P.File_T1;
     cd(pthf2);
  if xASL_exist( Fileoutput,'file') %rename T1 to T1_filled
     C = strsplit(namf2,'_'); %get subID, sesID and suffix for BIDS file naming
     subID = strsplit(C{1},'-');
     sesID = strsplit(C{2},'-');
     suffix = C{3};
     new_mFLAIR = fullfile(pthf2, BIDS_FileName(subID{2}, sesID{2},'ACPC','WMH','','filled', suffix, 'nii'));
     old_mFLAIR = fullfile(pthf2, Fileoutput);
     copyfile(old_mFLAIR, new_mFLAIR);
    
     %delete msub-001 output of ps_LST_lpa.m
     xASL_delete(old_mFLAIR);
 end
 
end