function [x] = xASL_init_FileSystem(x)
%xASL_init_FileSystem Contains definition of path/filenames within studies,
% per BIDS/COBIDAS. Can be used for each module/wrapper/function.
% Any changes here will modify the complete pipeline


%% ------------------------------------------------------------------------------------
%% Admin


x.P     = struct; % (re-)initiate P (otherwise the looping below gets too much)
x.P.rWMH_path = '';



%% create function here that creates all required dirs

CreateDir(x.D.T1CheckDir); 
% Put this dir creation in separate scripts
% Create derivative dir as subfolder in each subject folder, where to put output,
% keeping original files
% /dartel -> /Population folder, 
% standardize naming creation
%CreateDir(x.D.Derivatives); %% BIDS output
CreateDir(x.D.TissueVolumeDir);
CreateDir(x.D.CoregDir);
CreateDir(x.D.FlowFieldCheck);

CreateDir(x.D.FLAIR_CheckDir);
CreateDir(x.D.FLAIR_REGDIR);

CreateDir(x.D.ASLCheckDir);
CreateDir(x.D.MotionDir);
CreateDir(x.D.RawEPIdir);

CreateDir(x.D.M0CheckDir);  
CreateDir(x.D.M0regASLdir);

%% ------------------------------------------------------------------------------------
%% General definitions
 %x.D.SubDer         = fullfile(   x.D.Derivatives, ['Sub-' x.P.SubjectID],'anat');
      %CreateDir( x.D.SubDer);
      
if  isfield(x,'SUBJECTDIR')
      [dummy{1} x.P.SubjectID dummy{2}]       = fileparts(x.SUBJECTDIR);
      C = strsplit(x.P.SubjectID,'-');
      x.P.SubjectID = C{2}; %C{1} = 'sub' C{2} = '001'
end

%x.D.SubDer         = fullfile(   x.D.Derivatives, ['Sub-' x.P.SubjectID],'anat');
     %CreateDir( x.D.SubDer);

if  isfield(x,'SESSIONDIR')
      [Fpath    x.P.SessionID dummy{2}]       = fileparts(x.SESSIONDIR);
      C = strsplit(x.P.SessionID,'-');
      x.P.SessionID = C{2}; %C{1} = 'ses' C{2} = 'ASL01'
      [dummy{1} x.P.SubjectID dummy{2}]       = fileparts(Fpath);
      C = strsplit(x.P.SubjectID,'-');
      x.P.SubjectID = C{2}; %C{1} = 'sub' C{2} = '001' 
end




%% ------------------------------------------------------------------------------------
%% Default prefixes
x.P.STRUCT      = 'T1';
% c1T1 = GM
% c2T1 = WM
% yT1  = deformation field to common space

FileDef     = {'FLAIR' 'T1' 'R1' 'ASL4D' 'fMRI' 'ASL4D'  'M0' 'TT'}; % FileTypes in raw BIDS anat folder
Label       = {'GM' 'WM' 'WMH' 'CBF' 'c1T1' 'c2T1' 'c3T1' 'yT1' 'PV_pGM' 'PV_pWM' 'y_ASL' 'CBF' 'qCBF' 'qCBF4D' 'qCBF' 'PseudoCBF' 'PWI' 'PWI4D' 'mean_PWI' 'mean_control' 'SD' 'SNR' 'SD_control' 'SNR_control' 'SliceGradient' 'FoV' };
Prefix      = {'Filled' 'SEGM' 'PreSEGM' 'r' 'm' 's' 'mr' 'rm' 'rmr' 'rr' 'temp' 'rtemp' 'rSEGM' 'mask' 'BiasField' 'noSmooth' 'untreated' 'despiked' 'Clipped' 'extrapolated'}; % r=resample m=modulate s=smooth w=warp q=quantified p=probability % USE t for TEMP? replace w by r
Prefix_desc = {'Filled' 'Segmented' 'PreExistingSegmented' 'Resample' 'modulate' 'smooth' 'ModulateResample' 'ResampleModulate' 'rmr' 'rr' 'temp' 'ResampleTemp' 'ResampleSegmented' 'mask' 'BiasField' 'noSmooth' 'untreated' 'despiked' 'Clipped' 'extrapolated'};

% use "b" for backup & "o" for original, to reduce the number of suffixes
% need to define various ASL4D_session files still


for iFD=1:length(FileDef)
        x.P.(FileDef{iFD})  = FileDef{iFD};
end
     

% Add prefixes
for iPref=1:length(Prefix)
    x.P.([Prefix{iPref} FileDef{iFD}]) = [Prefix{iPref} FileDef{iFD}];          
end        
        
%% define derivatives folder path    
 
%% x.P.derivatives = [];
%% ------------------------------------------------------------------------------------
%% Subject-specific prefixes
if  isfield(x.P,'SubjectID')
    
    for iSess=1:length(x.SESSIONS)
        x.P.SessionDir{iSess}   = fullfile(x.SUBJECTDIR,x.SESSIONS{iSess});
    end

    if ~isfield(x.P,'SessionID')
        x.P.SessionID   = 'ASL01';
    end

   %% ------------------------------------------------------------------------------------
   %% File definitions
   % all files are stored in one folder (session - anat Directory) according to
   % BIDS
   for iSess=1:length(x.SESSIONS)
        anat_Path{iSess} = fullfile( x.P.SessionDir{iSess} ,'anat');
   end
   
for iSess=1:length(x.SESSIONS)    
   for iFD=1:length(FileDef)  
         
            % Create file & path
            x.P.(['File_' FileDef{iFD}])     =                            BIDS_FileName(x.P.SubjectID,x.P.SessionID, '','','','',FileDef{iFD}, 'nii');
            x.P.(['Path_' FileDef{iFD}])     = fullfile(anat_Path{iSess}, BIDS_FileName(x.P.SubjectID,x.P.SessionID, '','','','',FileDef{iFD}, 'nii'));
            x.P.(['Pop_Path_' FileDef{iFD}]) = fullfile(x.D.PopDir,       BIDS_FileName(x.P.SubjectID,x.P.SessionID, '','','','',FileDef{iFD}, 'nii'));
            
            % Add prefixes
            for iPref=1:length(Prefix)
                
                x.P.(['File_' Prefix{iPref} FileDef{iFD}])     =                      BIDS_FileName(x.P.SubjectID,x.P.SessionID,'', '','',Prefix_desc{iPref},FileDef{iFD}, 'nii');
                x.P.(['Path_' Prefix{iPref} FileDef{iFD}])     = fullfile(anat_Path{iSess},  BIDS_FileName(x.P.SubjectID,x.P.SessionID,'', '','',Prefix_desc{iPref},FileDef{iFD}, 'nii'));
                x.P.(['Pop_Path_' Prefix{iPref} FileDef{iFD}]) = fullfile(x.D.PopDir, BIDS_FileName(x.P.SubjectID,x.P.SessionID,'', '','',Prefix_desc{iPref},FileDef{iFD}, 'nii'));
               
            end   
            % -> add all possible paths in 1 go, with always SubjectID in the filename   
            
   end
   % Add BIDS Label (Tissue parts) key-value  
    for ilabel=1:length(Label)
         x.P.(['File_' Label{ilabel} ])     =                      BIDS_FileName(x.P.SubjectID,x.P.SessionID,'',  Label{ilabel},'','','', 'nii');
         x.P.(['Path_' Label{ilabel} ])     = fullfile(anat_Path{iSess},  BIDS_FileName(x.P.SubjectID,x.P.SessionID,'', Label{ilabel},'','','', 'nii'));
         x.P.(['Pop_Path_'  Label{ilabel} ]) = fullfile(x.D.PopDir, BIDS_FileName(x.P.SubjectID,x.P.SessionID,'', Label{ilabel},'','','', 'nii')); 
     end
   
      % Add BIDS Label (Tissue parts) & Describtion key-value  
     for ilabel=1:length(Label)
       for iPref=1:length(Prefix)
         x.P.(['File_' Label{ilabel} '_' Prefix{iPref}])     =                      BIDS_FileName(x.P.SubjectID,x.P.SessionID,'',  Label{ilabel},'',Prefix_desc{iPref},'', 'nii');
         x.P.(['Path_' Label{ilabel} '_' Prefix{iPref}])     = fullfile(anat_Path{iSess},  BIDS_FileName(x.P.SubjectID,x.P.SessionID,'', Label{ilabel},'',Prefix_desc{iPref},'', 'nii'));
         x.P.(['Pop_Path_'  Label{ilabel} '_' Prefix{iPref}]) = fullfile(x.D.PopDir, BIDS_FileName(x.P.SubjectID,x.P.SessionID,'', Label{ilabel},'',Prefix_desc{iPref},'', 'nii')); 
       end
     end
end




%% ------------------------------------------------------------------------------------------
%% Add sidecars

FieldList   = fieldnames(x.P);
for iL=12:length(FieldList) % here we should remove the prefixes
    
    % Add .mat sidecars
    CarSuf1   = {'_parms_mat' '_sn_mat' '_mat'}; % suffixes for sidecars (dots replaced by underscores)
    CarSuf2   = {'_parms.mat' '_sn.mat' '.mat'}; % suffixes for sidecars
    for iS=1:length(CarSuf1)
        x.P.([FieldList{iL} CarSuf1{iS}]) = [x.P.(FieldList{iL})(1:end-4) CarSuf2{iS}];
    end
end    
    
    


% -> add all possible paths in 1 go, with always SubjectID in the filename
% -> create transformation2BIDS script
        
% -> remove FileNames, PathNames are sufficient
% % rename c1 c2 -> pGM pWM
% 

 % remove some later, pruning, when not necessary anymore
% 


end