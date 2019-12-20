function cat_run_job1070(job,tpm,subj)
% run CAT 
% ______________________________________________________________________
%
% Initialization functions for the CAT preprocessing
%  * creation of the subfolder structure (if active)
%  * check of image resolution (avoid scans with very low resolution)
%  * interpolation 
%  * affine preprocessing (APP)
%    >> cat_run_job_APP_init
%    >> cat_run_job_APP_final
%  * affine registration
%  * initial SPM preprocessing
%
%   cat_run_job(job,tpm,subj)
% 
%   job  .. SPM job structure with main parameter
%   tpm  .. tissue probability map (hdr structure)
%   subj .. file name
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_run_job1070.m 1348 2018-08-02 13:01:57Z gaser $

%#ok<*WNOFF,*WNON>

    % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
    dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

    global cat_err_res; % for CAT error report

    stime  = clock; 
    stime0 = stime; % overall processing time

    % create subfolders if not exist
    pth = spm_fileparts(job.channel(1).vols{subj}); 
    if job.extopts.subfolders
    
      folders = {'mri','report'};
      for i=1:numel(folders)
        if ~exist(fullfile(pth,folders{i}),'dir')
          mkdir(fullfile(pth,folders{i}));
        end
      end
    
      if ~exist(fullfile(pth,'surf'),'dir') && job.output.surface
        mkdir(fullfile(pth,'surf'));
      end
    
      if ~exist(fullfile(pth,'label'),'dir') && job.output.ROI
        mkdir(fullfile(pth,'label'));
      end
      
      mrifolder    = 'mri';
      reportfolder = 'report';
    else
      mrifolder    = '';
      reportfolder = '';
    end
    
    % create subject-wise diagy file with the command-line output
    [pp,ff,ee,ex] = spm_fileparts(job.data{subj}); 
    catlog = fullfile(pth,reportfolder,['catlog_' ff '.txt']);
	warning('off','MATLAB:DELETE:Permission'); % ExploreASL hack
    if exist(catlog,'file'), delete(catlog); end % write every time a new file, turn this of to have an additional log file
	warning('on','MATLAB:DELETE:Permission'); % ExploreASL hack
    % diary(catlog); % ExploreASL fix - disable CAT12 outputs as we write everything to ExploreASL log file.
    
    % print current CAT release number and subject file
    [n,r] = cat_version;
    str  = sprintf('%s r%s: %d/%d',n,r,subj,numel(job.channel(1).vols));
    str2 = spm_str_manip(job.channel(1).vols{subj}(1:end-2),['a' num2str(70 - length(str))]);
    cat_io_cprintf([0.2 0.2 0.8],'\n%s\n%s: %s%s\n%s\n',...
          repmat('-',1,72),str,...
          repmat(' ',1,70 - length(str) - length(str2)),str2,...
          repmat('-',1,72));
    clear r str str2
    
    %  -----------------------------------------------------------------
    %  separation of full CAT preprocessing and SPM segmentation
    %  preprocessing (running DARTEL and PBT with SPM segmentation)
    %  -----------------------------------------------------------------
    [pp,ff,ee,ex] = spm_fileparts(job.data{subj}); 
    if exist(fullfile(pp,['c1' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,['c2' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,['c3' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,[ff(3:end) '_seg8.mat']),'file');
       
        job.data{subj}          = fullfile(pp,[ff ee]); 
        job.channel.vols{subj}  = fullfile(pp,[ff ee]); 

        % prepare SPM preprocessing structure 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
          images = char(images,job.channel(n).vols{subj});
        end

        obj.image    = spm_vol(images);
        obj.fwhm     = job.opts.fwhm;
        obj.biasreg  = cat(1,job.opts.biasreg);
        obj.biasfwhm = cat(1,job.opts.biasfwhm);
        obj.tpm      = tpm;
        obj.lkp      = [];
        spm_check_orientations(obj.image);
        
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end
        
        obj.reg      = job.opts.warpreg;
        obj.samp     = job.opts.samp;              
        cfname  = fullfile(pp,[ff ee]);
        ofname  = fullfile(pp,[ff(3:end) ee]); 
        nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
        copyfile(ofname,nfname); 

        Ysrc0    = single(spm_read_vols(obj.image)); 
        Ylesion  = single(isnan(Ysrc0) | isinf(Ysrc0) | Ysrc0==0); clear Ysrc0;
        
        res = load(fullfile(pp,[ff(3:end) '_seg8.mat']));
        job.channel(1).vols{subj}  = [nfname ex];
        job.channel(1).vols0{subj} = [ofname ex];
        res.image  = spm_vol([nfname ex]);
        res.image0 = spm_vol([ofname ex]);
        res.imagec = spm_vol([cfname ex]);
        res.spmpp  = 1; 
    else

        %  -----------------------------------------------------------------
        %  check resolution properties
        %  -----------------------------------------------------------------
        %  There were some images that should not be processed. So we have  
        %  to check for large slice thickness and low spatial resolution.
        %  -----------------------------------------------------------------
        for n=1:numel(job.channel) 
          V = spm_vol(job.channel(n).vols{subj});
          vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));

          if any(vx_vol>5)  % too thin slices
            error('CAT:cat_main:TooLowResolution', sprintf(...
                 ['Voxel resolution has to be better than 5 mm in any dimension \n' ...
                  'for reliable CAT preprocessing! \n' ...
                  'This image has a resolution %0.2fx%0.2fx%0.2f mm%s. '], ... 
                    vx_vol,char(179))); %#ok<SPERR>
          end
          if prod(vx_vol)>27  % too small voxel volume (smaller than 3x3x3 mm3)
            error('CAT:cat_main:TooHighVoxelVolume', ...
                 ['Voxel volume has to be smaller than 10 mm%s (around 3x3x3 mm%s) to \n' ...
                  'allow a reliable CAT preprocessing! \n' ...
                  'This image has a voxel volume of %0.2f mm%s. '], ...
                  char(179),char(179),prod(vx_vol),char(179));
          end
          if max(vx_vol)/min(vx_vol)>8 % anisotropy 
            error('CAT:cat_main:TooStrongIsotropy', sprintf(...
                 ['Voxel isotropy (max(vx_size)/min(vx_size)) has to be smaller than 8 to \n' ...
                  'allow a reliable CAT preprocessing! \n' ...
                  'This image has a resolution %0.2fx%0.2fx%0.2f mm%s and a isotropy of %0.2f. '], ...
                  vx_vol,char(179),max(vx_vol)/min(vx_vol))); %#ok<SPERR>
          end
        end

        % save original file name 
        for n=1:numel(job.channel) 
          job.channel(n).vols0{subj} = job.channel(n).vols{subj};
        end

        % always create the n*.nii image because of the real masking of the
        % T1 data for spm_preproc8 that include rewriting the image!
        for n=1:numel(job.channel) 
		
			fprintf('Initializing CAT pre-processing');%%% ExploreASL fix
		
          [pp,ff,ee] = spm_fileparts(job.channel(n).vols{subj}); 
          ofname  = fullfile(pp,[ff ee]); 
          nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
          if strcmp(ee,'.nii')
            copyfile(ofname,nfname); 
          elseif strcmp(ee,'.img')
            V = spm_vol(job.channel(n).vols{subj});
            Y = spm_read_vols(V);
            V.fname = nfname; 
            spm_write_vol(V,Y);
            clear Y; 
          end
          job.channel(n).vols{subj} = nfname;

          %% denoising
          if job.extopts.NCstr~=0
            NCstr.labels = {'none','full','light','medium','strong','heavy'};
            NCstr.values = {0 1 2 -inf 4 5}; 
            stime = cat_io_cmd(sprintf('SANLM denoising (%s)',...
              NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')}));
            cat_vol_sanlm(struct('data',nfname,'verb',0,'prefix','','NCstr',job.extopts.NCstr)); 
            fprintf('%5.0fs\n',etime(clock,stime));   
          end

          %% skull-stripping detection
          %  ------------------------------------------------------------
          %  Detect skull-stripping or defaceing because it strongly 
          %  affects SPM segmentation that expects gaussian distribution! 
          %  If a brain mask was used than we expect 
          %   - many zeros (50% for small background - 80-90% for large backgrounds)
          %   - a smaller volume because of missing skull (below 2500 cm3)
          %   - only one object (the masked regions)
          %   - only one background (not in every case?)
          %   - less variance of tissue intensity (only 3 brain classes)
          %  ------------------------------------------------------------
		  disp('Skull stripping detection');%%% ExploreASL fix
          VFn   = spm_vol(nfname); 
          YF    = spm_read_vols(VFn); 
          Oth   = cat_stat_nanmean(YF(YF(:)~=0 & YF(:)>cat_stat_nanmean(YF(:)))); 
          F0vol = cat_stat_nansum(YF(:)~=0) * prod(vx_vol) / 1000; 
          F0std = cat_stat_nanstd(YF(YF(:)>0.5*Oth & YF(:)>0)/Oth); 
          YFC = YF~=0; 
          if sum(YFC(:)>0)<numel(YFC)*0.9 && sum(YFC(:)>0)>numel(YFC)*0.1  % if there is a meanful background
            YFC = ~cat_vol_morph(YF~=0,'lc',1);                            % close noisy background
          end
          [YL,numo] = spm_bwlabel(double(YF~=0),26);  clear YL;            % number of objects
          [YL,numi] = spm_bwlabel(double(YFC==0),26); clear YL;            % number of background regions 
          skullstrippedpara = [sum(YF(:)==0)/numel(YF) numo numi F0vol F0std]; 
          skullstripped = ...
            skullstrippedpara(1)>0.5 && ...                     % many zeros
            skullstrippedpara(2)<5  && ...                      % only a few objects
            skullstrippedpara(3)<10 && ...                      % only a few background regions 
            F0vol<2500 && F0std<0.5;                                       % many zeros and not too big
          skullstripped = skullstripped || ...
            sum([skullstrippedpara(1)>0.8 F0vol<1500 F0std<0.4])>1; % or 2 extrem values
          clear YF YFC F0vol F0std numo numi; 

          %% Interpolation
          %  -----------------------------------------------------------------
          %  The interpolation can help reducing problems for morphological
          %  operations for low resolutions and strong isotropic images. 
          %  Especially for Dartel registration a native resolution larger than the Dartel 
          %  resolution helps to reduce normalization artifacts of the
          %  deformations. Furthermore, even if artifacts can be reduced by the final smoothing
          %  it is much better to avoid them.  

          % prepare header of resampled volume
          Vi        = spm_vol(job.channel(n).vols{subj}); 
          vx_vol    = sqrt(sum(Vi.mat(1:3,1:3).^2));
          vx_vol    = round(vx_vol*10^2)/10^2; % avoid small differences 

          % we have to look for the name of the field due to the GUI job struct generation! 
          restype   = char(fieldnames(job.extopts.restypes));
          switch restype
            case 'native'
              vx_voli  = vx_vol;
            case 'fixed', 
              vx_voli  = min(vx_vol ,job.extopts.restypes.(restype)(1) ./ ...
                         ((vx_vol > (job.extopts.restypes.(restype)(1)+job.extopts.restypes.(restype)(2)))+eps));
              vx_voli  = max(vx_voli,job.extopts.restypes.(restype)(1) .* ...
                         ( vx_vol < (job.extopts.restypes.(restype)(1)-job.extopts.restypes.(restype)(2))));
            case 'best'
              best_vx  = max( min(vx_vol) ,job.extopts.restypes.(restype)(1)); 
              vx_voli  = min(vx_vol ,best_vx ./ ((vx_vol > (best_vx + job.extopts.restypes.(restype)(2)))+eps));
            otherwise 
              error('cat_run_job:restype','Unknown resolution type ''%s''. Choose between ''fixed'',''native'', and ''best''.',restype)
          end

          % interpolation 
          if any( (vx_vol ~= vx_voli) )  
            stime = cat_io_cmd(sprintf('Internal resampling (%4.2fx%4.2fx%4.2fmm > %4.2fx%4.2fx%4.2fmm)',vx_vol,vx_voli));
           
            Vi        = rmfield(Vi,'private');
            imat      = spm_imatrix(Vi.mat); 
            Vi.dim    = round(Vi.dim .* vx_vol./vx_voli);
            imat(7:9) = vx_voli .* sign(imat(7:9));
            Vi.mat    = spm_matrix(imat);

            Vn = spm_vol(job.channel(n).vols{subj}); 
            Vn = rmfield(Vn,'private'); 
            cat_vol_imcalc(Vn,Vi,'i1',struct('interp',2,'verb',0));
            vx_vol = vx_voli;
          
            fprintf('%5.0fs\n',etime(clock,stime));     
          else
            vx_vol = sqrt(sum(Vi.mat(1:3,1:3).^2));
          end

          clear Vi Vn;
        end

        %  prepare SPM preprocessing structure 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
            images = char(images,job.channel(n).vols{subj});
        end

        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);

        obj.fwhm     = job.opts.fwhm;
        obj.biasreg  = cat(1,job.opts.biasreg);
        obj.biasfwhm = cat(1,job.opts.biasfwhm);
        obj.tpm      = tpm;        
        obj.reg      = job.opts.warpreg;
        obj.samp     = job.opts.samp;              
        obj.lkp      = [];
        
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end

        spm_check_orientations(obj.image);

        %% Initial affine registration.
        %  -----------------------------------------------------------------
        %  APP option with subparameter
        %  Skull-stripping is helpful for correcting affine registration of neonates and other species. 
        %  Bias correction is important for the affine registration.
        %  However, the first registation can fail and further control is required   
        % 
        %  bias = 0-5 = none, light, light threshold, light apply, fine apply (light=only for registration) 
        %  msk  = 0-4 = none, head msk, head hard, brain msk, brain hard (msk=mask only, hard=remove nonmsk)
        %  aff  = 0-1 = no affreg, affreg

        if ~strcmp(job.extopts.species,'human'), job.extopts.APP='nonhuman'; end

        switch job.extopts.APP
          case {0,'none'},     app.bias=0; app.msk=0; app.aff=1; % old default
          case {1,'light'},    app.bias=1; app.msk=1; app.aff=1; % affreg with BC; thresholding and head masking for SPM
          case {2,'medium'},   app.bias=2; app.msk=1; app.aff=1; % no-affreg; BC and head masking for SPM  
          case {3,'strong'},   app.bias=2; app.msk=1; app.aff=0; % no-affreg; BC and head masking for SPM  
          case {4,'heavy'},    app.bias=4; app.msk=3; app.aff=0; % no-affreg; BC and brain masking for SPM  
          case {5,'nonhuman'}, app.bias=4; app.msk=4; app.aff=0; % no-affreg; BC and brain masking for SPM  
          case {1070},         app.bias=1; app.msk=1; app.aff=1; % affreg with BC; thresholding and head masking for SPM
          otherwise
        end

        Affine  = eye(4);
        [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
        Pbt = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
        Pb  = char(job.extopts.brainmask);
        Pt1 = char(job.extopts.T1);
        
        if ~isempty(job.opts.affreg)      

            %% first affine registration (with APP)
            try 
                VG = spm_vol(Pt1);
            catch
                pause(rand(1))
                VG = spm_vol(Pt1);
            end
            VF = spm_vol(obj.image(1));

            % skull-stripping of the template
            if skullstripped
              % print a warning for all users because processing of
              % skull-stripped data is not standard!
              if job.extopts.APP==1 || job.extopts.APP==2
                cat_io_cprintf('warn',[...
                  'WARNING: Detected skull-stripped or strongly masked image. Skip APP. \n' ...
                  '         Use skull-stripped initial affine registration template and  \n' ...
                  '         TPM without head tissues (class 4 and 5)! \n']);
              else
                cat_io_cprintf('warn',[...
                  'WARNING: Detected skull-stripped or strongly masked image. \n' ...
                  '         Use skull-stripped initial affine registration template and \n' ...
                  '         TPM without head tissues (class 4 and 5)! \n']);
              end
              if job.extopts.verb>1
                cat_io_cprintf('warn',sprintf(...
                 ['           %0.2f%%%% zeros, %d object(s), %d background region(s) \n' ...
                  '           %4.0f cm%s, normalized SD of all tissues %0.2f \n'],...
                  skullstrippedpara(1:4),char(179),skullstrippedpara(5))); 
              end
  
              % skull-stripping of the template
              VB = spm_vol(Pb);
              [VB2,YB] = cat_vol_imcalc([VG,VB],Pbt,'i1 .* i2',struct('interp',3,'verb',0)); 
              VB2.dat(:,:,:) = eval(sprintf('%s(YB/max(YB(:))*255);',spm_type(VB2.dt))); 
              VB2.pinfo      = repmat([1;0],1,size(YB,3));
              VG             = cat_spm_smoothto8bit(VB2,0.5);
              clear VB2 YB; 
            end

            % Rescale images so that globals are better conditioned
            VF.pinfo(1:2,:) = VF.pinfo(1:2,:)/spm_global(VF);
            VG.pinfo(1:2,:) = VG.pinfo(1:2,:)/spm_global(VG);

            % APP step 1 rough bias correction 
            % --------------------------------------------------------------
            % Already for the rough initial affine registration a simple  
            % bias corrected and intensity scaled image is required, because
            % large head intensities can disturb the whole process.
            % --------------------------------------------------------------
            % ds('l2','',vx_vol,Ym, Yt + 2*Ybg,obj.image.private.dat(:,:,:)/WMth,Ym,60)
            if  app.bias ||  app.msk  
                stime = cat_io_cmd('APP: Rough bias correction'); 
                [Ym,Yt,Ybg,WMth,bias] = cat_run_job_APP_init1070(single(obj.image.private.dat(:,:,:)),vx_vol,job.extopts.verb);

                stime = cat_io_cmd('Coarse affine registration','','',1,stime); 
                
                % write data to VF
                VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
                VF.dat(:,:,:) = cat_vol_ctype(Ym * 200,'uint8'); 
                VF.pinfo      = repmat([1;0],1,size(Ym,3));
                clear WI; 

                % smoothing
                resa  = obj.samp*2; % definine smoothing by sample size
                VF1   = spm_smoothto8bit(VF,resa);
                VG1   = spm_smoothto8bit(VG,resa);

            else
                % standard approach with static resa value and no VG smoothing
                stime = cat_io_cmd('Coarse affine registration'); 
                resa  = 8;
                VF1   = spm_smoothto8bit(VF,resa);
                VG1   = VG; 
            end

            % prepare affine parameter 
            aflags     = struct('sep',obj.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1); 
            aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
            aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

            %% affine registration
            try
                spm_plot_convergence('Init','Coarse affine registration','Mean squared difference','Iteration');
            catch
                spm_chi2_plot('Init','Coarse affine registration','Mean squared difference','Iteration');
            end
            
            if app.aff
                warning off 
                try 
                  [Affine0, affscale]  = spm_affreg(VG1, VF1, aflags, eye(4)); Affine = Affine0; 
                catch
                  affscale = 0; 
                end
                if affscale>3 || affscale<0.5
                  stime  = cat_io_cmd('Coarse affine registration failed. Try fine affine registration.','','',1,stime);
                  Affine = eye(4); 
                end
                warning on
            else
              Affine = eye(4); affscale = 1;
            end

            %% APP step 2 - brainmasking and second tissue separated bias correction  
            %  ---------------------------------------------------------
            %  The second part of APP maps a brainmask to native space and 
            %  refines it by morphologic operations and region-growing to
            %  adapt for worse initial affine alignments. It is important
            %  that the mask covers the whole brain, whereas additional
            %  masked head is here less problematic.
            %  ---------------------------------------------------------
            %    ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,90)
            aflags.sep = obj.samp/2; 
            aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
            aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));
            if app.bias>2 || app.msk>2 
                %% apply (first affine) registration on the default brain mask
                VFa = VF; 
                if app.aff, VFa.mat = Affine * VF.mat; else Affine = eye(4); affscale = 1; end
                if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
                [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); Yb = Yb>0.5 & ~Ybg; 

                stime = cat_io_cmd('APP: Fine bias correction and skull-stripping','','',1,stime); 

                % fine APP
                [Ym,Yp0,Yb] = cat_run_job_APP_final(single(obj.image.private.dat(:,:,:)),...
                    Ym,Yb,Ybg,vx_vol,job.extopts.gcutstr,job.extopts.verb);
                stime = cat_io_cmd('Affine registration','','',1,stime); 

                %% smooth data
                VF.dat(:,:,:) =  cat_vol_ctype(Ym*200); 
                VF1 = spm_smoothto8bit(VF,aflags.sep);
                VG1 = spm_smoothto8bit(VG,aflags.sep);

                if 1 % brain masking for affine registration 
                  VB  = spm_vol(Pb);
                  Ybt = spm_read_vols(VB); 
                  VG1.dat(:,:,:) =  cat_vol_ctype(single(VG1.dat(:,:,:)) .* smooth3(Ybt));
                  VF1.dat(:,:,:) =  cat_vol_ctype(single(VF1.dat(:,:,:)) .* smooth3(Yb));
                end   

            elseif app.bias || app.msk 
                % smooth data
                stime = cat_io_cmd('Affine registration','','',1,stime); 
                VF.dat(:,:,:) =  cat_vol_ctype(Ym*200); 
                VF1 = spm_smoothto8bit(VF,aflags.sep);
                VG1 = spm_smoothto8bit(VG,aflags.sep);
            else
                % standard approach 
                stime = cat_io_cmd('Affine registration','','',1,stime); 
                VF1 = spm_smoothto8bit(VF,aflags.sep);
                VG1 = spm_smoothto8bit(VG,0.5); 
            end

            %% fine affine registration 
            if app.aff 
              try
                  spm_plot_convergence('Init','Affine registration','Mean squared difference','Iteration');
              catch
                  spm_chi2_plot('Init','Affine registration','Mean squared difference','Iteration');
              end
              warning off
              [Affine1,affscale1] = spm_affreg(VG1, VF1, aflags, Affine, affscale);  
              warning on
              if ~any(any(isnan(Affine1(1:3,:)))) && affscale>0.5 && affscale<3, Affine = Affine1; end
            end
            clear VG1 VF1
        else
          [Ym,Yt,Ybg,WMth,bias] = cat_run_job_APP_init1070(single(obj.image.private.dat(:,:,:)),vx_vol,job.extopts.verb);
        end

        
        
        %% Lesion masking as zero values of the orignal image (2018-06):
        %  We do not use NaN and -INF because (i) most images are only (u)int16
        %  and do not allow such values, (ii) NaN can be part of the background
        %  of resliced images, and (iii) multiple options are not required here. 
        %  Zero values can also occure by poor data scaling or processing in the 
        %  background but also by other (large) CSF regions and we have to remove  
        %  these regions later. 
        %  We further discussed to use a separate mask images but finally desided
        %  to keep this as simple as possible using no additional options!
        obj.image0 = spm_vol(job.channel(1).vols0{subj});
		% ExploreASL fix
        % Ysrc0      = spm_read_vols(obj.image0); 
        % Ylesion    = single(Ysrc0==0); clear Ysrc0; 
		if isfield(job.extopts,'xasl_lesion') && ~isempty(job.extopts.xasl_lesion) && (length(job.extopts.xasl_lesion{1}) > 1)
			Ysrc0 = spm_vol(job.extopts.xasl_lesion{1});
			Ysrc0 = spm_read_vols(Ysrc0);
			Ylesion    = single(Ysrc0>0.5); clear Ysrc0;
		else
			Ysrc0      = spm_read_vols(obj.image0);
			Ylesion    = single(Ysrc0==0); clear Ysrc0;
		end
		
        Ylesion(smooth3(Ylesion)<0.5)=0; % general denoising 
        if any( obj.image0.dim ~= obj.image.dim )
          mat      = obj.image0.mat \ obj.image.mat;
          Ylesion  = smooth3(Ylesion); 
          Ylesionr = zeros(obj.image.dim,'single'); 
          for i=1:obj.image.dim(3),
            Ylesionr(:,:,i) = single(spm_slice_vol(Ylesion,mat*spm_matrix([0 0 i]),obj.image.dim(1:2),[1,NaN]));
          end
          Ylesion = Ylesionr>0.5; clear Ylesionr;
        end
        if exist('Ybg','var'), Ylesion(Ybg)=0; end % denoising in background
        
        
        
        %% APP for spm_maff8
        %  optimize intensity range
        %  we have to rewrite the image, because SPM reads it again 
        if job.extopts.APP>0
            % WM threshold
            Ysrc = single(obj.image.private.dat(:,:,:)); 
            Ysrc(isnan(Ysrc) | isinf(Ysrc)) = min(Ysrc(:));


            if exist('Yb','var')
                Yb = cat_vol_morph(cat_vol_morph(Yb,'d',1),'lc',1);
                th = cat_stat_nanmean(Ysrc(Yb(:) & Ysrc(:)>cat_stat_nanmean(Ysrc(Yb(:))))) / ...
                     cat_stat_nanmean(Ym(Yb(:)   & Ym(:)>cat_stat_nanmean(Ym(Yb(:)))));
                if exist('WMth','var'), th = max(th,WMth); end
            else % only initial bias correction
                th = WMth;
            end
            bth = min( [ mean(single(Ysrc( Ybg(:)))) - 2*std(single(Ysrc( Ybg(:)))) , ...
                         mean(single(Ysrc(~Ybg(:)))) - 4*std(single(Ysrc(~Ybg(:)))) , ...
                         min(single(Ysrc(~Ybg(:))))]); 

            % add and write bias corrected (skull-stripped) image
            if app.bias<=1 % app.bias=1 is just a simple bias correction for affreg and maybe causes errors in the BWP cerebellum if used further! 
                Ymc = Ysrc; 
            else
                Ymc = Ym * abs(diff([bth,th])) + bth; % use the bias corrected image
            end

            % set variable and write image
            obj.image.dat(:,:,:)         = Ymc;  
            obj.image.pinfo              = repmat([255;0],1,size(Ymc,3));
            obj.image.private.dat(:,:,:) = Ymc; 

            obj.image.dt    = [spm_type('FLOAT32') spm_platform('bigend')];
            obj.image.pinfo = repmat([1;0],1,size(Ysrc,3));
            clear Ysrc; 
        end

        %  Fine Affine Registration with 3 mm sampling distance
        %  This does not work for non human data (or very small brains)
        stime = cat_io_cmd('SPM preprocessing 1 (estimate 1):','','',1,stime);
        if strcmp('human',job.extopts.species) 
            spm_plot_convergence('Init','Fine affine registration','Mean squared difference','Iteration');
            warning off 
            Affine2 = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.tpm,Affine ,job.opts.affreg,20); 
            Affine3 = spm_maff8(obj.image(1),obj.samp,obj.fwhm,       obj.tpm,Affine2,job.opts.affreg,20);
            warning on  
            if ~any(any(isnan(Affine3(1:3,:)))), Affine = Affine3; end
        end
        obj.Affine = Affine;

        % set original non-bias corrected image
        if job.extopts.APP==1
          obj.image = spm_vol(images);
        end
        cat_err_res.obj = obj; 

       
      

        if skullstripped 
          %% update number of SPM gaussian classes 
          Ybg = 1 - spm_read_vols(obj.tpm.V(1)) - spm_read_vols(obj.tpm.V(2)) - spm_read_vols(obj.tpm.V(3));
          obj.tpm.V(4).dat = Ybg; 
          obj.tpm.dat{4}   = Ybg; 
          obj.tpm.V(4).pinfo = repmat([1;0],1,size(Ybg,3));
          obj.tpm.dat(5:6) = []; 
          obj.tpm.V(5:6)   = []; 
          obj.tpm.bg1(4)   = obj.tpm.bg1(6);
          obj.tpm.bg2(4)   = obj.tpm.bg1(6);
          obj.tpm.bg1(5:6) = [];
          obj.tpm.bg2(5:6) = [];
            
          job.opts.ngaus = 3*ones(4,1); % this is more safe
          obj.lkp        = [];
          for k=1:numel(job.opts.ngaus)
            job.tissue(k).ngaus = job.opts.ngaus(k);
            obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
          end
        end

        %% SPM preprocessing 1
        %  ds('l2','a',0.5,Ym,Ybg,Ym,Ym,140);
        %  ds('l2','a',0.5,Ysrc/WMth,Yb,Ysrc/WMth,Yb,140);
        warning off 
        try 
          % inital estimate
          stime = cat_io_cmd('SPM preprocessing 1 (estimate 2):','','',job.extopts.verb-1,stime);
   
          if job.opts.redspmres==0 
            res = spm_preproc8(obj);
          else
            image1 = obj.image; 
            [obj.image,redspmres]  = cat_vol_resize(obj.image,'interpv',1);
            res = spm_preproc8(obj);
            res.image1 = image1; 
            clear reduce; 
          end
            
          % for non-skull-stripped brains use masked brains to get better estimates
          % esp. for brains with thinner skull or special defacing
          if ~skullstripped
            stime = cat_io_cmd('SPM preprocessing 1 (estimate skull-stripped):','','',job.extopts.verb-1,stime);
            % use dilated mask for spm_preproc8 because sometimes inital SPM segmentation
            % does not cover the whole brain for brains with thinner skull
            [Ym, Ycls] = cat_spm_preproc_write8(res,zeros(k,4),zeros(1,2),[0 0],1,1);
            Ym  = single(Ycls{1})/255 + single(Ycls{2})/255 + single(Ycls{3})/255;
            Yb  = (Ym > 0.5);
            Yb  = cat_vol_morph(cat_vol_morph(Yb,'lo'),'d',5/mean(vx_vol));
            res.biasreg   = obj.biasreg;
            res.biasfwhm  = obj.biasfwhm;
            res.reg       = obj.reg;
            res.samp      = obj.samp;
            res.tpm       = obj.tpm;
            res.fwhm      = obj.fwhm;
            res.msk       = res.image(1); 
            res.msk.pinfo = repmat([255;0],1,size(Yb,3));
            res.msk.dat(:,:,:) = Yb; % mask unused voxels!
            res = spm_preproc8(res);
  
            % final estimate without mask using parameters from previous run
            stime = cat_io_cmd('SPM preprocessing 1 (estimate skull-stripped):','','',job.extopts.verb-1,stime);
            res.biasreg   = obj.biasreg;
            res.biasfwhm  = obj.biasfwhm;
            res.reg       = obj.reg;
            res.samp      = obj.samp;
            res.tpm       = obj.tpm;
            res.fwhm      = obj.fwhm;
            res = spm_preproc8(res);
          end

        catch
            if any( (vx_vol ~= vx_voli) ) || ~strcmp(job.extopts.species,'human') 
                [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
                delete(fullfile(pp,[ff,ee]));
            end
            error('CAT:cat_run_job:spm_preproc8','Error in spm_preproc8. Check image and orientation. \n');
        end
        warning on 
        fprintf('%5.0fs\n',etime(clock,stime));   

        %% check contrast  
        clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
        Tgw = [cat_stat_nanmean(res.mn(res.lkp==1)) cat_stat_nanmean(res.mn(res.lkp==2))]; 
        Tth = [
          ... min(res.mn(res.lkp==6 & res.mg'>0.3)) ... % bg; ignore the background, because of MP2RGAGE, R1, and MT weighted images  
          max( min( clsint(3) ,  max(Tgw)+abs(diff(Tgw))) , min(Tgw)-abs(diff(Tgw)) ) ... % csf with limit for T2!
          clsint(1) ... gm
          clsint(2) ... wm 
          clsint(4) ... skull
          clsint(5) ... head tissue
          clsint(6) ... background
        ];
        
        res.Tth = Tth; 
        cat_err_res.res = res;   
        
        % inactive preprocessing of inverse images (PD/T2) 
        if job.extopts.INV==0 && any(diff(Tth(1:3))<=0)
          error('CAT:cat_main:BadImageProperties', ...
          ['CAT12 is designed to work only on highres T1 images.\n' ...
           'T2/PD preprocessing can be forced on your own risk by setting \n' ...
           '"cat12.extopts.INV=1" in the cat default file. If this was a highres \n' ...
           'T1 image then the initial segmentation might be failed, probably \n' ...
           'because of alignment problems (please check image orientation).']);    
        end

    end
    
    % updated tpm information for skull-stripped data should be available for cat_main
    if isfield(obj.tpm,'bg1')
      fname = res.tpm(1).fname;
      res.tpm       = obj.tpm;
      res.tpm(1).fname = fname;
    end
    spm_progress_bar('Clear');
            
    %% call main processing
    res.tpm    = tpm.V;
    res.stime  = stime0;
    res.catlog = catlog; 
    res.image0 = spm_vol(job.channel(1).vols0{subj}); 
    if exist('Ylesion','var'), res.Ylesion = Ylesion; else res.Ylesion = false(size(res.image.dim)); end; clear Ylesion;
    if exist('redspmres','var'); res.redspmres = redspmres; res.image1 = image1; end
    job.subj = subj; 
    cat_main(res,obj.tpm,job);
    
    % delete denoised/interpolated image
    [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
    if exist(fullfile(pp,[ff,ee]),'file'); 
      delete(fullfile(pp,[ff,ee]));
    end
    
%%
return
%=======================================================================

