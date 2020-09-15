function [Ym,Ybg,WMth,bias] = cat_run_job_APP_SPMinit(job,tpm,ppe,n,ofname,nfname,mrifolder,skullstripped)
%% APP bias correction (APP1 and APP2)
%  ------------------------------------------------------------
%  Bias correction is essential for stable affine registration. 
%  SPM further required Gaussian distributed data that is 
%  achieved by Smoothing in high resolution data and by 
%  additional noise in regions with many zeros typical in 
%  skull-stripped or defaced data. 
%
%  [Ym,Yt,Ybg,WMth,bias] 
%  ------------------------------------------------------------

  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  
  
  V = spm_vol(job.channel(n).vols{job.subj});
  vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
  
  [pp,ff] = spm_fileparts(ofname); 
  ppn     = spm_fileparts(nfname); 
  onfname = fullfile(ppn,['o' ff '.nii']);

  stime = cat_io_cmd('APPs bias correction');
  if job.extopts.verb>1, fprintf('\n'); end
  stime2 = cat_io_cmd('  Preparation','g5','',job.extopts.verb-1);
  %if debug, copyfile(ofname,nfname); end

  % SPM segmentation parameter
  preproc.channel.vols{1,1}   = nfname; 
  preproc.channel.biasreg     = min(0.01,max(0.0001,job.opts.biasreg));
  preproc.channel.biasfwhm    = min(90,max(30,job.opts.biasfwhm/2));
  preproc.channel.write       = [0 1];
  for ti=1:6
    preproc.tissue(ti).tpm    = {[tpm.V(ti).fname ',' num2str(ti)]};
    preproc.tissue(ti).ngaus  = job.opts.ngaus(ti);
    preproc.tissue(ti).native = [0 0];                           % native, dartel
    preproc.tissue(ti).warped = [0 0];                           % unmod, mod
  end
  preproc.warp.mrf            = 1;                               % 
  preproc.warp.cleanup        = 1;                               % this is faster and give better results!
  preproc.warp.reg            = [0 0.001 0.5 0.05 0.2];
  preproc.warp.affreg         = job.opts.affreg; 
  preproc.warp.write          = [0 0];
  %preproc.warp.fwhm           = 0;
  %preproc.warp.samp           = min(9,max(2,job.opts.samp*2));

  if skullstripped 
    %% update number of SPM gaussian classes 
    preproc.tpm = tpm; 
    Ybg = 1 - spm_read_vols(preproc.tpm.V(1)) - spm_read_vols(preproc.tpm.V(2)) - spm_read_vols(preproc.tpm.V(3));
    if 1
      for k=1:3
        preproc.tpm.dat{k}     = spm_read_vols(preproc.tpm.V(k));
        preproc.tpm.V(k).dt(1) = 64;
        preproc.tpm.V(k).dat   = double(preproc.tpm.dat{k});
        preproc.tpm.V(k).pinfo = repmat([1;0],1,size(Ybg,3));
      end
    end
    preproc.tissue(5:end)  = []; 
    preproc.tpm.V(4).dat   = Ybg;
    preproc.tpm.dat{4}     = Ybg; 
    preproc.tpm.V(4).pinfo = repmat([1;0],1,size(Ybg,3));
    preproc.tpm.V(4).dt(1) = 64;
    preproc.tpm.dat(5:6)   = []; 
    preproc.tpm.V(5:6)     = []; 
    preproc.tpm.bg1(4)     = preproc.tpm.bg1(6);
    preproc.tpm.bg2(4)     = preproc.tpm.bg1(6);
    preproc.tpm.bg1(5:6)   = [];
    preproc.tpm.bg2(5:6)   = [];
    preproc.tpm.V          = rmfield(preproc.tpm.V,'private');

    % tryed 3 peaks per class, but BG detection error require manual 
    % correction (set 0) that is simple with only one class  
    job.opts.ngaus = [([job.tissue(1:3).ngaus])';1]; % 3*ones(4,1);1; 
  end
  preproc.lkp        = [];
  for k=1:numel(job.opts.ngaus)
    preproc.ngaus(k)        = job.opts.ngaus(k);
    preproc.tissue(k).ngaus = job.opts.ngaus(k);
    preproc.lkp             = [preproc.lkp ones(1,job.tissue(k).ngaus)*k];
  end
  
  %% add noise in zero regions (skull-stripping / defacing)
  VF = spm_vol(nfname);
  YF = single(spm_read_vols(VF));
  % some average object intensity 
  Tthn = cat_stat_nanmean(YF(YF(:)>cat_stat_nanmean(YF(YF(:)~=0)) & YF(:)~=0));
  % limitation is required for division data (yv98_05mm_corrected_
  YF = min( Tthn*10, YF); 
  % smoothing for Gaussian distribution
  YF(isnan(YF(:)))=0;
  YFs = YF+0; spm_smooth(YFs,YFs,0.5./vx_vol);
  YM = abs( (YF-YFs)./max(eps,YFs));
  YM = min(1,smooth3(YM ./ min(0.2,max(YM(:)))));
  YF = YFs.*YM + YF.*(1-YM);
  Y0 = YF==0; 
  % add noise for Gaussian distribution in the background
  if ppe.affreg.skullstripped
    YF = YF + (Y0) .* (0.05*Tthn).*rand(size(YF));
  else
    YF(cat_vol_morph(Y0,'o'))=nan; 
  end
  % force floating point
  VF.dt(1) = 16; 
  VF = rmfield(VF,'private'); 
  if exist(VF.fname,'file'); delete(VF.fname); end
  spm_write_vol(VF,YF);  
  clear VF YF Tthn; 
            
  %% try SPM preprocessing
  %  ----------------------------------------------------------
  %  * if SPM failed we go on without bias correction 
  %  * further iterations on different resolution are usefull
  %    (the first run correct the worst problems and should 
  %     a correct registration of the next runs)
  %  * best results for strong corrections (fwhm 30 to 45 mm)
  %  * SPM preprocessing is very fast (and gives you results in 
  %    the template resolution?) but writing results in the 
  %    original resolution is very slow, so using 2 iterations
  %    is maybe optimal
  %  * much slower for biasfwhm<35 mm and samp<4.5 mm
  %  * the final operation is much slower (3 times) becauserun it
  %    required the estimation of the deformation field to
  %    write data in the orignal space
  %  * writing further outputs on disk will also cost more time
  %    >> further optimization is possible by avoiding
  %       temporary disk output 
  %  ----------------------------------------------------------
  if job.opts.biasfwhm<45                                   % strong (9 iterations)
    fwhmx  = 3:-1/4:1; 
    sampx  = 2:-1/(numel(fwhmx)-1):1.0; 
  elseif job.opts.biasfwhm>=45 && job.opts.biasfwhm<75      % medium (6 iterations) 
%    fwhmx  = 3:-1/2:1/2; 
%    sampx  = 2:-1/(numel(fwhmx)-1):1.0; 
    fwhmx  = [2,1]; %,1/2]; 
    sampx  = 2:-1/(numel(fwhmx)-1):1.0; 
  elseif job.opts.biasfwhm>=75                              % light (3 iteration)
    fwhmx  = 2:-3/4:1/2; 
    sampx  = 2:-1/(numel(fwhmx)-1):1.0; 
  end
  

  spmp0    = debug;  % job.extopts.verb>1; % for debugging: 0 - remove all APP data, 1 - save Ym, 2 - save Ym and Yp0
  optimize = 0;      % use low resolution (in mm) input image to increase speed >> limited speed advantage :/
                     % deformation requires a lot of time ... maybe its enough to use it only in case of ultra highres data
  if any(vx_vol<1.5), optimize = 1.5; end            

  preproc.Yclsout = false(1,6); 
  preproc.tol     = min(1e-2,job.opts.tol * 10); % 2 times less accurant than the final operation
  copyfile(nfname,onfname);
  
  if optimize>0 
    %% lower resolution to improve SPM processing time and use smoothing for denoising
    Vi        = spm_vol(nfname); 
    Vi        = rmfield(Vi,'private'); Vn = Vi;
    imat      = spm_imatrix(Vi.mat); 
    Vi.dim    = round(Vi.dim .* vx_vol./repmat(optimize,1,3));
    imat(7:9) = repmat(optimize,1,3) .* sign(imat(7:9));
    Vi.mat    = spm_matrix(imat);

    spm_smooth(nfname,nfname,repmat(0.75,1,3)); % denoising
    [Vi,Yi] = cat_vol_imcalc(Vn,Vi,'i1',struct('interp',5,'verb',0,'mask',-1));
    delete(nfname); spm_write_vol(Vi,Yi); 
    %if ~isinf(job.opts.biasstr), clear Yi; end
  end

  %%
  bias = zeros(size(sampx)); 
  for ix=1:numel(sampx)
    %% parameter update
    preproc.warp.samp         = min(9  ,max(1 ,job.opts.samp     * sampx(ix))); 
    preproc.channel.biasfwhm  = min(180,max(30,job.opts.biasfwhm * fwhmx(ix)));
    % preproc.warp.reg          = [0 0 0 0 0];ex

    stime2 = cat_io_cmd(sprintf('  SPM bias correction (samp: %0.2f mm, fwhm: %3.0f mm)',preproc.warp.samp,...
      preproc.channel.biasfwhm),'g5','',job.extopts.verb-1,stime2);

    %try
      % SPM bias correction
      warning off;  %#ok<WNOFF>
      if ix==numel(sampx) || bias(1) >= preproc.channel.biasfwhm 
        if job.extopts.APP==2 || spmp0>1, preproc.Yclsout = true(1,max(preproc.lkp)); end % 
        vout = cat_spm_preproc_run(preproc,'run'); res(ix) = vout.res;
      else
        vout = cat_spm_preproc_run(preproc,'run'); res(ix) = vout.res;
      end
      warning on;  %#ok<WNON>
%%
      Pmn    = fullfile(pp,mrifolder,['mn' ff '.nii']); 
      Pmn_ri = fullfile(pp,mrifolder,['mn' ff '_r' num2str(ix) '.nii']);
      Pmn_r0 = fullfile(pp,mrifolder,['mn' ff '_r0.nii']);

      % estimate bias strength based on the applied corrections
      % of the initial correction .. 
      % in case of updates, local biasfield strenght is maybe
      % better (only useful if strong changes are allowed)
      if ix==1 && exist('Yi','var')
        Vn = spm_vol(Pmn); Vn = rmfield(Vn,'private'); 
        Yn = spm_read_vols(Vn);
        bias(ix) = (1/cat_stat_nanstd(Yn(:)./Yi(:))) * 4; 
        %fprintf('bias=%5.0f mm ',bias(ix)); 
        bias(ix) = max(30,min(120,round(bias(ix) / 15) * 15));
        if ~debug, clear Yn; end
      else 
        bias(ix) = 0; 
      end

      % backup the bias corrected image
      if spmp0>0, copyfile(Pmn,Pmn_ri); end

      % backup for mixing
      if (ix==2 && numel(sampx)>2) || (ix==1 && numel(sampx)<=2) || bias(1) > preproc.channel.biasfwhm 
        copyfile(fullfile(pp,mrifolder,['mn' ff '.nii']),Pmn_r0); 
      end
      copyfile(fullfile(pp,mrifolder,['mn' ff '.nii']),nfname);

      % write segmentation 
      if spmp0>1 && exist('vout','var') && isfield(vout,'Ycls') 
        VF  = spm_vol(nfname); VF.fname = fullfile(pp,mrifolder,['p0n' ff '.nii']); VF.dt(1) = 2; 
        Yp0 = single(vout.Ycls{1})/255*2 + single(vout.Ycls{2})/255*3 + single(vout.Ycls{3})/255;
        spm_write_vol(VF,Yp0);  
      elseif spmp0==0 && ix==numel(sampx) % remove spm mat file
        delete(fullfile(pp,mrifolder,['n' ff '_seg8.mat']));
      end

      %% combine low and high frequency filted images 
      if ix==numel(sampx) %|| (bias(1) > preproc.channel.biasfwhm && ix>1)
        %%
        stime2 = cat_io_cmd('  Postprocessing','g5','',job.extopts.verb-1,stime2);
        if optimize
          % update Ym
          Vn  = spm_vol(nfname);  Vn = rmfield(Vn,'private'); 
          Vo  = spm_vol(onfname); Vo = rmfield(Vo,'private'); 
          [Vmx,vout.Ym] = cat_vol_imcalc([Vo,Vn],Vo,'i2',struct('interp',5,'verb',0,'mask',-1)); clear Vmx; %#ok<ASGLU>
          vout.Ym = single(vout.Ym); 

          %% remap Ym0
          Vn = spm_vol(Pmn_r0);  Vn = rmfield(Vn,'private'); 
          Vo = spm_vol(onfname); Vo = rmfield(Vo,'private'); 
          [Vmx,Ym0] = cat_vol_imcalc([Vo,Vn],Vo,'i2',struct('interp',5,'verb',0,'mask',-1)); clear Vmx; %#ok<ASGLU>
          Ym0 = single(Ym0); 

          %% replace interpolation boundary artefact
          %{
          Ybd     = ones(size(vout.Ym)); Ybd(3:end-2,3:end-2,3:end-2) = 0; 
          [D,I]   = cat_vbdist(single(~(isnan(vout.Ym) | Ybd)),true(size(vout.Ym))); clear D;  %#ok<ASGLU>
          vout.Ym = min(vout.Ym,vout.Ym(I));        
          Ym0     = Ym0(I);        
          clear I;
          %}

          %% update Ycls
          if isfield(vout,'Ycls')
            for i=1:numel(vout.Ycls)
              Vn  = spm_vol(nfname);  Vn = rmfield(Vn,'private');
              Vo  = spm_vol(onfname); Vo = rmfield(Vo,'private');
              Vn.pinfo = repmat([1;0],1,size(vout.Ycls{i},3));
              Vo.pinfo = repmat([1;0],1,Vo.dim(3));
              Vn.dt    = [2 0];
              Vo.dt    = [2 0]; 
              Vo.dat   = zeros(Vo.dim(1:3),'uint8');
              if ~isempty(vout.Ycls{i})
                Vn.dat   = vout.Ycls{i}; 
                [Vmx,vout.Ycls{i}]  = cat_vol_imcalc([Vo,Vn],Vo,'i2', ...
                  struct('interp',5,'verb',0,'mask',-1)); clear Vmx; %#ok<ASGLU>
                vout.Ycls{i} = cat_vol_ctype(vout.Ycls{i}); 
              end
            end
          end
        else
          Ym0 = spm_read_vols(spm_vol(Pmn_r0)); 
        
          %% replace interpolation boundary artefact
          %{
          Ybd   = ones(size(vout.Ym)); Ybd(3:end-2,3:end-2,3:end-2) = 0; 
          [D,I] = cat_vbdist(single(~(isnan(Ym0) | Ybd)),true(size(Ym0))); clear D;  %#ok<ASGLU>
          Ym0   = min(Ym0,Ym0(I)); clear I;
          %}
        end
        if debug, stime2 = cat_io_cmd(' ','g5','',job.extopts.verb-1,stime2); end

        if ~debug && exist(Pmn_r0,'file'), delete(Pmn_r0); end
        if ~debug && exist(Pmn_r0,'file'), delete(Pmn_r0); end

        %% mixing
        %  creation of segmentation takes a lot of time because
        %  of the deformations. So it is much faster to load a
        %  rought brain mask. 
        if isfield(vout,'Ycls')
          Yp0 = single(vout.Ycls{1})/255*2 + single(vout.Ycls{2})/255*3 + single(vout.Ycls{3})/255;
          YM2 = cat_vol_smooth3X(cat_vol_smooth3X(Yp0>0,16/mean(vx_vol))>0.95,10/mean(vx_vol)); 
          if ~debug, clear Yp0; end
        else
          Pb  = char(job.extopts.brainmask);
          Pbt = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
          VF  = spm_vol(onfname); 
          VFa = VF; %if job.extopts.APP~=5, VFa.mat = Affine * VF.mat; end
          [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1));  %#ok<ASGLU>
          Ybb = cat_vol_smooth3X(Yb>0.5,8/mean(vx_vol)); Ybb = Ybb./max(Ybb(:)); 
          YM2 = cat_vol_smooth3X(Ybb>0.95,8/mean(vx_vol)); YM2 = YM2./max(Ybb(:)); 
          if ~debug, clear Pb Pbt VFa Vmsk Yb Ybb; end
        end
        
        %% combine the low (~60 mm, Ym0) and high frequency correction (~30 mm, vout.Ym) 
        Yo  = spm_read_vols(spm_vol(onfname)); 
        
        % final bias field correction
        Yw  = vout.Ym.*(1-YM2) + (YM2).*Ym0;
        Yw  = Yo./ max(eps,Yw) .* (Yo~=0 & Yw~=0); 
        % correct undefined voxel and assure smoothness of the bias field 
        Yw  = cat_vol_approx(max(0,min(2,Yw)),'',vx_vol,2,struct('lfO',2));
        vout.Ym = Yo ./ max(eps,Yw);
        
        %% output variables
        %  the background class is maybe incorrect, so we use the minimum 
        %  of the most important values per class
        BGth = min(vout.res.mn(vout.res.mg(:) >= (0.5./preproc.ngaus(vout.res.lkp(:))')));
        WMth = mean(vout.res.mn(vout.res.lkp(:)==2) .* vout.res.mg(vout.res.lkp(:)==2)');
        vout.Ym(isnan(vout.Ym(:))) = min(vout.Ym(:)); 
        Ym   = (vout.Ym - BGth) ./ (WMth - BGth); 
        Ybg  = cat_vol_morph(cat_vol_morph(Ym<0.2,'ldo',2,vx_vol),'dc',3,vx_vol);
        
        %  The SPM bias correction can change the overall image intensity of 
        %  some images strongly (resulting in negative values!). Hence, we 
        %  have to rescale the intensity based on the first segmentation.
        BGth0   = min(res(1).mn(res(1).mg(:) >= (0.5./preproc.ngaus(res(1).lkp(:))')));
        WMth0   = mean(res(1).mn(res(1).lkp(:)==2) .* res(1).mg(res(1).lkp(:)==2)');
        PDT2    = res(1).mn(res(1).lkp(:)==1) < res(1).mn(res(1).lkp(:)==2); % inverse weighting > softer scalling 
        Yoc     = max(-0.1,min(3 + 7*PDT2 , Ym .* WMth0 + BGth0)); 
        
Yoc(isnan(Yo)) = nan; % not sure if this is good
        
        %%
        if ~debug, clear Yw Yo; end
        Vm = spm_vol(onfname); Vm.fname = nfname; Vm.dt(1) = 16;  
        spm_write_vol(Vm,Yoc); 

        % backup the bias corrected image
        if spmp0>0 
          copyfile(nfname,fullfile(pp,mrifolder,['mn' ff '_r' num2str(ix)+1 '.nii'])); 
        end
        if exist(Pmn_r0,'file'), delete(Pmn_r0); end
        %%
        break
      else
        movefile(fullfile(pp,mrifolder,['mn' ff '.nii']),nfname); 
      end
    try
    catch
      fprintf('\b\b\b\b\b\b\b\b\b(failed) ');   
      if exist(fullfile(pp,mrifolder,['mn' ff '.nii']),'file')
        delete(fullfile(pp,mrifolder,['mn' ff '.nii']));
      end
    end
    
    if 0
      %% just debugging
      Yp0  = single(vout.Ycls{1})/255*2 + single(vout.Ycls{2})/255*3 + single(vout.Ycls{3})/255;
      Tthx = max( [ median( vout.Ym(vout.Ycls{1}(:)>128)) , median( vout.Ym(vout.Ycls{2}(:)>128)) ]);
      ds('d2','a',vx_vol,Ym0/Tthx,Yp0/3,vout.Ym/Tthx,vout.Ym/Tthx,120)
      %%
      Tthx = mean( vout.Ym( YM2(:)>0.2 & vout.Ym(:)>mean(vout.Ym(YM2(:)>0.2))) );
      ds('d2','',vx_vol,Ym0/Tthx,YM2,vout.Ym/Tthx,Ym/Tthx,80)
    end
  end
  if debug, cat_io_cmd('','g5','',job.extopts.verb-1,stime); end
  Pmn = fullfile(pp,mrifolder,['mn' ff '.nii']); 
  if exist(Pmn,'file'), delete(Pmn); end

  
  
  %% try APPs preprocessing
  %  ---------------------------------------------------------
  %  SPM may failed in correction of high intensive gyri that 
  %  should be fixed by the APPs correction that is similar to
  %  the APPinit and APPfinal routines but used SPM as input 
  %  segmentation.
  %  However this can work very well in general, problems can 
  %  occure in subcortical structures and the cerebellum, 
  %  especially in special protocols. 
  %  It is also possible to skip this function if the estimated
  %  bias field is very small, but I am not sure if the is useful.
  %  ----------------------------------------------------------
  if exist('vout','var') && job.extopts.APP==2 && isfield(vout,'Ycls')
    stime2 = cat_io_cmd('  APP bias correction','g5','',job.extopts.verb-1,stime2);  %#ok<NASGU>

    try
      if isinf(job.opts.biasstr), catbias = 1 - (bias(1)-30)/60; else catbias = job.opts.biasstr; end
      [Yo,Ym] = cat_run_job_APP_SPMfinal(onfname,vout,vx_vol * (2 - strcmp(job.extopts.species,'human')),...
        job.extopts.verb-1,catbias); 
      spm_write_vol(spm_vol(nfname),Yo);  
     if ~debug, clear vout Yo Ym0 YM2; end  
      if spmp0
        copyfile(nfname,fullfile(pp,mrifolder,['mn' ff '_r' num2str(ix)+2 '.nii'])); 
      end
    catch
      fprintf('failed\n');
    end
  end
  
  if exist(onfname,'file'), delete(onfname); end
  fprintf('%5.0fs\n',etime(clock,stime));  




end