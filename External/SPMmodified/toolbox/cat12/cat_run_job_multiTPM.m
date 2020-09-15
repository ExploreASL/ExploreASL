function [Affine,tpm,res] = cat_run_job_multiTPM(job,obj,Affine,skullstripped,msk,acc)
% ______________________________________________________________________
% Tissue Probability Maps (TPM) play an important role for the affine
% registration and segmentation. Preprocessing of atypical subjects 
% (e.g. very old/young) with typical SPM TPM often lead to inoptimal 
% results or can fail completelly. Hence, it is important to use at 
% least specific TPMs in children. 
%
% This function focus on the determination of the best fitting TPMs  
% to create a subject specific TPM. It allows to create linear sub
% levels for age specific TPMs for e.g. very young subjects with 
% small head, healty adults, and older subjects with severe brain
% atrophy.  
%
%   [Affine,tpm,res] = ...
%     cat_run_job_multiTPM(job,obj[,Affine,skullstripped,msk,acc])
%  
%   job           .. CAT job structure
%   obj           .. SPM preprocessing structure 
%   Affine        .. affine in/output matrix
%   skullstripped .. use brain tissues only
%   msk           .. mask regions in the unified segmentation 
%   acc           .. accuracy (default=0, just for tests)
%                    0 - faster but less exact
%                    1 - slower but more exact
%
% This function is part of the CAT preprocessing (cat_main_job) and 
% utilize the TPM based affine registration spm_maffreg of SPM and the 
% unified segmentation. 
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_run_job_multiTPM.m 1435 2019-03-06 11:24:38Z dahnke $
% ______________________________________________________________________

%#ok<*WNOFF,*WNON>
  
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  
  stime0 = clock; 
  
  if ~exist('Affine','var'), Affine = eye(4); end
  if ~exist('skullstripped','var'), skullstripped = 0; end
  if ~exist('acc','var') % overwrite parameter for tests
    obj.samp = job.opts.samp + 1;       % higher res (default: 3 mm)
    obj.tol  = min(1e-2,job.opts.tol * 10); % faster pp (default: 1-e4;)
  else
    acc = min(1,max(0,acc)); 
    obj.samp = 4 - 2*acc;       % higher res (default: 3 mm)
    obj.tol  = 10^(-2 - 3*acc); % faster pp (default: 1-e4;)
  end
  if ~exist('msk','var'), msk = 0; end
  
  
  % mask probably masked/stripped voxels!
  if isfield(obj,'image0'), obj = rmfield(obj,'image0'); end
  obj.image = cat_vol_resize(obj.image,'interpv',1.5);
  if msk
    if isfield(obj,'msk')
      obj.msk = cat_vol_resize(obj.msk,'interpv',1.5);
    else
      obj.msk           = obj.image(1);
      obj.msk.dt        = [spm_type('uint8') spm_platform('bigend')];
      obj.msk.pinfo     = repmat([255;0],1,obj.image(1).dim(3));
      if numel(msk)==1
        obj.msk.dat     = uint8(spm_read_vols(obj.image)==0);
      else
        obj.msk.dat     = uint8(msk>0);
      end
      if isfield(obj.msk,'private'),   obj.msk   = rmfield(obj.msk,'private'); end
      if isfield(obj.image,'private'), obj.image = rmfield(obj.image,'private'); end
      obj.msk           = spm_smoothto8bit(obj.msk,0.1); 
    end
  else
    if isfield(obj,'msk'), obj = rmfield(obj,'msk'); end
  end

  
  if skullstripped>1
    Vi = obj.image; 
    if isfield(Vi,'dat'), Vi = rmfield(Vi,'dat'); end
    [Vmsk,Yb] = cat_vol_imcalc([Vi;obj.tpm.V(1:3)],'blub', ...
      'i2 + i3 + i4',struct('interp',5,'verb',0)); clear Vmsk;  %#ok<ASGLU>
    obj.image.dat = obj.image.dat .* (Yb>0.1);
    if ~debug, clear Yb; end
  end
  
  
  
  % In case of skull-stripped images it is better to use a more simple
  % model with only 4 classes with CSF, GM, WM and background. 3 peaks
  % per tissue class seamed to be more stable as only 1 peak!
  if skullstripped
    job.opts.ngaus = 3*ones(4,1); % this is more safe
    obj.lkp        = [];
    for k=1:numel(job.opts.ngaus)
      job.tissue(k).ngaus = job.opts.ngaus(k);
      obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
    end
  else 
    % simplyfied model? 
    %obj.lkp = 1:6;
  end
  
  
  % prepare variables
  mAffine = cell(1,numel(job.opts.tpm)); 
  weight  = zeros(1,numel(job.opts.tpm)); 
  
  
  % estimate affine registration to the current TPM
  %res = struct(); 
  for i=1:numel(job.opts.tpm)
    if numel(job.opts.tpm)>1
      spm_plot_convergence('Init',...
        sprintf('Fine affine registration (TPM %d)',i),...
        'Mean squared difference','Iteration');
      if i==1
        fprintf('\n'); 
        stime = cat_io_cmd(sprintf('  Evaluate TPM %d (%s)',i,...
          spm_str_manip(job.opts.tpm{i},'ra30')),'g5');
      else
        stime = cat_io_cmd(sprintf('  Evaluate TPM %d (%s)',i,...
          spm_str_manip(job.opts.tpm{i},'ra30')),'g5','',job.extopts.verb-1,stime);
      end
    else
      stime = clock; 
      spm_plot_convergence('Init','Fine affine registration',...
        'Mean squared difference','Iteration');
    end

    
    % select and load TPM
    [pp,ff,ee] = spm_fileparts(job.opts.tpm{i}); 
    job.opts.tpm{i} = fullfile(pp,[ff,ee]);
    obj.tpm = spm_load_priors8(job.opts.tpm{i});
    if all(numel(obj.tpm.V)~=[4 6])
      error('cat_run_job_multiTPM:badtpm',[...
        'The TPM has to include 4 (GM,WM,CSF,BG) or 6 (GM,WM,CSF,skull,head,BG)\n' ...
        'tissue classes! Your TPM has %d! Is this the right file?\n' ...
        '  %s'],numel(obj.tpm.V),job.opts.tpm{i});
    end
    
    % update number of SPM gaussian classes 
    if skullstripped 
      Yb2 = 1 - spm_read_vols(obj.tpm.V(1)) - spm_read_vols(obj.tpm.V(2)) - spm_read_vols(obj.tpm.V(3));
      if 1
        for k=1:3
          obj.tpm.dat{k}     = spm_read_vols(obj.tpm.V(k));
          obj.tpm.V(k).dt(1) = 64;
          obj.tpm.V(k).dat   = double(obj.tpm.dat{k});
          obj.tpm.V(k).pinfo = repmat([1;0],1,size(Yb2,3));
        end
      end

      obj.tpm.V(4).dat = Yb2;
      obj.tpm.dat{4}   = Yb2; 
      obj.tpm.V(4).pinfo = repmat([1;0],1,size(Yb2,3));
      obj.tpm.V(4).dt(1) = 64;
      obj.tpm.dat(5:6) = []; 
      obj.tpm.V(5:6)   = []; 
      obj.tpm.bg1(4)   = obj.tpm.bg1(6);
      obj.tpm.bg2(4)   = obj.tpm.bg1(6);
      obj.tpm.bg1(5:6) = [];
      obj.tpm.bg2(5:6) = [];
      obj.tpm.V = rmfield(obj.tpm.V,'private');
      if ~debug, clear Yb2; end
    end
    
    % affine registration to TPM with (rought) and without smoothing (fine)
    try
      warning off  
      mAffine{i} = spm_maff8(obj.image(1),...
        obj.samp*2,(obj.fwhm+1)*10,obj.tpm,Affine,job.opts.affreg,160*(1+acc)); 
      mAffine{i} = spm_maff8(obj.image(1),...
        obj.samp*2,(obj.fwhm+1)*5,obj.tpm,mAffine{i},job.opts.affreg,80*(1+acc)); 
      mAffine{i} = spm_maff8(obj.image(1),...
        obj.samp,obj.fwhm,obj.tpm,mAffine{i},job.opts.affreg,40*(1+acc));
      mAffine{i} = spm_maff8(obj.image(1),...
        obj.samp/2,obj.fwhm,obj.tpm,mAffine{i},job.opts.affreg,20*(1+acc));
      warning on 
    catch
      mAffine{i} = Affine; 
    end
    % futher checks and corrections??
    
    
    % call Unified Segmentation
    if ~any(isnan(mAffine{i}(:)))
      obj.Affine = mAffine{i};
      try
        res(i) = cat_spm_preproc8(obj);  %#ok<*AGROW>
      catch
        res(i).wp  = nan(1,6);
        res(i).lkp = obj.lkp;
        res(i).mn  = nan(1,numel(obj.lkp));
        res(i).mg  = nan(numel(obj.lkp),1);
        res(i).vr  = nan(1,1,numel(obj.lkp));
      end
    else
      res(i).wp  = nan(1,6); 
      res(i).lkp = obj.lkp;
      res(i).mn  = nan(1,numel(obj.lkp));
      res(i).mg  = nan(numel(obj.lkp),1);
      res(i).vr  = nan(1,1,numel(obj.lkp));
    end
    
    % estimate weighting 
    % The variable res(i).wp gives good information how good segmentation 
    % fits to the TPM. The first 3 classes should be relative similar and
    % have high values (>0.1), whereas class 4 to 6 should have low values 
    % (<0.15).
% old version ...
%     weight(i) = min(1,max(0,2 * ( mean(res(i).wp(1:2)) - std(res(i).wp(1:2)) + ...
%                                   mean(res(i).wp(4:6)) - std(res(i).wp(4:6))) + ...
%                                   max(0,0.5 - 2 * sum( ...
%                                     shiftdim(res(i).vr( res(i).lkp(:) < 7 )) ./ ...
%                                     res(i).mn( res(i).lkp(:) < 7 )' .* ...
%                                     res(i).mg( res(i).lkp(:) < 7 ) ) ) ));
    weight(i) = min(1,max(0,1 - sum( shiftdim(res(i).vr) ./ ...
      res(i).mn' .* res(i).mg ./ mean(res(i).mn(res(i).lkp(2)))) ));  

    if numel(job.opts.tpm)>1
      fprintf('%s (fits%4.0f%%)',sprintf(repmat('\b',1,12)),weight(i)*100);
    end
  end 

 
  %% Mix of different TPMs
  if numel(job.opts.tpm)>1 
   % select most fitting TPM based on the affine registration
    [tmp,weighti] = sort(weight,'descend'); clear tmp tpm;  %#ok<ASGLU>
    
    tpm{1} = spm_load_priors8(job.opts.tpm{weighti(1)});
    tpm{2} = spm_load_priors8(job.opts.tpm{weighti(2)});
    
    if skullstripped 
      for i=1:2
        % update number of SPM gaussian classes 
        Ybg = 1 - spm_read_vols(tpm{i}.V(1)) - spm_read_vols(tpm{i}.V(2)) - spm_read_vols(tpm{i}.V(3));
        tpm{i}.V(4).dat = Ybg; 
        tpm{i}.dat{4}   = Ybg; 
        tpm{i}.V(4).pinfo = repmat([1;0],1,size(Ybg,3));
        tpm{i}.dat(5:6) = []; 
        tpm{i}.V(5:6)   = []; 
        tpm{i}.bg1(4)   = tpm{i}.bg1(6);
        tpm{i}.bg2(4)   = tpm{i}.bg1(6);
        tpm{i}.bg1(5:6) = [];
        tpm{i}.bg2(5:6) = [];
      end
    end

    weightn = round(weight(weighti(1:2))/sum(weight(weighti(1:2)))*100)/100;
    
    if weightn(1)<0.9
      Affine = mAffine{end};
      stime = cat_io_cmd(sprintf('  Merge TPM %d and %d (%0.0f:%0.0f)',weighti(1:2),weightn*100),'g5','',job.extopts.verb-1,stime);
    
      obj.tpm = spm_load_priors8(job.opts.tpm{weighti(1)});
      for k1=1:min([numel(tpm{1}.dat),numel(tpm{2}.dat)]);
        obj.tpm.dat{k1} = double(tpm{1}.dat{k1} * weightn(1) + ...
                                 tpm{2}.dat{k1} * weightn(2)); 
      end
      clear tpm1 tpm2;
    
      % final affine regisration with individual TPM
      i = numel(job.opts.tpm) + 1; 
      mAffine{i} = spm_maff8(obj.image(1),...
        obj.samp,obj.fwhm,obj.tpm,mAffine{weighti(1)},job.opts.affreg,40);
      failedAffine = mean(abs(mAffine{end}(:) - mAffine{weighti(1)}(:))) < 0.2; 
      
      % Final evaluation to test the result
      if ~failedAffine
        res(i).Affine = mAffine{end};
        Affine  = mAffine{end}; 
      end

      res(i) = cat_spm_preproc8(obj);
      weight(i) = min(1,max(0,1 - sum( shiftdim(res(i).vr) ./ ...
        res(i).mn' .* res(i).mg ./ mean(res(i).mn(res(i).lkp(2)))) ));                              
      fprintf('%s (fits%4.0f%%)',sprintf(repmat('\b',1,12)),weight(i)*100);
        
      if max(weight(1:end-1)) > weight(i)  %failedAffine || 
         cat_io_cmd(sprintf('  Merging failed: Use only TPM %d!',weighti(1)),...
           'g5','',job.extopts.verb-1,stime);
         
         Affine  = mAffine{weighti(1)}; 
         obj.tpm = spm_load_priors8(job.opts.tpm{weighti(1)});
      else
        cat_io_cmd(' ','g5','',job.extopts.verb-1,stime);
      end
    else
      obj.tpm = job.opts.tpm{weighti(1)}; 
      Affine  = mAffine{weighti(1)}; 
      cat_io_cmd(sprintf('  Use only TPM %d.',weighti(1)),'g5','',job.extopts.verb-1,stime);
    end
  
    
    cat_io_cprintf('g5',sprintf('%5.0fs\n',etime(clock,stime0)));
    cat_io_cmd(' ','','',job.extopts.verb-1);
  
  else
    Affine  = mAffine{1}; 
  end
  
  
  tpm = obj.tpm;
    
end