function [Ym,Yt,Ybg,WMth,bias,Tth,pior] = cat_run_job_APP_init(Ysrco,vx_vol,opt)
%  _____________________________________________________________________
%  The rough bias correction is a subfunction of cat_run_rob.
% 
%  All tissues (low gradient areas) should have a similar intensity.
%  A strong smoothing of this approximation is essential to 
%  avoid anatomical filtering between WM and GM that can first 
%  be seen in overfitting of the subcortical structures.
%  However, this filtering will overcorrect head tissue with
%  a typical intensity around GM.
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id: cat_run_job_APP_init.m 1407 2018-12-15 19:58:05Z dahnke $

  if ~exist('opt','var'); opt = struct(); end
  def.APPstr = 0.5;  % strength of bias correction (smoothness of bias field; 0=smooth; 1=hard)
  def.APPred = 2.2;  % resolution limit for estimation  
  def.icall  = 1;    % test inverse processing  
  def.verb   = 1;    % verbose level
  def.iproc  = 0; 
  opt = cat_io_checkinopt(opt,def); 

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  
%    ds('l2','',0.5,Yo/WMth,Yg<0.2,Yo/WMth,Ym,80)
  Ysrcmin = min(Ysrco(~isinf(Ysrco(:))));
  Ysrco   = Ysrco - Ysrcmin;

  rf = 10^9; 
  if opt.verb, fprintf('\n'); end
  
  stime = cat_io_cmd('  Initialize','g5','',opt.verb);
  msize = 48; 

  %% initial noise reduction by resolution (for high resolution we expect more noise then in low resolution) 
  nsmooth = min( 1.6 , max( 0 , 2 - mean(vx_vol) ) * 0.8 ); 
  Ysrc    = Ysrco + 0; spm_smooth(Ysrc,Ysrc,nsmooth);
  [Ysrc,resT1] = cat_vol_resize(Ysrc,'reduceV',vx_vol,2,msize,'meanm'); 

  % correction for negative backgrounds (MT weighting)
  [Ym,BGth] = cat_stat_histth(Ysrc,99.99); BGth(2) = []; clear Ym;   %#ok<ASGLU>
  
  % initial thresholds
  % correction for negative backgrounds (MT weighting)
  [Ym,BGth] = cat_stat_histth(Ysrc,99.99); BGth(2) = []; clear Ym;   %#ok<ASGLU>
  
  Ysrc = Ysrc - BGth; Ysrco = Ysrco - BGth; BGth2 = BGth; 
  Yg   = cat_vol_grad(Ysrc,vx_vol) ./ max(eps,Ysrc); 
  Ygth = max(0.05,min(0.3,mean(Yg(Yg(:)>0))*0.75)); 

  
  Ybg0 = cat_vol_smooth3X(Yg<0.3 & Ysrc<cat_stat_nanmean([BGth cat_stat_nanmean(Ysrc(:))]),8)>0.05;
  Yw0  = Yg<0.3 & Ysrc>cat_stat_nanmean(  Ysrc(Yg(:)<0.3 & ~Ybg0(:))) & ~Ybg0; 
  WMth0 = roundx(single(cat_stat_nanmedian(Ysrc( Yw0(:) ))),rf); clear Ybg0 Yw0;
  BGth = max( min(Ysrc(:))*0.7 + 0.3*cat_stat_nanmean(Ysrc(:)) ,...
    cat_stat_nanmean(Ysrc(Ysrc(:)<cat_stat_nanmean(Ysrc(:))))); BGth = roundx(BGth,rf); 
  highBG = BGth>0.3; % MT / MP2RAGE
  if highBG 
    % In case of high intensity background we simply need the lowest 
    % intensity of the histogram without extrem outliers.
    Ysrcr = cat_vol_resize(Ysrc,'reduceV',vx_vol,2,24,'meanm'); 
    [Ysrcr,th] = cat_stat_histth(Ysrcr,99.99); %#ok<ASGLU>
    BGth = th(1);
    clear Ysrcr;
  end
  Ym   = (Ysrc - BGth) ./ (WMth0 - BGth);
  
  %% improved WM threshold
  Yg    = cat_vol_grad(Ym,resT1.vx_volr) ./ max(0.3,Ym); 
  Ydiv  = cat_vol_div(Ym,resT1.vx_volr/2) ./ (Ym+eps); Ydiv(Ydiv==0)=eps; % lower resolution is 8 times faster 
  if highBG 
    if opt.iproc
      Ymsk  = smooth3(Yg>0 & Yg<Ygth & Ym<0.8 & Ym<3)>0.5; 
    else
      Ymsk  = smooth3(Yg>0 & Yg<Ygth & Ym>0.2 & Ym<3)>0.5; 
    end
    WMth2 = roundx(single(cat_stat_nanmean( Ysrc( Ymsk(:) ) )),rf); if ~debug, clear WMth1 Ymsk, end
    Ym    = (Ysrc - BGth) ./ (WMth2 - BGth);
  end
  
  %% background
  stime = cat_io_cmd('  Estimate background','g5','',opt.verb,stime);
  if highBG
    % As far as our background have an intensity similar to the brain, we 
    % have to find the average intensity of the background create a mask.
    Ymr   = cat_vol_resize(Ysrc,'reduceV',vx_vol,3,32,'meanm'); 
    Ymsk  = false(size(Ymr)); Ymsk(3:end-2,3:end-2,3:end-2) = true;
    [HBGth,HBGsd] = kmeans3D(Ymr(Ymsk(:)),3); HBGsd = min(0.05, HBGsd);
    HBGth = HBGth(2); HBGsd = HBGsd(2); % use averge BG intensity
    Ybg   = Yg<0.1 & Ysrc>(HBGth - 2*HBGsd) & Ysrc<(HBGth + 2*HBGsd); 
    [Ybg,resT2] = cat_vol_resize(Ybg,'reduceV',vx_vol,2,32,'meanm'); 
    Ybg   = cat_vol_morph(Ybg>0.5,'ldo',2);
    Ybgi  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5; % inverse background = object
    
    %% estimate WM threshold
    Yw0   = Yg<0.2 & Ysrc>cat_stat_nanmean(  Ysrc(Yg(:)<0.3 & Ybgi(:))) & Ybgi; 
    Yw0   = cat_vol_morph(smooth3(Yw0)>0.5,'ldo'); 
    WMth2 = kmeans3D(Ysrc(Yw0(:)),1); 
    Ym    = (Ysrc - BGth) ./ (WMth2 - BGth);
    Ybgc  = Ybg;
  else
    Ybg   = Yg<(mean(Yg(Ym(:)~=0))) | isnan(Yg);
    % avoid edges on the image border
    bb = 1; Ybb1 = true(size(Ybg)); Ybb1(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
    bb = 4; Ybb2 = true(size(Ybg)); Ybb2(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
    Ybgc  = Ybb2 & cat_vol_morph(Ybg | (Ybb1 & Yg<0.5),'c',2);
    Ybgi = Ybg | Ybgc | smooth3(Yg./Ydiv > 100)<0.5;
    if ~debug, clear Ybb1 Ybb2 bb; end
  end 
  %% filling
  [Ybg,resT2] = cat_vol_resize(single(Ybgi),'reduceV',resT1.vx_volr,2,32,'meanm'); 
  Ybg  = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',4);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5;    
  zeroBG = cat_stat_nanmean(Ym(Ybg(:)>0))<0.2;
  if zeroBG
    BGth1 = cat_stat_nanmean(Ysrc(Ybg(:))); 
    BGth2 = cat_stat_nanmean(Ysrc(Ybg(:))) - cat_stat_nanmean(Ysrc(Ybg(:))) ;  
  else
    BGth1 = BGth; 
    BGth2 = BGth; 
  end
 %%
  if opt.iproc
    WMth3 = WMth2; % * roundx(single(cat_stat_nanmedian(Ym(Yg(:)>0 & Yg(:)<Ygth & Yg(:)./max(eps,abs(Ydiv(:)))<0.5 & ~Ybg(:) & ...
     %        Ym(:)<cat_stat_nanmean(Ym(Yg(:)<0.2 & ~Ybg(:) & Ym(:)<cat_stat_nanmean(Ym(:))))))),rf); if ~debug, clear WMth2, end
  else
    WMth3 = WMth2 * roundx(single(cat_stat_nanmedian(Ym(Yg(:)>0 & Yg(:)<Ygth & ~Ybg(:) & ...
             Ym(:)>cat_stat_nanmean(Ym(Yg(:)<Ygth & ~Ybg(:) & Ym(:)>cat_stat_nanmean(Ym(:))))))),rf);
     if ~debug, clear WMth2, end
  end
  
  %% initial noise reduction by resolution and by amount of variance in the background (multiply by 3 for 3 tissue classes)      
  Ywmi  = smooth3(cat_vol_morph( (Ym-Yg)>0.1 & (Ym+Yg)<1.2 & Yg<Ygth & abs(Ydiv)<0.2 , 'l'))>0.5; 
  Ymm   = Ym .* Ywmi; Ymm = cat_vol_resize(Ymm,'reduceV',resT1.vx_volr,resT1.vx_volr*3,16,'max');
  for i=1:1, Ymm = cat_vol_localstat(Ymm,Ymm>0,2,1); end;  
  bias0   = cat_stat_nanstd(Ymm(Ymm>0)) ./ cat_stat_nanmean(Ymm(Ymm>0)); clear Ymm;
  noise0  = cat_stat_nanmean( Yg( Ybg(:) ) ); 
  % smoothing parameter
  nsmooth = double(min( 1.6 , max( 0 , 2 - mean(vx_vol) ) * 0.8 ) *3*3*noise0 ); % noise reduction by resolution
  bsmooth = 6 + min(2,max(0,2 - (bias0 * 25))) * 2*(1-opt.APPstr) .* mean(resT1.vx_volr); 
  res     = cat_stat_nanmean(vx_vol) * max( 1 , min( opt.APPred , 0.1/max(0,bias0-0.05) )); 
  
  %%
  Ysrc  = Ysrco + 0; spm_smooth(Ysrc,Ysrc,nsmooth);
  [Ysrc,resT3] = cat_vol_resize(Ysrc,'reduceV',vx_vol,res,msize,'meanm'); 
  if any(resT1.vx_volr~=resT3.vx_volr)
    Ybg   = cat_vol_resize(single(Ybg) ,'dereduceV',resT1);    
    Ybgc  = cat_vol_resize(single(Ybgc),'dereduceV',resT1);
    Ybg   = cat_vol_resize(Ybg ,'reduceV',vx_vol,res,msize,'meanm')>0.5; 
    Ybgc  = cat_vol_resize(Ybgc,'reduceV',vx_vol,res,msize,'meanm')>0.5; 
  end
  Ym    = (Ysrc - BGth1) ./ (WMth3 - BGth1);
  Yg    = cat_vol_grad(Ym,resT3.vx_volr) ./ max(0.3,Ym); 
  Ydiv  = cat_vol_div(Ym,resT3.vx_volr/2) ./ (Ym+eps); Ydiv(Ydiv==0)=eps; % lower resolution is 8 times faster 
 
   
  %% first WM inhomogeneity with low tissue boundary (may include CSF > strong filtering for IXI175)
  stime = cat_io_cmd('  Initial correction','g5','',opt.verb,stime);
  Yms  = cat_vol_smooth3X( min(2 .* ~Ybg,Ym .* (Ydiv>-0.2) .* ~Ybg .* (Ym>0.1)),16*mean(vx_vol));     % this map is to avoid CSF in the mask!
  Yms  = (Yms ./ mean(Yms(~Ybg(:))));
  Yms  = cat_vol_smooth3X( min(Yms*1.5 .* ~Ybg,Ym .* ~Ybg),16*mean(vx_vol));
  Yms  = (Yms ./ mean(Yms(~Ybg(:))));
  if opt.iproc
    Yt   = Ym>max(0,Yms*0.4) & Ym<Yms*2 & Ym<(1+Yms*2) & Yg<Ygth*2 & Ydiv<0.2 & ~Ybg & Ym<1.2 & ...
           Ydiv>-0.6 & smooth3(Ym./(Yms+eps).*Yg.*Ydiv<-0.2)<0.3 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  else
    Yt   = Ym>max(0,Yms*0.3) & Ym<Yms*2 & Ym<(1+Yms*2) & Yg<0.5 & Ydiv<0.2 & ~Ybg & ...
         Ydiv>-0.6 & smooth3(Ym./(Yms+eps).*Yg.*Ydiv<-0.2)<0.3 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  end
  Yt   = Yt & cat_vol_morph(Yt & Yms>0.5 | (Yms>1.2),'l');       
  if ~highBG
    Ywi  = (Ym .* Yt) ./ max(eps,Yt);  
    if ~zeroBG 
      Ybg2 = Ybg(:) & Yg(:)<(cat_stat_nanmean(Yg(Ybg(:))) + 2*(cat_stat_nanstd(Yg(Ybg(:))))); 
      Ywi(Ybg2) = Ym(Ybg2) / cat_stat_nanmean(Ym(Ybg2(:))); clear Ybg2;
      %Ybg2 = (Yms<0.1 & Ybg) | smooth3(Yg./Ydiv > 1000)>0.5; 
    end
    [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
    for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
    for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
    Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
    Ywi  = cat_vol_smooth3X(Ywi,bsmooth); % highres data have may stronger inhomogeneities 
    Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
    Ym   = Ym./Ywi;
    WMt2 = roundx(cat_stat_nanmedian(Ym(Yg(:)<0.2 & Ym(:)>0.9)),rf); 
    Ywi  = Ywi * WMt2;
  else
    WMt2 = roundx(cat_stat_nanmedian(Ym(Yg(:)<0.2 & Ym(:)>0.9 & ~Ybg(:))),rf); 
    Ywi  = ones(size(Ym),'single') .* WMth3;
  end
  
  %% background update
  if zeroBG
    stime = cat_io_cmd('  Refine background','g5','',opt.verb,stime);
    Ybg  = ((Yg.*Ym)<cat_vol_smooth3X(Ym,2)*1.2) & Ym>0.2;
    Ybg  = Ybg & ~Ybgc;
    Ybg  = Ybg & ~isnan(Yg) & ~isnan(Ym); 
    [Ybg,resT2] = cat_vol_resize(single(Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
    Ybg  = Ybg>0.5;
    Ybg  = cat_vol_morph(Ybg,'lc',8);
    Ybg  = cat_vol_smooth3X(Ybg,2); 
    Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5 & (Ym<0.2 | isnan(Ym));
    Ybg  = cat_vol_morph(Ybg,'lo');
  end
  if ~debug, clear Ybgc Ybb1 Ybb3 bb; end
  
  %% second WM inhomogeneity with improved Yt with higher lower threshold (avoid CSF and less filtering)
  stime = cat_io_cmd('  Final correction','g5','',opt.verb,stime);
  Yt   = Ym>max(0,Yms*0.3) & (Ym-Yg)>0.2 & Ym<1.2 & Ym<Yms*2 & Yg<0.2 & Ydiv<0.2 & Ydiv>-0.6 & ...
         smooth3(Ym./(Yms+eps).*Yg.*Ydiv<-0.1)<0.1 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
       %%
  %Yt   = Yt | (~Ybg & Ym>0.1 & Ydiv./(Yg+eps)<0.5 & (Ym>0.3 & Yg>0.1 & Ydiv<0) | (~Ybg & Ym>0.6) & Ym<1.2 & Yg<0.1);
  Yt   = Yt & Ym>Yms*0.3 & Ym<Yms*1.2 & ~(-Ydiv.*Ym./max(eps,Yms)>0.3);
  Yt(smooth3(Yt)<0.5)=0;
  Yt   = Yt & cat_vol_morph(Yt & Yms>0.5 | (Yms>1.8),'l'); 
  %%
  if ~highBG
    Ywi2 = ( Ym .* Yt) ./ max(eps,Yt);
    % it would be nice to use futher regions, but as far as we did not know
    % their average intensity in relation to the bias field it not so easy
    if ~zeroBG
      Ybg2 = Ybg(:) & Yg(:)<(cat_stat_nanmean(Yg(Ybg(:))) + 2*(cat_stat_nanstd(Yg(Ybg(:))))); 
      Ywi2(Ybg2)   = Ym(Ybg2) / cat_stat_nanmean(Ym(Ybg2(:))); %clear Ybg2;
    else

      %Ybg2 = (Yms<0.1 & Ybg) | smooth3(Yg./Ydiv > 1000)>0.5; 
      %Yhht = -Ydiv.*Ym./Yms>0.2; 
      %Ywi2(Yhht)   = Ym(Yhht) / cat_stat_nanmean(Ym(Yhht(:))); %clear Yhht;
    end
    [Ywi2,resT2] = cat_vol_resize(Ywi2,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
    for i=1:1, Ywi2 = cat_vol_localstat(Ywi2,Ywi2>0,2,3); end % only one iteration!
    for i=1:4, Ywi2 = cat_vol_localstat(Ywi2,Ywi2>0,2,1); end
    Ywi2  = cat_vol_approx(Ywi2,'nn',resT2.vx_volr,2);
    Ywi2  = cat_vol_smooth3X(Ywi2,bsmooth); %.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
    Ywi2  = cat_vol_resize(Ywi2,'dereduceV',resT2);    
    Ywi   = Ywi2 .* Ywi; % both bias fields
    bias  = std(Ywi(:))/mean(Ywi(:)); 
  else
    Ywi  = ones(size(Ym),'single'); % .* WMt2;
    bias  = 0; 
  end
  
  %% BG inhomogeneity (important for normalization of the background noise)
  %[Ybc,Ygr,resT2] = cat_vol_resize({Ysrc./Ywi,Yg},'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*4,16,'meanm'); 
  %Ybc  = cat_vol_morph(Ybc<BGth/WMth*2 & Ygr<0.05,'lc',2);
  %Ybc  = cat_vol_resize(smooth3(Ybc),'dereduceV',resT2)>0.5; 
  %{
  if zeroBG
    stime = cat_io_cmd('  Background correction','g5','',opt.verb,stime);
    [Ybc,resT2] = cat_vol_resize(single(Ysrc .* Ybg),'reduceV',resT3.vx_volr,max(8,min(16,cat_stat_nanmean(resT3.vx_volr)*4)),8,'meanm'); 
    Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1);
    Ybc  = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the background! 
    Ybc  = cat_vol_smooth3X(Ybc,4);
    Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2) .* (~Ywi); 
  else
    Ybc  = ones(size(Ywi),'single')*BGth1;
  end
  %}
  
  %% back to original size
  %clear Ysrc; 
  stime = cat_io_cmd('  Final scaling','g5','',opt.verb,stime);
  Ywi       = cat_vol_resize(Ywi,'dereduceV',resT3);   
  Yg        = cat_vol_resize(Yg,'dereduceV',resT3); 
  Ydiv      = cat_vol_resize(Ydiv,'dereduceV',resT3); 
  Yms       = cat_vol_resize(Yms,'dereduceV',resT3); 
  Ymi       = cat_vol_resize(Ym,'dereduceV',resT3); 
  [Yt,Ybg]  = cat_vol_resize({single(Yt),single(Ybg)},'dereduceV',resT3); Yt = Yt>0.5; Ybg = Ybg>0.5;
  
  
  %% intensity normalization (Ybc is the average background noise)
  % in data with strong inhomogeneities (7T) the signal can trop below the noise level 
  if zeroBG, BGth = BGth2; else BGth = BGth1; end
  Ym   = ( (Ysrco - BGth) ./ (WMth3 - BGth)) ./ Ywi;
  noise = cat_stat_nanmean( Ym( Ybg(:) & (Yg(:)<0.5)) ); 
  
  if opt.iproc
    Ymw  = Ymi<Yms*0.5 & Yg<Ygth*(30*noise) & ~Ybg & Yg./abs(Ydiv)<min(3,max(0.5,15*noise)) & Yms>max(Yms(:))*0.5; 
  else
    Ymw  = Ymi<1.1 & Ymi>0.6 & Yg<Ygth*(30*noise) & ~Ybg & Yg./abs(Ydiv)<min(3,max(0.5,15*noise)) & Yms>max(Yms(:))*0.5; 
  end
  if sum(Ymw(:))>100000; Ymw = cat_vol_morph(Ymw,'lo',1); else Ymw = cat_vol_morph(Ymw,'l'); end
  Wth  = max(0.2,min(1.2,single(cat_stat_nanmedian(Ym(Ymw(:)))))); 
  if opt.iproc
    Ymw  = Ymi<Yms*0.5 & Yg<Ygth*(30*noise) & Ymi>Wth*0.3 & Ymi<Wth*1.5 & ~Ybg & Yg./abs(Ydiv)<min(3,max(0.5,15*noise)) & Yms>max(Yms(:))*0.5; 
  else
    Ymw  = Ymi<1.1 & Ymi>Wth*0.9 & Yg<Ygth*(30*noise) & ~Ybg & Yg./abs(Ydiv)<min(3,max(0.5,15*noise)) & Yms>max(Yms(:))*0.5; 
  end
  %if sum(Ymw(:))>100000; Ymw = cat_vol_morph(Ymw,'lo',1); else Ymw = cat_vol_morph(Ymw,'lo'); end
  [WIth3,WMv] = hist(Ym(Ymw(:)),1000);
  WItm = find(cumsum(WIth3)/sum(WIth3)>0.8,1,'first'); WItm = roundx(WMv(WItm),rf); 
  Ym   = Ym ./ WItm; 
  Yo   = ((Ysrco-BGth)/(WMth3-BGth))./WItm; 
  [WIth3,WMv] = hist(Yo(Ymw(:)),1000);
  WMto = find(cumsum(WIth3)/sum(WIth3)>0.7,1,'first'); WMto = roundx(WMv(WMto),rf); 
  Yo   = Yo / WMto;
  
  %% update WMth of the original image!
  Ysrco = Ysrco + Ysrcmin;
  [WIth3,WMv] = hist(Ysrco(Yg(:)<0.2 & Ym(:)>Wth*0.5 & Ym(:)<Wth*1.5 & ~Ybg(:)),1000);
  WMth = find(cumsum(WIth3)/sum(WIth3)>0.7,1,'first'); WMth = roundx(WMv(WMth),rf); 
  
  %% prepare intensity normalization by brain tissues
  [Ymr,Ytr,resT2] = cat_vol_resize({Ym,Yt},'reduceV',vx_vol,2,32,'meanm'); 
  Yb0r = cat_vol_morph(cat_vol_morph(Ytr>0.1,'dd',6,resT2.vx_volr),'ldc',10,resT2.vx_volr) & Ymr<1.2; 
  T3th = kmeans3D(Ymr(Yb0r(:)),5); T3th = T3th(1:2:5);
  clear Ymr Ytr Yb0r;

  
  %% check if correction was successful
  Ymx = Ymw; %cat_vol_morph(Ymw,'lo'); 
  rres = min(3,max(6,cat_stat_nanmean(vx_vol)*4)); 
  Ymm = Ym .* Ymx; Ymm = cat_vol_resize(Ymm,'reduceV',vx_vol,rres,msize,'meanm'); mm = Ymm(Ymm(:)>0);
  Yom = Yo .* Ymx; Yom = cat_vol_resize(Yom,'reduceV',vx_vol,rres,msize,'meanm'); mo = Yom(Yom(:)>0);
  Ymb = Ym .* Ymx; Ymb = cat_vol_resize(Ymb,'reduceV',vx_vol,rres,msize,'meanm'); bm = Ymb(Ymb(:)>0);
  Yob = Yo .* Ymx; Yob = cat_vol_resize(Yob,'reduceV',vx_vol,rres,msize,'meanm'); bo = Yob(Yob(:)>0);
  if opt.iproc
    Tth.biascorr = mean( (cat_stat_nanstd(mo)/cat_stat_nanmean(mo)) / (cat_stat_nanstd(mm)/cat_stat_nanmean(mm)) ); 
  else
    Tth.biascorr = mean([ (cat_stat_nanstd(mo)/cat_stat_nanmean(mo)) / (cat_stat_nanstd(mm)/cat_stat_nanmean(mm)) ...
                      (cat_stat_nanstd(bo)/cat_stat_nanmean(bo)) / (cat_stat_nanstd(bm)/cat_stat_nanmean(bm)) ]); 
  end
  if ~debug, clear Ymw Ymm Yom Ymb Yob; end 

  if opt.iproc
    Tth.inverse  = 1; 
  else
    Tth.inverse  = max(0,min(1,1-(Wth-0.6)/0.4));
  end
  Tth.Tmax     = WMth/WMth0;
  cat_io_cmd('','','',opt.verb,stime); 
  
  % PQA variable
  % - input/default parameter
  pior.ipara          = opt;            
  % - processing parameter
  pior.ppara.vx_vol   = resT3.vx_vol;   % voxelsize input image
  pior.ppara.vx_volr  = resT3.vx_volr;  % voxelsize for processing
  pior.ppara.bsmooth  = bsmooth;    % smoothing factor of bias field 
  pior.ppara.nsmooth  = nsmooth;        % smoothing factor to reduce noise (filtersize in voxel)
  pior.ppara.inverse  = Tth.inverse;    % contrast
  pior.ppara.zeroBG   = zeroBG;         % background
  pior.ppara.Tmax     = Tth.Tmax;       % 
  pior.ppara.BGth1    = BGth1; 
  pior.ppara.BGth2    = BGth2; 
  pior.ppara.WMth3    = WMth3; 
  % - test parameter
  pior.tpara.CVin     = mean([cat_stat_nanstd(mo)/cat_stat_nanmean(mo) cat_stat_nanstd(bo)/cat_stat_nanmean(bo)]);
  pior.tpara.CVout    = mean([cat_stat_nanstd(mm)/cat_stat_nanmean(mm) cat_stat_nanstd(bm)/cat_stat_nanmean(bm)]); 
  pior.tpara.CVR      = Tth.biascorr; 
  if opt.iproc
    pior.tpara.noise    = noise/18 * 125 / Wth;  % BWP correction
    pior.tpara.bias     = bias*6   * 250 / Wth;  % BWP correction
  else
    pior.tpara.noise    = noise/3 * 125 / Wth / (1+2*Tth.inverse);  % BWP correction
    pior.tpara.bias     = bias    * 250 / Wth / (1+2*Tth.inverse);  % BWP correction
  end
  if opt.verb
    cat_io_cprintf([0.5 0.5 0.5],sprintf('    APP-noise ~ %0.0f%%%%, APP-bias ~ %0.0f%%%% (BWP-level), APP-inv ~ %0.2f\n',...
      pior.tpara.noise,pior.tpara.bias,pior.ppara.inverse)); 
  end
  
  %% value greater 1 describe a successful correction
  if Tth.biascorr<1.01 || Tth.inverse>0.5
    
    % try inverse processing 
    if opt.icall && Tth.inverse %&& Tth.biascorr<1.05
      if opt.verb>1
        cat_io_cprintf('warn',sprintf('    T1-Bias correction failed (CJVR=%0.2f), try inverse correction: ',Tth.biascorr)); 
      end
      
      % second call with inverted image (due to maximum correction theorem)  
      Ysrcmax = max(Ysrco(~isinf(Ysrco(:)))); 
      opti = opt; opti.icall=0; opti.iproc=1; 
      [Ym2,Yt2,Ybg2,WMth3,bias2,Tth2,piorapp2] = cat_run_job_APP_init(Ysrcmax - Ysrco,vx_vol,opti); % no further recursion
    
      if Tth2.biascorr>=1.01
        Ym = (Ym2 - min(Ym2(:))) ./ (1 - min(Ym2(:))); Ym(Ybg2) = 0; % inverted map (T2/PD > T1
        Ym = Ym ./ mean(Ym(Ymx(:))); % well inverse map need another scaling otherwise we cut the intensity in affreg
        %Ym = (max(Ym2(Ybg2(:))) - Ym2) ./ ([max(Ym2(Ybg2(:))) - 1); % original map (same weighting as input)
        pior = piorapp2;  pior.ppara.inverse = Tth.inverse;
        Tth = Tth2;       Tth.inverse = pior.ppara.inverse;
        Yt = Yt2; Ybg = Ybg2; WMth = Ysrcmax - WMth3; bias = bias2; 
      else
        Ym = Yo;       
      end
    else
      if opt.verb>1
        cat_io_cprintf('warn',sprintf('    Bias correction failed use only scaling (CJVR=%0.2f). \n',Tth.biascorr)); 
      end
    
      Ym = Yo;
    end
  else    
    if opt.verb>1
      cat_io_cprintf([0 0.5 0],sprintf('    Bias correction successful (CJVR=%0.2f). \n',Tth.biascorr)); 
    end
  end
  
  %% intensity normalization
  [Ymr,Ytr,resT2] = cat_vol_resize({Ym,Yt},'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Yb0r = cat_vol_morph(cat_vol_morph(Ytr>0.1,'dd',6,resT2.vx_volr),'ldc',10,resT2.vx_volr) & Ymr<1.2; 
  T3th = kmeans3D(Ymr(Yb0r(:)),5); T3th = T3th(1:2:5);
  T3th2 = T3th; 
  if 1 % highBG
    T3thc = kmeans3D(Ymr(Yb0r(:) & ...
      Ymr(:)<cat_stat_nansum(T3th(1:2).*[0.8 0.2]) ),3); % close to minimum
    T3th2(1) = T3thc(1); T3th(1) = T3thc(1); 
    T3th2(2) = kmeans3D(Ymr(Yb0r(:) & ...
      Ymr(:)<cat_stat_nansum(T3th(2:3).*[0.8 0.2]) & ...
      Ymr(:)>cat_stat_nansum(T3th(1:2).*[0.8 0.2]) ),1);
    T3th2(3) = kmeans3D(Ymr(Yb0r(:) & ...
      Ymr(:)>cat_stat_nansum(T3th(2:3).*[0.5 0.5]) & ...
      Ymr(:)<(T3th(3) + diff(T3th(2:3))) ),1);
  end
  T3th = T3th2; 
  clear Ymr Ytr Yb0r;
  
  Tth.T3thx  = [0 T3th T3th(3)+diff(T3th(2:3))];
  Tth.T3th   = 0:1/3:4/3;
  Ym = cat_main_gintnormi(Ym/3,Tth);
  
  if opt.verb>1 && opt.icall, cat_io_cmd(' ','','',opt.verb); end
      
end
%=======================================================================
function r = roundx(r,rf)
  r(:) = round(r(:) * rf) / rf;
end
%=======================================================================
