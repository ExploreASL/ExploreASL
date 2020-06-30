function [Ym,Yt,Ybg,WMth,bias] = cat_run_job_APP_init1070(Ysrco,vx_vol,verb)
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
%  $Id: cat_run_job_APP_init1070.m 1414 2019-01-14 08:16:00Z dahnke $


%    ds('l2','',0.5,Yo/WMth,Yg<0.2,Yo/WMth,Ym,80)

  rf = 10^9; 
  bfsmoothness = 3; 
  if verb, fprintf('\n'); end
  
  stime = cat_io_cmd('  Initialize','g5','',verb);
  msize = 222; %round(222 ./ max(size(Ysrco).*vx_vol) .* min(size(Ysrco).*vx_vol));  

  [Ysrc,resT3] = cat_vol_resize(Ysrco,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*2),msize,'meanm'); 
  vx_vol = resT3.vx_volr; 
  
  %% correction for negative backgrounds (MT weighting)
  [Ym,BGth] = cat_stat_histth(Ysrc,99.99); BGth(2) = []; clear Ym;   %#ok<ASGLU>
  
  Ysrc = Ysrc - BGth; Ysrco = Ysrco - BGth; BGth2 = BGth; 
  Yg   = cat_vol_grad(Ysrc,vx_vol) ./ max(eps,Ysrc); 
  Ygth = max(0.05,min(0.2,mean(Yg(Yg(:)>0))*0.75));
  
  Ybg0 = cat_vol_smooth3X(Yg<Ygth & Ysrc<cat_stat_nanmean([BGth cat_stat_nanmean(Ysrc(:))]),8)>0.05;
  Yw0  = Yg<Ygth & Ysrc>cat_stat_nanmean(  Ysrc(Yg(:)<Ygth & ~Ybg0(:))) & ~Ybg0; 
  WMth = roundx(single(cat_stat_nanmedian(Ysrc( Yw0(:) ))),rf); clear Ybg0 Yw0;
  BGth = max( min(Ysrc(:))*0.7 + 0.3*cat_stat_nanmean(Ysrc(:)) ,...
    cat_stat_nanmean(Ysrc(Ysrc(:)<cat_stat_nanmean(Ysrc(:))))); BGth = roundx(BGth,rf); 
  highBG = BGth/WMth>0.3; % MT / MP2RAGE
  if highBG 
    % In case of high intensity background we simply need the lowest 
    % intensity of the histogram without extrem outliers.
    Ysrcr = cat_vol_resize(Ysrc,'reduceV',vx_vol,2,24,'meanm'); 
    [Ysrcr,th] = cat_stat_histth(Ysrcr,99.99); %#ok<ASGLU>
    BGth = th(1);
    clear Ysrcr;
  end
  Ym   = (Ysrc - BGth) ./ (WMth - BGth);
  
  Ydiv = cat_vol_div(Ym,vx_vol/2) ./ (Ym+eps); % lower resolution is 8 times faster 
  
  
  %% background
  stime = cat_io_cmd('  Estimate background','g5','',verb,stime);
  if highBG
    % As far as our background have an intensity similar to the brain, we 
    % have to find the average intensity of the background create a mask.
    Ymr   = cat_vol_resize(Ysrc,'reduceV',vx_vol,3,32,'meanm'); 
    Ymsk  = true(size(Ymr)); Ymsk(5:end-4,5:end-4,5:end-4) = false;
    [HBGth,HBGsd] = kmeans3D(Ymr(Ymsk(:)),1); 
    HBGsd = max(diff([BGth,WMth])*0.1, HBGsd); clear Ymr
    Ymsk  = true(size(Ym)); Ymsk(5:end-4,5:end-4,5:end-4) = false;
    Ybg   = (Yg<0.3 & Ysrc>(HBGth - 2*HBGsd) & Ysrc<(HBGth + 2*HBGsd)) | Ymsk; 
    [Ybg,resT2] = cat_vol_resize(smooth3(Ybg) ,'reduceV',vx_vol,2,32,'meanm'); 
    Ybg   = cat_vol_morph(cat_vol_morph(Ybg>0.5,'l',[100,400])<0.5,'lc',4)<0.5;
    Ybgid = cat_vbdist(single(Ybg)); Ybgid = Ybgid./max(Ybgid(Ybgid(:)<1000));
    Ybgi  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5; % inverse background = object
    Ybgid = cat_vol_resize(Ybgid,'dereduceV',resT2); % inverse background = object
    
    %% estimate WM threshold
    Yw0  = Yg<0.15 & Ysrc>cat_stat_nanmean( Ysrc(Yg(:)<0.15 & Ybgid(:)>0.5)) & Ybgid>0.8; 
    Yw0  = cat_vol_morph(smooth3(Yw0)>0.5,'l',[10,100])>0; 
    Yw0  = cat_vol_morph(smooth3(Yw0)>0.5,'l'); 
    WMth = kmeans3D(Ysrc(Yw0(:)),1); 
  else
    Ybgi = ((Yg.*Ym)<cat_vol_smooth3X(Ym,2)*1.2) & Ym>0.1; % inverse background = object
    Ybgi([1,end],:,:)=0; Ybgi(:,[1,end],:)=0; Ybgi(:,:,[1,end])=0; Ybgi(smooth3(Ybgi)<0.5)=0;
  end
  [Ybg,resT2] = cat_vol_resize(single(Ybgi),'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',8);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5; 
  if ~highBG
    BGth  = roundx(cat_stat_nanmean(Ysrc(Ybg(:))),rf);
    [Ybgr,resT2] = cat_vol_resize(smooth3(Ybg) ,'reduceV',vx_vol,2,32,'meanm'); 
    Ybgid = cat_vbdist(single(Ybgr)); Ybgid = Ybgid./max(Ybgid(Ybgid(:)<1000));
    Ybgid = cat_vol_resize(Ybgid,'dereduceV',resT2);
  end
  Ym   = (Ysrc - BGth) ./ (WMth - BGth);
  
  %% first WM inhomogeneity with low tissue boundary (may include CSF > strong filtering for IXI175)
  stime = cat_io_cmd('  Initial correction','g5','',verb,stime);
  Yms  = cat_vol_smooth3X( min(2 .* ~Ybg,Ym .* (Ydiv>-0.2) .* ~Ybg .* (Ym>0.1)),16*mean(vx_vol));     % this map is to avoid CSF in the mask!
  Yms  = (Yms ./ mean(Yms(~Ybg(:)))) * WMth;
  Yms  = cat_vol_smooth3X( min(Yms*1.5 .* ~Ybg,Ysrc .* ~Ybg),16*mean(vx_vol));
  Yms  = (Yms ./ mean(Yms(~Ybg(:)))) * WMth;
  if ~highBG
    Yt   = Ysrc>max(BGth,Yms*0.3) & Ysrc<Yms*2 & Ysrc<WMth*(1+Yms/WMth*2) & Yg<0.9 & Ydiv<0.2 & ...
           Ydiv>-0.6 & smooth3(Ysrc./(Yms+eps).*Yg.*Ydiv<-0.2)<0.3 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
    Ywi  = (Ysrc .* Yt) ./ max(eps,Yt);  
    [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
    for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
    for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
    Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
    Ywi  = cat_vol_smooth3X(Ywi,bfsmoothness.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
    Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);  
    Ybc  = (Ysrc - BGth) ./ Ywi;
    WMt2 = roundx(cat_stat_nanmedian(Ybc(Yg(:)<0.2 & Ym(:)>0.9)),rf); clear Ybc;
    Ywi  = Ywi * WMt2;
    Ym   = (Ysrc - BGth) ./ Ywi;
  else
    Yt   = Ysrc>max(BGth,Yms*0.3) & Ysrc<Yms*2 & Ysrc<WMth*(1+Yms/WMth*2) & Yg<0.9 & Ydiv<0.2 & Ybgid>0.5 & ...
           Ym>0.5 & Ym<1.1 & ...
           Ydiv>-0.6 & smooth3(Ysrc./(Yms+eps).*Yg.*Ydiv<-0.2)<0.3 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
    Yt   = cat_vol_morph(Yt,'ldo',1.5);
    Ywi  = ones(size(Ym),'single');
  end
  
  %% background update
  if ~highBG
    stime = cat_io_cmd('  Refine background','g5','',verb,stime);
    Ybg = ((Yg.*(Ym))<cat_vol_smooth3X(Ym,2)*1.2) & Ym>0.2; 
    Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg(smooth3(Ybg)<0.5)=0;
    [Ybg,resT2] = cat_vol_resize(single(Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
    Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg = Ybg>0.5;
    Ybg  = cat_vol_morph(Ybg,'ldc',8);
    Ybg  = cat_vol_smooth3X(Ybg,2); 
    Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5 & Ym<min(WMth*0.2,BGth*0.8+0.2*WMth);
    Ybg  = cat_vol_morph(Ybg,'lo');
  end

  %% second WM inhomogeneity with improved Yt with higher lower threshold (avoid CSF and less filtering)
  stime = cat_io_cmd('  Final correction','g5','',verb,stime);
  Yt   = Ysrc>max(BGth,Yms*0.3)  & Ym>0.9 & Ym<1.2  & Ybgid>0.5 & ...
         Ym<Yms/WMth*2 & Yg<0.9 & Ydiv<0.2 & Ydiv>-0.6 & ...
         smooth3(Ym.*Yg.*Ydiv<-0.1)<0.1 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  Yt   = Yt | (~Ybg & Ysrc>BGth/2 & Ysrc>Yms*0.5 & Ysrc<Yms*1.2 & Ydiv./(Yg+eps)<0.5 & ((Ym>0.3 & Yg>0.1 & Ydiv<0) | (~Ybg & Ysrc./(Ywi+eps)>0.6)) & Ysrc./(Ywi+eps)<1.2); 
  Yt(smooth3(Yt)<0.7)=0;
  if ~highBG
     Ywi  = (Ysrc .* Yt) ./ max(eps,Yt);  
    [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
    for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
    for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
    Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
    Ywi  = cat_vol_smooth3X(Ywi,bfsmoothness.*mean(vx_vol)); %.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
    Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
    Ybc  = (Ysrc - BGth) ./ Ywi;
    WMt2 = roundx(cat_stat_nanmedian(Ybc(Yg(:)<0.2 & Ym(:)>0.9)),rf); 
    Ywi  = Ywi * WMt2;
    Ym   = (Ysrc - BGth) ./ Ywi;
    bias = std(Ywi(:))/mean(Ywi(:)); 
  else 
    Ywi  = ones(size(Ym),'single') .* WMth;
    Yt   = cat_vol_morph(Yt,'lo');
    bias = 0;
  end
  
  %% BG inhomogeneity (important for normalization of the background noise)
  %[Ybc,Ygr,resT2] = cat_vol_resize({Ysrc./Ywi,Yg},'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*4,16,'meanm'); 
  %Ybc  = cat_vol_morph(Ybc<BGth/WMth*2 & Ygr<0.05,'lc',2);
  %Ybc  = cat_vol_resize(smooth3(Ybc),'dereduceV',resT2)>0.5; 
  if 0
    stime = cat_io_cmd('  Background correction','g5','',verb,stime);
    [Ybc,resT2] = cat_vol_resize(single(Ysrc .* Ybg),'reduceV',resT3.vx_volr,max(8,min(16,cat_stat_nanmean(resT3.vx_volr)*4)),16,'min'); 
    Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,2);
    Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1);
    %Ybc  = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the background! 
    Ybc  = cat_vol_smooth3X(Ybc,4);
    Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2) * WM; 
  end
  
  %% prepare intensity normalization by brain tissues
  [Ymr,Ytr,resT2] = cat_vol_resize({Ym,Yt},'reduceV',resT3.vx_volr,3,32,'meanm'); 
  Yb0r = cat_vol_morph(cat_vol_morph(Ytr>0.1,'dd',6,resT2.vx_volr),'ldc',10,resT2.vx_volr) & Ymr<1.2; 
  T3th = kmeans3D(Ymr(Yb0r(:)),5); T3th = T3th(1:2:5);
  clear Ymr Ytr Yb0r;

  %% back to original size
  stime = cat_io_cmd('  Final scaling','g5','',verb,stime);
  Ywi       = cat_vol_resize(Ywi,'dereduceV',resT3); 
  %Ybc       = cat_vol_resize(Ybc,'dereduceV',resT3); 
  Yg        = cat_vol_resize(Yg,'dereduceV',resT3); 
  [Yt,Ybg]  = cat_vol_resize({single(Yt),single(Ybg)},'dereduceV',resT3); Yt = Yt>0.5; Ybg = Ybg>0.5;
  Ysrc      = Ysrco; clear Ysrco;

  %% intensity normalization (Ybc is the average background noise)
  % in data with strong inhomogeneities (7T) the signal can trop below the noise level 
  Ym   = (Ysrc - BGth) ./ Ywi; %(Ywi - min(BGth,min(Ybc/2,Ywi/20))); 
  Wth  = single(cat_stat_nanmedian(Ym(Yg(:)<0.2 & Yt(:)))); 
  [WIth,WMv] = hist(Ym(Yg(:)<0.2 & Ym(:)>Wth*0.5 & Ym(:)<Wth*1.5 & ~Ybg(:)),0:0.01:2);
  WIth = find(cumsum(WIth)/sum(WIth)>0.8,1,'first'); WIth = roundx(WMv(WIth),rf); 
  Ym   = Ym ./ WIth; 
  % update WMth
  Ysrc = Ysrc + BGth2;
  [WIth,WMv] = hist(Ysrc(Yg(:)<0.2 & Ym(:)>Wth*0.5 & Ym(:)<Wth*1.5),1000);
  WMth = find(cumsum(WIth)/sum(WIth)>0.7,1,'first'); WMth = roundx(WMv(WMth),rf); 
  
  %% intensity normalization
  [Ymr,Ytr,resT2] = cat_vol_resize({Ym,Yt},'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Yb0r = cat_vol_morph(cat_vol_morph(Ytr>0.1,'dd',6,resT2.vx_volr),'ldc',10,resT2.vx_volr) & Ymr<1.2; 
  T3th = kmeans3D(Ymr(Yb0r(:)),5); T3th = T3th(1:2:5);
  T3th2 = T3th; 
  if 1 % highBG
    T3thc = kmeans3D(Ymr(Yb0r & Ymr>0.1 & ...
      cat_vol_morph(Ymr<cat_stat_nansum(T3th(1:2).*[0.8 0.2]) ,'de',1)),3); % close to minimum
    T3th2(1) = T3thc(1); T3th(1) = T3thc(1); 
    T3th2g = kmeans3D(Ymr(Yb0r(:) & ...
      Ymr(:)<cat_stat_nansum(T3th(2:3).*[0.8 0.2]) & ...
      Ymr(:)>cat_stat_nansum(T3th(1:2).*[0.8 0.2]) ),3);
    T3th2(2) = T3th2g(2);
    T3th2w = kmeans3D(Ymr(Yb0r & ...
      cat_vol_morph(Ymr>cat_stat_nansum(T3th(2:3).*[0.5 0.5]) & Ymr<(T3th(3) + diff(T3th(2:3))),'de',1) ),3);
    T3th2(3) = T3th2w(2); 
  end
  T3th = T3th2; 
  clear Ymr Ytr Yb0r;

  Tth.T3thx  = [0 T3th T3th(3)+diff(T3th(2:3))];
  Tth.T3th   = 0:1/3:4/3;
  Ym = cat_main_gintnormi(Ym/3,Tth);
  
  cat_io_cmd(' ','','',verb,stime); 
end
%=======================================================================
function r = roundx(r,rf)
  r(:) = round(r(:) * rf) / rf;
end
%=======================================================================
