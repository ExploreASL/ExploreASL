function [Yml,Ymg,Ycls,Ycls2,T3th] = cat_main_LAS1585(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,extopts,Tth) 
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
%
% Local Adaptive Segmentation (LAS):
%
% This version of the local adaptive intensity correction includes a  
% bias correction that based on a maximum filter for the WM and a mean
% filter of the GM to stabilize the correction in region with less WM.
%
% The extension based mostly on the assumption that the tissue next to 
% the CSF (and high divergence sulci) has to be WM (maximum, high 
% divergence) or GM. For each tissue a refined logical map is generated 
% and used to estimate the local intensity threshold.
%
% It is important to avoid high intensity blood vessels in the process, 
% because they will push down local WM and GM intensity. 
%
% There are further regionwise correction, e.g. , to avoid overfitting in 
% cerebellum, or adapt for age specific changes, e.g. enlarged ventricle.
%
% Based on this values a intensity transformation is used. Compared to 
% the global correciton this has to be done for each voxel. To save time
% only a rough linear transformation is used.
% ______________________________________________________________________
%
%   [Yml,Ycls,Ycls2,T3th] = ...
%     cat_main_LAS1585(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,PA,template)
%
%   Yml    .. local intensity correct image
%   Ycls   .. corrected SPM tissue class map
%   Ycls2  .. ?
%   T3th   .. tissue thresholds of CSF, GM, and WM in Ysrc
%
%   Ysrc   .. (bias corrected) T1 image
%   Ym     .. intensity corrected T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%   Ycls   .. SPM tissue class map
%   Yb0    .. brain mask
%   Yy     .. deformation map
%   res    .. SPM segmentation structure
%   vx_vol .. voxel dimensions
%   PA     .. CAT atlas map
%   template .. ?
% ______________________________________________________________________
% 
% internal maps:
%
%   Yg   .. gradient map   - edges between tissues
%   Ydiv .. divergence map - sulci, gyris pattern, and blood vessels
%   Yp0  .. label map      - tissue classes (BG=0,CSF=1,GM=2,WM=3) 
%
%   Ysw  .. save WM tissue map
%   Ybv  .. blood vessel map
%   Ycp  .. CSF / background area for distances estimation
%   Ycd  .. CSF / background distance
%   Ycm  .. CSF
%   Ygm  .. GM
%   Ywm  .. WM 
%   Yvt  .. WM next to the ventricle map 
%   Ygmt .. cortical thickness map
%   Ypp  .. cortical percentage position map
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_LAS.m 1036 2016-10-18 14:26:32Z dahnke $

  % set this variable to 1 for simpler debuging without reduceBrain
  % function (that normally save half of processing time)
  debug   = extopts.debug; % debug = 1;
  verb    = extopts.verb-1;
  vxv     = 1/ mean(vx_vol);
  dsize   = size(Ysrc);
  NS      = @(Ys,s) Ys==s | Ys==s+1;             % function to ignore brain hemisphere coding
  LASstr  = max(eps,min(1,extopts.LASstr));      % LAS strenght (for GM/WM threshold)3
  LAB     = extopts.LAB;                         % atlas labels
  cleanupstr  = min(1,max(0,extopts.gcutstr));   % required to avoid critical regions
  cleanupdist = min(3,max(1,1 + 2*cleanupstr));
    
  
  
%% ---------------------------------------------------------------------
%  First, we have to optimize the segments using further information that 
%  SPM do not use, such as the gradient, divergence and distance maps. 
%  The gradient map (average of the first derivate of the T1 map) is an 
%  edge map and independent of the image intensity. It helps to avoid PVE 
%  regions and meninges. 
%  The divergence (second derivate of the T1 map) help to identfiy sulcal
%  and gyral pattern and therefore to find WM and CSF regions for furhter 
%  corrections and to avoid meninges and blood vessels. 
%  Furhtermore, special assumption can be used. 
%  The first one is the maximum property of the WM in T1 data that allows
%  using of a maxim filter for the GM/WM region. 
%  The second is the relative stable estimation of CSF/BG that allows to 
%  estimat a distance map. Because, most regions have a thin layer of 
%  GM around the WM we can avoid overestimation of the WM by the other 
%  maps (especially the divergence). 
%  ---------------------------------------------------------------------
  fprintf('\n');
  stime = cat_io_cmd('  Prepare maps','g5','',verb); dispc=1;

  
  % brain segmentation can be restricted to the brain to save time 
  [Ym,Yb,BB] = cat_vol_resize({Ym,Yb0},'reduceBrain',vx_vol,round(10/mean(vx_vol)),Yb0);
  Yclsr=cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceBrain',vx_vol,BB.BB); end
  
  
  % help maps to detect edges (Yg) and sulci/gyris (Ydiv)
  Yg    = cat_vol_grad(Ym,vx_vol);                                                  % mean gradient map
  Ydiv  = cat_vol_div(max(0.33,Ym),vx_vol);                                         % divergence map 
  Yp0   = single(Yclsr{1})/255*2 + single(Yclsr{2})/255*3 + single(Yclsr{3})/255;   % tissue label map
  Yb    = smooth3(Yb | (cat_vol_morph(Yb,'d',2*vxv) & Ym<0.8 & Yg<0.3 & Ym>0 & Yp0>0.2))>0.5;   % increase brain mask, for missing GM 
  
  
  
  %% adding of atlas information (for subcortical structures)
  %  -------------------------------------------------------------------
  stime = cat_io_cmd('  Prepare partitions','g5','',verb,stime); dispc=dispc+1;

  % map atlas to RAW space
  for i=1:5
    try
      Vl1A = spm_vol(extopts.cat12atlas{1});
      break
    catch 
      % read error in parallel processing
      pause(1)
    end
  end
  Yl1  = cat_vol_ctype(round(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
  Yl1  = reshape(Yl1,dsize);

  % load WM of the TPM or Dartel/Shooting Template for WMHs
  Vtemplate = spm_vol(extopts.templates{end}); 
  Ywtpm = cat_vol_ctype(spm_sample_vol(Vtemplate(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
  Ywtpm = reshape(Ywtpm,dsize); spm_smooth(Ywtpm,Ywtpm,2*vxv);
  Ywtpm = single(Ywtpm)/255;

  % brain segmentation can be restricted to the brain to save time  
  Yl1    = cat_vol_resize(Yl1  ,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
  Ywtpm  = cat_vol_resize(Ywtpm,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
  if ~debug, clear Yy; end
  
  
  %% adaption of the LASstr depending on average basal values 
  LASmod = min(2,max(0,mean((Ym( NS(Yl1,LAB.BG) & Yg<0.1 & Ydiv>-0.05  & Yclsr{1}>4)) - 2/3) * 8));
  LASstr  = min(1,max(0.05,LASstr * LASmod)); clear LASmod                 % adaption by local BG variation
  LASfs   = 1 / max(0.05,LASstr);                                          % smoothing filter strength 
  LASi    = min(8,round(LASfs));                                           % smoothing interation (limited)
   
  
  %% GM thickness (Ygmt) and percentage possition map (Ypp) estimation
  %  -------------------------------------------------------------------
  %  The Ypp and Ygmt maps are used to refine the GM especially to correct
  %  highly myelinated GM regions.
  %  Interpolation look a little bit better, but I am not sure if it is
  %  necessary, so we do this only for lower/average resolutions. 
  %  It is further unclear who to handle subcortical regions ...
  if 0
    PBTinterpol = 0*mean(vx_vol)>0.9;
    tic
    if PBTinterpol
      Yli     = interp3(Yl1,1,'nearest'); 
      Ymi     = interp3(Ym,1); 
      Ybi     = interp3(single(Yb),1)>0.5;
      Ywtpmi  = interp3(Ywtpm,1); 
      Ydivi   = interp3(Ydiv,1);
      vx_voli = vx_vol/2; 
    else
      Yli     = Yl1; 
      Ymi     = cat_vol_smooth3X(Ym,0.5/mean(vx_vol)); 
      Ybi     = Yb; 
      Ywtpmi  = Ywtpm; 
      Ydivi   = Ydiv; 
      vx_voli = vx_vol; 
    end

    % PBT thickness and percentage position estimation 
    Ybgc  = cat_vol_smooth3X(Ymi>0.7 & Ymi<0.95 & (NS(Yli,LAB.BG) | NS(Yli,LAB.TH)),2/mean(vx_voli)); % correction for subcortical structures
    Ybgc  = (Ymi-2/3) .* (Ymi>2/3 & Ybgc>0.3);
    Ycsfd = cat_vbdist(2-Ybi-Ymi,(Ymi-Ybgc)<2.5/3);
    Ywmh  = cat_vol_morph((((Ymi-Ybgc+((Ywtpmi - Ydivi*0.1)>0.98)).*Ybi)*3-2)>0.5,'lc',1/mean(vx_voli)); % WM hyperintensity 
    clear Ywtpmi Ydivi;
    Ywmd  = cat_vbdist(((Ymi-Ybgc+Ywmh).*Ybi)*3-2,Ybi);
    Ygmt  = cat_vol_pbtp((Ymi-Ybgc+Ywmh).*Ybi*3,Ywmd,Ycsfd); 
    for i=1:1, Ygmt = cat_vol_localstat(Ygmt,Ygmt>0,1,1); end
    Ypp = zeros(size(Ybi)); 
    Ypp(Ygmt>0) = min(Ybi(Ygmt>0),min(Ycsfd(Ygmt>0),Ygmt(Ygmt>0)-Ywmd(Ygmt>0))./max(eps,Ygmt(Ygmt>0))); 
    Ypp((Ymi>5/6 & Ybi & ~(NS(Yli,LAB.BG) | NS(Yli,LAB.TH))) | (Ymi>0.95 & Ybi)) = 1;
    Ypp((((Ymi-Ybgc+Ywmh).*Ybi)*3-2)>0.5)=1;
    Ypp = cat_vol_smooth3X(Ypp,0.5/mean(vx_voli));

    % thickness correction and back to original resolution
    Ygmt = Ygmt*mean(vx_vol); Ygmt(Ygmt>10 | isnan(Ygmt)) = 0;
    clear Ycsfd Ywmd Ymi Ybi Yli vx_voli
    if PBTinterpol
      Ygmt  = cat_vol_resize(Ygmt,'reduceV',1,2,10,'meanm');
      Ypp   = cat_vol_resize(Ypp,'reduceV',1,2,10,'meanm');
      Ywmh  = cat_vol_resize(Ywmh,'reduceV',1,2,10,'meanm');
    end
    [D,I] = cat_vbdist(single(Ygmt>eps),Yp0>0); Ygmt = Ygmt(I); clear D I;  

    % subcortical structures
    %Ypp( NS(Yl1,LAB.BG) & Ym>2.1/3 & Ym<2.9/3 )=0.5;
    %Ypp   = cat_vol_median3(Ypp,Ypp<0.1 & Ym>=1.9/3 & Yb,true(size(Ym))); 

    Ywmh = Ywmh.*0;
  else
    [Ygmt,Ypp] = cat_vol_pbt1585( (Yp0 + (Ym*3 .* (Yp0>0)))/2 );
    Ygmt = Ygmt*mean(vx_vol); Ygmt(Ygmt>10 | isnan(Ygmt)) = 0;
    [D,I] = cat_vbdist(single(Ygmt>eps),Yp0>0); Ygmt = Ygmt(I); clear D I;  
  end
  
  %% helping segments
  %  -------------------------------------------------------------------
  stime = cat_io_cmd('  Prepare segments','g5','',verb,stime); dispc=dispc+1;
  % Ybb .. don't trust SPM to much by using Yp0 because it may miss some areas! Shood be better now with MRF.
  Ybb = cat_vol_morph((Yb & Ym>1.5/3 & Ydiv<0.05) | Yp0>1.5,'lo',vxv);
  
  % Ysw .. save WM and blood vessels mpas
  % Ybv .. possible blood vessels
  Ysw = cat_vol_morph(Yclsr{2}>128 & (min(1,Ym)-Ydiv)<1.5,'lc',vxv*2) & (Ym-Ydiv)>5/6; % 1.2 
  Ybv = ((min(1,Ym) - Ydiv + Yg)>2.0 | (Yclsr{5}>16 & Ym<0.6 & Yclsr{1}<192)) & ...
        ~cat_vol_morph(Ysw,'d',1) & Ym>0.2; 
      
  % Ycp .. for CSF/BG distance initialization 
  Ycp = (Yclsr{3}>240 & Ydiv>0 & Yp0<1.1 & Ym<0.5) | ...                   % typcial CSF
        (Yclsr{5}>8 & Yclsr{2}<32 & Ym<0.6 & Ydiv>0) | ...                 % venes
        ((Ym-Ydiv/4<0.4) & Yclsr{3}>4 & Yclsr{3}>16) | ...                 % sulcal CSF
        (single(Yclsr{6})+single(Yclsr{5})+single(Yclsr{4}))>192 | ...     % save non-csf 
        ~cat_vol_morph(Ybb,'lc',5) | ...                                   % add background
        Ym<0.3 | Ypp<min(0.5,0.1/Ygmt & Yl1<3);                            % but do not trust the brain mask!
  Ywd = cat_vbdist(single(Yp0>2.5),Yp0>0.5,vx_vol);                        % WM distance for skelelton
  Ycp(smooth3(Ycp)>0.4)=1;                                                 % remove some meninges
  Ycd = cat_vbdist(single(Ycp),~Ycp,vx_vol);                               % real CSF distance 
  Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Yclsr{3}>4 & Yclsr{3}>1) = ...            % correction for sulci 
    min(Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Yclsr{3}>4 & Yclsr{3}>1),1.5);      % maybe a second distance estimation???
  % we need to remove strong edge regions, because here is no GM layer between CSF and WM ???  
  % Yb  = cat_vol_morph(~Ycp | (Ycls{3}>128),'lc',1);
  Ybd = cat_vbdist(single(~Yb),Yb,vx_vol);
  Yvt = (Yg+abs(Ydiv))>0.4 & smooth3(single(Yclsr{1})/255)<0.5 & Ybd>20 & ...
    cat_vol_morph(Yclsr{3}>8,'d',vxv) & cat_vol_morph(Yclsr{2}>8,'d',vxv); 
  Yvt = smooth3(Yvt)>0.7;
  Yvt = smooth3(Yvt)>0.2;
  
  
  
  %% final tissue maps:  Ycm = CSF, Ygm = GM, Ywm = WM 
  %  -------------------------------------------------------------------
  Ysc = Ycp & Yb & Yclsr{3}>192 & ~Ybv & Ym<0.45 & Yg<0.1;
  Ycm = Ycp & Yb & Yclsr{3}>192 & ~Ybv & (Yb | Ym>1/6) & Ym<0.45 & Yg<0.25 & Ym>0; % & Ydiv>-0.05;
  %Ycm = Ycm | (Yb & (Ym-max(0,Ydiv))<0.5); 
  Ywm = (Ysw | Yclsr{2}>252 | ((Ycd-Ydiv)>2 & Ydiv<0 & (Ym-Ydiv)>(0.9+LASstr*0.05) & Yb) | ... % save WM 
        ((Ycd-Ydiv.*Ycd)>4 & (Ydiv<-0.01) & Yb & Ym>0.5 & Ybd<20 & Ycd>2) | ...
        (mean(cat(4,Ym,Yp0/3),4)-Ydiv)>0.9 & Ym<1.1 & Yp0>2.1) & ...
        ... ((Ycd-Ydiv*5)>3 & (Ydiv<-0.01 & (Yg + max(0,0.05-Ycd/100))<0.1) & Yb & Ym>0.4 & Ybd<20 & Ycd>2.5) ) & ... % further WM
        ~Ybv & Yb & Ybd>1 & (Ycd>1.0 | (Yvt & Yp0>2.9) | (mean(cat(4,Ym,Yp0/3),4)-Ydiv)>0.9 & Ym<1.1 & Yp0>2.1) & (Yg+Ydiv<(Ybd/50) | (Ydiv-Ym)<-1.2); % Ybd/800 + Ycd/50
  Ywm = Ywm & ~(Ygmt<3 & Ycd<3 & (Ym-Ydiv)<0.90);   
  Ywms = smooth3(Ywm); Ywm(Ywms>0.75)=1; clear Ywms;
  %%
  Ygm = ~Yvt & Ybb & ~Ybv & ~Ywm & ~Ycm & Ycd>0.5 & (Ym-Ydiv-max(0,2-Ycd)/10)<0.9 & ... (Ym+Ydiv)>0.5 & ... ~Ysk & 
        (Yclsr{1}>4 | (Ym>0.7 & Yclsr{3}>64) | Ycd<(Ym+Ydiv)*3 ) & Ypp>0.2 & ...
        smooth3(Yg>(Ybd/800) & Yclsr{2}<240 )>0.6; % avoid GM next to hard boundies in the middle of the brain
  Ygx = Ybb & ~Ycm & Ym>1/3 & Ym<2.8/3 & Yp0<2.5 & Yg<0.4 & (Ym-Ydiv)>1/3 & (Ym-Ydiv)<1; Ygx(smooth3(Ygx)<0.5) = 0;
  Ygm = Ygm | Ygx; clear Ygx;
  Ygm = Ygm | (Ym>1.5/3 & Ym<2.8/3 & Yp0<2.5 & ~Ycm & Ybb);
  Ygm = Ygm | (Ygmt>eps & Ygmt<8 & Ycd<mean(Ygmt(Ygmt>0))*(1.5-Ym) & Ym<0.98 & Yl1 & ~Yvt & Yl1<3 & Yp0<2.5 & Yp0>1.5);
  Ygm = Ygm | (Yp0>2 & Yp0<2.5);
  %%
  %Ygm = Ygm & Ypp>0.5/max(1,Ygmt*2) & ~((Ypp-Ygmt/10)>0.2 & Ygmt>max(2,min(3,mean(Ygmt(Ygmt>0))*0.75))); % & ~Ywmh;
  Ygm  = Ygm & (Ygmt>2 | Ypp>0.2); 
  Ygms = smooth3(Ygm); Ygm(Ygms<0.25)=0; Ygm(Ygms>0.75 & Ypp>0.2)=1; clear Ygms;
   
  %% Ygw = Ygm & smooth3(Ywm)<0.1 & smooth3(Ycm)<0.4 & Ycd>0 & Ycd<2 & Ydiv<0.4 & Ydiv>-0.3 & Yg<0.1; %& (Ydiv>-0.4 | Ycd>1.5)
  if ~debug, clear Ybv  Ycp; end %Ycd

  if debug>1
    try %#ok<TRYNC> 
      [pth,nam] = spm_fileparts(res.image0(1).fname); tpmci=0;
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'prepeaks'));
      save(tmpmat);
    end
  end



  %% ------------------------------------------------------------------
  % SPM GM segmentation can be affected by inhomogeneities and some GM
  % is missclassified as CSF/GM (Ycls{5}). But for some regions we can 
  % trust these information more
  % ------------------------------------------------------------------
  Ybd  = cat_vbdist(single(~Yb),Yb,vx_vol);
  Ycbp = cat_vbdist(single(NS(Yl1,LAB.CB)),Yb,vx_vol);                    % next to the cerebellum
  Ycbn = cat_vbdist(single(~NS(Yl1,LAB.CB)),Yb,vx_vol);                   % not to deep in the cerebellum
  Ylhp = cat_vbdist(single(mod(Yl1,2)==1 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the left hemisphere 
  Yrhp = cat_vbdist(single(mod(Yl1,2)==0 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the righ hemishpere
  Ybv2 = Yclsr{5}>2 & Ym<0.7 & Ym>0.3 & Yb & (... 
         ((Ylhp+Ybd/2)<cleanupdist*6 & (Yrhp+Ybd/2)<cleanupdist*6) | ... % between the hemispheres next to skull                 
         ((Ycbp+Ybd/2)<cleanupdist*8 & (Ycbn+Ybd/2)<cleanupdist*8));     % between cerebrum and cerebellum next to hull
  Ybv2 = smooth3(Ybv2)>0.5;
  Ybvv = (Ym-max(0,6-abs(Ycbp-6))/50)<0.6 & Ym>0.4 & Yb & Ycbp<8 & Ycbp>1;

  % subcortical map refinements
  THth = 0.8 - LASstr*0.6; %0.5; % lower more thalamus
  YTH = NS(Yl1,LAB.TH) | (cat_vol_morph(NS(Yl1,LAB.TH),'d',3) & Ym>0.5 & Yclsr{1}>128);
  Ytd = cat_vbdist(single(Ym<0.45),YTH | NS(Yl1,LAB.BG),vx_vol); Ytd(Ytd>2^16)=0; % CSF distance in the TH
  Yxd = cat_vbdist(single(NS(Yl1,LAB.BG)),YTH,vx_vol); Yxd(Yxd>2^16)=0; % BG distance in the TH
  %Yyd = cat_vbdist(single(NS(Yl1,LAB.TH)),NS(Yl1,LAB.BG),vx_vol); Yyd(Yyd>2^16)=0; % TH distance in the BG
  Yss = NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH); 
  Yss = Yss | (cat_vol_morph(Yss,'d',vxv*2) & Ym>2.25/3 & Ym<2.75/3 & Ydiv>-0.01); % add ihger tissue around mask
  Yss = Yss | (cat_vol_morph(Yss,'d',vxv*3) &  NS(Yl1,LAB.VT) & Yp0>1.5 & Yp0<2.3); % add lower tissue around mask
  Yss = Yss & Yp0>1.5 & (Yp0<2.75 | (Ym<(2.5+LASstr*0.45)/3 & Ydiv>-0.05)); % by intensity
  Yss = Yss | ((Yxd./max(eps,Ytd+Yxd))>THth/2 & (Yp0<2.75 | (Ym<(2.75+LASstr*0.20)/3 & Ydiv>-0.05))); 
  % save TH by distances - for overcorrected images
  Yss = cat_vol_morph(Yss,'o');
  Ynw = (Yxd./max(eps,Ytd+Yxd))>THth/2 | (NS(Yl1,LAB.BG) & Ydiv>-0.01);
  if ~debug, clear Ytd Yxd ; end
  % increase CSF roi
  Yvt = cat_vol_morph( (NS(Yl1,LAB.VT) | cat_vol_morph(Ycm,'o',3) ) ...
    & Ycm & ~NS(Yl1,LAB.BG) & ~NS(Yl1,LAB.TH) & Ybd>30,'d',vxv*3) & ~Yss; % ventricle roi to avoid PVE GM between WM and CSF
  Ycx = (NS(Yl1,LAB.CB) & ((Ym-Ydiv)<0.55 | Yclsr{3}>128)) | (((Ym-Ydiv)<0.45 &  Yclsr{3}>8)| Yclsr{3}>240);
  % in the crebellum tissue can be differentated by div etc.
  Ycwm = NS(Yl1,LAB.CB) & (Ym-Ydiv*4)>5/6 & Ycd>3 & Yg>0.05;
  Yccm = NS(Yl1,LAB.CB) & Ydiv>0.02 & Ym<1/2 & Yg>0.05;
  Ybwm = (Ym-Ydiv*4)>0.9 & Ycd>3 & Yg>0.05; %Ydiv<-0.04 & Ym>0.75 & Ycd>3;
  Ybcm = Ydiv>0.04 & Ym<0.55 & Yg>0.05;
  % correction 1 of tissue maps
  Ywmtpm = (Ywtpm.*Ym.*(1-Yg-Ydiv).*cat_vol_morph(NS(Yl1,1).*Ybd/5,'e',1))>0.6; % no WM hyperintensities in GM!
  Ygm = Ygm | (Yss & ~Yvt & ~Ycx & ~Ybv2 & ~Ycwm & ~(Yccm | Ybcm));
  Ygm = Ygm & ~Ywmtpm & ~Ybvv; % no WMH area
  Ygm = Ygm & ~Yvt; % & ~Ywmh; 
  Ywm = (Ywm & ~Yss & ~Ybv2  & ~Ynw) | Ycwm | Ybwm; clear Ybwm; %& ~NS(Yl1,LAB.BG)
  Ywmtpm(smooth3(Ywmtpm & Ym<11/12)<0.5)=0;
  Ywm = Ywm & ~Ywmtpm & ~Ybvv & ~Yss; % no WM area
  Ycm = Ycm | ( (Ycx | Yccm | Ybcm) & Yg<0.2 & Ym>0 & Ydiv>-0.05 & Ym<0.3 & Yb ) | Ybvv;
  if ~debug, clear Ycwm Yccm Ycd; end
  % mapping of the brainstem to the WM (well there were some small GM
  % structures, but the should not effect the local segmentation to much.
  Ybs = cat_vol_morph(NS(Yl1,LAB.BS) & Ym<1.2 & Ym>0.9 & Yp0>2.5,'c',2*vxv) & Ym<1.2 & Ym>0.9 & Yp0>1.5;
  Ygm = (Ygm & ~Ybs & ~Ybv2 & ~Ywm) | Yss;
  Ywm = Ywm | (Ybs & Ym<1.1 & Ym>0.9 & Yp0>1.5) ; 
  if ~debug, clear Ycx; end

  %% Parahippocampal Gyrus for surface reconstruction
  % noch nicht n?tig
  %{
  Yphcg = Ywtpm>0.1 & Ydiv<0.02 & Ym>2.2/3 & Ym<3.5/3 &...
    cat_vol_morph(NS(Yl1,LAB.PH),'d',3/mean(vx_vol)) & ...
    ~cat_vol_morph(NS(Yl1,LAB.BS) | NS(Yl1,LAB.CB),'d',1/mean(vx_vol)); 
  Yphcg = cat_vol_morph(Yphcg & mod(Yl1,2)==0,'l') | ...
          cat_vol_morph(Yphcg & mod(Yl1,2)==1,'l');
  %}
  
  
  %% back to original resolution for full bias field estimation
  [Ycm,Ygm,Ywm]      = cat_vol_resize({Ycm,Ygm,Ywm},'dereduceBrain',BB); % ,Ygw
  [Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2] = cat_vol_resize({Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2},'dereduceBrain',BB);
  [Ym,Yp0,Yl1] = cat_vol_resize({Ym,Yp0,Yl1},'dereduceBrain',BB);
  Ybd = cat_vol_resize(Ybd,'dereduceBrain',BB);
  %Yphcg = cat_vol_resize(Yphcg,'dereduceBrain',BB);
  [Yg,Ydiv] = cat_vol_resize({Yg,Ydiv},'dereduceBrain',BB);
  [Ygmt,Ypp] = cat_vol_resize({Ygmt,Ypp},'dereduceBrain',BB);
  [Ywd] = cat_vol_resize(Ywd,'dereduceBrain',BB); % Ysk
  [Ycd] = cat_vol_resize(Ycd,'dereduceBrain',BB); % Ysk
  clear Yclso Ybv;
  
  if debug>1
    try %#ok<TRYNC> % windows requires this... i don't know why ... maybe the file size
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'prepeaks'));
      save(tmpmat);
    end
  end
  
  
%% --------------------------------------------------------------------- 
%  Now, we can estimate the local peaks 
%  ---------------------------------------------------------------------
  % Estimation of the local WM threshold with "corrected" GM voxels to
  % avoid overfitting (see BWP cerebellum). 
  % CSF is problematic in high contrast or skull-stripped image should 
  % not be used here, or in GM peak estimation
  mres = 1.1; 
  stime = cat_io_cmd('  Estimate local tissue thresholds','g5','',verb,stime); dispc=dispc+1;
  Ysrcm = cat_vol_median3(Ysrc.*Ywm,Ywm,Ywm); 
  rf    = [10^5 10^4];
  T3th3 = max(1,min(10^6,rf(2) / (round(T3th(3)*rf(1))/rf(1))));
  Ysrcm = round(Ysrcm*T3th3)/T3th3;
  Ygw2 = Ycls{1}>128 & Ym>2/3-0.04 & Ym<2/3+0.04 & Ygm .*Ydiv>0.01;
  Ygw2 = Ygw2 | (Ycls{1}>128 & Yg<0.05 & abs(Ydiv)<0.05 & ~Ywm & Ym<3/4); % large stable GM areas - like the BWP cerebellum
  Ygw3 = Ycls{3}>128 & Yg<0.05 & ~Ywm & ~Ygm & Ywd<3; 
  Ygw3(smooth3(Ygw3)<0.5)=0;
  [Yi,resT2] = cat_vol_resize(Ysrcm,'reduceV',vx_vol,mres,32,'max'); % maximum reduction for the WM
  %%
  if cat_stat_nanmean(Ym(Ygw3))>0.1, % not in images with to low CSF intensity (error in skull-stripped)
    Ygi = cat_vol_resize(Ysrc.*Ygw2*T3th(3)/mean(Ysrc(Ygw2(:))) + ...
      Ysrc.*Ygw3*T3th(3)/mean(Ysrc(Ygw3(:))),'reduceV',vx_vol,mres,32,'meanm'); clear Ygw2; % mean for other tissues
  else
    % mean for other tissues
    Ygi = cat_vol_resize(Ysrc.*Ygw2*T3th(3)/mean(Ysrc(Ygw2(:))),'reduceV',vx_vol,mres,32,'meanm'); clear Ygw2; 
  end
  for xi=1:2*LASi, Ygi = cat_vol_localstat(Ygi,Ygi>0,2,1); end; Ygi(smooth3(Ygi>0)<0.3)=0;
  Yi = cat_vol_localstat(Yi,Yi>0,1,3); % one maximum for stabilization of small WM structures
  Yi(Yi==0 & Ygi>0)=Ygi(Yi==0 & Ygi>0);
  for xi=1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end % no maximum here!
  Yi = cat_vol_approx1585(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,LASfs); 
  Ylab{2} = max(eps,cat_vol_resize(Yi,'dereduceV',resT2)); 
 % Ylab{2} = Ylab{2} .* mean( [median(Ysrc(Ysw(:))./Ylab{2}(Ysw(:))),1] ); 
  if debug==0; clear Ysw; end

  %% update GM tissue map
  %Ybb = cat_vol_morph((Yb & Ysrc./Ylab{2}<(T3th(1) + 0.25*diff(T3th(1:2))) & Ydiv<0.05) | Yp0>1.5,'lo',1);
  Ygm(Ysrc./Ylab{2}>(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3))=0; % correct GM mean(T3th([2:3,3]))/T3th(3) 
  Ygm(Ysrc./Ylab{2}<(T3th(2) + 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ysrc./Ylab{2}<(T3th(2) - 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ydiv<0.3 & Ydiv>-0.3 & Ybb & ~Ywm & ~Yvt & ~Ybv2 & Ycls{1}>48)=1;
  Ygm = Ygm | (Ygmt>eps & Ygmt<8 & Ycd<mean(Ygmt(Ygmt>0))*(1.5-Ym) & NS(Yl1,1) & ~Ybs);
  Ygm = Ygm & Ypp>0.5/max(1,Ygmt) & ~((Ypp-Ygmt/10)>0.2 & Ygmt>max(2,min(3,mean(Ygmt(Ygmt>0))*0.75)));
  Ywmd2 = cat_vbdist(single(Ywm),Yb);
  Ygx = Ywmd2-Ym+Ydiv>0.5 & Ym+0.5-Ydiv-Yg-Ywmd2/10>1/3 & ~Ybv2 & ... low intensity tissue
    ~(Ym-min(0.2,Yg+Ywmd2/10-Ydiv)<1/4) & Yg<Ylab{2}/T3th(3)*0.3 & Ysrc<Ylab{2}*0.9; % no real csf
  Ygx(smooth3(Ygx)<0.5)=0; 
  Ygm = ~Ywm & (Ygm | Ygx); % correct gm (intensity based)
  Ygx = (single(Ycls{1})/255 - abs(Ydiv) + min(0,Ydiv) - Yg)>0.5 & ~Ywm;
  Ygx(smooth3(Ygx)<0.5)=0; 
  Ygm = Ygm | Ygx; % correct gm (spm based)
  %%
  Ycm = ~Ygm & ~Ywm & ~Ybv2 & Yg<0.6 & (Ycm | (Yb & (Ysrc./Ylab{2})<((T3th(1)*0.5 + 0.5*T3th(2))/T3th(3)))); 
  %
  Ycp = (Ycls{2}<128 & Ydiv>0 & Yp0<2.1 & Ysrc./Ylab{2}<mean(T3th(1)/T3th(2))) | Ycm | ...      % typcial CSF
        (Ycls{5}>32 & Ycls{2}<32 & Ysrc./Ylab{2}<T3th(2)/T3th(3) & Ydiv>0) | ...                % venes
        ((Ym-Ydiv<0.4) & Ycls{3}>4 & Ycls{3}>16 & Ysrc./Ylab{2}<mean(T3th(2)/T3th(3))) | ...    % sulcal CSF
        (single(Ycls{6})+single(Ycls{5})+single(Ycls{4}))>192 | ...                             % save non-csf 
        Ysrc./Ylab{2}<T3th(1)/T3th(3);                                                          % but do not trust the brain mask!
  Ycp(smooth3(Ycp)>0.4)=1;                                                                      % remove some meninges
  Ycp(Ypp<(0.2+0.1*max(0,5-Ygmt)) & ~Yss)=1; 
  Ycd = cat_vbdist(single(Ycp),~Ycp,vx_vol);  
  %%
  Ygm = Ygm & ~Ycm & ~Ywm & Ywd<5; % & ~Ywmh; %  & ~Ybvv  & ~Ysk
  Ygm = Ygm | (NS(Yl1,1) & Ybd<20 & (Ycd-Ydiv)<2 & Ycls{1}>0 & ~Ycm & Ybb & Ym>0.6 & Yg<max(0.5,1-Ybd/30)); 
  Ygm = Ygm & (Yg<0.1 | Ysrc./Ylab{2}<(T3th(2)*1/3+2/3*T3th(3))/T3th(3)); % outer high intensity GM
  Ygm = (Ygm | Yss) & ~Ycm & cat_vol_morph(~Ybs | ~Yvt,'e');
  %%
  %Ygm = Ygm & ~Yvt & Ypp>0.5/max(1,Ygmt*2) & ~((Ypp-Ygmt/10)>0.2 & Ygmt>max(2,min(3,mean(Ygmt(Ygmt>0))*0.75)));
  Ybb = cat_vol_morph(smooth3(Ygm | Ywm | Yp0>1.5 | (Ym>1.2/3 & Ym<3.1/3 & Yb))>0.6,'lo',min(1,vxv)); % clear Yp0 Yvt
  Ygm(~Ybb)=0; Ygm(smooth3(Ygm)<0.3)=0;
  Ygm(smooth3(Ygm)>0.4 & Ysrc./Ylab{2}>mean(T3th(1)/T3th(2)) & Ysrc./Ylab{2}<(T3th(2)*0.1+0.9*T3th(3)) & Ypp>0.2)=1;
  Ygm = (Ygm & Ypp>0.2 & Ywd<3) | (Ygm & Ywd>=2);
  %%
  if debug==0; clear Ybb Ybd Yvt Ybvv Ycp Ycd Yl1 Yss; end %Ydiv Yg
  %Ygm = Ygm & smooth3(smooth3(Ypp)>0.1)>0.7 & Ym<0.90; 
  
  %% GM
%  Yi = (Ysrc./Ylab{2} .* Ygm , Ysrc./Ylab{2}.*Ywm.*Ybs.*(T3th(2)+0.8*diff(T3th(2:3)))/T3th(3));
%  Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); 
  Yi = Ysrc./Ylab{2} .* Ygm; 
  Yi = round(Yi*rf(2))/rf(2);
 % Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); % ????
  Yi(Ybs)  = Ysrc(Ybs)./Ylab{2}(Ybs)   .* T3th(2)/T3th(3); 
  Yi = cat_vol_median3(Yi,Yi>0.5,Yi>0.5); 
%  Ycmx = smooth3(Ycm & Ysrc<(T3th(1)*0.8+T3th(2)*0.2))>0.9; Tcmx = mean(Ysrc(Ycmx(:))./Ylab{2}(Ycmx(:)))*T3th(3);
%  Yi(Ycmx) = Ysrc(Ycmx)./Ylab{2}(Ycmx)  .* T3th(2)/Tcmx; 
  %Yii =  Ysrc./Ylab{2} .* Ycm * T3th(2) / cat_stat_nanmedian(Ysrc(Ycm(:))); 
  [Yi,Yii,resT2] = cat_vol_resize({Yi,Ylab{2}/T3th(3)},'reduceV',vx_vol,1,32,'meanm');
  for xi=1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,3,1); end
  Yi = cat_vol_approx1585(Yi,'nh',resT2.vx_volr,2); 
  Yi = min(Yi,Yii*(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3));
  Yi = cat_vol_smooth3X(Yi,LASfs); 
  Ylab{1} = cat_vol_resize(Yi,'dereduceV',resT2).*Ylab{2};   
  %Ylab{1}(Ygm) = Ysrc(Ygm); Ylab{1} = cat_vol_smooth3X(Ylab{1},LASfs); % can lead to overfitting
  Ycm = (single(Ycls{3})/255 - Yg*4 + abs(Ydiv)*2)>0.5 &  Ysrc<(Ylab{1}*mean(T3th([1,1:2]))/T3th(2));
  %Ycm & Ysrc<(Ylab{1}*mean(T3th(1:2))/T3th(2)) & Yg<0.1;
  Ycm(smooth3(Ycm)<0.5)=0;
  Ycm(Yb & cat_vol_morph(Ysrc<mean(T3th(1:2)),'o'))=1;
  
  %% CSF & BG 
  Ynb = cat_vol_morph(smooth3(Ycls{6})>128 | (~Yb & Yg<=0.001),'e',4*vxv); 
  [Yx,Yc,resT2] = cat_vol_resize({round(Ysrc./Ylab{2} .* Ynb .* ((Ysrc./Ylab{2})<0.1) * rf(2))/rf(2),...
    round(Ysrc./Ylab{2} .* (smooth3(Ycm | Ysc)>0.5) * rf(2))/rf(2)},'reduceV',vx_vol,8,16,'min');% only pure CSF !!!
  Yx(Yc>0)=0; Yc(Yx>0)=0;
  meanYx = min(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  meanYc = max(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  stdYbc = mean([std(Yc(Yc(:)>0)),std(Yx(Yx(:)>0))]);
  %Yx = min(max(meanYx/2,Yx),min(meanYx*4,meanYc/2));
  %Yc = min(max(meanYc/2,Yx),meanYc/2);
  Yxa = cat_vol_approx1585(Yx ,'nh',resT2.vx_volr,16); %+(Yb>0).*stdYbc  + Yc.*meanYb/max(eps,meanYc)
  Yca = cat_vol_approx1585(Yc + min(max( meanYx + stdYbc , meanYc - stdYbc ),...
    Yx.*meanYc/max(eps,meanYx)),'nh',resT2.vx_volr,16); % + Yb.*meanYc/max(eps,meanYb)
  Yca = Yca*0.7 + 0.3*max(mean(Yca(:)),T3th(1)/T3th(3));
  %%
  Yxa = cat_vol_smooth3X(Yxa,LASfs*2); 
  Yca = cat_vol_smooth3X(Yca,LASfs*2); 
  Ylab{3} = cat_vol_smooth3X(cat_vol_resize(Yca,'dereduceV',resT2).*Ylab{2},LASfs*2);  
  Ylab{6} = cat_vol_smooth3X(cat_vol_resize(Yxa,'dereduceV',resT2).*Ylab{2},LASfs*2);
  clear Yxa Yca Yx Yc Y %Ydiv
  
  %% local intensity modification of the original image
  % --------------------------------------------------------------------
  Yml = zeros(size(Ysrc));  
  Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc-Ylab{2}) ./ max(eps,Ylab{2}-Ylab{1})) );
  Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc-Ylab{1}) ./ max(eps,Ylab{2}-Ylab{1})) );
  Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc-Ylab{3}) ./ max(eps,Ylab{1}-Ylab{3})) );
  Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc-Ylab{6}) ./ max(eps,Ylab{3}-Ylab{6})) );
  Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
  
  
  %% global
  Ymg = max(eps,Ysrc./Ylab{2});
  Ymg = cat_main_gintnorm(Ymg*Tth.T3th(5),Tth); 
  
  %%
  if debug>1
    try %#ok<TRYNC> % windows requires this... i don't know why
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'postpeaks'));
      save(tmpmat);
    end
  end
  
  %% fill up CSF in the case of a skull stripped image 
  if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.3))
    YM   = cat_vol_morph(Yb,'d'); 
    Ymls = smooth3(max(Yml,YM*0.5));
    Yml(YM & Yml<0.5)=Ymls(YM & Yml<0.5); 
    clear Ymls YM
  end
  
  
  %% class correction and second logical class map Ycls2
  Ynwm = Ywm & ~Ygm & Yml/3>0.95 & Yml/3<1.3;
  Ynwm = Ynwm | (smooth3(Ywm)>0.6 & Yml/3>5/6); Ynwm(smooth3(Ynwm)<0.5)=0;
  Yngm = Ygm & ~Ywm & Yml/3<0.95; Yngm(smooth3(Yngm)<0.5)=0;
  Yncm = ~Ygm & ~Ywm & ((Yml/3)>1/6 | Ycls{3}>128) & (Yml/3)<0.5 & Yb;
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) + (Ynwm & ~Yngm & Yp0>=1.5)*256 - (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Ynwm & ~Yngm & Yp0>=1.5)*256 + (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) - ((Ynwm | Yngm) & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) + (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls2 = {Yngm,Ynwm,Yncm};
  clear Yngm Ynwm Yncm;
  
  %%
  Yml = Yml/3;
  cat_io_cmd('','','',verb,stime);
 % if debug
 %   cat_io_cmd(' ','','',verb,stime); 
 % else
 %   cat_io_cmd(' ','','',verb,stime);   
 %   cat_io_cmd('cleanup',dispc,'',verb);
 % end

end

