function [Yml,Ymg,Ycls,Ycls2,T3th] = ...
  cat_main_LAS2(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,extopts,Tth) 
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
%
% Local Adaptive Segmentation (LAS):
%
% This version of the local adaptive intensity correction includes a  
% bias correction that is based on a maximum filter for the WM, a mean
% filter of GM and a mean filter of the head tissue (if available).
% For each tissue a refined logical map is generated to filter the local
% tissue intensity and approximate the local tissue intensity in the whole 
% volume. Based on these values an intensity transformation is used. 
% Compared to the global correction this has to be done for each voxel. 
% To save time only a rough linear transformation is used. Finally, a 
% second NLM-filter is used.
% It is important to avoid high intensity blood vessels in the process  
% because they will push down local WM and GM intensity. 
% 
% ______________________________________________________________________
%
%   [Yml,Ymg,Ycls,Ycls2,T3th] = ...
%     cat_main_LAS(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,extopts,Tth) 
%
%   Yml     .. local  intensity normalized image
%   Ymg     .. global intensity normalized image
%   Ycls    .. corrected SPM tissue class map
%   Ycls2   .. ?
%   T3th    .. tissue thresholds of CSF, GM, and WM in Ysrc
% 
%   Ysrc    .. (bias corrected) T1 image
%   Ycls    .. SPM tissue class map
%   Ym      .. intensity corrected T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%   Yb0     .. brain mask
%   Yy      .. deformation map
%   T3th    .. intensity thresholds from global intensity normalization
%   res     .. SPM segmentation structure
%   vx_vol  .. voxel dimensions
%   extopts ..
%   Tth     ..
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_LAS2.m 1523 2019-11-21 23:12:24Z gaser $



% ______________________________________________________________________
% 
% internal maps:
%
%   Yg   .. gradient map   - edges between tissues
%   Ydiv .. divergence map - sulci/gyri pattern, and blood vessels
%   Yp0  .. label map      - tissue classes (BG=0,CSF=1,GM=2,WM=3) 
%
%   Ysw  .. save WM tissue map
%   Ybv  .. blood vessel map
%   Ycp  .. CSF / background area for distances estimation
%   Ycd  .. CSF / background distance map
%   Ywd  .. WM distance map
%
%   Ycm  .. CSF
%   Ygm  .. GM
%   Ywm  .. WM 
%   Yvt  .. WM next to the ventricle map 
% ______________________________________________________________________
%
% Development / TODO: 
% * replace old calls of cat_vol_morph by distance based operations with 
%   resolution parameter 
% ______________________________________________________________________



  % The reduction to 1 mm is not really good for the Ycls map. If this 
  % is required to support faster and memory saving preprocessing also  
  % for ultra-high resolution data further test are necessary!
  %def.uhrlim     = 0.7; 
  %extopts        = cat_io_checkinopt(extopts,def); 
  extopts.uhrlim = 0.2; % no reduction for >0.4 mm 
  
  
  % set this variable to 1 for simpler debugging without reduceBrain
  % function (that normally save half of processing time)
  verb        = extopts.verb - 1;
  vxv         = 1 / mean(vx_vol);                % normalization of voxel size (mostly for old calls of cat_vol_morph) 
  dsize       = size(Ysrc);
  NS          = @(Ys,s) Ys==s | Ys==s+1;         % function to ignore brain hemisphere coding
  LASstr      = max(eps,min(1,extopts.LASstr));  % LAS strength (for GM/WM threshold) - manual correction based on R1109 (2017/02)
  LAB         = extopts.LAB;                     % atlas labels
  LABl1       = 1;                               % use atlas map
  cleanupstr  = min(1,max(0,extopts.gcutstr));   % required to avoid critical regions (only used in case of atlas maps)
  cleanupdist = min(3,max(1,1 + 2*cleanupstr));
    
  
  % set debug = 1 and do not clear temporary variables if there is a breakpoint in this file 
  dbs = dbstatus; debug = 0; 
  for dbsi = 1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_main_LAS'); debug = 1; break; end; end
  
  
  
  
%% ---------------------------------------------------------------------
%  First, we have to optimize the segments using further information that 
%  SPM do not use, such as the gradient, divergence and distance maps. 
%  (1) The gradient map (average of the first derivate of the T1 map) is  
%      an edge map and independent of the image intensity. It helps to   
%      avoid PVE regions and meninges. 
%  (2) The divergence (second derivate of the T1 map) help to identify 
%      sulcal and gyral pattern and therefore to find WM and CSF regions   
%      for further corrections and to avoid meninges and blood vessels. 
%  (3) Distance maps can help to describe the cortex by its expected 
%      thickness or to separate structures close to the skull or deep 
%      in the brain. 
%
%  Furthermore, special assumption can be used. 
%  The first is that the maximum property of the WM in T1 data that 
%  allows using of a maxim filter for the GM/WM region. 
%  The second is the relative stable estimation of CSF/BG that allows to 
%  estimate a distance map. Because, most regions have a thin layer of 
%  GM around the WM, we can avoid overestimation of the WM by the other 
%  maps (especially the divergence). 
%  ---------------------------------------------------------------------

  fprintf('\n');
  stime = cat_io_cmd('  Prepare maps','g5','',verb); 
  
  
  %  general resolution limitation
  %  -------------------------------------------------------------------
  %  Using a lot of maps make this function memory intensive that is
  %  critical for high resolution data. Furthermore, it can be expected
  %  that the full resolution is not required here. However, reducing 
  %  the tissue class map can lead to other problems. Hence, the uhrlim
  %  parameter is maybe inactive (see above).
  %  -------------------------------------------------------------------
  if any( vx_vol < extopts.uhrlim/2 )
    % store old data that is needed in full resolution
    Yclso2 = Ycls; 
    Ysrco2 = Ysrc;
    
    % reduce all input maps (Ysrc, Ycls, Yy, Ym, Yb0)
    [Ysrc,resT0] = cat_vol_resize( Ysrc , 'reduceV' , vx_vol , extopts.uhrlim , 64 ); 
    for ci = 1:numel(Ycls)
      Ycls{ci} = cat_vol_resize( Ycls{ci}, 'reduceV' , vx_vol , extopts.uhrlim , 64 ); 
    end
    Yy2 = ones([size(Ysrc),4],'single'); 
    for ci = 1:numel(Yy)
      Yy2(:,:,:,ci) = cat_vol_resize( Yy(:,:,:,ci) , 'reduceV' , vx_vol , extopts.uhrlim , 64 );
    end
    Yy  = Yy2; clear Yy2; 
    Ym  = cat_vol_resize( Ym          , 'reduceV' , vx_vol , extopts.uhrlim , 64 ); 
    Yb0 = cat_vol_resize( single(Yb0) , 'reduceV' , vx_vol , extopts.uhrlim , 64 )>0.5; 
    
    vx_vol  = resT0.vx_volr; 
    vxv     = 1/ mean(resT0.vx_volr);
    dsize   = size(Ysrc);
  end
  
  
  % brain bounding box to save processing time 
  Yclso   = Ycls; 
  [Ym,BB] = cat_vol_resize( Ym  , 'reduceBrain' , vx_vol , round(10/mean(vx_vol)) , Yb0 );
  Yb      = cat_vol_resize( Yb0 , 'reduceBrain' , vx_vol , round(10/mean(vx_vol)) , Yb0 ); clear Yb0; 
  for i = 1:6, Ycls{i} = cat_vol_resize(Ycls{i} , 'reduceBrain' , vx_vol , BB.BB); end
  
  
  % helping maps (Yg = mean gradient = edge) and divergence (CSF or WM skeleton)
  Yg    = cat_vol_grad( Ym , vx_vol );
  Ydiv  = cat_vol_div( max(0.33,Ym) , vx_vol );
  Yp0   = single( Ycls{1} )/255*2 + single( Ycls{2} )/255*3 + single( Ycls{3} )/255;  % create label map
  Yb    = smooth3( Yb | cat_vol_morph(Yb,'d',2*vxv) & Ym<0.8 & Yg<0.3 & Ym>0 )>0.5;   % increase brain mask, for missing GM 
  
  
  
  
  %% adding of atlas and template information (for subcortical structures)
  %  -------------------------------------------------------------------
%%%  20180801 - Why not use the partitioning before LAS to benefit by the 
  %             better ventricle, WMH and cerebellum maps?
  %           - Because the partitioning benefit in subcortical regions 
  %             from LAS. >> call it twice?
  %  -------------------------------------------------------------------
  if LABl1 
    stime = cat_io_cmd('  Prepare partitions','g5','',verb,stime); 

    % map atlas to RAW space
    for i=1:5
      try
        Vl1A = spm_vol(extopts.cat12atlas{1});
        break
      catch 
        % avoid read error in parallel processing
        pause(rand(1))
      end
    end
    Yl1  = cat_vol_ctype(round(spm_sample_vol(Vl1A,...
      double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
    Yl1  = reshape(Yl1,dsize);
    
    
    % load WM of the Dartel/Shooting Template for WMHs - use uint8 to save memory 
    % Ywtpm .. WM template map
    Vtemplate = spm_vol(extopts.templates{end}); 
    Ywtpm  = cat_vol_ctype(spm_sample_vol(Vtemplate(2),...
      double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
    if debug==0, clear Yy; end
    Ywtpm  = single(reshape(Ywtpm,dsize)); 
    spm_smooth(Ywtpm,Ywtpm,2*vxv); 
    Ywtpm = cat_vol_ctype(Ywtpm); 
    
    
    % apply boundary box for brain mask
    Yl1    = cat_vol_resize(Yl1  ,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
    Ywtpm  = cat_vol_resize(Ywtpm,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
    % The correction of image resolution is not required because of the adaption of the Yy?
    
    
    % do not reduce LASstr 
    LASmod = min(2,max(1,mean((Ym( NS(Yl1,LAB.BG) & Yg<0.1 & Ydiv>-0.05  & Ycls{1}>4)) - 2/3) * 8)); 
  else
    LASmod = 1; %#ok<UNRCH>
  end
  
  
  % adaption of the LASstr depending on average basal values 
  LASstr  = min(1,max(0.01,LASstr * LASmod));   % adaption by local BG variation
  LASfs   = 1 / max(0.01,LASstr);               % smoothing filter strength 
  LASi    = min(8,round(LASfs));                % smoothing iteration (limited)
   
  
  
  
  %% helping segments
  %  -------------------------------------------------------------------
  %  We will now try to identify the CSF and WM regions those finally  
  %  describe the GM area.
  %  * brain mask (Ybb), blood vessel maps (Ybv), 
  %  * distance maps from CSF (Ycd), WM (Ywd), and brain mask (Ybd)
  %  * tissue maps for CSF (Ycp >> Ycm), GM (Ygm), and WM (Ysw >> Ywm)
  %  * ventricle map (Yvt)
  %  -------------------------------------------------------------------
  stime = cat_io_cmd(sprintf('  Prepare segments (LASmod = %0.2f)',LASmod),'g5','',verb,stime); 
  
  % Ybb .. brain mask 
  %        Don't trust SPM too much by using Yp0 because it may miss some GM areas!
  %        20180801 - The brain mask should be correct now. 
  Ybb = cat_vol_morph((Yb & Ym>1.5/3 & Ydiv<0.05) | Yp0>1.5,'lo',vxv);
  if debug==0, clear Yb; end
  
  
  % Ysw .. save WM 
  % Ybv .. possible blood vessels
  Ysw1 = cat_vol_morph(Ycls{2}>128 & (min(1,Ym)-Ydiv)<1.5,'lc',vxv*2) & (Ym - Ydiv)>5/6; 
  Ybv  = ((min(1,Ym) - Ydiv + Yg)>2.0 | (Ycls{5}>16 & Ym<0.6 & Ycls{1}<192)) & ...
         ~cat_vol_morph(Ysw1,'d',1) & Ym>0.2;    
  
       
  % Ycp = mask for CSF/BG distance initialization 
  Ycp1 = Ycls{3}>240 & Ydiv>0 & Yp0<1.1 & Ym<0.5 ; % typical CSF
  Ycp2 = Ycls{5}>8 & Ycls{2}<32 & Ym<0.6 & Ydiv>0; % venes
  Ycp3 = (Ym-Ydiv/4<0.4) & Ycls{3}>4 & Ycls{3}>16; % sulcal CSF
  Ycp4 = (single(Ycls{6}) + single(Ycls{4}))>192;  % non-CSF .. class 5 with error for negative t1 values
  Ycp5 = ~cat_vol_morph(Ybb,'lc',5);               % add background
  Ycp6 = Ym<0.3;                                   % but do not trust the brain mask!
  Ycp  = Ycp1 | Ycp2 | Ycp3 | Ycp4 | Ycp5 | Ycp6;  % combine maps 
  Ycp(smooth3(Ycp)>0.4) = 1;                       % remove some meninges
  if debug==0; clear Ycp1 Ycp2 Ycp3 Ycp4 Ycp5 Ycp6; end
    
  
  % Ywd .. WM distance map    
  % Ysk .. divergence skeleton to improve CSF map 
  %        deactivated due to problems in non-cortical structures
  Ywd = cat_vbdist(single(Yp0>2.5),Yp0>0.5,vx_vol);                     % WM distance for skeleton
  %Ysk = cat_vol_div(min(5,Ywd),2);                                     % divergence skeleton
  %Ysk = (Ym + min(0,Ysk))<0.2;                                         % binary divergence skeleton
  %Ycp = Ycp | Ysk; if debug==0, clear Ysk; end                           % combine skeleton with CSF map  
  
  % Ycd .. CSF/BG distance map
  Ycd = cat_vbdist(single(Ycp),~Ycp,vx_vol);                            % real CSF distance 
  Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1) =  ...          % correction for sulci 
    min(Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1),1.5);     % ... maybe another distance estimation?
  
  
  % we need to remove strong edge regions, because here is no GM layer between CSF and WM ???
  % Yb  .. brain mask 
  % Ybd .. skull distance
%%% 20180801 - What is the difference of Yb2 to Ybb and the (previously replaced) Yb? 
  Yb2 = cat_vol_morph(~Ycp | (Ycls{3}>128),'lc',1);                     
  Ybd = cat_vbdist(single(~Yb2),Yb2,vx_vol); 
  
  
  % Yvt .. Ventricle map WITHOUT atlas data as large deep CSF areas
  %        need the ventricle to avoid PVE GM between WM and CSF.
  %        There is an update for atlas data.
  Yvt = (Yg + abs(Ydiv))>0.4 & smooth3(single(Ycls{1})/255)<0.5 & Ybd>20 & ...
    cat_vol_morph(Ycls{3}>8,'d',vxv) & cat_vol_morph(Ycls{2}>8,'d',vxv); 
  Yvt = smooth3(Yvt)>0.7;
  Yvt = smooth3(Yvt)>0.2;
  
  
  
  
  %% final tissue maps:  Ycm = CSF, Ygm = GM, Ywm = WM 
  %  -------------------------------------------------------------------
  
  % CSF:
  Ycm = Ycp & Yb2 & Ycls{3}>192 & ~Ybv & Ym<0.45 & Yg<0.25 & Ym>0 ; 
  %Ycm = Ycm | (Yb2 & (Ym - max(0,Ydiv))<0.5); % adding of divergence information was not save
  if debug==0, clear Ycp; end 

  
  %% WM:
  %  Different definitions of the possible WM (Ysw,Ysw2,Ysw3,Ysw4) 
  %  and a variable (Ygw) to avoid non WM areas. 
  %    Ysw1 .. defined above 
  %    Ysw2 .. addition WM map that further include CSF distance (Ycd) and
  %            divergence information (Ydiv) and use a flexible intensity 
  %    Ysw3 .. similar to Ysw2 with a skull-near and CSF distance criteria
  %            to reconstruct WM gyri
  %    Ysw4 .. was remove long ago
  %    Ygw  .. map to of regions that we want to avoid in all possible WM
  %            definitions Ysw*
  Ysw2 = Yb2 & (Ycd - Ydiv)>2        & Ydiv<0     & Ym>(0.9 + LASstr * 0.05); % general WM
  Ysw3 = Yb2 & (Ycd - Ydiv .* Ycd)>4 & Ydiv<-0.01 & Ym>0.5  & Ybd<20 & Ycd>2; % addition WM gyri close to the cortex 
  %Ysw4 = ( (Ycd - Ydiv*5)>3 & Yb2 & Ym>0.4 & Ybd<20 & Ycd>2.5) & ...         % further WM (rmoved long ago)
  %         (Ydiv<-0.01 & (Yg + max(0,0.05 - Ycd/100))<0.1 );                                  
  
  Ygw  = ~Ybv & Yb2 & Ybd>1 & (Ycd>1.0 | (Yvt & Yp0>2.9)) & ...
         (Yg + Ydiv<(Ybd/50) | (Ydiv - Ym)<-1); 
  
  Ywm  = ( Ycls{2}>252 | Ysw1 | Ysw2 | Ysw3 ) & Ygw; 
  if debug==0, clear Ysw Ysw2 Ysw3 Ygw; end 

  
  %% GM:
  %  Different criteria (Ygm1,Ygm2,Ygm3) to define the possible GM area 
  %  whereas Ygm4 describe general GM limitation.
  %  Ygm1 .. upper intensity limitation
  %  Ygm2 .. lower intensity limitation
  %  Ygm3 .. avoid GM next to hard boundaries in the middle of the brain
  Ygm1 = ( Ym - Ydiv - max(0,2 - Ycd)/10 )<0.9;                     % upper limit
  Ygm2 = Ycls{1}>4 | ( Ym>0.7 & Ycls{3}>64 ) | Ycd<(Ym + Ydiv)*3;   % lower limit
  Ygm3 = smooth3(Yg>(Ybd/800) & Ycls{2}<240 )>0.6;                  % avoid GM-WM PVE 
  Ygm  = Ybb & ~Yvt & ~Ybv & ~Ywm & ~Ycm & Ycd>0.5 & Ygm1 & Ygm2 & Ygm3;  
  if debug==0, clear Ybv Yvt Ygm1 Ygm2 Ygm3; end 
  
  %  Ygm4 .. general GM limitations
  Ygm4 = Ybb & ~Ycm & ~Ywm & Ym>1/3 & Ym<2.8/3 & ...
    Yg<0.4 & (Ym - Ydiv)>1/3 & (Ym - Ydiv)<1; 
  Ygm4( smooth3(Ygm4)<0.5 ) = 0; % remove small dots
  
  % combine possible area with limitation map
  Ygm  = Ygm | Ygm4;
  if debug==0, clear Ygm4; end
  Ygm  = Ygm | (Ym>1.5/3 & Ym<2.8/3 & ~Ywm & ~Ycm & Ybb);
  Ygm(smooth3(Ygm)<0.25) = 0;
  if debug==0, clear Ybv Ycp; end 


  
  
  
  %% further maps that require the atlas 
  if LABl1   
    % SPM GM segmentation can be affected by inhomogeneities and some GM
    % areas are miss-classified as CSF/GM (Ycls{5}) when the affine 
    % registration was not optimal. But for some regions we can trust
    % these information.
    
    
    
    %% regions for cleanup of meninges and blood vessels 
    %  Ycbp .. close to the cerebellum
    %  Ycbn .. not to deep in the cerebellum
    %  Ylhp / Yrhp .. GM next to the left/right hemisphere 
    %  Ybv2 .. close to the skull and between the hemispheres or between cerebrum and cerebellum 
    %  Ybvv .. sinus
    Ycbp = cat_vbdist(single( NS(Yl1,LAB.CB)),Yb2,vx_vol);                 % close to cerebellum
    Ycbn = cat_vbdist(single(~NS(Yl1,LAB.CB)),Yb2,vx_vol);                 % close to cerebellum surface
    Ylhp = cat_vbdist(single(mod(Yl1,2)==1 & Yb2 & Yl1>0),Yb2,vx_vol);     % GM next to the left  hemisphere 
    Yrhp = cat_vbdist(single(mod(Yl1,2)==0 & Yb2 & Yl1>0),Yb2,vx_vol);     % GM next to the right hemisphere
    Ybv2 = Ycls{5}>2 & Ym<0.7 & Ym>0.3 & Yb2 & (... 
           ((Ylhp+Ybd/2)<cleanupdist*6 & (Yrhp+Ybd/2)<cleanupdist*6) | ... % between the hemispheres next to skull                 
           ((Ycbp+Ybd/2)<cleanupdist*8 & (Ycbn+Ybd/2)<cleanupdist*8));     % between cerebrum and cerebellum next to hull
    if debug==0, clear Ycbn Ylhp Yrhp; end 
    Ybv2 = smooth3(Ybv2)>0.5;
    Ybvv = (Ym - max(0,6 - abs(Ycbp - 6))/50)<0.6 & Ym>0.4 & Yb2 & Ycbp<8 & Ycbp>1; % sinus
    if debug==0, clear Ycbp; end 
    
    
    
    %% refinements of subcortical regions
    % Thalamus (TH): 
    % YTH .. enlarged atlas ROI
    % Ytd .. CSF distance in the thalamus YTH
    % Yxd .. distance to the basal-ganglia (BGL) in the thalamus YTH
    THth = 0.8 - LASstr*0.6; %0.5; % lower more thalamus
    YTH  = NS(Yl1,LAB.TH) | (cat_vol_morph(NS(Yl1,LAB.TH),'d',3) & Ym>0.5 & Ycls{1}>128);
    Ytd  = cat_vbdist(single(Ym<0.45),YTH | NS(Yl1,LAB.BG),vx_vol); Ytd(Ytd>2^16)=0; % CSF distance in the TH
    Yxd  = cat_vbdist(single(NS(Yl1,LAB.BG)),YTH,vx_vol); Yxd(Yxd>2^16)=0;           % BGL distance in the TH
    if debug==0, clear YTH; end
    
    
    % Yss .. (enlarged) subcortical structures
    % Ynw .. no WM ??
    Yss = NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH); 
    Yss = Yss | (cat_vol_morph(Yss,'d',vxv*2) & Ym>2.25/3 & Ym<2.75/3 & Ydiv>-0.01);  % add higher tissue around mask
    Yss = Yss | (cat_vol_morph(Yss,'d',vxv*3) & NS(Yl1,LAB.VT) & Yp0>1.5 & Yp0<2.3);  % add lower  tissue around mask
    Yss = Yss & Yp0>1.5 & (Yp0<2.75 | (Ym<(2.5 + LASstr * 0.45)/3 & Ydiv>-0.05));     % limit the enlarged map by intensity
    Yss = Yss | ((Yxd./max(eps,Ytd + Yxd))>THth/2 & ... % save TH by distances - for overcorrected images
                 (Yp0<2.75 | (Ym<(2.75 + LASstr * 0.20)/3 & Ydiv>-0.05)));  
    Yss = cat_vol_morph(Yss,'o'); 
    
    Ynw = (Yxd./max(eps,Ytd + Yxd))>THth/2 | (NS(Yl1,LAB.BG) & Ydiv>-0.01);
    if debug==0, clear Ytd Yxd ; end
    
    
    % increase CSF ROI
    % Yvt .. improve ventricle ROI with atlas information
%%%%% 20180801 - Why not use the ventricle maps from the partitioning? 
    %          - Because the partitioning is after LAS
    %          - Why not use a fast partitioning before LAS?
    %            Why not use the code from the partitioning? 
    %            Is a good ventricle segmentation here not important?
    Yvt = cat_vol_morph( (NS(Yl1,LAB.VT) | cat_vol_morph(Ycm,'o',3) ) ...
      & Ycm & ~NS(Yl1,LAB.BG) & ~NS(Yl1,LAB.TH) & Ybd>30,'d',vxv*3) & ~Yss; 
    
    
    
    %% meningeal corrections 
    %  Ycx  .. cerebellar and other CSF
    %  Ycwm .. cerebellar WM
    %  Yccm .. cerebellar CSF
    %  Ybwm .. brain WM
    %  Ybcm .. brain CSF
    Ycx = ( NS(Yl1,LAB.CB) & ( (Ym - Ydiv)<0.55 | Ycls{3}>128 ) ) | ...
          ( ( (Ym - Ydiv)<0.45 &  Ycls{3}>8 ) | Ycls{3}>240 );
    % in the cerebellum tissue can be differentiated by div etc.
    Ycwm = NS(Yl1,LAB.CB) & (Ym - Ydiv*4)>5/6 & Ycd>3 & Yg>0.05;
    Yccm = NS(Yl1,LAB.CB) & Ydiv>0.02 & Ym<1/2 & Yg>0.05;
    Ybwm = (Ym - Ydiv*4)>0.9 & Ycd>3 & Yg>0.05; 
    Ybcm = Ydiv>0.04 & Ym<0.55 & Yg>0.05;
    
    
    
    %% correction 1 of tissue maps
    % Ywmtpm .. refined WM template map (original Ywtpm)
    Ywmtpm = ( single( Ywtpm )/255 .* Ym .* (1 - Yg - Ydiv) .* ... 
      cat_vol_morph( NS(Yl1,1) .* Ybd/5 ,'e',1) )>0.6; % no WM hyperintensities in GM!
    if debug==0, clear Ywtpm; end
    
    
    % refine the GM map
    Ygm = Ygm | (Yss & ~Yvt & ~Ycx & ~Ybv2 & ~Ycwm & ~(Yccm | Ybcm));
    Ygm = Ygm & ~Ywmtpm & ~Ybvv; % no WMH area
    
    
    % refine the WM map 
    Ywm = Ywm & ~Yss & ~Ybv2 & ~Ynw;
    Ywm = Ywm | Ycwm | Ybwm;
    Ywmtpm(smooth3(Ywmtpm & Ym<11/12)<0.5) = 0; % 20180801 - Why here?
    Ywm = Ywm & ~Ywmtpm & ~Ybvv & ~Yss; % no WM area
    
    
    Ycm = Ycm | ( (Ycx | Yccm | Ybcm) & Yg<0.2 & Ym>0 & Ydiv>-0.05 & Ym<0.3 & Yb2 ) | Ybvv;
    if debug==0, clear Ycwm Yccm Ycd Ybvv Ycx; end
    
    
    % Mapping of the brain stem to the WM (well there were some small GM
    % structures, but the should not effect the local segmentation to much.
    % Ybs .. brain stem mask
    Ybs = Ym<1.2 & Ym>0.9 & Yp0>1.5 & ...
      cat_vol_morph(NS(Yl1,LAB.BS) & Ym<1.2 & Ym>0.9 & Yp0>2.5,'c',2*vxv);
    Ygm = (Ygm & ~Ybs & ~Ybv2 & ~Ywm) | Yss;
    Ywm = Ywm | (Ybs & Ym<1.1 & Ym>0.9 & Yp0>1.5) ; 
  end
  
  
  
  
  
  
  
  
  %% back to original resolution for full bias field estimation
  %  -------------------------------------------------------------------
  Ycls = Yclso; clear Yclso;
  Ycm  = cat_vol_resize(Ycm  , 'dereduceBrain' , BB); 
  Ygm  = cat_vol_resize(Ygm  , 'dereduceBrain' , BB); 
  Ywm  = cat_vol_resize(Ywm  , 'dereduceBrain' , BB); 
  
  Yvt  = cat_vol_resize(Yvt  , 'dereduceBrain' , BB); 
  Yb2   = cat_vol_resize(Yb2 , 'dereduceBrain' , BB); 
  Yss  = cat_vol_resize(Yss  , 'dereduceBrain' , BB); 
  Ybb  = cat_vol_resize(Ybb  , 'dereduceBrain' , BB); 
  Ybs  = cat_vol_resize(Ybs  , 'dereduceBrain' , BB); 
  Ybv2 = cat_vol_resize(Ybv2 , 'dereduceBrain' , BB); 
  
  Ym   = cat_vol_resize(Ym   , 'dereduceBrain' , BB); 
  Yp0  = cat_vol_resize(Yp0  , 'dereduceBrain' , BB); 
  Yl1  = cat_vol_resize(Yl1  , 'dereduceBrain' , BB); 
  Ybd  = cat_vol_resize(Ybd  , 'dereduceBrain' , BB);
  Yg   = cat_vol_resize(Yg   , 'dereduceBrain' , BB);
  Ydiv = cat_vol_resize(Ydiv , 'dereduceBrain' , BB);
  Ywd  = cat_vol_resize(Ywd  , 'dereduceBrain' , BB); 
 
  
  % correction for negative values 
  [Ysrcx,thx] = cat_stat_histth(Ysrc,99); clear Ysrcx;  %#ok<ASGLU>
  srcmin = thx(1); 
  Ysrc = Ysrc - srcmin; 
  T3th = T3th - srcmin; 
  Tthc = Tth; Tthc.T3th = Tth.T3th - srcmin;
  if exist('Ysrco2','var'),Ysrco2 = Ysrco2 - srcmin; end
  
  
  
  
  
  
  
%% --------------------------------------------------------------------- 
%  Now, we can estimate the local peaks 
%  ---------------------------------------------------------------------
%  Estimation of the local WM intensity with help of intensity corrected
%  GM, CSF and head tissues to avoid overfitting (eg. BWP cerebellum). 
%  CSF is problematic in high contrast or skull-stripped images and  
%  should not be used here or in the GM peak estimation.
%  ---------------------------------------------------------------------
  mres  = 1.1; 
  stime = cat_io_cmd('  Estimate local tissue thresholds (WM)','g5','',verb,stime); 
  Ysrcm = cat_vol_median3(Ysrc.*Ywm,Ywm,Ywm); 
  rf    = [10^5 10^4];
  T3th3 = max(1,min(10^6,rf(2) / (round(T3th(3)*rf(1))/rf(1))));
  Ysrcm = round(Ysrcm*T3th3)/T3th3;
  
  
  % Ygwg .. large homogen GM regions (e.g. subcortical structures or in the BWP cerebellum)
  % Ygwc .. large homogen CSF regions
  Ygwg  = Ycls{1}>128 & Ym>(2/3 - 0.04) & Ym<(2/3 + 0.04) & Ygm .*Ydiv>0.01;
  Ygwg  = Ygwg | (Ycls{1}>128 & Yg<0.05 & abs(Ydiv)<0.05 & ~Ywm & Ym<3/4); 
  Ygwc  = Ycls{3}>128 & Yg<0.05 & ~Ywm & ~Ygm & Ywd<3; 
  Ygwc(smooth3(Ygwc)<0.5) = 0;
  if debug==0, clear Ywd; end
  
  
  % also use normalized GM tissues 
  % Ygi .. tissue-corrected intensity map (normalized for WM)
  Ygi = Ysrc .* Ygwg * T3th(3)/mean(Ysrc(Ygwg(:)));
  if cat_stat_nanmean(Ym(Ygwc))>0.1  
    % use normalized CSF tissue only in images with clear CSF intensity 
    % (to avoid errors in skull-stripped data)
    Ygi = Ygi + Ysrc .* Ygwc * T3th(3)/mean(Ysrc(Ygwc(:))); 
  end
  
  
  %% limit image resolution for fast processing
  [Yi,resT2] = cat_vol_resize( Ysrcm ,'reduceV',vx_vol,mres,32,'max');    % maximum reduction for the WM area
  Ygi        = cat_vol_resize( Ygi   ,'reduceV',vx_vol,mres,32,'meanm');  % mean for other (intensity normalized) tissues 
  for xi = 1:2*LASi, Ygi = cat_vol_localstat(Ygi,Ygi>0,2,1); end;         % intensity smoothing in the GM area
  Ygi(smooth3(Ygi>0)<0.3) = 0;                                            % remove small dots (meninges)
  Yi = cat_vol_localstat(Yi,Yi>0,1,3);                        % maximum filtering to stabilize small WM structures
  Yi(Yi==0 & Ygi>0) = Ygi(Yi==0 & Ygi>0);                     % combine the WM and GM map
  for xi = 1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end % intensity smoothing in both regions no maximum here!
  
  
  % Add head tissue if it is available. 
  if cat_stat_nanmean(Ym(Ygwc))>0.1 && cat_stat_nanmean(Ysrc(Ycls{6}(:)>128))>T3th(1)
    % definion of head tissue
    Ygwh = Ycls{6}>240 & ~Ygwc & ~Ygwg & Yg<0.5 & abs(Ydiv)<0.1 & ...
      Ysrc>min( res.mn(res.lkp==6) * 0.5 ) & Ysrc<max( res.mn(res.lkp==6) * 1.5 ); 
    % go to low resolution
    [Ygi,resTB] = cat_vol_resize(Ysrc .* Ygwh * T3th(3)/mean(Ysrc(Ygwh(:))),'reduceV',vx_vol,mres,32,'meanm');
    % additional approximation
    Ygi = cat_vol_approx(Ygi,'nh',resTB.vx_volr,2); Ygi = cat_vol_smooth3X(Ygi,2 * LASfs);
    Ygi = Ygwh .* max(eps,cat_vol_resize(Ygi,'dereduceV',resTB)); 
    Yi(Yi==0) = Ygi(Yi==0); 
    if debug==0; clear Ygwh; end
  end
  if debug==0, clear Ygwg Ygwc; end
  
  
  % final approximation of the WM inhomogeneity
  Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); 
  Yi = cat_vol_smooth3X(Yi,2 * LASfs); 
  Ylab{2} = max(eps,cat_vol_resize(Yi,'dereduceV',resT2)); 
% Ylab{2} = Ylab{2} .* mean( [median(Ysrc(Ysw(:))./Ylab{2}(Ysw(:))),1] ); 
  if debug==0; clear Ysw Yi; end

  
  
  
  
  
  %% update GM tissue map after the WM bias correction
  %  --------------------------------------------------------------------
  Ygm(Ysrc./Ylab{2}>(T3th(2) + 0.90 * diff(T3th(2:3)))/T3th(3)) = 0; 
  Ygm(Ysrc./Ylab{2}<(T3th(2) + 0.75 * diff(T3th(2:3)))/T3th(3) & ...
      Ysrc./Ylab{2}<(T3th(2) - 0.75 * diff(T3th(2:3)))/T3th(3) & ...
      Ydiv<0.3 & Ydiv>-0.3 & Ybb & ~Ywm & ~Yvt & ~Ybv2 & Ycls{1}>48) = 1;
  Ywd2 = cat_vbdist(single(Ywm),Yb2); % update WM distance
  Ygmn = (Ywd2 - Ym + Ydiv)>0.5 & ~Ybv2 & ...
    (Ym + 0.5 - Ydiv - Yg - Ywd2/10)>1/3 & ... low intensity tissue
    ~(Ym - min(0.2,Yg + Ywd2/10 - Ydiv)<1/4) & ...
    Yg<Ylab{2}/T3th(3)*0.3 & Ysrc<Ylab{2}*0.9; % no real csf
  Ygmn(smooth3(Ygmn)<0.5) = 0; 
  Ygm  = ~Ywm & (Ygm | Ygmn); % correct gm (intensity based)
  Ygmn = (single(Ycls{1})/255 - abs(Ydiv) + min(0,Ydiv) - Yg)>0.5 & ~Ywm;
  Ygmn(smooth3(Ygmn)<0.5) = 0; 
  Ygm = Ygm | Ygmn; % correct gm (SPM based)
  if debug==0; clear Ygm4; end
  
  
  
  
  
  
  %% update maps
%%%  ... more details  

  %  Update CSF maps.
  %  Ycm  .. updated CSF map
  %  Ycp1 .. typical CSF
  %  Ycp2 .. venes
  %  Ycp3 .. sulcal CSF
  %  Ycp4 .. save non-CSF
  %  Ycp5 .. but do not trust the brain mask!
  %  Ycpn .. updated CSF 2 map
  %  Ycd2 .. updated CSF distance 
  Ycm  = Ycm | (Yb2 & (Ysrc./Ylab{2})<((T3th(1)*0.5 + 0.5*T3th(2))/T3th(3))); 
  Ycm  = ~Ygm & ~Ywm & ~Ybv2 & Yg<0.6 & Ycm; 
  Ycp1 = Ycls{2}<128 & Ydiv>0 & Yp0<2.1 & Ysrc./Ylab{2}<mean(T3th(1)/T3th(2));            % typical CSF
  Ycp2 = Ycls{5}>32 & Ycls{2}<32 & Ysrc./Ylab{2}<T3th(2)/T3th(3) & Ydiv>0;                % venes
  Ycp3 = (Ym - Ydiv<0.4) & Ycls{3}>4 & Ycls{3}>16 & Ysrc./Ylab{2}<mean(T3th(2)/T3th(3));  % sulcal CSF
  Ycp4 = (single(Ycls{6}) + single(Ycls{5}) + single(Ycls{4}))>192;                       % save non-CSF 
  Ycp5 = Ysrc./Ylab{2}<T3th(1)/T3th(3);                                                   % but do not trust the brain mask!
  Ycpn = Ycp1 | Ycm | Ycp2 | Ycp3 | Ycp4 | Ycp5;
   if debug==0; clear Ycp1 Ycp2 Ycp2 Ycp3 Ycp4; end
  Ycpn( smooth3(Ycpn)>0.4 ) = 1; % remove some meninges
  Ycd2 = cat_vbdist(single(Ycpn),~Ycpn,vx_vol);  
 
  
  %% update GM map
  Ygm = Ygm & ~Ycm & ~Ywm & Ywd2<5; 
  if debug==0; clear Ywd2; end
  Ygm = Ygm | (NS(Yl1,1) & Ybd<20 & (Ycd2 - Ydiv)<2 & ...
    Ycls{1}>0 & ~Ycm & Ybb & Ym>0.6 & Yg<max(0.5,1-Ybd/30)); 
  if debug==0; clear Ybd Ycd2; end
  Ygm = Ygm & (Yg<0.1 | Ysrc./Ylab{2}<(T3th(2)*1/3+2/3*T3th(3))/T3th(3));  % outer high intensity GM
  % update by special atlas based maps
  if LABl1
    Ygm = (Ygm | Yss) & ~Ycm & cat_vol_morph(~Ybs | ~Yvt,'e');
  end
  if debug==0; clear Yvt Yss; end
  Ybb = cat_vol_morph(smooth3(Ygm | Ywm | Yp0>1.5 | ...
    (Ym>1.2/3 & Ym<3.1/3 & Yb2))>0.6,'lo',min(1,vxv)); 
  if debug==0, clear Yp0 Yvt; end
  Ygm(~Ybb) = 0; 
  Ygm(smooth3(Ygm)<0.3) = 0;
  Ygm(smooth3(Ygm)>0.4 & Ysrc./Ylab{2}>mean(T3th(1)/T3th(2)) & ...
    Ysrc./Ylab{2}<(T3th(2)*0.2+0.8*T3th(3))) = 1;
  if debug==0; clear Ybb Yl1; end 

  
  
  
  
  %% GM
  %  .. more details
  
  %  Yi .. main GM intensity map
  stime = cat_io_cmd('  Estimate local tissue thresholds (GM)','g5','',verb,stime); 
  Yi = Ysrc./Ylab{2} .* Ygm; 
  Yi = round(Yi * rf(2))/rf(2);
  Yi(Ybs) = Ysrc(Ybs)./Ylab{2}(Ybs) .* T3th(2)/T3th(3); 
  if debug==0; clear Ybs; end
  Yi = cat_vol_median3(Yi,Yi>0.5,Yi>0.5); % reduce artifacts
  
  % add CSF and CSF/GM areas (venes) to avoid overfitting
  % Ycmx .. venes / sinus 
  % Tcmx .. intensity correction 
  Ycmx = smooth3( Ycm & Ysrc<(T3th(1)*0.8 + T3th(2)*0.2) )>0.9; 
  Tcmx = mean( Ysrc(Ycmx(:)) ./ Ylab{2}(Ycmx(:)) ) * T3th(3);
  Yi(Ycmx) = Ysrc(Ycmx)./Ylab{2}(Ycmx)  .* T3th(2)/Tcmx; 
  if debug==0; clear Ycm Ycmx; end
  
  %% reduce data to lower resolution 
  [Yi,resT2] = cat_vol_resize( Yi              ,'reduceV',vx_vol,1,32,'meanm');
  Yii        = cat_vol_resize( Ylab{2}/T3th(3) ,'reduceV',vx_vol,1,32,'meanm');
  Ybgx       = cat_vol_resize( Ycls{6}>240     ,'reduceV',vx_vol,1,32,'meanm');
  for xi = 1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,3,1); end     % local smoothing
  Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2);                   % approximation of the GM bias 
  Yi = min(Yi,Yii*(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3));      % limit upper values to avoid wholes in the WM
  Yi(Ybgx) = Yii(Ybgx) * cat_stat_nanmean(Yi(~Ybgx(:)));          % use WM bias field in the background  
  if debug==0; clear Yii Ybgx; end
  Yi = cat_vol_smooth3X(Yi,LASfs);                                % final low res smoothing
  Ylab{1} = cat_vol_resize(Yi,'dereduceV',resT2).*Ylab{2};        % back to original resolution
  if debug==0; clear Yi; end
  
  
  % update CSF map
  Ycm = (single(Ycls{3})/255 - Yg*4 + abs(Ydiv)*2)>0.5 & ... 
    Ysrc<(Ylab{1}*mean(T3th([1,1:2]))/T3th(2));
  if debug==0; clear Ydiv; end
  Ycm(smooth3(Ycm)<0.5)=0;
  Ycm(Yb2 & cat_vol_morph(Ysrc<mean(T3th(1:2)),'o'))=1;
  
  
  
  
  
  %% CSF & BG 
  %  The differentiation of CSF and Background is not allway good. In some
  %  images with low intensity (especially for defacing/skull-stripping)
  %  the CSF can be similar to the background. 
  stime = cat_io_cmd('  Estimate local tissue thresholds (CSF/BG)','g5','',verb,stime); 
  Ynb   = smooth3( Ycls{6})>128 & Yg<0.05 & Ysrc<T3th(1); 
  Ynb   = cat_vol_morph(Ynb,'e',4*vxv); 
  if exist('Ygwh','var')
    Ynb = Ygwh & Ynb; 
  end  
  if debug==0, clear Ycpn; end
  Yc = round(Ysrc./Ylab{2} .* (smooth3(Ycm)>0.5) * rf(2))/rf(2); 
  Yx = round(Ysrc./Ylab{2} .* Ynb * rf(2))/rf(2);
  
  %% The correction should be very smooth (and fast) for CSF and we can use 
  %  much lower resolution. 
  [Yc,resT2] = cat_vol_resize(Yc,'reduceV',vx_vol,8,16,'min'); % only pure CSF !!!
  Yx         = cat_vol_resize(Yx,'reduceV',vx_vol,8,16,'mean');
  Yx(Yx<1) = 0; % no boundary voxel
  
  % CSF intensity show be higher than background intensity 
  Yx(Yc>0)=0; Yc(Yx>0)=0;
  if cat_stat_nanmean(Yx(Yx~=0)) < cat_stat_nanmean(Yc(Yc~=0))
    meanYx = min(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
    meanYc = max(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  else
    meanYc = min(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
    meanYx = max(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  end  
  stdYbc = mean([std(Yc(Yc(:)>0)),std(Yx(Yx(:)>0))]);
  
  % guaranty small difference between backgound and CSF intensity 
  Yxa = cat_vol_approx(Yx ,'nh',resT2.vx_volr,16); 
  Yca = cat_vol_approx(Yc + min(max( meanYx + stdYbc , meanYc - stdYbc ),...
    Yx .* meanYc/max(eps,meanYx)),'nh',resT2.vx_volr,16); 
  if debug==0, clear Yc Yx; end
  Yca = Yca * 0.7 + 0.3 * max( mean(Yca(:)) , T3th(1)/T3th(3) );
  % smoothing 
  Yxa = cat_vol_smooth3X( Yxa , LASfs * 2 ); 
  Yca = cat_vol_smooth3X( Yca , LASfs * 2 );
  % back to original resolution, combination with WM bias, and final smoothing
  Yca = cat_vol_resize(Yca,'dereduceV',resT2) .* Ylab{2};
  Yxa = cat_vol_resize(Yxa,'dereduceV',resT2) .* Ylab{2};
  Ylab{3} = cat_vol_smooth3X( Yca , LASfs * 2 );  
  Ylab{6} = cat_vol_smooth3X( Yxa , LASfs * 2 );
  if debug==0, clear Yxa Yca; end
  
  
  
  %% restore original resolution
  if exist('resT0','var')
    Ycls = Yclso2; clear Yclso2; 
    Ysrc = Ysrco2; clear Ysrco2; 
    for i=[1 2 3 6], Ylab{i}  = cat_vol_resize (Ylab{i} , 'dereduceV' , resT0 ); end
    Yb2 = cat_vol_resize( single(Yb2) , 'dereduceV' , resT0 )>0.5;
  end
  
  
  
  %% local intensity modification of the original image
  % --------------------------------------------------------------------
  cat_io_cmd('  Intensity transformation','g5','',verb,stime); 
  Yml = zeros(size(Ysrc)); 
  Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc - Ylab{2}) ./ max(eps,Ylab{2} - Ylab{3})) );
  Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc - Ylab{1}) ./ max(eps,Ylab{2} - Ylab{1})) );
  Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc - Ylab{3}) ./ max(eps,Ylab{1} - Ylab{3})) );
  Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc - Ylab{6}) ./ max(eps,Ylab{3} - Ylab{6})) );
  Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
  
  
  
  %% update of the global intensity normalized map
  if Tth.T3th(Tth.T3thx==1) == Tth.T3thx(Tth.T3thx==1) % invers intensities (T2/PD) 
    Ymg = max(eps,(Ysrc + srcmin)./Ylab{2}); 
  else
    Ymg = (Ysrc + srcmin)./max(eps,(Ylab{2} + srcmin)); 
    Ymg = Ymg * Tth.T3th(Tth.T3thx==3)/(Tth.T3thx(Tth.T3thx==3)/3);
  end
  clear Ylab Ysrc
  Ymg = cat_main_gintnorm(Ymg,Tth); 
  
 
  
  %% fill up CSF in the case of a skull stripped image 
  if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.3))
    YM   = cat_vol_morph(Yb2,'d'); 
    Ymls = smooth3(max(Yml,YM*0.5));
    Yml(YM & Yml<0.5) = Ymls(YM & Yml<0.5); 
    clear Ymls YM
  end
  clear res
  
  
  
  %% interpolate maps if an resolution adaption was applied
  if exist('resT0','var')
    Ywm = cat_vol_resize( single(Ywm) , 'dereduceV' , resT0 )>0.5;
    Ygm = cat_vol_resize( single(Ygm) , 'dereduceV' , resT0 )>0.5;
  end
  
  
  
  %% class correction and definition of the second logical class map Ycls2
  Ynwm = Ywm & ~Ygm & Yml/3>0.95 & Yml/3<1.3;
  Ynwm = Ynwm | (smooth3(Ywm)>0.6 & Yml/3>5/6); Ynwm(smooth3(Ynwm)<0.5)=0;
  Yngm =  Ygm & ~Ywm & Yml/3<0.95; Yngm(smooth3(Yngm)<0.5)=0;
  Yncm = ~Ygm & ~Ywm & ((Yml/3)>1/6 | Ycls{3}>128) & (Yml/3)<0.5 & Yb2;
  if debug==0, clear Ywm Ygm; end
  
  Yp0     = single( Ycls{1} )/255*2 + single( Ycls{2} )/255*3 + single( Ycls{3} )/255;  % recreate label map
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) + (Ynwm & ~Yngm & Yp0>=1.5)*256 - (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Ynwm & ~Yngm & Yp0>=1.5)*256 + (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) - ((Ynwm | Yngm) & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) + (Yb2 & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Yb2 & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) - (Yb2 & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls2 = {Yngm,Ynwm,Yncm};
  if debug==0, clear Yngm Ynwm Yncm Yb2; end
  
  
  
  %% final correction of Yml similar to Ymg 
  Yml = Yml/3;

end

