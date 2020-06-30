function [Ygmt,Ypp,Ymf,Ywmd,Ycsfdc] = cat_vol_pbt3(Ymf,opt)
% ______________________________________________________________________
%
% Cortical thickness and surface position estimation. 
%
% The PBT3 approach used the main concepts of the PBT approach in principle.
% In contrast to the classical PBT, it tries to focus on the GM intensities
% and not only the boundaries to describe the cortical ribbon to avoid the
% "racing line" that cut gyri and sulci. 
%
% Especially low quality data and subjects with thin structures (CSF in 
% sucli or WM in gyri) have badly described GM-WM and GM-CSF boundaries 
% and a central line around the average GM intensity works much better. 
% E.g. for a label map with CSF=1 and WM=3, an average GM threshold creates
% a very nice central intensity surface (CIS). Because this surface runs 
% somewhere in the GM and not directly usable but it helps to divide the GM 
% into two regions - one with lower intensities that needs a classical 
% suclus reconstruction and one with higher intensities that required a 
% gyrus reconstruction. PBT can now be used to estimate the thickness of the
% lower and the upper segment and map the data over the whole GM segment to
% sum up to the full cortical thickness. This also allows to create a real
% central surface.  
%
%   [Ygmt,Ypp,Ywmd,Ycsfd] = cat_vol_pbt(Ymf,opt)
%  
%   Ygmt ..      GM thickness map 
%   Ypp ..       percentage position map
%   Ywmd ..      WM distance map
% 
%   Ymf ..       tissue segment image or better the noise, bias, and 
%                intensity corrected 
%
%   opt ..       structure with options and parameters 
%    .resV       voxel resolution (only isotropic; given by cat_surf_createCS)
%    .enlarge    use larger boundary better projection (correct thickness)
%                (default = 1)
%    .verb       verbose level (default = cat_get_defaults('extopts.verb')-1)
%    .cb         cerebellar contrast enhancement (given by cat_surf_createCS) 
%    .fast       fast distance estimation for CIS 
% ______________________________________________________________________
%
%   Dahnke, R; Yotter R; Gaser C.
%   Cortical thickness and central surface estimation.
%   NeuroImage 65 (2013) 226-248.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_vol_pbt3.m 1520 2019-11-19 21:52:09Z dahnke $ 

%#ok<*ASGLU>

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

% default variables and check/set function  
  if ~exist('opt','var'), opt=struct(); end

  def.resV        = 1;
  def.enlarge     = 1; % use larger boundary better projection (correct thickness)
  def.debug       = cat_get_defaults('extopts.verb')>2;
  def.verb        = cat_get_defaults('extopts.verb')-1;
  def.cb          = 0; 
  def.fast        = 1;
  def.sharpening  = 1;
  def.wmbiascorr  = 0;
  def.gmbiascorr  = 0; % for myelination
  opt             = cat_io_checkinopt(opt,def);
  opt.resV        = mean(opt.resV);
  
  
  Ymfo=Ymf; 
  %% Create central intensity surface (CIS):
  %  Due to noise and artifacts we need to optimize the GM segment. 
  %  I tried to use some different smoothing filters but this create s
  if opt.verb, fprintf('\n'); end;  if debug, Ymf = Ymfo; tic; end
  stime = cat_io_cmd('    Create central intensity surface: ','g5','',opt.verb); 
  
  % cleanup - why not in createCS ??
  Ymf = Ymf .* cat_vol_morph( cat_vol_morph( Ymf>1.2,'ldo', 3 / opt.resV / (opt.cb+1) ), 'd', 0.5/opt.resV);
  Ymf = Ymf .* cat_vol_morph( cat_vol_morph( Ymf>1.5,'ldo', 2 / opt.resV / (opt.cb+1) ), 'd',  0.5/opt.resV);
  Ymf = Ymf .* cat_vol_morph( cat_vol_morph( Ymf>1.9 - 0.2*opt.cb,'ldo', 2 / opt.resV / (opt.cb+1) ), 'd',2 / opt.resV * (opt.cb+1) );
  Ymf = max(1,Ymf);
  
  Ymf = max( Ymf , smooth3(cat_vol_morph( cat_vol_morph( Ymf<2.9,'ldo', 3 / opt.resV / (opt.cb+1)), 'd',  0.5/opt.resV)<0.5) * 3);
  
 
  
  %% Remove outlier and artifacts:
  % This was already done in cat_surf_createCS and further smoothing may 
  % remove important anatomical information. 
  Ycis = Ymf; 
  Ycis = cat_vol_median3(Ycis,Ycis>1.5 & Ycis<2.5,true(size(Ycis)),0.05);
  
  % estimate diverence map for sharpening
  if opt.sharpening 
    if debug, tic; Ycis = Ymf;  end
    
    Ydiv = cat_vol_div( Ycis , opt.resV , 1 / (opt.cb+1) ); dth = [-0.5 0.8];
    Ydiv(Ydiv<dth(1)) = min(0,2*dth(1) - Ydiv(Ydiv<dth(1)));
    Ydiv(Ydiv>dth(2)) = max(0,2*dth(2) - Ydiv(Ydiv>dth(2)));
    
    if opt.cb
      Ymask = cat_vol_smooth3X(cat_vol_morph(cat_vol_morph(Ycis>1,'c',4),'e',4),2); 
      Yg = cat_vol_div( Ycis , opt.resV , 1 / (opt.cb+1));
      
      % Sharpening based on bias correction and the divergence map:
      % first bias correction (medium frequency)
      YMR = (Ycis .* ((( -Ydiv - Yg + 1).* Ymask).^3)>2); Ymw = Ymf .* YMR; 
      for i=1:round(2/opt.resV), Ymw = cat_vol_localstat(Ymw,YMR,1,3); end 
      Ymw(cat_vol_smooth3X(Ymf,16)<=1) = 3; 
      if ~debug, clear YMR; end
      Ymw = cat_vol_approx(Ymw/3,'nh',opt.resV,2); 
      Ycis = Ycis ./ Ymw;     
      % sharpening % high frequency)
      Ycis = Ycis .* ( -Ydiv ./ max(0.5,abs( min(3,max(1,Ycis))/1.5-0.75)) /2 + 1 ).^2; 
      Ycis = max(1,min(3,Ycis));
      Ymf  = Ymf ./ Ymw;
      if ~debug, clear Ymw; end
     
    else
      % sharpening % high frequency)
      Ycis = Ycis .* ( -Ydiv ./ max(0.5,abs( min(3,max(1,Ycis))/1.5-0.75)) /8 + 1 ).^2; 
      Ycis = max(1,min(3,Ycis));
      Ycis(Ymf<2 & Ycis>Ymf) = Ymf(Ymf<2 & Ycis>Ymf);
    end
    
    if debug, toc; end
  end

  

  
  %% Estimation of boundary maps:
  %  Use boundaries that are close to the average GM intensity do get more
  %  anatomical information. Although we can not further use these distance
  %  maps, we need high accuracy to obtain correct results that support 
  %  asymmetric structures. Hence, we have to use eidist rather then vbdist
  %  although this is slower.
  %  Even though, this is just an initial mid GM surface, it have to support 
  %  a good differentiation between the areas where we want a gyrus/sulcus
  %  reconstruction!

  % bias correction ###### this is not reqired here #######
  if opt.wmbiascorr
    if debug, tic; end
    YMR = (Ycis - Ydiv*5)>2.3 & Ycis>2.3 & Ycis<3.3 & Ydiv<0.05; YMR(smooth3(YMR)<0.45)=0;
    [YMR,resV] = cat_vol_resize(YMR  ,'reduceV', opt.resV, min(opt.resV*2, 1.4), 32, 'meanm');
    Ycisr      = cat_vol_resize(Ycis ,'reduceV', opt.resV, min(opt.resV*2, 1.4), 32, 'meanm');
    YMR = cat_vol_morph(YMR,'ldo',((1 - opt.cb)/2) / mean(resV.vx_volr) );
    YMR = cat_vol_morph(YMR,'ldc',0.5 / mean(resV.vx_volr) );
    Ymw = Ycisr .* YMR; 
    for i=1:round(3 / mean(resV.vx_volr) ), Ymw = cat_vol_localstat(Ymw,YMR,1,3); end
    for i=1:round(2 / mean(resV.vx_volr) ), Ymw = cat_vol_localstat(Ymw,YMR,1,1); end
    if ~debug, clear YMR; end
    Ymw(cat_vol_morph(Ycisr./Ymw < 2.5/3 , 'ldc',0.5/opt.resV)) = 0; clear Ycisr; 
    Ymw = cat_vol_approx(Ymw,'nh',opt.resV,2); 
    Ymw = cat_vol_resize(Ymw ,'dereduceV',resV); 
    Ycis = Ycis ./ Ymw * 3; 
    if ~debug, clear Ymw; end  
    if debug, toc; end
  end

  
  
  
  %% GM mean intensity for myelination 
  if opt.gmbiascorr && ~opt.cb
    if debug, tic; end
    Ycisth = [1.5 2.8]; % higher upper boundary for myelinated areas 
    Ymc = Ycis .* (Ycis>Ycisth(1) & Ycis<Ycisth(2) & cat_vbdist(single(Ycis<1.75))<5/opt.resV & Ydiv<0.5);
    Ypve = Ymc - cat_vol_localstat(max(2,Ymc),Ymc>0,1,2); 
    Ymc = Ymc .* (Ypve<0.7);
    Ymc = Ymc .* cat_vol_morph( cat_vol_morph( Ymc>0,'do', 1.5 ) , 'l' , [10 0.1]) ;
    [Ymc,resV] = cat_vol_resize(Ymc ,'reduceV', opt.resV, min(opt.resV*2, 1.4), 32, 'meanm');
    for i=1:round(2/opt.resV), Ymc = cat_vol_localstat(Ymc,Ymc>0,1,1); end
    Ymc = cat_vol_approx(Ymc,'nh',opt.resV,2);        % ######### noisy? staerker filtern und dafuer stabiler? #######
    Ymc = cat_vol_resize(Ymc ,'dereduceV',resV); 
    if debug, toc; end
  else
    Ymc = 2; 
  end
  
  
  
  
  %% intensity normalization 
  if 0 & opt.cb, tb = [1.9 2.2]; else, tb = [1.9 2.2]; end % 1.9 2.8 
  Ypp = zeros(size(Ymf)); 
  Ypp = Ypp + ( (Ymf>=tb(2) ) );
  Ypp = Ypp + ( (Ymf>=Ymc & Ymf<tb(2) ) .* (0.5 + 0.5 * (Ymf - Ymc  ) ./ max(eps,tb(2) - Ymc)) );
  Ypp = Ypp + ( (Ymf<Ymc              ) .* (      0.5 * (Ymf - tb(1)) ./ max(eps,Ymc - tb(1))) );
  Ypp(isnan(Ypp) | Ypp<0)=0; Ypp(Ypp>1)=1;
  Ypp = cat_vol_median3(Ypp,Ypp>0 & Ypp<1,true(size(Ymf)),0.05);
  %Ypp = 1 - Ypp; %cat_vol_smooth3X(Ypp,0.5); 
  
   [Yt,Ypp1] = cat_vol_pbt2(1 + 2*Ypp,struct('method','pbt2','resV',opt.resV,'vmat',opt.vmat)); % avoid underestimated thickness in gyri
  
  % Estimate wider maps of the major tissues to create the mid surface: 
  %  WM
  if debug, tic; end
  Ywm = Ypp>0.95;
  Ywm(smooth3(Ywm)<0.3) = 0; 
  Ywm(smooth3(Ywm)>0.8) = 1; 
  Ywm = cat_vol_morph(Ywm,'l');
  Ywm = ~cat_vol_morph(~Ywm,'l'); 
  
  %  CSF
  Ycm = Ypp<0.05; %(Ycisc-min(0,Ydiv*2))<2.1 & (Ycis-min(0,Ydiv*2.1))<2 & Ymc<2.3; if ~debug, clear Ymc; end  
  Ycm(smooth3(Ycm)>0.6) = 1; Ycm(smooth3(Ycm)<0.4) = 0; 
  Ycm = ~cat_vol_morph(~Ycm,'lo',0.5/opt.resV) & ~Ywm & smooth3(Ycis)<2; 
  

  
  % use Laplace-filter to remove blood-vessels 
  % (thin structures lose more engergy) 
  %{
  Ylt = (1 - Ycm + Ywm)/2; 
  Ylt = cat_vol_laplace3R(single(Ylt),Ylt>0 & Ylt<1,0.05);
  Ycm = Ycm | smooth3(Ylt<0.1.*Ymf)>0.25;
  
  % update Ywm distance map
  Ywm = (Ywm & smooth3(Ylt<0.1.*Ymf)<0.25);
  Ymf = min(2 + 1.5*smooth3(Ylt .* max(0,Ymf-2)),Ymf); %smooth3(abs(Ylt - (Ymf/3))>2/3));  
  %}
  
  if debug, toc; else, clear Ycisc; end
 
 
  
  %%
  % Speed maps for asymmetric mapping and distance estimation 
  % ######### puh - eidist 120s vs. vbdist 10s ##############
  if debug, tic; end
  if opt.fast
    Ycd  = cat_vbdist(single(Ycm),Ywm<1) - 0.5; 
    Ywd  = cat_vbdist(single(Ywm),Ycm<1) - 0.5; 
  else
    Ycm  = single(Ycm); Ycm(Ywm>0) = nan; 
    Ywm  = single(Ywm); Ywm(Ycm>0) = nan;
    F    = max(eps,max(eps,smooth3(min(1,(Ycis-1).^4)))); 
    Ywd  = cat_vol_eidist(Ywm,F); Ywd(isnan(Ywd))=0; %if ~debug, clear Ywm; end
    F    = smooth3(max(eps,smooth3(max(eps,min(1,(3-Ycis)).^4))));
    Ycd  = cat_vol_eidist(Ycm,F); Ycd(isnan(Ycd))=0; %if ~debug, clear Ycm; end
  end
  if debug, toc; end
 
  
  % PBT mapping ###### this are just 10 sekunden ####
  %  ~30s
  if debug, tic; end % we need doubles (copies) for pbtv!  
  Ygt1 = cat_vol_pbtp( single(2 - (Ycm>0) + (Ywm>0) + (Ycm>20)) , Ywd , Ycd );
  Ygt1 = min(Ygt1, Ycd + Ywd); % if opt.fast, Ygt1(Ygt1==0) = Ycd(Ygt1==0) + Ywd(Ygt1==0); end  % bei fast ging hier was schief wegen Ycis ...
  for i=1:round(5 / mean(opt.resV) ), Ygt1 = cat_vol_localstat(Ygt1,Ygt1>0,1,1); end
  if debug, toc; else, clear Ycd; end 
  %
  % Create volume/distance-based percentage position map (Yppv/Yppd): 
  % The Yppv map is strongly affected by (interpolation) artifacts and bad 
  % structures (broken gyri) need a lot of smoothing to avoid step artifacts. 
  % However, it support a meander line with wider gyri and thinner sulci.
  % The distance metric creates a smoother Ypp than the volume metric that 
  % helps to reduce artifacts but also show less anatomical details and 
  % use the "racing line". Hence, we combine both maps. 
  % This map has to run perfectly somewhere between the white and pial  
  % boundary and should not contain artifacts
  if debug, tic; end
  Ypp = (Ygt1 - Ywd) ./ max(eps,Ygt1); Ypp(Ywm>0)=1; 
  Ypp = cat_vol_median3(Ypp,Ypp>0.1 & Ypp<0.9,true(size(Ypp)),0.1); 
  Ypp(cat_vol_morph(Ypp>=0.5,'ldo',max(0.5,(1.5 - opt.cb)/(2*opt.resV)))==0 & Ypp>=0.5) = 0.49; %# ##############
  Ypp(Ypp<0.5 & ~cat_vol_morph(Ypp<0.5,'l')) = 0.51;
  if ~debug, clear Yv1 Yv2 Ywd; else, toc; end
  if ~debug
    clear Ygt1 Ygt; 
  end
  Ypp = 1 - Ypp; 

  
  
  %% distance from the central intensity surface (CIS)
  %  ----------------------------------------------------------------------
  %  Distance from the CIS to the CSF and WM. For the CSF case the removal 
  %  of meninges and blood-vessels is maybe important and depends on the 
  %  distance estimation. vbdist require removal, whereas eidist may works 
  %  better without and a good speed map F.
  stime = cat_io_cmd('    GM distance estimation: ','g5','',opt.verb,stime); if debug, tic; end
  
  % CIS boundary and area of distance estimation:
  YM = min(1,max(0,((cat_vol_smooth3X(Ypp,0.3/opt.resV) .* Ymf) - mean(Ymf(Ypp(:)>0.2 & Ypp(:)<0.8).*Ypp(Ypp(:)>0.2 & Ypp(:)<0.8)))/2 + 0.5));         % create a PVE bounary for the CIS
  YM(cat_vol_morph(YM<0.5,'e',0.5/opt.resV)) = 0;                          % remove lower PVE voxel 
  YM(cat_vol_morph(YM>0.5,'e',0.5/opt.resV)) = 1;                          % remove upper PVE voxel
  YM(~cat_vol_morph(YM<0.99,'l',0.5/opt.resV)) = 1;                        % remove small wholes
  YMM = ~cat_vol_morph(YM>0.1,'lo',0.25/opt.resV) & YM>0; YM(YMM) = 0;     % remove small objects (meninges) 
  
  %  CIS distance to the CSF and correction:
  %  We assume that the Eucledean and Eikonal are roughly the same. If the 
  %  Eikonal distance is much greater than it typically describes meninges
  %  or blood-vessels.
  % enlarge the estimation region to improve the thickness estimation - needs later correction! 
  YM( smooth3(cat_vol_morph( smooth3(~cat_vol_morph(Ymf>1.5,'lo',1/opt.resV))>0.3 | ...
    YMM ,'e',opt.enlarge))>0.5 ) = nan; clear YMM;
  F = smooth3(max(eps,max(eps,min(1,max(0.5,Ymf-0.5)/1.55)))); F(isnan(YM)) = nan;             % speed map for eikonal distance
  [Ycid,t,Yd] = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug); clear t;    
  YM(Yd>1.4*Ycid & Ycid>1) = nan; Ycid(Yd>1.4*Ycid & Ycid>1) = nan; clear Yd;   % correct for meninges
  if opt.enlarge
    % correc distance values behind the boundary
    YM2   = max(YM,max(0,min(1,Ymf-1))); YM2(isnan(F)) = nan; 
    Ycidc = cat_vol_eidist(YM2,F,[1 1 1],1,1,0,opt.debug); 
    Ycid  = Ycid - Ycidc; clear Ycidc YM2;                            
  end  
  if debug, toc; end

  %  CIS distance to the CSF and correction:
  %  As far as we can use the YM boundary description this part is much shorter. 
  if debug, tic; end
  F   = smooth3(max(eps,max(eps,min(1,(3-Ymf)).^2)));                      % speed map for eikonal distance
  YMM = cat_vol_morph(Ymf>2.5,'e',opt.enlarge) | isnan(Ymf);               % change YMM to have the deep WM
  YMI = YM; YMI(isnan(YMI))=0; YMI = 1 - YMI; YMI(YMM & YM>0.9) = nan; F(isnan(YM)) = nan;       % create inverse boundary map
  [Ycid2,t,Yd] = cat_vol_eidist(YMI,F,[1 1 1],1,1,0,opt.debug); clear t YMM;  
  Ycid2(Yd>1.4*Ycid2 & Ycid2>1) = nan; clear Yd;
  if opt.enlarge
    % correct distance values behind the boundary
    YM2   = max(YM,max(0,min(1,3-Ymf))); YM2(isnan(F)) = nan;
    Ycidc = cat_vol_eidist(YM2,F,[1 1 1],1,1,0,opt.debug); 
    Ycid2  = Ycid2 - Ycidc; clear Ycidc;                            
  end
  clear F;
   
  % combine both distance maps
  Ycid = Ycid - Ycid2; YnGM = isnan(Ycid); Ycid(isnan(Ycid)) = 0; clear Ycid2;
  Ycid = cat_vol_median3(Ycid,Ycid>0,Ycid>0,0.1); 
  if debug, toc; end

  
  %% intensity normalization 
  Ymp = zeros(size(Ymf)); tb = [1.0 3.0];
  Ymp = Ymp + ( (Ymf>=tb(2) ) .* 2);
  Ymp = Ymp + ( (Ymf>=Ymc & Ymf<tb(2) ) .* (1 + (Ymf - Ymc  ) ./ max(eps,tb(2) - Ymc)) );
  Ymp = Ymp + ( (Ymf<Ymc              ) .* (    (Ymf - tb(1)) ./ max(eps,Ymc - tb(1))) );
  Ymp(isnan(Ymp) | Ymp<0)=0; Ymp(Ymp>2)=2;
  Ymp = cat_vol_median3(Ymp,Ymp>0 & Ymp<2,true(size(Ycis)),0.05);
  Ymp = cat_vol_smooth3X(Ymp,0.5); 
  
  
  %% CSF and WM distance
  stime = cat_io_cmd('    CSF and WM distance estimation: ','g5','',opt.verb,stime); if debug, tic; end
  F     = smooth3(max(eps,max(eps,min(1,(2-Ymp)).^4)));                    % speed map for eikonal distance
  YM3   = max( max(0,min(1,1-Ymp)) , isnan(YM) );  % 2-Ymf
  YM3(YM3==0 & YnGM) = nan; F(isnan(YM3))=nan; if ~debug, clear Ycid2; end
  Ycsfd = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); Ycsfd(Ycsfd>100) = 0; 
  if opt.enlarge
    % correc distance values behind the boundary
    YM3    = max(YM3,max(0,min(1,2-Ymp))); % 3-Ymf
    YM3(cat_vol_morph(YM3>0.1,'e')) = 1; YM3(cat_vol_morph(YM3<0.9,'e')) = 0; 
    YM3((YM3==0 & YnGM) | isnan(F)) = nan; 
    Ycsfdc = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); 
    Ycsfd  = Ycsfd - Ycsfdc; clear Ycsfdc;                          
  end
  % WM distance 
  F     = smooth3(max(eps,max(eps,min(1,(Ymp).^4)))); %Ymf-1
  YM3   = max(max(0,min(1,(Ymp-1))) ,  isnan(YMI) ); YM3(YM3==0 & YnGM) = nan; F(isnan(YM3))=nan; %Ymf-2
  Ywmd  = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); 
  if opt.enlarge
    % correct distance values behind the boundary
    YM3   = max(YM3,max(0,min(1,Ymp))); YM3((YM3==0 & YnGM) | isnan(F)) = nan; %Ymf-1
    Ywmdc = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); 
    Ywmd  = Ywmd - Ywmdc; clear Ywmdc;                             
  end
  if debug, toc; else, clear YM3 F; end
  
  
  %% projection-based thickness estimation 
  stime = cat_io_cmd('    Thickness projection: ','g5','',opt.verb,stime); if debug, tic; end
  Ygmt1 = abs(cat_vol_pbtp( round(4 - Ymf) , -Ycid , Ywmd  )); 
  Ygmt2 = abs(cat_vol_pbtp( round(Ymf) ,  Ycid , Ycsfd ));
  Ygmtp = abs(cat_vol_pbtp( round(Ymf) ,  Ywmd , Ycsfd ));
  Ygmtd = Ycsfd + Ywmd;

  %% median filter
  fs    = 0; %min(1,max(0,1.5 - opt.resV));
  YM    = Ymf>1.5 & Ymf<2.5; 
  Ygmt1 = Ygmt1*(1-fs) + fs*cat_vol_median3(Ygmt1,YM,Ygmt1>0,0.05);
  Ygmt2 = Ygmt2*(1-fs) + fs*cat_vol_median3(Ygmt2,YM,Ygmt2>0,0.05);
  Ygmtp = Ygmtp*(1-fs) + fs*cat_vol_median3(Ygmtp,YM,Ygmtp>0,0.05);
  Ygmtd = Ygmtd*(1-fs) + fs*cat_vol_median3(Ygmtd,YM,Ygmtd>0,0.05);

  %% need this only for the next filter step to avoid crossing in sulci
  Ygmt    = min( Ygmt2 + Ygmt1 , min(Ygmtp,Ygmtd) );
  [D,I]   = cat_vbdist(single(Ygmt>0),Ymf>1.5 & Ymf<2.5); Ygmt = Ygmt(I); clear D I;   
  YM      = Ygmt>0; 
  Ypp     = zeros(size(Ymf),'single'); Ypp(Ymf>=2.5) = 1;
  Ypp(YM) = min(Ycsfd(YM),Ygmt(YM) - Ywmd(YM)) ./ (Ygmt(YM) + eps); Ypp(Ypp>2) = 0;
  clear Ygmt
  
  %% local filtering
  Ygmt1s = cat_vol_localstat(Ygmt1,YM & Ypp>0.1,1,1); Ygmt1(Ygmt1s>0) = Ygmt1(Ygmt1s>0)*(1-fs) + fs*Ygmt1s(Ygmt1s>0); clear Ygmt1s
  Ygmt2s = cat_vol_localstat(Ygmt2,YM & Ypp>0.1,1,1); Ygmt2(Ygmt2s>0) = Ygmt2(Ygmt2s>0)*(1-fs) + fs*Ygmt2s(Ygmt2s>0); clear Ygmt2s
  Ygmtps = cat_vol_localstat(Ygmtp,YM & Ypp>0.1,1,1); Ygmtp(Ygmtps>0) = Ygmtp(Ygmtps>0)*(1-fs) + fs*Ygmtps(Ygmtps>0); clear Ygmtps
  Ygmtds = cat_vol_localstat(Ygmtd,YM & Ypp>0.1,1,1); Ygmtd(Ygmtds>0) = Ygmtd(Ygmtds>0)*(1-fs) + fs*Ygmtds(Ygmtds>0); clear Ygmtds
 
  %% avoid meninges !
  Ygmt  = min( Ygmt2 + Ygmt1 , min(Ygmtp,Ygmtd) );
  %%
  if ~debug, clear Ygmt2 Ygmtp Ygmtd; end 
  
  % median filter to remove outliers
  Ygmt = cat_vol_median3(Ygmt,Ygmt>0,Ygmt>0,0.05);

  % filter result
  Ygmts = cat_vol_localstat(Ygmt,(Ygmt>1 | Ypp>0.1) & Ygmt>0 & (Ygmt>1 | Ymf>1.8),1,1); Ygmt(Ygmts>0) = Ygmts(Ygmts>0);
  clear Ypp; 
  
  
  %% Create percentage position map Ypp
  Ycd = min(Ycsfd,Ygmt - min(Ywmd,Ygmt1 + Ycid) );
  Ypp = max(0, min(1 , ( Ycd .* (Ygmt>0 & Ycsfd>0) ) ./ max(eps,Ygmt) )); 
  Ypp(Ypp==0 & Ymf>=2.5 & Ywmd<=0) = 1; 
  Ypp(Ypp==0 & Ymf>=2.5 & Ywmd<1)  = 0.99; 
  Ypp(Ypp==0 & Ymf>=2.0 & Ywmd<2)  = 0.985; 
  Ypp(Ypp==0 & Ymf>=2.5)           = 0.95; 
  
  Ypp = Ypp*(1-fs) + fs*cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & (Ypp>0.1 | Ypp==0),true(size(Ypp)),0.1);
  Ypp = Ypp*(1-fs) + fs*cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & (Ypp>0.1 | Ypp==0),true(size(Ypp)),0.1);
  
  %Ypp = cat_vol_laplace3R(Ypp,Ypp>0.05 & Ypp<0.95,0.4);
  %Ypp = cat_vol_laplace3R(Ypp,Ypp>0.35 & Ypp<0.65,0.4);
  Ypp(Ypp>0.5 & ~cat_vol_morph(Ypp>0.5,'ldo',0.50/opt.resV)) = 0.49;
  Ypp(Ypp<0.5 & ~cat_vol_morph(Ypp<0.5,'l'  ,0.25/opt.resV)) = 0.51;
  if debug, toc; else, clear Ygmt1 Ycid Ycsfd; end

  
 %% Final corrections for position map with removing of non brain objects.
  % ds('d2','',1,Ymf/3,Ywmd/3,Ygmt/5,Ypp,70)
  
  % Final corrections for thickness map with thickness limit of 10 mm. 
  % Resolution correction of the thickness map after all other operations, 
  % because PBT actually works only with the voxel-distance (isotropic 1 mm)
  Ygmt = Ygmt * opt.resV; 
  Ygmt(Ygmt>10) = 10; 
  if debug, toc; end
  if exist('Ycsfdc','var'), Ycsfdc = Ycsfdc*opt.resV; end
  if exist('Ywmd','var'), Ywmd = Ywmd*opt.resV; end
  cat_io_cmd(' ','g5','',opt.verb,stime);
  if opt.debug, cat_io_cmd(' ','','',opt.debug,stime2); end
 
  [D,I] = cat_vbdist(single(Ygmt>0),cat_vol_morph(Ymf<1.5 | Ymf>2.5,'e',2)); Ygmt = Ygmt(I); clear D I;  
   
  if 0
    %% just display - use sd to rotate the surface
    [D,I] = cat_vbdist(single(Ygmt>0),cat_vol_morph(Ypp<0.5,'d') & cat_vol_morph(Ypp>0.5,'d')); Ygmt = Ygmt(I); clear D I; 
    sd = 1; cat_surf_render2('Disp',isosurface(shiftdim(Ypp,sd),0.5,shiftdim(Ygmt,sd)));
  end
  
  
end
