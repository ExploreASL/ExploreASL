function [Ygmt,Ypp,Ywmd] = cat_vol_pbt3(Ymf,opt)
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
% $Id: cat_vol_pbt3.m 1311 2018-04-26 08:16:03Z dahnke $ 

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

% default variables and check/set function  
  if ~exist('opt','var'), opt=struct(); end

  def.resV        = 1;
  def.enlarge     = 2; % use larger boundary better projection (correct thickness)
  def.debug       = cat_get_defaults('extopts.verb')>2;
  def.verb        = cat_get_defaults('extopts.verb')-1;
  def.cb          = 0; 
  def.fast        = 0; 
  opt             = cat_io_checkinopt(opt,def);
  opt.resV        = mean(opt.resV);
  
  %% Create central intensity surface (CIS):
  %  Due to noise and artifacts we need to optimize the GM segment. 
  %  I tried to use some different smoothing filters but this create s
  if opt.verb, fprintf('\n'); end
  stime = cat_io_cmd('    Create central intensity surface: ','g5','',opt.verb); 
  
  % Remove outlier and artifacts:
  % This was already done in cat_surf_createCS and further smoothing may 
  % remove important anatomical information. 
  Ycis = Ymf; 
  Ycis = cat_vol_median3(Ycis,Ycis>1.5 & Ycis<2.5,true(size(Ycis)),0.05);

  %% estimate diverence map for sharpening
  Ydiv = cat_vol_div( Ycis , opt.resV , 1 ); dth = [-0.5 0.8];
  Ydiv(Ydiv<dth(1)) = min(0,2*dth(1) - Ydiv(Ydiv<dth(1)));
  Ydiv(Ydiv>dth(2)) = max(0,2*dth(2) - Ydiv(Ydiv>dth(2)));
  %%
  if opt.cb
    %% Sharpening based on bias correction and the divergence map:
    if debug, tic; Ycis = Ymf;  end
    % first bias correction (medium frequency)
    YMR = (Ycis .* ( -Ydiv+1).^3)>2; Ymw = Ymf .* YMR; 
    for i=1:round(2/opt.resV), Ymw = cat_vol_localstat(Ymw,YMR,1,3); end; 
    Ymw(Ymf<=1) = 3; 
    if ~debug, clear YMR; end
    Ymw = cat_vol_approx(Ymw/3,'nh',opt.resV,2); 
    Ycis = Ycis ./ Ymw;     
    if ~debug, clear Ymw; end
    % sharpening % high frequency)
    Ycis = Ycis .* ( -Ydiv ./ max(0.5,abs( min(3,max(1,Ycis))/1.5-0.75)) /2 + 1 ).^2; 
    Ycis = max(1,min(3,Ycis));
  else
    % sharpening % high frequency)
    Ycis = Ycis .* ( -Ydiv ./ max(0.5,abs( min(3,max(1,Ycis))/1.5-0.75)) /8 + 1 ).^2; 
    Ycis = max(1,min(3,Ycis));
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
    
  if debug, tic; end
  YMR = (Ycis-Ydiv*5)>2.3 & Ycis>2.3 & Ycis<3.3 & Ydiv<0.05; YMR(smooth3(YMR)<0.45)=0;
  YMR = cat_vol_morph(YMR,'ldo',((1-opt.cb)/2)/opt.resV);
  YMR = cat_vol_morph(YMR,'ldc',0.5/opt.resV);
  Ymw = Ycis .* YMR; 
  for i=1:round(3/opt.resV), Ymw = cat_vol_localstat(Ymw,YMR,1,3); end
  for i=1:round(2/opt.resV), Ymw = cat_vol_localstat(Ymw,YMR,1,1); end
  if ~debug, clear YMR; end
  Ymw(cat_vol_morph(Ycis./Ymw < 2.5/3 , 'ldc',0.5/opt.resV)) = 0;  
  Ymw = cat_vol_approx(Ymw,'nh',opt.resV,2); 
  Ycisc = Ycis ./ Ymw * 3; 
  if ~debug, clear Ymw; end  
  Ymc = Ycisc .* (Ycis>1.1 & Ycisc<2.9);
  for i=1:round(3/opt.resV), Ymc = cat_vol_localstat(Ymc,Ycis>1.1 & Ycisc<2.9,1,2); end
  for i=1:round(2/opt.resV), Ymc = cat_vol_localstat(Ymc,Ycis>1.1 & Ycisc<2.9,1,1); end
  Ymc(Ymc>2.4 | Ymc==0) = 0;  
  
  
  %% Estimate wider maps of the major tissues to create the mid surface: 
  %  WM
  Ywm = (Ycisc>2.4) | (Ycisc>2.3 & Ymc==0);
  Ywm(smooth3(Ywm)<0.3) = 0; 
  Ywm = cat_vol_morph(Ywm,'l');
  Ywm = ~cat_vol_morph(~Ywm,'l'); 
  %  CSF
  Ycm = (Ycisc-min(0,Ydiv*2))<2.1 & (Ycis-min(0,Ydiv*2.1))<2 & Ymc<2.2; if ~debug, clear Ymc; end  
  Ycm(smooth3(Ycm)>0.6) = 1; Ycm(smooth3(Ycm)<0.4) = 0; 
  Ycm = ~cat_vol_morph(~Ycm,'lo',0.5/opt.resV) & ~Ywm & smooth3(Ycis)<2; 
  Ywd = cat_vbdist(single(Ywm),Ycm<1); Ywd(Ywd>100)=0;
  YMR = Ywd > cat_stat_nanmedian(Ywd(Ywd(:)>0)) + 2*cat_stat_nanstd(Ywd(Ywd(:)>0)); 
  Ycm(YMR) = 1; Ywm(YMR) = 0; Ywd(YMR) = 0; clear YMR; 
  if debug, toc; else clear Ycisc; end
 
  % Speed maps for asymmetric mapping and distance estimation 
  if debug, tic; end
  if opt.fast
    Ycd  = cat_vbdist(single(Ycm),Ywm<1); 
  else
    Ycm  = single(Ycm); Ycm(Ywm>0) = nan; 
    Ywm  = single(Ywm); Ywm(Ycm>0) = nan;
    F    = max(eps,max(eps,smooth3(min(1,(Ycis-1).^2)))); 
    Ywd  = cat_vol_eidist(Ywm,F); Ywd(isnan(Ywd))=0; %if ~debug, clear Ywm; end
    F    = smooth3(max(eps,smooth3(max(eps,min(1,(3-Ycis)).^2))));
    Ycd  = cat_vol_eidist(Ycm,F); Ycd(isnan(Ycd))=0; %if ~debug, clear Ycm; end
  end
  if debug, toc; end
 
  
  % PBT mapping 
  %  ~30s
  if debug, tic; end
  Ywdc = Ywd+0; Ycdc = Ycd+0; if ~debug, clear Ycd; end % we need doubles (copies) for pbtv!  
  [Ygt1,Yv1,Yv2] = cat_vol_pbtv( Ycis , Ywdc , Ycdc ); clear Ywdc Ycdc; 
  if debug, toc; end
  
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
  Yppv = Yv1 ./ max(eps,Yv1 + Yv2); Yppv(Ywm>0)=1; Yppv(Ycm>0)=0; 
  Yppd = (Ygt1 - Ywd) ./ max(eps,Ygt1); Yppd(Ywm>0)=1; Yppd(Ycm>0)=0; 
  Yppv = cat_vol_median3(Yppv,Yppv>0.1 & Yppv<0.9,true(size(Yppv)),0.1); 
  Yppd = cat_vol_median3(Yppd,Yppd>0.1 & Yppd<0.9,true(size(Yppd)),0.1); 
  Ypp  = Yppv*0.5 + 0.5*Yppd; if ~debug, clear Yppv Yppd Ywd; end
  Ypp(cat_vol_morph(Ypp>=0.5,'ldo',max(0.5,(1.5 - opt.cb)/(2*opt.resV)))==0 & Ypp>=0.5) = 0.49;
  Ypp(Ypp<0.5 & ~cat_vol_morph(Ypp<0.5,'l')) = 0.51;
  if ~debug, clear Yv1 Yv2; else toc; end
  
  if ~debug
    clear Ygt1 Ygt; 
  elseif 0
    %% just display - use sd to rotate the surface
    [D,I] = cat_vbdist(single(Ygt1>0),cat_vol_morph(Ypp<0.5,'d') & cat_vol_morph(Ypp>0.5,'d')); Ygt1 = Ygt1(I); clear D I; 
    sd = 0; cat_surf_render2('Disp',isosurface(shiftdim(YM),0.5,shiftdim(cat_stat_histth(Ymf,0.99),sd)));
  end
  if debug, toc; end
  
  
  %% distance from the central intensity surface (CIS)
  %  ----------------------------------------------------------------------
  %  Distance from the CIS to the CSF and WM. For the CSF case the removal 
  %  of meninges and blood-vessels is maybe important and depends on the 
  %  distance estimation. vbdist require removal, whereas eidist may works 
  %  better without and a good speed map F.
  stime = cat_io_cmd('    GM distance estimation: ','g5','',opt.verb,stime);if debug, tic; end
  
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
  F = smooth3(max(eps,max(eps,min(1,max(0.5,Ymf-0.5)/1.55))));             % speed map for eikonal distance
  [Ycid,t,Yd] = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug); clear t;     %#ok<ASGLU>
  YM(Yd>1.4*Ycid & Ycid>1) = nan; Ycid(Yd>1.4*Ycid & Ycid>1) = nan; clear Yd;   % correct for meninges
  if opt.enlarge
    % correc distance values behind the boundary
    YM2   = max(YM,max(0,min(1,Ymf-1))); 
    Ycidc = cat_vol_eidist(YM2,F,[1 1 1],1,1,0,opt.debug); 
    Ycid  = Ycid - Ycidc; clear Ycidc;                            
  end  
  
  
  %  CIS distance to the CSF and correction:
  %  As far as we can use the YM boundary description this part is much shorter. 
  F   = smooth3(max(eps,max(eps,min(1,(3-Ymf)).^2)));                      % speed map for eikonal distance
  YMM = cat_vol_morph(Ymf>2.5,'e',opt.enlarge) | isnan(Ymf);               % change YMM to have the deep WM
  YMI = YM; YMI(isnan(YMI))=0; YMI = 1 - YMI; YMI(YMM & YM>0.9) = nan;     % create inverse boundary map
  [Ycid2,t,Yd] = cat_vol_eidist(YMI,F,[1 1 1],1,1,0,opt.debug); clear t YMM;   %#ok<ASGLU>
  Ycid2(Yd>1.4*Ycid2 & Ycid2>1) = nan; clear Yd;
  if opt.enlarge
    % correct distance values behind the boundary
    YM2   = max(YM,max(0,min(1,3-Ymf))); 
    Ycidc = cat_vol_eidist(YM2,F,[1 1 1],1,1,0,opt.debug); 
    Ycid2  = Ycid2 - Ycidc; clear Ycidc;                            
  end
  clear F;
   
  % combine both distance maps
  Ycid = Ycid - Ycid2; YnGM = isnan(Ycid); Ycid(isnan(Ycid)) = 0; clear Ycid2;
  Ycid = cat_vol_median3(Ycid,Ycid>0,Ycid>0,0.1); 
  if debug, toc; end

  
  %% CSF and WM distance
  stime = cat_io_cmd('    CSF and WM distance estimation: ','g5','',opt.verb,stime); if debug, tic; end
  F     = smooth3(max(eps,max(eps,min(1,(3-Ymf)).^2)));                    % speed map for eikonal distance
  Ycid2 = Ycid; Ycid2(isnan(Ycid))=0; Ycid2 = max(0,min(4,Ycid2 +1 - cat_stat_nanmean(Ycid2(Ycid2(:)>0 & Ymf(:)>1.3 & Ymf(:)<1.6)))); 
  YM3   = max( max(0,min(1,max(2-Ymf,Ycid2/4))) , smooth3(isnan(YM)) ); YM3(YM3==0 & YnGM) = nan;
  Ycsfd = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); 
  if opt.enlarge
    % correc distance values behind the boundary
    YM3    = max(YM3,max(0,min(1,3-Ymf))); 
    YM3(cat_vol_morph(YM3>0.1,'e')) = 1; YM3(cat_vol_morph(YM3<0.9,'e')) = 0; 
    YM3(YM3==0 & YnGM) = nan;
    Ycsfdc = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); 
    Ycsfd  = Ycsfd - Ycsfdc; clear Ycsfdc;                          
  end
  % WM distance 
  F     = max(eps,max(eps,min(1,(Ymf-1).^1.5))); 
  YM3   = max(max(0,min(1,(Ymf-2))) ,  smooth3(isnan(YMI)) ); YM3(YM3==0 & YnGM) = nan;
  Ywmd  = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); 
  if opt.enlarge
    % correct distance values behind the boundary
    YM3   = max(YM3,max(0,min(1,Ymf-1))); YM3(YM3==0 & YnGM) = nan;
    Ywmdc = cat_vol_eidist(YM3,F,[1 1 1],1,1,0,opt.debug); 
    Ywmd  = Ywmd - Ywmdc; clear Ywmdc;                             
  end
  if debug, toc; else clear YM3 F; end
  
  
  %% projection-based thickness estimation 
  stime = cat_io_cmd('    Thickness projection: ','g5','',opt.verb,stime); if debug, tic; end
 % Ycidx = Ycid+1000; Ycidx(isnan(Ycidx))=0; Ycidx = cat_vol_localstat(Ycidx,Ycidx>0,1,1); Ycidx(Ycidx>0)=Ycidx(Ycidx>0)-1000;
  Ygmt1 = cat_vol_pbtv( 4 - Ymf , -Ycid , Ywmd  ); 
 % Ycidx = Ycid+1000; Ycidx(isnan(Ycidx))=0; Ycidx = cat_vol_localstat(Ycidx,Ycidx>0,1,1); Ycidx(Ycidx>0)=Ycidx(Ycidx>0)-1000; 
  Ygmt2 = cat_vol_pbtv(     Ymf ,  Ycid , Ycsfd );
  Ygmt1(Ygmt1<0 | Ygmt2<0) = 0; Ygmt2(Ygmt1<0 | Ygmt2<0) = 0; % remove enlarged areas
  Ygmt  = (Ygmt2 + Ygmt1); Ygmt(YnGM)=0; 
  
  %% Create percentage position map Ypp
  Ypp = max(0, min(1 , (Ygmt - Ygmt1 - Ycid.*(Ygmt>0)) ./ max(eps,Ygmt) ));  Ypp(Ymf>=2.5 & Ywmd<0.5) = 1; 
  Ypp = cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & Ypp>0.1,true(size(Ypp)),0.1);
  Ypp = cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & Ypp>0.1,true(size(Ypp)),0.1);
  Ypp(Ypp>0.5 & ~cat_vol_morph(Ypp>0.5,'ldo',0.5/opt.resV)) = 0.49;
  Ypp(Ypp<0.5 & ~cat_vol_morph(Ypp<0.5,'l',0.25/opt.resV)) = 0.51;
  if debug, toc; end
  
  %% final corrections
  Ygmt = Ygmt * opt.resV; 
  Ygmt(Ygmt>10) = 10; 
  if debug, toc; end
  cat_io_cmd(' ','g5','',opt.verb,stime);
  if opt.debug, cat_io_cmd(' ','','',opt.debug,stime2); end

  if 0
    %% just display - use sd to rotate the surface
    [D,I] = cat_vbdist(single(Ygmt>0),cat_vol_morph(Ypp<0.5,'d') & cat_vol_morph(Ypp>0.5,'d')); Ygmt = Ygmt(I); clear D I; 
    sd = 1; cat_surf_render2('Disp',isosurface(shiftdim(Ypp,sd),0.5,shiftdim(Ygmt,sd)));
  end
  
  
end
