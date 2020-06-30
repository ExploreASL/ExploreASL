function [Ygmt,Ypp,Ymf,Ywmd,Ywmdc] = cat_vol_pbt2(Ymf,opt)
% ______________________________________________________________________
%
% Cortical thickness and surface position estimation. 
% 
%   [Ygmt,Ypp,Ywmd,Ycsfd] = cat_vol_pbt(Ymf,opt)
%  
%   Ygmt:      GM thickness map 
%   Ypp:       percentage position map
%   Ywmd:      WM distance map
%   Ycsfd:     CSF distance map
% 
%   Ymf:       tissue segment image or better the noise, bias, and 
%              intensity corrected 
%
%   opt.resV   voxel resolution (only isotropic)
%   opt.method choose of method {'pbt2x','pbt2'} with default=pbt2x as 
%              the method that is described in the paper.
%   opt.pbtlas GM intensity correction to reduce myelination effects
%              (added 201908)
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
%
% ______________________________________________________________________
% $Id: cat_vol_pbt2.m 1552 2020-01-17 10:19:24Z dahnke $ 


% default variables and check/set function  
  if ~exist('opt','var'), opt=struct(); end

  def.resV      = 1;
  def.dmethod   = 'eidist';
  def.method    = 'pbt2x';  % pbt is worse ... just for tests!
  def.debug     = cat_get_defaults('extopts.verb')>2;
  def.verb      = cat_get_defaults('extopts.verb')-1;
  def.pbtlas    = 1; % no/yes
  opt           = cat_io_checkinopt(opt,def);
  opt.resV      = mean(opt.resV);
  
  minfdist = 2; 
  
  if 0
    % RD 201803: 
    % Remove blood vessels & meninges
    % This block is too aggressive and removes gyral peaks (thickness overestimation).
    % However, some correction is required ...
    Ymx = cat_vol_morph(Ymf>2.5,'l');  Ymf(Ymf>2.5 & ~Ymx)=2;   clear Ymx
    Ymx = cat_vol_morph(Ymf>2.2,'lo'); Ymf(Ymf>2.2 & ~Ymx)=2.1; clear Ymx
    Ymx = cat_vol_morph(Ymf>1.5,'lo'); Ymf(Ymf>1.5 & ~Ymx)=1.2; clear Ymx
  end


  %% Distance maps
  if (sum(round(Ymf(:))==Ymf(:)) / numel(Ymf))>0.9, bin=1; else bin=0; end
  
  %  WM distance 
  %  Estimate WM distance Ywmd and the outer CSF distance Ycsfdc to correct
  %  the values in CSF area are to limit the Ywmd to the maximum value that 
  %  is possible within the cortex.
  %  The increasement of this area allow a more accurate and robust projection. 
  %  cat_vol_eidist used speed map to align voxel to the closer gyrus
  %  that is not required for the correction map.
  %
  %  RD 201803:
  %  The speed map weighting "max(0.5,min(1,Ymf/2))" is not strong enough to  
  %  support asymmetric structures. The map "max(eps,min(1,((Ymf-1)/1.1).^4))"  
  %  works much better but it leads to much higher thickness results (eg. in 
  %  the Insula).
  
  newspeedmapF = 1; 


  if opt.verb, fprintf('\n'); end; stime2=clock;
  YMM = cat_vol_morph(Ymf<1.5,'e',1) | isnan(Ymf);
  switch opt.dmethod
    case 'eidist' 
      % [D,I] = vbm_vol_eidist(B,L,[vx_vol,euclid,csf,setnan,verb])
      
      if newspeedmapF
        F = max(eps,min(1,((Ymf-1)/1.1).^4)); 
      else 
        F = max(0.5,min(1,Ymf/2)); % R1218
      end
      
      
      if opt.pbtlas == 0
      % Simple pure intensity based model that show strong variation in 
      % myelinated regions.
        YM  = max(0,min(1,(Ymf-2))); YM(YMM) = nan; 
        stime = cat_io_cmd('    WM distance: ','g5','',opt.verb); 
        
      else
      % New local intensity optimization to remove myelination effects
      % that works in a similar fasion as pbtlas but only focuses on the 
      % GM/WM boundary.
      % This model was implemented after the results of Shahrzad Kharabian
      % Masouleh (Juelich) that showed strong underestimations and increased
      % variance in CAT compared to FS and CIVET.
      % RD 201908
      
        stime = cat_io_cmd('    WM boundar optimization: ','g5','',opt.verb);
        vx_vol_org = 1; %max(1,sqrt(sum(opt.vmat(1:3,1:3).^2))));
        
        
        %% WM intensity:
        %  The result should be a map with an average value of 1 and no
        %  anatomical structures (maybe some differencence for regional
        %  folding independent pattern).
        %  Important to avoid gyral overestimations, however too small
        %  filter size / low number of iterations scan also lead to errors.
        [Ymfr,resT2] = cat_vol_resize(Ymf,'reduceV',opt.resV,mean(vx_vol_org),32,'meanm'); 
        YM2max = max(0,min(1,Ymfr-2)); YM2max(YM2max<0.0) = 0; % lower threshold > higher filter range/more iterations
        YM2max = YM2max .* cat_vol_morph(YM2max>0 ,'l',[10 0.1 ]);   
        for ix=1:3 %round(2/mean(resT2.vx_volr)) % iterations are faster/smoother than single calls with higher values
          YM2max = cat_vol_localstat(YM2max,YM2max>0.0,min(10,2.0/mean(resT2.vx_volr)),3); % higher treshold and maximum filter 
          YM2max = cat_vol_localstat(YM2max,YM2max>0.0,min(10,1.0/mean(resT2.vx_volr)),1); % mean filter for smooth data
        end
        YM2max  = cat_vol_approx(YM2max,'nh',resT2.vx_volr,2); % use (slow) high resolutions 
        YM2max  = cat_vol_resize(YM2max,'dereduceV',resT2);

 %  stime = cat_io_cmd(sprintf('    WM distance (GM - %0.2f): ',mean(vx_vol_org)),'g5','',opt.verb,stime);      
        %% GM intensity:
        %  find the higher intesity GM average and avoid the CSF PVE 
        lth = 1.8; %1.8; % 1.8 - 2.2
        [Ymfr,YM2maxr,resT2] = cat_vol_resize({Ymf,YM2max},'reduceV',opt.resV,mean(vx_vol_org),32,'meanm');
        
        % This is a simple tissue thickness model from a more central GM layer 
        % (currently 2.1, possible range 1.5 to 2.2) to avoid corrections 
        %% in thicker areas (mostly gyris). 
        Ygd   = cat_vbdist(2.1 - Ymfr); % estimate distance map to central/WM surface, lower thresholds are also possible (range 1.5 to 2.2, default 2.1) 
        Ygdt  = cat_vol_pbtp2(max(2,min(3,4-Ymfr)),Ygd,inf(size(Ygd),'single')) * mean(resT2.vx_volr);
        Ygdt  = cat_vol_median3(Ygdt,Ygdt>0.01,Ygdt>0.01);                    
        Ygdt  = cat_vol_localstat(Ygdt,Ygdt>0.1,1/mean(resT2.vx_volr),1);    

        Ywd   = cat_vbdist(2.5 - Ymfr);  
        Ywdt  = cat_vol_pbtp2(max(2,min(3,4-Ymfr)),Ywd,inf(size(Ywd),'single')) * mean(resT2.vx_volr);
        Ywdt  = cat_vol_median3(Ywdt,Ywdt>0.01,Ywdt>0.01);                    
        Ywdt  = cat_vol_localstat(Ywdt,Ywdt>0.1,1/mean(resT2.vx_volr),1);    

        
        %%
        % estimate local GM average intensity 
        YM2min   = max(0,min(1.5,Ymfr - lth)) .* ...
                      (Ymfr<( (YM2maxr+2) .* (2.5 + min(0.45,Ygdt/8) )) /3) .* ...
                      (Ymfr<( (YM2maxr+2) .* (2.5 + min(0.45,Ywdt/8) )) /3) .* ...
                      (Ywd<4*mean(resT2.vx_volr)) .*...
                      (cat_vbdist(single(Ymfr./YM2maxr<2))<3.5/mean(resT2.vx_volr)); %clear YM2maxr Ywdt 
                    
        % remove small dots
        YM2min = YM2min .* cat_vol_morph(YM2min>0 ,'l',[100 0.01]);   
        YM2min(smooth3(YM2min>0)<0.4) = 0;
        YM2min = cat_vol_median3(YM2min,YM2min>0,YM2min>0); 
        % avoid PVE
        YM2min2 = cat_vol_localstat(YM2min,YM2min>0,1,2);
        YM = Ymfr>2 & YM2min==0; YM = cat_vol_morph(YM,'dd',1) & ~YM; YM2min( YM ) = YM2min2( YM );
        YM2min2 = cat_vol_localstat(YM2min,YM2min>0,1,3);
        YM = Ymfr<2 & YM2min==0; YM = cat_vol_morph(YM,'dd',1) & ~YM; YM2min( YM ) = YM2min2( YM ); 
        %clear Ymfr YM

        % remove further atypical values
        YM2mina = YM2min; 
        for ix=1:10
          YM2mina = cat_vol_localstat(YM2mina,YM2mina>0,2/mean(resT2.vx_volr),1);
        end
        YM2min  = YM2min .* ((YM2mina - YM2min)>-0.2 & (YM2mina - YM2min)<0.3); clear YM2mina;
        
        % remove PVE 
        YM2min = cat_vol_median3(YM2min,YM2min>0); 
        
        % further local filtering
        for ix=1:10 %round(5/mean(resT2.vx_volr)) 
          YM2min = cat_vol_localstat(YM2min,YM2min>0,1/mean(resT2.vx_volr),1);
        end
        YM2min(Ywd>6*mean(resT2.vx_volr)) = 0.01; 
        
        YM2min  = cat_vol_approx(YM2min,'nh',resT2.vx_volr,2); % use (slow) high resolutions
        YM2min  = cat_vol_resize(YM2min,'dereduceV',resT2);
        Ywd     = cat_vol_resize(Ywd,'dereduceV',resT2);
        YM2min  = YM2min - (2 - lth);
 %stime = cat_io_cmd(sprintf('    WM distance (GM - %0.2f): ',mean(vx_vol_org)),'g5','',opt.verb,stime);              
        %% create intensity optimized boundary image
        YM    = max(0,min(1,((Ymf - (YM2min + 2)) ./ (YM2max - YM2min)  ))); 
        %% use only fat areas
        YMo  = max(0,min(1,(Ymf-2)));
        YMc  = YMo - YM;
        YMc(YMc>0 & YMc<0.15) = YMc(YMc>0 & YMc<0.15);
        YMc(YMc>0 & YMc<0.10) = YMc(YMc>0 & YMc<0.10)/2;
        YMc(YMc>0 & YMc<0.05) = YMc(YMc>0 & YMc<0.05)/4;
        YMc(YMc<0.05) = 0;
        YM   = YMo  - YMc; 
        clear YMo YMc YMcs;
      
        % distance corrections to stabilize thin structures
        YM    = max(YM,min(1,Ywd/mean(resT2.vx_volr)/2 - 2)); % L4 distance based corrections
        YMd   = 0.75 - YM; YMd(YM>0.9) = nan; 
        YM    = max(YM,min(1,cat_vol_eidist(YMd,ones(size(YM),'single'),[1 1 1],1,1,0,opt.debug) - 0.5));  
   %stime = cat_io_cmd(sprintf('    WM distance (GM - %0.2f): ',mean(vx_vol_org)),'g5','',opt.verb,stime);            
        %%
        for ix=1:3, YM(YM>0.5 & smooth3(YM>0.5)<0.2) = 0.45; end
        YM(cat_vol_morph(YM>0.9 ,'lc')) = 1; 
        YM(YMM) = nan;
        clear YM2max YM2min
        
        stime = cat_io_cmd('    WM distance: ','g5','',opt.verb,stime); 
      end
      
      
      %% final distance estimation
      YMwm   = YM;
      Ywmd   = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug);
      
      if newspeedmapF
        F = max(eps,min(1,((Ymf-1)/1.1).^4)); 
      else 
        F = max(1.0,min(1,Ymf/2)); % R1218 - no speed difference!
      end
      YM  = max(0,min(1,(Ymf-1))); YM(YMM) = nan; Ywmdc = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug); 
      
      clear F;
      
    case 'vbdist'
      stime = cat_io_cmd('    WM distance: ','g5','',opt.verb); stime2=stime;
      YM  = max(0,min(1,(Ymf-2))); Ywmd   = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
      YM  = max(0,min(1,(Ymf-1))); Ywmdc = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
  end
  clear YMM; 
  if ~bin
    % limit the distance values outside the GM/CSF boudary to the distance possible in the GM
    YM  = Ywmd>minfdist & Ymf<=1.5; Ywmd(YM) = Ywmd(YM) - Ywmdc(YM); Ywmd(isinf(Ywmd)) = 0; clear Ycsfdc;
    % smoothing of distance values inside the GM
    %YM  = Ywmd>minfdist & Ymf> 1.5; YwmdM = Ywmd; YwmdM = cat_vol_localstat(YwmdM,YM,1,1); Ywmd(YM) = YwmdM(YM);
    % smoothing of distance values outside the GM
    YM  = Ywmd>minfdist & Ymf<=1.5; YwmdM = Ywmd; for i=1:1, YwmdM = cat_vol_localstat(YwmdM,YM,1,1); end; Ywmd(YM) = YwmdM(YM);
    % reducing outliers in the GM/CSF area
    YM  = Ywmd>minfdist & Ymf< 2.0; YwmdM = Ywmd; YwmdM = cat_vol_median3(YwmdM,YM,YM); Ywmd(YM) = YwmdM(YM); clear YwmdM YM;
  end
  
  minfdist = 2/opt.resV; 
  %  CSF distance
  %  Similar to the WM distance, but keep in mind that this map is
  %  incorrect in blurred sulci that is handled by PBT
  stime = cat_io_cmd('    CSF distance: ','g5','',opt.verb,stime);
  if exist('YMwm','var')
    YMM = cat_vol_morph(Ymf<1.5,'e',1) | cat_vol_morph(YMwm>0.5,'e',1) | isnan(Ymf); % this was dilate???
  else
    YMM = cat_vol_morph(Ymf<1.5,'e',1) | cat_vol_morph(Ymf>2.5,'e',1) | isnan(Ymf); % this was dilate???
  end
  switch opt.dmethod
    case 'eidist'
      if newspeedmapF
        F = max(eps,min(1,(4-Ymf)/2).^2);
      else
        F = max(0.5,min(1,(4-Ymf)/2)); % R1218
      end
      YM  = max(0,min(1,(2-Ymf)));   YM(YMM) = nan; Ycsfd = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug); 
      if newspeedmapF
        F = max(eps,min(1,(4-Ymf)/2).^2);
      else
        F = max(1,min(1,(4-Ymf)/2)); % R1218
      end
      if exist('YMwm','var')
        YM  = 1 - YMwm; YM(YMM) = nan; Ywmdc = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug);  Ywmdx = Ywmdc;
      else
        YM  = max(0,min(1,(3-Ymf)));   YM(YMM) = nan; Ywmdc = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug);
        YM  = max(0,min(1,(2.7-Ymf))); YM(YMM) = nan; Ywmdx = cat_vol_eidist(YM,F,[1 1 1],1,1,0,opt.debug)+0.3;
      end
      clear F;
    case 'vbdist'
      YM  = max(0,min(1,(2-Ymf)));   Ycsfd = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
      YM  = max(0,min(1,(3-Ymf)));   Ywmdc = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
      YM  = max(0,min(1,(2.7-Ymf))); Ywmdx = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.2); 
  end
  Ywmdc = min(Ywmdc,Ywmdx);
  clear YMM Ywmdx;
  if ~bin
    if exist('YMwm','var')
      YM = Ycsfd>minfdist & YMwm>=0.5; Ycsfd(YM) = Ycsfd(YM) - Ywmdc(YM); Ycsfd(isinf(-Ycsfd)) = 0; clear Ywmdc;
      %YM = Ycsfd>minfdist & YMwm< 0.5; YcsfdM = Ycsfd; YcsfdM = cat_vol_localstat(YcsfdM,YM,1,1); Ycsfd(YM) = YcsfdM(YM);
      YM = Ycsfd>minfdist & YMwm>=0.5; YcsfdM = Ycsfd; for i=1:1, YcsfdM = cat_vol_localstat(YcsfdM,YM,1,1); end; Ycsfd(YM) = YcsfdM(YM);
    else
      YM = Ycsfd>minfdist & Ymf>=2.5; Ycsfd(YM) = Ycsfd(YM) - Ywmdc(YM); Ycsfd(isinf(-Ycsfd)) = 0; clear Ywmdc;
      %YM = Ycsfd>minfdist & Ymf< 2.5; YcsfdM = Ycsfd; YcsfdM = cat_vol_localstat(YcsfdM,YM,1,1); Ycsfd(YM) = YcsfdM(YM);
      YM = Ycsfd>minfdist & Ymf>=2.5; YcsfdM = Ycsfd; for i=1:1, YcsfdM = cat_vol_localstat(YcsfdM,YM,1,1); end; Ycsfd(YM) = YcsfdM(YM);
    end
    YM = Ycsfd>minfdist & Ymf> 2.0; YcsfdM = Ycsfd;  YcsfdM = cat_vol_median3(YcsfdM,YM,YM); Ycsfd(YM) = YcsfdM(YM); clear YcsfdM YM;
  end  


  %% PBT thickness mapping 
  %  PBT is the default thickness estimation, but PBT2x is the optimized
  %  version that usex both sulci and gyri refinements, because not only 
  %  thin sulci can be blurred. PBT2x is also the method that is
  %  described in the paper.
  %  PBTv is new version that uses the volume rather than the distance. 
  %  Although this works in principle, this is biased by interpolation 
  %  artifacts and not-optimal WMD mapping.  
  iter = 0;%round(0.5/mean(opt.resV));
  fs   = 0.5;
  if strcmp(opt.method,'pbtv')  
    Ywmdo = Ywmd+0; 
    
    %% estimate cortical thickness and map the local volumes
    stime = cat_io_cmd('    PBTV thickness: ','g5','',opt.verb,stime);
    [Ygmt,Yv1,Yv2,Ypp] = cat_vol_pbtv(Ymf,Ywmd,Ycsfd);
    Ygmts = Ygmt; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,(Ygmt>1 | Ypp>0.1) & Ygmt>0 & (Ygmt>1 | Ymf>1.8),1,1); end; Ygmt(Ygmts>0) = Ygmts(Ygmts>0); 
  
    % volume-based PP map
    Ypp = ((Ygmt - Ywmd)./(Ygmt + eps) + (Ymf>2.5))*0.5 + 0.5*min( 1, Yv1./(Yv1 + Yv2 + eps) + (Ymf>2.5)); 
    
    
  elseif strcmp(opt.method,'pbt2x') ||  strcmp(opt.method,'pbt2x2')  
    % Estimation of the cortical thickness with sulcus (Ygmt1) and gyri 
    % correction (Ygmt2) to create the final thickness as the minimum map
    % of both.
    stime = cat_io_cmd('    PBT2x thickness: ','g5','',opt.verb,stime);
    
    % estimate thickness with PBT approach
    if opt.pbtlas, Ymfo=Ymf; Ymf = single(1 + 2*((Ywmd>0 & Ycsfd>0) | Ymfo>2)) - (Ywmd>0 & Ycsfd>0);  end
    Ygmt1 = cat_vol_pbtp2((Ymf),Ywmd,Ycsfd);  
    Ygmt2 = cat_vol_pbtp2((4-Ymf),Ycsfd,Ywmd);
    
    % Error handling 
    % For some unkown reasons the sulcus reconstruction of cat_vol_pbtp failed in some cases (not directly reproducable).      
    % Reprocessing is solving this problem, but further investigation of cat_vol_pbtp.cpp would be good (RD 20190811).
    % Maybe it depends on the initialization of the regions, e.g., using Ymf without rounding and incorrect boundary seams to increase the problems.  
    % Now use cat_vol_pbtp2.cpp (RD20200111).
    mask   = @(Y) Y(:)>0 & Y(:)<1000000; 
    rerun = 0; rerunlim = 3; 
    while rerun <= rerunlim && isnan( mean( Ygmt1(mask(Ygmt1))) ) || mean( Ygmt1(mask(Ygmt1)))>100
      Ygmt1 = cat_vol_pbtp2((Ymf),Ywmd,Ycsfd);  
      rerun = rerun + 1; 
      pause(rand*3);
    end
    if rerun == rerunlim && isnan( mean( Ygmt1(mask(Ygmt1))) ) || mean( Ygmt1(mask(Ygmt1)))>100
      error('cat_vol_pbtp2:bad_mapping1','Untypcial values in PBT thickness mapping detected. ');
    end
    rerun = 0; 
    while rerun <= rerunlim && isnan( mean( Ygmt2(mask(Ygmt2))) ) || mean( Ygmt2(mask(Ygmt2)))>100
      Ygmt2 = cat_vol_pbtp2((4-Ymf),Ycsfd,Ywmd);
      rerun = rerun + 1; 
      pause(rand*3);
    end
    if rerun == rerunlim && rerunlim && isnan( mean( Ygmt2(mask(Ygmt2))) ) || mean( Ygmt2(mask(Ygmt2)))>100
      error('cat_vol_pbtp2:bad_mapping2','Untypcial values in PBT thickness mapping detected. ');
    end
    
    %% avoid meninges !
    update_WMD =  strcmp(opt.method,'pbt2x2'); 
  
    Ygmt1 = min(Ygmt1,Ycsfd+Ywmd);
    if update_WMD % eigentlich muss das hier an sein ...
      Ygmt2 = min(Ygmt2 + (Ygmt2>0) .* 2,Ycsfd+Ywmd);  
    else
      Ygmt2 = min(Ygmt2,Ycsfd+Ywmd);  
    end
    
    % median filter to remove extreme outliers 
    YM    = Ymf>1.5 & Ymf<2.5; 
    Ygmt1 = cat_vol_median3(Ygmt1,YM,Ygmt1>0,0.2);
    Ygmt2 = cat_vol_median3(Ygmt2,YM,Ygmt2>0,0.2);
    
    %% estimation of Ypp for further GM filtering without sulcul blurring
    if update_WMD
      [Ygmt,Yi] = min(cat(4,Ygmt1,Ygmt2),[],4);
      Ywmdc   = (Ywmd.*(Yi==1) + (Yi==2).*(Ygmt2 - Ycsfd)); 
      Ywmdc   = cat_vol_median3(Ywmdc,Yi==2,Yi>0,0.05);
      Ywmdcs  = cat_vol_localstat(Ywmdc,Yi==2,2,1); Ywmdc(Ywmdcs>0) = Ywmdc(Ywmdcs>0); clear Ycsfdcs
      Ywmd    = Ywmdc; 
    else
      Ygmt    = min(cat(4,Ygmt1,Ygmt2),[],4);
    end
    
    Ygmt    = Ygmt*(1-fs) + fs*cat_vol_median3(Ygmt,Ygmt>0,Ygmt>0,0.05);
    Ypp     = zeros(size(Ymf),'single'); Ypp(Ymf>=2.5) = 1;  YM      = Ygmt>0; 
    Ypp(YM) = min(Ycsfd(YM),Ygmt(YM) - Ywmd(YM)) ./ (Ygmt(YM) + eps); Ypp(Ypp>2) = 0;
    Ygmts   = cat_vol_localstat(Ygmt,Ygmt>0 & Ypp>0.1,1,1); Ygmt(Ygmts>0) = Ygmt(Ygmts>0)*(1-fs) + fs*Ygmts(Ygmts>0); clear Ygmtps
    %%
    Ypp(Ypp==0 & Ymf>=2.5 & Ywmd<=0) = 1; 
    Ypp(Ypp==0 & Ymf>=2.5 & Ywmd<1)  = 0.99; 
    Ypp(Ypp==0 & Ymf>=2.0 & Ywmd<2)  = 0.985; 
    Ypp(Ypp==0 & Ymf>=2.5)           = 0.95; 

    Ypp = Ypp*(1-fs) + fs*cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & (Ypp>0.1 | Ypp==0),true(size(Ypp)),0.1);
    Ypp = Ypp*(1-fs) + fs*cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & (Ypp>0.1 | Ypp==0),true(size(Ypp)),0.1);

    Ypp(Ypp>0.5 & ~cat_vol_morph(Ypp>0.5,'ldo',0.50/opt.resV)) = 0.49;
    Ypp(Ypp<0.5 & ~cat_vol_morph(Ypp<0.5,'l'  ,0.25/opt.resV)) = 0.51;

 
  else
    % Estimation of thickness map Ygmt and percentual position map Ypp.
    stime = cat_io_cmd('    PBT2 thickness: ','g5','',opt.verb,stime);
    
    % estimate thickness with PBT approach
    if opt.pbtlas, Ymfo=Ymf; Ymf = single(1 + 2*((Ywmd>0 & Ycsfd>0) | Ymfo>2)) - (Ywmd>0 & Ycsfd>0);  end
    
    % estimate thickness with PBT approach
    if opt.pbtlas, Ymfo=Ymf; Ymf = single(1 + 2*((Ywmd>0 & Ycsfd>0) | Ymfo>2)) - (Ywmd>0 & Ycsfd>0);  end
    Ygmt    = cat_vol_pbtp2( (Ymf) ,Ywmd,Ycsfd);   
    Ygmt    = Ygmt*(1-fs) + fs*cat_vol_median3(Ygmt,Ygmt>0,Ygmt>0,0.05);
    Ypp     = zeros(size(Ymf),'single'); Ypp(Ymf>=2.5) = 1;  YM      = Ygmt>0; 
    Ypp(YM) = min(Ycsfd(YM),Ygmt(YM) - Ywmd(YM)) ./ (Ygmt(YM) + eps); Ypp(Ypp>2) = 0;
    Ygmts   = cat_vol_localstat(Ygmt,Ygmt>0 & Ypp>0.1,1,1); Ygmt(Ygmts>0) = Ygmt(Ygmts>0)*(1-fs) + fs*Ygmts(Ygmts>0); clear Ygmtps
    
    % Error handling 
    % For some unkown reasons the sulcus reconstruction of cat_vol_pbtp failed in some cases (not directly reproducable).      
    % Reprocessing is solving this problem, but further investigation of cat_vol_pbtp.cpp would be good (RD 20190811).
    % Maybe it depends on the initialization of the regions, e.g., using Ymf without rounding and incorrect boundary seams to increase the problems.  
    % Now use cat_vol_pbtp2.cpp (RD20200111).
    mask   = @(Y) Y(:)>0 & Y(:)<1000000; 
    rerun = 0; rerunlim = 3; 
    while rerun <= rerunlim && isnan( mean( Ygmt(mask(Ygmt))) ) || mean( Ygmt(mask(Ygmt)))>100
      Ygmt = cat_vol_pbtp2((Ymf),Ywmd,Ycsfd);  
      rerun = rerun + 1; 
      pause(rand*3);
    end
    if rerun == rerunlim && isnan( mean( Ygmt(mask(Ygmt))) ) || mean( Ygmt(mask(Ygmt)))>100
      error('cat_vol_pbtp2:bad_mapping','Untypcial values in PBT thickness mapping detected. ');
    end
    
    % avoid meninges !
    Ygmt = min(Ygmt,Ycsfd+Ywmd);
    
    % median filter to remove outliers
    Ygmt = cat_vol_median3(Ygmt,Ygmt>0,Ygmt>0);
    
    %YM      = Ymf>=1.5 & Ymf<2.5; Ypp = zeros(size(Ymf),'single'); Ypp(Ymf>=2.5) = 1;
    %Ypp(YM) = min(Ycsfd(YM),Ygmt(YM) - Ywmd(YM)) ./ (Ygmt(YM) + eps); Ypp(Ypp>2) = 0;
    
    Ypp(Ypp==0 & Ymf>=2.5 & Ywmd<=0) = 1; 
    Ypp(Ypp==0 & Ymf>=2.5 & Ywmd<1)  = 0.99; 
    Ypp(Ypp==0 & Ymf>=2.0 & Ywmd<2)  = 0.985; 
    Ypp(Ypp==0 & Ymf>=2.5)           = 0.95; 

    Ypp = Ypp*(1-fs) + fs*cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & (Ypp>0.1 | Ypp==0),true(size(Ypp)),0.1);
    Ypp = Ypp*(1-fs) + fs*cat_vol_median3(Ypp,Ymf>1.1 & Ymf<2.9 & (Ypp>0.1 | Ypp==0),true(size(Ypp)),0.1);

    Ypp(Ypp>0.5 & ~cat_vol_morph(Ypp>0.5,'ldo',0.50/opt.resV)) = 0.49;
    Ypp(Ypp<0.5 & ~cat_vol_morph(Ypp<0.5,'l'  ,0.25/opt.resV)) = 0.51;

    
    % filter result
    if exist('Ymfo','var'); Ymf=Ymfo; end
    Ygmts = Ygmt; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,(Ygmt>1 | Ypp>0.1) & Ygmt>0 & (Ygmt>1 | Ymf>1.8),1,1); end; Ygmt(Ygmts>0) = Ygmts(Ygmts>0); 
    
    % filter result
    Ygmts = Ygmt; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,(Ygmt>1 | Ypp>0.1) & Ygmts>0 & (Ygmt>1 | Ymf>1.8),1,1); end; Ygmt(Ygmts>0) = Ygmts(Ygmts>0);
  end
  
  
 
  
  
  
  %% Estimation of a mixed percentual possion map Ypp.
  if nargout>3 || ~strcmp(opt.method,'pbtv')  
    YM  = Ymf>=1.5 & Ymf<2.5 & Ygmt>eps;
    Ywmdc = Ycsfd; Ywmdc(YM) = min(Ycsfd(YM),Ygmt(YM) - Ywmd(YM)); 
  end
  if ~strcmp(opt.method,'pbtv')  
    Ypp = zeros(size(Ymf),'single'); Ypp(Ymf>=2.5)=1;
    Ypp(YM) = Ywmdc(YM) ./ (Ygmt(YM) + eps); 
    Ypp(Ypp>2) = 0;
  end
 
 % YM  = (Ygmt<=opt.resV & Ywmd<=opt.resV & Ygmt>0); Ypp(YM) = (Ymf(YM)-1)/2 - 0.2; % correction of voxel with thickness below voxel resolution
  Ypp(isnan(Ypp)) = 0; 
  Ypp(Ypp<0) = 0; 
  
  %% Final corrections for position map with removing of non brain objects.
  % ds('d2','',1,Ymf/3,Ywmd/3,Ygmt/5,Ypp,70)
  
  % Final corrections for thickness map with thickness limit of 10 mm. 
  % Resolution correction of the thickness map after all other operations, 
  % because PBT actually works only with the voxel-distance (isotropic 1 mm)
  Ygmt = Ygmt*opt.resV; 
  Ygmt(Ygmt>10) = 10; 
  if exist('Ycsfdc','var'), Ywmdc = Ywmdc*opt.resV; end
  if exist('Ywmd','var'), Ywmd = Ywmd*opt.resV; end
  
  cat_io_cmd(' ','g5','',opt.verb,stime);
  if opt.debug, cat_io_cmd(' ','','',opt.debug,stime2); end
end