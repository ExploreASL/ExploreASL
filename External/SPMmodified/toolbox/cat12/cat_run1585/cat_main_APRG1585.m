function [Yb,Ym0,Yg,Ydiv] = cat_main_APRG1585(Ysrc,P,res,T3th)
% ______________________________________________________________________
%  Skull-stripping subfunction APRG (adaptive probability region-growing)
%  of cat_main_updateSPM.
%
%  [Yb,Ym0,Yg,Ydiv] = cat_main_updateSPM1585(Ysrc,P,res)
%
%   Ysrc  .. original input images
%   P     .. tissue segmenation (4D)
%   res   .. SPM preprocessing structure
%   Yb    .. binary brain mask
%   Ym0   .. probability brain mask
%   Yg    .. absolute gradient map (eg. for tissue edges)
%   Ydiv  .. divergence maps (eg. for blood vessels)
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_APRG.m 1435 2019-03-06 11:24:38Z dahnke $

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  % adaptive probability region-growing
  if debug, Po=P; end

  vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
  clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;

  %% tissue treshholds depending on MR modality
  if T3th(1) < T3th(3) % T1
    cth = min(res.mn(res.lkp==3));
    cth = min(cth,T3th(1));
    cth = min(cth,cat_stat_nanmedian(Ysrc( cat_vol_morph( smooth3(P(:,:,:,3))>200 & ...
            Ysrc<sum(T3th(1:2).*[0.8 0.2]) ,'de',2,vx_vol) ) ));
  else
    cth = max(res.mn(res.lkp==3));
    cth = max(cth,T3th(1));
    cth = max(cth,cat_stat_nanmedian(Ysrc( cat_vol_morph( smooth3(P(:,:,:,3))>200 & ...
            Ysrc>sum(T3th(1:2).*[0.8 0.2]),'de',2,vx_vol) ) ));
  end  
  if max(res.lkp)==4
    bth = 0;
  else
    bth = min( res.mn(res.lkp==4) );
  end
  if T3th(1) < T3th(3) % T1: CSF<GM<WM
    tth(1,:)  = [ mean(T3th(1:2))                     mean(T3th(2:3))                ]; % GM
    tth(2,:)  = [ mean(T3th(2:3))                     T3th(3) + 0.25*diff(T3th(2:3)) ]; % WM
    tth(3,:)  = [ T3th(1) - 0.25*diff([bth,T3th(1)])  sum(T3th(1:2) .* [0.75 0.25])  ]; % CSF
  elseif T3th(1) > T3th(3) % T2/PD: WM<GM<CSF ... not tested
    tth(1,:)  = [ mean(T3th(2:3))                 mean(T3th(1:2))                ];
    tth(2,:)  = [ T3th(3) - 0.25*diff(T3th(2:3))  mean(T3th(2:3))                ];
    tth(3,:)  = [ sum(T3th(1:2) .* [0.75 0.25])   sum(T3th(1:2) .* [0.75 0.25]) ];
  else % other contrast
    tth(1,:)  = [ T3th(1) - std(T3th(1:2))        T3th(1) + std(T3th(1:2)) ];
    tth(2,:)  = [ T3th(2) - std(T3th(1:2))        T3th(1) + std(T3th(2:3)) ];
    tth(3,:)  = [ T3th(3) - std(T3th(2:3))        T3th(1) + std(T3th(2:3)) ];
  end


  %% CSF mask
  %  Yc .. Combination of the CSF probability map and intensity map to 
  %        avoid meninges (especially in older subjectes) at the outer
  %        boundary. Due to failed registration we directly use the 
  %        brain mask Yb in the center of the brain, ie. we allow brain 
  %        tissue in the CSF far from the skull. 
  if T3th(1) < T3th(3)
    Yc   = single(P(:,:,:,3))/255 .* ...
      max(0,min(1,1 - ( max(0, Ysrc - cth) / abs( mean(res.mn(res.lkp(:)==1).*res.mg(res.lkp(:)==1)' ) - cth ) + ...
                        max(0,-Ysrc + cth) / abs( mean(res.mn(res.lkp(:)==4).*res.mg(res.lkp(:)==4)' ) - cth ) ) ));
  else
    Yc   = single(P(:,:,:,3))/255 .* ...
      max(0,min(1,1 - ( max(0,Ysrc - cth) / abs( min(res.mn(res.lkp==3)) - mean(res.mn(res.lkp==1)) )) ));
  end
  Ycg  = (single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 - Yc) .* ...
    max(0,min(1,1 - (0.5*abs(Ysrc - T3th(2)) / abs(diff(T3th(1:2))) ) ));


  % improved brain mask by region growing
  %  Yb .. improving the brain mask is necessary in case of missing
  %        structures (e.g. in children) or failed registration where 
  %        the TPM does not fit well and the CSF include brain tissue
  %        simply by chance.
  BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
  Yg   = cat_vol_grad((Ysrc - BGth)/diff([BGth,T3th(3)]),vx_vol);
  BVth = abs(diff(T3th(1:2:3)) / abs(T3th(3)) * 3);   % avoid blood vessels (high gradients) 
  RGth = abs(diff(T3th(2:3))   / abs(T3th(3)) * 0.1); % region growing threshold




  %% initial brain mask by region-growing
% Todo: T2/PD, skull-stripped
  %  as CSF/GM+GM+WM without blood vessels (Yg<0.5) 
  Yb   = min(1,single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255 + Ycg/2); 
  Yb   = cat_vol_median3(Yb,Yb>0,true(size(Yb)));
  if T3th(1) < T3th(3) % T1 
    Yb   = min(1,Yb + cat_vol_morph(smooth3(sum(P(:,:,:,1:3),4))>253,'ldc',1,vx_vol));
  else % T2/PD ... more tolerant
    Yb   = min(1,Yb + cat_vol_morph(smooth3(sum(P(:,:,:,1:3),4))>200,'ldc',1,vx_vol));
  end
  Yb   = cat_vol_morph(Yb>0.8,'ldo',1.9,vx_vol);
  Yb   = cat_vol_morph(Yb,'ldc',1.9,vx_vol);
  
  
  %% mask for region growing and WM-GM region growing
  Yb2  = single(cat_vol_morph(Yb,'de',1.9,vx_vol)); 
  if T3th(1) < T3th(3) % T1 
    Yh   = (Yb2<0.5) & (Ysrc<sum(T3th(2:3).*[0.75 0.25]) | Ysrc>(T3th(3)*1.2) | Yg>BVth);
  else
    Yh   = (Yb2<0.5) & (Ysrc<mean([T3th(3),BGth]) | Ysrc>sum(T3th(2:3).*[0.75 0.25]) | Yg>BVth);
  end
  Yh   = cat_vol_morph(Yh,'dc',1,vx_vol); 
  Yh   = cat_vol_morph(Yh,'de',1,vx_vol); Yb2(Yh) = nan; if ~debug, clear Yh; end
  if T3th(1) < T3th(3) % T1 
    [Yb2,YD] = cat_vol_downcut(Yb2,Ysrc/T3th(3),RGth/2); clear Yb2; %#ok<ASGLU>
  else
    [Yb2,YD] = cat_vol_downcut(Yb2,1 - Ysrc/T3th(3),RGth/2); clear Yb2; %#ok<ASGLU>
  end
  Yb(YD<400/mean(vx_vol)) = 1; clear YD; 
  Yb(smooth3(Yb)<0.5) = 0; 
  Yb   = cat_vol_morph(Yb,'ldo',1.9,vx_vol);
  Yb   = cat_vol_morph(Yb,'ldc');
  

  %% GM-CSF region
  Yb2  = single(cat_vol_morph(Yb,'de',1.9,vx_vol)); 
  if T3th(1) < T3th(3)
    Yh   = (Yb2<0.5) & (Ysrc<sum(T3th(1:2).*[0.9 0.1]) | sum(P(:,:,:,4:6),4)>250 | ...
            Ysrc>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth); 
  else
    Yh   = (Yb2<0.5) & (Ysrc>sum(T3th(1:2).*[0.9 0.1]) | sum(P(:,:,:,4:6),4)>250 | ...
            Ysrc<(T3th(3) - sum(T3th(2:3).*[0.5 0.5])) | Yg>BVth); 
  end
  Yh   = cat_vol_morph(Yh,'dc',1) | cat_vol_morph(~Yb,'de',10,vx_vol); 
  Yh   = cat_vol_morph(Yh,'de',1,vx_vol);  Yb2(Yh) = nan; if ~debug, clear Yh; end
  if T3th(1) < T3th(3) % T1 
    [Yb2,YD] = cat_vol_downcut(Yb2,Ysrc/T3th(3),-RGth); clear Yb2; %#ok<ASGLU>
  else
    [Yb2,YD] = cat_vol_downcut(Yb2,1 - Ysrc/T3th(3),-RGth); clear Yb2; %#ok<ASGLU>
  end
  Yb(YD<400/mean(vx_vol)) = 1; clear YD; 
  Yb(smooth3(Yb)<0.5) = 0; Yb(smooth3(Yb)>0.5) = 1; 
  Yb   = cat_vol_morph(Yb,'ldo',1.9,vx_vol);
  Yb   = cat_vol_morph(Yb,'ldc');



  %% CSF mask 2 with Yb
  %  Yc .. Combination of the CSF probability map and intensity map to 
  %        avoid meninges (especially in older subjectes) at the outer
  %        boundary. Due to failed registration we directly use the 
  %        brain mask Yb in the center of the brain, ie. we allow brain 
  %        tissue in the CSF far from the skull. 
  cth2 = min(cth,cat_stat_nanmedian( Ysrc( smooth3( ...
    cat_vol_morph(P(:,:,:,3)>200 & Yb & Yg<0.2,'de',1.5))>0.9 ) ));
  if T3th(1) < T3th(3)
    Yc   = single(P(:,:,:,3) )/255 .* ...
      max(Yb & Yg<BVth & Ysrc>tth(3,1) & Ysrc<tth(3,2), ...
      min(1,1 - ( max(0, Ysrc - cth2) / (2 * abs( mean(res.mn(res.lkp(:)==1).*res.mg(res.lkp(:)==1)' - cth2) )) + ...
                  max(0,-Ysrc + cth2) / (2 * abs( mean(res.mn(res.lkp(:)==4).*res.mg(res.lkp(:)==4)' - cth2) )) ) ));
  else
    Yc   = single(P(:,:,:,3))/255 .* ...
      max(Yb & Yg<BVth & Ysrc>tth(3,1) & Ysrc<tth(3,2), ...
      min(1,1 - ( min(0,Ysrc - cth2) / (2 * abs( mean(res.mn(res.lkp==1)) - cth2) ) ) ));
  end
  % single(P(:,:,:,3))/255 .* ...
  %  max(Yb & Ysrc>cth2,min(1,1 - ( abs(Ysrc - cth2) / ...
  %  (2 * abs( mean(res.mn(res.lkp==1)) - cth2) ) ) ));




  %% update GM map (include tissue that was previously labeled as head) 
  for i=1:3
    % smaller mask in case of WM
    if i==2, Ybt = cat_vol_morph(Yb,'de',1.5,vx_vol); else Ybt = Yb; end
    % smooth mask in case of GM
    if i==1
      Ytmp = single(P(:,:,:,4)) .* smooth3(Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth); 
    else %if i==2
      Ytmp = single(P(:,:,:,4)) .*        (Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth); 
    %{
    elseif i==3
      % add CSF around mask only if the intensity fit very well
      Ytmp = single(P(:,:,:,4)) .* Yc .* ...
        max(Yb & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth, ...
        min(1,1 - ( abs(Ysrc - cth2) / (4 * abs( mean(res.mn(res.lkp==1)) - cth2) ) ) ));
    %}
   end
    % tissue transfer
    P(:,:,:,i) = cat_vol_ctype(single(P(:,:,:,i)) + Ytmp);
    P(:,:,:,4) = cat_vol_ctype(single(P(:,:,:,4)) - Ytmp);
  end
  % transfer tissue in case of bad TPM matching
  for i=[1,2,4]
    % smaller mask in case of WM
    if i==2, Ybt = cat_vol_morph(Yb,'de',1.5,vx_vol); else Ybt = Yb; end
    % smooth mask in case of GM
    if i==1
      Ytmp = single(P(:,:,:,3)) .* smooth3(Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth);
    elseif i==2
      Ytmp = single(P(:,:,:,3)) .*        (Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth); 
    elseif i==4
      Ytmp = single(P(:,:,:,3)) .* (1 - min(1,Ybt + Yc)); 
    end
    % tissue transfer
    P(:,:,:,i) = cat_vol_ctype(single(P(:,:,:,i)) + Ytmp);
    P(:,:,:,3) = cat_vol_ctype(single(P(:,:,:,3)) - Ytmp);
  end


  %% qa control parameter?
  % SPMsegfit = mean(P(Yb(:))-Po(Yb(:))); 


  %% create brain level set map
  %  Ym .. combination of brain tissue and CSF that is further corrected
  %        for noise (median) and smoothness (Laplace) an finally 
  %        threshholded 
  Ym  = min(1,Yc + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255 + Yb);
  Ym  = cat_vol_median3(Ym,Ym>0 & Ym<1);  % remove noise 
  % region-growing 
  Ym2 = Ym; Ym2(Ym2==0)=nan;
  [Ym2,YD] = cat_vol_downcut(single(Ym2>0.99),Ym2,0.01); clear Ym2; %#ok<ASGLU>
  Ym(YD>400/mean(vx_vol))=0; clear YD; 
  Ym(cat_vol_morph(Ym>0.95,'ldc',1)) = 1; 
  Ym(cat_vol_morph(Yb,'e') & Ym<0.9 & Yc<0.25) = 0;
  Ym  = Ym .* cat_vol_morph(Ym>0.5,'ldo',2);  % remove extensions (BV, eye)
  Ym = cat_vol_laplace3R(Ym,Ym>0.1 & Ym<0.9,0.2); % smooth mask
  Ym = cat_vol_laplace3R(Ym,Ym<0.25,0.2); 
  Ym(cat_vol_morph(Yb,'e') & Ym<0.9 & Yc<0.25) = 0;

  
  %% cutting parameter
  %  This is a nice parameter to control the CSF masking.
  %  The default value of 0.5 which works best for validation data.
  %  However, by setting the value to 1 we can also use a surface to find the 
  %  optimal value, which is currently not stable enough.
  cutstr    = 0.5; 
  cutstrs   = linspace(0.05,0.95,4); % 0.05,0.35,0.65,0.95]; 
  cutstrval = nan(1,4); 
  if debug, cutstrsa = zeros(0,8); end
  Ysrc2 = (Ysrc>T3th(1)) .* (abs(Ysrc - T3th(1))/(T3th(2) - T3th(1))) + ...
          (Ysrc<T3th(1)) .* (abs(Ysrc - T3th(1))/(T3th(1) - BGth)) ;
  Ysrc2 = smooth3(Ysrc2);
  
  % automatical estimation of cutting parameter is currently not stable enough!
  if cutstr == 1 % auto
    for l=1:5
      for i=1:numel(cutstrs)
        if isnan( cutstrval(i) )
          S = isosurface(Ym,cutstrs(i),Ysrc2); 
          cutstrval(i) = cutstrs(i)/10 + cat_stat_nanmean(S.facevertexcdata.^2).^0.5; 
          %... % litte offset to get more CSF% + std(S.facevertexcdata);
        end
      end
      [tmp,cutstrid] = sort(cutstrval); clear tmp; %#ok<ASGLU>
      if debug, cutstrsa  = [cutstrsa; cutstrs, cutstrval]; end %#ok<AGROW>
      cutstrs   = linspace(cutstrs(max(1,cutstrid(1)-1)),cutstrs(min(4,cutstrid(1)+1)),4);
      cutstrval = [cutstrval(max(1,cutstrid(1)-1)),nan,nan,cutstrval(min(4,cutstrid(1)+1))];

    end
    cutstr = cutstrs(cutstrid(1));
  end


  %% normalize this map depending on the cutstr parameter 
  Yb  = cat_vol_morph(cat_vol_morph(Ym > cutstr,'lo'),'lc',2);
  Yb  = cat_vol_morph(Yb,'e') | (Ym>0.9) | (Yb & Yc>0.5);
  Yb(smooth3(Yb)<0.5)=0;
  Ybb = cat_vol_ctype( max(0,min(1,(Ym - cutstr)/(1-cutstr))) * 256); 
  Ym0 = Ybb; 

  
  %% estimate gradient (edge) and divergence maps
  [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
  Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
  Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
  Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
  Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
end