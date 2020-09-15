function [P,res,stime2] = cat_main_kamap(Ysrc,Ycls,Yy,tpm,job,res,vx_vol,stime2)
% ______________________________________________________________________
%
% (initial) k-means AMAP segmentation
%
% This is an alternative pipeline in case of failed SPM brain tissue
% classification in datasets with abnormal anatomy, i.e. superlarge 
% ventricle. However, SPM is used for head tissue classification and
%  bias correction. It is (only) called from the cat_main function. 
%
%   [P,res,stime2] = cat_main_amap(Ysrc,Ycls,job,res,stime2)
%
%   Ysrc    .. cell structure of 3D probability maps of the segmentation 
%              [GM,WM,CSF,SK,HD,BG] 
%   P       .. 4D probability map of the segmentation
%   res     .. SPM segmentation structure
%   vx_vol  .. voxel_resolution
%   stime2  .. just a time stamp
% ______________________________________________________________________
%
% Robert Dahnke, Christian Gaser 
% Structural Brain Mapping Group
% University Jena 
% ______________________________________________________________________
% $Id: cat_main_kamap.m 1435 2019-03-06 11:24:38Z dahnke $ 


% ______________________________________________________________________
% Possible developments:
% (1) Update for PD/T2 weightings?
%     Yes, in future!
% (2) Add head class model? 
%     Maybe later.
% (-) More complex skull-stripping? 
%     No, this is given by the input segmenation  
% ______________________________________________________________________


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs  = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  
  stime2 = cat_io_cmd('  K-means AMAP bias-correction','g5','',job.extopts.verb-1,stime2);

  
  % Create brain mask based on the the TPM classes
  % cleanup with brain mask - required for ngaus [1 1 2 4 3 2] and R1/MP2Rage like data 
  d   = res.image(1).dim(1:3);
  YbA = zeros(d,'single');
  Vb  = tpm.V(1); Vb.pinfo(3) = 0; Vb.dt=16; 
  Vb.dat = single(exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3})); 
  for z=1:d(3)
    YbA(:,:,z) = spm_sample_vol(Vb,double(Yy(:,:,z,1)),double(Yy(:,:,z,2)),double(Yy(:,:,z,3)),1); 
  end
  if round(max(YbA(:))/Vb.pinfo(1)), YbA=YbA>0.1*Vb.pinfo(1); else YbA=YbA>0.1; end
  % add some distance around brainmask (important for bias!)
  YbA = YbA | cat_vol_morph(YbA & (Ycls{1} + Ycls{2})>4 ,'dd',2.4,vx_vol);
  
  
  
  %% get some low intensity values (backgroudn or CSF)
  %  don't trust SPM background value in case of large ventricles
  [Ysrc2,th] = cat_stat_histth(Ysrc); clear Ysrc2;
  T4th = [ th(1) ... 
    sum(res.mn(res.lkp==3).*res.mg(res.lkp==3)') ...
    sum(res.mn(res.lkp==1).*res.mg(res.lkp==1)') ...
    sum(res.mn(res.lkp==2).*res.mg(res.lkp==2)') ...
    ];
  % add CSF that was missclassified as background
  if  T4th(3) < T4th(4) % T1 
    newCSF = YbA & Ysrc>sum(T4th(1:2).*[2/3 1/3]) & Ysrc<sum(T4th(2:3).*[0.5 0.5]);
  else % T2 / PD
    newCSF = YbA & Ysrc>sum(T4th(2:3).*[0.5 0.5]);
  end
  % modify tissue classes
  if numel(Ycls)>5
    newCSF = newCSF .* (single(Ycls{4}) + single(Ycls{5}) + single(Ycls{6})); 
  else
    newCSF = newCSF .* single(Ycls{4});
  end
  newCSF = smooth3(newCSF/255)>0.5;
  newCSF = cat_vol_morph(newCSF,'do',1.9);
  % transfer CSF like voxel from class 4-6 to CSF class 3
  for ci=4:numel(Ycls), Ycls{ci} = cat_vol_ctype( single(Ycls{ci}) .* single(~newCSF)); end
  Ycls{3} = cat_vol_ctype( single(Ycls{3}) + single(newCSF)*255 );
  
  
  
  %% simple skull-stripping
  prob = cat_main_clean_gwc(cat(4,Ycls{1},Ycls{2},Ycls{3}),1,0); % no new cleanup here!
  Ycls{1} = prob(:,:,:,1); Ycls{2} = prob(:,:,:,2); Ycls{3} = prob(:,:,:,3); clear prob;
  Yb = cat_vol_ctype( single(Ycls{1} + Ycls{2} + Ycls{3}) + 255*newCSF );
  
  
  % In case of skull-stripping with have to combine SPM segmentation and 
  % the possible skull-stripping information 
  if (max(res.lkp)==4), Yb = Yb .* uint8(Ysrc>0); end % if skullstripped
  Yb = cat_vol_morph( smooth3(Yb)>64 , 'ldc', 8,vx_vol);
  Yb = cat_vol_morph( Yb , 'ldo', 8,vx_vol);
 
  
  %% intensity normalization
  Ym = Ysrc; Ym(~cat_vol_morph(Yb,'e',3)) = nan; 
  [Ymic,th] = cat_stat_histth(Ym,0.98); 
  Ym(isnan(Ym) | Ym>(th(2) + diff(th)*2)) = 0;
  [T3th,T3sd,T3md] = kmeans3D(Ymic((Ym(:)>0)),5); clear Ymic;  %#ok<ASGLU>
  if T3md(end)<0.08, T3th(end)=T3th(end-1); end
  Ym = (Ysrc - th(1)) ./ (T3th(end) - th(1)); 
  Tth.T3thx  = [T3th(1)-diff(T3th(1:2:3)) T3th(1:2:5) T3th(5)+diff(T3th(3:2:5))];
  Tth.T3th   = 0:1/3:4/3;
  Ymi = cat_main_gintnormi(Ysrc/3,Tth);
  
  %% update CSF regions for skull-stripping of the cleanup and update skull-stripping
  Yhd  = cat_vbdist(single(Yb<0.5)); Yhd = Yhd ./ max(Yhd(Yhd <100));  
  Ycls{3}(~newCSF & single(Ycls{3})/255 .* abs(Ymi*3 - 1) .* cat_vol_morph(Ycls{1}==0,'e')>0.3 & Yhd<0.2) = 0; 
  Yb = cat_vol_ctype( single(Ycls{1} + Ycls{2} + Ycls{3}) + 255*newCSF );
  %%
  Yb = cat_vol_morph( smooth3(Yb)>64 , 'ldo', 1,vx_vol);
  Yb = cat_vol_morph( Yb , 'ldc', 8,vx_vol);
  Yb = cat_vol_morph( Yb , 'ldo', 8,vx_vol);
  if ~debug, clear newCSF; end
  
  %%
  Ymsk = cat_vol_morph(Yb & Ymi>0.45 & Ymi<0.95,'do',1); 
  tx = kmeans3D(Ysrc(Ymsk(:)),1); T3th(3) = mean(tx);
  Ymsk = cat_vol_morph(Yb & Ymi>0.9 & Ymi<1.1,'do',1.5); 
  tx = kmeans3D(Ysrc(Ymsk(:)),1); T3th(5) = mean(tx);
  Tth.T3thx  = [T3th(1)-diff(T3th(1:2:3)) T3th(1:2:5) T3th(5)+diff(T3th(3:2:5))];
  Ymi = cat_main_gintnormi(Ysrc/3,Tth);
 
  %Yb(Ym>1.2) = 0;
  %clear th; 

  
  %% bias correction for white matter (Yw) and ventricular CSF areas (Yv)
  Yw  = Yb & Ymi>0.9 & Ymi<1.5; Yw(smooth3(Yw)<0.6) = 0; Yw(smooth3(Yw)<0.6) = 0;
  %Yg  = Ym<0.95 & Ym>0.45 & Yb; Yg(smooth3(Yg)<0.4) = 0; Yg(smooth3(Yg)<0.6) = 0; % not used?
  Yv  = Yb & Ymi<0.4 & Yb; Yv  = cat_vol_morph(Yv,'e',4); 
  
  %% estimate value vth to mix CSF and WM information
  Yvw = cat_vol_smooth3X(Yv,6)>0.05 & cat_vol_morph(Yw,'e',2); 
  Ywv = cat_vol_smooth3X(Yw,6)>0.05 & cat_vol_morph(Yv,'e',2); 
  if sum(Yvw(:))>10 && sum(Yvw(:))>10
    vth = cat_stat_nanmedian(Ysrc(Yvw(:))) ./ cat_stat_nanmedian(Ysrc(Ywv(:)));
  elseif  sum(Yvw(:))>0 && sum(Yvw(:))>0
    vth = cat_stat_nanmedian(Ysrc(Yv(:))) ./ cat_stat_nanmedian(Ysrc(Yw(:)));
  else
    vth = 1; 
  end
  if debug==0, clear Ywv Yvw; end

  
  %% some precorrections and bias field approximation 
  Yi   = cat_vol_localstat(Ysrc .* Yw,Yw,1,3);
  Yiv  = cat_vol_localstat(Ysrc .* Yv,Yv,1,1);
  Yi   = Yi + Yiv * vth; if debug==0, clear Yiv; end
  Yi(Yb==0) = T3th(5);
  Yi   = cat_vol_approx(Yi,'nn',vx_vol,4); 
  Yi   = cat_vol_smooth3X(Yi,4); 
  Ysr2 = Ysrc ./ (Yi ./ median(Yi(Yw(:))));
  Ymi  = Ysr2 ./ Yi  .* Yb;
  clear Ysr2 Yi; 
  
  
  %% intensity normalisation
  [T3th,T3sd,T3md] = kmeans3D(Ymi(Yb(:)),5);  %#ok<ASGLU>
  Tth.T3thx  = [T3th(1)-diff(T3th(1:2:3)) T3th(1:2:5) T3th(5)+diff(T3th(3:2:5))];
  Tth.T3th   = 0:1/3:4/3;
  Ymi2 = cat_main_gintnormi(Ymi/3,Tth) .* Yb;
  %%
  Ymsk = cat_vol_morph(Yb & Ymi2>0.45 & Ymi2<0.95,'do',1); 
  tx = kmeans3D(Ymi2(Ymsk(:)),1); T3th(3) = mean(tx);
  Ymsk = cat_vol_morph(Yb & Ymi2>0.9 & Ymi2<1.1,'do',1.5); 
  tx = kmeans3D(Ymi2(Ymsk(:)),1); T3th(5) = mean(tx);
  clear Ymi2;
  
  Tth.T3thx  = [T3th(1)-diff(T3th(1:2:3)) T3th(1:2:5) T3th(5)+diff(T3th(3:2:5))];
  Ymi = cat_main_gintnormi(Ymi/3,Tth) .* Yb;
  
  
  
  %% remove of high intensity structures
  Ymi = cat_vol_median3(Ymi/median(Ymi(Yw(:))),Yb,Yb,0.2)*median(Ymi(Yw(:)));
  if debug==0, clear Ym Yi Yw Yv; end


  %%  similar to later AMAP call
  %   -------------------------------------------------------------------

  % correct for harder brain mask to avoid meninges in the segmentation
  Ymib = min(1.5,Ymi); Ymib(~Yb | Ymi>1.2) = 0; 
  rf = 10^4; Ymib = round(Ymib*rf)/rf;
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));

  fprintf('prepare data for segmentation:  ');
  %  prepare data for segmentation
  % more direct method ... a little bit more WM, less CSF
  Yp0  = (Ymi * 3); Yp0(~Yb | Ymi>1.2) = 0; 
  Yp0c = cat_vol_ctype(255*cat(4,Yp0toC(Yp0,2),Yp0toC(Yp0,3),Yp0toC(Yp0,1))); 
  Yp0c = cat_main_clean_gwc(Yp0c,1,1); % new cleanup here!
  Yp0  = uint8(single(Yp0c(:,:,:,3))/255 + 2*single(Yp0c(:,:,:,1))/255 + ...
    3*single(Yp0c(:,:,:,2))/255); clear Yp0c;
  %Yp0 = cat_vol_ctype(max(1,min(3,single(Yb) + (Ymi>T3th(2)/T3th(5)) + ...
  %  (Ymi>sum(T3th(4:5).*[0.6 0.4])/T3th(5))))); Yp0(~Yb) = 0;
  
  % use index to speed up and save memory
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

  % Yb source image because Amap needs a skull stripped image
  % set Yp0b and source inside outside Yb to 0
  Yp0b = Yp0(indx,indy,indz);  %#ok<NASGU>
  Ymib = Ymib(indx,indy,indz); 

  % adaptive mrf noise 
  if job.extopts.mrf>=1 || job.extopts.mrf<0; 
    % estimate noise
    [Yw,Yg] = cat_vol_resize({Ymi.*(Ycls{1}>240),Ymi.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.15,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * 0.5; 
    clear Yn Yg
  end

  % display something
  stime2 = cat_io_cmd(sprintf('  Amap using k-means segmentations (MRF filter strength %0.2f)',...
    job.extopts.mrf),'g5','',job.extopts.verb-1,stime2);

  
  %% Amap parameters  - default sub=16 caused errors with highres data!
  % don't use bias_fwhm, because the Amap bias correction is not that efficient
  % and also changes intensity values
  Ymib = double(Ymib); n_iters = 200; sub = round(128/min(vx_vol)); %#ok<NASGU>
  n_classes = 3; pve = 5; bias_fwhm = 120; init_kmeans = 1;  %#ok<NASGU>
  if job.extopts.mrf~=0, iters_icm = 50; else iters_icm = 0; end %#ok<NASGU>

  
  % do segmentation and rep
  evalc(['prob = cat_amap(Ymib, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ' ...
    'job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);']);
  clear Ymib Yp0b;
 
  
  % reorder probability maps according to spm order
  prob = prob(:,:,:,[2 3 1]);  

  
  % cleanup
  prob = cat_main_clean_gwc(prob,1,1); % and here!

  
  % update probability maps
  P = zeros([size(Ycls{1}) numel(Ycls)],'uint8');
  for i=1:numel(Ycls), P(:,:,:,i) = Ycls{i}; end
  for i=1:3
     P(:,:,:,i) = 0; P(indx,indy,indz,i) = prob(:,:,:,i);
  end
  clear prob;
  Ys = cat_vol_ctype(255 - sum(P(:,:,:,1:3),4));
  for i=4:size(P,4)
    P(:,:,:,i) = P(:,:,:,i) .* Ys; 
  end
  clear Ys;

  sP = (sum(single(P),4) + eps)/255;
  for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end

  
  % update SPM segmentation information 
  for i=1:3
    [res.mn(res.lkp==i),tmp,res.mg(res.lkp==i)] = ...
      kmeans3D(Ysrc(P(:,:,:,i)>64),sum(res.lkp==i)); %#ok<ASGLU>
  end

end
