function [Ysrc,Ycls,Yb,Yb0,job,res,T3th,stime2] = cat_main_updateSPM(Ysrc,P,Yy,tpm,job,res,stime,stime2)
% ______________________________________________________________________
%  Update SPM preprocessing. 
%  Subfunction of cat_main.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_updateSPM.m 1561 2020-02-04 15:49:34Z gaser $


  %dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  global cat_err_res; 

  clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
  
  [pth,nam] = spm_fileparts(res.image0(1).fname); %#ok<ASGLU> % original   
  
  % voxel size parameter
  vx_vol  = sqrt(sum(res.image(1).mat(1:3,1:3).^2));    % voxel size of the processed image
  vx_volp = prod(vx_vol)/1000;
  
  d = res.image(1).dim(1:3);

  
  stime2 = cat_io_cmd('  Update Segmentation','g5','',job.extopts.verb-1,stime2); 
  

  % Create brain mask based on the the TPM classes
  % cleanup with brain mask - required for ngaus [1 1 2 4 3 2] and R1/MP2Rage like data 
  YbA = zeros(d,'single');
  Vb = tpm.V(1); Vb.pinfo(3) = 0; Vb.dt=16; 
  Vb.dat = single(exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3})); 
  for z=1:d(3)
    YbA(:,:,z) = spm_sample_vol(Vb,double(Yy(:,:,z,1)),double(Yy(:,:,z,2)),double(Yy(:,:,z,3)),1); 
  end
  if round(max(YbA(:))/Vb.pinfo(1)), YbA=YbA>0.1*Vb.pinfo(1); else YbA=YbA>0.1; end
  % add some distance around brainmask (important for bias!)
  YbA = YbA | cat_vol_morph(YbA & sum(P(:,:,:,1:2),4)>4 ,'dd',2.4,vx_vol);
  
  
  % transfer tissue outside the brain mask to head  ... 
  % RD 201807: I am not sure if this is a good idea. Please test this with children! 
  for i=1:3
    P(:,:,:,4) = cat_vol_ctype(single(P(:,:,:,4)) + single(P(:,:,:,i)) .* single(~YbA)); 
    P(:,:,:,i) = cat_vol_ctype(single(P(:,:,:,i)) .* single(YbA)); 
  end
  clear YbA;
  
  
  % Cleanup for high resolution data
  % Alghough the old cleanup is very slow for high resolution data, the   
  % reduction of image resolution removes spatial segmentation information. 
  if job.opts.redspmres==0 % already done in case of redspmres
    if max(vx_vol)<1.5 && mean(vx_vol)<1.3
      for i=1:size(P,4), [Pc1(:,:,:,i),RR] = cat_vol_resize(P(:,:,:,i),'reduceV',vx_vol,job.extopts.uhrlim,32); end %#ok<AGROW>
      Pc1 = cat_main_clean_gwc(Pc1,1);
      for i=1:size(P,4), P(:,:,:,i)   = cat_vol_resize(Pc1(:,:,:,i),'dereduceV',RR); end 
      clear Pc1 Pc2;
    end
  end

  
  % garantee probability 
  sP = (sum(single(P),4)+eps)/255;
  for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
  clear sP;

  
  % Use median for WM threshold estimation to avoid problems in case of WMHs!
  WMth = double(max( clsint(2) , cat_stat_nanmedian(Ysrc(P(:,:,:,2)>192)) )); 
  if clsint(3)>clsint(2) % invers
    CMth = clsint(3); 
  else
    CMth = min( [  clsint(1) - diff([clsint(1),WMth]) , clsint(3) ]);
  end
  T3th = [ CMth , clsint(1) , WMth];


  %% Some error handling
  %    ds('l2','',vx_vol,Ysrc./WMth,Yp0>0.3,Ysrc./WMth,Yp0,80)
  Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
    res.YlesionFull = res.Ylesion; % ExploreASL fix
    res.Ylesion = cat_vol_ctype( single(res.Ylesion) .* (Yp0>0.2) ); 
    for k=1:size(P,4), Yl = P(:,:,:,k); Yl(res.Ylesion>0.5) = 0; P(:,:,:,k) = Yl; end  
    Yl = P(:,:,:,3); Yl(res.Ylesion>0.5) = 255; P(:,:,:,3) = Yl; clear Yl; 
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  end
  if sum(Yp0(:)>0.3)<100 
    % this error often depends on a failed affine registration, where SPM
    % have to find the brain in the head or background
    BGth  = min(cat_stat_nanmean(Ysrc( P(:,:,:,end)>128 )),clsint(6));
    HDHth = clsint(5);
    HDLth = clsint(4);
    clsvol = nan(1,size(P,4)); for ci=1:size(P,4), Yct = P(:,:,:,ci)>128; clsvol(ci) = sum(Yct(:))*vx_volp; end; clear Yct; 
    if size(P,4)==6
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Please check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM','HDH','HDL', ...
          [ clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol)*100, clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th,HDHth,HDLth]));  %#ok<SPERR>
    elseif size(P,4)==4 % skull-stripped
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Please check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM', ...
          [ clsvol([4 3 1 2])/cat_stat_nansum(clsvol)*100, clsvol([4 3 1 2])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th]));  %#ok<SPERR>
    else
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ['Empty Segmentation: ' ...
           'Possibly the affine registration failed. Please check image orientation.\n']); 
    end
  end

  
  
  %%
  Yp0(smooth3(cat_vol_morph(Yp0>0.3,'lo'))<0.5)=0; % not 1/6 because some ADNI scans have large "CSF" areas in the background 
  Yp0     = Yp0 .* cat_vol_morph(Yp0 & (Ysrc>WMth*0.05),'lc',2);
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));

  
  % values are only used if errors occur
  cat_err_res.init.T3th = T3th; 
  cat_err_res.init.subjectmeasures.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)), ... CSF
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)), ... GM 
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3)), ... WM
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),4))];  % WMH
  cat_err_res.init.subjectmeasures.vol_TIV     =  sum(cat_err_res.init.subjectmeasures.vol_abs_CGW); 
  cat_err_res.init.subjectmeasures.vol_rel_CGW =  cat_err_res.init.subjectmeasures.vol_abs_CGW ./ ...
                                                  cat_err_res.init.subjectmeasures.vol_TIV;
  [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
  cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);
  clear Yp0; 

  % ### This can not be reached because the mask field is removed by SPM! ###
  if isfield(res,'msk') 
    Ybg = ~res.msk.dat; 
    P4  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg==0) + single(P(:,:,:,4)) .* (Ybg<1) ); % remove air in head
    P5  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg==0) + single(P(:,:,:,5)) .* (Ybg<1) ); % remove air in head
    P6  = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg==1) + single(P(:,:,:,6)) .* (Ybg>0) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6 Ybg; 
  end
  

  
  
  %% Skull-Stripping
  %  ----------------------------------------------------------------------
  %  Update Skull-Stripping 1
  %  ----------------------------------------------------------------------
  stime2 = cat_io_cmd('  Update Skull-Stripping','g5','',job.extopts.verb-1,stime2); 
  if size(P,4)==4 % skull-stripped
    [Yb,Ybb,Yg,Ydiv,P] = cat_main_updateSPM_skullstriped(Ysrc,P,res,vx_vol,T3th);
  elseif job.extopts.gcutstr==0 
    [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th);
  elseif job.extopts.gcutstr==2
    [Yb,Ybb,Yg,Ydiv] = cat_main_APRG(Ysrc,P,res,T3th);
  else
    [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcutold(Ysrc,P,res,vx_vol,T3th);
  end

  
  
  
  
  %% save brainmask using SPM12 segmentations for later use
  if ~exist('Ym0','var'),
    Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
  end
  Yb0 = (Ym0 > min(0.5,max(0.25, job.extopts.gcutstr))); clear Ym0
  Yb0 = cat_vol_morph(cat_vol_morph(Yb0,'lo'),'c');


  
  
%%
  stime2 = cat_io_cmd('  Update probability maps','g5','',job.extopts.verb-1,stime2);
  if ~(job.extopts.INV && any(sign(diff(T3th))==-1))
    %% Update probability maps
    disp('Update probabiliy maps');
    % background vs. head - important for noisy backgrounds such as in MT weighting
    if size(P,4)==4 % skull-stripped
      Ybg = ~Yb;
    else
      if sum(sum(sum(P(:,:,:,6)>240 & Ysrc<cat_stat_nanmean(T3th(1:2)))))>10000
        Ybg = P(:,:,:,6); 
        [Ybgr,Ysrcr,resT2] = cat_vol_resize({Ybg,Ysrc},'reduceV',vx_vol,2,32); 
        Ybgrth = max(cat_stat_nanmean(Ysrcr(Ybgr(:)>128)) + 2*std(Ysrcr(Ybgr(:)>128)),T3th(1));
        Ybgr = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ybgr>128,'d') & Ysrcr<Ybgrth,'lo',1),'lc',1);
        Ybg  = cat_vol_resize(cat_vol_smooth3X(Ybgr,1),'dereduceV',resT2); 
        clear Ysrcr Ybgr; 
      else
        Ybg = ~Yb;
      end
    end
    P4   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg==0) + single(P(:,:,:,4)) .* (Ybg<1) ); % remove air in head
    P5   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg==0) + single(P(:,:,:,5)) .* (Ybg<1) ); % remove air in head
    P6   = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg==1) + single(P(:,:,:,6)) .* (Ybg>0) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6;

    
    %% correct probability maps to 100% 
    disp('correct probability maps to 100%');
    sumP = cat_vol_ctype(255 - sum(P(:,:,:,1:6),4));
    P(:,:,:,1) = P(:,:,:,1) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>cat_stat_nanmean(T3th(1:2)) & Ysrc<cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,2) = P(:,:,:,2) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>=cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,3) = P(:,:,:,3) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc<=cat_stat_nanmean(T3th(1:2)));
    P(:,:,:,4) = P(:,:,:,4) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc<T3th(2));
    P(:,:,:,5) = P(:,:,:,5) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc>=T3th(2));
    P(:,:,:,6) = P(:,:,:,6) + sumP .* uint8( Ybg>=0.5 & ~Yb );
    clear Ybg sumP;

    
    %% head to WM 
    % Undercorrection of strong inhomogeneities in high field scans 
    % (>1.5T) can cause missalignments of the template and therefore 
    % miss classifications of the tissues that finally avoid further 
    % corrections in by LAS. 
    % Typically the alginment failed in this cases because the high 
    % intensities next to the head that were counted as head and not
    % corrected by SPM.
    % e.g. HR075, Magdeburg7T, SRS_SRS_Jena_DaRo81_T1_20150320-191509_MPR-08mm-G2-bw330-nbc.nii, ...
    disp('head to WM correction');
    Ywm = single(P(:,:,:,2)>128 & Yg<0.3 & Ydiv<0.03); Ywm(Ybb<128 | (P(:,:,:,1)>128 & abs(Ysrc/T3th(3)-2/3)<1/3) | Ydiv>0.03) = nan;
    [Ywm1,YD] = cat_vol_downcut(Ywm,1-Ysrc/T3th(3),0.02); Yb(isnan(Yb))=0; Ywm(YD<300)=1; Ywm(isnan(Ywm))=0; clear Ywm1 YD; %#ok<ASGLU>
    Ywmc = uint8(smooth3(Ywm)>0.7);
    Ygmc = uint8(cat_vol_morph(Ywmc,'d',2) & ~Ywmc & Ydiv>0 & Yb & cat_vol_smooth3X(Yb,8)<0.9);
    P(:,:,:,[1,3:6]) = P(:,:,:,[1,3:6]) .* repmat(1-Ywmc,[1,1,1,5]);
    P(:,:,:,2:6)     = P(:,:,:,2:6)     .* repmat(1-Ygmc,[1,1,1,5]);
    P(:,:,:,1)       = max(P(:,:,:,1),255*Ygmc);
    P(:,:,:,2)       = max(P(:,:,:,2),255*Ywmc);
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ygmc Ywmc; 

    
    %% head to GM ... important for children
    disp('head to GM correction');
    [Ywmr,Ybr,resT2] = cat_vol_resize({Ywm,Yb},'reduceV',vx_vol,2,32); 
    Ygm = cat_vol_morph(Ywmr>0.5,'d',3) & (cat_vol_morph(~Ybr,'d',3) | cat_vol_morph(Ybr,'d',1)); clear Ybr Ywmr;  % close to the head
    Ygm = cat_vol_resize(single(Ygm),'dereduceV',resT2)>0.5;
    Ygm = Ygm & Yp0<2/3 & Yb & Yg<cat_stat_nanmean(Yg(P(:,:,:,1)>64)) & Ydiv<cat_stat_nanmean(Ydiv(P(:,:,:,1)>64)); % add GM with low SPM prob ... 
    Ygm = Ygm & (Ysrc>cat_stat_nansum(T3th(1:2).*[0.5 0.5])) & (Ysrc<cat_stat_nansum(T3th(2:3).*[0.2 0.8])); % but good intensity
    Ygm(smooth3(Ygm)<0.5)=0; 
    clear Yg Ydiv;
    Ygm = uint8(Ygm); 
    P(:,:,:,5) = P(:,:,:,5) .* (1-Ygm);
    P(:,:,:,3) = P(:,:,:,3) .* (1-Ygm);
    P(:,:,:,2) = P(:,:,:,2) .* (1-Ygm);
    P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) + 255*single(Ygm));
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ywm Ygm;

    
    %% remove brain tissues outside the brainmask ...
    disp('remove brain tissues outside the brainmask ...');
    %  tissues > skull (within the brainmask)
    Yhdc = uint8(smooth3( Ysrc/T3th(3).*(Ybb>cat_vol_ctype(0.2*255)) - Yp0 )>0.5); 
    sumP = sum(P(:,:,:,1:3),4); 
    P(:,:,:,4)   =  cat_vol_ctype( single(P(:,:,:,4)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc<T3th(2)));
    P(:,:,:,5)   =  cat_vol_ctype( single(P(:,:,:,5)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc>=T3th(2)));
    P(:,:,:,1:3) =  P(:,:,:,1:3) .* repmat(uint8(~(Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ),[1,1,1,3]);
    clear sumP Yp0 Yhdc; 
  end
  clear Ybb;
  
  
  

  %% MRF
  % Used spm_mrf help and tested the probability TPM map for Q without good results. 
  stime = cat_io_cmd('spm_MRF');
  nmrf_its = 0; % 10 iterations better to get full probability in thin GM areas 
  spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
  fprintf('%s\n','MRF iterations:   '); % ExploreASL hack
  if isfield(res,'mg'), Kb = max(res.lkp); else Kb = size(res.intensity(1).lik,2); end
  G   = ones([Kb,1],'single');
  vx2 = single(sum(res.image(1).mat(1:3,1:3).^2));
  % P = zeros([d(1:3),Kb],'uint8');
  % P = spm_mrf(P,Q,G,vx2); % init: transfer data from Q to P 
  if 0
    %% use TPM as Q
    Q = zeros(size(P),'uint8');
    for di=1:6
      xASL_TrackProgress(di,6); % ExploreASL hack for counting on screen  
      vol = cat_vol_ctype(spm_sample_vol(tpm.V(di),...
        double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
      Q(:,:,:,di) = reshape(vol,d);
    end
  end
  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
  for iter=1:nmrf_its,
      xASL_TrackProgress(iter,nmrf_its); % ExploreASL hack for counting on screen
      P = spm_mrf(P,single(P),G,vx2); % spm_mrf(P,Q,G,vx2);
      spm_progress_bar('set',iter);
  end

  % update segmentation for error report
  Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
  cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);
 
  spm_progress_bar('clear');
  for k1=1:size(P,4)
      xASL_TrackProgress(k1,size(P,4));
      Ycls{k1} = P(:,:,:,k1); %#ok<AGROW>
  end
  clear Q P q q1 Coef b cr N lkp n wp M k1


  if job.extopts.verb>2
    % save information for debugging and OS test
    % input variables + bias corrected, bias field, class image
    % strong differences in bias fields can be the result of different 
    % registration > check 'res.image.mat' and 'res.Affine'
    [pth,nam] = spm_fileparts(res.image0(1).fname); 
    tpmci  = 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postbias'));
    save(tmpmat,'res','tpm','job','Ysrc','Ybf','Ycls');
  end
 
  clear Ybf
  stime2 = cat_io_cmd(' ','g5','',job.extopts.verb-1,stime2); 
  fprintf('%5.0fs\n',etime(clock,stime));



end
function [Yb,Ybb,Yg,Ydiv,P] = cat_main_updateSPM_skullstriped(Ysrc,P,res,vx_vol,T3th)
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    Yb   = Yp0>=0.5/3 & Ysrc>0; 
    Ybb  = cat_vol_ctype(Yb)*255; 

    P(:,:,:,6) = P(:,:,:,4); 
    P(:,:,:,4) = zeros(size(Ysrc),'uint8');
    P(:,:,:,5) = zeros(size(Ysrc),'uint8'); 
    res.lkp = [res.lkp 5 6];
    res.mn  = [res.mn(1:end-1),0,0,0];
    res.mg  = [res.mg(1:end-1);1;1;1];
    res.vr(1,1,numel(res.lkp)-1:numel(res.lkp)) = 0;
     
    [Ysrcb,BB] = cat_vol_resize(Ysrc,'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3); clear Yp0;
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
end
function [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th)
    % brain mask
    disp('Resizing brain mask');
    Ym   = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
    Yb   = (Ym > 0.5);
    Yb   = cat_vol_morph(cat_vol_morph(Yb,'lo'),'c');
    Ybb  = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*256); 

    [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
end
function [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcutold(Ysrc,P,res,vx_vol,T3th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1 only > remove in future if gcut is removed too!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
    voli   = @(v) (v ./ (pi * 4./3)).^(1/3);   % volume > radius
    Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));

    % old skull-stripping
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    brad = voli(sum(Yp0(:)>0.5).*prod(vx_vol)/1000); 
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    %Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
    Yg   = cat_vol_grad((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ydiv = cat_vol_div((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ybo  = cat_vol_morph(cat_vol_morph(Yp0>0.3,'lc',2),'d',brad/2/mean(vx_vol)); 
    BVth = diff(T3th(1:2:3))/abs(T3th(3))*1.5; 
    RGth = diff(T3th(2:3))/abs(T3th(3))*0.1; 
    Yb   = single(cat_vol_morph((Yp0>1.9/3) | (Ybo & Ysrcb>mean(T3th(2)) & ...
           Ysrcb<T3th(3)*1.5 & Yg<0.5),'lo',max(0,0.6/mean(vx_vol)))); 
    
    %% region-growing GM 1
    Yb(~Yb & (~Ybo | Ysrcb<cat_stat_nanmean(T3th(2)) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU> 
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    
    %% region-growing GM 2
    Yb(~Yb & (~Ybo | Ysrcb<T3th(1) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/2); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU>
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    
    %% region-growing GM 3
    Yb(~Yb & (~Ybo | Ysrcb<mean([BGth,T3th(1)]) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan; clear Ybo;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/10); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU>
    Yb(smooth3(Yb)<0.5)=0; Yb(Yp0toC(Yp0*3,1)>0.9 & Yg<0.3 & Ysrcb>BGth & Ysrcb<T3th(2)) = 1; 
    
    %% ventricle closing
    [Ybr,Ymr,resT2] = cat_vol_resize({Yb>0,Ysrcb/T3th(3)},'reduceV',vx_vol,2,32); clear Ysrcb
    Ybr = Ybr | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6)); clear Ymr;  % large ventricle closing
    Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
    Yb  = Yb | cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.7; clear Ybr
    Yb  = smooth3(Yb)>0.5; 
    Ybb = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*255); 
    Yb   = cat_vol_resize(Yb   , 'dereduceBrain' , BB);
    Ybb  = cat_vol_resize(Ybb  , 'dereduceBrain' , BB);
    Yg   = cat_vol_resize(Yg   , 'dereduceBrain' , BB);
    Ydiv = cat_vol_resize(Ydiv , 'dereduceBrain' , BB);
    clear Ysrcb Ybo;
end
