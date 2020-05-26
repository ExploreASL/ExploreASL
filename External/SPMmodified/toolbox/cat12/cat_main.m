function Ycls = cat_main(res,tpm,job)
% Write out CAT preprocessed data
%
% FORMAT Ycls = cat_main(res,tpm,job)
%
% based on John Ashburners version of
% spm_preproc_write8.m 2531 2008-12-05 18:59:26Z john $
%
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_main.m 1361 2018-08-31 08:44:57Z gaser $

%#ok<*ASGLU>

% if there is a breakpoint in this file set debug=1 and do not clear temporary variables
dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

% this limits ultra high resolution data, i.e. images below ~0.4 mm are reduced to ~0.7mm!
% used in cat_main_partvol, cat_main_gcut, cat_main_LAS
def.extopts.uhrlim = 0.7 * 2; % default 0.7*2 that reduce images below 0.7 mm
def.extopts.xasl_quality = 1; % default for xasl_quality == 1
job.extopts = cat_io_checkinopt(job.extopts,def.extopts);
clear def;

debug=1; %%% ExploreASL fix
global cat_err_res; % for CAT error report

tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;

% EXPLOREASL HACK
CurrDir=pwd;
if exist('./mri','dir')
    cd('mri');
end
if xASL_exist('nT1.nii')
    tIM = xASL_io_Nifti2Im('nT1.nii');
    tIM(isnan(tIM)) = 0;
    xASL_io_SaveNifti('nT1.nii', 'nT1.nii', tIM, [], false);
end
cd(CurrDir);

%% complete job structure
defr.ppe = struct();
res = cat_io_checkinopt(res,defr);

cat_io_cmd('Initiate CAT12 main settings'); % ExploreASL fix

% Deprecated, now using the field job.extopts.xasl_quality
%Fpath   = xASL_fileparts(char(job.data)); % ExploreASL hack
%if  exist(fullfile(Fpath,'LowQ.status'),'file')
%    LowQ = true;
%else
%    LowQ = false;
%end


def.cati            = 0;
def.color.error     = [0.8 0.0 0.0];
def.color.warning   = [0.0 0.0 1.0];
def.color.warning   = [0.8 0.9 0.3];
def.color.highlight = [0.2 0.2 0.8];
job = cat_io_checkinopt(job,def);

% definition of subfolders
if job.extopts.subfolders
  mrifolder     = 'mri';
  reportfolder  = 'report';
  labelfolder   = 'label';
else
  mrifolder     = '';
  reportfolder  = '';
  labelfolder   = '';
end

M1   = tpm.M;

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
bb = job.extopts.bb;
vx = job.extopts.vox(1);
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
bb(1,:) = vx.*round(bb(1,:)./vx);
bb(2,:) = vx.*round(bb(2,:)./vx);
res.bb = bb;
clear vx vx1 bb1

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end

N = numel(res.image);
if N > 1
  warning('CAT12:noMultiChannel',...
    'CAT12 does not support multiple channels. Only the first channel will be used.');
end

do_dartel = 1 + (job.extopts.regstr(1)~=0);      % always use dartel (do_dartel=1) or shooting (do_dartel=2) normalization
if do_dartel
  need_dartel = any(job.output.warps) || ...
    job.output.bias.warped || ...
    job.output.label.warped || ...
    any(any(tc(:,[4 5 6]))) || job.output.jacobian.warped || ...
    job.output.ROI || ...
    any([job.output.atlas.warped]) || ...
    numel(job.extopts.regstr)>1 || ...
    numel(job.extopts.vox)>1;
  if ~need_dartel
    do_dartel = 0;
  end
end
if do_dartel<2, job.extopts.templates = job.extopts.darteltpms; else job.extopts.templates = job.extopts.shootingtpms; end % for LAS

%%% ExploreASL fix
if job.extopts.xasl_disabledartel % DARTEL disabling request
	do_dartel = 0;
	fprintf('%s\n','DARTEL disabled');
else
	do_dartel = 1;
	fprintf('%s\n','DARTEL enabled');
end

res.do_dartel = do_dartel;


stime = cat_io_cmd('SPM preprocessing 2 (write)');
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');%%% ExploreASL fix
if ~isfield(res,'spmpp')
  if job.extopts.verb>1, fprintf('\n'); end
  stime2 = cat_io_cmd('  Write Segmentation','g5','',job.extopts.verb-1);
end
[Ysrc,Ycls,Yy] = cat_spm_preproc_write8(res,zeros(max(res.lkp),4),zeros(1,2),[0 0],0,0);
if isfield(res,'redspmres')
  % Update Ycls: cleanup on original data
  Yb = Ycls{1} + Ycls{2} + Ycls{3};
  for i=1:numel(Ycls), [Pc(:,:,:,i),BB] = cat_vol_resize(Ycls{i},'reduceBrain',repmat(job.opts.redspmres,1,3),2,Yb); end %#ok<AGROW>
    Pc = clean_gwc(Pc,1);
  for i=1:numel(Ycls), Ycls{i} = cat_vol_resize(Pc(:,:,:,i),'dereduceBrain',BB); end; clear Pc Yb;
  for ci=1:numel(Ycls)
    Ycls{ci} = cat_vol_ctype(cat_vol_resize(Ycls{ci},'deinterp',res.redspmres,'linear'));
  end

  % Update Yy:
  Yy2 = zeros([res.redspmres.sizeO 3],'single');
  for ci=1:size(Yy,4)
    Yy2(:,:,:,ci) = cat_vol_ctype(cat_vol_resize(Yy(:,:,:,ci),'deinterp',res.redspmres,'linear'));
  end
  Yy   = Yy2; clear Yy2;

  % Update Ysrc:
  Ysrc = cat_vol_resize(Ysrc,'deinterp',res.redspmres,'cubic');
  Ybf  = res.image1.dat ./ Ysrc;
  Ybf  = cat_vol_approx(Ybf .* (Ysrc~=0 & Ybf>0.25 & Ybf<1.5),'nn',1,8);
  Ysrc = res.image1.dat ./ Ybf; clear Ybf;
  res.image = res.image1;
  res  = rmfield(res,'image1');
end

% remove noise/interpolation prefix
VT  = res.image(1);  % denoised/interpolated n*.nii
VT0 = res.image0(1); % original
fname0 = VT0.fname;
[pth,nam] = spm_fileparts(VT0.fname);

% voxel size parameter
vx_vol  = sqrt(sum(VT.mat(1:3,1:3).^2));    % voxel size of the processed image
vx_volr = sqrt(sum(VT0.mat(1:3,1:3).^2));   % voxel size of the original image
vx_volp = prod(vx_vol)/1000;
voli    = @(v) (v ./ (pi * 4./3)).^(1/3);   % volume > radius

% delete old xml file
oldxml = fullfile(pth,reportfolder,['cat_' nam '.xml']);
if exist(oldxml,'file'), delete(oldxml); end
clear oldxml

d = VT.dim(1:3);
P = zeros([size(Ycls{1}) numel(Ycls)],'uint8');
for i=1:numel(Ycls), P(:,:,:,i) = Ycls{i}; end
clear Ycls;

if ~isfield(res,'spmpp')
  stime2 = cat_io_cmd('  Update Segmentation','g5','',job.extopts.verb-1,stime2);
  % cleanup with brain mask - required for ngaus [1 1 2 4 3 2] and R1/MP2Rage like data

  % tpm brain mask
  YbA = false(d);
  Vb = tpm.V(1); Vb.pinfo(3) = 0; Vb.dat = single(exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3}));

  for z=1:d(3)
    YbA(:,:,z) = spm_sample_vol(Vb,double(Yy(:,:,z,1)),double(Yy(:,:,z,2)),double(Yy(:,:,z,3)),1)>0.1;
  end

  % add some distance around brainmask (important for bias!)
  YbA = YbA | cat_vol_morph(YbA & sum(P(:,:,:,1:2),4)>4 ,'dd',2.4,vx_vol);
  % transfer tissue outside the brain mask to head  ...
  % RD 201807: I am not sure if this is a good idea. Please test this with children!
  for i=1:3
    P(:,:,:,4) = P(:,:,:,4) + P(:,:,:,i) .* uint8(~YbA);
    P(:,:,:,i) = P(:,:,:,i) .* uint8(YbA);
  end
  clear YbA;

  % cleanup for high resolution data
  % Alghough the old cleanup is very slow for high resolution data, the
  % reduction of image resolution removes spatial segmenation information.
  if job.opts.redspmres==0 % already done in case of redspmres
    if max(vx_vol)<1.5 && mean(vx_vol)<1.3
		fprintf('Cleanup for high resolution data:   ');
      for i=1:size(P,4), [Pc1(:,:,:,i),RR] = cat_vol_resize(P(:,:,:,i)  ,'reduceV',vx_vol,job.extopts.uhrlim,32); end %#ok<AGROW>
      for i=1:size(P,4), [Pc2(:,:,:,i),BB] = cat_vol_resize(Pc1(:,:,:,i),'reduceBrain',vx_vol,2,sum(Pc1,4)); end %#ok<AGROW>
      Pc2 = clean_gwc(Pc2,1);
      for i=1:size(P,4), Pc1(:,:,:,i) = cat_vol_resize(Pc2(:,:,:,i),'dereduceBrain',BB); end
      for i=1:size(P,4), P(:,:,:,i)   = cat_vol_resize(Pc1(:,:,:,i),'dereduceV',RR); end
      clear Pc1 Pc2;
    end
  end

  % garantee probability
  sP = (sum(single(P),4)+eps)/255;
  for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
  clear sP;

  % median in case of WMHs!
  WMth = double(max( clsint(2) , cat_stat_nanmedian(Ysrc(P(:,:,:,2)>192)) ));
  if clsint(3)>clsint(2) % invers
    CMth = clsint(3);
  else
    CMth = min( [  clsint(1) - diff([clsint(1),WMth]) , clsint(3) ]);
  end
  T3th = [ CMth , clsint(1) , WMth];


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
          'Possibly the affine registration failed. Pleace check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM','HDH','HDL', ...
          [ clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol)*100, clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th,HDHth,HDLth]));  %#ok<SPERR>
    elseif size(P,4)==4 % skull-stripped
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Pleace check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM', ...
          [ clsvol([4 3 1 2])/cat_stat_nansum(clsvol)*100, clsvol([4 3 1 2])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th]));  %#ok<SPERR>
    else
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ['Empty Segmentation: ' ...
           'Possibly the affine registration failed. Pleace check image orientation.\n']);
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

  %%
 stime2 = cat_io_cmd('  Update Skull-Stripping','g5','',job.extopts.verb-1,stime2);
  skullstripped = max(res.lkp) == 4;

  if skullstripped % skull-stripped
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    Yb   = Yp0>=0.5/3;
    Ybb  = cat_vol_ctype(Yb)*255;

    P(:,:,:,6) = P(:,:,:,4);
    P(:,:,:,4) = zeros(size(Yp0),'uint8');
    P(:,:,:,5) = zeros(size(Yp0),'uint8');
    res.lkp = [res.lkp 5 6];
    res.mn  = [res.mn(1:end-1),0,0,0];
    res.mg  = [res.mg(1:end-1);1;1;1];
    res.vr(1,1,numel(res.lkp)-1:numel(res.lkp)) = 0;

    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  elseif (job.extopts.experimental || (job.extopts.INV && any(sign(diff(T3th))==-1))) && job.extopts.gcutstr>0
   % (sum( abs( (Ysrc(:)==0) - (Yp0(:)<0.5) ) ) / sum(Ysrc(:)==0)) < 0.1  || ...
  % use gcut2

    %   brad = voli(sum(Yp0(:)>0).*prod(vx_vol)/1000);
    %   noise = nanstd(Ysrc(cat_vol_morph(Yp0>0.8,'o') & Yp0>0.99)/diff(T3th(2:3)));

    Vl1 = spm_vol(job.extopts.cat12atlas{1});
    Yl1 = cat_vol_ctype(spm_sample_vol(Vl1,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
    Yl1 = reshape(Yl1,size(Ysrc)); [D,I] = cat_vbdist(single(Yl1>0)); Yl1 = Yl1(I);

    %%
    clear D I
    Ybg = P(:,:,:,6)>128;
    T3ths = [min(min(min(single(Ysrc(P(:,:,:,6)>128))))),...
             min( cat_stat_nanmean(Ysrc(Ybg(:))) + 2*cat_stat_nanstd(Ysrc(Ybg(:))) , ...
              mean([cat_stat_nanmean(Ysrc(Ybg(:))),min(T3th)])), ...
             T3th, T3th(3) + cat_stat_nanmean(diff(T3th))]; clear Ybg;
    T3thx = [0,0.05,1,2,3,4];
    if T3th(1)>T3th(3), T3thx = [0,0.05,1,2,3,2]; T3ths(end) = T3ths(2); end;
    [T3ths,si] = sort(T3ths);
    T3thx      = T3thx(si);
    Ym = Ysrc+0;
    for i=numel(T3ths):-1:2
      M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3ths(end);
    Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));
    Ym = Ym / 3;
    clear M;

    for k1=1:3, Ycls{k1} = P(:,:,:,k1); end
    Yb = cat_main_gcut(Ym,Yp0>0.1,Ycls,Yl1,false(size(Ym)),vx_vol,...
      struct('gcutstr',0.1,'verb',0,'LAB',job.extopts.LAB,'LASstr',0,'red',1));
    clear Ycls;
    Ybb  = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*255);

    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);

    Yb   = smooth3(Yb)>0.5;
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  elseif job.extopts.gcutstr==0 || job.extopts.gcutstr==2 || job.extopts.gcutstr==3
    % brain mask
	disp('Resizing brain mask');
    Ym  = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
    Yb   = (Ym > 0.5);
    Yb = cat_vol_morph(cat_vol_morph(Yb,'lo'),'c');
    Ybb  = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*256);

    [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  else
    % old skull-stripping
    brad = voli(sum(Yp0(:)>0.5).*prod(vx_vol)/1000);
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    %Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
    Yg   = cat_vol_grad((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ydiv = cat_vol_div((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ybo  = cat_vol_morph(cat_vol_morph(Yp0>0.3,'lc',2),'d',brad/2/mean(vx_vol));
    BVth = diff(T3th(1:2:3))/abs(T3th(3))*1.5;
    RGth = diff(T3th(2:3))/abs(T3th(3))*0.1;
    Yb   = single(cat_vol_morph((Yp0>1.9/3) | (Ybo & Ysrcb>mean(T3th(2)) & Ysrcb<T3th(3)*1.5 & Yg<0.5),'lo',max(0,0.6/mean(vx_vol))));
    %% region-growing GM 1
    Yb(~Yb & (~Ybo | Ysrcb<cat_stat_nanmean(T3th(2)) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD;
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    %% region-growing GM 2
    Yb(~Yb & (~Ybo | Ysrcb<T3th(1) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/2); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD;
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    %% region-growing GM 3
    Yb(~Yb & (~Ybo | Ysrcb<mean([BGth,T3th(1)]) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan; clear Ybo;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/10); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD;
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
    clear Ybo;
  end

  clear Ysrcb Yp0;

  % save brainmask using SPM12 segmentations for later use
  Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
  Yb0 = (Ym0 > min(0.5,max(0.25, job.extopts.gcutstr))); clear Ym0
  Yb0 = cat_vol_morph(cat_vol_morph(Yb0,'lo'),'c');

  %%
  stime2 = cat_io_cmd('  Update probability maps','g5','',job.extopts.verb-1,stime2);
  if ~(job.extopts.INV && any(sign(diff(T3th))==-1))
    disp('Update probabiliy maps');
	%% Update probability maps
    % background vs. head - important for noisy backgrounds such as in MT weighting
    if job.extopts.gcutstr==0
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
	disp('head to WM correction');
    % e.g. HR075, Magdeburg7T, SRS_SRS_Jena_DaRo81_T1_20150320-191509_MPR-08mm-G2-bw330-nbc.nii, ...
    Ywm = single(P(:,:,:,2)>128 & Yg<0.3 & Ydiv<0.03); Ywm(Ybb<128 | (P(:,:,:,1)>128 & abs(Ysrc/T3th(3)-2/3)<1/3) | Ydiv>0.03) = nan;
    [Ywm1,YD] = cat_vol_downcut(Ywm,1-Ysrc/T3th(3),0.02); Yb(isnan(Yb))=0; Ywm(YD<300)=1; Ywm(isnan(Ywm))=0; clear Ywm1 YD;
    Ywmc = uint8(smooth3(Ywm)>0.7);
    %%
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
    %%
    Ygm = uint8(Ygm);
    P(:,:,:,5) = P(:,:,:,5) .* (1-Ygm);
    P(:,:,:,3) = P(:,:,:,3) .* (1-Ygm);
    P(:,:,:,2) = P(:,:,:,2) .* (1-Ygm);
    P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) + 255*single(Ygm));
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ywm Ygm;

    %% remove brain tissues outside the brainmask ...
	disp('remove brain tissues outside the brainmask ...');
    % tissues > skull (within the brainmask)
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
  nmrf_its = 0; % 10 interations better to get full probability in thin GM areas
  spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
  fprintf('%s\n','MRF iterations:   '); % ExploreASL hack
  G   = ones([Kb,1],'single');
  vx2 = single(sum(VT.mat(1:3,1:3).^2));
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

  spm_progress_bar('clear');
  for k1=1:size(P,4)
	  xASL_TrackProgress(k1,size(P,4));
      Ycls{k1} = P(:,:,:,k1);
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



  %% ---------------------------------------------------------------------
  %  check if the brain hits a image boundary
  %  ---------------------------------------------------------------------

  %Yib1 = cat_vol_morph(true(size(Yb)),'e',1);
  %Yib3 = cat_vol_morph(true(size(Yb)),'e',2);
  %Ybbs = sum(Yb & Yib);



  %% ---------------------------------------------------------------------
  %  Global (and local) intensity normalization and partioning
  %  ---------------------------------------------------------------------
  %  Global and local intensity corrections are the basis of most of the
  %  following functions. The global normalization based on the SPM tissue
  %  thresholds (res.mn) and were used anyway. For strong differences
  %  (mostly by the CSF) the median will used, because it is e.g. more
  %  stable. This will cause a warning by the cat_main_gintnorm.
  %
  %  The local adaptive segmentation include a further bias correction
  %  and a global  and local intensity correction. The local intensity
  %  correction refines the tissue maps to aproximate the local tissue
  %  peaks of WM (maximum-based), GM, and CSF.
  %
  %  If you want to see intermediate steps of the processing use the "ds"
  %  function:
  %    ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,80)
  %  that display 4 images (unterlay, overlay, image1, image2) for one
  %  slice. The images were scaled in a range of 0 to 1. The overlay
  %  allows up to 20 colors
  %
  %  ---------------------------------------------------------------------
  stime = cat_io_cmd('Global intensity correction');
  if debug
    Ym   = Ysrc / T3th(3); %#ok<NASGU> % only WM scaling
    Yp0  = (single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255)/3; %#ok<NASGU> % label map
  end
  if any( min(vx_vol*2,1.4)./vx_vol >= 2 )
    %%
    Ysrcr = cat_vol_resize(Ysrc,'reduceV',vx_vol,min(vx_vol*2,1.4),32,'meanm');
    Ybr   = cat_vol_resize(single(Yb),'reduceV',vx_vol,min(vx_vol*2,1.4),32,'meanm')>0.5;
    Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,min(vx_vol*2,1.4),32); end
    [Ymr,Ybr,T3th,Tth,job.inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrcr,Yclsr,Ybr,vx_vol,res,Yy,job.extopts);
    clear Ymr Ybr Ysrcr Yclsr;
    Ym = cat_main_gintnorm(Ysrc,Tth);
  else
    [Ym,Yb,T3th,Tth,job.inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res,Yy,job.extopts);;
  end

  % update in inverse case ... required for LAS
  if job.inv_weighting
%    Ysrc = Ym * Tth.T3th(5); Tth.T3th = Tth.T3thx * Tth.T3th(5);
    if T3th(1)>T3th(3) && T3th(2)<T3th(3) && T3th(1)>T3th(2)
      Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
      %ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,90)
      Yb2  = cat_vol_morph(Yp0>0.5,'lc',2);
      prob = cat(4,cat_vol_ctype(Yb2.*Yp0toC(Ym*3,2)*255),...
                   cat_vol_ctype(Yb2.*Yp0toC(Ym*3,3)*255),...
                   cat_vol_ctype(Yb2.*Yp0toC(min(3,Ym*3),1)*255));

      for i=1:size(prob,4), [prob2(:,:,:,i),BB] = cat_vol_resize(prob(:,:,:,i),'reduceBrain',vx_vol,2,sum(prob,4)); end %#ok<AGROW>
      prob2 = clean_gwc(prob2,1);
      for i=1:size(prob,4), prob(:,:,:,i) = cat_vol_resize(prob2(:,:,:,i),'dereduceBrain',BB); end; clear prob2;

      for ci=1:3, Ycls{ci} = prob(:,:,:,ci); end;
      %job.extopts.mrf = 0.3;
      clear prob;
    end

    Ysrc = Ym; Tth.T3thx(3:5) = 1/3:1/3:1; Tth.T3th = Tth.T3thx; T3th = 1/3:1/3:1;

  end
  fprintf('%5.0fs\n',etime(clock,stime));

  if job.extopts.verb>2
    tpmci  = tpmci + 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postgintnorm'));
    save(tmpmat,'Ysrc','Ycls','Ym','Yb','T3th','vx_vol');
  end


  %% After the intensity scaling and with correct information about the
  % variance of the tissue, a further harder noise correction is meaningful.
  % Finally, a stronger NLM-filter is better than a strong MRF filter!

  if job.extopts.NCstr~=0
    NCstr.labels = {'none','full','light','medium','strong','heavy'};
    NCstr.values = {0 1 2 -inf 4 5};
    stime = cat_io_cmd(sprintf('SANLM denoising after intensity normalization (%s)',...
      NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')}));

    % filter only within the brain mask for speed up
    [Yms,Ybr,BB] = cat_vol_resize({Ym,Yb},'reduceBrain',vx_vol,round(2/cat_stat_nanmean(vx_vol)),Yb); Ybr = Ybr>0.5;
    Yms = cat_vol_sanlm(struct('data',res.image0.fname,'verb',0,'NCstr',job.extopts.NCstr),res.image,1,Yms);
    Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = Yms .* Ybr + ...
      Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr);
    clear Yms Ybr BB;

    if job.inv_weighting
      Ysrc = Ym;
    else
      Ysrc = cat_main_gintnormi(Ym,Tth);
    end

    fprintf('%5.0fs\n',etime(clock,stime));

  end




  %% prepared for improved partitioning - RD20170320, RD20180416
  %  Update the initial SPM normalization by a fast version of Shooting
  %  to improve the skull-stripping, the partitioning and LAS.
  %  We need stong deformations in the ventricle for the partitioning
  %  but low deformations for the skull-stripping. Moreover, it has to
  %  be really fast > low resolution (3 mm) and less iterations.
  %  The mapping has to be done for the TPM resolution, but we have to
  %  use the Shooting template for mapping rather then the TPM because
  %  of the cat12 atlas map.
  if job.extopts.WMHC || job.extopts.SLC || job.extopts.gcutstr==2 || job.extopts.gcutstr==3
    stime = cat_io_cmd(sprintf('Fast Shooting registration'),'','',job.extopts.verb);

    job2 = job;

    if  job2.extopts.xasl_quality % ExploreASL hack fewer iterations to speed up
		job2.extopts.reg.nits = 16;     % more iterations
	else
        job2.extopts.reg.nits = 2;      % fewer iterations
    end

    job2.extopts.regstr   = 15;     % low resolution
    job2.extopts.verb     = debug;  % do not display process (people would may get confused)
    job2.extopts.vox      = abs(res.tpm(1).mat(1));  % TPM resolution to replace old Yy
    job2.extopts.shootingtpms(3:end) = [];             % remove high templates, we only need low frequency corrections
    res2 = res;
    res2.do_dartel        = 2;      % use shooting
    if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
      %[trans,res.ppe.reginitp] = cat_main_registration(job2,res2,Ycls(1:2),Yy,tpm.M,res.Ylesion);
	  [trans,res.ppe.reginitp] = cat_main_registration(job2,res2,Ycls(1:2),Yy,tpm.M,res.YlesionFull); % ExploreASL fix
    else
      [trans,res.ppe.reginitp] = cat_main_registration(job2,res2,Ycls(1:2),Yy,tpm.M);
    end
    %%
    Yy2  = trans.warped.y;
    %Ydtw = trans.jc.dt2;
    if ~debug, clear trans job2 res2; end



    % Shooting did not include areas outside of the boundary box
    % >> add to cat_main_registration?
    Ybd = true(size(Ym)); Ybd(3:end-2,3:end-2,3:end-2) = 0; Ybd(~isnan(Yy2(:,:,:,1))) = 0; Yy2(isnan(Yy2))=0;
    for k1=1:3
      Yy2(:,:,:,k1) = Yy(:,:,:,k1) .* Ybd + Yy2(:,:,:,k1) .* (1-Ybd);
      Yy2(:,:,:,k1) = cat_vol_approx(Yy2(:,:,:,k1),'nn',vx_vol,3);
    end
    if ~debug, Yy = Yy2; end

    %% preparte mapping in case of shooting templates with other resolutions
    Vb2    = spm_vol(job.extopts.shootingtpms{2});
    amat   = Vb2(2).mat \ M1;
    if any(any(amat~=eye(4)));
      yn     = numel(Yy2);
      p      = ones([4,yn/3],'single');
      p(1,:) = Yy2(1:yn/3);
      p(2,:) = Yy2(yn/3+1:yn/3*2);
      p(3,:) = Yy2(yn/3*2+1:yn);
      p      = amat * p;

      Yy2 = zeros(size(Yy2),'single');
      Yy2(1:yn/3)        = p(1,:);
      Yy2(yn/3+1:yn/3*2) = p(2,:);
      Yy2(yn/3*2+1:yn)   = p(3,:);
      clear p;
    end

    %% Used (flipped) Jacobian determinant to detect lesions as strongly deformed regions.
    % I mapped the normalized determinant to individual space because I was not
    % able to create the determinant from the inverse deformation field.
    % Furthermore this idea works only in single cases yet and is therefore not activated yet.
    % I will remove this block and related variables in October 2018 if it still not working.
    if 0
      Vdtw = Vb2(1); Vdtw.dt(1) = 16; Vdtw.pinfo(1) = 1; Vdtw.pinfo(3) = 0; Vdtw.dat = Ydtw;
      Vdtw.mat = trans.warped.M1; Vdtw.dim = size(Ydtw);
      if isfield(Vdtw,'private'), Vdtw = rmfield(Vdtw,'private'); end
      Ydt  = spm_sample_vol(Vdtw,double(Yy2(:,:,:,1)),double(Yy2(:,:,:,2)),double(Yy2(:,:,:,3)),5);
      Ydt  = reshape(Ydt,size(Ym));
      %%
      mati = spm_imatrix(Vdtw.mat); mati([1,7]) = -mati([1,7]); Vdtw.mat = spm_matrix(mati);
      Vdtw.dat = flipud(Ydtw);
      Ydti = spm_sample_vol(Vdtw,double(Yy2(:,:,:,1)),double(Yy2(:,:,:,2)),double(Yy2(:,:,:,3)),5);
      Ydti = reshape(Ydti,size(Ym));
    end

    %%
    clear Ydtw;
    if job.extopts.gcutstr==2 || job.extopts.gcutstr==3
    %% Update brain mask
    stime = cat_io_cmd(sprintf('SPM+ Skull-Stripping'),'','',job.extopts.verb,stime);

      % load template tissue maps
      Yg  = spm_sample_vol(Vb2(1),double(Yy2(:,:,:,1)),double(Yy2(:,:,:,2)),double(Yy2(:,:,:,3)),5);
      Yw  = spm_sample_vol(Vb2(2),double(Yy2(:,:,:,1)),double(Yy2(:,:,:,2)),double(Yy2(:,:,:,3)),5);
      if ~debug, clear Yy2; end
      Yg  = reshape(Yg,size(Ym));
      Yw  = reshape(Yw,size(Ym));

      if debug, Ybo2 = Yb; Yclso=Ycls; end

      %% update brainmask
      Yp0 = (single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255)/3;
      ith = 0.2; % lower values will result in wider mask
      Yb  = ( cat_vol_smooth3X(Yg + Yw,2/mean(vx_vol)) .* ... smooth brain template mask
              (Yp0toC(Ym*3,2) + Yp0toC(Ym*3,3) - Yp0toC(Ym*3,4)/2 - Yp0toC(Ym*3,1)/max(2,noise*10)) )>ith & ... individual intensities
              Ym>0.5*(Yg+Yw) & (Ym<0.8 | Ym<2*(Yg+Yw)) &  ... major boundaries
              Ym>0.3 & Ym<1.2 & cat_vol_smooth3X(Yp0,2/mean(vx_vol))>0.05/3; % GM and WM
      Yb = Yb & ~(Ym>1 & cat_vol_smooth3X(Yp0)<1.1/3);                % try to avoid previously removed head tissue
      Yb = Yb | Ycls{2}>128 | (Ycls{1}>64 & Yp0toC(Ym*3,max(1.5,min(2.2,2-noise*2)))>0);    % add save WM and GM
      Yb = ~cat_vol_morph( ~Yb, 'l', [2,0.5]);                        % close voxel wholes (missed high intensity WM voxels)
      Yb = cat_vol_morph( Yb, 'l', [2,0.5]);                          % remove small voxel (head tissue)
      Yb = Yb | (Ycls{1}/255>0.2 & cat_vol_morph( Yb, 'ldc', max(0.5,min(2,noise*10))));
      Yb = cat_vol_morph( Yb, 'lo', 1);                               % remove larger parts of head tissue
      Yb = Yb | smooth3( (Ycls{3} + Ycls{1})/255>0.4 & Ym>0.3/3 & Ym<1.5/3)>0.6;  % add CSF
      Yb = cat_vol_morph( cat_vol_morph( Yb, 'ldo', 2 , vx_vol), 'ldc', 3, vx_vol); % remove large wholes and head tissue
      % very small dilation to get PVE voxel (better to include menignes than to miss GM)
      Yb = Yb | (cat_vol_morph(Yb,'dd',min(2,1 + noise*10)) & Ym<3/3 & Ym>0.9/3);
      Yb = Yb | (cat_vol_morph(Yb,'dd',min(2,1 + noise*10)) & Ym<2/3 & Ym>0.7/3);
      Yb = smooth3(Yb)>0.5; % make it nice
      if ~debug, clear Yg Yw; end

      %% update segmentation (add GM)
      Ygd = Yb .* min(1,Yp0toC(Ym*3,2)*1.2) - single(Ycls{1})/255;
      Ygd(smooth3(Ygd)<0.2 | Yp0>2/3) = 0;
      Ycls{1} = Ycls{1} + cat_vol_ctype(Ygd * 255);
      Ycls{3} = Ycls{3} - cat_vol_ctype(Ygd * 255);
      if ~debug, clear Yp0 Ygd; end
      Yclss = zeros(size(Ym),'uint8');
      for ci=1:numel(Ycls), if ~isempty(Ycls{ci}), Yclss = Yclss + Ycls{ci}; end; end
      for ci=1:numel(Ycls), if ~isempty(Ycls{ci}), Yclss = cat_vol_ctype(single(Ycls{ci}) ./ single(Yclss) * 255); end; end
      if ~debug, clear Yclss; end
    end
    fprintf('%5.0fs\n',etime(clock,stime));
  end



  %% Local Intensity Correction
  Ymo = Ym;
  if job.extopts.LASstr>0
    if job.extopts.LASstr>1
      extoptsLAS2 = job.extopts;
      extoptsLAS2.LASstr = extoptsLAS2.LASstr-1;
      stime = cat_io_cmd(sprintf('Local adaptive segmentation 2 (LASstr=%0.2f)',extoptsLAS2.LASstr));
      [Ymi,Ym,Ycls] = cat_main_LASs(Ysrc,Ycls,Ym,Yb,Yy,Tth,res,vx_vol,extoptsLAS2); % use Yclsi after cat_vol_partvol
    else
      stime = cat_io_cmd(sprintf('Local adaptive segmentation (LASstr=%0.2f)',job.extopts.LASstr));
      [Ymi,Ym,Ycls] = cat_main_LAS(Ysrc,Ycls,Ym,Yb,Yy,T3th,res,vx_vol,job.extopts,Tth);
    end
    fprintf('%5.0fs\n',etime(clock,stime));

    if job.extopts.NCstr~=0
      %% noise correction of the local normalized image Ymi, whereas only small changes are expected in Ym by the WM bias correction
      stime = cat_io_cmd(sprintf('  SANLM denoising after LAS (%s)',...
        NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')}),'g5');

      [Ymis,Ymior,BB]  = cat_vol_resize({Ymi,Ymo},'reduceBrain',vx_vol,round(2/mean(vx_vol)),Yb);
      Ymis = cat_vol_sanlm(struct('data',res.image0.fname,'verb',0,'NCstr',job.extopts.NCstr),res.image,1,Ymis);

      Yc = abs(Ymis - Ymior); Yc = Yc * 6 * min(2,max(0,abs(job.extopts.NCstr)));
      spm_smooth(Yc,Yc,2./vx_vol); Yc = max(0,min(1,Yc)); clear Ymior;
      % mix original and noise corrected image and go back to original resolution
      Ybr = Yb(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6));
      Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = ...
        Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr) + ...
        (1-Yc) .* Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* Ybr + ...
        Yc .* Ymis .* Ybr;

      % extreme background denoising to remove wholes?
      Ymis = cat_vol_median3(Ymi,Ymi>0 & Ymi<0.4,Ymi<0.4); Ymi = Ymi.*max(0.1,Ymi>0.4) + Ymis.*min(0.9,Ymi<=0.4);
      Ymis = cat_vol_median3(Ym,Ym>0 & Ym<0.4,Ym<0.4); Ym = Ym.*max(0.1,Ym>0.4) + Ymis.*min(0.9,Ym<=0.4);

      clear Ymis;
    end

    cat_io_cmd(' ','','',job.extopts.verb,stime);
    fprintf('%5.0fs\n',etime(clock,stime));
  else
    Ymi = Ym;
  end
  if ~debug; clear Ysrc ; end


  if job.extopts.verb>2
    tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postLAS'));
    save(tmpmat,'Ysrc','Ycls','Ymi','Yb','T3th','vx_vol');
  end

  %  ---------------------------------------------------------------------
  %  Partitioning:
  %  ---------------------------------------------------------------------
  %  For most of the following adaptions further knowledge of special
  %  regions is helpfull. Also Ymi is maybe still a little bit inhomogen
  %  the alignment should work. Only strong inhomogenities can cause
  %  problems, especially for the blood vessel detection.
  %  But for bias correction the ROIs are important too, to avoid over
  %  corrections in special regions like the cerbellum and subcortex.
  %  ---------------------------------------------------------------------
  stime = cat_io_cmd('ROI segmentation (partitioning)');
  if job.extopts.SLC
    if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
      [Yl1,Ycls,YMF] = cat_vol_partvol(Ymi,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise,job,res.Ylesion); %,Ydt,Ydti);
      fprintf('%5.0fs\n',etime(clock,stime));
    else
      [Yl1,Ycls,YMF] = cat_vol_partvol(Ymi,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise,job,false(size(Ym)));
      fprintf('%5.0fs\n',etime(clock,stime));
      if isfield(res,'Ylesion') && sum(res.Ylesion(:)==0) && job.extopts.SLC==1
        cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main_SLC_noExpDef','SLC is set for manual lesions corection but no lesions were found!');
        fprintf('\n');
      end
    end
  else
    [Yl1,Ycls,YMF] = cat_vol_partvol(Ymi,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise,job,false(size(Ym)));
    fprintf('%5.0fs\n',etime(clock,stime));
    if job.extopts.expertgui && isfield(res,'Ylesion') && sum(res.Ylesion(:))>1000
      cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main_SLC_noExpDef',sprintf(['SLC is deactivated but there are %0.2f cm' ...
          char(179) ' of voxels with zero value inside the brain!'],prod(vx_vol) .* sum(res.Ylesion(:)) / 1000 ));
      fprintf('\n');
    end
  end
  if ~debug; clear YBG Ycr Ydt; end


  %  ---------------------------------------------------------------------
  %  Blood Vessel Correction
  %  ---------------------------------------------------------------------
  %  Blood vessel correction has to be done before the segmentation to
  %  remove high frequency strutures and avoid missclassifications.
  %  Problems can occure for strong biased images, because the partioning
  %  has to be done before bias correction.
  %  Of course we only want to do this for highres T1 data!
  %  ---------------------------------------------------------------------
  NS = @(Ys,s) Ys==s | Ys==s+1;
  if job.extopts.BVCstr && ~job.inv_weighting && all(vx_vol<2);
    stime = cat_io_cmd(sprintf('Blood vessel correction (BVCstr=%0.2f)',job.extopts.BVCstr));

    Ybv  = cat_vol_smooth3X(cat_vol_smooth3X( ...
      NS(Yl1,7) .* (Ymi*3 - (1.5-job.extopts.BVCstr)),0.3).^4,0.1)/3;

    % correct src images
    Ymi   = max(0,Ymi - Ybv*2/3);
    Ymi   = cat_vol_median3(Ymi,cat_vol_morph(Ybv>0.5,'dilate'));
    Ymis  = cat_vol_smooth3X(Ymi); Ymi(Ybv>0.5) = Ymis(Ybv>0.5); clear Ymis;

    % update classes
    Ycls{1} = min(Ycls{1},cat_vol_ctype(255 - Ybv*127));
    Ycls{2} = min(Ycls{2},cat_vol_ctype(255 - Ybv*127));
    Ycls{3} = max(Ycls{3},cat_vol_ctype(127*Ybv));

    fprintf('%5.0fs\n',etime(clock,stime));
    clear Ybv p0;
  end





  %% ---------------------------------------------------------------------
  %  Segmentation part
  %  ---------------------------------------------------------------------
  %  Now, it is time for skull-stripping (gcut,morph), AMAP tissue
  %  segmentation, and further tissue corrections (cleanup,LAS,finalmask).
  %  ---------------------------------------------------------------------


  %  -------------------------------------------------------------------
  %  skull-stipping
  %  -------------------------------------------------------------------
  %  For skull-stripping gcut is used in general, but a simple and very
  %  old function is still available as backup solution.
  %  Futhermore, both parts prepare the initial segmentation map for the
  %  AMAP function.
  %  -------------------------------------------------------------------
  if job.extopts.gcutstr>0 && job.extopts.gcutstr~=2
    %  -----------------------------------------------------------------
    %  gcut+: skull-stripping using graph-cut
    %  -----------------------------------------------------------------
    try
      stime = cat_io_cmd(sprintf('Skull-stripping using graph-cut (gcutstr=%0.2f)',job.extopts.gcutstr));
      [Yb,Yl1] = cat_main_gcut(Ymo,Yb,Ycls,Yl1,YMF,vx_vol,job.extopts);

      % extend gcut brainmask by brainmask derived from SPM12 segmentations if necessary
      if ~job.inv_weighting, Yb = Yb | Yb0; end

      fprintf('%5.0fs\n',etime(clock,stime));
    catch %#ok<CTCH>
      fprintf('\n'); cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_gcut:err99','Unknown error in cat_main_gcut. Use old brainmask.'); fprintf('\n');
      job.extopts.gcutstr = 99;
    end
  end

  % correct mask for skull-stripped images
  if skullstripped
    Yb = Yb.*(spm_read_vols(res.image0(1)) > 0);
  end

  %% -------------------------------------------------------------------
  %  AMAP segmentation
  %  -------------------------------------------------------------------
  %  Most corrections were done before and the AMAP routine is used with
  %  a low level of iterations and no further bias correction, because
  %  some images get tile artifacts.
  %  -------------------------------------------------------------------

  % correct for harder brain mask to avoid meninges in the segmentation
  Ymib = Ymi; Ymib(~Yb) = 0;
  rf = 10^4; Ymib = round(Ymib*rf)/rf;

  fprintf('prepare data for segmentation:  ');
  %  prepare data for segmentation
  if 1
    %% classic approach, consider the WMH!
    Kb2 = 4;
    cls2 = zeros([d(1:2) Kb2]);
    Yp0  = zeros(d,'uint8');
    for i=1:d(3)
        for k1 = 1:Kb2, cls2(:,:,k1) = Ycls{k1}(:,:,i); end
        % find maximum for reordered segmentations
        [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb2]),[],3);
        k1ind = [1 2 3 1 0 0 1 0];
        for k1 = 1:Kb2
          Yp0(:,:,i) = Yp0(:,:,i) + cat_vol_ctype((maxind == k1) .* (maxi~=0) * k1ind(k1) .* Yb(:,:,i));
        end
    end
    if ~debug, clear maxi maxind Kb k1 cls2; end
  else
    % more direct method ... a little bit more WM, less CSF
    % Yp0 = uint8(max(Yb,min(3,round(Ymi*3)))); Yp0(~Yb) = 0;
  end


  % use index to speed up and save memory
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

  % Yb source image because Amap needs a skull stripped image
  % set Yp0b and source inside outside Yb to 0
  Yp0b = Yp0(indx,indy,indz);
  Ymib = Ymib(indx,indy,indz);


  % remove non-brain tissue with a smooth mask and set values inside the
  % brain at least to CSF to avoid wholes for images with CSF==BG.
  if job.extopts.LASstr>0
    Ywmstd = cat_vol_localstat(single(Ymib),Yp0b==3,1,4);
    CSFnoise(1) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol);
    Ywmstd = cat_vol_localstat(cat_vol_resize(single(Ymib),'reduceV',vx_vol,vx_vol*2,16,'meanm'),...
      cat_vol_resize(Yp0==3,'reduceV',vx_vol,vx_vol*2,16,'meanm')>0.5,1,4);
    CSFnoise(2) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol);
    Ycsf = double(0.33 * Yb(indx,indy,indz)); spm_smooth(Ycsf,Ycsf,0.6*vx_vol);
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),0.5) * max(0.005,min(0.2,CSFnoise(1)/4)); % high-frequency noise
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),1.0) * max(0.005,min(0.2,CSFnoise(2)*1)); % high-frequency noise
    Ymib = max(Ycsf*0.8 .* cat_vol_smooth3X(Ycsf>0,2),Ymib);
    clear Ycsf;
    % Yb is needed for surface reconstruction
  %     if ~job.output.surface, clear Yb; end
  end

  % adaptive mrf noise
  if job.extopts.mrf>=1 || job.extopts.mrf<0;
    % estimate noise
    [Yw,Yg] = cat_vol_resize({Ymi.*(Ycls{1}>240),Ymi.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.15,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * 0.5;
    clear Yn Yg
  end

  % display something
  stime = cat_io_cmd(sprintf('Amap using initial SPM12 segmentations (MRF filter strength %0.2f)',job.extopts.mrf));

  % Amap parameters  - default sub=16 caused errors with highres data!
  % don't use bias_fwhm, because the Amap bias correction is not that efficient and also changes
  % intensity values
  Ymib = double(Ymib); n_iters = 50; sub = round(32/min(vx_vol)); n_classes = 3; pve = 5; bias_fwhm = 0; init_kmeans = 0;  %#ok<NASGU>
  if job.extopts.mrf~=0, iters_icm = 50; else iters_icm = 0; end %#ok<NASGU>

  % remove noisy background for kmeans
  if init_kmeans, Ymib(Ymib<0.1) = 0; end

  % do segmentation
  amapres = evalc(['prob = cat_amap(Ymib, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ' ...
    'job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);']);
  fprintf('%5.0fs\n',etime(clock,stime));

  % analyse segmentation ... the input Ym is normalized an the tissue peaks should be around [1/3 2/3 3/3]
  amapres = textscan(amapres,'%s'); amapres = amapres{1};
  th{1}   = cell2mat(textscan(amapres{11},'%f*%f'));
  th{2}   = cell2mat(textscan(amapres{12},'%f*%f'));
  th{3}   = cell2mat(textscan(amapres{13},'%f*%f'));


  if job.extopts.verb>1
    fprintf('    AMAP peaks: [CSF,GM,WM] = [%0.2f%s%0.2f,%0.2f%s%0.2f,%0.2f%s%0.2f]\n',...
      th{1}(1),char(177),th{1}(2),th{2}(1),char(177),th{2}(2),th{3}(1),char(177),th{3}(2));
  end
  if th{1}(1)<0 || th{1}(1)>0.6 || th{2}(1)<0.5 || th{2}(1)>0.9 || th{3}(1)<0.95-th{3}(2) || th{3}(1)>1.1
    warning('cat_main:amap',['AMAP estimated untypical tissue peaks that point to an \n' ...
                           'error in the preprocessing bevor the AMAP segmentation. ']);
  end
  % reorder probability maps according to spm order
  clear Yp0b Ymib;
  prob = prob(:,:,:,[2 3 1]);
  clear vol Ymib

  % finally use brainmask before cleanup that was derived from SPM12 segmentations and additionally include
  % areas where GM from Amap > GM from SPM12. This will result in a brainmask where GM areas
  % hopefully are all included and not cut
  if job.extopts.gcutstr >= 0.5 && ~job.inv_weighting
    Yb0(indx,indy,indz) = Yb0(indx,indy,indz) | ((prob(:,:,:,1) > 0) & Yb(indx,indy,indz)); % & ~Ycls{1}(indx,indy,indz));
    for i=1:3
      prob(:,:,:,i) = prob(:,:,:,i).*uint8(Yb0(indx,indy,indz));
    end
  end

  %% -------------------------------------------------------------------
  %  final cleanup
  %  There is one major parameter to control the strength of the cleanup.
  %  As far as the cleanup has a strong relation to the skull-stripping,
  %  cleanupstr is controlled by the gcutstr.
  %     Yp0ox = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255; Yp0o = zeros(d,'single'); Yp0o(indx,indy,indz) = Yp0ox;
  %     Yp0   = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b;
  %  -------------------------------------------------------------------
  if job.extopts.cleanupstr>0 && job.extopts.cleanupstr<=1
    %%
    for i=1:size(prob,4), [prob2(:,:,:,i),BB] = cat_vol_resize(prob(:,:,:,i),'reduceBrain',vx_vol,2,sum(prob,4)); end %#ok<AGROW>
    prob2 = clean_gwc(prob2,round(job.extopts.cleanupstr*4)); % old cleanup
    for i=1:size(prob,4), prob(:,:,:,i) = cat_vol_resize(prob2(:,:,:,i),'dereduceBrain',BB); end; clear prob2


    [Ycls,Yp0b] = cat_main_cleanup(Ycls,prob,Yl1(indx,indy,indz),Ymo(indx,indy,indz),job.extopts,job.inv_weighting,vx_volr,indx,indy,indz);
  else
    if 0% job.extopts.cleanupstr == 2 % old cleanup for tests
      stime = cat_io_cmd('Old cleanup');

      for i=1:size(prob,4), [prob2(:,:,:,i),BB] = cat_vol_resize(prob(:,:,:,i),'reduceBrain',vx_vol,2,sum(prob,4)); end %#ok<AGROW>
      prob2 = clean_gwc(prob2,round(job.extopts.cleanupstr*4)); % old cleanup
      for i=1:size(prob,4), prob(:,:,:,i) = cat_vol_resize(prob2(:,:,:,i),'dereduceBrain',BB); end; clear prob2

      fprintf('%5.0fs\n',etime(clock,stime));
    end
    for i=1:3
       Ycls{i}(:) = 0; Ycls{i}(indx,indy,indz) = prob(:,:,:,i);
    end
    Yp0b = Yb(indx,indy,indz);
  end;
  if ~debug; clear Ymo; end
  clear prob




  %% -------------------------------------------------------------------
  %  Correction of WM hyperintensities
  %  -------------------------------------------------------------------
  %  The correction of WMH should be important for a correct normalization.
  %  It is only important to close the mayor WMH structures, and further
  %  closing can lead to problems with small gyri. So keep it simple here
  %  and maybe add further refinements in the partitioning function.
  %  -------------------------------------------------------------------
  LAB  = job.extopts.LAB;
  Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b;
  Ywmhrel = single(Ycls{1})/255 .* NS(Yl1,23);
  qa.subjectmeasures.WMH_abs    = sum(Ywmhrel(:));                                            % absolute WMH volume without PVE
  qa.subjectmeasures.WMH_rel    = 100*qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(0.5/3*255));   % relative WMH volume to TIV without PVE
  qa.subjectmeasures.WMH_WM_rel = 100*qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(2.5/3*255));   % relative WMH volume to WM without PVE
  qa.subjectmeasures.WMH_abs    = prod(vx_vol)/1000 * qa.subjectmeasures.WMH_abs;             % absolute WMH volume without PVE in cm^3
  clear Ywmhrel Yp0




  %% correction for normalization [and final segmentation]
  if ( (job.extopts.WMHC && job.extopts.WMHCstr>0) || job.extopts.SLC) && ~job.inv_weighting
    % display something
    %{
    if job.extopts.WMHC==1
      cat_io_cmd(sprintf('Internal WMH correction for spatial normalization')); % (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    elseif job.extopts.WMHC>1
      cat_io_cmd(sprintf('Permanent WMH correction')); % (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    end
    fprintf('\n');
    if job.extopts.SLC==1
      cat_io_cmd('Internal stroke lesion correction for spatial normalization');
    elseif job.extopts.SLC>1
      cat_io_cmd('Permanent stroke lesion correction');
    end
    fprintf('\n');
    %}

    % prepare correction map
    Ynwmh = NS(Yl1,LAB.TH) | NS(Yl1,LAB.BG) | NS(Yl1,LAB.HC) | NS(Yl1,LAB.CB) | NS(Yl1,LAB.BS);
    Ynwmh = cat_vol_morph(cat_vol_morph( Ynwmh, 'dd', 8 , vx_vol),'dc',12 , vx_vol) & ...
            ~cat_vol_morph( NS(Yl1,LAB.VT), 'dd', 4 , vx_vol);
    Ywmh  = cat_vol_morph( Ycls{7}>0, 'dd', 1.5);
    Ywmh  = Ycls{7}>0 | (~Ynwmh & (Ycls{2}==255 | ...
            cat_vol_morph( cat_vol_morph(Ycls{2}>128 | Ywmh,'ldc',1) ,'de' , 1.5)));
    Ywmh  = Ywmh .* cat_vol_smooth3X(Ywmh,0.5); % smooth inside


    %% transfer tissue from GM and CSF to WMH
    if job.extopts.SLC>0
      % WMHs and lesions
      if job.extopts.SLC==1
        Yls     = res.Ylesion;
        Ycls{8} = cat_vol_ctype( Yls*255 );
      elseif job.extopts.SLC==2
        Yls     = NS(Yl1,LAB.LE)>0.5 | res.Ylesion;
        Ycls{8} = cat_vol_ctype( Yls  .* single(Ycls{1})  +  Yls  .* single(Ycls{3})  + 255*single(res.Ylesion) );
      end
      Ycls{7} = cat_vol_ctype( Ywmh .* single(Ycls{1})  +  Ywmh .* single(Ycls{3}));
      Ycls{1} = cat_vol_ctype( single(Ycls{1}) .* (1 - Ywmh - single(Yls)) );
      Ycls{3} = cat_vol_ctype( single(Ycls{3}) .* (1 - Ywmh - single(Yls)) );
    else
      % only WMHS
      Ycls{7} = cat_vol_ctype( Ywmh .* single(Ycls{1})  +  Ywmh .* single(Ycls{3}));
      Ycls{1} = cat_vol_ctype( single(Ycls{1}) .* (1 - Ywmh) );
      Ycls{3} = cat_vol_ctype( single(Ycls{3}) .* (1 - Ywmh) );
    end
    if ~debug, clear Ynwmh Ywmh Yls; end


    % different types of WMH correction as GM, WM or extra class
    % different types of lesion correction as CSF or extra class
    if numel(Ycls)>7 && job.extopts.SLC==1
      if job.extopts.WMHC<2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*2/5 + single(Ycls{8})*1/5,'uint8');
      elseif job.extopts.WMHC==2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*3/5 + single(Ycls{8})*1/5,'uint8');
      elseif job.extopts.WMHC==3
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*4/5 + single(Ycls{8})*1/5,'uint8');
      end
    elseif numel(Ycls)>7 && job.extopts.SLC==2
      if job.extopts.WMHC<2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*2/5 + single(Ycls{8})*1.5/5,'uint8');
      elseif job.extopts.WMHC==2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*3/5 + single(Ycls{8})*1.5/5,'uint8');
      elseif job.extopts.WMHC==3
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*4/5 + single(Ycls{8})*1.5/5,'uint8');
      end
    else % no stroke lesion handling
      if job.extopts.WMHC<2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*2/5,'uint8');
      elseif job.extopts.WMHC==2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*3/5,'uint8');
      elseif job.extopts.WMHC==3
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*4/5,'uint8');
      end
    end

  else
    Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5,'uint8');

    if qa.subjectmeasures.WMH_rel>3 || qa.subjectmeasures.WMH_WM_rel>5 % #% of the TIV or the WM are affected
      cat_warnings = cat_io_addwarning(cat_warnings,...
        'MATLAB:SPM:CAT:cat_main:uncorrectedWMH',...
        sprintf('Uncorrected WM lesions greater (%2.2f%%%%%%%% of the WM)!\\n',qa.subjectmeasures.WMH_rel));
    end
  end
  Yp0b = Yp0b(indx,indy,indz);
  clear Yclsb;

  if job.extopts.verb>2
    tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'preDartel'));
    save(tmpmat,'Yp0','Ycls','Ymi','T3th','vx_vol','Yl1');
    clear Yp0;
  end


  % clear last 3 tissue classes to save memory
  % please do not try to write out these segmentations because class 4-6 are form SPM12
  % and class 1-3 from CAT12 and these are completely different segmentation approaches
  %if all(cell2mat(struct2cell(job.output.TPMC)')==0)
  %  for i=6:-1:4, Ycls(i)=[]; end
  %end

  %% ---------------------------------------------------------------------
  %  Affine registration:
  %  Use the segmentation result (p0 label image) for the final affine
  %  registration to a pseudo T1 image of the TPM. The dartel tempalte
  %  can not be used here, due to special TPM image for chrildren or other
  %  species!
  %  ---------------------------------------------------------------------
  if 0
    % parameter
    aflags = struct('sep',job.opts.samp,'regtype',job.opts.affreg,'WG',[],'WF',[],'globnorm',0);

    % VG template
    cid = [2 3 1];
    VG = tpm.V(1); VG.dat = zeros(VG.dim,'uint8'); VG.dt = [2 0]; VG.pinfo(3) = 0;
    for ci=1:3
      Yt = spm_read_vols(tpm.V(ci));
      VG.dat = VG.dat + cat_vol_ctype(Yt * 255 * cid(ci)/3,'uint8'); clear Yt
    end

    % VF image
    VF = VT; VF.fname = [tempname '.nii']; VF = rmfield(VF,'private'); VF.pinfo(3)=0;
    VF.dat = zeros(d,'uint8'); VF.dt = [2 0];
    VF.dat(indx,indy,indz) = Yp0b;

    % just brain mask - this improves affine registration but lead to worse
    % resultes in dartel normalization .. retry later without CSF
    % VG.dat = uint8(VG.dat>85)*255; VF.dat = uint8(VF.dat>85)*255;

    % smooth source with job.opts.samp mm
    VF = spm_smoothto8bit(VF,job.opts.samp/2); % smoothing, because of the TPM smoothness!

    %% affreg registration for one tissue segment

    % deactivated on 20160405 because something failed in the logitudinal process
    warning off;
    Affine = spm_affreg(VG, VF, aflags, res.Affine);
    warning on;
    rf=10^6; Affine    = fix(Affine * rf)/rf;
  else
    Affine = res.Affine;
  end

  res.Affine0        = res.Affine;
  res.Affine         = Affine;

  clear VG VF cid %tpm

else
%% SPM segmentation input
%  ------------------------------------------------------------------------
%  Here, DARTEL and PBT processing is prepared.
%  We simply use the SPM segmentation as it is without further modelling of
%  a PVE or other refinements.
%  ------------------------------------------------------------------------
  job.extopts.WMHC = 0;
  job.extopts.SLC  = 0;

  NCstr.labels = {'none','full','light','medium','strong','heavy'};
  NCstr.values = {0 1 2 -inf 4 5};

  % here we need the c1 filename
  VT0 = res.imagec(1);
  [pth,nam] = spm_fileparts(VT0.fname);

  vx_vol              = sqrt(sum(VT.mat(1:3,1:3).^2));          % voxel size of the processed image
  cat_warnings        = struct('identifier',{},'message',{});   % warning structure from cat_main_gintnorm
  NS                  = @(Ys,s) Ys==s | Ys==s+1;                % for side independent atlas labels

  % QA WMH values required by cat_vol_qa later
  qa.subjectmeasures.WMH_abs    = nan;  % absolute WMH volume without PVE
  qa.subjectmeasures.WMH_rel    = nan;  % relative WMH volume to TIV without PVE
  qa.subjectmeasures.WMH_WM_rel = nan;  % relative WMH volume to WM without PVE
  qa.subjectmeasures.WMH_abs    = nan;  % absolute WMH volume without PVE in cm^3

  % load SPM segments
  [pp,ff,ee] = spm_fileparts(res.image0(1).fname);
  Ycls{1} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c1' ff ee])))*255);
  Ycls{2} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c2' ff ee])))*255);
  Ycls{3} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c3' ff ee])))*255);

  % create (resized) label map and brainmask
  Yp0  = single(Ycls{3})/5 + single(Ycls{1})/5*2 + single(Ycls{2})/5*3;
  Yb   = Yp0>0.5;

  % load original images and get tissue thresholds
  Ysrc = spm_read_vols(spm_vol(fullfile(pp,[ff ee])));
  WMth = double(max(clsint(2),...
           cat_stat_nanmedian(cat_stat_nanmedian(cat_stat_nanmedian(Ysrc(Ycls{2}>192))))));
  T3th = [ min([  clsint(1) - diff([clsint(1),WMth]) ,clsint(3)]) , clsint(2) , WMth];
  if T3th(3)<T3th(2) % inverse weighting allowed
    job.inv_weighting   = 1;
  else
    job.inv_weighting   = 0;
  end
  if ~debug, clear Ysrc; end


  % the intensity normalized images are here represented by the segmentation
  Ym   = Yp0/255;
  Ymi  = Yp0/255;


  % low resolution Yp0b
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));
  Yp0b = Yp0(indx,indy,indz);
  if ~debug, clear Yp0 Yb; end

  % load atlas map and prepare filling mask YMF
  % compared to CAT default processing, we have here the DARTEL mapping, but no individual refinement
  Vl1 = spm_vol(job.extopts.cat12atlas{1});
  Yl1 = cat_vol_ctype(spm_sample_vol(Vl1,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
  Yl1 = reshape(Yl1,size(Ym)); [D,I] = cat_vbdist(single(Yl1>0)); Yl1 = Yl1(I);
  YMF = NS(Yl1,job.extopts.LAB.VT) | NS(Yl1,job.extopts.LAB.BG) | NS(Yl1,job.extopts.LAB.BG);

  fprintf('%5.0fs\n',etime(clock,stime));
end



%% ---------------------------------------------------------------------
%  Spatial Registration with Dartel or Shooting
%  ---------------------------------------------------------------------
  Yclsd = Ycls(1:2); % use only GM and WM for deformation
  if job.extopts.WMHC>0 && numel(Ycls)>6
    Yclsd{2} = cat_vol_ctype(min(255,single(Ycls{2}) + single(Ycls{7}))); % set WMHs as WM in some cases
  end

  if (~isempty(job.extopts.xasl_lesion{1})) || (job.extopts.SLC && isfield(res,'Ylesion') && sum(res.Ylesion(:)>0))
    % lesion detection in the original space with the original data
    LSstr   = 0.5;
    Yvt     = cat_vol_morph( NS(Yl1,job.extopts.LAB.VT),'do',4,vx_vol);      % open to get lesions close to the ventricle
    Yvt     = cat_vol_morph( Yvt ,'dd',4,vx_vol);                            % add some voxels for smoothness
    res.Ylesion = cat_vol_ctype( single(res.Ylesion) .* (1 - (Yvt & Ym>0.9 & Ym<1.1) ));
	res.YlesionFull = cat_vol_ctype( single(res.YlesionFull) .* (1 - (Yvt & Ym>0.9 & Ym<1.1) ));
    if ~debug, clear Yvt Ybgvt Ybgn;  end
    % add lesion of automatic lesion estimation? - in development
    if job.extopts.WMHC>3
      res.Ylesion = cat_vol_ctype( single(res.Ylesion) + ...
        255* smooth3( Ym<1.5/3 & cat_vol_morph(NS(Yl1,job.extopts.LAB.LE),'dd',4*(1-LSstr))) );
	  res.YlesionFull = cat_vol_ctype( single(res.YlesionFull) + ...
        255* smooth3( Ym<1.5/3 & cat_vol_morph(NS(Yl1,job.extopts.LAB.LE),'dd',4*(1-LSstr))) );
    end
    %Ylesions = cat_vol_smooth3X(single(res.Ylesion)/255,4); % final smoothing to have soft boundaries
	Ylesions = cat_vol_smooth3X(single(res.YlesionFull),1.5); % final smoothing to have soft boundaries
  else
	Ylesions = zeros(size(Ym),'single');
  end

























  %%% Cost Function Masking -> if Lesion_*.nii are found, these are masked out.
  %%% Otherwise, this part leaves the segmentation untouched
  %%% It only stores an original segmentation image before DARTEL
  xASL_im_SaveOriginal4CAT(Yclsd, job.channel.vols0{1});
  %Yclsd         = xASL_im_LesionRemoval4CAT(Yclsd, job.channel.vols0{1});
  %Ycls          = xASL_im_LesionRemoval4CAT(Ycls, job.channel.vols0{1}); %%% ExploreASL fix
  % LesionImOut         = repmat(LesionImOut,[1 1 1 3]);
  % trans.warped.y(LesionImOut)     = Yy(LesionImOut); % within the lesion mask, replace DARTEL flow field by SPM affine+DCT flow field
  %%%
  %%% ANOTHER EXPLOREASL HACK: remove NaNs from transformations before DARTEL
  if exist('Yy', 'var')
     Yy = single(xASL_im_FillNaNs(Yy, 1, job.extopts.xasl_quality, [1.5 1.5 1.5])); % SPM doesn't like double format Yy
  end

  %% call Dartel/Shooting registration
  if ~job.extopts.xasl_quality % ExploreASL hack fewer iterations/lower resolution to speed up
      job.extopts.reg.nits = 4; % less iterations
      job.extopts.regstr = 15; % low resolution
  end

  [trans,res.ppe.reg] = cat_main_registration(job,res,Yclsd,Yy,tpm.M,Ylesions);
  clear Yclsd Ylesions;
  if ~do_dartel
    if job.extopts.regstr == 0
      fprintf('Dartel registration is not required.\n');
    else
      fprintf('Shooting registration is not required.\n');
    end
  end

  % EXPLOREASL HACK
  if isfield(trans,'warped')
      trans.warped.y = single(xASL_im_FillNaNs(trans.warped.y, 3, job.extopts.xasl_quality, [1.5 1.5 1.5])); % SPM doesn't like double format here
  end


  %% update WMHs?
  if numel(Ycls)>6
    if isfield(trans,'warped')
      %% load template
      VwmA = spm_vol([job.extopts.templates{end},',2']);
      %VgmA = spm_vol([job.extopts.templates{end},',1']);
      if any( VwmA.dim ~= trans.warped.odim )
        % interpolation
        yn = numel(trans.warped.y);
        p  = ones([4,yn/3],'single');
        p(1,:) = trans.warped.y(1:yn/3);
        p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
        p(3,:) = trans.warped.y(yn/3*2+1:yn);
        amat   = VwmA.mat \ trans.warped.M1;
        p      = amat(1:3,:) * p;

        Yy = zeros([res.image(1).dim(1:3),3],'single');
        Yy(1:yn/3)        = p(1,:);
        Yy(yn/3+1:yn/3*2) = p(2,:);
        Yy(yn/3*2+1:yn)   = p(3,:);

        Yy = double(Yy);
      else
        Yy = double(trans.warped.y);
      end
      YwmA = single( spm_sample_vol( VwmA ,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),1)); YwmA = reshape(YwmA,size(Ym));
      %YgmA = single( spm_sample_vol( VgmA ,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),1)); YgmA = reshape(YgmA,size(Ym));
    else
      %% load template
      VwmA = spm_vol([job.extopts.templates{end},',2']);
      %VgmA = spm_vol([job.extopts.templates{end},',1']);
      if any( VwmA.dim ~= size(Yy(:,:,:,1)) )
        % interpolation
        yn = numel(trans.atlas.Yy);
        p  = ones([4,yn/3],'single');
        p(1,:) = trans.atlas.Yy(1:yn/3);
        p(2,:) = trans.atlas.Yy(yn/3+1:yn/3*2);
        p(3,:) = trans.atlas.Yy(yn/3*2+1:yn);
        amat   = VwmA.mat \ M1;
        p      = amat(1:3,:) * p;

        Yy = zeros([res.image(1).dim(1:3),3],'single');
        Yy(1:yn/3)        = p(1,:);
        Yy(yn/3+1:yn/3*2) = p(2,:);
        Yy(yn/3*2+1:yn)   = p(3,:);

        Yy = double(Yy);
      else
        Yy = double(trans.warped.y);
      end
      YwmA = single( spm_sample_vol( VwmA ,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),1)); YwmA = reshape(YwmA,size(Ym));
      %YgmA = single( spm_sample_vol( VgmA ,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),1)); YgmA = reshape(YgmA,size(Ym));
    end
    %%
    Yclst   = cat_vol_ctype( single(Ycls{7}) .* YwmA );
    Ycls{2} = Ycls{2} + (Ycls{7} - Yclst);
    Ycls{7} = Yclst;

    %Yp0   = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*5/255;
    %Yclst = cat_vol_morph( cat_vol_morph( (YwmA + YgmA)>0.9 & Ym<1.25/3,'do',1) , 'dd', 4) & Ym<1.25/3;
  end

%%  --------------------------------------------------------------------
%  write results
%  ---------------------------------------------------------------------

if job.extopts.WMHC<2 && numel(Ycls)>6
  Ycls{1} = Ycls{1} + Ycls{7}; % WMH as GM
elseif job.extopts.WMHC==2
  Ycls{2} = Ycls{2} + Ycls{7}; % WMH as WM
elseif job.extopts.WMHC>=3
  Ycls{2} = Ycls{2} + Ycls{7}; % WMH as own class
end



stime = cat_io_cmd('Write result maps');
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*5/255;

% bias, noise and global corrected without masking for subject space and with masking for other spaces
cat_io_writenii(VT0,Ym,mrifolder,'m', ...Dartel
  'bias and noise corrected, global intensity normalized','uint16',[0,0.0001], ...
  min([1 0 2],[job.output.bias.native job.output.bias.warped job.output.bias.dartel]),trans);
cat_io_writenii(VT0,Ym.*(Yp0>0.1),mrifolder,'m', ...
  'bias and noise corrected, global intensity normalized (masked due to normalization)','uint16',[0,0.0001], ...
  min([0 1 0],[job.output.bias.native job.output.bias.warped job.output.bias.dartel]),trans);

% bias, noise and local intensity corrected without masking for subject space and with masking for other spaces
cat_io_writenii(VT0,Ymi,mrifolder,'mi', ...
  'bias and noise corrected, local intensity normalized','uint16',[0,0.0001], ...
  min([1 0 2],[job.output.las.native job.output.las.warped job.output.las.dartel]),trans);
cat_io_writenii(VT0,Ymi.*(Yp0>0.1),mrifolder,'mi', ...
  'bias and noise corrected, local intensity normalized (masked due to normalization)','uint16',[0,0.0001], ...
  min([0 1 0],[job.output.las.native job.output.las.warped job.output.las.dartel]),trans);

% Yp0
cat_io_writenii(VT0,Yp0,mrifolder,'p0','label map','uint8',[0,5/255],job.output.label,trans);
clear Yp0;

% partitioning
cat_io_writenii(VT0,Yl1,mrifolder,'a0','brain atlas map for major structures and sides',...
  'uint8',[0,1],job.output.atlas,trans);

% class maps 1-3
fn = {'GM','WM','CSF','head','head','background'};
for clsi=1:3
  cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
    min([1 1 0 3],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
    job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
  cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
    min([0 0 3 0],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
    job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
end

% write WMH class maps
if job.extopts.WMHC>=3 && job.extopts.WMHCstr>0 && ~job.inv_weighting;
  cat_io_writenii(VT0,single(Ycls{7})/255,mrifolder,'p7','WMH tissue map','uint8',[0,1/255],...
    min([1 1 0 3],[job.output.WMH.native job.output.WMH.warped ...
    job.output.WMH.mod job.output.WMH.dartel]),trans); % 1 0 0 0
  cat_io_writenii(VT0,single(Ycls{7})/255,mrifolder,'p7','WMH tissue map','uint16',[0,1/255],...
    min([0 0 3 0],[job.output.WMH.native job.output.WMH.warped ...
    job.output.WMH.mod job.output.WMH.dartel]),trans); % 0 1 2 2


elseif any([job.output.WMH.native job.output.WMH.warped ...
    job.output.WMH.mod job.output.WMH.dartel])
  % should write, but can not ...
  if job.extopts.WMHC<3
    cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_WMHC:output','Cannot write WMH images because no seperate class is used (set WMHC>=3).');
  elseif job.extopts.WMHCstr==0
    cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_WMHC:output','Cannot write WMH images because WMHCstr is 0.');
  elseif job.inv_weighting
    cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_WMHC:output','Cannot write WMH images because inverse contrast was detected.');
  end
  %
end
% export lesions
if job.extopts.SLC>0 && numel(Ycls)>7
  cat_io_writenii(VT0,single(Ycls{8}),mrifolder,'p8','stroke lesion map','uint8',[0,1/255],...
    min([1 1 0 3],[job.output.SL.native job.output.SL.warped ...
    job.output.SL.mod job.output.SL.dartel]),trans); % 1 0 0 0
  cat_io_writenii(VT0,single(Ycls{8}),mrifolder,'p8','stroke lesion map','uint16',[0,1/255],...
    min([0 0 3 0],[job.output.WMH.native job.output.SL.warped ...
    job.output.SL.mod job.output.SL.dartel]),trans); % 1 0 0 0
end

% developer - intensity scaled tissue classe maps
% ----------------------------------------------------------------------
% The strong normalization of the T1 data can directly be used as tissue
% segmentation. The Ymi images is scaled to get similar maps for each
% tissue class, with good visible differences in the sulci.
job.output.intsegments = job.extopts.experimental;
if job.output.intsegments
  if (any(tc(:)) || job.extopts.WMHC==3 && job.extopts.WMHCstr>0 && ~job.inv_weighting);

    % intensity scaled tissue maps
    Yclsi = cell(1,3);
    for clsi=1:3
      clsid = [2 3 1];
      Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*5/255;
      Yclsi{clsi} = Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),clsid(clsi));
      switch clsi
        case 1
          Yclsi{clsi} = Yclsi{clsi} .* (Ycls{clsi}>0) + ...
            Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),3) .* (Ycls{2}==0 & (Ycls{1}>0));
        case 2
          Yclsi{clsi} = Yclsi{clsi} .* (Ycls{clsi}>0) +  ...
            Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2) .* (Ycls{2}==255 & ~cat_vol_morph(Ycls{2}>192,'e'));
          Ywmhp = cat_vol_ctype( Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2) .* cat_vol_morph(Ycls{2}>192,'e') * 255);
        case 3
          Yclsi{clsi} = Yclsi{clsi} + ...
            (Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2)) .* (Ycls{1}==0 & (Ycls{3}>0)) + ...
            (Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),3)) .* (Ycls{2}==0 & (Ycls{3}>0));
      end
      Yclsi{clsi} = cat_vol_ctype(Yclsi{clsi} * 255);
    end
    clear Yp0;

    % class maps 1-3
    % Yclss = single(Yclsi{1} + Yclsi{2} + Yclsi{3} + Ywmhp + Ywmh)/255; % + single(Ywmh)/255;
    for clsi=1:3
      cat_io_writenii(VT0,single(Yclsi{clsi})/255,mrifolder,sprintf('pi%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
        min([1 1 0 3],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
        job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
      cat_io_writenii(VT0,single(Yclsi{clsi})/255,mrifolder,sprintf('pi%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
        min([0 0 2 0],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
        job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
    end
    clear Yclsi;

    % write WMH class maps
    if job.extopts.WMHC>=3 && job.extopts.WMHCstr>0 && ~job.inv_weighting;
      cat_io_writenii(VT0,(single(Ywmhp) + single(Ywmh))/255,mrifolder,...
        'pi7','WMH tissue map','uint8',[0,1/255],...
        min([1 1 0 2],[job.output.WMH.native job.output.WMH.warped ...
        job.output.WMH.mod job.output.WMH.dartel]),trans); % 1 0 0 0
      cat_io_writenii(VT0,(single(Ywmhp) + single(Ywmh))/255,mrifolder,...
        'pi7','WMH tissue map','uint16',[0,1/255],...
        min([0 0 2 0],[job.output.WMH.native job.output.WMH.warped ...
        job.output.WMH.mod job.output.WMH.dartel]),trans); % 0 1 2 2
    end
    clear Ywmhp;
  end
end
% ----------------------------------------------------------------------


% classe maps 4-6 (for full TPM/template creation, e.g. for apes)
if any(cell2mat(struct2cell(job.output.TPMC)'))
  for clsi=4:6
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
      min([1 1 0 3],[job.output.TPMC.native job.output.TPMC.warped ...
      job.output.TPMC.mod job.output.TPMC.dartel]),trans);
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
      min([0 0 3 0],[job.output.TPMC.native job.output.TPMC.warped ...
      job.output.TPMC.mod job.output.TPMC.dartel]),trans);
  end
end
%clear cls clsi fn Ycls; % we need this maps later for the ROIs

% write jacobian determinant
if job.output.jacobian.warped
  %%
  if do_dartel==2 % shooting
    dt = trans.jc.dt2;
    dx = 10; % smaller values are more accurate, but large look better;
    [D,I] = cat_vbdist(single(~(isnan(dt) | dt<0 | dt>100) )); D=min(1,D/min(dx,max(D(:))));
    dt = dt(I); dt = dt .* ((1-D) + D); dt(isnan(dt))=1;
    %dt = 1/max(eps,dt);
    clear y0 D I
  else %dartel
    [y0, dt] = spm_dartel_integrate(reshape(trans.jc.u,[trans.warped.odim(1:3) 1 3]),[1 0], 6);
    clear y0
  end
  N      = nifti;
  N.dat  = file_array(fullfile(pth,mrifolder,['wj_', nam, '.nii']),trans.warped.odim(1:3),...
             [spm_type('float32') spm_platform('bigend')],0,1,0);
  N.mat  = trans.warped.M1;
  N.mat0 = trans.warped.M1;
  N.descrip = ['Jacobian' VT0.descrip];
  create(N);
  N.dat(:,:,:) = dt;
end


%%
if job.output.ROI || any(cell2mat(struct2cell(job.output.atlas)'))
  % get atlases
  FAF = job.extopts.atlas;
  FA  = {}; fai = 1;
  AN  = fieldnames(job.output.atlases);
  for ai = 1:numel(AN)
    fafi = find(cellfun('isempty',strfind(FAF(:,1),[AN{ai} '.']))==0);
    if ~isempty(fafi) && job.output.atlases.(AN{ai}), FA(fai,:) = FAF(fafi,:); fai = fai+1; end %#ok<AGROW>
  end
end

if exist('FA','var')
  if isempty(FA)
    % deactivate output
    FN = job.output.atlas;
    for ai = 1:numel(AN)
      job.output.atlas.(FN{ai}) = 0;
    end
  else
    % get atlas resolution
    % we sort the atlases to reduce data resampling
    VA = spm_vol(char(FA(:,1)));
    for ai=1:numel(VA), VAvx_vol(ai,:) = sqrt(sum(VA(ai).mat(1:3,1:3).^2)); end   %#ok<AGROW>
    [VAs,VAi] = sortrows(VAvx_vol);
    FA = FA(VAi,:); VA = VA(VAi,:); VAvx_vol = VAvx_vol(VAi,:); %clear VA;
  end
end

%% write atlas output
if any(cell2mat(struct2cell(job.output.atlas)'))
  for ai=1:size(FA,1)
    [px,atlas] = fileparts(FA{ai,1});

    % map atlas in native space
    Vlai = spm_vol(FA{ai,1});
    if any( Vlai.dim ~= trans.warped.odim )
      % interpolation
      Vlai = spm_vol(FA{ai,1});
      yn = numel(trans.warped.y);
      p  = ones([4,yn/3],'single');
      p(1,:) = trans.warped.y(1:yn/3);
      p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
      p(3,:) = trans.warped.y(yn/3*2+1:yn);
      amat   = Vlai.mat \ trans.warped.M1;
      p      = amat(1:3,:) * p;

      Yy = zeros([res.image(1).dim(1:3),3],'single');
      Yy(1:yn/3)        = p(1,:);
      Yy(yn/3+1:yn/3*2) = p(2,:);
      Yy(yn/3*2+1:yn)   = p(3,:);

      Yy = double(Yy);
    else
      Yy = double(trans.warped.y);
    end
    Ylai = cat_vol_ctype(spm_sample_vol(Vlai,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),0));
    Ylai = reshape(Ylai(:),trans.native.Vi.dim);
    if ~debug, clear Yy; end

    % write map (mri as tissue subforder and mri_atals as ROI subfolder)
    if isempty(mrifolder), amrifolder = ''; else amrifolder = 'mri_atlas'; end
    cat_io_writenii(VT0,Ylai,amrifolder,[atlas '_'],[atlas ' original'],...
      'uint8',[0,1],job.output.atlas,trans);
    if ~debug, clear Vlai Ylai; end
  end
end


%% deformations y - dartel > subject
if job.output.warps(1)
	disp('Saving the subject-DARTEL transformation as NIfTI');
    % ExploreASL hack: Here, create flow field without DARTEL
    Yy        = spm_diffeo('invdef',trans.atlas.Yy,[121 145 121],eye(4),res.image.mat); % assuming 1.5 mm MNI
    %Yy        = spm_diffeo('invdef',trans.warped.y,trans.warped.odim,eye(4),trans.warped.M0);
    N         = nifti;

	N.dat     = file_array(fullfile(pth,mrifolder,['y_', nam, '_withoutDARTEL.nii']),[121,145,121,1,3],'float32',0,1,0);
	N.mat     = M1;
	N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
	N.dat(:,:,:,:,:) = reshape(Yy,[121,145,121,1,3]);
    if ~debug, clear Yy; end
	clear N;

	% ExploreASL hack, here create flow field + DARTEL
    if  do_dartel
        Yy        = spm_diffeo('invdef',trans.warped.y,trans.warped.odim,eye(4),trans.warped.M0);
        N         = nifti;
        N.dat     = file_array(fullfile(pth,mrifolder,['y_', nam, '.nii']),[trans.warped.odim(1:3),1,3],'float32',0,1,0);
        N.mat     = trans.warped.M1;
        N.mat0    = trans.warped.M1;
        N.descrip = 'Deformation';
        create(N);
        N.dat(:,:,:,:,:) = reshape(Yy,[trans.warped.odim,1,3]);
        clear N
    end

end

%% deformation iy - normalized > subject
if job.output.warps(2)
	disp('Saving the DARTEL-subject transformation as NIfTI (i.e. inverse)');
  if any(trans.native.Vo.dim~=trans.native.Vi.dim)
    %%
    vx_voli  = sqrt(sum(trans.native.Vi.mat(1:3,1:3).^2));
    vx_volo  = sqrt(sum(trans.native.Vo.mat(1:3,1:3).^2));
    eyev = eye(4); eyev([1 6 11]) = eyev([1 6 11]) .* vx_volo./vx_voli;
    Yy2  = zeros([trans.native.Vo.dim 1 3],'single');
    for k1=1:3
      for i=1:trans.native.Vo.dim(3),
        Yy2(:,:,i,:,k1) = trans.warped.M1(k1,4) + trans.warped.M1(k1,k1) * ...
          single(spm_slice_vol(trans.warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]), ...
          trans.native.Vo.dim(1:2),[1,NaN])); % adapt for res
      end
    end
  else
    yn = numel(trans.warped.y);
    p  = ones([4,yn/3],'single');
    p(1,:) = trans.warped.y(1:yn/3);
    p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
    p(3,:) = trans.warped.y(yn/3*2+1:yn);
    p      = trans.warped.M1(1:3,:) * p;

    Yy2 = zeros([res.image(1).dim(1:3),1,3],'single');
    Yy2(1:yn/3)        = p(1,:);
    Yy2(yn/3+1:yn/3*2) = p(2,:);
    Yy2(yn/3*2+1:yn)   = p(3,:);
  end
  clear p;

  % f2 = spm_diffeo('resize', f1, dim)
  % write new output
  Ndef      = nifti;
  Ndef.dat  = file_array(fullfile(pth,mrifolder,['iy_', nam, '.nii']),[res.image0(1).dim(1:3),1,3],...
              [spm_type('float32') spm_platform('bigend')],0,1,0);
  Ndef.mat  = res.image0(1).mat;
  Ndef.mat0 = res.image0(1).mat;
  Ndef.descrip = 'Inverse Deformation';
  create(Ndef);
  Ndef.dat(:,:,:,:,:) = Yy2;

  if ~debug, clear Yy2; end
end
fprintf('%5.0fs\n',etime(clock,stime));



%% ---------------------------------------------------------------------
%  surface creation and thickness estimation
%  ---------------------------------------------------------------------
if (job.output.surface || any( [job.output.ct.native job.output.ct.warped job.output.ct.dartel] )) && ...
   ~(job.output.surface==9 && job.output.ROI==0 && ~any( [job.output.ct.native job.output.ct.warped job.output.ct.dartel] ))
    % ... not required, if only thickness but no output
  stime = cat_io_cmd('Surface and thickness estimation');;

  % specify WM/CSF width/depth/thickness estimation
  if job.output.surface>10
    job.output.surface=job.output.surface-10;
    WMT = 1;
  else
    WMT = 0;
  end

  if job.extopts.experimental || job.extopts.expertgui==2
    WMT = 1;
  end

  % specify surface
  switch job.output.surface
    case 1, surf = {'lh','rh'};
    case 2, surf = {'lh','rh','lc','rc'};
    case 3, surf = {'lh'};
    case 4, surf = {'rh'};
    % fast surface reconstruction without simple spherical mapping
    case 5, surf = {'lhfst','rhfst'};
    case 6, surf = {'lhfst','rhfst','lcfst','rcfst'};
    % fast surface reconstruction with simple spherical mapping
    case 7, surf = {'lhsfst','rhsfst'};
    case 8, surf = {'lhsfst','rhsfst','lcsfst','rcsfst'};
    case 9, surf = {'lhv','rhv'}; %,'lcv','rcv'};
  end

  if ~job.output.surface && any( [job.output.ct.native job.output.ct.warped job.output.ct.dartel] )
    surf = {'lhv','rhv'}; %,'lcv','rcv'};
  end

  if job.output.surface>4 && job.output.surface~=9 % fast but not volume
    job.extopts.pbtres = max(0.8,min([((min(vx_vol)^3)/2)^(1/3) 1.0]));
  end

  % brain masking and correction of blood vessels
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*5/255;

  % surface creation and thickness estimation
  if 1
    %% using the Ymi map
    Ymix = Ymi .* cat_vol_morph(Yp0>0.5,'d');
    if job.extopts.pbtres==99
    % development block
      smeth = [3 1];
      sres  = [1 0.5];
      for smi = 1:numel(smeth)
        for sresi = 1:numel(sres)
          %%
          if smeth(smi)==1, pbtmethod = 'pbt2xf'; elseif smeth(smi)==3, pbtmethod = 'pbt3'; end

          cat_io_cprintf('blue',sprintf('\nPBT Test99 - surf_%s_%0.2f\n',pbtmethod,sres(sresi)));
          surf = {'lhfst'}; %,'lcfst','rhfst','rcfst'};

          [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Ymix,Yl1,Yp0/3,YMF,...
          struct('pbtmethod',pbtmethod,'interpV',sres(sresi),'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
          'verb',job.extopts.verb,'WMT',WMT));
        end
      end
    else
      pbtmethod = 'pbt2x';
      [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Ymix,Yl1,Yp0/3,YMF,...
        struct('pbtmethod',pbtmethod,'interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
        'verb',job.extopts.verb,'WMT',WMT));
    end
  else
  %% using the label image
    if job.extopts.pbtres==99
    % development block
      smeth = [3 1];
      sres  = [1 0.5];
      for smi = 1:numel(smeth)
        for sresi = 1:numel(sres)
          %%
          if smeth(smi)==1, pbtmethod = 'pbt2xf'; elseif smeth(smi)==3, pbtmethod = 'pbt3'; end

          cat_io_cprintf('blue',sprintf('\nPBT Test99 - surf_%s_%0.2f\n',pbtmethod,sres(sresi)));
          surf = {'lhfst'}; %,'lcfst','rhfst','rcfst'};

          [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Yp0/3,Yl1,Yp0/3,YMF,...
              struct('pbtmethod',pbtmethod,'interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
              'verb',job.extopts.verb,'WMT',WMT));
        end
      end
    else
      pbtmethod = 'pbt2x';
      [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Yp0/3,Yl1,Yp0/3,YMF,...
          struct('pbtmethod',pbtmethod,'interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
          'verb',job.extopts.verb,'WMT',WMT));
    end
  end

  % save Euler characteristics (absolute value)
  qa.subjectmeasures.EC_abs = EC;
  qa.subjectmeasures.defect_size = defect_size;

  if exist('S','var') && numel(fieldnames(S))==0 && isempty(Psurf)
    clear S Psurf;
  end

  % thickness map
  if isfield(job.output,'ct')
    cat_io_writenii(VT0,Yth1,mrifolder,'ct','cortical thickness map','uint16',...
      [0,0.0001],job.output.ct,trans,single(Ycls{1})/255,0.1);
  end

  if ~debug, clear Yp0; end

  cat_io_cmd('Surface and thickness estimation');
  fprintf('%5.0fs\n',etime(clock,stime));
  if ~debug; clear YMF; end
  if ~debug && ~job.output.ROI && job.output.surface, clear Yth1; end
else
  if ~debug; clear Ymi; end
end


%% ---------------------------------------------------------------------
%  ROI Partitioning
%  ---------------------------------------------------------------------
%  This part estimated indivudal measurements for different ROIs.
%  The ROIs are described in the CAT normalized space and there are to
%  ways to estimate them - (1) in subject space, and (2) in normalized
%  space. Estimation in normalized space is more direct an avoid further
%  transformations. The way over the subject space have the advantage
%  that indivdiual anatomical refinients are possible, but the this has
%  to be done and evalutated for each atlas.
%  ---------------------------------------------------------------------
if job.output.ROI
  stime = cat_io_cmd('ROI estimation');
  if job.extopts.verb, fprintf('\n'); end;

  for ai=1:size(FA,1)
    %%
    if ai==1 || any(VAvx_vol(ai,:)~=VAvx_vol(ai-1,:))
      % resampe data in atlas resolution for the first time or if the atlas resolution change

      Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*5/255;

      % map data to actual template space
      if ai==1
        stime2  = cat_io_cmd('  Data mapping to normalized atlas space','g5','', job.extopts.verb-1);
      else
        stime2  = cat_io_cmd('  Data mapping to normalized atlas space','g5','', job.extopts.verb-1,stime2);
      end
      transw      = trans.warped;                     % dartel/shooting deformation data
      transw.odim = VA(ai).dim;                       % adaption for atlas image size
      transw.ress = job.extopts.vox(1)./VAvx_vol(ai,:);  % adaption for atlas sampling resolution

      wYp0     = cat_vol_ROInorm(Yp0,transw,1,0,FA);
      wYcls    = cat_vol_ROInorm(Ycls,transw,1,1,FA);

      if exist('Ywmh','var')
        wYcls(7) = cat_vol_ctype(cat_vol_ROInorm({single(Ywmh)},transw,1,1,FA));
      end

      % correction for voxel size of the orignal image
      for ci=1:numel(wYcls)
        if ~isempty(wYcls{ci})
          wYcls{ci} = wYcls{ci} * prod(vx_vol);
        end
      end
      if debug
        %% check voxel size of an atlas map that labels the whole brain
        fprintf('%8s %8s %8s\n','GM','WM','CSF')
        fprintf('%8.2f %8.2f %8.2f\n', [cat_stat_nansum(wYcls{1}(:)),cat_stat_nansum(wYcls{2}(:)),cat_stat_nansum(wYcls{3}(:))]  / 1000);
        fprintf('%8.2f %8.2f %8.2f\n', [cat_stat_nansum(Ycls{1}(:)),cat_stat_nansum(Ycls{2}(:)),cat_stat_nansum(Ycls{3}(:))] * prod(vx_vol) / 1000 / 255);
      end

      %wYm      = cat_vol_ROInorm(Ym,transw,1,0,job.extopts.atlas);         % intensity

      if exist('Yth1','var')
        % ROI based thickness of all GM voxels per ROI
        Yth1x = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5) = nan;
        wYth1 = cat_vol_ROInorm(Yth1x,transw,1,0,FA);
        wYth1(wYth1==0) = nan;
        clear Yth1x;
      end
    end

    % ds('l2','',1.5,wYm,round(wYp0),wYm,single(wYa)/50 .* (wYp0<2.5),70)

    [px,atlas] = fileparts(FA{ai,1});
    stime2 = cat_io_cmd(sprintf('  ROI estimation of ''%s'' atlas',atlas),'g5','', job.extopts.verb-1,stime2);

    % map atlas to actual template space
    transa      = trans.warped;
    transa.M1   = VA(ai).mat;
    transa.odim = transw.odim;
    wYa   = cat_vol_ROInorm([],transa,ai,0,FA);


    %% extract ROI data
    csv   = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],FA{ai,3},FA);  % volume

    % thickness
    if exist('Yth1','var')
    % For thickness we want to avoid values in noncortical regions such as
    % the ventricles or regions with relative low GM volume or high CSF volume.
      csv  = cat_vol_ROIestimate(wYp0,wYa,wYth1,ai,'T',csv,{'gm'},job.extopts.atlas); %.*wYmim
      % correct for ventricular regions that use the 'Ven' keyword.
      ven  = find(cellfun('isempty',strfind( csv(:,2) , 'Ven'))==0);
      csv(ven,end) = {nan};
      % correct for regions with relative low GM (<10%) or high CSF volume (>50%).
      csvf = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],{'csf','gm','wm'},FA);
      vola = [nan,nan,nan;cell2mat(csvf(2:end,3:end))];
      volr = vola ./ repmat(sum(vola,2),1,3);
      csv(volr(:,2)<0.1 | volr(:,2)>0.5,end) = {nan};
    end

    % xml-export one file for all (this is a structure)
    ROI.(atlas) = csv;

  end

  % write results
  catROI = cat_roi_fun('csvtab2xmlroi',ROI);
  cat_io_xml(fullfile(pth,labelfolder,['catROI_' nam '.xml']),catROI,'write');

  cat_io_cmd(' ','g5','',job.extopts.verb,stime2);
  cat_io_cmd('','n','',1,stime);
end
if ~debug, clear wYp0 wYcls wYv trans; end


%%  estimate volumes and TIV
qa.subjectmeasures.vol_abs_CGW = [
  prod(vx_vol)/1000/255 .* sum(Ycls{3}(:)), ... CSF
  prod(vx_vol)/1000/255 .* sum(Ycls{1}(:)), ... GM
  prod(vx_vol)/1000/255 .* sum(Ycls{2}(:)) 0 0]; %, ... WM
if numel(Ycls)>6, qa.subjectmeasures.vol_abs_CGW(4) = prod(vx_vol)/1000/255 .* sum(Ycls{7}(:)); end % WMHs qa.subjectmeasures.WMH_abs
if numel(Ycls)>7, qa.subjectmeasures.vol_abs_CGW(5) = prod(vx_vol)/1000/255 .* sum(Ycls{8}(:)); end % SL
qa.subjectmeasures.vol_TIV     =  sum(qa.subjectmeasures.vol_abs_CGW);
qa.subjectmeasures.vol_rel_CGW =  qa.subjectmeasures.vol_abs_CGW ./ qa.subjectmeasures.vol_TIV;
if ~debug, clear Ycls; end

%% ---------------------------------------------------------------------
%  XML-report and Quality Assurance
%  ---------------------------------------------------------------------
stime = cat_io_cmd('Quality check'); job.stime = stime;
Yp0   = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*5; Yp0(Yp0>3.1) = nan;  %qa2=qa;
qa    = cat_vol_qa('cat12',Yp0,fname0,Ym,res,cat_warnings,job.extopts.species, ...
          struct('write_csv',0,'write_xml',1,'method','cat12','job',job,'qa',qa));

% WMH updates? ... has to be done within cat_vol_qa?!
%qa.subjectmeasures.vol_abs_CGW(2) = qa.subjectmeasures.vol_abs_CGW(2) - qa2.subjectmeasures.WMH_abs;
%qa.subjectmeasures.vol_abs_CGW(4) = qa2.subjectmeasures.WMH_abs;
%qa.subjectmeasures.vol_rel_CGW    = qa.subjectmeasures.vol_abs_CGW ./ qa.subjectmeasures.vol_TIV;
if job.output.surface && exist('S','var')
  % metadata
  if isfield(S,'lh') && isfield(S.lh,'th1'), th=S.lh.th1; else th=[]; end;
  if isfield(S,'rh') && isfield(S.rh,'th1'), th=[th; S.rh.th1]; end
  qa.subjectmeasures.dist_thickness{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th;
  if job.extopts.expertgui>1
    if isfield(S,'lh') && isfield(S.lh,'th2'), th=S.lh.th2; else th=[]; end;
    if isfield(S,'rh') && isfield(S.lh,'th2'), th=[th; S.rh.th2]; end
    qa.subjectmeasures.dist_gyruswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th;
    if isfield(S,'lh') && isfield(S.lh,'th3'), th=S.lh.th3; else th=[]; end;
    if isfield(S,'rh') && isfield(S.lh,'th3'), th=[th; S.rh.th3]; end
    qa.subjectmeasures.dist_sulcuswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th;
  end

  qam = cat_stat_marks('eval',job.cati,qa,'cat12');

  cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
    ... 'subjectratings',qam.subjectmeasures, ... not ready
    'subjectmeasures',qa.subjectmeasures,'ppe',res.ppe),'write+'); % here we have to use the write+!
elseif exist('Yth1','var')
  qa.subjectmeasures.dist_thickness{1} = [cat_stat_nanmean(Yth1(Yth1(:)>1)) cat_stat_nanstd(Yth1(Yth1(:)>1))];

  qam = cat_stat_marks('eval',job.cati,qa,'cat12');

  cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
    ... 'subjectratings',qam.subjectmeasures, ... not ready
    'subjectmeasures',qa.subjectmeasures,'ppe',res.ppe),'write+'); % here we have to use the write+!
end
clear Yo Yp0 qas Yth1;
fprintf('%5.0fs\n',etime(clock,stime));


if job.extopts.print
  %% display and print result if possible
  %  ---------------------------------------------------------------------
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);

    warning off; %#ok<WNOFF> % there is a div by 0 warning in spm_orthviews in linux


    %mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
    marks2str = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
    mark2rps    = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
    grades      = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
    mark2grad   = @(mark) grades{min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3)))};


    % CAT GUI parameter:
    % --------------------------------------------------------------------
    str = [];

    % 1 line: Matlab, SPM12, CAT12 version number and GUI and experimental mode
    str = [str struct('name', 'Version: Matlab / SPM12 / CAT12:','value',...
      sprintf('%s / %s / %s (%s)',qa.software.version_matlab,qa.software.version_spm,qa.software.version_cat,qa.software.revision_cat))];
    if     job.extopts.expertgui==1, str(end).value = [str(end).value '\bf\color[rgb]{0 0.2 1}e'];
    elseif job.extopts.expertgui==2, str(end).value = [str(end).value '\bf\color[rgb]{0 0.2 1}d'];
    end
    if job.extopts.experimental, str(end).value = [str(end).value '\bf\color[rgb]{0 0.2 1}x']; end

    % 2 lines: TPM, Template, Normalization method with voxel size
    str = [str struct('name', 'Tissue Probability Map:','value',strrep(spm_str_manip(res.tpm(1).fname,'k40d'),'_','\_'))];
    if do_dartel
      if job.extopts.regstr==0 % Dartel
        str = [str struct('name', 'Dartel Registration to: ',...
                          'value',strrep(spm_str_manip(job.extopts.darteltpm{1},'k40d'),'_','\_'))];
      elseif job.extopts.regstr==4 % Dartel
        str = [str struct('name', 'Shooting Registration to: ',...
                          'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
      else
        if job.extopts.expertgui==0
          str = [str struct('name','Optimized Shooting Registration to:',...
                            'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
        else
          str = [str struct('name', sprintf('Optimized Shooting Registration (regstr:%s) to :',sprintf('%g ',job.extopts.regstr)),...
                            'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
        end
      end
    end

    % 1 line 1: Affreg
    str = [str struct('name', 'affreg:','value',sprintf('%s',job.opts.affreg))];
    % 1 line 2: APP
    APPstr = {'none','light','full','APPi','APPf','animal'}; APPstr{1071} = 'APP';
    str(end).name  = [str(end).name(1:end-1) ' / APP '];
    str(end).value = [str(end).value sprintf(' / %s',APPstr{job.extopts.APP+1})];
    % 1 line 3: biasstr / biasreg+biasfwhm
    if job.opts.biasstr>0
      biasstr = {'ultralight','light','medium','strong','heavy'};
      str(end).name  = [str(end).name(1:end-1) ' / biasstr '];
      str(end).value = [str(end).value sprintf(' / %s',biasstr{round(job.opts.biasstr*4)+1})];
      if job.extopts.expertgui % add the value
        str(end).value = [str(end).value sprintf('(%0.2f;breg:%0.2f;bfwhm:%0.2f)',job.opts.biasstr,job.opts.biasreg,job.opts.biasfwhm)];
      end
    else
      str(end).name  = [str(end).name(1:end-1) ' / biasreg / biasfwhm'];
      str(end).value = [str(end).value sprintf(' / %0.2f / %0.2f',job.opts.biasreg,job.opts.biasfwhm)];
    end


    % 1 line: adaptive noise parameter ( MRFstr + SANLM + NCstr )
    defstr  = {'none','ultralight','light','medium','strong','heavy',... sanlm vs. isarnlm
               'ultralight+','ultralight+','light+','medium+','strong+','heavy+'};
    defstrm = @(x) defstr{ round(max(0,min(2,x))*4) + 1 + (x>0) + (x>1)};
    str = [str struct('name', 'Noise reduction:','value','')];
    if job.extopts.NCstr==0
      if job.extopts.mrf==0
        str(end).value = 'no noise correction';
      else
        if job.extopts.expertgui==0
          str(end).value = 'MRF';
        else
          str(end).value = sprintf('MRF(%0.2f)',job.extopts.mrf);
        end
      end
    else
      str(end).value = sprintf('SANLM(%s)',NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')});
    end

    if job.extopts.NCstr~=0 && job.extopts.mrf
      if job.extopts.expertgui==0
        str(end).value = '+MRF';
      else
        str(end).value = sprintf('+MRF(%0.2f)',job.extopts.mrf);
      end
    end


    % 1 line(s): LASstr / GCUTstr / CLEANUPstr
    str(end).name  = 'LASstr / GCUTstr / CLEANUPstr:';
    if ~job.extopts.expertgui
      str(end).value = sprintf('%s / %s / %s',defstrm(job.extopts.LASstr),defstrm(job.extopts.gcutstr),defstrm(job.extopts.cleanupstr));
    else
      str(end).value = sprintf('%s(%0.2f) / %s(%0.2f) / %s(%0.2f)',defstrm(job.extopts.LASstr),job.extopts.LASstr,...
        defstrm(job.extopts.gcutstr),job.extopts.gcutstr,defstrm(job.extopts.cleanupstr),job.extopts.cleanupstr);
    end
    restype = char(fieldnames(job.extopts.restypes));
    if job.extopts.expertgui
      str = [str struct('name', 'WMHC / SLC / restype:','value',...
             sprintf('%d / %d / %s',...
            job.extopts.WMHC,job.extopts.SLC,restype))];
    else
      str = [str struct('name', 'restype:','value',sprintf('%s',restype))];
    end
    if ~strcmp('native',restype)
      str(end).value = [str(end).value sprintf('(%0.2f %0.2f)',job.extopts.restypes.(restype))];
    end;

    % line 8: surfae parameter
    if job.output.surface
      str = [str struct('name', 'Voxel resolution (original > internal > PBT; vox):',...
             'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s > %4.2f mm%s; %4.2f mm  ', ...
             qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179),job.extopts.pbtres,char(179),job.extopts.vox(1)))];
    else
      str = [str struct('name', 'Voxel resolution (original > intern; vox):',...
             'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s; %4.2f mm', ...
             qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179),job.extopts.vox(1)))];
    end
    % str = [str struct('name', 'Norm. voxel size:','value',sprintf('%0.2f mm',job.extopts.vox))]; % does not work yet


    % Image Quality measures:
    % --------------------------------------------------------------------
    str2 =       struct('name', '\bfImage and Preprocessing Quality:','value','');
    str2 = [str2 struct('name',' Resolution:','value',marks2str(qa.qualityratings.res_RMS,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.res_RMS),mark2grad(qa.qualityratings.res_RMS))))];
    str2 = [str2 struct('name',' Noise:','value',marks2str(qa.qualityratings.NCR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.NCR),mark2grad(qa.qualityratings.NCR))))];
    str2 = [str2 struct('name',' Bias:','value',marks2str(qa.qualityratings.ICR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.ICR),mark2grad(qa.qualityratings.ICR))))]; % not important and more confusing
    str2 = [str2 struct('name','\bf Weighted average (IQR):','value',marks2str(qa.qualityratings.IQR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR))))];
    if isfield(qa.subjectmeasures,'EC_abs') && ~isempty(qa.subjectmeasures.EC_abs)
      str2 = [str2 struct('name',' Surface Euler number:','value',sprintf('%d', qa.subjectmeasures.EC_abs))];
    end
    if isfield(qa.subjectmeasures,'defect_size') && ~isempty(qa.subjectmeasures.defect_size)
      str2 = [str2 struct('name',' Size of topology defects:','value',sprintf('%d', qa.subjectmeasures.defect_size))];
    end

    % Subject Measures
    % --------------------------------------------------------------------
    % Volume measures

    % header
    str3 = struct('name', '\bfVolumes:','value',sprintf('%5s %5s %5s ','CSF','GM','WM'));
    if job.extopts.WMHC>1, str3(end).value = [str3(end).value sprintf('%5s ','WMH')]; end
    if job.extopts.SLC>0,  str3(end).value = [str3(end).value sprintf('%5s ','SL')];  end

    % absolute volumes
    str3 = [str3 struct('name', ' Absolute volume:','value',sprintf('%5.0f %5.0f %5.0f ', ...
            qa.subjectmeasures.vol_abs_CGW(1:3)))];
    if job.extopts.WMHC>1,  str3(end).value = [str3(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_abs_CGW(4))]; end;
    if job.extopts.SLC>0,   str3(end).value = [str3(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_abs_CGW(5))]; end;
    str3(end).value = [str3(end).value 'cm' char(179)];

    % relative volumes
    str3 = [str3 struct('name', ' Relative volume:','value',sprintf('%5.1f %5.1f %5.1f ', ...
            qa.subjectmeasures.vol_rel_CGW(1:3)*100))];
    if job.extopts.WMHC>1,  str3(end).value = [str3(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_rel_CGW(4)*100)]; end;
    if job.extopts.SLC>0,   str3(end).value = [str3(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_rel_CGW(5)*100)]; end;
    str3(end).value = [str3(end).value '%'];

    Ywmhrel = NS(Yl1,23);
    Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*5;
    WMH_rel    = 100*sum(Ywmhrel(:)) / sum(Yp0(:)>0.5);
    WMH_WM_rel = 100*sum(Ywmhrel(:)) / sum(Yp0(:)>2.5);
    if WMH_rel>3 || WMH_WM_rel>5
      str3(end).value = [str3(end).value sprintf('\\bf\\color[rgb]{0.8 0 0} WMHs!')];
    end
    clear Yp0;

    str3 = [str3 struct('name', ' TIV:','value', sprintf(['%0.0f cm' char(179)],qa.subjectmeasures.vol_TIV))];


    % Surface measures - Thickness, (Curvature, Depth, ...)
    if isfield(qa.subjectmeasures,'dist_thickness') && ~isempty(qa.subjectmeasures.dist_thickness)
      str3 = [str3 struct('name', '\bfThickness:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_thickness{1}(1),177,qa.subjectmeasures.dist_thickness{1}(2)))];
      if isfield(qa.subjectmeasures,'dist_gyruswidth') && ~isnan(qa.subjectmeasures.dist_gyruswidth{1}(1))
        str3 = [str3 struct('name', '\bfGyruswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
               qa.subjectmeasures.dist_gyruswidth{1}(1),177,qa.subjectmeasures.dist_gyruswidth{1}(2)))];
      end
      if isfield(qa.subjectmeasures,'dist_sulcuswidth') && ~isnan(qa.subjectmeasures.dist_sulcuswidth{1}(1))
        str3 = [str3 struct('name', '\bfSulcuswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
               qa.subjectmeasures.dist_sulcuswidth{1}(1),177,qa.subjectmeasures.dist_sulcuswidth{1}(2)))];
      end
    end

    % Preprocessing Time
    str2 = [str2 struct('name','\bfProcessing time:','value',sprintf('%02.0f:%02.0f min', ...
      floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60)))];

    % Warnings
    if numel(cat_warnings)>0 && job.extopts.expertgui>0
      str2 = [str2 struct('name', '','value','')];
      str2 = [str2 struct('name', '\bfWarnings:','value','')];
      for wi=1:numel(cat_warnings)
        shorter = cat_warnings(wi).identifier;
        % remove leading MATLAB, SPM or CAT elements
        dots    = max([min(strfind(shorter,'MATLAB')+7), ...
                       min(strfind(shorter,'SPM')+4), ...
                       min(strfind(shorter,'CAT')+4)]);
        if ~isempty(dots), shorter = shorter(dots:end); end
        % limit lenght of the string and replace critical character
        shorter = spm_str_manip(shorter,'l40');
        shorter = marks2str(4,shorter);
        shorter = strrep(shorter,'_','\_');
        str2    = [str2 struct('name',shorter,'value','')];  %#ok<AGROW>
      end
    end

    % adding one space for correct printing of bold fonts
    for si=1:numel(str)
      str(si).name   = [str(si).name  '  '];  str(si).value  = [str(si).value  '  '];
    end
    for si=1:numel(str2)
      str2(si).name  = [str2(si).name '  '];  str2(si).value = [str2(si).value '  '];
    end
    for si=1:numel(str3)
      str3(si).name  = [str3(si).name '  '];  str3(si).value = [str3(si).value '  '];
    end


    %
    fg = spm_figure('FindWin','Graphics');
    set(0,'CurrentFigure',fg)
	%%% Hacking for xASL use here, assume always parallel option on
    if isempty(fg)
      %if job.nproc, fg = spm_figure('Create','Graphics','visible','off'); else fg = spm_figure('Create','Graphics'); end;
	  fg = spm_figure('Create','Graphics','visible','off');
    else
      if job.nproc, set(fg,'visible','off'); end
    end
    set(fg,'windowstyle','normal');
    spm_figure('Clear','Graphics');
    switch computer
      case {'PCWIN','PCWIN64'}, fontsize = 8;
      case {'GLNXA','GLNXA64'}, fontsize = 8;
      case {'MACI','MACI64'},   fontsize = 9.5;
      otherwise,                fontsize = 9.5;
    end
    ax=axes('Position',[0.01 0.75 0.98 0.24],'Visible','off','Parent',fg);

    text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image0(1).fname,'k60d') '       '],...
      'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','none','Parent',ax);

    cm = job.extopts.colormap;

    % check colormap name
    switch lower(cm)
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter',...
          'gray','bone','copper','pink','bcgwhw','bcgwhn'}
      otherwise
        cat_io_cprintf(job.color.warning,'WARNING:Unknown Colormap - use default.\n');
        cm = 'gray';
    end

    % SPM_orthviews seams to allow only 60 values
    % It further requires a modified colormaps with lower values that the
    % colorscale and small adaption for the values.
    surfcolors = 128;
    switch lower(cm)
      case {'bcgwhw','bcgwhn'} % cat colormaps with larger range
        ytick       = [1,5:5:60];
        yticklabel  = {' BG',' ',' CSF',' CGM',' GM',' GWM',' WM',' ',' ',' ',' ',' ',' BV / HD '};
        yticklabelo = {' BG',' ','    ','    ','   ','     ',' avg WM  ',' ',' ',' ',' ',' ',' BV / HD '};
        yticklabeli = {' BG',' ','    ','    ','   ','  ','  ',' ',' ',' ',' ',' ',' BV / HD '};
        %colormap(cat_io_colormaps(cm,60));
        cmap = [cat_io_colormaps([cm 'ov'],60);flipud(cat_io_colormaps([cm 'ov'],60));jet(surfcolors)];
        cmmax = 2;
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
        ytick       = [1 20 40 60];
        yticklabel  = {' BG',' CSF',' GM',' WM'};
        yticklabelo = {' BG','    ','   ',' WM'};
        yticklabeli = {' BG','    ','   ','   '};
        cmap = [eval(sprintf('%s(60)',cm));flipud(eval(sprintf('%s(60)',cm)));jet(surfcolors)];
        cmmax = 1;
    end
    colormap(cmap);
    spm_orthviews('Redraw');

    htext = zeros(5,2,2);
    for i=1:size(str,2)   % main parameter
      htext(1,i,1) = text(0.01,0.98-(0.055*i), str(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
      htext(1,i,2) = text(0.51,0.98-(0.055*i), str(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end
    for i=1:size(str2,2)  % qa-measurements
      htext(2,i,1) = text(0.01,0.40-(0.055*i), str2(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
      htext(2,i,2) = text(0.25,0.40-(0.055*i), str2(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end
    % qa-scala
    %htext(5,1,1) = text(0.01,0.45-(0.055*(i+2)),str4(1).name,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    for i=1:size(str3,2)  % subject-measurements
      htext(3,i,1) = text(0.51,0.40-(0.055*i), str3(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
      htext(3,i,2) = text(0.70,0.40-(0.055*i), str3(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end



    pos = [0.01 0.38 0.48 0.36; 0.51 0.38 0.48 0.36; ...
           0.01 0.01 0.48 0.36; 0.51 0.01 0.48 0.36];
    spm_orthviews('Reset');




    % BB box is not optimal for all images
    disptype = 'affine';
    switch disptype
      case 'affine'
        dispmat = res.Affine;
        spm_orthviews('BB', job.extopts.bb*0.95 );
      case 'ridid'
        % this does not work so good... AC has a little offset ...
        aff = spm_imatrix(res.Affine);  scale = aff(7:9);
        spm_orthviews('BB', job.extopts.bb ./ mean(scale));
        dispmat = R;
    end


    % Yo - original image in original space
    % using of SPM peak values didn't work in some cases (5-10%), so we have to load the image and estimate the WM intensity
    try %#ok<TRYNC>
      Yo  = single(VT.private.dat(:,:,:));
    end
    Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*5;
    if exist('Yo','var')
      if job.inv_weighting
        WMth = min([...
          cat_stat_nanmedian(Yo(Yp0(:)>0.8 & Yp0(:)<1.2))*2,...
          cat_stat_nanmedian(Yo(Yp0(:)>1.8 & Yp0(:)<2.2))*1.5,...
          ]);
        T1txt = '*.nii (Original PD/T2)';
      else
        WMth = cat_stat_nanmedian(Yo(Yp0(:)>2.8 & Yp0(:)<3.2)); clear Yo;
        T1txt = '*.nii (Original T1)';
      end
      if ~debug, clear Yo; end

      if isfield(res,'spmpp')
        VT0x = res.image0(1);
      else
        VT0x = VT0;
      end
      VT0x.mat = dispmat * VT0x.mat;
      hho = spm_orthviews('Image',VT0x,pos(1,:));
      spm_orthviews('Caption',hho,{T1txt},'FontSize',fontsize,'FontWeight','Bold');
      spm_orthviews('window',hho,[0 WMth*cmmax]); caxis([0,2]);
      cc{1} = axes('Position',[pos(1,1) + 0.30 0.37 0.02 0.15],'Parent',fg);
      try, image(cc{1},(60:-1:1)'); end

      if job.inv_weighting
        set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabeli),'XTickLabel','','XTick',[],'TickLength',[0 0],...
          'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
      else
        set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabelo),'XTickLabel','','XTick',[],'TickLength',[0 0],...
          'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
      end
    else
      cat_io_cprintf('warn','WARNING: Can''t display original file "%s"!\n',VT.fname);
    end


    % Ym - normalized image in original space
    if ~isfield(res,'spmpp')
      %%
      Vm        = res.image(1);
      Vm.fname  = '';
      Vm.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
      Vm.dat(:,:,:) = single(Ym);
      Vm.pinfo  = repmat([1;0],1,size(Ym,3));
      Vm.mat    = dispmat * Vm.mat;
      hhm = spm_orthviews('Image',Vm,pos(2,:));
      spm_orthviews('Caption',hhm,{'m*.nii (Int. Norm.)'},'FontSize',fontsize,'FontWeight','Bold');
      spm_orthviews('window',hhm,[0 cmmax]); caxis([0,2]);
      cc{2} = axes('Position',[pos(2,1) + 0.30 0.37 0.02 0.15],'Parent',fg);
      try, image(cc{2},(60:-1:1)'); end
      set(cc{2},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
        'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
    end

    % Yo - segmentation in original space
    VO        = res.image(1);
    VO.fname  = '';
    VO.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
    VO.dat(:,:,:) = single(Yp0/3);
    VO.pinfo  = repmat([1;0],1,size(Yp0,3));
    VO.mat    = dispmat * VO.mat;
    hhp0 = spm_orthviews('Image',VO,pos(3,:)); if ~debug, clear Yp0; end
    spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hhp0,[0 cmmax]); caxis([0,2]);
    cc{3} = axes('Position',[pos(3,1) + 0.30 0.02 0.02 0.15],'Parent',fg);
    try, image(cc{3},(60:-1:1)'); end
    set(cc{3},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
      'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
    spm_orthviews('Reposition',[0 0 0]);


    % surface
    if job.extopts.print>1
      if exist('Psurf','var')
        try
          hCS = subplot('Position',[0.50 0.05 0.55 0.30],'visible','off');
          hSD = cat_surf_display(struct('data',Psurf(1).Pthick,'readsurf',0,'expert',2,...
            'multisurf',job.output.surface,'view','s',...
            'parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));
          colormap(cmap);  set(hSD{1}.colourbar,'visible','off');
          cc{3} = axes('Position',[0.63 0.02 0.3 0.01],'Parent',fg); image(cc{3},(121:1:120+surfcolors));
          set(cc{3},'XTick',1:(surfcolors-1)/6:surfcolors,'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
            'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
        catch
          cat_io_cprintf('warn','WARNING: Can''t display surface!\n',VT.fname);
        end
      end
    end



    %% print subject report file as standard PDF/PNG/... file
    job.imgprint.type  = 'pdf';
    job.imgprint.dpi   = 600;
    job.imgprint.fdpi  = @(x) ['-r' num2str(x)];
    job.imgprint.ftype = @(x) ['-d' num2str(x)];
    job.imgprint.fname     = fullfile(pth,reportfolder,['catreport_' nam '.' job.imgprint.type]);

    fgold.PaperPositionMode = get(fg,'PaperPositionMode');
    fgold.PaperPosition     = get(fg,'PaperPosition');
    fgold.resize            = get(fg,'resize');

    % it is necessary to change some figure properties especialy the fontsizes
    set(fg,'PaperPositionMode','auto','resize','on','PaperPosition',[0 0 1 1]);
    for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize*0.8); end; end
    for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize*0.8); end;
    warning off
    print(fg, job.imgprint.ftype(job.imgprint.type), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fname);
    warning on
    for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize); end; end
    for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize); end;
    set(fg,'PaperPositionMode',fgold.PaperPositionMode,'resize',fgold.resize,'PaperPosition',fgold.PaperPosition);
    try
      fprintf('Print ''Graphics'' figure to: \n  %s\n',job.imgprint.fname);% windows error?
    end

    %% reset colormap to the simple SPM like gray60 colormap
    if exist('hSD','var')
      % if there is a surface than we have to use the gray colormap also here
      % because the colorbar change!
      try %#ok<TRYNC>
        cat_surf_render2('ColourMap',hSD{1}.axis,gray(128));
        cat_surf_render2('Clim',hSD{1}.axis,[0 6]);
        axes(cc{3}); image(cc{3},0:60);
        set(cc{3},'XTick',max(1,0:10:60),'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
          'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
      end
    end

    WMfactor = 4/3; cmap = gray(60); colormap(cmap); caxis([0,numel(cmap)]);

    % new colorscale
    if job.inv_weighting
      % if exist('hho' ,'var'), spm_orthviews('window',hho ,[0 T3th(3)*WMfactor]); end
    else
      if exist('hho' ,'var'), spm_orthviews('window',hho ,[0 T3th(3)*WMfactor]); end
    end
    if exist('hhm' ,'var'), spm_orthviews('window',hhm ,[0 WMfactor]); end
    if exist('hhp0','var'), spm_orthviews('window',hhp0,[0 WMfactor]); end

    warning on;  %#ok<WNON>


    %% command window output
    fprintf('\n%s',repmat('-',1,72));
    fprintf(1,'\nCAT preprocessing takes %0.0f minute(s) and %0.0f second(s).\n', ...
      floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60));
    cat_io_cprintf(color(QMC,qa.qualityratings.IQR), sprintf('Image Quality Rating (IQR):  %5.2f%%%% (%s)\n',...
      mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR)));

    % print subfolders
    if job.extopts.subfolders
      fprintf('Segmentations are saved in %s%s%s\n',pth,filesep,'mri');
      fprintf('Reports are saved in %s%s%s\n',pth,filesep,'report');
      if job.output.ROI
        fprintf('Labels are saved in %s%s%s\n',pth,filesep,'label');
      end
      if job.output.surface && exist('S','var')
        fprintf('Surface measurements are saved in %s%s%s\n',pth,filesep,'surf');
      end
    end

    fprintf('%s\n\n',repmat('-',1,72));

    % finish diary entry of "../report/cmdln_*.txt"
    % read diagry and add the command-line output to the *.xml and *.mat file
    % diary off; % ExploreASL fix - disable CAT12 outputs as we write everything to ExploreASL log file.
    try
      fid=fopen(res.catlog);
      txt=fread(fid,200000,'uint8=>char');
      fclose(fid);
      txt2=textscan(txt,'%s','Delimiter','');
      cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
        'catlog',txt2),'write+'); % here we have to use the write+!
    end

    spm_progress_bar('Clear');
    clear C c Ymi Ymf Ym

end
%%
return;
%=======================================================================

%=======================================================================
function wYv = cat_vol_ROInorm(Yv,warped,ai,mod,FA)
% ----------------------------------------------------------------------
% normalized space:
% ----------------------------------------------------------------------
% for normalized space no further adaptions are available, but
% a masking based on the tissue map can be used
% ----------------------------------------------------------------------

  % load mask (and complete undefined parts)
  if isempty(Yv)
    % no input - load atlas

    if ~exist(FA{ai,1},'file')
      error('cat:cat_main:missAtlas','Miss cat atlas-file ''%s''!',FA{ai,1});
    end
    % try multiple times, because of read error in parallel processing
    for i=1:5
      try
        wVv = spm_vol(FA{ai,1});
        wYv = spm_read_vols(wVv);
        break
      catch
        pause(0.5)
      end
    end

    % resample atlas, if the atlas resolution differs from the actual template resolution
    if 0  %wVv.mat(1) ~= warped.M1(1)
      wVv2 = wVv; wVv2.mat = warped.M1; wVv2.dim = warped.odim;
      [t,wYv] = cat_vol_imcalc(wVv,wVv2,'i1',struct('interp',0,'verb',0));
    end
    wYv = cat_vol_ctype(wYv,wVv(1).private.dat.dtype);
  else
    % map image to atlas space
    for yi=1:numel(warped.ress), warped.y(:,:,:,yi) = warped.y(:,:,:,yi) * warped.ress(yi); end
    if mod==0
      old = 0;

      Vlai = spm_vol(FA{ai,1});
      vx_vol_Vlai   = sqrt(sum(Vlai.mat(1:3,1:3).^2));
      vx_vol_Vdef   = sqrt(sum(warped.M1(1:3,1:3).^2));

      if old % this did not work for the thickness map?
        [wYv,w] = spm_diffeo('push',Yv,warped.y,warped.odim(1:3)); spm_field('boundary',1);
        wYv = spm_field(w,wYv,[sqrt(sum(warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]);
      elseif any( vx_vol_Vdef > vx_vol_Vlai*2)
        %% increase resolution
        fc   = ceil(vx_vol_Vdef / vx_vol_Vlai);
        ddim = (size(warped.y)-1)*fc+1; ddim(4)=[];
        eyev = eye(4); eyev(1:end-1) = eyev(1:end-1) * 1/fc;
        Yy   = zeros([ddim 3],'single');
        for k1=1:3
          for i=1:ddim(3),
            Yy(:,:,i,k1) = single(spm_slice_vol(warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]),ddim(1:2),[1,NaN])); % adapt for res
          end
        end
        Yvi  = zeros(ddim,'single');
        for i=1:ddim(3),
          Yvi(:,:,i) = single(spm_slice_vol(Yv(:,:,:),eyev*spm_matrix([0 0 i]),ddim(1:2),[1,NaN])); % adapt for res
        end

        [wYv,w]  = spm_diffeo('push',Yvi,Yy,warped.odim(1:3));
        % divide by jacdet to get unmodulated data
        wYv = wYv./(w+0.001);


      else
        [wYv,w]  = spm_diffeo('push',Yv,warped.y,warped.odim(1:3));
        % divide by jacdet to get unmodulated data
        wYv = wYv./(w+0.001);
      end
    elseif mod==1 && iscell(Yv) % tissue case
      nicemapping = 1;
      if nicemapping
        % Modulation using spm_diffeo and push introduces aliasing artefacts,
        % thus we use the def2det function of the inverted deformations to obtain the old and
        % in my view a more appropriate jacobian determinant
        % The 2nd reason to use the old modulation is compatibility with cat_vol_defs.m
        Yy = spm_diffeo('invdef',warped.y,warped.odim,eye(4),warped.M0);
        w  = spm_diffeo('def2det',Yy)/det(warped.M0(1:3,1:3)); clear Yy;
        % ensure that jacobian det is positive (no clue why some times the sign is switched)
        if mean(w(~isnan(w))) < 0, w = -w; end
        w(:,:,[1 end]) = NaN; w(:,[1 end],:) = NaN; w([1 end],:,:) = NaN;
      end

      wYv = cell(1,numel(Yv));
      for i=1:numel(Yv)
        if nicemapping
          [wYv{i},w2] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3));
          % divide by jacdet to get unmodulated data
          wYv{i} = wYv{i}./(w2+0.001);
          wYv{i} = wYv{i} .* w;
        else
          % simple push
          [wYv{i},w] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3));
        end
      end
    else
      error('unknown case');
    end
    %spm_smooth(wYv,wYv,[1 1 1]);
  end


return
%=======================================================================

%=======================================================================
function csv = cat_vol_ROIestimate(Yp0,Ya,Yv,ai,name,csv,tissue,FA)
% ----------------------------------------------------------------------
% estimate values
% ----------------------------------------------------------------------


% load atlas-csv-file

  [pp,ff] = fileparts(FA{ai,1});
  csvf = fullfile(pp,[ff '.csv']);

  if isempty(csv)
    if exist(csvf,'file')
      csv = cat_io_csv(csvf,'','',struct('delimiter',';'));
    else
      csv = [num2cell((1:max(Ya(:)))') ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')]) ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')])];
    end

    % remove empty rows and prepare structure names
    if size(csv,2)>2, csv(:,3:end)=[]; end
    for ri=size(csv,1):-1:1
      if isempty(csv{ri,1}) || isempty(csv{ri,2});
        csv(ri,:)=[];
      elseif csv{ri,1}==0
        csv(ri,:)=[];
      end
    end
  end
  name = genvarname(strrep(strrep(name,'-','_'),' ','_'));



  %% volume case
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
  % other maps with masks
  for ti=1:numel(tissue)
    switch name(1)
      case 'V' % volume
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
        for ri=2:size(csv,1)
          switch lower(tissue{ti})
            case 'csf',   Ymm=single(Yv{3}) .* single(Ya==csv{ri,1});
            case 'gm',    Ymm=single(Yv{1}) .* single(Ya==csv{ri,1});
            case 'wm',    Ymm=single(Yv{2}) .* single(Ya==csv{ri,1});
            case 'wmh',   Ymm=single(Yv{2}) .* single(Ya==csv{ri,1});
            case 'brain', Ymm=single(Yv{1} + Yv{2} + Yv{3}) .* single(Ya==csv{ri,1});
            case '',      Ymm=single(Ya==csv{ri,1});
          end
          csv{ri,end} = 1/1000 * cat_stat_nansum(Ymm(:));
        end
      case 'c'
        return
      otherwise %
        csv{1,end+1} = strrep([name tissue{ti}],'Tgm','ct');  %#ok<AGROW>
        switch lower(tissue{ti})
          case 'csf',   Ymm=Yp0toC(Yp0,1);
          case 'gm',    Ymm=Yp0toC(Yp0,2);
          case 'wm',    Ymm=Yp0toC(Yp0,3);
          case 'wmh',   Ymm=Yp0toC(Yp0,4);
          case 'brain', Ymm=Yp0>0.5;
          case '',      Ymm=true(size(Yp0));
        end
        for ri=2:size(csv,1)
          csv{ri,end} = cat_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
        end
    end
  end

return
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================

%=======================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = bsxfun(@minus,cr,mn(:,k)')*inv(chol(vr(:,:,k))); %#ok<MINV>
    p(:,k) = amp*exp(-0.5*sum(d.*d,2)) + eps;
end
return;
%==========================================================================

%==========================================================================
% function [P] = clean_gwc(P,level)
%==========================================================================
function [P] = clean_gwc(P,level)
if nargin<4, level = 1; end

b    = P(:,:,:,2);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level>1, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter  = 32;
niter2 = 32;
spm_progress_bar('Init',niter+niter2,'Extracting Brain','Iterations completed');
fprintf('Clean Up 1:   ');
for j=1:niter
	xASL_TrackProgress(j,niter);
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = single(P(:,:,i,1));
        wp       = single(P(:,:,i,2));
        bp       = single(b(:,:,i))/255;
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = cat_vol_ctype(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
    spm_progress_bar('Set',j);
end

% Also clean up the CSF.
fprintf('\n');
fprintf('Also clean up the CSF:   ');
if niter2 > 0,
    c = b;
    for j=1:niter2
		xASL_TrackProgress(j,niter2);
        for i=1:size(b,3)
            gp       = single(P(:,:,i,1));
            wp       = single(P(:,:,i,2));
            cp       = single(P(:,:,i,3));
            bp       = single(c(:,:,i))/255;
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = cat_vol_ctype(round(bp));
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j+niter);
    end
end

th = 0.05;
fprintf('\n');
fprintf('Some slice-work:   ');
for i=1:size(b,3)
	xASL_TrackProgress(i,size(b,3));
    slices = cell(1,size(P,4));
    for k1=1:size(P,4),
        slices{k1} = single(P(:,:,i,k1))/255;
    end
    bp        = single(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{1}+slices{2}))>th;
    slices{1} = slices{1}.*bp;
    slices{2} = slices{2}.*bp;

    if niter2>0,
        cp        = single(c(:,:,i))/255;
        cp        = ((cp>th).*(slices{1}+slices{2}+slices{3}))>th;
        slices{3} = slices{3}.*cp;
    end
    if numel(slices)>=5
      slices{5} = slices{5}+1e-4; % Add something to the soft tissue class
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4),
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4),
        P(:,:,i,k1) = cat_vol_ctype(round(slices{k1}./tot*255));
    end
end

spm_progress_bar('Clear');
return;
%==========================================================================
