function [Yth1,S,Psurf,EC,defect_size,res] = cat_surf_createCS(V,V0,Ym,Ya,YMF,opt)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth1,S,Psurf,EC]=cat_surf_createCS(V,V0,Ym,Ya,YMF,opt)
%
% Yth1  = thickness map
% S     = structure with surfaces, like the left hemishere, that contains
%        vertices, faces, GM thickness (th1), and the transformation to
%        map to nifti space (vmat) and back (vmati).
% Psurf = name of surface files
% EC    = Euler characteristics
% defect_size = size of topology defects
% V     = spm_vol-structure of internally interpolated image
% V0    = spm_vol-structure of original image
% Ym    = the (local) intensity, noise, and bias corrected T1 image
% Ya    = the atlas map with the ROIs for left and right hemispheres
%        (this is generated with cat_vol_partvol)
% Yp0   = label image for surface deformation
% YMF   = a logical map with the area that has to be filled
%        (this is generated with cat_vol_partvol)
%   
% opt.surf       = {'lh','rh'[,'lc','rc']} - side
%    .reduceCS   = 100000 - number of faces
%
% Options set by cat_defaults.m
%    .interpV    = 0.5    - mm-resolution for thickness estimation
% 
% Here we used the intensity normalized image Ym, rather that the Yp0
% image, because it has more information about sulci that we need 
% especially for asymmetrical sulci.
% Furthermore, all non-cortical regions and blood vessels were removed 
% (for left and right surface). Blood vessels (with high contrast) can 
% lead to strong error in the topology correction. Higher resolution 
% also helps to reduce artifacts.
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_surf_createCS.m 1561 2020-02-04 15:49:34Z gaser $ 

% Turn off gifti data warning in gifti/subsref (line 45)
%   Warning: A value of class "int32" was indexed with no subscripts specified. 
%            Currently the result of this operation is the indexed value itself, 
%            but in a future release, it will be an error. 
warning('off','MATLAB:subscripting:noSubscriptsSpecified');
cstime = clock;

  % variables to tranfer from MATLAB to image coordinates used by loadSurf and saveSurf subfunctions
  global vmat vmati mati
  
  % surface evaluation paramter 
  res = struct('euler_characteristic',nan,'defect_size_promile',nan,'lh',struct(),'rh',struct()); 

%#ok<*AGROW>
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  % set defaults
  vx_vol  = sqrt(sum(V.mat(1:3,1:3).^2));   % further interpolation based on internal resolution 
  vx_vol0 = sqrt(sum(V0.mat(1:3,1:3).^2));  % final surface resolution based on original image resolution
  if ~exist('opt','var'), opt=struct(); end
  def.verb      = max(2,2 + cat_get_defaults('extopts.expertgui')); 
  def.surf      = {'lh','rh'};
  
  % reducepatch has some issues with self intersections and should only be used for "fast" option
  def.reduceCS = 0;
  
  def.vdist     = max(1,mean(vx_vol0)); % distance between vertices ... at least 1 mm ?
  def.LAB       = cat_get_defaults('extopts.LAB');  
  def.SPM       = 0; 
  def.pbtmethod = 'pbt2x';
  def.WMT       = 0; % WM/CSF width/depth/thickness
  def.sharpenCB = 0; % in development
  def.thick_measure      = 1; % Tfs: Freesurfer method using mean(Tnear1,Tnear2)
  def.thick_limit        = 5; % 5mm upper limit for thickness (same limit as used in Freesurfer)
  def.extract_pial_white = 0; % Estimate pial and white matter surface (in development and very slow!)
  def.new_release        = 0; % developer flag to test new functionality for new release (currently not used)
  def.pbtlas             = 0;
  def.interpV            = 0.5;
  def.add_parahipp       = cat_get_defaults('extopts.add_parahipp');
  def.scale_cortex       = cat_get_defaults('extopts.scale_cortex');
  def.close_parahipp     = cat_get_defaults('extopts.close_parahipp');
  
  opt            = cat_io_updateStruct(def,opt);
  opt.fast       = any(~cellfun('isempty',strfind(opt.surf,'fst'))) + any(~cellfun('isempty',strfind(opt.surf,'sfst')));
  opt.vol        = any(~cellfun('isempty',strfind(opt.surf,'v')));
  opt.interpV    = max(0.1,min([opt.interpV,1.5]));
  opt.interpVold = opt.interpV; 
  opt.surf       = cat_io_strrep(opt.surf,{'sfst','fst','v'},'');

  % check for self-intersections during surface refinement with CAT_SurfDeform
  if opt.extract_pial_white
    force_no_selfintersections = 1;
  else
    force_no_selfintersections = 0;
  end
  
  if opt.fast==2, opt.reduceCS = 40000; end
  if opt.fast
    opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k'); 
  else
    opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  end

  Psurf = struct(); 

  % correction for 'n' prefix for noise corrected and/or interpolated files
  [pp,ff]   = spm_fileparts(V.fname);

  if cat_get_defaults('extopts.subfolders')
    if strcmp(opt.pbtmethod,'pbt3')
      surffolder = sprintf('surf_%s_%0.2f',opt.pbtmethod,opt.interpV);
    elseif strcmp(opt.pbtmethod,'pbt2xf')
      opt.pbtmethod = 'pbt2x';
      surffolder = sprintf('surf_%s_%0.2f',opt.pbtmethod,opt.interpV);
    else
      surffolder = 'surf';
    end
    mrifolder = 'mri';
    pp = spm_str_manip(pp,'h'); % remove 'mri' in pathname that already exists
    if ~exist(fullfile(pp,surffolder),'dir'), mkdir(fullfile(pp,surffolder)); end
  else
    surffolder = '';
    mrifolder = '';
  end

  if ff(1)=='n'
    if (exist(fullfile(pp,[ff(2:end) '.nii']), 'file')) || (exist(fullfile(pp,[ff(2:end) '.img']), 'file'))
      ff = ff(2:end);
    end
  end

  % function to estimate the number of interactions of the surface deformation: d=distance in mm and a=accuracy 
  QMC    = cat_io_colormaps('marks+',17);
  color  = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  rate   = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
 
  
  %% get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
    
  % noise reduction for higher resolutions (>=1 mm full correction, 1.5 mm as lower limit)
  % (added 20160920 ~R1010 due to servere sulcus reconstruction problems with 1.5 Tesla data)
  Yms = Ym + 0; cat_sanlm(Yms,3,1);
  %noise = std(Yms(Yms(:)>0) - Ym(Yms(:)>0)); % more selective filtering?
  %vx_vol = [0.5;0.75;1;1.25;1.5;2]; [vx_vol
  %min(1,max(0,3-2*mean(vx_vol,2))) min(1,max(0,1-mean(vx_vol,2))/2) 0.5*min(1,max(0,1.5-mean(vx_vol,2)))] % filter test 
  mf  = min(1,max(0,3-2*mean(vx_vol,2))); 
  Ym  = mf * Yms  +  (1-mf) * Ym;
  clear Yms;
   
  % filling
  Ymf  = max(Ym,min(1,YMF)); 
  Ymfs = cat_vol_smooth3X(Ymf,1); 
  Ytmp = cat_vol_morph(YMF,'d',3) & Ymfs>2.3/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs; 
  Ymf = Ymf*3;
  
  %% reduction of artifact, blood vessel, and meninges next to the cortex
  % (are often visible as very thin structures that were added to the WM 
  % or removed from the brain)
  if ~opt.SPM
    Ydiv  = cat_vol_div(Ymf,vx_vol); 
    Ycsfd = cat_vbdist(single(Ymf<1.5),Ymf>1,vx_vol);
    Yctd  = cat_vbdist(single(Ymf<0.5),Ymf>0,vx_vol); 
    Ysroi = Ymf>2  &  Yctd<10  & Ycsfd>0 & Ycsfd<2 & ...
            cat_vol_morph(~NS(Ya,opt.LAB.HC) & ~NS(Ya,opt.LAB.HI) & ...
              ~NS(Ya,opt.LAB.PH) & ~NS(Ya,opt.LAB.VT),'erode',4); 
    Ybv   = cat_vol_morph(Ymf+Ydiv./max(1,Ymf)>3.5,'d') & Ymf>2; 
    Ymf(Ybv) = 1.4; 
    Ymfs  = cat_vol_median3(Ymf,Ysroi | Ybv,Ymf>eps & ~Ybv,0.1); % median filter
    %%
    Ymf   = mf * Ymfs  +  (1-mf) * Ymf;

    %% closing of small WMHs and blood vessels
    %vols = [sum(round(Ymf(:))==1 & Ya(:)>0) sum(round(Ymf(:))==2)  sum(round(Ymf(:))==3)] / sum(round(Ymf(:))>0); 
    %volt = min(1,max(0,mean([ (vols(1)-0.20)*5  (1 - max(0,min(0.3,vols(3)-0.2))*10) ]))); 
    %Ywmh = cat_vol_morph(Ymf>max(2.2,2.5 - 0.3*volt),'lc',volt); 
    %Ymf  = max(Ymf,smooth3(Ywmh)*2.9); 
  
    % gaussian filter? ... only in tissue regions
    %Ymfs = cat_vol_smooth3X(max(1,Ymf),0.5*min(1,max(0,1.5-mean(vx_vol)))); 
    %Ymf(Ymf>1) = Ymfs(Ymf>1);
  end
  if ~debug, clear Ysroi Ymfs Yctd Ybv Ymfs; end
  
  %% sharpening of thin structures (gyri and sulci)
  % WARNING: this will change cortical thickness!
  if ~opt.SPM && opt.sharpenCB
    Ydiv = cat_vol_div(Ymf); %Ydivl  = cat_vol_div(Ymf,vx_vol); 
    Ywmd = cat_vbdist(single(Ymf>2.5),Ymf>1,vx_vol);
    if 0
      %% divergence based
      %  this works in principle but gyral crones and sulcal values are
      %  overestimated ... need limit
      Ymsk = (NS(Ya,opt.LAB.CB) & ((Ymf<2.8 & Ymf>2.0          ) | (Ymf<1.9 & Ymf>1.2         )) ) | ... sulci and gyri in the cerebellum 
             (NS(Ya,opt.LAB.CT) & ((Ymf<2.8 & Ymf>2.0 & Ycsfd>3) | (Ymf<1.9 & Ymf>1.2 & Ywmd>3)) ) | ... distant gyri and sulci in the cerebrum
             (NS(Ya,opt.LAB.PH) & ((Ymf<2.8 & Ymf>2.0 & Ycsfd>3) | (Ymf<1.9 & Ymf>1.2 & Ywmd>3)) );
      Ymf  = min(3,max( min(1,Ymf) , Ymf - (abs(Ydivl) .* Ydiv) .* Ymsk));
    end
    
    if 1
      %% biascorrection based
      % WM 
      Ymsk = ((NS(Ya,opt.LAB.CB) | YMF) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01) ) ) | ...                     % sulci and gyri in the cerebellum 
             (NS(Ya,opt.LAB.PH) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01) ) ) | ...                             % hippocampal gyri
             (NS(Ya,opt.LAB.CT) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01 & Ycsfd>cat_stat_nanmean(Ycsfd(Ycsfd(:)>0 & Ycsfd(:)<100)) )*1.0) ); % distant gyri and sulci in the cerebrum
      Yi   = cat_vol_localstat(Ymf,Ymsk,1,3);
      % GM
      Ymsk = (NS(Ya,opt.LAB.CB) & ( Ymf>1.9 & Ymf<2.2 & Ycsfd>0 & Ydiv>-0.05) ) | ...                   % sulci and gyri in the cerebellum 
             (NS(Ya,opt.LAB.PH) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 ) ) | ...                               % hippocampal gyri
             (NS(Ya,opt.LAB.CT) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 & Ywmd>cat_stat_nanmean(Ywmd(Ywmd(:)>0 & Ywmd(:)<100))*0.2 ) );   % distant gyri and sulci in the cerebrum
      Yi   = Yi + cat_vol_localstat(Ymf,Yi==0 & Ymsk,1,1)/2*3;
      Yi   = cat_vol_localstat(Yi,Yi>0,1,3);
      Yi   = cat_vol_localstat(Yi,Yi>0,1,1); 
      if ~debug, clear Ywmd; end
      %% CSF - instable and not required
      %Ymsk = NS(Ya,opt.LAB.VT) & Ymf>=0.5 & Ymf<1.5;                               % sulci and gyri in the cerebellum 
      %Yi  = Yi + cat_vol_localstat(Ymf,Yi==0 & Ymsk,1,3)*3;
      %%
      Ywi = cat_vol_approx(Yi,'nn',1,2,struct('lfO',2)); 
      
      %%
      Ymf = Ymf./Ywi * 3; 
      if ~debug, clear Ywi Yi; end
    end
    if ~debug, clear Ymsk; end
  end
  if ~debug, clear Ydiv Ycsfd; end
  
  Yth1 = zeros(size(Ymf),'single'); 
  if opt.WMT > 1
    Ywd  = zeros(size(Ymf),'single'); 
    Ycd  = zeros(size(Ymf),'single'); 
  end
  
  [D,I] = cat_vbdist(single(Ya>0)); Ya = Ya(I); % for sides
  
  % use sum of EC's and defect sizes for all surfaces, thus set values initially to 0
  EC            = 0;
  defect_size   = 0;
  defect_area   = 0; 
  defect_number = 0; 


  for si=1:numel(opt.surf)
   
    % surface filenames
    Praw       = fullfile(pp,surffolder,sprintf('%s.central.nofix.%s.gii',opt.surf{si},ff));    % raw
    Psphere0   = fullfile(pp,surffolder,sprintf('%s.sphere.nofix.%s.gii',opt.surf{si},ff));     % sphere.nofix
    Pcentral   = fullfile(pp,surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % central
    Pcentralr  = fullfile(pp,surffolder,sprintf('%s.central.resampled.%s.gii',opt.surf{si},ff));          % central
    Ppial      = fullfile(pp,surffolder,sprintf('%s.pial.%s.gii',opt.surf{si},ff));             % pial (GM/CSF)
    Pwhite     = fullfile(pp,surffolder,sprintf('%s.white.%s.gii',opt.surf{si},ff));            % white (WM/GM)
    Pthick     = fullfile(pp,surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % FS thickness / GM depth
    Ppbt       = fullfile(pp,surffolder,sprintf('%s.pbt.%s',opt.surf{si},ff));                  % PBT thickness / GM depth
    Pmask      = fullfile(pp,surffolder,sprintf('%s.mask.%s',opt.surf{si},ff));                 % mask
    Ptemp      = fullfile(pp,surffolder,sprintf('%s.temp.%s',opt.surf{si},ff));                 % temporary file
    Pgwo       = fullfile(pp,surffolder,sprintf('%s.depthWMo.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgw        = fullfile(pp,surffolder,sprintf('%s.depthGWM.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgww       = fullfile(pp,surffolder,sprintf('%s.depthWM.%s',opt.surf{si},ff));              % gyrus witdh of the WM / WM depth
    Pgwwg      = fullfile(pp,surffolder,sprintf('%s.depthWMg.%s',opt.surf{si},ff));             % gyrus witdh of the WM / WM depth
    Psw        = fullfile(pp,surffolder,sprintf('%s.depthCSF.%s',opt.surf{si},ff));             % sulcus width / CSF depth / sulcal span
    Pdefects0  = fullfile(pp,surffolder,sprintf('%s.defects.%s',opt.surf{si},ff));              % defects temporary file
    Pdefects   = fullfile(pp,surffolder,sprintf('%s.defects.%s.gii',opt.surf{si},ff));          % defects
    Psphere    = fullfile(pp,surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    Pspherereg = fullfile(pp,surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    Pfsavg     = fullfile(opt.fsavgDir, sprintf('%s.central.freesurfer.gii',opt.surf{si}));     % fsaverage central
    Pfsavgsph  = fullfile(opt.fsavgDir, sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));      % fsaverage sphere    
    Pfsavgmask = fullfile(opt.fsavgDir, sprintf('%s.mask',opt.surf{si}));                       % fsaverage mask    
    
    surffile = {'Praw','Psphere0','Pcentral','Pthick','Ppbt','Pgw','Pgww','Psw',...
      'Pdefects0','Pdefects','Psphere','Pspherereg','Pfsavg','Pfsavgsph','Pwhite','Ppial'};
    for sfi=1:numel(surffile)
      eval(sprintf('Psurf(si).%s = %s;',surffile{sfi},surffile{sfi})); 
    end
        
    % reduce for object area
    switch opt.surf{si}
      case {'lh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==1); Yside = mod(Ya,2)==1; 
      case {'rh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==0); Yside = mod(Ya,2)==0;  
      case {'lc'},  Ymfs = Ymf .* (Ya>0) .*   NS(Ya,opt.LAB.CB).* (mod(Ya,2)==1); Yside = mod(Ya,2)==1; 
      case {'rc'},  Ymfs = Ymf .* (Ya>0) .*   NS(Ya,opt.LAB.CB).* (mod(Ya,2)==0); Yside = mod(Ya,2)==0; 
    end 
    
    switch opt.surf{si}
      case {'lh','rh'}, opt.interpV = opt.interpVold; 
      case {'lc','rc'}, opt.interpV = opt.interpVold / 2 ; 
    end 
    
    % check for cerebellar hemis
    iscerebellum = strcmp(opt.surf{si},'lc') || strcmp(opt.surf{si},'rc');
    
    % scaling factor for reducing patches and refinement for cerebellar hemis 2..4 according to voxel size
    % or 1 for cerebrum
    scale_cerebellum  = 1 + (iscerebellum*max(1,min(3,1/mean(vx_vol,2))));
    
    % get dilated mask of gyrus parahippocampalis and hippocampus of both sides
    if ~iscerebellum
      mask_parahipp = cat_vol_morph(NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.HC),'d',6);
    end
    
    %% thickness estimation
    if si==1, fprintf('\n'); end
    
    switch opt.fast
      case 2, fprintf('%s - fast with registration:\n',opt.surf{si});
      case 1, fprintf('%s - fast without registration:\n',opt.surf{si});
      case 0, fprintf('%s:\n',opt.surf{si});
    end
    
    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1'))); stimet =stime;
    
    % removing background (smoothing to remove artifacts)
    switch opt.surf{si}
      case {'lh','rh'},  [Ymfs,Ysidei,mask_parahipp,BB] = cat_vol_resize({Ymfs,Yside,mask_parahipp},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
      case {'lc','rc'},  [Ymfs,Ysidei,BB] = cat_vol_resize({Ymfs,Yside},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
    end
    interpBB = BB; interpBB.interpV = opt.interpV; 
     
    imethod         = 'cubic'; % cubic should be better in general - however, linear is better for small thickness (version?)
    [Ymfs,resI]     = cat_vol_resize(max(1,Ymfs),'interp',V,opt.interpV,imethod);                  % interpolate volume
    Ysidei          = cat_vol_resize(single(Ysidei>0.5),'interp',V,opt.interpV,imethod)>0.5;       % interpolate volume (small dilatation)
    
    if ~iscerebellum
      mask_parahipp   = cat_vol_resize(mask_parahipp,'interp',V,opt.interpV)>0.5;          % interpolate volume
    end 
    
    Ymfs = min(3,max(1,Ymfs));

    %% pbt calculation
    if strcmp(opt.pbtmethod,'pbt3')
      [Yth1i,Yppi] = cat_vol_pbt3(Ymfs,struct('method',opt.pbtmethod,'cb',iscerebellum,'resV',opt.interpV,'vmat',V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1])); % avoid underestimated thickness in gyri
    else
      [Yth1i,Yppi] = cat_vol_pbt(Ymfs,struct('method',opt.pbtmethod,'resV',opt.interpV,'vmat',V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'pbtlas',opt.pbtlas)); % avoid underestimated thickness in gyri
    end  
    %%
    Yth1i(Yth1i>10)=0; Yppi(isnan(Yppi))=0;  
    [D,I] = cat_vbdist(Yth1i,Ysidei); Yth1i = Yth1i(I); clear D I Ysidei;   % add further values around the cortex
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI); %clear Yth1i;             % back to original resolution
    Yth1t = cat_vol_resize(Yth1t,'dereduceBrain',BB);                       % adding background
    Yth1  = max(Yth1,Yth1t .* Yside);                                       % save on main image
    clear Yth1t;
    %fprintf('%5.0fs\n',etime(clock,stime)); 
    
    if opt.vol
      S = struct(); Psurf = '';
      fprintf('%5.0fs\n',etime(clock,stime)); 
      continue
    end
    
    %% PBT estimation of the gyrus and sulcus width 
    if opt.WMT > 1 
      %% gyrus width / WM depth
      %  For the WM depth estimation it is better to use the L4 boundary
      %  and correct later for thickness, because the WM is very thin in
      %  gyral regions and will cause bad values. 
      %  On the other side we do not want the whole filled block of the 
      %  Yppi map and so we have to mix both the original WM map and the
      %  Yppi map. 
      %  As far as there is no thickness in pure WM regions there will
      %  be no correction. 
      %
      %    figure, isosurface(smooth3(Yppi),0.5,Yth1i), axis equal off
      stime = cat_io_cmd('  WM depth estimation');
      [Yar,Ymr,BB] = cat_vol_resize({Ya,Ym},'reduceBrain',vx_vol,BB.BB);    % removing background
      Yar   = uint8(cat_vol_resize(Yar,'interp',V,opt.interpV,'nearest'));  % interpolate volume
      Ymr   = cat_vol_resize(Ymr,'interp',V,opt.interpV);                   % interpolate volume
      switch opt.surf{si}
        case {'lh'}, 
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==1);
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==1)); 
        case {'rh'},
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==0);    
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==0)); 
        case {'lc'}, Ymr = Ymr .* (Yar>0) .* NS(Yar,3) .* (mod(Yar,2)==1);
        case {'rc'}, Ymr = Ymr .* (Yar>0) .* NS(Yar,3) .* (mod(Yar,2)==0);
      end 
     % clear Yar; 
      %%
      Yppis = Yppi .* (1-Ynw) + max(0,min(1,Ymr*3-2)) .* Ynw;                         % adding real WM map 
      Ywdt  = cat_vol_eidist(1-Yppis,ones(size(Yppis),'single'));                     % estimate distance map to central/WM surface
      Ywdt  = cat_vol_pbtp(max(2,4-Ymfs),Ywdt,inf(size(Ywdt),'single'))*opt.interpV;
      [D,I] = cat_vbdist(single(Ywdt>0.01),Yppis>0); Ywdt = Ywdt(I); clear D I Yppis; % add further values around the cortex
      Ywdt  = cat_vol_median3(Ywdt,Ywdt>0.01,Ywdt>0.01);                    
      Ywdt = cat_vol_localstat(Ywdt,Ywdt>0.1,1,1);     % smoothing
      Ywdt  = cat_vol_resize(Ywdt,'deinterp',resI);                                   % back to original resolution
      Ywdt  = cat_vol_resize(Ywdt,'dereduceBrain',BB);                                % adding background
      Ywd   = max(Ywd,Ywdt); 
      clear Ywdt;
      
      %% sulcus width / CSF depth
      %  for the CSF depth we cannot use the origal data, because of
      %  sulcal blurring, but we got the PP map at half distance and
      %  correct later for half thickness
      fprintf('%5.0fs\n',etime(clock,stime)); 
      stime = cat_io_cmd('  CSF depth estimation');
      YM    = single(smooth3(cat_vol_morph(Ymr<0.1,'o',4))<0.5); YM(YM==0)=nan;       % smooth CSF/background-skull boundary 
      Yppis = Yppi .* ((Ymr+0.25)>Yppi) + min(1,Ymr*3-1) .* ((Ymr+0.25)<=Yppi);       % we want also CSF within the ventricle (for tests)
      Ycdt  = cat_vol_eidist(Yppis,YM);                                               % distance to the cental/CSF-GM boundary
      Ycdt  = cat_vol_pbtp(max(2,Ymfs),Ycdt,inf(size(Ycdt),'single'))*opt.interpV; Ycdt(isnan(Ycdt))=0;
      [D,I] = cat_vbdist(single(Ycdt>0),Yppis>0 & Yppis<3); Ycdt = Ycdt(I); clear D I Yppis; % add further values around the cortex
      Ycdt  = cat_vol_median3(Ycdt,Ycdt>0.01,Ycdt>0.01);                              % median filtering
      Ycdt = cat_vol_localstat(Ycdt,Ycdt>0.1,1,1);                                    % smoothing
      Ycdt  = cat_vol_resize(Ycdt,'deinterp',resI);                                   % back to original resolution
      Ycdt  = cat_vol_resize(Ycdt,'dereduceBrain',BB);                                % adding background
      Ycd   = max(Ycd,Ycdt); 
      clear Ycdt;
      %fprintf('%5.0fs\n',etime(clock,stime));
      clear Ymr;
    end
    
    if debug, Yppio=Yppi; end
    fprintf('%5.0fs\n',etime(clock,stime));
    
    %% Replace isolated voxels and holes in Ypp by its median value
    
    % indicate isolated holes and replace by median of the neighbors
    Yppi(Yppi<0.35 & ~cat_vol_morph(Yppi<1,'l'))=1;  % close major wholes in the WM 
    Ymsk = Yppi==0 & cat_vol_morph(Yppi>0.9,'d',1); % filter small wholes close to the WM
    Yppi = cat_vol_median3(single(Yppi),Ymsk,~Ymsk); 
    
    %% indicate isolated objects and replace by median of the neighbors
    Yppi(Yppi>0.65 & cat_vol_morph(Yppi==0,'l'))=0;
    Ymsk = Yppi>0.95 & cat_vol_morph(Yppi<0.1,'d',1); 
    Yppi = cat_vol_median3(single(Yppi),Ymsk,~Ymsk);
    if ~debug, clear Ymsk; end
    
    %% Write Ypp for final deformation
    %  Write Yppi file with 1 mm resolution for the final deformation, 
    %  because CAT_DeformSurf achieved better results using that resolution
    Yppt = cat_vol_resize(Yppi,'deinterp',resI);                        % back to original resolution
    Yppt = cat_vol_resize(Yppt,'dereduceBrain',BB);                     % adding of background
    
    % scale Yppt so that backgrounds remains 0 and WM 1, but cortical band is 
    % now in the range of 0.1..0.9
    if opt.extract_pial_white
      indi = find((Yppt>0) & (Yppt<0.99999));
      Yppt(indi) = 0.1 + (0.8*Yppt(indi));
    end
    Vpp  = cat_io_writenii(V,Yppt,'',sprintf('%s.pp',opt.surf{si}),'percentage position map','uint8',[0,1/255],[1 0 0 0]);
    
    Vpp1 = Vpp; 
    Vpp1.fname    = fullfile(pp,mrifolder,sprintf('%s.pp1%s.nii',opt.surf{si},ff));
    vmat2         = spm_imatrix(Vpp1.mat);
    Vpp1.dim(1:3) = round(Vpp1.dim .* abs(vmat2(7:9)*(1 + iscerebellum)));   % use double resolution in case of cerebellum
    vmat2(7:9)    = sign(vmat2(7:9)).*[1 1 1]/(1 + iscerebellum);            % use double resolution in case of cerebellum
    Vpp1.mat      = spm_matrix(vmat2);

    Vpp1 = spm_create_vol(Vpp1); 
    for x3 = 1:Vpp1.dim(3),
      M    = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1]) * inv(Vpp1.mat) * Vpp.mat); %#ok<MINV>
      v    = spm_slice_vol(Vpp,M,Vpp1.dim(1:2),1);       
      Vpp1 = spm_write_plane(Vpp1,v,x3);
    end;
    clear M v x3; 

    %% surface coordinate transformations
    if opt.verb>2
      stime = clock; cat_io_cprintf('g5',sprintf('  Create topology opt. surface '));
    else
      stime = cat_io_cmd('  Create initial surface','g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
    end
    vmatBBV = spm_imatrix(V.mat);

     % surface coordinate transformations that are used in the "saveCS" and "loadCS" functions  
    mati  = spm_imatrix(V.mat); 
    vmat  = V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    vmati = inv([vmat; 0 0 0 1]); vmati(4,:) = [];    

    %% transformation matrix
    matIBB            = spm_imatrix(V.mat * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    matIBB(7:9)       = sign( mati(7:9)) .* repmat( opt.interpV , 1 , 3); 
    Smat.matlabi_mm   = V.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    Smat.matlabIBB_mm = spm_matrix(matIBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % .* mati(7:9)'
    
    
    % smooth mask to have smooth border
    if ~iscerebellum
      mask_parahipp_smoothed = zeros(size(mask_parahipp));
      spm_smooth(double(mask_parahipp),mask_parahipp_smoothed,[8 8 8]);
    end 

    % parameter for isosurface of Yppi
    th_initial = 0.5;

    ind0   = find(Yppi<=0);
    Yppisc = opt.scale_cortex*Yppi;
    
    if ~iscerebellum
      Yppisc  = Yppisc + opt.add_parahipp/opt.scale_cortex*mask_parahipp_smoothed;
    end
    Yppisc(ind0) = 0;
    clear ind0;

    % optionally apply closing inside mask for parahippocampal gyrus to get rid of the holes that lead to large cuts in gyri
    % after topology correction
    if opt.close_parahipp && ~iscerebellum
      tmp = cat_vol_morph(Yppisc,'labclose',1);
      Yppisc(mask_parahipp) = tmp(mask_parahipp);
    end

    if opt.reduceCS>0
      % apply voxel-based topology correction only for smaller defects < 30 voxel
      [tmp,CS] = cat_vol_genus0opt(Yppisc,th_initial,15 * (1-iscerebellum),debug);
      
      % correction for the boundary box used within the surface creation process 
      CS = cat_surf_fun('smat',CS,Smat.matlabIBB_mm); % translate to mm coordinates 
      [Yvxdef,defect_number0] = spm_bwlabel( double(abs(tmp - (Yppisc>0.5))>0) ); clear tmp
    else
      % if no mesh reduction is selected use lower-scaled Yppt with original voxel size
      Yppt = cat_vol_resize(Yppisc,'deinterp',resI);                        % back to original resolution
      Yppt = cat_vol_resize(Yppt,'dereduceBrain',BB);                     % adding of background

      % apply voxel-based topology correction only for smaller defects < 30 voxel
      [tmp,CS] = cat_vol_genus0opt(Yppt,th_initial,15 * (1-iscerebellum),debug);
      CS   = cat_surf_fun('smat',CS,Smat.matlabi_mm); % translate to mm coordinates 
      [Yvxdef,defect_number0] = spm_bwlabel( double(abs(tmp - (Yppt>0.5))>0) ); clear tmp
    end
    EC0            = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
    vdefects       = cat_surf_fun('isocolors',Yvxdef,CS.vertices)>0; clear Yvxdef; 
    defect_size0   = sum(vdefects > 0) / length(vdefects) * 100; % percent
    defect_area0   = sum(vdefects > 0) / length(vdefects) .* ...
      sum(cat_surf_fun('area',CS)) / opt.interpV / 100; % cm2
    if opt.verb>2
      cat_io_cprintf('g5',sprintf('( SC/EC/DN/DS = %0.1f/',opt.scale_cortex));
      cat_io_cprintf( color( rate( abs( EC0 - 2 ) , 0 ,100 * (1+4*iscerebellum) )) ,sprintf('%d/',EC0));
      cat_io_cprintf( color( rate( defect_number0 , 0 ,100 * (1+4*iscerebellum) )) ,sprintf('%d/',defect_number0));
      cat_io_cprintf( color( rate( defect_size0   , 1 , 10 * (1+4*iscerebellum) )) ,sprintf('%0.2f%%%%' ,defect_size0));
      cat_io_cprintf('g5',' )');
      fprintf(repmat(' ',1,max(0,14 - numel(sprintf('%d/%d/%0.2f. )',EC0,defect_number0,defect_size0))))); 
    end
    clear Yppisc;
    
    if opt.verb>2 
      fprintf(txt); 
      fprintf('%s %4.0fs\n',repmat(' ',1,66),etime(clock,stime)); 
    end
      
    % correct the number of vertices depending on the number of major objects
    if opt.reduceCS>0 
      CS = reducepatch(CS,opt.reduceCS * scale_cerebellum); % adaption for cerebellum
      if opt.verb>2     
        stime = cat_io_cmd(sprintf('  Reduce surface to %d faces:',size(CS.faces,1)),'g5','',opt.verb);
      elseif opt.verb>0
        stime = cat_io_cmd(sprintf('  Reduce surface to %d faces:',size(CS.faces,1)),'g5','',opt.verb,stime);
      end
    end
    
    
    %% transform coordinates 
    if opt.fast==1
      %%
      save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),Pcentral);    

      % remove some unconnected meshes
      cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Pcentral,Pcentral); % CAT_SeparatePolygon works here
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      
      % deform initial surface to central surface
      % CAT_DeformSurf "vol" "activity_file?|none" nx ny nz "inputmesh" "outputmesh" "originalposition|none" maxdist n_modls  
      %'CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" ' ... "vol" "activity_file?|none" nx ny nz "inputmesh" "outputmesh"
      %              'none 0 1 -1 .1 ' ...                  "originalposition|none"   maxdist  n_modls  up_to_n_points  model_weight
      %               'avg -0.1 0.1 ' ...                    "model_file...|avg|none"  mincurv  maxcurv 
      %               '.2 .1 5 0 ' ...                       fract_step  max_step  max_search_distance  degrees_continuity  
      %               '"%g" "%g" n ' ...                     min_isovalue  max_isovalue  +/-/n 
      %               '0 0 0 ' ...                           gradient_threshold  angle tolerance  
      %               '150 0.01 0.0 0'], ...                 max_iterations movement_threshold  stop_treshold force_no_selfintersections
                     
      th = 0.5;
      cmds = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" ' ... "vol" "activity_file?|none" nx ny nz "inputmesh" "outputmesh"
                     'none  0  1  -1  .1 ' ...               "originalposition|none"   maxdist  n_modls  up_to_n_points  model_weight
                     'avg  -0.1  0.1 ' ...                   "model_file...|avg|none"  mincurv  maxcurv 
                     '.2  .1  2  0 ' ...                     fract_step  max_step  max_search_distance  degrees_continuity  
                     '"%g"  "%g"  n ' ...                    min_isovalue  max_isovalue  +/-/n 
                     '0  0  0 ' ...                          gradient_threshold  angle  tolerance  
                     '10  0.03  0.0 0'], ...                 max_iterations movement_threshold  stop_threshold force_no_selfintersections
                      Vpp.fname,Pcentral,Pcentral,th,th);
      [ST, RS] = cat_system(cmds); cat_check_system_output(ST,RS,opt.verb-2);
    
      % load surf and project thickness
      CS = gifti(Pcentral);
      % ignore this warning writing gifti with int32 (eg. cat_surf_createCS:580 > gifti/subsref:45)
      warning off MATLAB:subscripting:noSubscriptsSpecified
      facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS,Smat.matlabIBB_mm); 
      cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
      
      % map WM and CSF width data (corrected by thickness)
      if opt.WMT > 1
        %%
        facevertexcdata2  = cat_surf_fun('isocolors',Ywd,CS,Smat.matlabi_mm);
        facevertexcdata2c = max(eps,facevertexcdata2 - facevertexcdata/2);
        cat_io_FreeSurfer('write_surf_data',Pgwo,facevertexcdata2c); % gyrus width WM only
        facevertexcdata2c = correctWMdepth(CS,facevertexcdata2c,100,0.2);
        cat_io_FreeSurfer('write_surf_data',Pgww,facevertexcdata2c); % gyrus width WM only
        facevertexcdata3c = facevertexcdata2c + facevertexcdata; % );
        cat_io_FreeSurfer('write_surf_data',Pgw,facevertexcdata3c); % gyrus width (WM and GM)
        facevertexcdata4 = estimateWMdepthgradient(CS,facevertexcdata2c);
        cat_io_FreeSurfer('write_surf_data',Pgwwg,facevertexcdata4); % gyrus width WM only > gradient
        % smooth resampled values
        try
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Pgwwg,3,Pgwwg);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
        end
        %%
        %clear facevertexcdata2 facevertexcdata2c facevertexcdata3c facevertexcdata4; 
        % just a test ... problem with other species ...
        %norm = sum(Ymf(:)>0.5) / prod(vx_vol) / 1000 / 1400;
        %norm = mean([2 1 1].*diff([min(CS.vertices);max(CS.vertices)])); 
        %norm = mean([2 1 1].*std(CS.vertices)); % maybe the hull surface is better...

        facevertexcdata3 = cat_surf_fun('isocolors',Ycd,CS,Smat.matlabi_mm);  
        facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
        cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3);
      end
      
      % distance between linked surfaces 
      %Tlink = @(S1,S2) sum( [ sum( abs( S1.vertices(:,1:2) - S2.vertices(:,1:2) ).^2 , 2 ).^0.5  abs(S1.vertices(:,3) - S2.vertices(:,3)) ].^2 , 2 ).^0.5;
              
      % save datastructure
      S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata);
      if opt.WMT > 1
        setfield(S.(opt.surf{si}),'th2',nan(size(facevertexcdata)));
        setfield(S.(opt.surf{si}),'th3',nan(size(facevertexcdata)));
      end
      
      if ~debug
        delete(Vpp.fname);
        delete(Vpp1.fname);
      end
      fprintf('%5.0fs\n',etime(clock,stime)); 
      
      % estimate Euler characteristics: EC = #vertices + #faces - #edges
      EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
      EC  = EC + abs(EC0);
      
      % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
      if opt.thick_measure == 1
        cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      
        % apply upper thickness limit
        facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick);  
        facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
        cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
      else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
        copyfile(Ppbt,Pthick);
      end
      
      
      %% intensity based evaluation
      CS = gifti(Pcentral);
      % ignore this warning writing gifti with int32 (eg. cat_surf_createCS:580 > gifti/subsref:45)
      warning off MATLAB:subscripting:noSubscriptsSpecified
      try
        facevertexcdata1 = cat_io_FreeSurfer('read_surf_data',Ppbt);
      catch
        facevertexcdata1 = cat_surf_fun('isocolors',Yth1i,CS,Smat.matlabIBB_mm); 
      end
      cat_surf_fun('saveico',CS,facevertexcdata1,Pcentral,'fast',Ymfs,Smat.matlabIBB_mm);
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',CS,facevertexcdata1,Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);

      %%
      clear CS
      continue
    end
    save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),Praw,'Base64Binary');
    
    if opt.reduceCS>0 
      % after reducepatch many triangles have very large area which causes isses for resampling
      % RefineMesh adds triangles in those areas
      if opt.fast
        cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 0',Praw,Praw,4 * opt.vdist / scale_cerebellum); % adaption for cerebellum
      else  
        cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 0',Praw,Praw,2 * opt.vdist / scale_cerebellum); % adaption for cerebellum
      end
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    
      % remove some unconnected meshes
      cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw); % CAT_SeparatePolygon works here
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    end

    % spherical surface mapping 1 of the uncorrected surface for topology correction
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Praw,Psphere0);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

    % estimate size of topology defects 
    cmd = sprintf('CAT_MarkDefects "%s" "%s" "%s"',Praw,Psphere0,Pdefects0); 
    [ST, RS] = cat_system(cmd);
    sdefects       = cat_io_FreeSurfer('read_surf_data',Pdefects0); delete(Pdefects0);  
    defect_number0 = defect_number0 + ceil( max(sdefects )); 
    defect_size0   = defect_size0   + sum(sdefects > 0) / length(sdefects) * 100; % percent
    defect_area0   = defect_area0   + sum(sdefects > 0) / length(sdefects) .* ...
      sum(cat_surf_fun('area',CS)) / opt.interpV / 100; clear defects; % cm2
    % estimate Euler characteristics: EC = #vertices + #faces - #edges
    EC0            = (EC0-2) + ( size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1) - 2) + 2;
    EC             = EC + abs(EC0 - 2) + 2; % -2 is the correction for the sphere
    defect_size    = defect_size   + defect_size0;
    defect_area    = defect_area   + defect_area0;
    defect_number  = defect_number + defect_number0;
    
    %% topology correction and surface refinement 
    stime = cat_io_cmd('  Topology correction and surface refinement:','g5','',opt.verb,stime); 
    if opt.verb>2, fprintf('\n'); end
    if opt.fast
      cmd = sprintf('CAT_FixTopology -lim 128 -bw 512 -n 40960 -refine_length %g "%s" "%s" "%s"',4 * opt.vdist / scale_cerebellum,Praw,Psphere0,Pcentral);
      try
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      catch
        cmd = sprintf('CAT_FixTopology -lim 128 -bw 512 -n 81920 -refine_length %g "%s" "%s" "%s"',2 * opt.vdist / scale_cerebellum,Praw,Psphere0,Pcentral);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      end  
    else
      cmd = sprintf('CAT_FixTopology -lim 128 -bw 512 -n 81920 -refine_length %g "%s" "%s" "%s"',2  * opt.vdist / scale_cerebellum,Praw,Psphere0,Pcentral);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    end
    
    % read final surface and map thickness data
    CS = gifti(Pcentral);
    % ignore this warning writing gifti with int32 (eg. cat_surf_createCS:580 > gifti/subsref:45)
    warning off MATLAB:subscripting:noSubscriptsSpecified
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS,Smat.matlabIBB_mm); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
  
    % final correction of central surface in highly folded areas with high mean curvature with weight of 0.7
    stime = cat_io_cmd('  Correction of central surface in highly folded areas','g5','',opt.verb,stime);
    cmd = sprintf(['CAT_Central2Pial -equivolume -weight 0.7 "%s" "%s" "%s" 0.2'], ...
                       Pcentral,Ppbt,Pcentral);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
  
    % we need some refinement because some vertices are too large to be deformed with high accuracy
    if opt.fast
      cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 0',Pcentral,Pcentral,4 * opt.vdist / scale_cerebellum); % adaption for cerebellum
    else
      cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 1',Pcentral,Pcentral,2 * opt.vdist / scale_cerebellum); % adaption for cerebellum
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    
    % surface refinement by surface deformation based on the PP map
    stime = cat_io_cmd('  Refine central surface','g5','',opt.verb,stime);
    th = 0.5;
    if opt.fast
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'avg -0.1 0.1 .2 .1 5 0 "%g" "%g" n 0 0 0 50 0.02 0.0 %d'], ...
                     Vpp.fname,Pcentral,Pcentral,th,th,0);
    else
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'avg -0.1 0.1 .2 .1 5 0 "%g" "%g" n 0 0 0 150 0.01 0.0 %d'], ...
                     Vpp.fname,Pcentral,Pcentral,th,th,force_no_selfintersections);
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

    % need some more refinement because some vertices are distorted after CAT_DeformSurf
    if opt.fast
      cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 0',Pcentral,Pcentral,4 * opt.vdist / scale_cerebellum); % adaption for cerebellum
    else
      cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 1',Pcentral,Pcentral,1.6 * opt.vdist / scale_cerebellum); % adaption for cerebellum
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

    if opt.fast
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'flat -0.15 0.15 .5 .1 5 0 "%g" "%g" n 0 0 0 25 0.01 0.0 %d'], ...
                     Vpp.fname,Pcentral,Pcentral,th,th,0);
    else
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.15 0.15 .1 .1 5 0 "%g" "%g" n 0 0 0 150 0.01 0.0 %d'], ...
                     Vpp.fname,Pcentral,Pcentral,th,th,force_no_selfintersections);
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

    % read final surface and map thickness data
    CS = gifti(Pcentral);
    % ignore this warning writing gifti with int32 (eg. cat_surf_createCS:580 > gifti/subsref:45)
    warning off MATLAB:subscripting:noSubscriptsSpecified
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS,Smat.matlabIBB_mm); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);

    % final correction of central surface in highly folded areas with high mean curvature
    if ~opt.fast
      stime = cat_io_cmd('  Correction of central surface in highly folded areas 2','g5','',opt.verb,stime);
      cmd = sprintf(['CAT_BlurSurfHK "%s" "%s" 1'], Pcentral,Pcentral);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      cmd = sprintf(['CAT_Central2Pial -equivolume -weight 0.4 "%s" "%s" "%s" 0'], ...
                       Pcentral,Ppbt,Pcentral);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    end

    % final correction of cortical thickness using pial and WM surface
    if opt.extract_pial_white && ~opt.fast
  
      % estimation of pial surface
      th2 = 0.1; % GM/CSF border
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 300 0.001 0.0 %d'], ...
                     Vpp.fname,Pcentral,Ppial,th2,th2,force_no_selfintersections);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

      % if deformation stopped earlier call CAT_Central2Pial to get closer to the pial surface
      if ~isempty(strfind(RS,'Stopped after'))
        stime = cat_io_cmd('  Estimation of pial surface','g5','',opt.verb,stime);
        cmd = sprintf(['CAT_Central2Pial -check_intersect "%s" "%s" "%s" 0.3'], ...
                       Pcentral,Ppbt,Ppial);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
        
        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 300 0.001 0.0 %d'], ...
                     Vpp.fname,Ppial,Ppial,th2,th2,force_no_selfintersections);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      end

      stime = cat_io_cmd('  Correction of pial surface in highly folded areas','g5','',opt.verb,stime);
      cmd = sprintf(['CAT_Central2Pial -equivolume -weight 0.5 "%s" "%s" "%s" 0'], ...
                       Ppial,Ppbt,Ppial);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

      % estimation of white matter surface
      th2 = 0.9; % GM/WM border
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 300 0.001 0.0 %d'], ...
                     Vpp.fname,Pcentral,Pwhite,th2,th2,force_no_selfintersections);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

      % if deformation stopped earlier call CAT_Central2Pial to get closer to the white matter surface
      if ~isempty(strfind(RS,'Stopped after'))
        stime = cat_io_cmd('  Estimation of white matter surface','g5','',opt.verb,stime);
        cmd = sprintf(['CAT_Central2Pial -check_intersect "%s" "%s" "%s" -0.3'], ...
                       Pcentral,Ppbt,Pwhite);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
        
        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 300 0.001 0.0 %d'], ...
                     Vpp.fname,Pwhite,Pwhite,th2,th2,force_no_selfintersections);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      end

      stime = cat_io_cmd('  Correction of white matter surface in highly folded areas','g5','',opt.verb,stime);
      cmd = sprintf(['CAT_Central2Pial -equivolume -weight 0.5 "%s" "%s" "%s" 0'], ...
                       Pwhite,Ppbt,Pwhite);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

      % update central surface as average between white and pial surface
      cmd = sprintf(['CAT_AverageSurfaces -avg "%s" "%s" "%s"'], ...
                       Pcentral,Pwhite,Ppial);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      
      updateThickness = 1;
      if updateThickness % more testing necessary to also correct thickness
        % correction of cortical thickness
        cmd = sprintf('CAT_Hausdorff  "%s" "%s" "%s"',Pwhite,Ppial,Ppbt);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      end
      
      saveiconame = 'extract_pial_white';
    else
      updateThickness = 0;
      saveiconame = 'default';
    end
    
    
    
    %% just a shortcut for manual tests 
    writedebug = 0; %cat_get_defaults('extopts.expertgui')==2;
    % intensity based evaluation
    CS = gifti(Pcentral);
    % ignore this warning writing gifti with int32 (eg. cat_surf_createCS:580 > gifti/subsref:45)
    warning off MATLAB:subscripting:noSubscriptsSpecified
    if updateThickness
      facevertexcdata1 = cat_io_FreeSurfer('read_surf_data',Ppbt);
    else
      facevertexcdata1 = cat_surf_fun('isocolors',Yth1i,CS,Smat.matlabIBB_mm); 
    end
    fprintf('%5.0fs',etime(clock,stime)); stime = []; 
    if writedebug
      cat_surf_fun('saveico',CS,facevertexcdata1,Pcentral,saveiconame,Ymfs);
    else
      fprintf('\n');
    end
    res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',CS,facevertexcdata1,Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2,cat_get_defaults('extopts.expertgui')>1);

    
    

    %% spherical surface mapping 2 of corrected surface
    stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
    if opt.fast
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Pcentral,Psphere);
    else
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 10',Pcentral,Psphere);
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    
    % spherical registration to fsaverage template
    stime = cat_io_cmd('  Spherical registration','g5','',opt.verb,stime);
    if opt.fast
      cmd = sprintf(['CAT_WarpSurf -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s" ' ...
        '-size 256 128 -loop 1 -steps 1 -runs 1 -v -fwhm 10 -fwhm-surf 20 -lmreg 0.01'],...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    else
      cmd = sprintf('CAT_WarpSurf -steps 2 -avg -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"', ...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    
    
    % set thickness values to zero for masked area (use inverse transformation to map mask)
    % does not work properly for all data...
    if 0
      stime = cat_io_cmd('  Correct thickness','g5','',opt.verb,stime);
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"', ...
        Pfsavg,Pfsavgsph,Pspherereg,Ptemp,Pfsavgmask,Pmask);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);      
      resampled_mask = cat_io_FreeSurfer('read_surf_data',Pmask);
      
      % set thickness to 0 for masked area and write thickness data
      facevertexcdata(resampled_mask < 0.5) = 0;
      cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);  
      delete(Pmask)
      delete(Ptemp);
    end
    
    % map WM and CSF width data (corrected by thickness)
    if opt.WMT > 1
      %%
      facevertexcdata2  = cat_surf_fun('isocolors',Ywd,CS,Smat.matlabi_mm);
      facevertexcdata2c = max(eps,facevertexcdata2 - facevertexcdata/2);
      cat_io_FreeSurfer('write_surf_data',Pgwo,facevertexcdata2c); % gyrus width WM only
      facevertexcdata2c = correctWMdepth(CS,facevertexcdata2c,100,0.2);
      cat_io_FreeSurfer('write_surf_data',Pgww,facevertexcdata2c); % gyrus width WM only
      facevertexcdata3c = facevertexcdata2c + facevertexcdata; % );
      cat_io_FreeSurfer('write_surf_data',Pgw,facevertexcdata3c); % gyrus width (WM and GM)
      facevertexcdata4 = estimateWMdepthgradient(CS,facevertexcdata2c);
      cat_io_FreeSurfer('write_surf_data',Pgwwg,facevertexcdata4); % gyrus width WM only > gradient
      % smooth resampled values
      try
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Pgwwg,3,Pgwwg);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      end
      %%
      %clear facevertexcdata2 facevertexcdata2c facevertexcdata3c facevertexcdata4; 
      % just a test ... problem with other species ...
      %norm = sum(Ymf(:)>0.5) / prod(vx_vol) / 1000 / 1400;
      %norm = mean([2 1 1].*diff([min(CS.vertices);max(CS.vertices)])); 
      %norm = mean([2 1 1].*std(CS.vertices)); % maybe the hull surface is better...
 
      facevertexcdata3 = cat_surf_fun('isocolors',Ycd,CS,Smat.matlabi_mm); 
      facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
      cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3);
    end
    fprintf('%5.0fs\n',etime(clock,stime)); 
    
    if opt.verb>1
      cat_io_cprintf( 'g5', sprintf('    Euler char. / def. number / def. size: '));
      cat_io_cprintf( color( rate(  EC0 - 2        , 0 , 2 * 50 * (1+9*iscerebellum)) ) , sprintf('%0.0f / '   , EC0 ) );
      cat_io_cprintf( color( rate(  defect_number0 , 0 , 2 * 50 * (1+9*iscerebellum)) ) , sprintf('%0.0f / '   , defect_number0 ) );
      cat_io_cprintf( color( rate(  defect_size0   , 0 , 2 * 5  * (1+9*iscerebellum)) ) , sprintf('%0.2f%%%% ' , defect_size0 ) );
      fprintf('\n');
    end
    %fprintf('  Surface Euler number:                  %d\n',EC0);
    %fprintf('  Number of defects:                     %0.0f\n',defect_number);
    %fprintf('  Overall size of topology defects:      %0.2f%% (~%0.2f cm%s)',...
    %  defect_size0,defect_area0,char(178));
    
    if 0 %writedebug 
    % This part is not highly relevant for the individual surface reconstruction 
    % but it can help to test and optimize the spatial registration. 
    
      % filenames for resmapling
      Presamp   = fullfile(pp,surffolder,sprintf('%s.tmp.resampled.%s'    ,opt.surf{si},ff));  
      Ppbtr     = fullfile(pp,surffolder,sprintf('%s.pbt.resampled.%s'    ,opt.surf{si},ff));  
      Ppbtr_gii = [Ppbtr '.gii'];
      
      % resample values using warped sphere 
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavgsph,Presamp,Ppbt,Ppbtr);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      
      if 0 
        % resample surface using warped sphere with better surface quality (using Spherical harmonics)
        % ###
        % deactivated because the resampling of the surface alone leads to displacements of the textures (RD20190927)!
        % ###
        cmd = sprintf('CAT_ResampleSphericalSurfSPH -n 327680 "%s" "%s" "%s"',Pcentral,Pspherereg,Presamp);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);

        % resample surface according to freesurfer sphere
        cmd = sprintf('CAT_ResampleSurf "%s" NULL "%s" "%s"',Presamp,Pfsavgsph,Presamp);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2); 
      end
      
      % add values to resampled surf and save as gifti
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Ppbtr,Ppbtr_gii); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2); 
      if exist(Ppbtr,'file'), delete(Ppbtr); end

      % remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
      [pp2,ff2,ex2] = spm_fileparts(Ppbtr_gii);
      g = gifti(Ppbtr_gii);
      g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
      save(g, Ppbtr_gii, 'Base64Binary');
      
      % intensity based evaluation
      CS1 = gifti(Ppbtr_gii);
      CSr = struct('vertices',CS1.vertices,'faces',CS1.faces,'cdata',CS1.cdata,'vmat',vmat,'mati',mati); 
      CSr.vertices = (vmati*[CSr.vertices' ; ones(1,size(CSr.vertices,1))])';
      if mati(7)<0, CSr.faces = [CSr.faces(:,1) CSr.faces(:,3) CSr.faces(:,2)]; end
      warning off MATLAB:subscripting:noSubscriptsSpecified
      cat_surf_fun('saveico',CSr,CSr.cdata,Pcentralr,[saveiconame '_resampled']);
      res.(opt.surf{si}).createCS_resampled = cat_surf_fun('evalCS',CSr,CSr.cdata,Ymfs,Yppi,Pcentralr);
      clear CSr CS1
    end
    if ~isfield( res.(opt.surf{si}),'createCS_final')
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',CS1,cat_io_FreeSurfer('read_surf_data',Ppbt),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2,cat_get_defaults('extopts.expertgui')>1);
    else 
      fprintf('\n'); 
    end
    %clear Yppi; 
    
    % visualize a side
    % csp=patch(CS); view(3), camlight, lighting phong, axis equal off; set(csp,'facecolor','interp','edgecolor','none')
    
    % create output structure
    warning off MATLAB:subscripting:noSubscriptsSpecified
    S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata);
    if opt.WMT > 1
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th2',facevertexcdata2);
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th3',facevertexcdata3);
    end
    clear Yth1i

    % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
    if opt.thick_measure == 1
      % not ready yet
      if 0
%      if opt.extract_pial_white && ~opt.fast % use white and pial surfaces
        cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',Pwhite,Ppial,Pthick);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      else % use central surface and thickness
        cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
      end
      
      % apply upper thickness limit
			facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick);  
			facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
			cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
    else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
      copyfile(Ppbt,Pthick);
    end
    
    
    %% average final values
    FNres = fieldnames( res.(opt.surf{si}).createCS_final );
    for fnr = 1:numel(FNres)
      if ~isfield(res,'final') || ~isfield(res.final,FNres{fnr})
        res.final.(FNres{fnr}) = res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      else
        res.final.(FNres{fnr}) = res.final.(FNres{fnr}) + res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      end
    end
    if isfield(res.(opt.surf{si}),'createCS_resampled') 
      FNres = fieldnames( res.(opt.surf{si}).createCS_resampled );
      for fnr = 1:numel(FNres)
        if ~isfield(res,'createCS_resampled') || ~isfield(res.createCS_resampled,FNres{fnr}) 
          res.resampled.(FNres{fnr}) = res.(opt.surf{si}).createCS_resampled.(FNres{fnr}) / numel(opt.surf);
        else
          res.resampled.(FNres{fnr}) = res.resampled.(FNres{fnr}) + res.(opt.surf{si}).createCS_resampled.(FNres{fnr}) / numel(opt.surf);
        end
      end
    end
    
    
    % we have to delete the original faces, because they have a different number of vertices after
    % CAT_FixTopology!
    delete(Praw);  
    if opt.verb > 2
      delete(Pdefects0);  
    end
    delete(Psphere0);
    delete(Vpp.fname);
    delete(Vpp1.fname);
    clear CS
    
    % create white and central surfaces
    if cat_get_defaults('extopts.expertgui') == 2
      cat_surf_fun('white',Pcentral);
      cat_surf_fun('pial',Pcentral);
    end
  end  
  
  % calculate mean EC and defect size for all surfaces
  mnth = []; sdth = []; mnRMSE_Ypp = []; mnRMSE_Ym = []; sdRMSE_Ym = []; sdRMSE_Ypp = []; 
  SIw = []; SIp = []; SIwa = []; SIpa = []; 
  for si=1:numel(opt.surf)
    if any(strcmp(opt.surf{si},{'lh','rh'}))
      mnth        = [ mnth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(1) ]; 
      sdth        = [ sdth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(2) ]; 
      mnRMSE_Ym   = [ mnRMSE_Ym   mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4 ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ]) ];
      sdRMSE_Ym   = [ sdRMSE_Ym  std([...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4 ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ]) ];
      mnRMSE_Ypp  = [ mnRMSE_Ypp  mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_central ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ]) ];
      sdRMSE_Ypp  = [ sdRMSE_Ypp  std([...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_central ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ]) ];
      if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
        SIw     = [ SIw  res.(opt.surf{si}).createCS_final.white_self_interections ]; 
        SIp     = [ SIp  res.(opt.surf{si}).createCS_final.pial_self_interections  ]; 
        SIwa    = [ SIwa res.(opt.surf{si}).createCS_final.white_self_interection_area ]; 
        SIpa    = [ SIpa res.(opt.surf{si}).createCS_final.pial_self_interection_area  ]; 
      end
    end
  end

  EC            = EC / numel(opt.surf);
  defect_area   = defect_area / numel(opt.surf);
  defect_size   = defect_size / numel(opt.surf);
  defect_number = defect_number / numel(opt.surf);
  
  % final res structure
  res.Smat        = Smat; 
  res.EC          = EC; 
  res.defect_size = defect_size;
  res.defect_area = defect_area;
  res.defects     = defect_number;
  res.RMSE_Ym     = mean(mnRMSE_Ym);
  res.RMSE_Ypp    = mean(mnRMSE_Ypp);
  if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
    res.self_intersections      = mean([SIw,SIp]);
    res.self_intersections_area = mean([SIwa,SIpa]);
  end
    
  
  if opt.verb && ~opt.vol  
    % display some evaluation 
    % - For normal use we limited the surface measures.  
    % - Surface intensity would be interesting as cortical measure similar to thickness (also age dependent).
    %   Especially the outer surface will describe the sulcal blurring in childeren. 
    %   But the mixing of surface quality and anatomical features is problematic. 
    % - The position value describes how good the tranformation of the PBT map into a surface worked. 
    %   Also the position values depend on age. Children have worste pial values due to sulcal blurring but
    %   the white surface is may effected by aging, e.g. by WMHs.
    % - However, for both intensity and position some (average) maps would be also interesting. 
    %   Especially, some Kappa similar measure that describes the differences to the Ym or Ypp would be nice.
    % - What does the Euler chararteristic say?  Wouldn't the defect number more useful for users? 
    if any(~cellfun('isempty',strfind(opt.surf,'cb'))), cbtxt = 'cerebral '; else cbtxt = ''; end
    fprintf('Final %ssurface processing results: \n', cbtxt);
      
    if cat_get_defaults('extopts.expertgui')
    % color output currently only for expert ...
      fprintf('  Average thickness:                     ');
      cat_io_cprintf( color( rate( abs( mean(mnth) - 2.5 ) , 0 , 2.0 )) , sprintf('%0.4f'  , mean(mnth) ) );  fprintf(' %s ',native2unicode(177, 'latin1'));
      cat_io_cprintf( color( rate( abs( mean(sdth) - 0.5 ) , 0 , 1.0 )) , sprintf('%0.4f mm\n', mean(sdth) ) );
  
      fprintf('  Surface intensity / position RMSE:     ');
      cat_io_cprintf( color( rate( mean(mnRMSE_Ym)  , 0.05 , 0.3 ) ) , sprintf('%0.4f / ', mean(mnRMSE_Ym) ) );
      cat_io_cprintf( color( rate( mean(mnRMSE_Ypp) , 0.05 , 0.3 ) ) , sprintf('%0.4f\n', mean(mnRMSE_Ypp) ) );
    
      if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
        fprintf('  Pial/white self-intersections:         ');
        cat_io_cprintf( color( rate(  mean([SIw,SIp]) , 0 , 20 ) ) , sprintf('%0.2f%%%% (%0.2f mm%s)\n'  , mean([SIw,SIp]) , mean([SIwa,SIpa]) , char(178) ) );
      end
      
      fprintf('  Euler char. / def. number / def. size: ');
      cat_io_cprintf( color( rate(  EC - 2        , 0 , 100 * (1+9*iscerebellum)) ) , sprintf('%0.1f / '   , EC ) );
      cat_io_cprintf( color( rate(  defect_number , 0 , 100 * (1+9*iscerebellum)) ) , sprintf('%0.1f / '   , defect_number ) );
      cat_io_cprintf( color( rate(  defect_size   , 0 ,  10 * (1+9*iscerebellum)) ) , sprintf('%0.2f%%%% ' , defect_size ) );
      fprintf('\n');
    else
      fprintf('  Average thickness:                     %0.4f %s %0.4f mm\n' , mean(mnth), native2unicode(177, 'latin1'), mean(sdth));
      %fprintf('  Surface intensity / position RMSE:     %0.4f & %0.4f\n'    , mean(sdRMSE_Ym) ,mean(sdRMSE_Ypp) );
      fprintf('  Euler characteristic / defect size:    %0d / %0.2f%%%% \n'  , EC, defect_size);
    end
    
    for si=1:numel(Psurf)
      fprintf('  Display thickness: %s\n',spm_file(Psurf(si).Pthick,'link','cat_surf_display(''%s'')'));
    end
  end
end

function varargout = cat_vol_genus0opt(Yo,th,limit,debug)
% cat_vol_genus0opt: Voxel-based topology optimization and surface creation 
%   The correction of large defects is often not optimal and this function
%   uses only small corrections. 
% 
%    [Yc,S] = cat_vol_genus0vol(Yo[,limit,debug])
%  
%    Yc    .. corrected volume 
%    Yo    .. original volume 
%    S     .. surface 
%    th    .. threshold for creating surface
%    limit .. maximum number of voxels to correct a defect (default = 30)
%    debug .. print details.  
%

  if nargin < 2, th = 0.5; end
  if nargin < 3, limit = 30; end
  if nargin < 4, debug = 0; end
  
  Yc = Yo; nooptimization = limit==0;  %#ok<NASGU>
  if limit==0
    % use all corrections
    if nargout>1
      txt = evalc(sprintf('[Yc,S.faces,S.vertices] = cat_vol_genus0(Yo,th,nooptimization);'));
    else
      txt = evalc(sprintf('Yc = cat_vol_genus0(Yo,th,nooptimization);'));
    end
    
    if debug, fprintf(txt); end
  else
    % use only some corrections
    txt = evalc(sprintf('Yc = cat_vol_genus0(Yo,th,nooptimization);'));
    
    % remove larger corrections
    Yvxcorr = abs(Yc - (Yo > th))>0; 
    Yvxdef  = spm_bwlabel( double( Yvxcorr ) ); clear Yppiscrc; 
    Yvxdef  = cat_vol_morph(Yvxdef,'l',[inf limit]) > 0; % large corrections that we remove 
    
    if debug
      fprintf(txt); 
      fprintf('  Number of voxels of genus-topocorr: %d\n  Finally used corrections:  %0.2f%%\n', ...
        sum(Yvxcorr(:)) , 100 * sum(Yvxcorr(:) & ~Yvxdef(:)) / sum(Yvxcorr(:)) );
    end
    
    Yc = Yc & ~Yvxdef; 
  
    % final surface creation without correction
    if nargout>1
      evalc(sprintf('[Yt,S.faces,S.vertices] = cat_vol_genus0( single(Yc) ,th,1);')); 
    end
  
  end
  
  varargout{1} = Yc; 
  if nargout>1, varargout{2} = S; end
end

%=======================================================================
function [cdata,i] = correctWMdepth(CS,cdata,iter,lengthfactor)
% ______________________________________________________________________
% Correct deep WM depth values that does not fit to the local thickness 
% of the local gyri.
% 
% lengthfactor should be between 0.2 and 0.4
% ______________________________________________________________________

  if ~exist('lengthfactor','var'), lengthfactor = 1/3; end
  if ~exist('iter','var'), iter = 100; end

  %%
  SV  = CS.vertices;                                                          % Surface Vertices 
  SE  = unique([CS.faces(:,1:2);CS.faces(:,2:3);CS.faces(:,3:-2:1)],'rows');  % Surface Edges
  SEv = single(diff(cat(3,SV(SE(:,1),:),SV(SE(:,2),:)),1,3));                 % Surface Edge Vector
  SEL = sum(SEv.^2,2).^0.5;                                                   % Surface Edge Length  
  clear SEv

  
  %%
  i=0; cdatac = cdata+1; pc = 1; oc = 0; 
  while i<iter && pc~=oc; 
  %%
    pc = sum( abs(cdata - cdatac)>0.05 ); 
    i=i+1; cdatac = cdata;
    
    M  = (cdatac(SE(:,1)) - SEL(SE(:,1))*lengthfactor ) > cdatac(SE(:,2)); 
    cdata(SE(M,1)) = cdatac(SE(M,2)) + SEL(SE(M,1))*lengthfactor; 
    M  = (cdata(SE(:,2)) - SEL(SE(:,2))*lengthfactor ) > cdatac(SE(:,1));
    cdata(SE(M,2)) = cdatac(SE(M,1)) + SEL(SE(M,1))*lengthfactor; 
    oc = sum( abs(cdata - cdatac)>0.05 );
    
    %fprintf('%d - %8.2f - %d\n',i,sum( abs(cdata - cdatac)>0.05 ),pc~=oc)
    
  end
  
end
%=======================================================================
function V = isocolors2(R,V,opt)
% ______________________________________________________________________
% calculates an interpolated value of a vertex in R  
% We have to calculate everything with double, thus larger images will 
% cause memory issues.
% ______________________________________________________________________
  
  if isempty(V), return; end
  if ndims(R)~=3,  error('MATLAB:isocolor2:dimsR','Only 2 or 3 dimensional input of R.'); end
  if ~exist('opt','var'), opt=struct(); end
  
  def.interp = 'linear';
  opt = cat_io_checkinopt(opt,def);
  
  if  isa(R,'double'), R = single(R); end
  if ~isa(V,'double'), V = double(V); VD=0; else VD=1; end
  
  nV   = size(V,1);
  ndim = size(V,2);
  
  switch opt.interp
    case 'nearest'
      V = max(1,min(round(V),repmat(size(R),nV,1))); 
      V = R(sub2ind(size(R),V(:,2),V(:,1),V(:,3)));
    case 'linear'
      nb  = repmat(shiftdim(double([0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]'),-1),nV,1);  
      enb = repmat(shiftdim((ones(8,1,'double')*[size(R,2),size(R,1),size(R,3)])',-1),nV,1);  

      % calculate the weight of a neigbor (volume of the other corner) and
      w8b = reshape(repmat(V,1,2^ndim),[nV,ndim,2^ndim]); clear V;
      % if the streamline is near the boundary of the image you could be out of range if you add 1 
      n8b = min(floor(w8b) + nb,enb); clear enb
      n8b = max(n8b,1);
      w8b = flipdim(prod(abs(n8b - w8b),2),3);        

      % multiply this with the intensity value of R
      V = sum(R(sub2ind(size(R),n8b(:,2,:),n8b(:,1,:),n8b(:,3,:))) .* w8b,3);
  end  
  if ~VD, V = single(V); end
end
    
%=======================================================================
function cdata = estimateWMdepthgradient(CS,cdata)
% ______________________________________________________________________
% Estimates the maximum local gradient of a surface. 
% Major use is the WM depth that grows with increasing sulcal depth. 
% It measures the amount of WM behind the cortex, but more relevant is
% the amount of WM fibers that this region will add to the WM depth. 
% The width of the street next to a house gives not the connectivity of
% this house, but the width of the entrance does!
% This measure can be improved by furhter information of sulcal depth.
% ______________________________________________________________________

  %%
  SV  = CS.vertices;                                                          % Surface Vertices 
  SE  = unique([CS.faces(:,1:2);CS.faces(:,2:3);CS.faces(:,3:-2:1)],'rows');  % Surface Edges
  SEv = single(diff(cat(3,SV(SE(:,1),:),SV(SE(:,2),:)),1,3));                 % Surface Edge Vector
  SEL = sum(SEv.^2,2).^0.5;                                                   % Surface Edge Length  
  clear SEv

  
  %%
  cdata_l = inf(size(cdata),'single'); 
  cdata_h = zeros(size(cdata),'single'); 
  for i=1:size(SE,1)
    val = (cdata(SE(i,2)) - cdata(SE(i,1)))*SEL(SE(i,1));
    cdata_l(SE(i,1)) = min([cdata_l(SE(i,1)),val]);
    cdata_h(SE(i,1)) = max([cdata_h(SE(i,2)),val]);
  end
  cdata = cdata_h - cdata_l; 
end
              
