function [Yth,S,Psurf,EC,defect_size,res] = cat_surf_createCS2(V,V0,Ym,Ya,YMF,Ytemplate,opt)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth1,S,Psurf,EC]=cat_surf_createCS(V,V0,Ym,Ya,YMF,Ytemplate,opt)
%
% Yth1   .. thickness map
% S      .. structure with surfaces, like the left hemisphere, that contains
%           vertices, faces, GM thickness (th1)
% Psurf  .. name of surface files
% EC     .. Euler characteristics
% defect_size .. size of topology defects
% res    .. intermediate and final surface creation information
% V      .. spm_vol-structure of internally interpolated image
% V0     .. spm_vol-structure of original image
% Ym     .. the (local) intensity, noise, and bias corrected T1 image
% Ya     .. the atlas map with the ROIs for left and right hemispheres
%           (this is generated with cat_vol_partvol)
% YMF    .. a logical map with the area that has to be filled
%           (this is generated with cat_vol_partvol)
% Ytemplate .. Shooting template to improve cerebellar surface
%              reconstruction
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
% $Id: cat_surf_createCS2.m 1561 2020-02-04 15:49:34Z gaser $ 

%#ok<*AGROW,*STREMP,*ASGLU,*SFLD,*STFLD>


  % Turn off gifti data warning in gifti/subsref (line 45)
  %   Warning: A value of class "int32" was indexed with no subscripts specified. 
  %            Currently the result of this operation is the indexed value itself, 
  %            but in a future release, it will be an error. 
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  cstime = clock; 
  
  
  % set debuging variable
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  
  % set defaults
  if ~exist('opt','var'), opt = struct(); end                 % create variable if not exist
  vx_vol        = sqrt(sum(V.mat(1:3,1:3).^2));               % further interpolation based on internal resolution 
  def.verb      = cat_get_defaults('extopts.expertgui');      % 0-none, 1-minimal, 2-default, 3-details, 4-debug
  def.surf      = {'lh','rh'};                                % surface reconstruction setting with {'lh','rh','rc','lc'} and
                                                              % additional fast option 'fst' (e.g., 'lhfst') without topo.-corr. and sphere-reg. 
  
  % reducepatch has some issues with self intersections and should only be used for "fast" option
  % There is a new SPM approach spm_mesh_reduce that is maybe more robust. 
  % Higher resolution is at least required for animal preprocessing that is given by cat_main.
  def.LAB                 = cat_get_defaults('extopts.LAB');  % brain regions 
  def.SPM                 = 0;                                % surface-reconstration based on SPM segmentation input (see cat_main)
  def.pbtlas              = 0;                                % myelination correction option (in development - not working correctly in all data, RD201907)  
  def.pbtmethod           = 'pbt2x';                          % projection-based thickness (PBT) estimation ('pbt2x' (with minimum setting) or 'pbt2')
  def.sharpen             = 1;                                % sharpening function (in development, RD2019)
  def.sharpenCB           = 1;                                % sharpening function for the cerebellum (in development, RD2017-2019)
  def.thick_measure       = 1;                                % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS))
  def.thick_limit         = 5;                                % 5mm upper limit for thickness (same limit as used in Freesurfer)
  def.collcorr            = 5;                                % correction of surface collisions: 0 - none; 1 - CG (extract_pial_white); 2 - CAT_SI, 3 - PBT, 4 - PBT+, 5 - PBT+ & CAT_SI 
  def.surf_measures       = 1;                                % 0 - none, 1 - only thickness, 2 - expert maps (myelin,defects), 3 - developer (WMT,CSFT, ...),
                                                              % 4 - debug output, 5 - debug extended (more substeps and mex output)
  %def.new_release         = 0;                                % developer flag to test new functionality for new release (currently not used)
  %def.WMT                 = 0;                                % pbt-based WM/CSF width/depth/thickness estimation 
  def.fsavgDir            = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.add_parahipp        = cat_get_defaults('extopts.add_parahipp');
  def.scale_cortex        = cat_get_defaults('extopts.scale_cortex');
  def.close_parahipp      = cat_get_defaults('extopts.close_parahipp');
  def.reduce_mesh         = 1;  % 0 - surface creation on PBT resolution, no mesh-reduction (very slow) 
                                % 1 - optimal resolution depending on final mesh resolution, no mesh-reduction 
                                % 2 - internal resolution, no mesh-reduction (slow for highres data) 
                                % 3/4 - SPM/MATLAB reduce on initial surface - there seams to be a bug in the c-function that kills matlab 
                                % 5/6 - SPM/MATALB reduce on initial and final surface   
                                % 7 - call matlab reduce in external matlab
  def.outputpp.native     = 0;  % output of Ypp map for cortical orientation in EEG/MEG 
  def.outputpp.warped     = 0;
  def.outputpp.dartel     = 0;
  
  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.fast                = any(~cellfun('isempty',strfind(opt.surf,'sfst'))) + ...
                            any(~cellfun('isempty',strfind(opt.surf,'fst'))); % fast registration without topo.-corr. and sphere-reg.
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness estimation  
  opt.surf                = cat_io_strrep(opt.surf,{'sfst','fst','v'},'');    % after definition of the 'fast' and 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,1.5]));                  % general limitation of the PBT resolution
  opt.extract_pial_white          = opt.collcorr==1;                          % estimate pial and white matter surface (in development and very slow!)
  opt.force_no_selfintersections  = opt.extract_pial_white;                   % exact estimation requires this setting

  if any(~cellfun('isempty',strfind(opt.surf,'sfst')))
    opt.vdist    = 27/3; 
    opt.pbtres   = 0.8;
  end
  if opt.fast
    opt.collcorr = 3;
  end
  
  % isosurface threshold
  th_initial = 0.5;

  % distance between vertices that can be set directly by "vdist" or indirectly by "interpV"  
  % - surface should have more than 80k faces to look nice, whereas more than 400k does not improve the visual quality 
  % - controlled by power function to avoid a quadratic grow of the number of faces
  % - vdisto = [4 2 1 0.5] => [100k 200k 400k 800k] faces
  %max( 0.5 , min( 2 , min( opt.interpV , mean(vx_vol0) ))); % use square to use sqrt in general  
  if ~isfield(opt,'vdist') || opt.vdist == 0, opt.vdist  = 4/3; end
  opt.vdisto = opt.vdist; 
  opt.vdist  = sqrt(opt.vdist * 2); % here we use the sqrt to support linear mesh resolution increase (default input is [4 2 1 0.5])
  
  % Another parameter to control runtime is the accuracy of the surface
  % deformation. As far as we primary adapt the mesh resolution above, it 
  % is useful to use sqrt values rather than linear values to avoid square
  % processing times for higher quality levels. Otherwise, we can simple 
  % avoid changes here ... so we can define this parameter utilizing the 
  % vdist parameter by simply divide it by 100.
  %def.surfaccuracy        = 0.01; % no adaption here otherwise processing time will not simply double 
  def.surfaccuracy = opt.vdist / 200; 
  def.reduceCS     = (300000 * sqrt(4/3 * 2) ) ./ opt.vdist; % to test ... fprintf('%g ',300000 * sqrt(1.3*2) ./ (( [4 2 1.3 1 0.5] * 2).^0.5))
  opt              = cat_io_updateStruct(def,opt);
  
  if opt.surf_measures > 4, opt.verb = 3; end
  
  % function to estimate the number of interactions of the surface deformation: d=distance in mm and a=accuracy 
  moveth = @(d,a) [ round(d / a) , a ]; 
  QMC    = cat_io_colormaps('marks+',17);
  color  = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  rate   = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
   
  
  % some interal overview for developers
  if opt.verb>2
    fprintf('\nSurface reconstruction:              %s\n',....
      sprintf('%s',char( cellfun(@(x) [x ' '],opt.surf,'UniformOutput',0) )')); 
    fprintf('  PBT resolution:                    %0.3f\n',opt.interpV);
    fprintf('  lower face limit:                  %g\n',opt.reduceCS);
    fprintf('  maximal vertex distance:           %0.3f mm\n',opt.vdist);
    fprintf('  optimization stepsize:             %0.3f mm\n',opt.surfaccuracy);
    fprintf('  collcorr / reduce_mesh:            %d / %d',20 + opt.collcorr,opt.reduce_mesh);
  end  
  
  
  % correction for 'n' prefix for noise corrected and/or interpolated files
  [pp,ff]   = spm_fileparts(V.fname);

  if cat_get_defaults('extopts.subfolders')
    surffolder = 'surf'; 
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

  
  
  %% get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
    
  % noise reduction for higher resolutions (>=1 mm full correction, 1.5 mm as lower limit)
  % (added 20160920 ~R1010 due to severe sulcus reconstruction problems with 1.5 Tesla data)
  Yms = Ym + 0; cat_sanlm(Yms,3,1);
  mf  = min(1,max(0,3-2*mean(vx_vol,2))); 
  Ym  = mf * Yms  +  (1-mf) * Ym;
  clear Yms;
   
  % filling
  Ymf  = max(Ym,min(1,YMF & ~NS(Ya,opt.LAB.HC) & ~( cat_vol_morph( NS(Ya,opt.LAB.HC),'d',2) & NS(Ya,opt.LAB.TH) ))); 
  Ymfs = cat_vol_smooth3X(Ymf,1); 
  Ytmp = cat_vol_morph(YMF,'d',3) & Ymfs>2.3/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs; 
  Ymf = Ymf * 3;
  
  % removing fine WM structures in the hippocampus area to reduce topological and geometrical defects (added RD20190912)
  % use erode to reduce probability of cutting other gyri
  HCmask = cat_vol_morph( NS(Ya,opt.LAB.HC) , 'de', 1.5, vx_vol); 
  Ymf( HCmask ) =  min(2,Ymf( HCmask )); clear HCmask; 

  % surface output and evaluation paramter 
  Psurf = struct(); 
  res   = struct('euler_characteristic',nan,'defect_size',nan,'lh',struct(),'rh',struct()); 
  % initialize WM/CSF thickness/weidth/depth maps
  Yth   = zeros(size(Ymf),'single'); 
  Ypp   = zeros(size(Ymf),'single'); 
  if opt.surf_measures > 3 
    Ywd = zeros(size(Ymf),'single'); 
    Ycd = zeros(size(Ymf),'single'); 
  end
  
  
  
  %% reduction of artifact, blood vessel, and meninges next to the cortex in SPM segmentations 
  %  (are often visible as very thin structures that were added to the WM 
  %  or removed from the brain)
  if ~opt.SPM
    Ydiv  = cat_vol_div(Ymf,vx_vol); 
    Ycsfd = cat_vbdist(single(Ymf<1.8),Ymf>1,vx_vol);
    Yctd  = cat_vbdist(single(Ymf<0.5),Ymf>0,vx_vol); 
    Ysroi = Ymf>2  &  Yctd<10  & Ycsfd>0 & Ycsfd<2.0 & ...
            cat_vol_morph(~NS(Ya,opt.LAB.HC) & ~NS(Ya,opt.LAB.HI) & ...
              ~NS(Ya,opt.LAB.PH) & ~NS(Ya,opt.LAB.VT),'erode',4); 
    if ~debug, clear Ycsfd Yctd; end 
    Ybv   = cat_vol_morph(Ymf+Ydiv./max(1,Ymf)>3.5,'d') & Ymf>2; 
    Ymf(Ybv) = 1.4; 
    Ymfs  = cat_vol_median3(Ymf,Ysroi | Ybv,Ymf>eps & ~Ybv,0.1); % median filter
    Ymf   = min(Ymf, mf * Ymfs  +  (1-mf) * Ymf);
    clear Ysroi Ydiv Ybv

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

  
  
  if 0
    %% Blood vessel correction 
    %  It would be much faster to do the blood vessel correction on the 
    %  internal rather the interpolated resolution (see below)
    %  (not tested now RD201912)
    
    if opt.verb>1, stime = cat_io_cmd('  Blood vessel correction'); end
    Ylt = (1 - (Ymf<1.9) + cat_vol_morph(cat_vol_morph(Ymf>2.2,'lo'),'d'))/2; 
    Ylt = cat_vol_laplace3R(single(Ylt),Ylt>0 & Ylt<1,0.001);

    % update Ywm distance map
    Ymskc = (Ylt<0.8 & Ymf>2) | Ymf<=1.5;
    Ymskc = smooth3(Ymskc); 

    % apply correction and filter modified area
    Ymf  = Ymf.*(1-Ymskc) + Ymskc.*max(1,min(Ymf, 4.2 - cat_vol_smooth3X(Ymf,1))); 
    Ymf  = min(Ymf,cat_vol_median3(Ymf,Ymskc>0));

    if ~debug, clear Ylt Ymskc; end
    if opt.verb>1, fprintf('%5.0fs\n',etime(clock,stime)); end
  end
  
  
  
  
  if 1 
    %% Amygdala Hippocampus smoothing 
    %  We use a median filter to remove the nice details of the hippocampus
    %  that will cause topology errors and self-intersections. 
    %  Currently, I have no CAT ROI for Amygdala - but it would be more 
    %  robust to filter (simple smoothing) this region strongly because 
    %  the "random" details especially in good data introduce more variance.
    %  (RD 201912)
    
    Ymsk = ~(NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.BS) );
    Ymsk = Ymf>0 & cat_vol_morph( NS(Ya,opt.LAB.HC) , 'd' , 3 , vx_vol ) & Ymsk; 
    Ymf  = cat_vol_median3( Ymf , Ymsk ); 
    Ymf  = cat_vol_median3( Ymf , Ymsk ); 
    if ~debug, clear Ymsk; end
    
    % further cleanup
    Ymsk = ~(NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.BS) );
    Ymsk = Ymf>0 & NS(Ya,opt.LAB.HC) & Ymsk; 
    Ymsk = smooth3(Ymsk); 
    Ymf  = min(Ymf,3-Ymsk); 
    if ~debug, clear Ymsk; end
  end
  
  
  
  if opt.sharpen %strcmp(opt.pbtmethod,'pbt2')
    %% Sharpening
    %  This function works quite good in the cerebellum and it allows to
    %  stabilize thin structures avoiding thickness overestimations if it 
    %  used moderatly.  Abuse can increase problems by local artefacts or 
    %  unwanted details that finally cause strong local unterestimation. 
    %  (RD 201911)
    
    gmv = sum(round(Ymf(:))==2) / sum(round(Ymf(:))==3); gmvm = max(0,min(1,1 / gmv)); 
    for i=1:gmv, Ymf = Ymf.*(gmvm) +  (1 - gmvm).*max(1,min(3, Ymf - smooth3(cat_vol_median3(Ymf,Ymf>1,Ymf>1) - Ymf) )); end  
  end
  
  
  
  
  %% Sharpening of thin structures in the cerebellum (gyri and sulci)
  %  Processing of the cerebellum needs a complete updated with a simplified
  %  reconstruction rather than the maximum level of detail. However, the 
  %  main branches are quite similare and their is a lot of variance in 
  %  aging and between MR protocols. 
  %  (RD 2015-201912)
  if ~opt.SPM && opt.sharpenCB && any(~cellfun('isempty',strfind(opt.surf,'cb')))
    if opt.verb>2, fprintf('\n'); stime = cat_io_cmd('  Sharpen cerebellum'); end
    %Ytemplate = Ytemplate * 3; 
    
    Ywmd = cat_vbdist(single(Ymf>2.5),Ymf>1,vx_vol);

    %% normalized Ym divergence
    Ydiv = cat_vol_div(max(2,Ymf)); %Ydivl  = cat_vol_div(Ymf,vx_vol); 
    Ydc  = cat_vol_localstat(abs(Ydiv),Ymf>0,1,3);
    Ydc  = cat_vol_localstat(Ydc,Ymf>0,1,1);
    Ydiv = Ydiv ./ Ydc; clear Ydc;
    Ydiv( abs(Ydiv) > 1.5 ) = 0;
    
    %% normalized Ytemplate divergence
    if ~isempty(Ytemplate)
      Ydivt = cat_vol_div(max(2,Ytemplate));
      Ydct  = cat_vol_localstat(abs(Ydivt),Ytemplate>0,1,3);
      Ydct  = cat_vol_localstat(Ydct,Ymf>0,1,1);
      Ydivt = Ydivt ./ Ydct; clear Ydct;
      Ydivt( abs(Ydivt) > 1.5 ) = 0;
    else
      Ytemplate = Ymf; 
      Ydivt     = Ydiv; 
    end
    
    %% bias-correction based
    % WM 
    Ymsk = ((cat_vol_morph(NS(Ya,opt.LAB.CB),'e',3) | YMF) & ( (Ym-Ydiv).*(Ytemplate/3-Ydivt) )>2/3 ) |  ...
           (NS(Ya,opt.LAB.PH) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01) ) ) | ...                  % hippocampal gyri
           (NS(Ya,opt.LAB.CT) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01 & ...
              Ycsfd>cat_stat_nanmean(Ycsfd(Ycsfd(:)>0 & Ycsfd(:)<100)) )*1.0) );            % distant gyri and sulci in the cerebrum
    Yi   = cat_vol_localstat(Ymf,Ymsk,1,3);
    % GM
    Ymsk = (NS(Ya,opt.LAB.CB) & ( Ymf>1.9 & Ymf<2.2 & Ycsfd>0 & Ydiv>-0.05) ) | ...         % sulci and gyri in the cerebellum 
           (NS(Ya,opt.LAB.PH) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 ) ) | ...                     % hippocampal gyri
           (NS(Ya,opt.LAB.CT) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 & ...
              Ywmd>cat_stat_nanmean(Ywmd(Ywmd(:)>0 & Ywmd(:)<100))*0.2 ) );                 % distant gyri and sulci in the cerebrum
    Yi   = Yi + cat_vol_localstat(Ymf,Yi==0 & Ymsk,1,1)/2*3;%& ( cat_vol_morph(Yi==0,'e') & Ymf>2.2)
    Yi   = cat_vol_localstat(Yi,Yi>0,1,3);
    Yi   = cat_vol_localstat(Yi,Yi>0,1,1); 
    if ~debug, clear Ywmd Ymsk Ydiv; end
    % CSF - instable and not required
    Ywi = cat_vol_approx(Yi,'nn',1,1,struct('lfO',2)); clear Yi
    
    %% only cerebellum
    Ycb = cat_vol_smooth3X( NS(Ya,opt.LAB.CB) , 2 ); 
    Ywi = 3*ones(size(Ywi),'single').*(1-Ycb) + Ycb.*Ywi; clear Ycb 
    Ymf = Ymf./Ywi * 3; 
    
    % denoising result
    Ycb = Ymf .* NS(Ya,opt.LAB.CB); cat_sanlm(Ycb,3,1); Ymf(NS(Ya,opt.LAB.CB)) = Ycb(NS(Ya,opt.LAB.CB)); clear Ycb
    
    
    %% sharpening (RD 201912)
    Ycb = (NS(Ya,opt.LAB.CB)>0.5) .* max(0,min(1,min(2,max(-1,(Ymf/3).^2 - 0.1*Ydiv) .* max(-1,Ytemplate/3 - 0.02*Ydivt - 0.03*Ydiv )*3 - 1)/3 + 2/3)); 
    for i=1:3, Ycb = max(0,min(1, Ycb - smooth3(cat_vol_median3(Ycb,Ycb>0,Ycb>0) - Ycb) )); end
    Ycb = min(1,Ycb); 
    cat_sanlm(Ycb,3,1); 
    
    %% final mixing
    Ymsk = cat_vol_smooth3X(NS(Ya,opt.LAB.CB) & Ycb<1.9/3,0.5); 
    Ycb  = Ycb.*(1-Ymsk) + Ymsk.*Ymf/3;
    Ycs  = NS(Ya,opt.LAB.CB) .* cat_vol_smooth3X( NS(Ya,opt.LAB.CB) , 4 ); 
    Ymf  = Ymf.*(1-Ycs) + Ycs.*Ycb*3; clear Ycs 
    
    if ~debug, clear Ywi Yi Ycb; end
    if ~debug, clear Ymsk; end
    if opt.verb>2, fprintf('%5.0fs\n',etime(clock,stime)); end
  end
  if ~debug, clear Ydiv Ycsfd; end
  
  
  
  % complete atlas map 
  [D,I] = cat_vbdist(single(Ya>0)); Ya = Ya(I); clear D;  
  
  
  % use sum of EC's and defect sizes for all surfaces, thus set values initially to 0
  EC            = 0;
  defect_size   = 0;
  defect_area   = 0; 
  defect_number = 0; 

  
  % main loop for each surface structure 
  for si=1:numel(opt.surf)
   
    % surface filenames
    Pm         = fullfile(pp,mrifolder, sprintf('m%s',ff));    % raw
    Praw       = fullfile(pp,surffolder,sprintf('%s.central.nofix.%s.gii',opt.surf{si},ff));    % raw
    Psphere0   = fullfile(pp,surffolder,sprintf('%s.sphere.nofix.%s.gii',opt.surf{si},ff));     % sphere.nofix
    Pcentral   = fullfile(pp,surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % central
    Pcentralr  = fullfile(pp,surffolder,sprintf('%s.central.resampled.%s.gii',opt.surf{si},ff));% central
    Player4    = fullfile(pp,surffolder,sprintf('%s.layer4.%s.gii',opt.surf{si},ff));           % layer4
    PintL4     = fullfile(pp,surffolder,sprintf('%s.intlayer4.%s',opt.surf{si},ff));            % layer4 intensity
    Ppial      = fullfile(pp,surffolder,sprintf('%s.pial.%s.gii',opt.surf{si},ff));             % pial (GM/CSF)
    Pwhite     = fullfile(pp,surffolder,sprintf('%s.white.%s.gii',opt.surf{si},ff));            % white (WM/GM)
    Pthick     = fullfile(pp,surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % FS thickness / GM depth
    Ppbt       = fullfile(pp,surffolder,sprintf('%s.pbt.%s',opt.surf{si},ff));                  % PBT thickness / GM depth
    Pgwo       = fullfile(pp,surffolder,sprintf('%s.depthWMo.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgw        = fullfile(pp,surffolder,sprintf('%s.depthGWM.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgww       = fullfile(pp,surffolder,sprintf('%s.depthWM.%s',opt.surf{si},ff));              % gyrus witdh of the WM / WM depth
    Pgwwg      = fullfile(pp,surffolder,sprintf('%s.depthWMg.%s',opt.surf{si},ff));             % gyrus witdh of the WM / WM depth
    Psw        = fullfile(pp,surffolder,sprintf('%s.depthCSF.%s',opt.surf{si},ff));             % sulcus width / CSF depth / sulcal span
    Pdefects0  = fullfile(pp,surffolder,sprintf('%s.defects0.%s',opt.surf{si},ff));             % defects temporary file
    Pdefects   = fullfile(pp,surffolder,sprintf('%s.defects.%s',opt.surf{si},ff));              % defects
    Psphere    = fullfile(pp,surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    Pspherereg = fullfile(pp,surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    Pfsavg     = fullfile(opt.fsavgDir, sprintf('%s.central.freesurfer.gii',opt.surf{si}));     % fsaverage central
    Pfsavgsph  = fullfile(opt.fsavgDir, sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));      % fsaverage sphere    
    
    
    % add the variables defined in "surffile" to the "Psurf" output variable
    surffile = {'Praw','Psphere0','Pcentral','Pthick','Ppbt','Pgw','Pgww','Psw',...
      'Pdefects0','Pdefects','Psphere','Pspherereg','Pfsavg','Pfsavgsph','Pwhite','Ppial'};
    for sfi=1:numel(surffile)
      eval(sprintf('Psurf(si).%s = %s;',surffile{sfi},surffile{sfi})); 
    end
    
    
    %% reduce for object area
    switch opt.surf{si}
      case {'lh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==1); Yside = mod(Ya,2)==1; 
      case {'rh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==0); Yside = mod(Ya,2)==0;  
      case {'cb'},  Ymfs = Ymf .* (Ya>0) .*   NS(Ya,opt.LAB.CB); Yside = true(size(Ymfs)); % new full cerebellum reconstruction
    end 
   
    % check for cerebellar hemis 
    iscerebellum = strcmp(opt.surf{si},'cb');
    
    
    if ~iscerebellum
      mask_parahipp = NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.HC);
    end
    
    % print something
    if si==1, fprintf('\n'); end
    switch opt.fast
      case 2, fprintf('%s - fast with registration:\n',opt.surf{si});
      case 1, fprintf('%s - fast without registration:\n',opt.surf{si});
      case 0, fprintf('%s:\n',opt.surf{si});
    end
    
    % removing background (smoothing to remove artifacts)
    switch opt.surf{si}
      case {'lh','rh'},  [Ymfs,Ysidei,mask_parahipp,BB] = cat_vol_resize({Ymfs,Yside,mask_parahipp},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
      case {'cb'},       [Ymfs,Ysidei,BB] = cat_vol_resize({Ymfs,Ysidei},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
    end
    
    % interpolation 
    imethod         = 'cubic'; % cubic should be better in general - however, linear is better for small thickness (version?)
    [Ymfs,resI]     = cat_vol_resize(max(1,Ymfs),'interp',V,opt.interpV,imethod);                  % interpolate volume
    Ysidei           = cat_vol_resize(Ysidei,'interp',V,opt.interpV,imethod)>0.5;                    % interpolate volume (small dilatation)
    
    if ~iscerebellum
      % get dilated mask of gyrus parahippocampalis and hippocampus of both sides
      mask_parahipp = smooth3( cat_vol_morph(mask_parahipp,'dd',6) );
      mask_parahipp = cat_vol_resize(mask_parahipp,'interp',V,opt.interpV)>0.5;          % interpolate volume
    end 
    
    % PVE with background will lead to a light underestimation?
    Ymfs = min(3,max(1,Ymfs));

    
    
    %% cleanup
    if 1
    %% new cleanup (201911) ... very slow ... ~45 seconds
    %  * this support closer values for the different resolutions of the thickness phantom 
    %  * possibly to aggressive
    %  * better on original resolution 
    %  * use of Yp0 map? 
    %  * region to cerebellum showed problems at least in the SIMON2018 test-retest dataset
  
      if opt.verb>1, stime = cat_io_cmd('  Cleanup'); end
      Ymfs = Ymfs .* (cat_vol_morph( cat_vol_morph( Ymfs>1.2,'ldo', 3.0 / opt.interpV / (iscerebellum+1) ), 'd', 0.5/opt.interpV));
      Ymfs = Ymfs .* (cat_vol_morph( cat_vol_morph( Ymfs>1.5,'ldo', 2.5 / opt.interpV / (iscerebellum+1) ), 'd', 0.5/opt.interpV));
      Ymfs = Ymfs .* (cat_vol_morph( cat_vol_morph( Ymfs>1.9 - 0.2*(iscerebellum+1),'ldo', 2 / opt.interpV / (iscerebellum+1) ), 'd',2 / opt.interpV * (iscerebellum+1) ));
      Ymfs = max(1,Ymfs);

      if opt.verb>1, fprintf('%5.0fs\n',etime(clock,stime)); end
    end
    
    
    
    if opt.sharpen
      %% new sharpening function also for the cerebrum (201912)
      %  * this function clearly helps in case of blurred WM and CSF
      %    structures but it also strongly change the thickness values! 
      gmv = sum(round(Ymfs(:))==2) / sum(round(Ymfs(:))==3); gmvm = max(0,min(1,1 / gmv)); 
      for i=1:gmv, Ymfs = Ymfs.*(gmvm) +  (1 - gmvm).*max(1,min(3, Ymfs - smooth3(cat_vol_median3(Ymfs,Ymfs>1,Ymfs>1) - Ymfs) )); end     
    end
    if opt.sharpenCB
      %% new sharpening function also for the cerebrum (201912)
      %  * this function clearly helps in case of blurred WM and CSF
      %    structures but it also strongly change the thickness values!
      %gmv = max(3,min(9,6 * sum(round(Ymfs(:))==2) / sum(round(Ymfs(:))==3))); gmvm = max(0,min(1,1 / gmv)); 
      gmv = max(1,min(3,2 * sum(round(Ymfs(:))==2) / sum(round(Ymfs(:))==3))); gmvm = max(0,min(1,1 / gmv)); 
      for i=1:gmv, Ymfs = Ymfs.*(gmvm) +  (1 - gmvm).*max(1,min(3, Ymfs - smooth3(cat_vol_median3(Ymfs,Ymfs>1,Ymfs>1) - Ymfs) )); end     
    end  
    
    
    
    % blood vessel correction
    if ~iscerebellum
    %% new blood vessel correction (201911)
    %  * use Laplace-filter to remove blood-vessels (thin structures lose more energy) 
    %  * works quite well but it should fit better in blood vessel correction of cat_main  
    %  * optimized for speed ~15 seconds
    %  * relative save 
      
      if opt.verb>1, stime = cat_io_cmd('  Blood vessel correction'); end
      Ylt = (1 - (Ymfs<1.9) + cat_vol_morph(cat_vol_morph(Ymfs>2.2,'lo'),'d'))/2; 
      Ylt = cat_vol_laplace3R(single(Ylt),Ylt>0 & Ylt<1,0.01);

      % update Ywm distance map
      Ymskc = (Ylt<0.5 & Ymfs>2) | Ymfs<=1.5;
      Ymskc = smooth3(single(Ymskc)); 
      
      % apply correction and filter modified area
      Ymfs  = Ymfs.*(1-Ymskc) + Ymskc.*max(1,min(Ymfs, 4.2 - cat_vol_smooth3X(Ymfs,1))); 
      Ymfs  = cat_vol_median3(Ymfs,Ymskc>0);
      
      if ~debug, clear Ylt Ymskc; end
      if opt.verb>1, fprintf('%5.0fs\n',etime(clock,stime)); end
    end
    
    

    
    
    %% pbt calculation
    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1'))); stimet =stime;
    if strcmp(opt.pbtmethod,'pbt3') %|| opt.collcorr>0
      [Yth1i,Yppi] = cat_vol_pbt3(Ymfs,struct('method',opt.pbtmethod,'cb',iscerebellum,'resV',opt.interpV,...
        'vmat',V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1])); % avoid underestimated thickness in gyri
    else %opt.pbtmethod
      [Yth1i,Yppi,Ymfs] = cat_vol_pbt(Ymfs,struct('method',opt.pbtmethod,'resV',opt.interpV,'vmat',...
        V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'pbtlas',opt.pbtlas)); % avoid underestimated thickness in gyri
    end  
    
    
    % back to internal resolution (at least for thickness)
    Yth1i(Yth1i>10) = 0; Yppi(isnan(Yppi)) = 0;                            % general thickness limit  
    [D,I] = cat_vbdist(Yth1i,Ysidei); Yth1i = Yth1i(I); clear D I;         % add further values around the cortex
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI);                         % back to original resolution
    Yth1t = cat_vol_resize(Yth1t,'dereduceBrain',BB);                      % adding background
    Yth  = max(Yth,Yth1t .* Yside);                                      % save on main image
%    clear Yth1t Ysidei;
fprintf('%26.22f | %26.22f | %26.22f\n',std(Ymfs(:)),std(Yth1i(:)),std(Yppi(:))); 

    
    %% PBT estimation of the gyrus and sulcus width 
    if opt.surf_measures > 3 
      %% gyrus width / WM depth
      %  For the WM depth estimation it is better to use the L4 boundary
      %  and correct later for thickness, because the WM is very thin in
      %  gyral regions and will cause bad values. 
      %  On the other side we do not want the whole filled block of the 
      %  Yppi map and so we have to mix both the original WM map and the
      %  Yppi map. 
      %  As far as there is no thickness in pure WM regions there will
      %  be no correction. 
      
      stime = cat_io_cmd('  WM depth estimation','g5','',opt.verb-1,stime); 
      [Yar,Ymr] = cat_vol_resize({Ya,Ym},'reduceBrain',vx_vol,BB.BB);       % removing background
      Yar   = uint8(cat_vol_resize(Yar,'interp',V,opt.interpV,'nearest'));  % interpolate volume
      Ymr   = cat_vol_resize(Ymr,'interp',V,opt.interpV);                   % interpolate volume
      switch opt.surf{si}
        case {'lh'} 
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==1);
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==1)); 
        case {'rh'}
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==0);    
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==0)); 
        case {'cb'}
          Ymr = Ymr .* (Yar>0) .* NS(Yar,3);
          Ynw = true(size(Ymr)); 
      end 
      clear Yar; 
     
      %
      Yppis = Yppi .* (1-Ynw) + max(0,min(1,Ymr*3-2)) .* Ynw; clear Ynw;              % adding real WM map 
      Ywdt  = cat_vol_eidist(1-Yppis,ones(size(Yppis),'single'));                     % estimate distance map to central/WM surface
      Ywdt  = cat_vol_pbtp(max(2,4-Ymfs),Ywdt,inf(size(Ywdt),'single'))*opt.interpV;
      [D,I] = cat_vbdist(single(Ywdt>0.01),Yppis>0); Ywdt = Ywdt(I); clear D I Yppis; % add further values around the cortex
      Ywdt  = cat_vol_median3(Ywdt,Ywdt>0.01,Ywdt>0.01);                    
      Ywdt  = cat_vol_localstat(Ywdt,Ywdt>0.1,1,1);                                   % smoothing
      Ywdt  = cat_vol_resize(Ywdt,'deinterp',resI);                                   % back to original resolution
      Ywdt  = cat_vol_resize(Ywdt,'dereduceBrain',BB);                                % adding background
      Ywd   = max(Ywd,Ywdt); 
      clear Ywdt;
      
      % sulcus width / CSF depth
      %  for the CSF depth we cannot use the original data, because of
      %  sulcal blurring, but we got the PP map at half distance and
      %  correct later for half thickness
      stime = cat_io_cmd('  CSF depth estimation','g5','',opt.verb-1,stime); 
      YM    = single(smooth3(cat_vol_morph(Ymr<0.1,'o',4))<0.5); YM(YM==0)=nan;       % smooth CSF/background-skull boundary 
      Yppis = Yppi .* ((Ymr+0.25)>Yppi) + min(1,Ymr*3-1) .* ((Ymr+0.25)<=Yppi);       % we want also CSF within the ventricle (for tests)
      Ycdt  = cat_vol_eidist(Yppis,YM);                                               % distance to the cental/CSF-GM boundary
      Ycdt  = cat_vol_pbtp(max(2,Ymfs),Ycdt,inf(size(Ycdt),'single'))*opt.interpV; Ycdt(isnan(Ycdt))=0;
      [D,I] = cat_vbdist(single(Ycdt>0),Yppis>0 & Yppis<3); Ycdt = Ycdt(I); clear D I Yppis; % add further values around the cortex
      Ycdt  = cat_vol_median3(Ycdt,Ycdt>0.01,Ycdt>0.01);                              % median filtering
      Ycdt  = cat_vol_localstat(Ycdt,Ycdt>0.1,1,1);                                   % smoothing
      Ycdt  = cat_vol_resize(Ycdt,'deinterp',resI);                                   % back to original resolution
      Ycdt  = cat_vol_resize(Ycdt,'dereduceBrain',BB);                                % adding background
      Ycd   = max(Ycd,Ycdt); 
      clear Ycdt;
      
      clear Ymr;
    end
    %if debug, Yppio=Yppi; end
    fprintf('%5.0fs\n',etime(clock,stime));
    
    
    
    % Replace isolated voxels and holes in Ypp by its median value
    
    % indicate isolated holes and replace by median of the neighbors
    Yppi(Yppi<0.35 & ~cat_vol_morph(Yppi<1,'l'))=1;  % close major wholes in the WM 
    Ymsk = Yppi==0 & cat_vol_morph(Yppi>0.9,'d',1); % filter small wholes close to the WM
    Yppi = cat_vol_median3(single(Yppi),Ymsk,~Ymsk); 
    
    % indicate isolated objects and replace by median of the neighbors
    Yppi(Yppi>0.65 & cat_vol_morph(Yppi==0,'l'))=0;
    Ymsk = Yppi>0.95 & cat_vol_morph(Yppi<0.1,'d',1); 
    Yppi = cat_vol_median3(single(Yppi),Ymsk,~Ymsk);
    if ~debug, clear Ymsk; end
    
    %  Write Ypp for final deformation
    %  Ytt has the interal resolution (e.g. <1 mm) and is only used temparary 
    %  Write Yppi file with 1 mm resolution for the final deformation, 
    %  because CAT_DeformSurf achieved better results using that resolution
    Yppt = cat_vol_resize(Yppi,'deinterp',resI);                        % back to original resolution
    Yppt = cat_vol_resize(Yppt,'dereduceBrain',BB);                     % adding of background

    % scale Yppt so that backgrounds remains 0 and WM 1, but cortical band is now in the range of 0.1..0.9
    if opt.extract_pial_white
      % mask hemispheres and regions
      indi = find((Yppt>0) & (Yppt<0.99999));
      Yppt(indi) = 0.1 + (0.8*Yppt(indi));
      clear indi
    end
    Ypp = max(Ypp,Yppt .* Yside); 
    
    % this is the map that we want to keep in the original image resolution
    Vpp  = cat_io_writenii(V0,Ypp,'mri',sprintf('%s.pp',opt.surf{si}) ,'percentage position map','uint8',[0,1/255],...
      min([1 1 2],[1 opt.outputpp.warped opt.outputpp.dartel]),opt.trans);
    Vpp  = Vpp(1); 
    
    % just an internal map with 1 mm resolution 
    Vppt          = cat_io_writenii(V,Yppt,'',sprintf('%s.ppt',opt.surf{si}) ,'percentage position map','uint8',[0,1/255],[1 0 0]); clear Yppt
    Vpp1          = Vppt(1);
    %{
    defres        = 1 * ones(1,3); % not tested 
    Vpp1.fname    = fullfile(pp,mrifolder,['pp1' ff '.nii']);
    vmat2         = spm_imatrix(Vpp1.mat);
    Vpp1.dim(1:3) = round(Vpp1.dim .* round(abs(vmat2(7:9)*(1 + iscerebellum))./defres) );   % use double resolution in case of cerebellum
    vmat2(7:9)    = sign(vmat2(7:9)) .* defres./(1 + iscerebellum);            % use double resolution in case of cerebellum
    Vpp1.mat      = spm_matrix(vmat2);
    Vpp1 = spm_create_vol(Vpp1); 
    for x3 = 1:Vpp1.dim(3)
      M    = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1]) * inv(Vpp1.mat) * Vppt.mat); %#ok<MINV>
      v    = spm_slice_vol(Vppt,M,Vpp1.dim(1:2),1);       
      Vpp1 = spm_write_plane(Vpp1,v,x3);
    end
    if exist(Vppt.fname,'file'), delete(Vppt.fname); end
    %}
    clear M v x3 Vppt; 

    if opt.vol
      S = struct(); Psurf = '';
      if opt.verb<2, fprintf('%5.0fs',etime(clock,stime)); end
      continue; 
    end
    
    
    
    
    %% surface creation 
    %  --------------------------------------------------------------------
    %  Surface create should be at 0.5 mm to support a useful description
    %  in even narrow sulci, i.e., for 1 mm thickness the sulci would be 
    %  2 mm wide were a 1 mm percentage position map is strongly limited.
    %  The surface is reconstructed by marching cubes on a binary version 
    %  of the PBT radial position map Yppi to use a simple voxel-based
    %  topology correction to avoid small defects and incorrect
    %  triangulation (e.g., matlab isosurface function). Next, the
    %  complexity of the surface is reduced by the spm_mesh_reduce function
    %  (more accurate and stable than the matlab reduce function?) to
    %  improve the performance of the following main topology correction. 
    %  However, the reduction can partially lead to very large faces that
    %  need an additional refine. 
    %  Voxelbased topology optimization is important for the hippo. gyrus. 
    %  --------------------------------------------------------------------
    if opt.verb==1, fprintf('%s %4.0fs\n',repmat(' ',1,66),etime(clock,stimet)); end
    nosurfopt = iscerebellum; 
    if opt.verb>1 
      stime = clock;
      switch sprintf('%d%d',nosurfopt,opt.reduce_mesh==1 || opt.reduce_mesh==2)
        case '00'
          cat_io_cprintf('g5',sprintf('  Create topo. opt. HR surface '));
        case '01'
          cat_io_cprintf('g5',sprintf('  Create topo./res. opt. surface '));
        case '10'
          cat_io_cprintf('g5',sprintf('  Create intial HR surface '));
        case '11'
          cat_io_cprintf('g5',sprintf('  Create res. opt. surface '));
      end  
    else
      stime = cat_io_cmd('  Create initial surface','g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
    end
    
    % surface coordinate transformation matrix
    matI              = spm_imatrix(V.mat); 
    matI(7:9)         = sign( matI(7:9))   .* repmat( opt.interpV , 1 , 3); 
    matiBB            = spm_imatrix(V.mat   * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    matIBB            = matiBB; 
    matIBB(7:9)       = sign( matiBB(7:9)) .* repmat( opt.interpV , 1 , 3); 
    Smat.matlabi_mm   = V.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];                % CAT internal space
    Smat.matlabI_mm   = spm_matrix(matI) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % PBT interpolated space
    Smat.matlabIBB_mm = spm_matrix(matIBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated
    Smat.matlabiBB_mm = spm_matrix(matiBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated
 
    
  
    %% scaling correction to reduce topology errors in the parahippocampus
    % ### not tested for cerebellum yet (RD201909) ###
    if ~iscerebellum
      scale_cortex = opt.scale_cortex; %/0.7 * 0.9;
      
      ind0   = find(Yppi<=0);
      Yppisc = scale_cortex * Yppi;
      
      % smooth mask to have smooth border
      mask_parahipp_smoothed = zeros(size(mask_parahipp));
      spm_smooth(double(mask_parahipp),mask_parahipp_smoothed,[4 4 4] / opt.interpV);
      Yppisc  = Yppisc + opt.add_parahipp / scale_cortex * mask_parahipp_smoothed;
      Yppisc(ind0) = 0; clear ind0; %#ok<FNDSB>
      
      % optionally apply closing inside mask for parahippocampal gyrus to get rid 
      % of the holes that lead to large cuts in gyri after topology correction
      if opt.close_parahipp
        tmp = cat_vol_morph(Yppisc,'labclose',1);
        Yppisc(mask_parahipp) = tmp(mask_parahipp); 
        if ~debug, clear tmp; end 
      end
      clear mask_parahipp; 
    else
      scale_cortex = 0.7;
      Yppisc = scale_cortex * Yppi;
    end
    
    if 1
      %% thickness depending cortical scaling - this seams to work but need further tests (RD201911)
      Yth1i  = cat_vol_localstat(Yth1i,Yth1i>0,2,2);
      Yts    = cat_vol_approx(Yth1i,2);  
      Yts    = 1 + max(-0.5,min(0.5,(Yts - mean(Yth1i(:))) / (2 * mean(Yth1i(:)))  )); 
      Yts    = Yts .* (1-mask_parahipp_smoothed) + mask_parahipp_smoothed;  % no thickness adaption in the hippocampus!  
      Yppisc = max(0.55 .* (Yppi>=1),min(1.5, Yppisc .* Yts )); % factor 1 is Ypp at 0.5 > limit 0.55 and 1.5 
      scale_cortex = scale_cortex * median( Yts(:) );
    end
    clear mask_parahipp_smoothed;
    
    % we use another scaling and can therefore update this map
    Yppisc = Yppisc .* Ymfs/2; 
    Yppisc(smooth3(Yppisc)<0.3) = 0; 
    Yppisc(smooth3(Yppisc)>0.7) = 1; 
    Yppisc( Yppisc<0.5 & ~cat_vol_morph(Yppisc<0.5,'l')) = 1;  % close major wholes in the WM 
    Yppisc( Yppisc>0.5 & ~cat_vol_morph(Yppisc>0.5,'l')) = 0;  % remove small dots

   
    %%
    % Marching cubes surface creation and correction for the boundary box 
    % used within the surface creation process.  It is better to use the 
    % full voxel resolution in combination with surface reduction rather   
    % then using lower voxel resolutions! Moreover, it was NOT necessary to 
    % use surface deformation before the surface reduction!
    % #####
    %   I am not sure if the topologoy correction is optimal.
    %   Moreover, a correction should also change the Yppi to avoid self-intersections. 
    %   Maybe a smooth adaption similar to "mask_parahipp_smoothed" can be used here. 
    %   However, this is quite complex and I miss the time go on ... 
    %   RD201911
    % #####
    if iscerebellum
      %% region-growing
      % Ylt = single( cat_vol_morph(Yppi<0.1,'l') +  2 * cat_vol_morph(Yppi>0.9,'l') );
      Ylt = single( cat_vol_morph(Ymfs<1.9,'l') +  2 * cat_vol_morph(Ymfs>2.5,'l') );
      [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.05); Ylt(D>10) = 0;  % Yppi is to similiar
      [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.10); Ylt(D>10) = 0; 
      [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.20); Ylt(D>10) = 0; 
      [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.50); Ylt(D>10) = 0; 
      Ylts    = smooth3(Ylt);
      if ~debug, clear Ylt; end

      %% relative position map between core and hull
      %  to control opening and closing operations to reduce topology defects
      %  the topology correction is critical part where we need a quite
      %  robust inner core to avoid superlarge wholes 
      
      % the hull is quite simple in general 
      Ycbh  = cat_vol_morph(Yppi>0.1 & Ymfs>1.9,'lc',4); 
      % the core is more complicated because we have assure that something
      % is there (at least the highest 20% of the hull distance)
      Ycbhd = cat_vbdist( single( ~Ycbh ) ); 
      Ycbhd = Ycbhd/max(Ycbhd(Ycbhd<1000)); 
      Ycbc  = cat_vol_smooth3X(Yppi + Ymfs/3 + (Ycbhd*0.8) ,2)>2; 
      Ycbc  = cat_vol_smooth3X( cat_vol_morph( cat_vol_morph( cat_vol_morph( Ycbc ,'o',1) ,'l',[2 0.2]), 'lc',2),2) >0.5;
      %% more over we can use the Laplace filter for subtile openen and closing
      if debug, tic; end
      Yltw  = (Yppi>0.6 | cat_vol_morph(Yppi>0.9,'dc',1) | Ycbc)/2 + Ycbc/2; 
      Yltw  = cat_vol_laplace3R(single(Yltw),Yltw>0 & Yltw<1,0.01);
      Yltw  = max(Yltw,Ycbc | cat_vol_morph( cat_vol_morph( Yltw>0.25 , 'do' , 1.5 ) , 'l' , [2 0.2])); 
      Yltw  = cat_vol_laplace3R(single(Yltw),Yltw>0 & Yltw<1,0.001);
      if debug, toc; end
      %%
      if debug, tic; end
      Yltc  = Ycbh/2 + (Yltw>0.002)/2;
      Yltc  = cat_vol_laplace3R(single(Yltc),Yltc>0 & Yltc<1,0.01);
      Yltw(Yltc>0.99) = 1; 
      if debug, toc; end
      %%
      Ycbc  = cat_vol_smooth3X( cat_vol_smooth3X( Yltw , 4)> 0.6 & Yltw>0.5 ,2)>.5 ;
      Yppisc = max(0.55 .* (Yppi>=1),min(1.5, Yppisc .* Yts )); % factor 1 is Ypp at 0.5 > limit 0.55 and 1.5 
      if debug, clear Yts; end
      
      %Ycbc  = cat_vol_morph( Ycbc , 'd') & Yppi>0.5;
      %% distance mapping
      if 1
        Ycbhd = cat_vbdist( single( ~Ycbh ) , ~Ycbc ); 
        Ycbcd = cat_vbdist( single(  Ycbc ) ,  Ycbh ); 
      else
        Ycbhd = single( ~Ycbh ); Ycbhd(  Ycbc ) = nan; [Ytmp,Ycbhd] = cat_vol_downcut( Ycbhd , 4-Ymfs , 0.5);
        Ycbcd = single(  Ycbc ); Ycbcd( ~Ycbh ) = nan; [Ytmp,Ycbcd] = cat_vol_downcut( Ycbcd , Ymfs   , 0.5);
      end
      Ycbpp = min(Ycbh,max(Ycbc,Ycbhd ./ max(eps,Ycbhd+Ycbcd))); 
      Ycbpp = cat_vol_localstat(Ycbpp,Ycbpp>0,1,1); 
      Ycbcd = Ycbcd ./ mean(Ycbcd(Ycbcd(:)>0 & Ycbcd(:)<1000));
      Ycbhd = Ycbhd ./ mean(Ycbhd(Ycbhd(:)>0 & Ycbhd(:)<1000));
      if ~debug, clear Ycbc Ycbh Ycbhd; end
      
      %% 
      Ycbth  = @(lth,pth,mth,cth) max(Ylt,Ylts)>(1.5*lth) & (Yppi - ( 0.5 - Ycbpp )/2 )>pth & Ymfs>mth & ...
                 ( (cth>=0).*(Ycbpp>cth)  | (cth<0).*(Ycbpp<-cth) ); 
      Ycbth2 = @(lth,pth,mth,cth) max(Ylt,Ylts)>(1.5*lth) & (Yppi - ( 0.5 - Ycbpp )/2 )>pth & Ymfs>mth & ...
                 ( (cth>=0).*(Ycbcd<cth)  | (cth<0).*(Ycbcd>-cth) ); 
      if ~debug, clear Ylt Ycbpp Ycbcd; end
      if 0
        % this is to complex ...
        Ycbm = Ycbpp>0.99 | ...
               Ycbth(1,0.95,2.5,0.0) | ...
               Ycbth(1,0.9,2.1,1/3) | ...
               Ycbth(1,0.95,2.5,1/3) | ...
               Ycbth(0.8,0.05,2,2/3) | ...
               (cat_vol_morph( Ycbth(1,0.9,2.2,-1/3) ,'do',2) ) | ...
               (cat_vol_morph(Ycbth2(1,0.95,2.7,2/3) ,'dc',3) & Ycbth2(0,0.3,1.7,1/2)) | ... 
               (cat_vol_morph( Ycbth(1,0.9,2.2,0)   ,'dc',2) & Ycbth(0,0.2,2,2/3)); 
        %%     
        Ycbm(smooth3(Ycbm)<0.3) = 0; Ycbm = single(Ycbm); %Ycbm(Ycbpp==0) = nan;  
        [Ycbm,D] = cat_vol_simgrow(Ycbm,Ymfs+Yppi,0.01); Ycbm(D>2) = 0; 
        Ycbm = Ycbm | cat_vol_morph( Ycbm ,'dc',1) & (Ycbth(1,0.5,2.4,0) | Ycbth(1,0.5,2.1,0.5)); 
        Ycbm(smooth3(Ycbm)<0.5) = 0; Ycbm(smooth3(Ycbm)>0.5) = 1; Ycbm(smooth3(Ycbm)>0.3 & Ycbth(1,0,2,0.8)) = 1;
        Ycbm = Ycbm .* Ycbth(0,0,0,0.5) + cat_vol_morph(Ycbm,'do',1) .* Ycbth(0,0,0,-0.5); 
        Ycbm = Ycbm .* Ycbth(0,0,0,0.3) + cat_vol_morph(Ycbm,'do',3) .* Ycbth(0,0,0,-0.2); 
        Ycbm = cat_vol_morph( Ycbm, 'l'); 
        Ycbm = single( Ycbm ); 
      else
        %%
        Ycbm   = Yppi .* min( 0.4 + Ycbhd/5, 0.5 + Ycbpp*0.60); % in sum larger than one to have save core structure
        if 1
          Ycbms  = smooth3(Ycbm);
          % opening by smoothing in outer regions
          Ycbppt = max(0,min(1,Ycbpp * 5)); 
          Ycbm   = Ycbms .* (1-Ycbppt) + Ycbm .* Ycbppt; 
          % closing by smoothing in inner regions
          Ycbppt = max(0,min(1,(1 - Ycbpp ) * 5));
          Ycbm   = Ycbms .* (1-Ycbppt) + Ycbm .* Ycbppt; 
          clear Ycbms; 
        end
      end
      %%
      
      
      %%
      Yppi05c = Ycbm; evalc(sprintf('clear CS; [Yppi05c,CS.faces,CS.vertices] = cat_vol_genus0(Ycbm,0.5,nosurfopt);')); % no_adjustment 
      [Yvxdef,defect_number0] = spm_bwlabel( double(abs(Yppi05c - (Yppisc>0.5))>0) ); clear Yppi05c;
    else
      % Main initial surface creation using cat_vol_genus0
      % cat_vol_genus0 uses a "simple" marching cube without use of isovalues
      % that is used in the MATLAB isosurface function. Our test showed that
      % the surface deformation allows the same or better accuracy and also
      % that the meshes of cat_vol_genus0 are more regular and also allow 
      % voxel-based topology optimization.  
      
      if opt.reduce_mesh==0 || opt.reduce_mesh>2 % full resolution  
        
        [Yppi05c,CS] = cat_vol_genus0opt(Yppisc,th_initial,15 * (1-iscerebellum),debug);
        
        %Yppi05c = Yppisc; evalc(sprintf('clear CS; [Yppi05c,CS.faces,CS.vertices] = cat_vol_genus0(Yppisc,0.5,nosurfopt);')); % no_adjustment
        [Yvxdef,defect_number0] = spm_bwlabel( double(abs(Yppi05c - (Yppisc>0.5))>0) ); clear Yppi05c;
      else % lower resolutions
        rf = min(1.,max(0.75,opt.vdist)); % because I use V
        VI = V; VI.mat = spm_matrix(matI); 

        %{
        Yppiscr = Yppisc;  
        if 0
          evalc(sprintf('Yppiscr = cat_vol_genus0(Yppisc,0.5,1);'));
        else
          evalc(sprintf('Yppiscr = cat_vol_genus0(Yppisc,0.5,nosurfopt);'));

          % remove larger corrections
          Yvxcorr = abs(Yppiscr - (Yppisc>0.5))>0; 
          Yvxdef  = spm_bwlabel( double( Yvxcorr ) ); clear Yppiscrc; 
          Yvxdef  = cat_vol_morph(Yvxdef,'l',[inf 30])>0; % only small corrections 
          if debug, fprintf('  # vx. of genus-topocorr: %d; final used:  %0.2f%%\n', ...
              sum(Yvxcorr(:)) , 100 * sum(Yvxcorr(:) & ~Yvxdef(:)) / sum(Yvxcorr(:)) ); end
          Yppiscr = Yppiscr & ~Yvxdef; 
          clear Yvxcorr
        end
%}
        Yppiscr = cat_vol_genus0opt(Yppisc,th_initial,15 * (1-iscerebellum),debug);
        
        if 1
          % optimized downsampling to avoid blurring of thin gyri/sulci 
          if 0
            Yppi_mn  = cat_vol_localstat( Yppi , Ymfs>0 , 1 , 2); 
            Yppi_mx  = cat_vol_localstat( Yppi , Ymfs>0 , 1 , 3); 

            Yppi_mn  = smooth3(Yppi_mn); 
            Yppi_mx  = smooth3(Yppi_mx); 
          end

          d = max(1,rf / opt.interpV / 2) * 1.5; 
          distmorph = 1; if distmorph, dm = 'd'; else, dm = ''; end
          Yppi_o  = cat_vol_morph(Yppiscr>0.5,[dm 'o'], d ); % rf / opt.interpV - 1 
          Yppi_c  = cat_vol_morph(Yppiscr<0.5,[dm 'o'], d ); 

          Yppi_o2 = cat_vol_morph(Yppiscr>0.5,[dm 'o'], d*2 ); % rf / opt.interpV - 1 
          Yppi_c2 = cat_vol_morph(Yppiscr<0.5,[dm 'o'], d*2 ); 

          Yppi_od = cat_vol_morph(Yppiscr>0.5 & ~Yppi_o & ~cat_vol_morph(Yppiscr<0.5 & ~Yppi_c2,[dm 'd'],d),[dm 'd'], d ); % (rf / opt.interpV - 1)/3 
          Yppi_cd = cat_vol_morph(Yppiscr<0.5 & ~Yppi_c & ~cat_vol_morph(Yppiscr>0.5 & ~Yppi_o2,[dm 'd'],d),[dm 'd'], d ); 
          if ~debug,  clear Yppi_c Yppi_c2 Yppi_o Yppi_o2; end 

          Yppi_od = smooth3(Yppi_od)>0.5; 
          Yppi_cd = smooth3(Yppi_cd)>0.5; 
          
          %%
          if exist('Yppi_mx','var')
            Yppi_gyri  = Yppi_mn>0.1 & Yppi_od & ~Yppi_cd;
            Yppi_sulci = Yppi_mx<0.9 & Yppi_cd & ~Yppi_od;

            Yppi_gyri  = Yppi_gyri  | (Yppi_mx>0.9 & Yppi>0.3 & Yppi_mn>0.5); 
            Yppi_sulci = Yppi_sulci | (Yppi_mn<0.1 & Yppi<0.7 & Yppi_mx<0.5); 
            if ~debug, clear Yppi_mn Yppi_mx; end 
          else
            Yppi_gyri  = Yppi_od & ~Yppi_cd;
            Yppi_sulci = Yppi_cd & ~Yppi_od;
          end
          
          if ~debug, clear Yppi_od Yppi_cd Yppiscrc Yppiscmn Yppiscmx; end

          Yppiscr  = min( ~(Yppi_sulci & ~Yppi_gyri) , max( Yppisc , Yppi_gyri & ~Yppi_sulci)); 
       
          if ~debug, clear Yppi_gyri Yppi_sulci; end 
        end
        
        if ~debug, clear Yppigyri Yppislci; end
        [Yppiscr,resL] = cat_vol_resize(Yppiscr,'interp',VI,rf);
        %evalc(sprintf('clear CS; [Yppi05c,CS.faces,CS.vertices] = cat_vol_genus0(Yppiscr,0.5,1);')); % no_adjustment
        [Yppi05c,CS] = cat_vol_genus0opt(Yppiscr,th_initial,5 * (1-iscerebellum),debug);
        [Yvxdef,defect_number0] = spm_bwlabel( double(abs(Yppi05c - (Yppiscr>0.5))>0) ); clear Yppi05c;
        Yvxdef = cat_vol_resize(cat_vol_morph(Yvxdef,'d'),'deinterp',resL); % #### this is not ideal and need refinement ###  
        if ~debug, clear Yppiscr; end
        
        CS.vertices = CS.vertices * rf/opt.interpV;
      end
    end
    %
    EC0            = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
    vdefects       = cat_surf_fun('isocolors',CS,cat_vol_morph(Yvxdef,'d'),Smat.matlabIBB_mm)>0; clear Yvxdef;
    defect_size0   = sum(vdefects > 0) / length(vdefects) * 100; % percent
    defect_area0   = sum(vdefects > 0) / length(vdefects) .* ...
      sum(cat_surf_fun('area',CS)) / opt.interpV / 100; % cm2
    if opt.verb>1
      cat_io_cprintf('g5',sprintf('( SC/EC/DN/DS = %0.2f/',scale_cortex));
      cat_io_cprintf( color( rate( abs( EC0 - 2 ) , 0 ,100 * (1+9*iscerebellum) )) ,sprintf('%d/',EC0));
      cat_io_cprintf( color( rate( defect_number0 , 0 ,100 * (1+9*iscerebellum) )) ,sprintf('%d/',defect_number0));
      cat_io_cprintf( color( rate( defect_size0   , 1 , 10 * (1+9*iscerebellum) )) ,sprintf('%0.2f%%%%' ,defect_size0));
      cat_io_cprintf('g5',' )');
      fprintf(repmat(' ',1,max(0,14 - numel(sprintf('%d/%d/%0.2f%%%% )',EC0,defect_number0,defect_size0))))); 
    end
    
    % translate to mm coordinates
    CS = cat_surf_fun('smat',CS,Smat.matlabIBB_mm);   
    if opt.surf_measures > 1, CSraw0 = CS; end       % need this map later to create a common defect map
    if ~debug, clear mask_parahipp_smoothed; end
    
    
    % reduce resolution with higher resolution for cerebellum and fast option
    %  ##########
    %  * Both the SPM as well as the MATLAB function crashed my MATLAB
    %    multiple times (unreproducible and fatal).
    %    However, I have no idea why this happen and if it only on my system
    %    or how I could avoid or catch it because it is not just a simple error. 
    %    > This also happens if I only use double.
    %    > It also happens on the server. 
    %  RD201911
    %  * use the same mesh resolution for the cerebellum for acceptable processing times. 
    %  ##########
    if opt.reduce_mesh>2 
      CS.vertices = double(CS.vertices); CS.faces = double(CS.faces); 
      if opt.reduce_mesh == 3 || opt.reduce_mesh == 5
        CS = spm_mesh_reduce(CS, 81920 / (1 + (opt.vdist>2)) * (1 + 0*iscerebellum) );
      elseif  opt.reduce_mesh == 4 || opt.reduce_mesh == 6
        CS = reducepatch(CS, 81920 / (1 + (opt.vdist>2)) * (1 + 0*iscerebellum) );
      elseif  opt.reduce_mesh == 7
        CS = cat_surf_fun('reduce',CS, 81920 / (1 + (opt.vdist>2)) * (1 + 0*iscerebellum) ); 
      end  
      % remove bad faces 
      CS = correctReducePatch(CS);
    end
    saveSurf(CS,Praw);

    % remove unconnected meshes
    cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
     
    % refine super-large faces with adaption for cerebellum and fast option
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,3 / ( 1 + (opt.fast==1)) ); % only deformation for fast pipeline
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm); 
fprintf('\nRAW: V=%d, MN(CT)=%0.20f, SD(CT)=%0.20f\n',size(CS.vertices,1),mean(facevertexcdata(:)),std(facevertexcdata(:)));    
res.(opt.surf{si}).createCS_0_initfast = cat_surf_fun('evalCS',CS,cat_surf_fun('isocolors',CS,Yth1i,Smat.matlabIBB_mm),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);
   
    
    % Create a smooth surface for the topology correction. 
    % It don't has to be perfect because it will replaced completelly!
    for li = 2.^(0:2) 
      cmds = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                     'avg  %0.3f  %0.3f .2  .1  2  0 "0.5"  "0.5"  n 0  0  0 %d %g  0.0 0'], ...          
                      Vpp1.fname,Praw,Praw,-0.1/li,0.1/li,moveth(2/li,opt.surfaccuracy*2/li));
    end
    [ST, RS] = cat_system(cmds); cat_check_system_output(ST,RS,opt.verb-3);
    
    % load surf and map thickness
    CS = loadSurf(Praw);
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS,Smat.matlabIBB_mm); 
fprintf('IS: V=%d, MN(CT)=%0.20f, SD(CT)=%0.20f\n',size(CS.vertices,1),mean(facevertexcdata(:)),std(facevertexcdata(:)));    
res.(opt.surf{si}).createCS_0_initfast = cat_surf_fun('evalCS',CS,cat_surf_fun('isocolors',CS,Yth1i,Smat.matlabIBB_mm),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);
           
    if opt.fast
    %% Fast processing without topology correction and spherical registration
    %  --------------------------------------------------------------------
    %  The one and only fast option that is equal to the init surface but 
    %  with collision correction. For fast surface and thickness outputs 
    %  for visual analysis. 
    %  --------------------------------------------------------------------
   
        % correction for surface collision of the IS and OS
        if opt.collcorr > 1      
          % final correction of central surface in highly folded areas with high mean curvature
          cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
          cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0',Praw,Ppbt,Player4);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
          Sl4 = loadSurf(Player4);
          Yl4 = cat_surf_fun('isocolors',Ymf,Sl4.vertices,Smat.matlabi_mm); delete(Player4), clear Sl4; 

          stime = cat_io_cmd(sprintf('  Correction of surface collisions:'),'g5','',opt.verb-1,stime); 
          [CS,facevertexcdata] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdata,Ymfs,Yppi,Yl4,struct('optimize',opt.collcorr - 3,'verb',opt.verb>1,'mat',Smat.matlabIBB_mm));
          collcorrstr = 'collcorr'; 
        else
          collcorrstr = '';
        end

        % save final data and data structure
        saveSurf(CS,Pcentral);
        cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
        S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata);

        % create white and central surfaces
        cat_surf_fun('white',Pcentral);
        cat_surf_fun('pial',Pcentral);

        % evaluate and save results
        fprintf('%5.0fs',etime(clock,stime)); stime = []; 
        if opt.surf_measures > 3 % final result
          cat_surf_fun('saveico',CS,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm),Pcentral,...
            sprintf('createCS_1_initfast_%s_pbtres%0.2fmm_vdist%0.2fmm',collcorrstr,opt.interpV,opt.vdist),Ymfs,Smat.matlabIBB_mm); 
        end   
        res.(opt.surf{si}).createCS_1_initfast = cat_surf_fun('evalCS',CS,cat_surf_fun('isocolors',CS,Yth1i,Smat.matlabIBB_mm),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);
        res.(opt.surf{si}).createCS_final      = res.(opt.surf{si}).createCS_1_initfast;
       
          
        %  map WM and CSF width data (corrected by thickness)
        %  cat_surf_parameters and removed here (RD201909)
        if opt.surf_measures > 3 
            facevertexcdata2  = cat_surf_fun('isocolors',Ywd,CS.vertices,Smat.matlabi_mm); 
            facevertexcdata2c = max(eps,facevertexcdata2 - facevertexcdata/2); clear facevertexcdata2
            cat_io_FreeSurfer('write_surf_data',Pgwo,facevertexcdata2c); % gyrus width WM only
            facevertexcdata2c = correctWMdepth(CS,facevertexcdata2c,100,0.2);
            cat_io_FreeSurfer('write_surf_data',Pgww,facevertexcdata2c); % gyrus width WM only
            facevertexcdata3c = facevertexcdata2c + facevertexcdata; % );
            cat_io_FreeSurfer('write_surf_data',Pgw,facevertexcdata3c); clear facevertexcdata3c; % gyrus width (WM and GM)
            facevertexcdata4 = estimateWMdepthgradient(CS,facevertexcdata2c); clear facevertexcdata2c
            cat_io_FreeSurfer('write_surf_data',Pgwwg,facevertexcdata4); clear facevertexcdata4 % gyrus width WM only > gradient
          
            % smooth resampled values
            try %#ok<TRYNC>
              cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Pgwwg,3,Pgwwg);
              [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
            end

            facevertexcdata3 = cat_surf_fun('isocolors',Ycd,CS.vertices,Smat.matlabi_mm);
            facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
            cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3); clear facevertexcdata3

            setfield(S.(opt.surf{si}),'th2',nan(size(facevertexcdata)));  
            setfield(S.(opt.surf{si}),'th3',nan(size(facevertexcdata)));
            clear facevertexcdata; 
        end

        if ~debug
          if exist(Vpp1.fname ,'file'), delete(Vpp1.fname); end
          if exist(Vpp.fname,'file') && ~opt.outputpp.native, delete(Vpp.fname); end
        end

        % estimate Euler characteristics: EC = #vertices + #faces - #edges
        EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
        if any( strcmp( {'lh','rh'} , opt.surf{si} ))
          EC  = EC + abs(EC0 - 2) + 2;
        end

        % estimate FreeSurfer thickness measure Tfs using mean(Tnear1,Tnear2)
        if opt.thick_measure == 1
          stime = cat_io_cmd('  Tfs thickness estimation:','g5','',opt.verb,stime);
          cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

          % apply upper thickness limit
          facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick);  
          facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
          cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
        else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
          copyfile(Ppbt,Pthick);
        end

        fprintf('%5.0fs\n',etime(clock,stime)); 
        
        clear CS
        continue
    end
    %% evaluate and save results
    if isempty(stime), stime = clock; end
    fprintf('%5.0fs',etime(clock,stime)); stime = []; 
    if opt.surf_measures > 4 % just a substep
      cat_surf_fun('saveico',CS,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm),Pcentral,sprintf('createCS_1_init_pbtres%0.2fmm_vdist%0.2fmm',opt.interpV,opt.vdist),Ymfs,Smat.matlabIBB_mm); 
    else
      fprintf('\n'); 
    end
    res.(opt.surf{si}).createCS_init = cat_surf_fun('evalCS',CS,cat_surf_fun('isocolors',CS,Yth1i,Smat.matlabIBB_mm),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);
    
    
    
        
    
    %% Topology correction and surface refinement
    %  --------------------------------------------------------------------
    %  This topology correction creates a completely new surface based on  
    %  spherical hormonic functions resulting in a relative unbalanced
    %  local resolution (i.e., oversampled in the insula) that is corrected
    %  in the next block.  However, this also means the resolution of the
    %  input surface don't have to be super high (see above). 
    %  --------------------------------------------------------------------
    stime = cat_io_cmd('  Topology correction:','g5','',opt.verb,stime); 
    
    % spherical surface mapping 1 of the uncorrected surface for topology correction
    % We do not need so much smoothing as for the final surface but the
    % cerebellum needs maybe more due to servere topology defects.
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',Praw,Psphere0,...
      5 + round( sqrt( size(CS.faces,1) / 10000 * (1 + 3*iscerebellum) ) )); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % estimate size of topology defects 
    cmd = sprintf('CAT_MarkDefects "%s" "%s" "%s"',Praw,Psphere0,Pdefects0); 
    [ST, RS] = cat_system(cmd);
    sdefects       = cat_io_FreeSurfer('read_surf_data',Pdefects0); delete(Pdefects0);  
    defect_number0 = defect_number0 + max(sdefects); 
    defect_size0   = defect_size0   + sum(sdefects > 0) / length(sdefects) * 100; % percent
    defect_area0   = defect_area0   + sum(sdefects > 0) / length(sdefects) .* ...
      sum(cat_surf_fun('area',CS)) / opt.interpV / 100; % cm2
    % estimate Euler characteristics: EC = #vertices + #faces - #edges
    EC0            = (EC0-2) + ( size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1) - 2) + 2;
    if any( strcmp( {'lh','rh'} , opt.surf{si} ))
      EC             = EC + abs(EC0 - 2) + 2; % -2 is the correction for the sphere
      defect_size    = defect_size   + defect_size0;
      defect_area    = defect_area   + defect_area0;
      defect_number  = defect_number + defect_number0;
    end    
    
    % topology correction and surface refinement
    % Higher -n will result in larger but still unbalanced meshes and the 
    % refine_lenght parameter is more important to obtain nice meshes.
    if opt.verb>3, fprintf('\n'); end
    cmd = sprintf('CAT_FixTopology -lim %d -bw %d -n %d -refine_length %g "%s" "%s" "%s"',...
    ...  512,1024, 81920, opt.vdist ,Praw,Psphere0,Pcentral);
      256 / (1 + iscerebellum),1024 / (1 + iscerebellum), 81920 * (1 + 0*3*iscerebellum), opt.vdist ,Praw,Psphere0,Pcentral); % avoid to long processing in cerebellum
    ...  512 * (1 + iscerebellum),1024 * (1 + iscerebellum), 81920 * (1 + 0*3*iscerebellum), opt.vdist ,Praw,Psphere0,Pcentral);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
CS = loadSurf(Pcentral); 
facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm); 
fprintf('TC: V=%d, MN(CT)=%0.20f, SD(CT)=%0.20f\n',size(CS.vertices,1),mean(facevertexcdata(:)),std(facevertexcdata(:)));    
res.(opt.surf{si}).createCS_0_initfast = cat_surf_fun('evalCS',CS,cat_surf_fun('isocolors',CS,Yth1i,Smat.matlabIBB_mm),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);


  
  
    %% Optimize topology corrected mesh: 
    %  --------------------------------------------------------------------
    %  Oversampling can lead to problems in the normal transformation to 
    %  obtain the inner and outer surface, especially in the Insula/Amygdala.
    %  The idea is to reduce the sampling of the mesh and then to refine 
    %  the mesh again. 
    %  However it is highly essential that most regions where not corrected
    %  and a specific masking of the Insula (relative small triangles in a
    %  specific area on one/all of the main cortical surfaces or flipping) 
    %  is possible useful.
    %  ... seams that this is working and it takes only a few seconds!
    %  --------------------------------------------------------------------
    % refinement - important for sulci .. here we need a lot of details with a similar resolution as the Insula 
    if opt.reduce_mesh > 4 % superinterpolation 
      stime = cat_io_cmd('  Surface optimization and refinement:','g5','',opt.verb,stime); 
      meshres = 0.8; % this will create a super resolution that has to be reduced or will result in long processing times 
    else
      stime = cat_io_cmd('  Surface refinement:','g5','',opt.verb,stime); 
      meshres = opt.vdist;
    end
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Pcentral,Pcentral, meshres ); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % surface refinement (this time even before reduction)
    for li = 2.^(0:1)
      cmds = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                     'avg  -0.1  0.1 .2  .1  %d  0 "0.5"  "0.5"  n 0  0  0 %d  %0.2f  0.0 0'], ...    
                      Vpp1.fname,Pcentral,Pcentral,1/li,moveth(0.4/li,opt.surfaccuracy * 4/li));
      [ST, RS] = cat_system(cmds); cat_check_system_output(ST,RS,opt.verb-3); % hier ging was schief
    end
    
    % reduce - as far as the Insula/Amygdala is not so heavily folded compared to sulci this region is first reduced 
    CS = loadSurf(Pcentral); 
    if opt.reduce_mesh > 4
      rfaces =  min( max( 81920 , opt.reduceCS/2 ) , 81920 * 4 );
      CS.vertices = double(CS.vertices); CS.faces = double(CS.faces); 
      if opt.reduce_mesh == 5
        CS = spm_mesh_reduce(struct('vertices',CS.vertices,'faces',CS.faces),rfaces); 
      elseif  opt.reduce_mesh == 6
        CS = reducepatch(struct('vertices',CS.vertices,'faces',CS.faces),rfaces); 
      elseif  opt.reduce_mesh == 7
        CS = cat_surf_fun('reduce',CS,rfaces); 
      end
      % remove bad faces 
      CS = correctReducePatch(CS);
    end
    saveSurf(CS,Pcentral); 

     % remove unconnected meshes
    cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Pcentral,Pcentral); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
   
    % refinement - guaranty our default resolution
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Pcentral,Pcentral,opt.vdist);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    
    % surface defomation for relaxation after reduction and refinement
    for li = 2.^(1:2)
      cmds = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                     'avg  -0.1  0.1 .2  .1  %d  0 "0.5"  "0.5"  n 0  0  0 %d  %0.2f  0.0 0'], ...    
                      Vpp1.fname,Pcentral,Pcentral,1/li,moveth(0.4/li,opt.surfaccuracy*4/li));
      [ST, RS] = cat_system(cmds); cat_check_system_output(ST,RS,opt.verb-3);
    end 
    % read final surface and map thickness data
    CS = loadSurf(Pcentral);
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
fprintf('SR1: V=%d, MN(CT)=%0.20f, SD(CT)=%0.20f\n',size(CS.vertices,1),mean(facevertexcdata(:)),std(facevertexcdata(:)));    
res.(opt.surf{si}).createCS_0_initfast = cat_surf_fun('evalCS',CS,cat_surf_fun('isocolors',CS,Yth1i,Smat.matlabIBB_mm),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);

    
    
    
    
    %% final correction of central surface in highly folded areas 
    %  with high mean curvature with weight of 0.7 and further refinement
    %  of the mesh and its vertices based on the position map 
    cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.7 "%s" "%s" "%s" 0',Pcentral,Ppbt,Pcentral);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % we need some refinement because some vertices are too large to be deformed with high accuracy
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 1',Pcentral,Pcentral,opt.vdist); % adaption for cerebellum
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % surface refinement by surface deformation based on the PP map
    for li = 2.^(0:1)
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'avg -0.1 0.1 .2 .1 %d 0 "0.5" "0.5" n 0 0 0 %d %0.2f 0.0 %d'], ...
                     Vpp1.fname,Pcentral,Pcentral,1/li,moveth(0.2,opt.surfaccuracy * 2/li),any(opt.collcorr==[1,5])); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    end
    
    % need some more refinement because some vertices are distorted after CAT_DeformSurf
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 1',Pcentral,Pcentral,opt.vdist); % adaption for cerebellum
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % final surface refinement
    for li = 2.^(0:1)
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'avg -0.05 0.05 .5 .1 %d 0 "0.5" "0.5" n 0 0 0 %d %0.2f 0.0 %d'], ...
                     Vpp1.fname,Pcentral,Pcentral,1/li,moveth(0.1/li,opt.surfaccuracy * 2/li),any(opt.collcorr==[1,5])); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    end

    % read final surface and map thickness data
    CS = loadSurf(Pcentral);
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
   
    % evaluate and save results
    fprintf('%5.0fs',etime(clock,stime)); stime = []; 
    if opt.surf_measures > 4 % just a substep
      cat_surf_fun('saveico',CS,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm),Pcentral,sprintf('createCS_2_refined_pbtres%0.2fmm_vdist%0.2fmm',opt.interpV,opt.vdist),Ymfs,Smat.matlabIBB_mm); 
    else
      fprintf('\n'); 
    end
    res.(opt.surf{si}).createCS_2_refine = cat_surf_fun('evalCS' ,CS,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm),...
      Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2,opt.collcorr==0 && cat_get_defaults('extopts.expertgui')>1);
    res.(opt.surf{si}).createCS_final    = res.(opt.surf{si}).createCS_2_refine; 
fprintf('SR2: V=%d, SD(CT)=%0.20f\n',size(CS.vertices,1),std(facevertexcdata(:)));    
    
    
    
    %% Collision correction by Delaunay triangularization
    if opt.collcorr > 1  
      stime = cat_io_cmd('  Correction of surface collisions:','g5','',opt.verb,stime); 

      % final correction of central surface in highly folded areas with high mean curvature
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0',Pcentral,Ppbt,Player4);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
      Sl4 = loadSurf(Player4);
      Yl4 = cat_surf_fun('isocolors',Ymf,Sl4,Smat.matlabi_mm); delete(Player4), clear Sl4; 
 
      
      %% call collision correction
      if debug, if exist('CSO','var'), CS = CSO; facevertexcdata = facevertexcdatao; else, CSO = CS; facevertexcdatao = facevertexcdata; end; stime2 = clock; else, stime2 = [];  end
      if opt.collcorr == 2
        [CS,facevertexcdata] = cat_surf_fun('collisionCorrectionRY',CS,facevertexcdata,Ymfs,struct('Pcs',Pcentral,'verb',opt.verb>1,'mat',Smat.matlabIBB_mm,'accuracy',1/2^6));
      elseif opt.collcorr == 3 || opt.collcorr == 4
        [CS,facevertexcdata] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdata,Ymfs,Yppi,Yl4,struct('optimize',opt.collcorr - 3,'verb',opt.verb>1,'mat',Smat.matlabIBB_mm));
      elseif opt.collcorr == 5
        [CS,facevertexcdata] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdata,Ymfs,Yppi,Yl4,struct('optimize',4,'verb',opt.verb>1,'mat',Smat.matlabIBB_mm));
        fprintf('\b\b');
        [CS,facevertexcdata] = cat_surf_fun('collisionCorrectionRY',CS,facevertexcdata,Ymfs,struct('Pcs',Pcentral,'verb',opt.verb>1,'mat',Smat.matlabIBB_mm,'accuracy',1/2^3));
      end
      if ~debug, clear Yl4; end
      
      % saveSurf(CS,Pcentral); 
      cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
      
      % evaluate and save results
      if opt.verb > 2, cat_io_cmd(' ','g5','',opt.verb);  end
      fprintf('%5.0fs',etime(clock,min([stime2;stime],[],1))); if ~debug, stime = []; end
      % final result ... test for self-intersections only in developer mode? 
      if opt.surf_measures > 2  ||  opt.verb > 2
        cat_surf_fun('saveico',CS,facevertexcdata,Pcentral,sprintf('createCS_3_collcorr_%0.2fmm_vdist%0.2fmm',opt.interpV,opt.vdist),Ymfs,Smat.matlabIBB_mm); 
      else
        fprintf('\n');
      end
      res.(opt.surf{si}).createCS_3_collcorr = cat_surf_fun('evalCS' ,CS,facevertexcdata,Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-1,cat_get_defaults('extopts.expertgui')>1);
      res.(opt.surf{si}).createCS_final      = res.(opt.surf{si}).createCS_3_collcorr; 
      
      %%
      clear Yl4 CSO; 
      
      
    elseif opt.collcorr == 1
    %% final correction of cortical thickness using pial and WM surface
    %  This does not work right now. Although, there are maybe no self
    %  intersections, the deformation to pial and white are still
    %  incomplete, resulting in high thickness underestimation! 
    
      % estimation of pial surface
      th2 = 0.1; % GM/CSF border
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp1.fname,Pcentral,Ppial,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % if deformation stopped earlier call CAT_Central2Pial to get closer to the pial surface
      if ~isempty(strfind(RS,'Stopped after'))
        stime = cat_io_cmd('  Estimation of pial surface','g5','',opt.verb,stime);
        cmd = sprintf('CAT_Central2Pial -check_intersect "%s" "%s" "%s" 0.3',Pcentral,Ppbt,Ppial);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp1.fname,Ppial,Ppial,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end

      stime = cat_io_cmd('  Correction of pial surface in highly folded areas','g5','',opt.verb,stime);
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.5 "%s" "%s" "%s" 0',Ppial,Ppbt,Ppial);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % estimation of white matter surface
      th2 = 0.9; % GM/WM border
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp1.fname,Pcentral,Pwhite,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % if deformation stopped earlier call CAT_Central2Pial to get closer to the white matter surface
      if ~isempty(strfind(RS,'Stopped after')) 
        stime = cat_io_cmd('  Estimation of white matter surface','g5','',opt.verb,stime);
        cmd = sprintf('CAT_Central2Pial -check_intersect "%s" "%s" "%s" -0.3', ...
                       Pcentral,Ppbt,Pwhite);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp1.fname,Pwhite,Pwhite,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end

      stime = cat_io_cmd('  Correction of white matter surface in highly folded areas','g5','',opt.verb,stime);
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.3 "%s" "%s" "%s" 0',Pwhite,Ppbt,Pwhite);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % update central surface as average between white and pial surface
      cmd = sprintf('CAT_AverageSurfaces -avg "%s" "%s" "%s"',Pcentral,Pwhite,Ppial);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % ####
      % More testing necessary to also correct thickness.
      % However, if Tnear is bad than also the boundary surfaces are
      % incorrect, ie. this pipeline does not work yet.
      % ####
      cmd = sprintf('CAT_Hausdorff  "%s" "%s" "%s"',Pwhite,Ppial,Ppbt);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      saveiconame = 'deform';
      facevertexcdata = cat_io_FreeSurfer('read_surf_data',Ppbt);
     
      % evaluate and save results
      CS = loadSurf(Pcentral); 
      fprintf('%5.0fs',etime(clock,stime)); stime = []; 
      if opt.surf_measures > 3 % final result
        cat_surf_fun('saveico',CS,facevertexcdata,Pcentral,...
          sprintf('debug_3_extract_pial_white_pbtres_%s_pbtres%0.2fmm_vdist%0.2fmm',saveiconame,opt.interpV,opt.vdist),Ym);  
      else
        fprintf('\n'); 
      end
      res.(opt.surf{si}).createCS_3_extract_pial_white = cat_surf_fun('evalCS' ,CS,facevertexcdata,Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,opt.verb-2);
      res.(opt.surf{si}).createCS_final                = res.(opt.surf{si}).createCS_3_extract_pial_white; 
    end

    
    % Further error measures: 
    % - Map local noise pattern?  
    %   >> Detection of movement artifacts?
    % - Map local bias field?     
    %   >> This finally codes the heat noise of the original image >> input by cat_main 
    %   >> Relevant for local interpretation of results, e.g. low signal (high noise) ares are less reliable. 
    % ... both aspects are more general QA topics and not so relevant here
    %
    % - Map collision error?      
    %   >> output of collision detection? - not independent
    %   >> only relevant for me? 
    % - Map distance error? 
    %   >> To what? RAW surface? This would include topology defect errors.
    %      Possibly mask by defect map? And then? 
    %   >> estimate distance of each GM voxels to its closes surface point
    %      . describes the number of voxels out of thickness range (missing GM voxels - but also artifacts) 
    %      estimate distance of each ~GM brain voxels within the local thickness
    %      . describes the number of voxels the are in the cortical surface ribbon but should not 
    %      ... the biggest problem are here the artifacts ...
    %      ... use defect map to avoid counting in some region 

    
    
    % map defects to the final surface 
    if opt.surf_measures > 1
      %%
      CSraw   = loadSurf(Praw); 
      
      Yvdef   = cat_surf_fun('surf2vol',struct('vertices',CSraw0.vertices,'faces',CSraw.faces),Ymfs>1.1,(vdefects)>0,'val',struct('mat',Smat.matlabIBB_mm,'verb',0));
      Ysdef   = cat_surf_fun('surf2vol',struct('vertices',CSraw.vertices ,'faces',CSraw.faces),Ymfs>1.1,(sdefects)>0,'val',struct('mat',Smat.matlabIBB_mm,'verb',0));
      
      Yvdef(isnan(Yvdef)) = 0;
      Ysdef(isnan(Ysdef)) = 0;

      %%
      Ydef = Yvdef + Ysdef; clear Ysdef; 
      [Ydef,ndef] = spm_bwlabel(double(Ydef)); 
      Ydefn = Ydef; 
      for i=1:ndef
        Ydefn(Ydef==i) = sum(Ydef(:)==i) / sum(Ydef(:)>0); 
      end
        
      %%
      defects = cat_surf_fun('isocolors',Ydefn,CS.vertices,Smat.matlabIBB_mm,'nearest'); 
      cat_io_FreeSurfer('write_surf_data',Pdefects,defects);
      clear defects Ydef Ydefn;
    end
    clear Yvxdef;
    

    
    
    % Test without surface registration - just a shortcut for manual tests! 
    if 0
      cat_io_cmd('  ','g5','',opt.verb,stime);  
      S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'vmat',vmat,'vmati',vmati,'mati',mati,'th1',facevertexcdata); 
      if si==numel(opt.surf) && si == 1
        cat_io_cmd('  ','g5','',opt.verb,cstime);
        sprintf('%5ds\n',round(etime(clock,cstime)));
      end
      continue
    end
    
    
    
    
    %% spherical surface mapping and registration of the final corrected surface
    %  use more iterations for higher surfaces (sqrt due to surface area)
    %  no extra rule for the cerebellum here because it is has the correct topology. 
    stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',Pcentral,Psphere,...
      5 + round( sqrt( size(CS.faces,1) / 10000 ) + 1 )); % 300k with value 10
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    
    % spherical registration to fsaverage template
    stime = cat_io_cmd('  Spherical registration','g5','',opt.verb,stime);
    if opt.vdist>2 % low quality 
      cmd = sprintf(['CAT_WarpSurf -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s" ' ...
        '-size 256 128 -loop 1 -steps 1 -runs 2'],...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    else
      cmd = sprintf('CAT_WarpSurf -steps 2 -avg -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"', ...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    
    % setting of thickness values to zero for masked area (by the inverse transformation)
    % does not work properly for all data and is moreover part of cat_surf_resample ...

    % display some evaluation 
    if opt.verb>1
      fprintf('%5.0fs\n',etime(clock,stime)); 
      fprintf('    Euler char. / def. number / def. size: ');
      cat_io_cprintf( color( rate(  EC0 - 2        , 0 , 100 * (1+4*iscerebellum)) ) , sprintf('%0.0f / '   , EC0 ) );
      cat_io_cprintf( color( rate(  defect_number0 , 0 , 100 * (1+4*iscerebellum)) ) , sprintf('%0.0f / '   , defect_number0 ) );
      cat_io_cprintf( color( rate(  defect_size0   , 0 ,  10 * (1+4*iscerebellum)) ) , sprintf('%0.2f%%%% ' , defect_size0 ) );
      fprintf('\n');
    end
    
    % evaluate and save results
    if opt.surf_measures > 4 
    % This part is not highly relevant for the individual surface reconstruction 
    % but it can help to test and optimize the spatial registration. 
      
      % filenames for resmapling
      Presamp   = fullfile(pp,surffolder,sprintf('%s.tmp.resampled.%s'    ,opt.surf{si},ff));  
      Ppbtr     = fullfile(pp,surffolder,sprintf('%s.pbt.resampled.%s'    ,opt.surf{si},ff));  
      Ppbtr_gii = [Ppbtr '.gii'];
      
      % resample values using warped sphere 
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavgsph,Presamp,Ppbt,Ppbtr);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      
      if 0 
        % resample surface using warped sphere with better surface quality (using Spherical harmonics)
        % ###
        % deactivated because the resampling of the surface alone leads to displacements of the textures (RD20190927)!
        % ###
        cmd = sprintf('CAT_ResampleSphericalSurfSPH -n 327680 "%s" "%s" "%s"',Pcentral,Pspherereg,Presamp);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

        % resample surface according to freesurfer sphere
        cmd = sprintf('CAT_ResampleSurf "%s" NULL "%s" "%s"',Presamp,Pfsavgsph,Presamp);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3); 
      end
      
      % add values to resampled surf and save as gifti
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Ppbtr,Ppbtr_gii); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3); 
      if exist(Ppbtr,'file'), delete(Ppbtr); end

      %% remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
      [pp2,ff2,ex2] = spm_fileparts(Ppbtr_gii); 
      g = gifti(Ppbtr_gii);
      g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
      save(g, Ppbtr_gii, 'Base64Binary');
      
      %% intensity based evaluation
      CSr = loadSurf(Ppbtr_gii); 
      CSr = struct('vertices',CSr.vertices,'faces',CSr.faces,'cdata',CSr.cdata);
      cat_surf_fun('saveico',CSr,CSr.cdata,Pcentralr,sprintf('createCS_4_resampled_pbtres%0.2fmm_vdist%0.2fmm',opt.interpV,opt.vdist),Ymfs,Smat.matlabIBB_mm); 
      res.(opt.surf{si}).createCS_resampled = cat_surf_fun('evalCS',CSr,CSr.cdata,Ymfs,Yppi,Pcentralr,Smat.matlabIBB_mm);
      clear CSr
    end
    if ~isfield( res.(opt.surf{si}),'createCS_final')
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',loadSurf(Pcentral),cat_io_FreeSurfer('read_surf_data',Ppbt),Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm);
    else 
      fprintf('\n'); 
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
        if isfield(res.(opt.surf{si}),'createCS_resampled') 
          if ~isfield(res,'createCS_resampled') || ~isfield(res.createCS_resampled,FNres{fnr}) 
            res.resampled.(FNres{fnr}) = res.(opt.surf{si}).createCS_resampled.(FNres{fnr}) / numel(opt.surf);
          else
            res.resampled.(FNres{fnr}) = res.resampled.(FNres{fnr}) + res.(opt.surf{si}).createCS_resampled.(FNres{fnr}) / numel(opt.surf);
          end
        end
      end
    end
    
    
    % create white and central surfaces
    cat_surf_fun('white',Pcentral);
    cat_surf_fun('pial',Pcentral);
    
    
    % write myelination map (Ypp intensity of layer 4)  
    if opt.surf_measures > 1 
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0',Pcentral,Ppbt,Player4);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
      L4  = gifti(Player4);
      L4v = cat_surf_fun('isocolors',Ymf,L4,Smat.matlabi_mm); clear L4
      cat_io_FreeSurfer('write_surf_data',PintL4,L4v); clear L4v
    end
    
    
    % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
    if opt.thick_measure == 1
      % not ready yet
      if 0
%      if opt.extract_pial_white && ~opt.fast % use white and pial surfaces
        cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',Pwhite,Ppial,Pthick);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      else % use central surface and thickness
        cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end
      
      % apply upper thickness limit
      facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick);  
      facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
      cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
    else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
      copyfile(Ppbt,Pthick);
    end
    
    
    
    
    %% WM and CSF thickness
    %  Will hopefully be improved in future and may become part of 
    %  cat_surf_parameters and removed here (RD201909)
    if opt.surf_measures > 3 
      % map WM and CSF width data (corrected by thickness)
      facevertexcdata2  = cat_surf_fun('isocolors',Ywd,CS.vertices,Smat.matlabi_mm); 
      facevertexcdata2c = max(eps,facevertexcdata2 - facevertexcdata/2);
      cat_io_FreeSurfer('write_surf_data',Pgwo,facevertexcdata2c); % gyrus width WM only
      facevertexcdata2c = correctWMdepth(CS,facevertexcdata2c,100,0.2);
      cat_io_FreeSurfer('write_surf_data',Pgww,facevertexcdata2c); % gyrus width WM only
      facevertexcdata3c = facevertexcdata2c + facevertexcdata; % );
      cat_io_FreeSurfer('write_surf_data',Pgw,facevertexcdata3c); clear facevertexcdata3c; % gyrus width (WM and GM)
      facevertexcdata4 = estimateWMdepthgradient(CS,facevertexcdata2c); clear facevertexcdata2c; 
      cat_io_FreeSurfer('write_surf_data',Pgwwg,facevertexcdata4); clear facevertexcdata4; % gyrus width WM only > gradient
      
      % smooth resampled values
      try %#ok<TRYNC>
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Pgwwg,3,Pgwwg);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end
      
      facevertexcdata3 = cat_surf_fun('isocolors',Ycd,CS.vertices,Smat.matlabi_mm); 
      facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
      cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3);
    end
    
    
    % create output structure
    S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata); clear facevertexcdata; 
    if opt.surf_measures > 3
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th2',facevertexcdata2); clear facevertexcdata2; 
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th3',facevertexcdata3); clear facevertexcdata3; 
    end
    clear Yth1i

    
    % we have to delete the original faces, because they have a different number of vertices after
    % CAT_FixTopology!
    if exist(Praw      ,'file'), delete(Praw); end
    if exist(Psphere0  ,'file'), delete(Psphere0); end
    if exist(Vpp1.fname ,'file'), delete(Vpp1.fname); end
    if exist(Vpp.fname,'file') && ~opt.outputpp.native, delete(Vpp.fname); end
    if opt.verb > 2 && exist(Pdefects0,'file'), delete(Pdefects0); end
    if opt.extract_pial_white && exist(Vyp0s.fname,'file'), delete(Vyp0s.fname); end
    clear CS

    % processing time per side for manual tests
    if si == numel(opt.surf) && si == 1
      %%
      fprintf('\n');
      cat_io_cmd('  ','g5','',opt.verb);
      fprintf('%5ds',round(etime(clock,cstime)));
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
    %   Especially the outer surface will describe the sulcal blurring in children. 
    %   But the mixing of surface quality and anatomical features is problematic. 
    % - The position value describes how good the transformation of the PBT map into a surface worked. 
    %   Also the position values depend on age. Children have worse pial values due to sulcal blurring but
    %   the white surface is may effected by aging, e.g., by WMHs.
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
      cat_io_cprintf( color( rate(  defect_size   , 0 , 10  * (1+9*iscerebellum)) ) , sprintf('%0.2f%%%% ' , defect_size ) );
      fprintf('\n');
    else
      fprintf('  Average thickness:                     %0.4f %s %0.4f mm\n' , mean(mnth), native2unicode(177, 'latin1'), mean(sdth));
      %fprintf('  Surface intensity / position RMSE:     %0.4f & %0.4f\n'      , mean(sdRMSE_Ym) ,mean(sdRMSE_Ypp) );
      fprintf('  Euler characteristic / defect size:    %0d / %0.2f%%%%\n'     , EC, defect_size);
    end
    
    for si=1:numel(Psurf)
      fprintf('  Display thickness: %s\n',spm_file(Psurf(si).Pthick,'link','cat_surf_display(''%s'')'));
    end
    
    
    %% surfaces in spm_orthview
    if exist(Pm,'file'), Po = Pm; else, Po = V0.fname; end
    
    Porthfiles = '{'; Porthcolor = '{'; Porthnames = '{';
    for si=1:numel(Psurf)
      Porthfiles = [ Porthfiles , sprintf('''%s'',''%s'',',Psurf(si).Ppial, Psurf(si).Pwhite )]; 
      Porthcolor = [ Porthcolor , '''-g'',''-r'',' ]; 
      Porthnames = [ Porthnames , '''white'',''pial'',' ];
    end
    Porthfiles = [ Porthfiles(1:end-1) '}']; 
    Porthcolor = [ Porthcolor(1:end-1) '}']; 
    Porthnames = [ Porthnames(1:end-1) '}']; 
      
    fprintf('    Show in orthview: %s\n',spm_file(Po ,'link',...
      sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Po,Porthcolor,Porthnames))) ;

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

function saveSurf(CS,P)
  save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),P,'Base64Binary');
end

function CS1 = loadSurf(P)
  CS = gifti(P);
  CS1.vertices = CS.vertices; CS1.faces = CS.faces; 
  if isfield(CS,'cdata'), CS1.cdata = CS.cdata; end
end

function CS = correctReducePatch(CS)
  % remove bad faces 
  badv = find( sum( spm_mesh_neighbours(CS)>0,2) == 2); 
  badf = []; for fi=1:numel(badv); [badfi,badfj] = find( CS.faces == badv(fi) );  badf = [badf;badfi]; end
  CS.faces(badf,:) = [];  
  for fi=numel(badv):-1:1; CS.faces(CS.faces > badv(fi)) = CS.faces(CS.faces > badv(fi)) - 1; end
  CS.vertices(badv,:) = []; 
end

%=======================================================================
function [cdata,i] = correctWMdepth(CS,cdata,iter,lengthfactor)
% ______________________________________________________________________
% Correct deep WM depth values that does not fit to the local thickness 
% of the local gyri.
% 
% length factor should be between 0.2 and 0.4
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
  while i<iter && pc~=oc 
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
function cdata = estimateWMdepthgradient(CS,cdata)
% ______________________________________________________________________
% Estimates the maximum local gradient of a surface. 
% Major use is the WM depth that grows with increasing sulcal depth. 
% It measures the amount of WM behind the cortex, but more relevant is
% the amount of WM fibers that this region will add to the WM depth. 
% The width of the street next to a house gives not the connectivity of
% this house, but the width of the entrance does!
% This measure can be improved by further information of sulcal depth.
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
              
