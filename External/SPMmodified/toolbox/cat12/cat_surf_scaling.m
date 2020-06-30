function rn = cat_surf_scaling(job)
%cat_surf_scaling: scaling to normalize surfaces by spherical properties
% 
% rn = cat_surf_scaling(job)
% 
% rn      .. linear scaling factor
% job     .. SPM job structure
%  .file  .. input file
%  .norm  .. scaling option
%  .fname .. outpout file name
%

  def.file  = {};
  def.norm  = 31;
  def.fname = {};
  job = cat_io_checkinopt(job,def);

  [pp,ff] = spm_fileparts(job.file); 
  
  S = gifti( job.file );

  if job.norm == 12 % affine 
    % RD20200211 - get XML information for affine normalization? 
    if strcmp(pp(end-3:end),'surf')
      reportdir = [pp(1:end-4) strrep(pp(end-3:end),'surf','report')]; 
    else
      reportdir = pp;
    end
    Pxml = fullfile( reportdir , cat_io_strrep( ['cat_' ff '.xml' ],{'lh.central.','rh.central.','cb.central.'},'') );

    X    = cat_io_xml( Pxml ); 
    mati = spm_imatrix( X.parameter.spm.Affine ) .* [0 0 0 0 0 0 1 1 1 0 0 0] - [ min( S.vertices )*1.1  0 0 0  0 0 0  0 0 0] ;
    mat  = spm_matrix(mati);

    vertices = mat * ([ S.vertices , ones( size(S.vertices,1) , 1) ])';
    vertices = vertices(1:3,:)';
    
  elseif job.norm == 1 || job.norm == 11
    % normalization as distance of all vertices to the center of mass (COM)
    % or to the vertices of the hull 
    if job.norm == 1
      COM = mean(S.vertices); 
    else % hull
      hf  = convhulln(double(S.vertices)); 
      COM = mean(S.vertices(unique( hf) ,:)); 
    end
    DCOM  = sum((S.vertices - repmat( COM , size(S.vertices,1), 1)).^2,2).^0.5; clear COM; 
    r     = mean(DCOM); clear DCOM; 
  
  elseif job.norm == 2 || job.norm == 21
    % normalization by the surface area
    if job.norm == 2
      SA  = cat_surf_fun('area',S); 
    else
      hf  = convhulln(double(S.vertices));
      SH  = struct('vertices',S.vertices,'faces',hf); 
      SA  = cat_surf_fun('area',SH); 
    end
    r  = sqrt( sum(SA) / (4 * pi) ); clear SA, % A_sphere = 4 * pi * r^2
  
  elseif job.norm == 31
    % normalization by hull volume 
    [hf,vol]  = convhulln(double(S.vertices)); clear hf;  %#ok<ASGLU>
    r  = ( vol / (4/3 * pi) )^(1/3); clear vol; % V_sphere = 4/3 * pi * r^3
    
  elseif job.norm == 32
    % normalization by volumebased TIV
    % RD20200211 - get XML information for affine normalization? 
    if strcmp(pp(end-3:end),'surf')
      reportdir = [pp(1:end-4) strrep(pp(end-3:end),'surf','report')]; 
    else
      reportdir = pp;
    end
    Pxml = fullfile( reportdir , cat_io_strrep( ['cat_' ff '.xml' ],{'lh.central.','rh.central.','cb.central.'},'') );

    X    = cat_io_xml( Pxml ); 
    r    = X;
    rn   = 60/r;
  end
  
  if job.norm == 12
  else
    rn = 60/abs(r);
    vertices = S.vertices * rn;          
  end
  
  % write data
  if ~isempty(job.fname)
    if exist([job.fname(1:end-3),'dat'],'file'), delete([job.fname(1:end-3),'dat']); end
    save( gifti(struct('faces',S.faces,'vertices',vertices)),job.fname,'Base64Binary');
  end
  
  
  % just a block of warnings for extrem normalization factors
  if     rn<0.25 || rn>2.00
    cat_io_cprintf('err' ,sprintf(['cat_surf_scale:smallSurface:Warning the normalisation factor is quite ' ... 
                                   'low (%0.2f), check for possible sampling artifcats in Toro''s GI.\n'],rn)); 
  elseif rn<0.50 || rn>1.50
    cat_io_cprintf('warn',sprintf(['cat_surf_scale:largeSurface:Warning the normalisation factor is quite ' ...
                                   'high (%0.2f), check for possible (ocipital) boundary artifacts in Toro''s GI.\n'],rn)); 
  end
end