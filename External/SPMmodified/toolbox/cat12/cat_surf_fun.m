function varargout = cat_surf_fun(action,S,varargin)
%cat_surf_fun Set of functions to modify and evaluate surfaces. 
%  Call without parameters (only action) will print the help of this action.
% 
%  varargout = cat_surf_fun(action,S,varargin)
% 
%  varargout, varargin .. variable input depending on the called action
%  action              .. string to call a subfuction 
%  S                   .. surface structure with vertices and faces 
%
%  Actions: 
%  * Distance estimation:
%    <a href="matlab:help cat_surf_fun>cat_surf_dist;">dist</a>    Estimate edge length of faces in S.
%              D = cat_surf_fun('dist',S);    
%    <a href="matlab:help cat_surf_fun>cat_surf_vdist;">vdist</a>   Estimate distance of each voxel to S.
%              D = cat_surf_fun('vdist',S);    
%    <a href="matlab:help cat_surf_fun>cat_surf_thickness;">tfs</a>     Estimate FreeSurfer thickness.
%              T = cat_surf_fun('Tfs',S,varargin{1}); 
%    <a href="matlab:help cat_surf_fun>cat_surf_thickness;">tmin</a>    Estimate minimum thickness.
%              T = cat_surf_fun('Tmin',S,varargin{1}); 
%    <a href="matlab:help cat_surf_fun>cat_surf_thickness;">tmax</a>    Estimate maximum thickness.
%              T = cat_surf_fun('Tmax',S,varargin{1}); 
%
%  * Data mapping:
%    <a href="matlab:help cat_surf_fun>cat_surf;">?</a>        Map face data to vertices. 
%                     V = cat_surf_fun('?',S,F);   
%    <a href="matlab:help cat_surf_fun>cat_surf_cdatamapping;">cdatamapping</a>   Map texture cdata from S2 to S1.
%                     C = cat_surf_fun('cdatamapping',S1,S2,cdata[,opt]); 
%    <a href="matlab:help cat_surf_fun>cat_surf_createEdgemap;">createEdgemap</a>  Creates mapping between two surfaces.  
%                     edgemap = cat_surf_fun('createEdgemap',S1,S2);
%    <a href="matlab:help cat_surf_fun>cat_surf_createEdgemap;">useEdgemap</a>     Apply mapping between two surfaces.  
%                     cdata2 = cat_surf_fun('useEdgemap',cdata,edgemap);
%
%  * Surface (data) rendering:
%    <a href="matlab:help cat_surf_fun>cat_surf_surf2vol;">surf2vol</a> Render surface (data) to a volume. 
%              V = cat_surf_fun('surf2vol',S,varargin{1}); 
%
%  * Surface modification:
%    <a href="matlab:help cat_surf_fun>cat_surf_hull;">hull</a>      Estimate hull surface. 
%              HS = cat_surf_fun('hull',S);      
%    <a href="matlab:help cat_surf_fun>cat_surf_core;">core</a>      Estimate core  surface.  
%              HS  = cat_surf_fun('core',S);    
%    <a href="matlab:help cat_surf_fun>cat_surf_createinneroutersurface;">inner</a>     Estimate inner surface.  
%              IS  = cat_surf_fun('inner',S,T);  
%    <a href="matlab:help cat_surf_fun>cat_surf_createinneroutersurface;">outer</a>     Estimate outer surface.  
%              OS  = cat_surf_fun('outer',S,T);  
%    <a href="matlab:help cat_surf_fun>cat_surf_saveICO;">saveICO</a>   Save multiple output surfaces & measures.  
%              cat_surf_saveICO(S,T,Pcentral,subdir,writeTfs,C)
%
%  * Other functions:
%    <a href="matlab:help cat_surf_fun>cat_suf_smoothtexture;">smoothcdata</a>   Smooth surface data.  
%              
% 
%
%  * Helping functions: 
%    <a href="matlab:help cat_surf_fun>cat_surf_area;">area</a>    Estimate face surface area. 
%              A = cat_surf_fun('area',S);  
%    <a href="matlab:help cat_surf_fun>cat_surf_graph2edge;">graph2edge</a>    Extract edges of a triangulation T. 
%              E = cat_surf_fun('graph2edge',T); 
%
%  * Test functions:
%    <a href="matlab:help cat_surf_fun>cat_surf_createEdgemap;">cdatamappingtst</a> 
%
% See also spm_mesh_* functions.
% _________________________________________________________________________
% Robert Dahnke, 2016-2019
% $Id: cat_surf_fun.m 1577 2020-03-09 17:36:03Z dahnke $ 

% See also spm_mesh_area, spm_mesh_borders,
%   spm_mesh_calc, spm_mesh_clusters, spm_mesh_contour, spm_mesh_curvature,
%   spm_mesh_detect, spm_mesh_distmtx, spm_mesh_edges, spm_mesh_euler, 
%   spm_mesh_geodesic, spm_mesh_get_lm, spm_mesh_inflate, spm_mesh_join, 
%   spm_mesh_label, spm_mesh_max, spm_mesh_neighbours, spm_mesh_normals, 
%   spm_mesh_polyhedron, spm_mesh_project, spm_mesh_reduce, spm_mesh_refine
%   spm_mesh_resels, spm_mesh_smooth, spm_mesh_sphere, spm_mesh_split, 
%   spm_mesh_to_grid, spm_mesh_transform, spm_mesh_utils.
%   spm_mesh_adjacency

%#ok<*ASGLU,*AGROW>

  switch lower(action)
    case 'normals'
      if nargin<2, help cat_surf_fun>cat_surf_normals; return; end
      varargout{1} = cat_surf_normals(S); 
    
    case 'dist'
      if nargin<2, help cat_surf_fun>cat_surf_dist; return; end
      varargout{1} = cat_surf_dist(S);
      
    case 'surf2surf'
      % create mapping between similar surfaces for area projection
      % using Delaunay (will be slow for large surface >1'000'000
      if nargin<2, help cat_surf_fun>cat_surf_surf2surf; return; end
      varargout{1} = cat_surf_surf2surf(S,varargin{1});
      
    case 'area'
       % simple area estimation of S 
       if nargin<2, help cat_surf_fun>cat_surf_area; return; end
       varargout{1} = cat_surf_area(S);
    
    case {'smoothcdata','smoothtexture'}
      % use CAT smoothing rather than spm_mesh_smooth 
      % differences? 
      if nargin<2, help cat_surf_fun>cat_surf_smoothtexture; return; end
      switch nargin
        case 2
          if isfield(S,'cdata'), varargout{1} = cat_surf_smoothtexture(S,S.cdata,1); end
        otherwise, varargout{1} = cat_surf_smoothtexture(S,varargin{:});
      end
      
    case 'maparea'
      % do the final area projection
      if nargin<2, help cat_surf_fun>cat_surf_maparea; return; end
      if nargout==1, varargout{1} = cat_surf_maparea(S,varargin{:}); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_maparea(S,varargin{:}); end
    
    case 'hull'
      % create a hull surface
      if nargin<2, help cat_surf_fun>cat_surf_hull; return; end
      if nargout==1, varargout{1} = cat_surf_hull(S); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_hull(S); end
    
    case 'core'
      % create a core surface
      if nargin<2, help cat_surf_fun>cat_surf_core; return; end
      if nargout==1, varargout{1} = cat_surf_core(S,varargin{:}); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_core(S,varargin{:}); end
    
    case {'tfs','tmin','tmax'}
      % estimate the distance between two linked? surfaces
      if nargin<2, help cat_surf_fun>cat_surf_thickness; return; end
      if numel(varargin)==1
        if nargout==1
          varargout{1} = cat_surf_thickness(action,S,varargin{1}); 
        else
          cat_surf_thickness(action,S,varargin{1});
        end
      else
        if nargout==1
          varargout{1} = cat_surf_thickness(action,S); 
        else
          cat_surf_thickness(action,S);
        end
      end
      
    case {'inner','outer','white','pial','innervar','outervar','whitevar','pialvar'}
    % create different cortical surfaces
      if nargin<2, help cat_surf_fun>cat_surf_GMboundarySurface; return; end
      if numel(varargin)==1
        switch nargout % surface & texture input
          case 0, cat_surf_GMboundarySurface(action,S,varargin{:});
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S,varargin{:}); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S,varargin{:}); 
        end
      else % file input
        switch nargout
          case 0, cat_surf_GMboundarySurface(action,S);
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S); 
        end
      end
      
    case 'disterr'
      switch nargout
        case 0, cat_surf_disterr(S,varargin{:});
        case 1, varargout{1} = cat_surf_disterr(S,varargin{:});
        case 2, [varargout{1},varargout{2}] = cat_surf_disterr(S,varargin{:});
      end
        
      
    case 'evalcs'
    % evaluation of the central and inner/outer surfaces
      if nargin<2, help cat_surf_fun>cat_surf_evalCS; return; end
      varargout{1} = cat_surf_evalCS(S,varargin{:});
    
    case 'createinneroutersurface'
    % not a good name ... this is for the Laplacian approach
      if nargin<2, help cat_surf_fun>cat_surf_createinneroutersurface; return; end
      cat_surf_createinneroutersurface(S,varargin{:});
    
    case 'show_orthview'
    % show cortical surfaces in the orthview window
      if nargin<2, help cat_surf_fun>show_orthview; return; end
      cat_surf_show_orthview(S,varargin{:});
    
    case 'saveico'
    % save the inner, central, and outer surfaces
      if nargin<2, help cat_surf_fun>show_orthview; return; end
      cat_surf_saveICO(S,varargin{:});
      
    case 'collisioncorrection'
    % Delaunay-based correction of surface self-intersections - not working  
      if nargin<2, help cat_surf_fun>show_orthview; return; end
      [varargout{1},varargout{2},varargout{3}] = cat_surf_collision_correction(S,varargin{:});

    case 'collisioncorrectionry'
    % CAT_SelfIntersect-based correction of surface self-intersections   
      if nargin<2, help cat_surf_fun>show_orthview; return; end
      [varargout{1},varargout{2}] = cat_surf_collision_correction_ry(S,varargin{:});
    
    case 'collisioncorrectionpbt'
    % correction of self-intersections based on the PBT PP map
      if nargin<2, help cat_surf_fun>show_orthview; return; end
      [varargout{1},varargout{2}] = cat_surf_collision_correction_pbt(S,varargin{:});
    
    case 'vdist'
    % ?
      if nargin<2, help cat_surf_fun>show_orthview; return; end
      [varargout{1},varargout{2}] = cat_surf_vdist(S,varargin);
    
    case 'surf2vol'
    % Render surface (data) into volume space. See also spm_mesh_to_grid. 
      if nargin<2, help cat_surf_fun>show_orthview; return; end
      if nargout>2
        if nargin>2
          [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S,varargin{:});
        else
          [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S);
        end
      elseif nargout>1
        if nargin>2
          [varargout{1},varargout{2}] = cat_surf_surf2vol(S,varargin{:});
        else
          [varargout{1},varargout{2}] = cat_surf_surf2vol(S);
        end
      else
        if nargin>2
          varargout{1} = cat_surf_surf2vol(S,varargin{:});
        else
          varargout{1} = cat_surf_surf2vol(S);
        end
      end
    case 'smat'
    % Apply matrix transformation. See also spm_mesh_transform.
      if nargin<2, help cat_surf_fun>cat_surf_mat; return; end
      varargout{1} = cat_surf_mat(S,varargin{:});
    
    case 'graph2edge'
    % extract edges from a give graph
      if nargin<2, help cat_surf_fun>cat_surf_edges; return; end
      switch nargout
        case 0, cat_surf_edges(S); 
        case 1, varargout{1} = cat_surf_edges(S); 
        case 2, [varargout{1},varargout{2}] = cat_surf_edges(S); 
      end
      
    case 'cdatamappingtst'
      cat_surf_cdatamappingtst;
    
    case 'createedgemap'
    % ???
      if nargin<2, help cat_surf_fun>cat_surf_surf2surf; return; end
      varargout{1} = cat_surf_surf2surf(S,varargin{:});
    
    case 'useedgemap'
    % ???
      if nargin<2, help cat_surf_fun>cat_surf_maparea; return; end
      varargout{1} = cat_surf_maparea(S,varargin{:});
    
    case 'gmv'
    % ???
      if nargin<2, help cat_surf_fun>cat_surf_gmv; return; end
      varargout{1} = cat_surf_gmv(S,varargin{:});
    
    case 'cdatamapping'
    % ???
      if nargin<2, help cat_surf_fun>cat_surf_cdatamapping; return; end
      if nargin<3, varargin{3} = ''; end
      if nargin<4, varargin{4} = struct(); end
      if nargout>1
        [varargout{1},varargout{2}] = cat_surf_cdatamapping(S,varargin{:});
      else
        varargout{1} = cat_surf_cdatamapping(S,varargin{:});
      end
    case 'reduce'
      varargout{1} = cat_surf_reduce(S,varargin{:});
    case 'isocolors'
    % Similar to MATLAB isocolor function but faster and support to use mat
    % files. See also spm_mesh_project. 
      if nargin<2, help cat_surf_fun>isocolors2; return; end
      varargout{1} = cat_surf_isocolors2(S,varargin{:});
    
    otherwise
      error('Unknow action "%s"! Check "help cat_surf_fun".\n',action); 
  end
    
end


function S = cat_surf_reduce(S,red)
  Ptemp = tempname; 
  Ptmpo = [Ptemp 'o.gii'];
  Ptmpn = [Ptemp 'n.gii'];
  
  % save surface
  save(gifti(struct('faces',S.faces,'vertices',S.vertices)),Ptmpo,'Base64Binary');
  clear Sgii 
  
  % call  
  for i=1:3
    %fprintf('.');
    if i>1, pause(5); end
    if ispc
      cmds = sprintf('set PATH=%s; start ',fullfile(matlabroot,'bin')); 
    else
      cmds = sprintf('%s/',fullfile(matlabroot,'bin'));
    end
    cmd = sprintf(['matlab -nojvm -nosplash -nodisplay -r ' ...
      '"try, %s; catch, disp(''Error''); end; exit; "'],...
      sprintf(['addpath(''%s''); S = gifti(''%s''); '...
        'S = reducepatch( struct(''vertices'',S.vertices,''faces'',S.faces) , %g); ' ...
        'save(gifti(S),''%s'',''Base64Binary''); '],fullfile(spm('dir'),'toolbox','cat12'),Ptmpo,red,Ptmpn));
    evalc('[ST, RS] = system([cmds,cmd]); cat_check_system_output(ST,RS,1);');
    
    if exist(Ptmpn,'file'); break; end
  end
  delete(Ptmpo)
  if ~exist(Ptmpn,'file'); 
    error('cat_surf_fun:reducemesh','Failed %d times to reduce surface resolution.\n%s',i);
  end
  
  Sgii = gifti(Ptmpn); delete(Ptmpn); 
  S = struct('vertices',Sgii.vertices,'faces',Sgii.faces);
end


%% area mapping concept
%  ------------------------------------------------------------------------
%  We need two functions, one that create the mapping between two very 
%  close surfaces (cat_surf_surf2surf) and another one (cat_surf_maparea) 
%  that finally performs the mapping. 
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/04
function edgemap = cat_surf_surf2surf(S1,S2,normalize)
%  ------------------------------------------------------------------------
%  Create mapping between two surface by nearest neighbor search on a 
%  Delaunay graph.
%   
%    edgemap = cat_surf_surf2surf(S1,S2,normalize)
%
%  see also <a href="matlab:help
%  cat_surf_surf>cat_surf_maparea">cat_surf_surf>cat_surf_maparea</>
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/04

  if ~exist('normalize','var'), normalize=1; end

  % 0. Normalize input
  if normalize
    S1.vertices = S1.vertices/max(S1.vertices(:));
    S2.vertices = S2.vertices/max(S2.vertices(:)) * 0.98; % need slight difference for Delaunay! 
  end
  
  % 1. Delaunay triangulation for fast search
  D1 = delaunayn(double(S1.vertices)); 
  D2 = delaunayn(double(S2.vertices)); 
 
  % 2. find minimal relation between the vertices of S1 and S2 as nearest neighbor
  E1(:,1) = 1:size(S1.vertices,1);
  E2(:,2) = 1:size(S2.vertices,1);
  E1(:,2) = dsearchn(double(S2.vertices),D2,double(S1.vertices)); 
  E2(:,1) = dsearchn(double(S1.vertices),D1,double(S2.vertices)); 
  E = unique([E1;E2],'rows'); 
  
  % 3. estimate edge length as weighting function 
  EL = sum( ( S1.vertices(E(:,1),:) - S2.vertices(E(:,2),:) ).^2 , 2) .^ 0.5; 
  
  % 4. Estimate the number of all by c-function
  %    - outgoing edges of S1 (connections from v element of S1) and
  %    - incoming edges of S2 (connections from v element of S2) 
  nE1 = zeros(max(E(:,1)),1,'single'); EL1 = zeros(max(E(:,1)),1,'single');
  nE2 = zeros(max(E(:,2)),1,'single'); EL2 = zeros(max(E(:,2)),1,'single');
  for i=1:size(E,1)
    nE1( E(i,1) ) =  nE1( E(i,1) ) + 1;
    EL1( E(i,1) ) =  EL1( E(i,1) ) + EL(i);
    nE2( E(i,2) ) =  nE2( E(i,2) ) + 1; 
    EL2( E(i,2) ) =  EL2( E(i,2) ) + EL(i);
  end
  
  % 5. Create a weighting function to project data from Si2St and St2Si.
  edgemap.edges     = E;
  edgemap.dist      = EL;
  edgemap.nvertices = [size(S1.vertices,1),size(S2.vertices,1)];
  edgemap.nforward  = 1  ./ nE1(E(:,1)); 
  edgemap.nbackward = 1  ./ nE2(E(:,2));
  edgemap.dforward  = EL ./ EL1(E(:,1)); 
  edgemap.dbackward = EL ./ EL2(E(:,2)); 
end

function varargout = cat_surf_maparea(varargin)
%  Apply graph-based mapping
%  ------------------------------------------------------------------------
%  use a c-function to process cat_surf_surf2surf mapping function 
%
%   cdata = cat_surf_maparea(cdatain,edgemap[,weighting])
%  
%   cdata     .. texture values at the output surface
%   edgemap   .. mapping structure between two surfaces
%   direction .. direction of cdata mapping if both surfaces have cdata 
%                with direction: 'forward' == '', 'backward' == 'invers'
%   weighting .. type of weighting: 
%                 'num'  .. by number of vertices
%                 'dist' .. by distance to the vertices (default)
%
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/04

  cdata     = varargin{1}; 
  edgemap = varargin{2};
  if nargin>2 
    dir = varargin{3}; 
  else
    dir = '';
  end
  switch dir
    case {'','forward'},        idir = 0; 
    case {'invers','backward'}, idir = 1;
    otherwise, error('Unkown mapping direction %s.\n',dir);
  end

  varargout{1} = cat_surf_edgemap(edgemap,cdata,idir);
  
  % filter with 1/mean(edgelegnth?
  % check for area?
end

function cdata2 = cat_surf_edgemap(edgemap,cdata,idir)
  if idir==0
    cdata2 = zeros(edgemap.nvertices(2),1,'single');
    for i=1:size(edgemap.edges,1)
      cdata2(edgemap.edges(i,2)) = cdata2(edgemap.edges(i,2)) + ...
        cdata(edgemap.edges(i,1)) * edgemap.dforward(i);
    end
  else
    cdata2 = zeros(edgemap.nvertices(1),1,'single');
    for i=1:size(edgemap.edges,1)
      cdata2(edgemap.edges(i,1)) =  cdata2(edgemap.edges(i,1)) + ...
        cdata2(edgemap.edges(i,2)) * edgemap.dbackward(i);
    end    
  end
end

function [IS,OS] = cat_surf_createinneroutersurface(S,T,Yp0)
  if ~exist('Yp0','var')
    % render surface
    
    
    Yp0 = 0; 
  end
  
  % call laplace 
  L = cat_surf_laplace(Yp0);
  
  % create streamlines
  IS.vertices = cat_surf_steams(L  ,T/2);
  OS.vertices = cat_surf_steams(1-L,T/2);
end

function VV = cat_surf_gmv(IS,OS)
% Estimate the volume between the given inner and outer boundary surface.
% 
%   VV = cat_surf_gmv(IS,OS)
%
%   VV .. volume between IS and OS
%   IS .. inner surface 
%   OS .. outer surface 
%
% Robert Dahnke 201904

  IV = IS.vertices*0.45 + 0.55*OS.vertices; 
  OV = OS.vertices*0.45 + 0.55*OS.vertices;

  % create Delaunay triangulation 
  D  = delaunayn(double([IV;OV])); clear IV OV;  
  
  % #############
  % classify and remove non GM tetraeder that have only WM points (ok) but
  % only GM? points ( this will not work for all gyri/sulci)
  % #############
  % - you maybe can use the centerpoint and then ...
  % - or you just use two verly lighty displaced surfaces, ie. 0.01 mm
  %   thickness that would give a Delaunay triangulation without the bad
  %   gyral effects!
  % #############
  DS = D>size(IS.vertices,1);
  D( sum(DS,2)==0 | sum(DS,2)==4 ,:)  = [];  
  clear DS;
  
  % estimate tetraeder volume
  DV = tetraedervolume(D,double([IS.vertices;OS.vertices]));
  
  %% map volume to faces
  VV = zeros(size(IS.vertices,1),1,'single'); 
  DF = D; DF(DF>size(IS.vertices,1)) = DF(DF>size(IS.vertices,1)) - size(IS.vertices,1); % to use the IS indices for mapping 
  for i=1:size(DV,1)
    for j=1:4
      VV(DF(i,j)) = VV(DF(i,j)) + DV( i ) / 4;  
    end
  end
 
end

function DV = tetraedervolume(D,V)
% estimate tetraeder volume by the Cayley-Menger determinant
% Robert Dahnke 201904

  % edgelength
  r = sum( ( V(D(:,1),:) - V(D(:,2),:) ).^2 , 2).^0.5; 
  p = sum( ( V(D(:,2),:) - V(D(:,3),:) ).^2 , 2).^0.5; 
  q = sum( ( V(D(:,3),:) - V(D(:,1),:) ).^2 , 2).^0.5; 
  a = sum( ( V(D(:,1),:) - V(D(:,4),:) ).^2 , 2).^0.5; 
  b = sum( ( V(D(:,2),:) - V(D(:,4),:) ).^2 , 2).^0.5; 
  c = sum( ( V(D(:,3),:) - V(D(:,4),:) ).^2 , 2).^0.5; 

  % volume
  DV = zeros(size(D,1),1);
  for i=1:size(D,1)
    DM = [ 0    r(i) q(i) a(i) 1;...
           r(i) 0    p(i) b(i) 1;...
           q(i) p(i) 0    c(i) 1;...
           a(i) b(i) c(i) 0    1;...
           1    1    1    1    0];
    DV(i) = sqrt(det( DM .^2 ) / 288);
  end
end

function varargout = cat_surf_GMboundarySurface(type,varargin)
%
%   ... = cat_surf_GMboundarySurface(type,Ps,Pth)
%
%   type .. projection direction
%           ['inner'|'outer'|'white'|'pial'] .. varargout is a filename
%           ['innervar'|'outervar'|'whitevar'|'pialvar'] 
%             .. varargout is a variable
%   Ps   .. filename of a given surface
%   Pth  .. thickness of the surface 
%

  if strfind(type,'var')
    varout=1; type = strrep(type,'var',''); 
  else
    varout=0;
  end
  switch type
    case {'white','inner'}, direction = -0.5;
    case {'pial' ,'outer'}, direction =  0.5;
  end
  
  if nargin>=2
    %% use filenames
    [pp,ff,ee] = spm_fileparts(varargin{1}); 
    
    if strcmp(ee,'')
      Praw = cat_io_FreeSurfer('fs2gii',varargin{1}); 
      Praw = Praw{1};
    else
      Praw   = varargin{1};
    end
    if nargin==3
      Pthick = varargin{2};
    else
      Pthick = cat_io_strrep(Praw,{'central','.gii'},{'pbt',''});
      if ~exist(Pthick,'file') 
        Pthick = cat_io_strrep(Praw,{'central','.gii'},{'thickness',''});
      end
      if ~exist(Pthick,'file') && exist(cat_io_strrep(Praw,'central','pbt'),'file')
        Pthick = cat_io_strrep(Praw,{'central','.gii'},{'pbt',''});
        movefile(cat_io_strrep(Praw,{'central'},{'pbt'}),Pthick);
      end
    end
    Ptype  = cat_io_strrep(Praw,'central',type);
    
    cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" %0.2f',Praw,Pthick,Ptype,direction); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,1);

    if strcmp(ee,'')
      Ptype = cat_io_FreeSurfer('gii2fs',Ptype); 
    end
    
    % filename
    if varout
      % load surface 
      varargout{1} = gifti(Ptype); 

      % delete temp files
      delete(Ptype);
    else
      varargout{1} = Ptype; 
    end
  else
    % write temp files ...
    Praw   = 'central.';
    Pthick = strrep(Praw,'central','pbt');
    Ptype  = strrep(Praw,'central',type);
   
    cmd = sprintf('CAT_Central2Pial "%s" "%s" %0.2f',Praw,Pthick,direction); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,1);
    
    % load surface 
    varargout{1} = gifti(Ptype); 
    
    % delete temp files
    delete(Praw,Pthick,Ptype);
  end
end

function cat_surf_cdatamappingtst
% ??? 
%  need at least some input
%
%


%% Testdata
   Psubcentral  = ['/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cg_vbm_defaults_template/template_NKI/'...
     'surf/lh.central.NKI_HC_NKI_1013090_T1_SD000000-RS00.gii'];
   PsubsphereA  = strrep(Psubcentral,'central','sphere.reg');              
   %Psubthick    = strrep(strrep(Psubcentral,'central','pbt'),'.gii','');               
   Psubthickres = strrep(strrep(Psubcentral,'central','thickness.resampled'),'lh.','s15mm.lh.'); 
   Psubtmp      = strrep(Psubcentral,'central','tmp'); 
   Pavgtmp      = strrep(strrep(Psubcentral,'central','tmp.resampled'),'lh.','s15mm.lh.'); 
 
   %Pavgcentral  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.freesurfer.gii'));
   PavgsphereA  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.sphere.freesurfer.gii'); 
   PavgDKT40    = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces','lh.aparc_DKT40JT.freesurfer.annot');
   
%% Test 1 - avg2sub - ok
   Ssub = gifti(PsubsphereA);
   Savg = gifti(PavgsphereA); 
   [vertices, label, colortable]  = cat_io_FreeSurfer('read_annotation',PavgDKT40); 
   Savg.cdata = label; 
   
   S3 = gifti(Psubcentral); 
   S3.cdata = cat_surf_fun('cdatamapping',Ssub,Savg,'nearest');
   save(gifti(S3),Psubtmp);
   
%% Test2 - sub2avg - ok
   Savg = gifti(PavgsphereA); 
   Ssub = gifti(PsubsphereA);
   %Ssub.cdata = cat_io_FreeSurfer('read_surf_data',Psubthick); 
   Ssub.cdata = cat_surf_fun('area',gifti(Psubcentral));
   
   S3 = gifti(Psubthickres); 
   mapping = {'directed'}; %,'undirected'}; %'nearest',
   for mi = 1:numel(mapping)
     S3.cdata  = cat_surf_fun('cdatamapping',Savg,Ssub,mapping{mi},1);
     S3.cdata  = spm_mesh_smooth(struct('vertices',S3.vertices,'faces',S3.faces),double(S3.cdata'),5);
     fprintf('mapping = %10s: A(sub) = %0.2f vs. A(avg) = %0.2f\n',mapping{mi},sum(Ssub.cdata(:)),sum(S3.cdata(:))); 
     save(gifti(S3),Pavgtmp); cat_surf_display(Pavgtmp)
   end
   
end

function varargout = cat_surf_cdatamapping(S1,S2,cdata,opt) 
% nearest connection between to surfaces
  if ischar(S1), S1 = gifti(S1); end
  if ischar(S2), S2 = gifti(S2); end
  if ischar(cdata)
    Pcdata = cdata;
    [pp,ff,ee] = spm_fileparts(cdata); 
    switch ee
      case '.annot'
        [vertices, cdata]  = cat_io_FreeSurfer('read_annotation',Pcdata); 
        clear vertices
      case '.gii'
        Scdata = gifti(S2); 
        if isfield(Scdata,'cdata')
          cdata = SX.cdata;
        else
          error('cat_surf_fun:cdatamapping:noTexture','No texture found in "%s"!\n',Pcdata);
        end
      otherwise
        cdata =  cat_io_FreeSurfer('read_surf_data',Pcdata);   
    end
  end
  
  if ~exist('cdata','var') || isempty(cdata)
    if isfield(S2,'cdata'), cdata = S2.cdata; end
  end
  
  if ~exist('opt','var'), opt = struct(); end
  def.method = 'nearest';
  def.verb   = 0; 
  def.smooth = 0; 
  opt        = cat_io_checkinopt(opt,def);
  
  if opt.verb, stime1 = cat_io_cmd(sprintf('Data-mapping (%s)',method)); fprintf('\n'); end
  
  % prepare vertices
  S1.vertices = S1.vertices ./ repmat(max(S1.vertices),size(S1.vertices,1),1)*1.1; % *100 
  S2.vertices = S2.vertices ./ repmat(max(S2.vertices),size(S2.vertices,1),1); 
  verticesS1  = double(S1.vertices - repmat(mean(S1.vertices),size(S1.vertices,1),1)); 
  verticesS2  = double(S2.vertices - repmat(mean(S2.vertices),size(S2.vertices,1),1)); 
  
  
  % estimate mapping
  switch opt.method
    case {'nearest'}
      [varargout{2},varargout{3}] = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
      varargout{1} = cdata(varargout{2}); 
    case {'undirected','directed'}
      %% use the surface as delauny graph
      switch opt.method 
        case 'directed'
          if opt.verb,  stime = cat_io_cmd('  Edge-Estimation (Nearest)','g5',''); end
          nextS2fromS1 = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
          nextS1fromS2 = dsearchn([verticesS1;inf(1,3)],double([S1.faces ones(size(S1.faces,1),1)*(size(S1.vertices,1)+1)]),verticesS2);
          tmp = nextS1fromS2; nextS1fromS2 = nextS2fromS1; nextS2fromS1 = tmp;
          nearestedges = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
          nearestedges = unique(nearestedges,'rows');
        case 'undirected'
          if opt.verb,  stime = cat_io_cmd('  Edge-Estimation (Delaunay','g5',''); end
          % nearest is required too
          nextS2fromS1 = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
          nextS1fromS2 = dsearchn([verticesS1;inf(1,3)],double([S1.faces ones(size(S1.faces,1),1)*(size(S1.vertices,1)+1)]),verticesS2);
          tmp = nextS1fromS2; nextS1fromS2 = nextS2fromS1; nextS2fromS1 = tmp;
          nearestedges  = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
          nearestedges1 = unique(nearestedges,'rows');
          % delauany
          triangulation = delaunayn([verticesS2;verticesS1]);              % delaunay triangulation
          nearestedges  = cat_surf_fun('graph2edge',triangulation);        % get edges 
          nearestedges(sum(nearestedges<=size(verticesS2,1),2)~=1,:)=[];   % only edges between S1 and S2
          nearestedges(:,2) = nearestedges(:,2) - size(verticesS2,1); 
          nearestedges = unique([nearestedges;nearestedges1],'rows');
      end
      if opt.verb, stime = cat_io_cmd('  Weighting','g5','',1,stime); end
      
      if 0
        %% my little testset
        nextS1fromS2 = [1; 1; 3; 4; 4; 4; 5; 5]; 
        nextS2fromS1 = [1; 3; 3; 5; 8; 8];
        cdata        = [1 1 1 1 1 1]';
        nearestedges = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
        nearestedges = unique(nearestedges,'rows');
      end
      
      
      %% simplify edges 1
      if 0
        % simpler, but much slower 
        nearestedges = [nearestedges, ones(size(nearestedges,1),1)]; % default weight
        [NeighborsS1,NidS1]  = hist(nearestedges(:,1),1:1:max(nearestedges(:,1)));
        for ni=NidS1(NeighborsS1>1)
          NumNi = nearestedges(:,1)==ni; 
          nearestedges(NumNi,3) =  nearestedges(NumNi,3) ./ sum(NumNi);
        end
      else
        % faster 
        %nearestedges = [nearestedges, ones(size(nearestedges,1),1)]; % default weight
        dist = sum( (S2.vertices(nearestedges(:,1),:) - S1.vertices(nearestedges(:,2),:)).^2 , 2) .^ 0.5; 
        nearestedges = [nearestedges, dist]; % default weight
        list = [1; find(nearestedges(1:end-1,1)~=nearestedges(2:end,1))+1; size(nearestedges,1)]; 
        for ni=1:numel(list)-1
          %nearestedges(list(ni):list(ni+1)-1,3) = nearestedges(list(ni):list(ni+1)-1,3) ./ (list(ni+1) - list(ni)); 
          nearestedges(list(ni):list(ni+1)-1,3) = nearestedges(list(ni):list(ni+1)-1,3) ./ sum(nearestedges(list(ni):list(ni+1)-1,3)); 
        end
      end
      if opt.verb, stime = cat_io_cmd('  Mapping','g5','',1,stime); end

      %%
      if 0
        % correct & simple, but very slow
        varargout{1} = zeros(1,max(nearestedges(:,2)));
        for ni=1:size(nearestedges,1)
          varargout{1}(nearestedges(ni,2)) = varargout{1}(nearestedges(ni,2)) + ...
            cdata(nearestedges(ni,1))' .* nearestedges(ni,3)';
        end
      else
        varargout{1} = zeros(1,max(nearestedges(:,2)));
        if 0
          list = [1; find(nearestedges(1:end-1,2)~=nearestedges(2:end,2))+1; size(nearestedges,1)+1]; 
          for ni=1:numel(list)-1
            varargout{1}(nearestedges(list(ni),2)) = varargout{1}(nearestedges(list(ni),2)) + ...
              sum(cdata(nearestedges(list(ni):list(ni+1)-1,1)) .*  nearestedges(list(ni):list(ni+1)-1,3));
          end
        else
          nearestedges2 = sortrows([nearestedges(:,2) nearestedges(:,1) nearestedges(:,3)]);  
          list = [1; find(nearestedges2(1:end-1,1)~=nearestedges2(2:end,1))+1; size(nearestedges2,1)+1]; 
          for ni=1:numel(list)-1
            varargout{1}(nearestedges2(list(ni),1)) = varargout{1}(nearestedges2(list(ni),1)) + ...
              sum(cdata(nearestedges2(list(ni):list(ni+1)-1,2)) .*  nearestedges2(list(ni):list(ni+1)-1,3));
          end
        end
      end
      if numel(varargout{1})<20, disp(varargout{1}); end
      if opt.verb, cat_io_cmd(' ','g5','',1,stime); end
  end
  
  % default smoothing???
  if opt.smooth
    varargout{1}  = spm_mesh_smooth(struct('vertices',S3.vertices,'faces',S3.faces),double(varargout{1}'),opt.smooth);
  end
  
  if isfield(opt,'fname')
    save(gifti(struct('vertices',S1.vertices,'faces',S1.faces,'cdata',varargout{1})),opt.fname); 
  end
  
  if opt.verb, cat_io_cmd('','','',1,stime1); end
end

function [E,uE] = cat_surf_edges(T)
% Extract edges of a given surface struture or its faces
%
%   [E,uE] = cat_surf_edges(T)
%
%   T .. surface structure or faces
%   E .. edges
%

  if isstruct(T) && isfield(T,'faces')
    T = T.faces;
  end

  T = sort(T,2); E = []; 
  for i=1:size(T,2)-1
    E = [E; T(:,[i i+1])]; 
  end
  [E,uE] = unique(E,'rows');
end

function D = cat_surf_dist(S)
% Estimate the distance between the vertices of the faces of S.
%
%   D = cat_surf_dist(S)
%
%   D = [c,a,b] = [d(AB),d(BC),d(CA)]
%

  D = [sum( (S.vertices(S.faces(:,1),:) - S.vertices(S.faces(:,2),:)).^2 , 2) .^ 0.5, ...
       sum( (S.vertices(S.faces(:,2),:) - S.vertices(S.faces(:,3),:)).^2 , 2) .^ 0.5, ...
       sum( (S.vertices(S.faces(:,3),:) - S.vertices(S.faces(:,1),:)).^2 , 2) .^ 0.5]; 
     
end

function [AV,AF] = cat_surf_area(S)
% Calculate surface area of the faces AF (Horonsche Form) and map it to the
% vertices AV.
%
%   [AV,AF] = cat_surf_area(S)
%
%   AV .. area per vertex (1/3 of all connected faces)
%   AF .. area per face
%   S  .. surface structure with vertices and faces
%

  % facearea (Horonsche Form)
  method = 2;
  if method == 1
    %%
    D = cat_surf_dist(S);
    facesp = sum(D,2) / 2;  % s = (a + b + c) / 2;
    AF = (facesp .* (facesp - D(:,1)) .* (facesp - D(:,2)) .* (facesp - D(:,3))).^0.5; % area=sqrt(s*(s-a)*(s-b)*(s-c));

    % numerical (to small point differences) and mapping problems (crossing of streamlines)
    % -> correction because this is theoretical not possible (Laplace field theory)
    AF(AF==0) = eps; % to small values
    AF = abs(AF);    % streamline-crossing

    AV = cat_surf_F2V(S,AF);
  elseif method == 2
    %%
    AF = spm_mesh_area(S,1);
    AV = cat_surf_F2V(S,AF);
  else
    %% edge points
    vn   = size(S.vertices,1);
    V1   = S.vertices(S.faces(:,1),:); 
    V2   = S.vertices(S.faces(:,2),:);
    V3   = S.vertices(S.faces(:,3),:);
    V12  = mean( cat( 3 , V1 , V2 ) , 3); 
    V13  = mean( cat( 3 , V1 , V3 ) , 3); 
    V23  = mean( cat( 3 , V2 , V3 ) , 3); 
    %V123 = circumcenter( struct( 'Points' , [ V1 ; V2 ; V3 ] , 'ConnectivityList' , [ (1:vn)' (1:vn)'+vn (1:vn)'+vn*2 ] )); % center of mass
    V123 = circlefit3d( V1 , V2 , V3 );

%%   
    SR11.faces = S.faces; SR11.vertices = [V1 V12 V123]; 
    SR12.faces = S.faces; SR12.vertices = [V1 V13 V123]; 
    SR21.faces = S.faces; SR21.vertices = [V2 V12 V123]; 
    SR22.faces = S.faces; SR22.vertices = [V2 V23 V123];
    SR31.faces = S.faces; SR31.vertices = [V3 V13 V123]; 
    SR32.faces = S.faces; SR32.vertices = [V3 V23 V123];
    
    %% outer circle point
    AF3 = [ abs(spm_mesh_area(SR11,1))' + abs(spm_mesh_area(SR12,1))' , ...
            abs(spm_mesh_area(SR21,1))' + abs(spm_mesh_area(SR22,1))' , ...
            abs(spm_mesh_area(SR31,1))' + abs(spm_mesh_area(SR32,1))'];
    
     AF = sum( AF3( S.faces ) , 2); 
     AV = cat_surf_F2V(S,AF); 
     
     % the COM is equal to divide by three and the estimation of the
     % circumcircle point was not so easy or fast possible as I hoped ...
     %   https://en.wikipedia.org/wiki/Triangle#Computing_the_area_of_a_triangle
     %   https://en.wikipedia.org/wiki/Circumscribed_circle
  end
   %%
%   AV = cat_mesh_smooth(S,AV,5); 

end

function data = cat_surf_F2V(S,odata)
%% mapping of facedata to vertices

  %A  = spm_mesh_distmtx(S,0);
  %A  = full(sum(A,2)); 
  %FA = A(S.faces);
  %FA = (1./FA.^16); 
  %FA = FA ./ repmat(sum(FA,2),1,3); 
  
  data   = zeros(size(S.vertices,1),1);
  [v,f]  = sort(S.faces(:)); 
  [f,fj] = ind2sub(size(S.faces),f);  
  far    = odata(f);
%  FA     = FA(f,:);
  for i=1:numel(v)
    data(v(i)) = data(v(i)) + far(i) / 3;% .* FA(i); %
  end

  %  data = data ./ vcount; %size(S.vertices,2); % Schwerpunkt... besser Voronoi, aber wie bei ner Oberflaeche im Raum???
  
end

function A = cat_surf_smoothtexture(S,A,smooth,Amax)
%  create smooth area texture files
%  ---------------------------------------------------------------------
  debug = 0;
  
  if ~exist('smooth','var'), smooth=1; end

  % temporary file names
  Pname  = tempname; 
  Pmesh  = [Pname 'mesh'];
  Parea  = [Pname 'area'];

  if exist('Amax','var');  A = min(A,Amax); end 
  
  % write surface and textures
  cat_io_FreeSurfer('write_surf',Pmesh,S);
  cat_io_FreeSurfer('write_surf_data',Parea,A);
  
  % smooth textures
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pmesh,Parea,smooth,Parea);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,debug);
  
  % load smoothed textures
  A  = cat_io_FreeSurfer('read_surf_data',Parea);
  
  % delete temporary file
  delete(Parea);
end

function [SH,V] = cat_surf_hull(S)
%  Create the hull surface SH and volume V of a given surface structure S.
%
%  [SH,V] = cat_surf_hull(S)
%
%  SH .. hull surface (different from S!)
%  V  .. volume of S
%  S  .. input surface
%

  % render surface points
  Vi = cat_surf_fun('surf2vol',S);
  
  % fill mesh
  V  = cat_vol_morph(Vi,'ldc',mean(size(Vi))/6); clear Vi; % closing 
  V  = cat_vol_smooth3X(V,2);    % smoothing
  SH = isosurface(V,0.4);        % create hull 
  V  = V>0.4;
  
  % final mesh operations
  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end

function PTN = cat_surf_thickness(action,PS,PT)
% Estimation of different cortical thickness metrics. 
%

  if ~exist('T','var')
    % create inner and outer surfaces
    PIS  = cat_surf_fun('inner',PS);  % estimate inner surface
    POS  = cat_surf_fun('outer',PS);  % estimate outer surface
  else
    % create inner and outer surfaces
    PIS  = cat_surf_fun('inner',PS,PT);  % estimate inner surface
    POS  = cat_surf_fun('outer',PS,PT);  % estimate outer surface
  end
  
  % load surfaces
  IS = gifti(PIS); 
  OS = gifti(POS);
  
  % edgemap
  % create mapping between 
  Pedgemap = cat_io_strrep(PS,{'.central.';'.gii'},{'.edgemapnative.';'.mat'});
  if 0%exist(Pedgemap,'file') % ... you have to test if central is older than the edgemap to use this 
    load(Pedgemap,'edgemap'); 
  else
    %%
    stime2  = clock;
    fprintf('  Estimate mapping for native surface');
    edgemap = cat_surf_surf2surf(IS,OS,0); 
    %edgemap.dist = sum ( (IS.vertices(edgemap.edges(:,1),:) - OS.vertices(edgemap.edges(:,2),:)).^2 , 2).^0.5;  
    %save(Pedgemap,'edgemap'); 
    fprintf(' takes %ds\n',round(etime(clock,stime2))); 
  end

  % create thickness metric mapping matrix
  switch lower(action)
    case {'tfs','tmin'}
      Tnear = inf(edgemap.nvertices(1),2,'single'); 
      for i=1:size(edgemap.edges,1)
        Tnear(edgemap.edges(i,1),1) = min( [ Tnear(edgemap.edges(i,1),1) edgemap.dist(i)  ] ) ;
        Tnear(edgemap.edges(i,2),2) = min( [ Tnear(edgemap.edges(i,2),2) edgemap.dist(i)  ] ) ;
      end
    case 'tmax'
      Tfar  = zeros(edgemap.nvertices(1),2,'single'); 
      for i=1:size(edgemap.edges,1)
        Tfar(edgemap.edges(i,1),1)  = max( [ Tfar(edgemap.edges(i,1),1) edgemap.dist(i)  ] ) ;
        Tfar(edgemap.edges(i,2),2)  = max( [ Tfar(edgemap.edges(i,2),2) edgemap.dist(i)  ] ) ;
      end
  end
    
  switch lower(action)
    case 'tfs'
      TN  = mean(Tnear,2);
      PTN = cat_io_strrep(PS,{'.central.';'.gii'},{'.thicknessfs.';''});
    case 'tmin'
      TN  = min(Tnear,[],2);
      PTN = cat_io_strrep(PS,{'.central.';'.gii'},{'.thicknessmin.';''});
    case 'tmax'
      TN  = max(Tfar,[],2);
      PTN = cat_io_strrep(PS,{'.central.';'.gii'},{'.thicknessmax.';''});
  end
  
  % save smoothed textures
  cat_io_FreeSurfer('write_surf_data',PTN,TN);
end

function [SH,V] = cat_surf_core(S,opt)
% _________________________________________________________________________
% Estimation of a core surface of central WM regions without gyri that is 
% the inverse of the hull surface in principle. 
% 
%   [SH,V] = cat_surf_core(S,opt)
%
%   SH      .. core surface 
%   V       .. volume map of the core surface 
%   S       .. input surface 
%   opt     .. parameter structure
%    .type  .. type of core creation 
%    .th    .. threshold for core creation 
%
% This function is in development!
% 
% This is much more complicated that the hull definition. So I will need 
% different types of core definitions. However, I first have to find one
% (or multiple) anatomical definitions.
%
%  For estimation I can use different techniques: 
%   * morphological operations 
%     > very inaccurate and error-prone 
%     > use of distance & smoothing functions 
%   * smoothing with/without boundaries
%   * anatomical information from volume or better surface atlas maps 
%   * use of other measures such as 
%     - thickness (no)
%     - sulcal depth or outward folding GI (maybe)
%     - curvature (not really)
%
%   * use of percentage scalings
%   * use of multiple threshold levels and averaging to avoid using only 
%     one threshold (=multiband) 
%   * definition as fractal dimension measure?
% _________________________________________________________________________

  def.type = 1; 
  def.th   = 0.15;
  opt = cat_io_checkinopt(opt,def); 
  
  
  % render surface points
  Vi = cat_surf_fun('surf2vol',S);
   
  %% break gyri
  if opt.type == 1
    %%
    Vd  = cat_vbdist(single(Vi<0.5));
    
    %%
    Vdn = Vd ./ max(Vd(:));
    Vdn = cat_vol_laplace3R(Vdn,Vdn>0 & Vdn<0.8 ,0.001);
    Vdn = min(opt.th + cat_vol_morph(Vdn>opt.th,'l'),Vdn); 
    V   = Vdn > opt.th; 
    
    %%
    SH  = isosurface(Vdn,opt.th);        % create hull 
  elseif opt.type == 2
    SiGI = S; SiGI.cdata = opt.iGI;
    V    = cat_surf_fun('surf2vol',SiGI,struct('pve',3));
    V    = V ./ max(V(:));
    SH   = isosurface(V,opt.th);        % create hull 
  else
    Vs = cat_vol_smooth3X(Vi,8);    % smoothing
    V  = ~cat_vol_morph(Vs<0.5,'ldc',min(size(Vi))/(6*1.5));  % opening 
    V  = cat_vol_smooth3X(V,6);    % smoothing
    SH = isosurface(V .* smooth3(Vi),0.6);        % create hull 
    V  = min(V>0.6,Vi==1);
    V  = cat_vol_smooth3X(V,4);    % smoothing
    V  = min(V,Vi);
    V  = cat_vol_laplace3R(V,Vi>0 & V<0.9,0.4);
  end   %clear Vi;

  % final mesh operations
  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end

function res = cat_surf_evalCS(CS,T,Ym,Ypp,Pcentral,mat,verb,estSI)
% _________________________________________________________________________
% cat_surf_evalCS in cat_surf_fun print out some values to characterize a
% cortical surface.
%
%   res = cat_surf_fun('evalCS',CS[,T,Ym,Yppt])
%   res = cat_surf_evalCS(CS[,T,Ym,Yppt])
%    
%   CS        .. central surface in world space (mm)
%   T         .. cortical thickness in world space (mm)
%   Ym        .. intensity normalized file with BG=0, CSF=1/3, GM=2/3, and WM=1
%   Ypp       .. percent position map 
%   Pcentral  .. number of classes for further thickness evaluation or a
%                given filename to detect specific thickness phantom rules
%   verb      .. print results (default = 1)
%   estSI     .. estimate self intersections (SI) .. slow! (default = 0)
%
%   res       .. structure with data fields of the printed values
% _________________________________________________________________________
% Robert Dahnke 201909

% - maybe also save and link (surface + hist) some files in future
% - the Layer4 handling with the global variables is horrible 
 
  QMC       = cat_io_colormaps('marks+',17);
  color     = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  rms       = @(x) mean( x.^2 ).^0.5;
  rate      = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
  
  if ~exist('verb','var'),  verb  = 1; end
  if ~exist('estSI','var'), estSI = 0; end
  
  M  = spm_mesh_smooth(CS);    % smoothing matrix
  if exist('T','var') % thickness in voxel space 
    N  = spm_mesh_normals(CS);
    
    VI = CS.vertices + N .* repmat(T / 2 ,1,3); % white surface
    VO = CS.vertices - N .* repmat(T / 2 ,1,3); % pial surface
  end
  
  if exist('Pcentral','var') && ischar(Pcentral)
    Player4   = strrep(Pcentral,'.central.','.layer4.'); 
    Pcentralx = strrep(Pcentral,'.central.','.centralx.'); 
    Player4x  = strrep(Pcentral,'.central.','.layer4x.'); 
    Ppbtx     = strrep(Pcentral(1:end-4),'.central.','.pbtx.'); 
  
    if ~isempty(mat); 
      save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),Pcentralx,'Base64Binary');
      cat_io_FreeSurfer('write_surf_data',Ppbtx,T);  
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0', ...
                     Pcentralx,Ppbtx,Player4x);
      try
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
        L4 = gifti(Player4x);
        uL4 = 1; 
        delete(Player4x);
      catch
        uL4 = 1; 
      end
      if exist(Pcentralx,'file'), delete(Pcentralx); end
      if exist(Ppbtx,'file'),     delete(Ppbtx); end
    else
      uL4 = 0;
    end
  else
    uL4 = 0;
  end
  
  
  
  
  %% Evaluation of local intensities
  %  Here we have to use the Layer 4 rather than the central surface.
  %  All values will depend on age!
  if exist('Ym','var')
    warning off MATLAB:subscripting:noSubscriptsSpecified  
    II = cat_surf_isocolors2(Ym,VI,mat);  
    IO = cat_surf_isocolors2(Ym,VO,mat); 
    % local adaption for GM intensity changes by myelination 
    IIs = single(spm_mesh_smooth(M,double(II),round(100 * sqrt(size(CS.faces,1)/180000)))); 
    IOs = single(spm_mesh_smooth(M,double(II),round(100 * sqrt(size(CS.faces,1)/180000)))); 
    % normalization
    II  = II./(IIs/mean(IIs)) - 2.5; clear IIs;
    IO  = IO./(IOs/mean(IOs)) - 1.5; clear IOs;
    % 
    
    if uL4
      ML  = spm_mesh_smooth(L4);    % smoothing matrix
      IC  = cat_surf_isocolors2(Ym,L4,mat);  
      ICs = single(spm_mesh_smooth(ML,double(IC),round(100 * sqrt(size(CS.faces,1)/180000)))); 
      IC  = IC./(ICs/mean(ICs)) - mean(ICs); clear ICs;
    end
    % divide by 2 because of the CSF-GM (1-2) and the GM-WM area (2-3) 
    % and to obtain a similar scaling as for the Ypp (also two segments and
    % we do not dived)
    II = II / 2; IO = IO / 2; IC = IC / 2;
    % output
    if verb
      fprintf('    Local intensity RMSE (lower=better): ')
      if uL4
        cat_io_cprintf( color( rate( mean( [rms(II),rms(IC),rms(IO)] ) , 0.05 , 0.30 )) , sprintf('%0.4f ',mean( [rms(II),rms(IC),rms(IO)] )) ); 
      else
        cat_io_cprintf( color( rate( mean( [rms(II),rms(IO)] ) , 0.05 , 0.30 )) , sprintf('%0.4f ',mean( [rms(II),rms(IO)] )) ); 
      end
      cat_io_cprintf( color( rate( rms(II) , 0.05 , 0.30 )) , sprintf('(IS=%0.4f,',rms(II)) ); 
      if uL4, cat_io_cprintf( color( rate( rms(IC) , 0.05 , 0.30 )) , sprintf('L4=%0.4f,',rms(IC)) ); end
      cat_io_cprintf( color( rate( rms(IO) , 0.05 , 0.30 )) , sprintf('OS=%0.4f)\n',rms(IO)) ); 
    end
    res.RMSE_Ym_white  = rms(II);
    if uL4, res.RMSE_Ym_layer4 = rms(IC); end
    res.RMSE_Ym_pial   = rms(IO);
    clear II IO; 
  end
  
  
  
  
  %% Evaluation of surface position values
  %  Here we can of course use the central surface
  %  This will be relative age independent.
  if exist('Ypp','var')
    II = cat_surf_isocolors2(Ypp,VI,mat);          
    IC = cat_surf_isocolors2(Ypp,CS,mat); 
    IO = cat_surf_isocolors2(Ypp,VO,mat);         
    II = II - 1.0;
    IC = IC - 0.5;
    IO = IO - 0.0;
    % output
    if verb
      fprintf('    Local position  RMSE (lower=better): '); 
      cat_io_cprintf( color( rate( mean( [rms(IC),rms(II),rms(IO)]) , 0.05 , 0.30 )) ,sprintf('%0.4f ',mean( [rms(IC),rms(II),rms(IO)] )) ); 
      cat_io_cprintf( color( rate( rms(II) , 0.05 , 0.30 )) , sprintf('(IS=%0.4f,' ,rms(II)) ); 
      cat_io_cprintf( color( rate( rms(IC) , 0.05 , 0.30 )) , sprintf('CS=%0.4f,'  ,rms(IC)) ); 
      cat_io_cprintf( color( rate( rms(IO) , 0.05 , 0.30 )) , sprintf('OS=%0.4f)\n',rms(IO)) ); 
    end
    res.RMSE_Ypp_white   = rms(II);
    res.RMSE_Ypp_pial    = rms(IO);
    res.RMSE_Ypp_central = rms(IC);
  end 
  
  
  
  
  %% CAT_SelfIntersect 
  %  This is very slow and we may want to keep the result. 
  if estSI && exist('Pcentral','var') && ischar(Pcentral)
    [pp,ff,ee] = spm_fileparts(Pcentral);
    
    Pwhite = fullfile(pp,strrep([ff ee],'central','whitex'));   
    Ppial  = fullfile(pp,strrep([ff ee],'central','pialx'));   
    Pselfw = fullfile(pp,strrep(ff,'central','whiteselfintersect'));
    Pselfp = fullfile(pp,strrep(ff,'central','pialselfintersect'));
    ePwhite = exist(Pwhite,'file'); 
    ePpial  = exist(Ppial, 'file');
    %ePselfw = exist(Pselfw,'file');
    %ePselfp = exist(Pselfp,'file');
    
    % save surfaces
    save(gifti(struct('faces',CS.faces,'vertices',VI)),Pwhite,'Base64Binary');
    save(gifti(struct('faces',CS.faces,'vertices',VO)),Ppial,'Base64Binary');

    % write self intersection maps
    cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Pwhite,Pselfw); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Ppial,Pselfp); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    
    selfw = cat_io_FreeSurfer('read_surf_data',Pselfw);
    selfp = cat_io_FreeSurfer('read_surf_data',Pselfp);

    area = cat_surf_area(CS);

    res.white_self_interection_area = sum((selfw(:)>0) .* area(:)) / 100;
    res.pial_self_interection_area  = sum((selfp(:)>0) .* area(:)) / 100;
    res.white_self_interections     = res.white_self_interection_area / sum(area(:)/100) * 100;
    res.pial_self_interections      = res.pial_self_interection_area  / sum(area(:)/100) * 100;

    if verb
      fprintf('    Self intersections (white,pial):     '); 
      cat_io_cprintf( color( rate( res.white_self_interections , 0 , 20 )) , ...
        sprintf('%0.2f%%%% (%0.2f cm%s) ',res.white_self_interections,res.white_self_interection_area,char(178))); 
      cat_io_cprintf( color( rate( res.pial_self_interections , 0 , 20 )) , ...
        sprintf('%0.2f%%%% (%0.2f cm%s)\n',res.pial_self_interections,res.pial_self_interection_area,char(178))); 
    end

    % delete temparary files
    if ~ePwhite, delete(Pwhite); end
    if ~ePpial,  delete(Ppial);  end 
    %if ~ePselfw, delete(Pselfw); end
    %if ~ePselfp, delete(Pselfp); end
  end
  
  
  
  
  %% thickness analysis
  if exist('T','var')
    if exist('Tclasses','var') && ~isempty(Pcentral)
      if ischar(Pcentral)
        if strfind(Pcentral,'dilated1.5-2.5mm')
          T = cat_stat_histth(T,0.95);
          Pcentral = 3;
        else
          Pcentral = 0;
        end
      end
      
      if Pcentral>0
        if Pcentral>7 || Pcentral<1
          warning('Tclasses has to be between 2 and 7.');
          Pcentral = min(7,max(3,Pcentral)); 
        end
        [mn,sd] = kmeans3D(T,Pcentral);
        if verb
          fprintf('    Thickness mean (%d class(es)):       ',Pcentral)
          fprintf('%7.4f',mn); fprintf('\n'); 
          fprintf('    Thickness std  (%d class(es)):       ',Pcentral)
          fprintf('%7.4f',sd); fprintf('\n');
        end
        res.thickness_mean_nclasses = mn;
        res.thickness_std_nclasses  = sd;
      end
    end
      
    res.thickness_mn_sd_md_mx = [mean(T),std(T),median(T),max(T)];
    if verb
      fprintf('    Thickness values:                    %0.4f%s%0.4f (md=%0.4f,mx=%0.4f)\n',...
        res.thickness_mn_sd_md_mx(1),native2unicode(177, 'latin1'),res.thickness_mn_sd_md_mx(2:end)); 
    end
  end
  
  
  
  
  % curvature analyis - not realy relevant 
  if 0
    C = abs(spm_mesh_curvature(CS));
    res.abscurv_mn_sd_md_mx = [mean(C),std(C),median(C),max(C)];
    if verb
      fprintf('    Abs mean curvature values:           %0.4f%s%0.4f (md=%0.4f,mx=%0.4f)\n',...
        res.abscurv_mn_sd_md_mx(1),native2unicode(177, 'latin1'),res.abscurv_mn_sd_md_mx(2:end)); 
    end
  end
  
  
  
  
  % surface values
  warning off MATLAB:subscripting:noSubscriptsSpecified  
  EC  = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
  res.euler_characteristic = EC; 
  if verb
    fprintf('    Faces / Euler:                       '); 
    cat_io_cprintf( color( rate( 1 - max(0,size(CS.faces,1)/300000) , 0 , 0.9 )) , sprintf('%d / ',size(CS.faces,1)));
    cat_io_cprintf( color( rate( abs(EC-2) , 0 , 30 )) , sprintf('%d',EC));
    fprintf('\n'); 
  end
end

function cat_surf_saveICO(S,Tpbt,Pcs,subdir,Pm,mat,writeTfs,writeSI,writeL4,writeInt)
% _________________________________________________________________________
% Save surface data for debugging:
% Creates and save the white and pial surfaces based on the displacement by
% the half thickness along the surface normals and use the inner and outer
% surfaces to create the layer4 surface.
% Saves also the thickness file (PBT by default).
%
% Writing of FreeSurfer thickness takes about 10s.
% Writing of self-intersection maps takes about 90s. 
% Writing of the Layer4 also take about 8s.
% Writing of the Layer4 also take about 6s if Ym is given else 90s.
%
%   cat_surf_saveICO(S,Tpbt,Pcs,subdir,Pm,writeTfs,writeSI,writeL4,writeInt)
%
%   S         .. central surface in voxel space
%   T         .. cortical thickness (in mm)
%   Pcs       .. central surface file name (with full path)
%   subdir    .. addition subdirectory in the standard surface directory
%   Pm        .. intensity information 
%   writeTfs  .. estimate FreeSurfer thickness metric
%   Pm        .. intensity file/volume to map data to the surfaces
%   writeTfs  .. write default thickness metric rather than given PBT 
%   writeSI   .. write self-intersection maps
%   writeL4   .. create L4 surface (and project data on it) 
%   writeInt  .. write some intensity maps
% _________________________________________________________________________
% Robert Dahnke 201908

  opt.verb = 0; 

  [pp,ff,ee] = spm_fileparts(Pcs);
  if ~exist('Pm','var'), Pm = ''; end
  if ~exist('writeTfs','var'), writeTfs = 0; end
  if ~exist('writeInt','var'), writeInt = 0; end
  if ~exist('writeSI','var'),  writeSI  = 0; end
  if ~exist('writeL4','var'),  writeL4  = 1; end
  if ~exist('subdir','var')
    subdir = '';
  else
    if ~exist(fullfile(pp,subdir),'dir')
      mkdir(fullfile(pp,subdir)); 
    end
  end
 
  %if isfield(S,'vmat') && isfield(S,'mati') 
  %  S.vertices = (S.vmat * [S.vertices' ; ones(1,size(S.vertices,1))])';
  %  if S.mati(7)<0, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; end
  %end
  
  
  % normalized surface normals
  N = spm_mesh_normals(S);                             
  %  matx = spm_imatrix( [S.vmat; 0 0 0 1] ); matx([1:3,9:12]) = 0; matx(6:9) = 1; % only rotate
  %  vmat = spm_matrix( matx ); vmat(4,:) = []; 
  %  N = -(vmat  * [N' ; ones(1,size(N,1))])';
  
  % inner and outer surface
  VI  = S.vertices + N .* repmat(Tpbt / 2,1,3); 
  VO  = S.vertices - N .* repmat(Tpbt / 2,1,3); 
 
  % surface filenames
  Pcentral = fullfile(pp,subdir,[ff ee]);   
  Pwhite   = fullfile(pp,subdir,strrep([ff ee],'central','white'));   
  Ppial    = fullfile(pp,subdir,strrep([ff ee],'central','pial'));   
  Pthick   = fullfile(pp,subdir,strrep(ff,'central','thickness'));   
  Ppbt     = fullfile(pp,subdir,strrep(ff,'central','pbt'));   
  PintIS   = fullfile(pp,subdir,strrep(ff,'central','intwhite'));
  PintOS   = fullfile(pp,subdir,strrep(ff,'central','intpial'));
  PintL4   = fullfile(pp,subdir,strrep(ff,'central','intlayer4'));
  Player4  = fullfile(pp,subdir,strrep([ff ee],'central','layer4'));   
  Pselfw   = fullfile(pp,subdir,strrep(ff,'central','whiteselfintersect'));
  Pselfp   = fullfile(pp,subdir,strrep(ff,'central','pialselfintersect'));

  % save surfaces
  save(gifti(struct('faces',S.faces,'vertices',S.vertices)),Pcentral,'Base64Binary');
  save(gifti(struct('faces',S.faces,'vertices',VI)),Pwhite,'Base64Binary');
  save(gifti(struct('faces',S.faces,'vertices',VO)),Ppial,'Base64Binary');

  % write self intersection maps
  if writeSI
    cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Pwhite,Pselfw); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Ppial,Pselfp); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
  end
  
  % save thickness
  cat_io_FreeSurfer('write_surf_data',Ppbt,Tpbt);
  if exist('writeTfs','var') && ~isempty(writeTfs) && writeTfs
    cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',Pwhite,Ppial,Pthick);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    fprintf('Display thickness: %s\n',spm_file(Pthick ,'link','cat_surf_display(''%s'')'));
  end
  
  % final correction of central surface in highly folded areas with high mean curvature
  if writeL4
    cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0', ...
                     Pcentral,Ppbt,Player4);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
  end
  
  % save intensities
  if isempty(Pm) && writeInt
    % volume filenames for spm_orthview
    sinfo = cat_surf_info(Pcentral);
    
    if cat_get_defaults('extopts.subfolders')
      Pm = fullfile(spm_str_manip(pp,'h'),'mri',['m' sinfo.name '.nii']); 
    else
      Pm = fullfile(pp,['m' sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = fullfile(spm_str_manip(pp,'hh'),'mri',['m' sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = fullfile(spm_str_manip(pp,'h'),[sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = ''; 
    end
  end
  
  if ~isnumeric( Pm ) && exist(Pm,'file')
    % use the file data ... slow????
    cmd = sprintf('CAT_3dVol2Surf -linear -steps 1 -start 0 -end 1 "%s" "%s" "%s"',Pwhite , Pm, PintIS);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    cmd = sprintf('CAT_3dVol2Surf -linear -steps 1 -start 0 -end 0 "%s" "%s" "%s"',Ppial  , Pm, PintOS);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    if writeL4
      cmd = sprintf('CAT_3dVol2Surf -linear -steps 1 -start 0 -end 0 "%s" "%s" "%s"',Player4, Pm, PintL4);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    end
  elseif ndims(Pm)==3
    int = cat_surf_isocolors2(Pm,VO,mat); cat_io_FreeSurfer('write_surf_data',PintOS,int);
    int = cat_surf_isocolors2(Pm,VI,mat); cat_io_FreeSurfer('write_surf_data',PintIS,int);
    if writeL4
      % use a given volume 
      SL = gifti(Player4);
      
      warning off MATLAB:subscripting:noSubscriptsSpecified
      %SL.vertices = (S.vmati * [SL.vertices' ; ones(1,size(SL.vertices,1))])';
      %if S.mati(7)<0,  SL.faces = [SL.faces(:,1) SL.faces(:,3) SL.faces(:,2)]; end
      %VL = SL.vertices;

      int = cat_surf_isocolors2(Pm,SL,mat); cat_io_FreeSurfer('write_surf_data',PintL4,int);
    end

    % define filename (same block as above)
    
    % volume filenames for spm_orthview
    sinfo = cat_surf_info(Pcentral);
    
    if cat_get_defaults('extopts.subfolders')
      Pm = fullfile(spm_str_manip(pp,'h'),'mri',['m' sinfo.name '.nii']); 
    else
      Pm = fullfile(pp,['m' sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = fullfile(spm_str_manip(pp,'hh'),'mri',['m' sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = fullfile(spm_str_manip(pp,'h'),[sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = ''; 
    end
    
  end
 
  % display something to click
  fprintf('\n    Display surface:  %s\n',spm_file(Ppbt  ,'link','cat_surf_display(''%s'')'));
  if ~isempty(Pm)
    fprintf('    Show in orthview: %s\n',spm_file(Pm ,'link',...
      [ sprintf('cat_surf_fun(''show_orthview'',{''%s'';''%s'';''%s''},',Pcentral,Ppial,Pwhite) '''%s'')']));
  end
  
  if 0
    if ndims(Pm)==3, Ym=Pm; else, Ym=spm_read_vols(spm_vol(Pm)); end
    res = cat_surf_evalCS(S,Tpbt,Ym,Ypp,Tclasses)
  end
end

function Sg = cat_surf_volgrad(varargin)
% _________________________________________________________________________
% This function estimates the local gradient in an image along the surface
% normals.
%
%  Sg = cat_surf_volgrad(S[,N],Y,mat[,d])
%
%  S      .. surface in voxel-space
%  N      .. surface normals in voxel-space
%  Y      .. voxel in voxel-space
%  mat    .. structure with resolution information 
%  d      .. distance along the surface-normal (default = 0.05)
% _________________________________________________________________________
% Robert Dahnke 201910

  d = 0.05;
  if isstruct(varargin{1})
    S = varargin{1};
    Y = varargin{2};
    mat = varargin{3};
    if nargin==4
      d = varargin{4};
    end
    
    N = spm_normals(S); 
    V = S.vertices; 
  else
    V = varargin{1}; 
    N = varargin{2}; 
    Y = varargin{3}; 
    mat = varargin{4}; 
  
    if nargin==5
      d = varargin{5};
    end
  end

  V1  = V - N*d;
  V2  = V + N*d;
  
  if exist('mat','var')
    Si1 = cat_surf_isocolors2(Y,V1,mat);
    Si2 = cat_surf_isocolors2(Y,V2,mat);
  else
    Si1 = cat_surf_isocolors2(Y,V1);
    Si2 = cat_surf_isocolors2(Y,V2);
  end
  
  Sg  = Si1 - Si2; 
end

function [S,Tn] = cat_surf_collision_correction_ry(S,T,Y,opt) 
% _________________________________________________________________________
% Correction of self-intersection (S) by iterative use of the 
% CAT_selfintersect function of Rachel Yotter / Christian Gaser.
% For typical 160k surfaces 60 to 90s are required per iteration and at
% least 5 iterations are required for reasonable quality (accuracy 1/2^4
% with aobut 0.25% SIs). 
%
%  [S,Tn] = cat_surf_collision_correction_ry(S,T,Y,opt) 
%
%  S          .. central surface in voxel space
%  T          .. thickness in world space 
%  Y          .. intensity normalized image (to estimate face normal orientation)
%  opt        .. structure with further options
%  .verb      .. print progress
%  .redterm   .. size of the reduction interval (default 2 = half interval) 
%  .accurarcy .. stop criteria (default 1/2^4 for about 5 iterations)
%  .iterfull  .. maximum number of iterations that is automatically defined 
%                by the accuracy 
%  .interBB   .. structure with resolution information 
%   .BB       .. not used here
%   .interpV  .. voxel-space resolution of S 
% _________________________________________________________________________
% Robert Dahnke 201910

  if ~exist('opt','var'), opt = struct(); end

  % defaults
  def.redterm  = 2;      % reduction term of the iteration, e.g. 2 = half resolution every iteration
  def.verb     = 1;      % print changes
  def.debug    = 2;      % use debugging
  def.mod      = 1;      % 1 - interval-based correction, 0 - linear correction   
  def.accuracy = 0.0625; % smallest used stepsize (the value should not be to small because there
                         % should also be a gab between the structures and a reasonable running time) 
  %def.interpBB = struct('BB',[],'interpV',1);
  def.mat      = []; 
  opt          = cat_io_checkinopt(opt,def); 
  opt.iterfull = find( 1 ./ opt.redterm.^(1:100)  < opt.accuracy , 1); % + 2 maximum number of iterations
  
  sf           = round( sqrt( size(S.faces,1) / 50000) ) * 1; % smoothing iterations depend on mesh size ### empirical value 
  M            = spm_mesh_smooth(S);                          % for spm_smoothing matrix

  % filenames
  [pp,ff,ee]   = spm_fileparts(opt.Pcs);
  Pwhite       = fullfile(pp,strrep([ff ee],'central','white'));   
  Ppial        = fullfile(pp,strrep([ff ee],'central','pial'));   
  Pselfw       = fullfile(pp,strrep(ff,'central','whiteselfintersect'));
  Pselfp       = fullfile(pp,strrep(ff,'central','pialselfintersect'));

  % inner and outer thickness seen from the central surface
  Tvxw  = T / 2;                    % white matter side thickness 
  Tvxp  = T / 2;                    % pial side thickness
  
  % detection and correction for flipped faces to have always the same normal direction
  flipped = cat_surf_checkNormalDir(S); %,T,Y,opt.interpBB);
  if flipped, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; S.mati(7) = - S.mati(7); end

  % simple surface smoothing
  smoothsurf = @(V,s) [ ...        
    spm_mesh_smooth(M,double(V(:,1)),s) , ...
    spm_mesh_smooth(M,double(V(:,2)),s) , ...
    spm_mesh_smooth(M,double(V(:,3)),s) ];
  
  Tn = T; 
  N    = spm_mesh_normals(S); 
  i    = 0; final = 0; 
  
  if opt.verb, fprintf('\n'); end
  while i < opt.iterfull
    i = i + 1;
   
    % In theory the search/correction size would be half each time (0.5 0.25 0.125 ... 2^i)
    % but due to the smoothing a slower sloop is used (1.8^i). 
    % The changes are limited by opt.accuracy for anatomical (min sulcus-width) and runtime reasons.
    % If optimization is used than the test starts again but with higher accuracy ( 0.125 ). 
    if opt.mod
      corrsize = max( opt.accuracy , 1 ./ opt.redterm.^(i+1) );
    else
      corrsize = opt.accuracy;
    end
    
    % outer and inner boundary
    VO = S.vertices - N .* repmat(Tvxp,1,3);
    VI = S.vertices + N .* repmat(Tvxw,1,3); 

    save(gifti(struct('faces',S.faces,'vertices',VI)),Pwhite,'Base64Binary');
    save(gifti(struct('faces',S.faces,'vertices',VO)),Ppial ,'Base64Binary');

    % call self-intersection function 
    cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Pwhite,Pselfw); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Ppial,Pselfp); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);

    selfw = cat_io_FreeSurfer('read_surf_data',Pselfw)>0; 
    selfp = cat_io_FreeSurfer('read_surf_data',Pselfp)>0; 
    
    
    % correction scheme that based on the original thickness
    % selfo self direction
    %   0     0      0
    %   1     0     -1 % not the last iteration!
    %   1     1      1
    %   0     1      1
    
    % at the beginning with have no old information, at the end we only want to reduce thickness 
    % here we use the original thickness value - test if new is better!
    if opt.mod~=1 || i==1 || final %corrsize <= opt.accuracy || 
      Twc = (T/2) .* corrsize .* selfw;
      Tpc = (T/2) .* corrsize .* selfp;
    else
      Twc = (T/2) .* corrsize .* ( (selfw | selfwo) - 2 .* ( selfwo & ~selfw ) );
      Tpc = (T/2) .* corrsize .* ( (selfp | selfpo) - 2 .* ( selfpo & ~selfp ) );
    end
    
    % correction in specified areas that also include a general
    % smoothness constrains of the cortical thickness
    Twc = single(spm_mesh_smooth(M,double(Twc) , sf)) * sf;  Tvxw = max(eps,Tvxw - Twc); 
    Tws = single(spm_mesh_smooth(M,double(Tvxw)  , sf/2));     Tvxw(Twc~=0) = Tws(Twc~=0);  clear Tws;
    Tpc = single(spm_mesh_smooth(M,double(Tpc) , sf)) * sf;  Tvxp = max(eps,Tvxp - Tpc);  
    Tps = single(spm_mesh_smooth(M,double(Tvxp)  , sf/2));     Tvxp(Tpc~=0) = Tps(Tpc~=0);  clear Tps;
    
    % save for next iteration
    selfwo = selfw;
    selfpo = selfp;
    
    % update thickness and surface
    VOC = S.vertices - N .* repmat(Tvxp,1,3); 
    VIC = S.vertices + N .* repmat(Tvxw,1,3); 
 
    % adaptive smoothing
    Tpc = Tpc./max(Tpc(:)); Tpc = single(spm_mesh_smooth(M,double(Tpc) , 1)); VOCSf = repmat(0.2 .* Tpc(Tpc>0),1,3); 
    Twc = Twc./max(Twc(:)); Twc = single(spm_mesh_smooth(M,double(Twc) , 1)); VICSf = repmat(0.5 .* Twc(Twc>0),1,3);

    % extra thickenss smoothing
    Sa   = cat_surf_fun('area',S); 
    Tsw  = repmat( max(0,min(1,max( (Tn-3)/6 , (Tn ./ (Sa * 4) - 3) / 6) )) ,1,3); clear Sa;
    VOCS = smoothsurf(VOC,2); VOC = VOC.*(1-Tsw) + Tsw.*VOCS;
    VICS = smoothsurf(VIC,2); VIC = VIC.*(1-Tsw) + Tsw.*VICS;
    clear VOCS VICS Tsw;
    
    for is=1
      VOCS = smoothsurf(VOC,1); VOC(Tpc>0,:) = VOC(Tpc>0,:).*(1-VOCSf) + VOCSf.*VOCS(Tpc>0,:);
      VICS = smoothsurf(VIC,1); VIC(Twc>0,:) = VIC(Twc>0,:).*(1-VICSf) + VICSf.*VICS(Twc>0,:);
      clear VOCS VICS Tsw;
    end
    
    % update thickness and central surface position
    Tn  = sum( (VIC - VOC).^2 , 2) .^ 0.5;
    Tvxw  = Tn/2; 
    Tvxp  = Tn/2; 
    S.vertices = mean(cat(3,VIC,VOC),3); 

    % update normals (only for smaller corrections otherwise this will introduce further errors!)
    if corrsize<0.25
      N   = spm_mesh_normals(S);   
    end
    
    % cleanup
    delete(Pwhite);
    delete(Ppial);
    delete(Pselfw);
    delete(Pselfp);
    
    % iteration
    if corrsize <= opt.accuracy, final = final + 1; end
    if opt.verb
      cat_io_cprintf('g5',sprintf( ...
        '    Step %2d (SS=%02.0f%%%%, SI=%5.2f%%%%, T=%4.2f%s%4.2f)\n',...
        i,corrsize*100, (sum(selfw>0)/2 + sum(selfp>0)/2) / numel(selfw) * 100,...
        mean(Tn),native2unicode(177, 'latin1'),std(Tn)));  
      %fprintf('/sprintf('%s',repmat('\b',1,73*2)));
    end
    if  (sum(selfw>0)/2 + sum(selfp>0)/2) / numel(selfw) * 100  <  opt.accuracy % if changes are below a specified relative level
      if final < 4 && (sum(selfw>0)/2 + sum(selfp>0)/2)>0                       % do some additional iterations if required
        final = final + 1; 
      else
        break
      end
    end
  end
  
  fprintf('\n');
  % final settings: back to world thickness in mm 
  if flipped, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; S.mati(7) = - S.mati(7); end
end

function flipped = cat_surf_checkNormalDir(S)
% estimate the orientation of the surface by adding the normals and
% testing which directions makes the surface great again :D
  N       = spm_mesh_normals(S); 
  Snin    = mean( abs(S.vertices(:) - N(:) ) ) > mean( abs(S.vertices(:) + N(:) ) ); 
  flipped = ~Snin;
end

function V = cat_surf_smooth(M,V,s,mode)
  if ~exist('s','var'), s = 1; end
  if ~exist('mode','var'), mode = 0; end

  smoothsurf = @(V,s) [ ...         % simple surface smoothing 
    spm_mesh_smooth(M,V(:,1),s) , ...
    spm_mesh_smooth(M,V(:,2),s) , ...
    spm_mesh_smooth(M,V(:,3),s) ];
  
  if isa(V,'single')
    singleV = 1;
    V = double(V); 
  else
    singleV = 0;
  end
  
  if mode == 0
    V = smoothsurf(V,s);
  else
    VS = smoothsurf(V,s);
      
    VD = sum( (VS - V).^2 , 2).^0.5;

    if mode<0
      V(VD<mode) = VS(VD<mode); 
    else
      V(VD>mode) = VS(VD>mode); 
    end  
  end
  
  if singleV
    single(V);
  end
end

function [S,Tn] = cat_surf_collision_correction_pbt(S,T,Y,Ypp,Yl4,opt)
% _________________________________________________________________________
% This function utilize the percentage position map Ypp map of the 
% projection-based thickness to detect and correct self-intersections (SI)
% of a central surface S (in voxel-space) and thickness T (in mm).
%
% [S,Tn] = cat_surf_collision_correction_pbt(S,T,Y,Ypp,Yl4,opt)
% 
%  S          .. central surface in voxel space 
%  T          .. thickness in world space (mm)
% ####
%  Y          .. intensity normalized image (resolution ????)
% ####
%  Ypp        .. percentage position map in voxel space
%  opt        .. structure with further options
%  .verb      .. print progress
%  .redterm   .. size of the reduction interval (default 2 = half interval) 
%  .accurarcy .. stop criteria (default 1/2^4 for about 5 iterations)
%  .iterfull  .. maximum number of iterations that is automatically defined 
%                by the accuracy 
%  .optimize  .. use position/intensity values for optimization (def = 1) 
%  .optmod    .. optimization mode (0 - Ypp, 1 - Ym, 2 - Ypp & Ym; def = 0)
%  .interBB   .. structure with resolution information 
%   .BB       .. not used here
%   .interpV  .. voxel-space resolution of S 
% _________________________________________________________________________
% Robert Dahnke 201910

  if ~exist('opt','var'), opt = struct(); end

  % defaults
  def.redterm  = 2; 
  def.verb     = 1;       % display information 
  def.debug    = 2;       % save results 
  def.relvert  = 0.01;    % percentual error of vertices with self-intersections (1.00 = 1%)  
  def.accuracy = 1/2^3;  % smallest used step-size (the value should not be to small   ... 0.00625
                          % because there should also be a gab between the structures and a reasonalbe running time) 
  def.optimize = 1;       % use position/intensity values for optimization  
  def.optmod   = 2;       % optimization mode (0 - Ypp, 1 - Ym, 2 - Ypp & Ym)
  def.mat      = [];
  %def.interpBB = struct('BB',[],'interpV',1);
  opt          = cat_io_checkinopt(opt,def); 
  opt.iteropt  = find( 1 ./ opt.redterm.^(1:100)  < opt.accuracy , 1); 
  opt.iterfull = round( opt.iteropt * (1 + opt.optimize) ); 
  
  sf           = round( sqrt( size(S.faces,1) / 50000) ); %2;  % ### empirical value 
  M            = spm_mesh_smooth(S);                           % for spm_smoothing matrix



  % inner and outer thickness seen from the central surface
  Tw = T / 2; 
  Tp = T / 2;
  
  % detection and correction for flipped faces to have always the same normal direction
  flipped = cat_surf_checkNormalDir(S); %,T,Y,opt.interpBB);
  if flipped, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; S.mati(7) = - S.mati(7); end

  smoothsurf = @(V,s) [ ...         % simple surface smoothing 
    spm_mesh_smooth(M,double(V(:,1)),s) , ...
    spm_mesh_smooth(M,double(V(:,2)),s) , ...
    spm_mesh_smooth(M,double(V(:,3)),s) ];
  
  if opt.optimize && opt.optmod
    % use smoothed layer4 values only once to avoid re-estimation 
    Yl4 = single(spm_mesh_smooth(M,double(Yl4),sf/4 * 50));
  end
  
  Tn  = T; meshnorm = 1; 
  if meshnorm
    N   = spm_mesh_normals(S);
    
  else
    %%
    V   = cat_surf_mat(S.vertices,opt.mat,1); 
    N   = isonormals(1-Ypp,V); %[V(:,2),V(:,1),V(:,3)]);
    % normalize
    Ns  = sum(N.^2,2).^.5;
    N   = N ./ repmat(Ns,1,3); %clear Ns; 
    N   = cat_surf_mat(N,opt.mat .* [[ones(3,3) zeros(3,1)];zeros(1,3) 1],0); 
    
    %%
    %NS  = spm_mesh_normals(S);
   % N(Ns<0.2) = NS(Ns<0.2);  
    
  end
  N   = smoothsurf(N,sf * 2); 
  Ns  = sum(N.^2,2).^.5;
  N   = N ./ repmat(Ns,1,3); %clear Ns; 
  
  i   = 0; final = 0; 
  SIO = 1; 
 
  if opt.verb, fprintf('\n'); end
  while i < opt.iterfull * 2
    i = i + 1;
   
    % In theory the search/correction size would be half each time (0.5 0.25 0.125 ... 2^i)
    % but due to the smoothing a slower sloop is used (1.8^i). 
    % The changes are limited by opt.accuracy for anatomical (min sulcus-width) and runtime reasons.
    % If optimization is used than the test starts again but with higher accuracy ( 0.125 ). 
    corrsize  = max( opt.accuracy , 1 ./ opt.redterm.^(i+2) );
   
    % outer and inner boundary
    VO = S.vertices - N .* repmat(Tp,1,3);
    VI = S.vertices + N .* repmat(Tw,1,3); 

    % estimate the gradient in the Ypp map of each cortical surface 
    Vg  = cat_surf_volgrad(S.vertices,N,Ypp,opt.mat);
    VOg = cat_surf_volgrad(VO,N,Ypp,opt.mat);
    VIg = cat_surf_volgrad(VI,N,Ypp,opt.mat);

    % if the gradients point in totally other directions that we had 
    % crossed an anatomical border and have an self intersection
    opt.alphaYpp = 20; % max(5 , min(30, 15 / T ))';
    selfp = ( spm_mesh_smooth( M , double( cat_surf_edgeangle( Vg , VOg )) ,1 ) > opt.alphaYpp ); % | (T>2 & angle( VIg , VOg ) > opt.alphaYpp*2 ); 
    %selfc = ( spm_mesh_smooth( M , double( angle( Vg , N   )) ,1 ) > opt.alphaYpp ); % | (T>2 & angle( VIg , VOg ) > opt.alphaYpp*2 ); 
    selfw = ( spm_mesh_smooth( M , double( cat_surf_edgeangle( Vg , VIg )) ,1 ) > opt.alphaYpp );
    
    % correction scheme that based on the original thickness
    % selfo self direction
    %   0     0      0
    %   1     0     -1 % not the last iteration!
    %   1     1      1
    %   0     1      1
    % 
    if i==1 || final % ||  corrsize <= opt.accuracy || i > opt.iterfull;  % at the beginning with have no old information, at the end we only want to reduce thickness 
      Twc = T/2 .* corrsize .* selfw;
      Tpc = T/2 .* corrsize .* selfp;
    else
      Twc = T/2 .* corrsize .* ( (selfw | selfwo) - 2 .* ( selfwo & ~selfw ) );
      Tpc = T/2 .* corrsize .* ( (selfp | selfpo) - 2 .* ( selfpo & ~selfp ) );
    end
    
    % optimization by local intensities
    if opt.optimize && i<opt.iteropt*3 %&& ~(corrsize <= opt.accuracy || final || i > opt.iterfull)  %
      
      if opt.optmod == 3
        % intensity model based on the Ym
        YI = cat_surf_isocolors2(Y,VI,opt.mat);
        YO = cat_surf_isocolors2(Y,VO,opt.mat);

        if opt.verb, fprintf('  YIC:%5.2f%s%0.2f, YOC:%5.2f%s%0.2f',mean(YI),native2unicode(177, 'latin1'),std(YI),mean(YO),native2unicode(177, 'latin1'),std(YO)); end 
      
        % to fast corrections can jump gyri/sulci 
        GWth = 3 * 0.5 + 0.5 * Yl4; 
        CGth = 1 * 0.7 + 0.3 * Yl4;

        % correction value 
        YI = max(-1, min(1, ( GWth - YI ) * max(0.02,0.8 .* max(0.5,1 - i / opt.iteropt) ) ));
        YO = max(-1, min(1, ( YO - CGth ) * max(0.02,0.8 .* max(0.5,1 - i / opt.iteropt) ) ));
      elseif opt.optmod == 2
        % intensity model based on the Ym
        YI = cat_surf_isocolors2(Ypp,VI,opt.mat); YI = cat_surf_isocolors2(Y,VI,opt.mat) + YI - (YI>eps) .* 0.95; 
        YO = cat_surf_isocolors2(Ypp,VO,opt.mat); YO = cat_surf_isocolors2(Y,VO,opt.mat) + YO - (YO>eps) .* 0.05;  

        if opt.verb, fprintf('  YIC:%5.2f%s%0.2f, YOC:%5.2f%s%0.2f',mean(YI),native2unicode(177, 'latin1'),std(YI),mean(YO),native2unicode(177, 'latin1'),std(YO)); end 
      
        % to fast corrections can jump gyri/sulci 
        GWth = 3 * 0.5 + 0.5 * Yl4; 
        CGth = 1 * 0.7 + 0.3 * Yl4;

        % correction value 
        YI = max(-1, min(1, ( GWth - YI ) * max(0.02,0.4 .* max(0.5,1 - i / opt.iteropt / 2) ) ));
        YO = max(-1, min(1, ( YO - CGth ) * max(0.02,0.4 .* max(0.5,1 - i / opt.iteropt / 2) ) ));
      else
        GWth = 0.98;
        CGth = 0.02;
        
        % position model based on the Ypp
        YI = cat_surf_isocolors2(Ypp,VI,opt.mat); 
        YO = cat_surf_isocolors2(Ypp,VO,opt.mat);  

        if 1+opt.verb>1, fprintf('  YIC:%5.2f%s%0.2f, YOC:%5.2f%s%0.2f',mean(YI),native2unicode(177, 'latin1'),std(YI),mean(YO),native2unicode(177, 'latin1'),std(YO)); end 

        % correction value 0.01, 0.02
        YI = max(-1, min(1, ( GWth - YI ) * max(0.02,0.2 .* (1 - i / opt.iteropt / 2 )) ));
        YO = max(-1, min(1, ( YO - CGth ) * max(0.02,0.2 .* (1 - i / opt.iteropt / 2 )) ));
      end  
      
      % test new inner surface position
      VIC   = S.vertices + N .* repmat( Tw - (Twc - YI) ,1,3);    % inner surface
      VOC   = S.vertices - N .* repmat( Tp - (Tpc - YO) ,1,3);    % outer surface 
      
      if opt.optmod == 3
        YIC   = cat_surf_isocolors2(Y,VIC,opt.mat);
        YOC   = cat_surf_isocolors2(Y,VOC,opt.mat);  
      elseif opt.optmod == 2   
        YIC = cat_surf_isocolors2(Ypp,VIC,opt.mat); YIC = cat_surf_isocolors2(Y,VIC,opt.mat) + YIC - (YIC>eps) .* 0.95; 
        YOC = cat_surf_isocolors2(Ypp,VOC,opt.mat); YOC = cat_surf_isocolors2(Y,VOC,opt.mat) + YOC - (YOC>eps) .* 0.05;  
      else
        YIC   = cat_surf_isocolors2(Ypp,VIC,opt.mat);
        YOC   = cat_surf_isocolors2(Ypp,VOC,opt.mat);  
      end
      
      YppI  = cat_surf_isocolors2(Ypp,VI,opt.mat);
      YppIC = cat_surf_isocolors2(Ypp,VIC,opt.mat);  
      VIg   = cat_surf_volgrad(VIC,N,Ypp,opt.mat);
      YI    = YI .* ( Twc==0 & ...
        abs( GWth - YI ) > abs( GWth - YIC )  & ...
        cat_surf_edgeangle( Vg , VIg ) < opt.alphaYpp & ...
        ( YppIC > YppI | YppI>0.98) & (YppIC<0.98 | YIC<2.9) & YIC<2.95 ); % 0.98 & 2.9
      clear YppIC YIC YppI VIg; 
      
      % test new outer surface position
      YppO  = cat_surf_isocolors2(Ypp,VO,opt.mat);  
      YppOC = cat_surf_isocolors2(Ypp,VOC,opt.mat);  
      VOg   = cat_surf_volgrad(VOC,N,Ypp,opt.mat);
      if opt.optmod 
        YO  = YO .* YppOC .* ( Tpc==0 & ...
          abs( YO - CGth ) > abs( YOC - CGth ) & ...
          cat_surf_edgeangle( Vg , VOg ) < opt.alphaYpp & ...
          YppOC < YppO & YppOC>0.02 & YOC>1.50); % 0.01 & 1.5 
      else
        YO  = YO .* (0.5+0.5*YppOC) .* ( Tpc==0 & ...
          cat_surf_edgeangle( Vg , VOg ) < opt.alphaYpp & ...
          YppOC>0.02 );
      end
      clear YppOC YOC YppO VOg; 
      
      if opt.verb>1, fprintf('  YIC:%5.2f%s%0.2f, YOC:%5.2f%s%0.2f',mean(YI),native2unicode(177, 'latin1'),std(YI),mean(YO),native2unicode(177, 'latin1'),std(YO)); end
      % YO = YO ./ max(0.8,TN/mean(TN(:))); % less correction if thicker than the average - this does not work!  
      
      % smooth data 
      YO = single(spm_mesh_smooth(M,double(YO), sf )) * sf/2; 
      YI = single(spm_mesh_smooth(M,double(YI), sf )) * sf/2; 
    
      % add to other correction map
      Tpc = Tpc - YO; clear YO; 
      Twc = Twc - YI; clear YI; 
    end

    % das hilft hier tatsaechlich
    sulciGyriWidth = 0.05; % war 0.2 aber sulci zu breit
    Twc(Twc>0) = Twc(Twc>0) + sulciGyriWidth * 2;
    Tpc(Tpc>0) = Tpc(Tpc>0) + sulciGyriWidth;
    
    % correction in specified areas that also include a general
    % smoothness constrains of the cortical thickness
    % correction in specified areas that also include a general
    % smoothness constrains of the cortical thickness
    Twc = single(spm_mesh_smooth(M,double(Twc) , sf)) * sf;  Tw = max(eps,Tw - Twc);  
    Tws = single(spm_mesh_smooth(M,double(Tw)  , sf/2));     Tw(Twc~=0) = Tws(Twc~=0);  clear Tws;
    Tpc = single(spm_mesh_smooth(M,double(Tpc) , sf)) * sf;  Tp = max(eps,Tp - Tpc);  
    Tps = single(spm_mesh_smooth(M,double(Tp)  , sf/2));     Tp(Tpc~=0) = Tps(Tpc~=0);  clear Tps;
        
    selfwo = selfw;
    selfpo = selfp;
    
    % update thickness and surface
    VOC = S.vertices - N .* repmat(Tp,1,3); 
    VIC = S.vertices + N .* repmat(Tw,1,3); 
    
    % edge flip
    if 1 % this has a large effect
      sm   = sf; % smoothness - 1 is not enougth
      fa   = 60; % error angle (worst case is 180 )
      E    = spm_mesh_edges(S);
      V    = S.vertices;

      VCOA = cat_surf_edgeangle( V(E(:,1),:) - V(E(:,2),:) , VOC(E(:,1),:) - VOC(E(:,2),:) ); 
      VCOC = VOC(E(VCOA>fa,1),:)/2 + VOC(E(VCOA>fa,2),:)/2; 
      VOC(E(VCOA>fa,1),:) = VCOC; VOC(E(VCOA>fa,2),:) = VCOC; clear VCOC;
      VTPO = T*0; VTPO(E(VCOA>fa,1)) = 1; VTPO(E(VCOA>fa,2)) = 1; clear VCOA;
      VTPO = repmat(min(1,max(0,single(spm_mesh_smooth(M,double(VTPO) , 1*sm )) * 1*sm)),1,3);
      VOCS = cat_surf_smooth(M,VOC,sm); VOC = VOC.*(1-VTPO) + VTPO.*VOCS;
      clear VCOC VTPM;

      
      VCIA = cat_surf_edgeangle( V(E(:,1),:) - V(E(:,2),:) , VIC(E(:,1),:) - VIC(E(:,2),:) ); 
      VCIC = VIC(E(VCIA>fa,1),:)/2 + VIC(E(VCIA>fa,2),:)/2; 
      VIC(E(VCIA>fa,1),:) = VCIC; VIC(E(VCIA>fa,2),:) = VCIC; clear VCIC;
      VTPI = T*0; VTPI(E(VCIA>fa,1)) = 1; VTPI(E(VCIA>fa,2)) = 1; clear VCIA;
      VTPI = repmat(min(1,max(0,single(spm_mesh_smooth(M,double(VTPI) , 2*sm))) * 2*sm),1,3);
      VICS = cat_surf_smooth(M,VIC,sm); VIC = VIC.*(1-VTPI) + VTPI.*VICS;
      clear VCIC VTPM;
      clear E V;
    else
      VTPO = zeros(size(T),'single'); 
      VTPI = zeros(size(T),'single');
    end
    
   % adaptive smoothing
  %  Tpc = abs(Tpc); Twc = abs(Twc);
    Tf  = max(1,min(5,Tn)); inorm = ( opt.iteropt*3 - i ) / (opt.iteropt*3); 
    Tpc = Tpc./max(Tpc(:)); Tpc = single(spm_mesh_smooth(M,double(Tpc) , 1 )) * 1/2 .* Tf; VOCSf = max(0,min(1,repmat(0.2 .* inorm .* abs(Tpc(Tpc>0 | VTPO(:,1)>0)),1,3))); 
    Twc = Twc./max(Twc(:)); Twc = single(spm_mesh_smooth(M,double(Twc) , 1 )) * 2/2 .* Tf; VICSf = max(0,min(1,repmat(0.8 .* inorm .* abs(Twc(Twc>0 | VTPI(:,1)>0)),1,3)));

    for is = 1
      VOCS = smoothsurf(VOC,1); VOC(Tpc>0 | VTPO(:,1)>0,:) = VOC(Tpc>0 | VTPO(:,1)>0,:).*(1-VOCSf) + VOCSf.*VOCS(Tpc>0 | VTPO(:,1)>0,:);
      VICS = smoothsurf(VIC,1); VIC(Twc>0 | VTPI(:,1)>0,:) = VIC(Twc>0 | VTPI(:,1)>0,:).*(1-VICSf) + VICSf.*VICS(Twc>0 | VTPI(:,1)>0,:);
    end
    clear VTPI VTPO

    % extra thickenss smoothing
    if opt.optimize
      Sa   = cat_surf_fun('area',S); 
      Tsw  = repmat( max(0,min(1,max( (Tn-3)/6 , (Tn ./ (Sa * 4) - 3) / 6) )) ,1,3) * 0.5; % * min(1,0.02 / opt.iterfull); 
      VOCS = smoothsurf(VOC,sf); VOC = VOC.*(1-Tsw) + Tsw.*VOCS;
      VICS = smoothsurf(VIC,sf); VIC = VIC.*(1-Tsw) + Tsw.*VICS;
    end
    
    % full smooting 
    if opt.optimize  
      VOCSf = repmat(Tf * 0.02 / opt.iterfull,1,3);  VOCS = smoothsurf(VOC,sf/2); VOC = VOC.*(1-VOCSf) + VOCSf.*VOCS; clear VOCSf;
      VICSf = repmat(Tf * 0.05 / opt.iterfull,1,3);  VICS = smoothsurf(VIC,sf/2); VIC = VIC.*(1-VICSf) + VICSf.*VICS; clear VICSf;
    end
    clear VICS VOCS 
    
    % outlier smoothing 
    VOC = cat_surf_smooth(M,VOC,sf,1);
    VIC = cat_surf_smooth(M,VIC,sf,1);
    
    % remove outlier
    VOCc = spm_mesh_curvature( struct('vertices',VOC,'faces',S.faces) ); 
    VICc = spm_mesh_curvature( struct('vertices',VIC,'faces',S.faces) ); 
    
    if mod(i,5)==0 || final
      Tn   = max(0.01,sum( (VIC - VOC).^2 , 2) .^ 0.5 );
      Tns  = spm_mesh_smooth(M,Tn,sf/2);
      Tnss = spm_mesh_smooth(M,Tn,sf/2 * 5);
      Tnm  = abs( Tns - T )>0.5/(final+1) & (VOCc>60 & VICc>60); clear Tns;  
      Tn( Tnm ) = Tnss( Tnm ); clear Tnss; 
      VIC( Tnm ,:)  = S.vertices( Tnm ,:) + N( Tnm ,:) .* repmat( Tn( Tnm )/2 ,1,3);    % inner surface
      VOC( Tnm ,:)  = S.vertices( Tnm ,:) - N( Tnm ,:) .* repmat( Tn( Tnm )/2 ,1,3);    % outer surface 
      clear Tnm; 
    end
    
    % final new thickness
    Tn  = max(0.01,sum( (VIC - VOC).^2 , 2) .^ 0.5 ); % -  0.0002 * ((i<opt.iterfull*2) + 2*(opt.optimize | i<opt.iterfull)); %single(spm_mesh_smooth(M, double(,1))); %* (Tpc>0 | Twc>0)
    Tw  = Tn/2; 
    Tp  = Tn/2; 
    S.vertices = mean(cat(3,VIC,VOC),3); 

    % update normals
    if corrsize<=0.25
      if meshnorm
        N   = spm_mesh_normals(S);
      else
        V   = cat_surf_mat(S.vertices,opt.mat,1); 
        N   = isonormals(1-Ypp,V); %[V(:,2),V(:,1),V(:,3)]);
        % normalize
        Ns  = sum(N.^2,2).^.5;
        N   = N ./ repmat(Ns,1,3); %clear Ns; 
        N   = cat_surf_mat(N,opt.mat .* [[ones(3,3) zeros(3,1)];zeros(1,3) 1],0); 

        %NS  = spm_mesh_normals(S);
        %N(Ns<0.2) = NS(Ns<0.2);  
      end
      N   = smoothsurf(N,sf * 2); 
      Ns  = sum(N.^2,2).^.5;
      N   = N ./ repmat(Ns,1,3); %clear Ns; 
    end
    
    SI = (sum(selfw>0)/2 + sum(selfp>0)/2) / numel(selfw) * 100; 
    % iteration
    if corrsize <= opt.accuracy, final = final + 1; end
    if opt.verb
      cat_io_cprintf('g5',sprintf( ...
        '    Step %2d (SS=%02.0f%%%%, SI=%5.2f%%%%, T=%4.2f%s%4.2f)',...
        i,corrsize*100, SI, mean(Tn),native2unicode(177, 'latin1'),std(Tn)));  
      fprintf('\n')
      %fprintf(sprintf('%s',repmat('\b',1,73*2)));
    end
    if  ((sum(selfw>0)/2 + sum(selfp>0)/2) / numel(selfw) * 100  <  opt.accuracy) || (i>opt.iterfull && SI<0.1 && abs(SI - SIO)<0.005) % if changes are below a specified relative level
      if final < 4 && (sum(selfw>0)/2 + sum(selfp>0)/2)>0 && i>opt.iterfull  && abs(SI - SIO)>0.005                      % do some additional iterations if required
        final = final + 1; 
      else
        break
      end
    end
    SIO = SI;
  end
  
  fprintf('\n');
  if flipped, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; S.mati(7) = - S.mati(7); end
  
end 
function [Sve,Svde] = cat_surf_disterr(S,T,mat,Yp0,tol)

  % render volume
  % [Yp,Yt,vmat1,vmat1i] = cat_surf_surf2vol(S,Y,T,type,opt)
  Yp0S = cat_surf_surf2vol(S,Yp0,T,'seg',struct('mat',mat));
  Ycs  = cat_surf_surf2vol(S,Yp0,T,'p0',struct('mat',mat));
  
  Ydiff = Yp0S - Yp0;
  
  % Delaunay triangulation 
  D  = delaunayn(S);
  
  % select relevant points for projection
  PI = find( abs(Ydiff) > tol & Ycs>0 ); 
  PO = find( abs(Ydiff) > tol & Ycs<0 );
  
  % project points
  [DI,II] = dnear(D,PI); 
  [DO,IO] = dnear(D,PI); 
  
  % remove to close points
  
  
  % mapping
  for di = 1:numel(DI)
    Sve(  II(di) )  = Sve(  II(di) ) + Ydiff( PI(di) );             % just the volume 
    Svde( II(di) )  = Svde( II(di) ) + DI(di) * Ydiff( PI(di) );    % weighted by distance 
  end
end
function [SN,TN,E] = cat_surf_collision_correction(S,T,Y,Ypp,Yl4,opt) 
% _________________________________________________________________________
% Delaunay based collision detection:
% 1) Correction for local curvature
% 2) Creation of a Delaunay graph 
% 3) Differentiation of intra (edges between surface points) and inter
%    surface edges (e.g., edges between to opposite gyri or sulci) and 
%    removal of intra-surface edges. 
% 4) Use of the inter-surfaces edges to detect collision by normal 
%    transformations of the half thickness to obtain the inner and outer 
%    surfaces.
% 5) Further correction of possible flips by normal transformation
%
% This is a prototype that allows correction of the worst things but not 
% all collisions. 
%
%   [SN,TN,E] = cat_surf_collision_correction(S,T,Y[,debug,E,Pcs])
%
%   SN      .. new surface
%   TN      .. new thickness
%   S       .. original surface
%   T       .. original thickness
%   Y       .. segmentation map or intensity normalized images 
%              for intra/inter surface edge definition
%   Yl4     .. layer4 intensity surface map
%   opt     .. parameter structure
%    .debug .. option to write the un- and corrected cortical surfaces in a
%              subdirectory
%    .verb  ..
%    .E     .. Delaunay edge map (from previous run or empty matrix as input)
%    .Pcs   .. central surface file name to write debugging files
%    ... experimental settings
%    .smoothCSinput  ..
%    .model          ..
%    .PVEcorr        ..
%    .slowdown       ..
%    .
% _________________________________________________________________________
% Robert Dahnke 201909


% -------------------------------------------------------------------------
% Todo:
% - full support of parameters & recall 
% - support of different correction models 
%   (e.g. only collision vs. intensity correction)
% - full documentation and detailed comments
% - helping boundary surfaces? 
%   > partially implemented 
%   > very slow +20-60s for each just for creation  
% - stable subset/list (also as internal/external error measure
% - use of gradients and divergence rather simple intensity information
% - use of Ypp
% - optimization of the layer 4 (layer concept)
% - face-flipping correction that can not be handled by Delaunay because of
%   its neighbor limits
% - improved evaluation concept
% - improved validation concept
% - fast mapping c-function for edge to surf that combine multiple values
%   given by an index map by different functions (mean,min,max,std)
% - surface filter sub-function to remove outlier
% - triangle height rather than edge distance (or combination)
% -------------------------------------------------------------------------

  if ~exist('opt','var'), opt = struct(); end 

  % default variables 
  def.Pcs               = '';     % filename to write debugging output data
  def.debug             = 1;      % debugging output vs. memory optimization 
  def.verb              = 1;      % display debugging information 
  def.E                 = [];     % inter-surface Delaunay edges (of a previous run) 
  def.boundarySurfaces  = 0;      % use inner and outer boundary surface to improve the Delaunay graph 
  def.smoothCSinput     = 0;      % smooth the input CS for more stable Delaunay triangulation in case of locally oversampled surfaces 
  def.PVEcorr           = 1;      % correction of PVE values for 2 boundaries in on voxel (experimental)
  def.slowdown          = 1;      % slowdown may stabilize the process over the iterations   
  def.model             = 2;      % 0 - only collcorr, 1 - only intopt, 2 - both  
  def.vx_vol            = 1; 
  opt                   = cat_io_checkinopt(opt,def); clear def;
  def.write             = opt.debug & ~isempty(opt.Pcs);
  opt                   = cat_io_checkinopt(opt,def); 
  
  
  % helping smoothing functions for data and surfaces
  M   = spm_mesh_smooth(S);         % for spm_smoothing matrix
  rms = @(x) mean( x.^2 ) .^ 0.5;   % for error handling of mad vertices
  smoothsurf = @(V,s) [ ...         % simple surface smoothing 
    spm_mesh_smooth(M,double(V(:,1)),s) , ...
    spm_mesh_smooth(M,double(V(:,2)),s) , ...
    spm_mesh_smooth(M,double(V(:,3)),s) ];
  
  
  % detection and correction for flipped faces to have always the same normal direction
  lim = 1:round(size(S.vertices,1)/1000):size(S.vertices,1); 
  N   = spm_mesh_normals(S);   
  VOl = S.vertices(lim,:) - N(lim,:) .* repmat(T(lim)/2,1,3); 
  VIl = S.vertices(lim,:) + N(lim,:) .* repmat(T(lim)/2,1,3); 
  YOl = cat_surf_isocolors2(Y,VOl); 
  YIl = cat_surf_isocolors2(Y,VIl); 
  flipped = mean(YOl) > mean(YIl); 
  clear N VOl VIl YOl YIl lim; 
  if flipped, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; S.mati(7) = - S.mati(7); end
  
  
  % larger surface need more smoothing to avoid triangulation problems 
  sf = round( sqrt( size(S.faces,1) / 50000) );  % ### empirical value 
  if max(Y(:))<1.5, Y = Y.*2+1; else, Y = max(1,Y); end              
  if opt.debug, fprintf('\n'); end
  if opt.write, cat_surf_saveICO(S,T,mat,Pcs,sprintf('pre_collcorr_%0.0fk',round( size(S.faces,1)/1000 / 10) * 10 ),0); end
  stime = cat_io_cmd(sprintf('    Delaunay triangulation of %d vertices (sf=%d):',size(S.vertices,1),sf),'g5','',opt.debug); 

  
  

  
  
  %% Creation of the inter surface edges based on a Delaunay graph 
  %  ----------------------------------------------------------------------
  %  There is a short cut to apply further iterations without processing
  %  the graph again. 
  %  ----------------------------------------------------------------------
  if isfield(opt,'E') && isempty(opt.E)
    
    % Early versions used a smoothed surface to reduce problems due to
    % artifacts. However, the improved input meshed (createCS2 pipeline)
    % do not need smoothing and surface smoothing can not be combined 
    % with helping surfaces (inner and outer surface points). 
    VS = double(S.vertices); 
    if opt.smoothCSinput
      % Surface smoothing as loop to correct for outlier due to incorrect surfaces.
      % Using the smoothing directly create some extrem large spikes - don't know why (RD20190912).
      % The smoothing is not required in newer version (RD20190922)
      for i = 1:opt.smoothCSinput                           
        VSi = smoothsurf(VS,2); 
        VM  = rms(VSi - VS)<2; 
        VS(VM,:) = VSi(VM,:); 
      end
      VS = smoothsurf(VS,1);
    end
    if ~opt.debug, clear VSi VM; end

    
    % helping boundary surface - (uncorrected) inner or/and outer surface 
    % this was slow (20-60 seconds) and did not work so simple/fast ... need further work
    % maybe just use low resolution surfaces (1 mm)
    if opt.boundarySurfaces
      if opt.smoothCSinput, error('Initial surface smoothing "smoothCSinput" can not be combined with "boundarySurfaces".\n'), end
      VB = S.vertices;                             
      if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
        UOS = isosurface(Y,1.5); 
        VS  = [ VS ; UOS.vertices ]; 
        VB  = [ VB ; UOS.vertices ]; 
        if ~opt.debug, clear UOS; end 
      end
      if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
        UIS = isosurface(Y,1.5); 
        VS  = [ VS ; UIS.vertices ]; 
        VB  = [ VB ; UIS.vertices ]; 
        if ~opt.debug, clear UIS; end 
      end
    end
    
    
    % Delaunay graph
    D  = single(delaunayn( VS )); 
    if ~opt.debug, clear VS; end            

    % decompose delaunay graph into its edges
    E  = uint32(cat_surf_edges(D));   
    nE = size(E,1); 
    if ~opt.debug, clear D; end
    if opt.debug, cat_io_cprintf('g5',sprintf('%5.0fs\n',etime(clock,stime))); end


    % separate helping boundary surfaces
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      Etmep = sum( E>( size(S.vertices,1) + any(opt.boundarySurfaces==[1,3]).*size(UOS.vertices,1) ) , 2 )>0; 
      EUIS  = E( Etmep , : ); 
      EUIS( sum( EUIS>numel(S.vertices) , 2 )~=1, :) = []; % remove all edges that are not between the CS and the UIS
      EUIS  = sort(EUIS,2); 
      E( Etmep , : ) = [];
    end
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      Etmep = sum( E>( size(S.vertices,1) ) , 2 )>0; 
      EUOS  = E( Etmep , : ); 
      EUOS( sum( EUOS>numel(S.vertices) , 2 )~=1, :) = []; % remove all edges that are not between the CS and the UOS
      EUOS  = sort(EUOS,2); 
      E( Etmep , : ) = []; 
    end
    

    %% Remove intra-surface edges 
    %  --------------------------------------------------------------------
    %  If we remove too much then the correction will not work.
    %  If we do not remove enough then it will add sulci in regions without sulci

    V  = S.vertices;

    % remove edge that we know from the surface - super save
    stime = clock; 
    EF = uint32(cat_surf_edges(S.faces));           
    E  = setdiff(E,EF,'rows'); clear EF; 
    
    % remove edges between neigbors of each point - relative save
    % get neighbor matrix
    [NE,MED] = spm_mesh_neighbours(M); 
    nNE = size(NE,1);

    % extra element that link on itself and replace the 0 in NE that does not allow matrix indexing 
    NIL = nNE + 1; 
    NE(NIL,:)  = NIL * ones(1,size(NE,2)); NE(NE==0) = nNE+1;
    
    % further levels
    % use higher levels only for large surfaces (use sqrt to compensate area grow factor)
    nlevel = 2; % max(2,round( 1 + sqrt( ceil( size(S.faces,1) / 300000 ) ))); 
    for nli = 1:nlevel % nice idear but not working yet
      for ni = 1:size(NE,2)
        NEN = sum( NE == repmat( NE(NE(:,ni),ni) , 1 , size(NE,2) ),2)>0; 
        NE  = [ NE  min( NIL , NE(:,ni)  + NIL*NEN) ];  %#ok<AGROW>
        MED = [ MED MED(:,ni) ];  %#ok<AGROW>
      end
      
      % sort entries
      [NE,NEsi] = sort(NE,2); MED = MED(NEsi); clear NEsi
      NILi = min([ size(NE,2) , find( sum(NE == NIL,1) >= size(NE,1)*0.5 , 1, 'first') - 1]);
      NE  = NE(:,1:NILi); 
      MED = MED(:,1:NILi);  
    end
    NE(NE==NIL) = 0; 

    % remove edges from the neigbor list
    for i=2:size(NE,2)
      E = setdiff(E,[NE(:,1) NE(:,i)],'rows'); 
      E = setdiff(E,[NE(:,i) NE(:,1)],'rows'); 
    end
    clear NE
    if opt.debug
      cat_io_cprintf('g5',sprintf('    remove edges by surface (l%d):%8d > %9d (%0.2f%%%%) %9.0fs',...
        nlevel,nE,size(E,1),size(E,1)./nE,etime(clock,stime)));
    else
      clear nE
    end

    if 1
        % remove edge by distance - this is not clear but it helps
        stime = clock; 
        LE  = sum( (V(E(:,1),:) - V(E(:,2),:)).^2 , 2) .^ 0.5; % length of edge
        DE  = min( LE > max(0.5,min(1,max(MED(:)))) , min( LE - T(E(:,1))*0.33 , LE - T(E(:,2))*0.33 ));
        NEd = abs(DE); clear LE DE MED;

        % remove by angle .. sum(NEa)./numel(NEa), figure, hist( S1alpha, -180:1:180)
        %   N(S)alpha  .. angle between the (smoothed) normals of each edge
        %                 (~0? = surface edge; ~180? between surface edge)       
        %   S[12]alpha .. angle between the edge and the first normal    
        %                 (~0?/~180? = surface edge; ~90? = between surface edge)
        % 
        N  = spm_mesh_normals(S);                 
        NS = N; for i=1:80*sf, NSS = smoothsurf(NS,1); NM = rms(NS - NSS)<0.5; NS(NM,:) = NSS(NM,:); end 
        Nalpha  = [cat_surf_edgeangle(NS(E(:,1),:), NS(E(:,2),:)), ...
                   cat_surf_edgeangle(NS(E(:,2),:), NS(E(:,1),:))]; clear NS
        SNalpha = [cat_surf_edgeangle(N(E(:,1),:),  V(E(:,1),:) - V(E(:,2),:)), ...
                   cat_surf_edgeangle(N(E(:,2),:),  V(E(:,2),:) - V(E(:,1),:))]; 
        NEna    = mean(Nalpha/180,2); clear Nalpha                       % figure, hist( NEna , 0:0.01:1);
        NEsa    = (abs(90  - SNalpha)/90  + abs(90  - SNalpha)/90)/2;    % figure, hist( NEsa , 0:0.01:1);
        clear SNalpha; 

        % remove by intensity given by the centroids of the edges
        VC  = cat_surf_centroid(V,E); 
        IC  = cat_surf_isocolors2(Y,VC); clear VC;
        % outer surface intensity
        VO  = V - N .* repmat(T/2,1,3); 
        VOC = cat_surf_centroid(VO,E); 
        IO  = cat_surf_isocolors2(Y,VOC); clear VOC VO; 
        % inner surface intensity
        VI  = V + N .* repmat(T/2,1,3) + 0.1; % GM/WM  
        VIC = cat_surf_centroid(VI,E); 
        II  = cat_surf_isocolors2(Y,VIC); clear VIC VI; 
        VI  = V + N .* repmat(T/2,1,3) + 0.5; % save WM  
        VIC = cat_surf_centroid(VI,E); 
        II  = max(II,cat_surf_isocolors2(Y,VIC)); clear VIC VI; % use max to get WM value 
        VI  = V + N .* repmat(T/2,1,3) + 1.0; % supersave WM  
        VIC = cat_surf_centroid(VI,E); 
        II  = max(II,cat_surf_isocolors2(Y,VIC)); clear VIC VI; % use max to get WM value 
        % combine all intensities 
        NEi = 1 - min(1,max(abs(diff([II IC IO],1,2)),[],2)); 
        %ET  = mean([II IC IO],2)>2.25; % edge classification 
        if ~opt.debug, clear II IC IO; end

        % combine all measures by product to remove many things
        % I also though about an adaptive threshold but it is not so easy ...
        NE = prod( [NEd NEi NEna*2 NEsa*2] ,2); % 1.75 % larger values > remove less
        NE = NE < .05; %05; %max(eps,mean(NE) - 1*std(NE)); % smaller values > remove less
        E (NE,:) = []; %if exist('ET','var'), ET(NE) = []; end
        if opt.debug
          cat_io_cprintf('g5',sprintf('\n    remove edges by intensity:            > %9d (%0.2f%%%%) %9.0fs',...
            size(E,1),size(E,1)./nE,etime(clock,stime))); stime = clock;
        else
          clear NE NEd NEi NEna NEsa
        end
    else
        N  = spm_mesh_normals(S);   
    end
    
    
    %fprintf('\nsf = %0.2f',sf);
  else
    N  = spm_mesh_normals(S);   
  end
  
 
  if opt.debug
    cat_io_cprintf('g5','\n    Prepare Optimization:'); stime = clock; 
  end
  
  %% updated measures
  SNalpha = [cat_surf_edgeangle(N(E(:,1),:),  V(E(:,1),:) - V(E(:,2),:)), ...
             cat_surf_edgeangle(N(E(:,2),:),  V(E(:,2),:) - V(E(:,1),:))]; 

  VC  = cat_surf_centroid(V,E); 
  IC  = isocolors(Y,VC); clear VC;

  OE  = min(1,(IC<2.15)); %(min(SNalpha,[],2)<90) + 
  IE  = min(1,(IC>2.15));%(max(SNalpha,[],2)>90) + 
  if 0
    %% map outer or inner edges
    %VLE = inf(size(T),'single');
    VLE = zeros(size(T),'single');
    for ni=1:size(E,1)
      VLE(E(ni,1)) = VLE(E(ni,1)) + IE(ni); 
    end
  end
  
  
  if 0
    % avoid PBT overestimation in gyri (well thickness is correct but
    % measures non-linear/non-orthogonal)
    TN = single(spm_mesh_smooth(M,double(T), sf * 20 )); 
    T  = min(T,TN); 
  end
  
  TN = T; SN = S; TCsum = 0; %#ok<NASGU>
  TCsumo = inf; TNold = inf;
  maxiter  = 5;  % main number of iterations 
  maxiter2 = 5;  % limit of adapting the mixing model
  
  Yl4 = single(spm_mesh_smooth(M,double(Yl4),sf/4 * 100));
  
  % I did not manage to use curvate ...
  %C   = spm_mesh_curvature(S); 
  %C   = spm_mesh_smooth(M,C,1);

  % PVE doubleside correction: 
  % If a voxel contain 38% GM and 62% CSF and has one boundary, it is approximately at the 38% position of the voxel.  
  % If the same voxel contain two boundaries, the each boundary is approximately at the 19% position of that voxel.    
  % Hence, I try to measure the filling effect in regions of two boundaries by the local minimum/maximum to estimate 
  % double the PVE effect (like a sharpening).
  % However, this is relative slow ...
  if 0 %opt.PVEcorr 
    Ypvec = cat_vol_localstat(max(1,cat_vol_localstat(min(2,max(1,Y)),Y>1,2,3)),Y>1,2,2);
    Ypvew = cat_vol_localstat(max(2,cat_vol_localstat(min(3,max(2,Y)),Y>1,2,2)),Y>1,2,3);
    Y = max(1,min(3,Y - ((max(Ypvec,Y))-Y) + (Y-(min(Ypvew,Y)))));
  end
  
  if opt.debug
    stime = cat_io_cmd('  Optimize surface:','g5','',opt.verb,stime); fprintf('\n');
  end
      
  %% Iterative correction routine
  for j=1:maxiter+1
    V   = single(SN.vertices);

    % update surface normales  
    N   = spm_mesh_normals(SN); 
    
    % inner and outer surface
    VO  = V - N .* repmat(TN/2,1,3);
    VI  = V + N .* repmat(TN/2,1,3);

    
    % First correction step that works but also could be improved.
    % ---------------------------------------------------------------------
    % Complex side specific correction by the inter-surface edges, that 
    % used the angle between the edges and normals to define edges within 
    % a suclus (outer) or within a gyrus (inner). 
    % In general, only inter-surface edges are expected here, those
    % distance describes the maximal local thickness.  We also add some 
    % sulcus-width to avoid collisions but the effect will be small due 
    % to the smoothing.
    % There are problems that not all points have a inter-surface edge, 
    % so it is necessary to smooth to include unconnected neighbors  
    % LEC, LEOC, and LEIC represent the distance error by collisions
    % of each edge.
    % ---------------------------------------------------------------------

    
    % edgelength of the central, inner, and outer surface
    % ###
    %   The edgelength is just the simples measure - the high of the
    %   tetraeder would be more excact.
    %% ###
    LE  = sum( (V(E(:,1),:)  - V(E(:,2),:)).^2  , 2) .^ 0.5;          % distance between the central surface (~ thickness/2 + thickness/2)
    LEO = sum( (VO(E(:,1),:) - VO(E(:,2),:)).^2 , 2) .^ 0.5;          % distance between the outer surface (~?minimal/maximal distance) 
    LEI = sum( (VI(E(:,1),:) - VI(E(:,2),:)).^2 , 2) .^ 0.5;          % distance between the inner surface (~ minimal/maximal distance)
    %% angle correction
    if 0
      %%
      Nalpha  = [angle(V(E(:,1),:), N(E(:,1),:)), angle(V(E(:,2),:), N(E(:,2),:))];
      Nalphas = 1; %(Nalpha>90) * 2 - 1;
      Nalpha  = min( Nalpha , abs( 180 - Nalpha )) .* Nalphas;
      LE      = min( repmat( LE ,1,2) .* cosd(Nalpha) , [] , 2); 
%%
      NalphaO  = [angle(VO(E(:,1),:), N(E(:,1),:)), angle(VO(E(:,2),:), N(E(:,2),:))];
      NalphaOs = (NalphaO>90) * 2 - 1;
      NalphaO  = min( NalphaO , abs( 180 - NalphaO )) .* NalphaOs;
      LEO      = min( repmat( LEO ,1,2) .* cosd(NalphaO)  .* (NalphaO<60),[],2); clear NalphaO;
    
      NalphaI  = [angle(VI(E(:,1),:), N(E(:,1),:)), angle(VI(E(:,2),:), N(E(:,2),:))];
      NalphaIs = (NalphaI>90) * 2 - 1;
      NalphaI  = min( NalphaI , abs( 180 - NalphaI )) .* NalphaIs;
      LEI      = min( repmat( LEI ,1,2) .* cosd(NalphaI),[],2); clear NalphaI;
    
      if 0
        %%
        %VLE = inf(size(T),'single');
        VLE = zeros(size(T),'single');
        for ni=1:size(E,1)
          %VLE(E(ni,1)) = min( VLE(E(ni,1)) , LE(ni) ./ OE(ni) );
          VLE(E(ni,1)) = VLE(E(ni,1)) + IE(ni); 
        end
      end
    end
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      LEUOS = sum( (VB(EUOS(:,1),:)  - VB(EUOS(:,2),:)).^2  , 2) .^ 0.5;  % distance b
    end
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      LEUIS = sum( (VB(EUIS(:,1),:)  - VB(EUIS(:,2),:)).^2  , 2) .^ 0.5;
    end
%TNalpha = sum( [ min( TN(E(:,1))/2, TN(E(:,1))/2 .* cosd(Nalpha(E(:,1))) )  ....* (Nalpha(E(:,1))<60) )  ...
%                 min( TN(E(:,2))/2, TN(E(:,2))/2 .* cosd(Nalpha(E(:,2))) ) ] ,2); %.* (Nalpha(E(:,2))<60) ) ] , 2); 
TNalpha = sum( [ TN(E(:,1))/2, TN(E(:,2))/2 ] , 2 ); %.* (Nalpha(E(:,2))<60) ) ] , 2); 

    LE = max(0,LE - 0.02); % minimum sulcus/gyrusweidth 

    % estimate error for each Delaunay edge 
    % (sum local thickness and sulcuswidth vs. the length of the edge) 
    %sulcuswidth = 0.0; % worse results with additional width
    LECP = max(0, LE - TNalpha); %( TN(E(:,1))/2 + TN(E(:,2))/2 + 0.02 ) ) / 2; % - sulcuswidth )); 
    LEC  = max(-inf, TNalpha - LE); %(TN(E(:,1))/2 + TN(E(:,2))/2) - (LE - 0.02) ) / 2; % - sulcuswidth )); 
    LEOC = LEO .* max(-inf, 0.02 - LEO) / 2; %clear LEO; % minimum distance between points (rare spaecial case) 
    LEIC = LEI .* max(-inf, 0.02 - LEI) / 2; %clear LEI; % minimum distance between points (rare spaecial case)
    TNP  = repmat(TN,3,1); %2 + any(opt.boundarySurfaces == [1,3]) + any(opt.boundarySurfaces == [2,3]) ,1);
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      LEUOC = max(0, LEUOS - ( TNP(EUOS(:,1))/2 ) );
    end
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      LEUIC = max(0, LEUIS - ( TNP(EUIS(:,1))/2 ) );
    end
    
    
    
    %% map the Delaunay edge correction to the vertices (simple maximum)
    % ###
    %   You may (also) use some intensity information here! ... added
    %   Moreover, a loop is very slow and the estimation of a mapping
    %   would be better! But how? ... partially implemented
    % ###
    %{
            EOid = TN*0; EIid = TN*0; 
            for ni=1:size(E,1)
              EOid(E(ni,1)) = EOid(E(ni,1)) .* , LEC(ni) .* OE(ni), LEOC(ni)]); 
              EOid(E(ni,2)) = max([EOid(E(ni,2)), LEC(ni) .* OE(ni), LEOC(ni)]); 
              EIid(E(ni,1)) = max([TIC(E(ni,1)), LEC(ni) .* IE(ni), LEIC(ni)]); 
              EIid(E(ni,2)) = max([TIC(E(ni,2)), LEC(ni) .* IE(ni), LEIC(ni)]); 
            end
    %}
    OE = min(1,(min(SNalpha,[],2)<60) + IC<2.15);
    IE = min(1,(max(SNalpha,[],2)>60) + IC>2.15); 
        
    DN  = TN*0;
    TOC = TN*0; TIC = TN*0; TOCP = TN*0; TICP = TN*0; %PVE_LEOC = TN*0; PVE_LEIC = TN*0;
    app = 0; 
    if app == 1
      for ni=1:size(E,1)
       % OE = min(1,(min(SNalpha(ni,:))<60) + IC(ni)<2.15);
       % IE = min(1,(max(SNalpha(ni,:))>60) + IC(ni)>2.15); 
        
        %{
        PVE_LEOC(E(ni,1)) = max(0,opt.vx_vol - LEOC(ni)); 
        PVE_LEOC(E(ni,2)) = max(0,opt.vx_vol - LEOC(ni)); 
        PVE_LEIC(E(ni,1)) = max(0,opt.vx_vol - LEIC(ni)); 
        PVE_LEIC(E(ni,2)) = max(0,opt.vx_vol - LEIC(ni)); 
        %}
        
      %  DN(E(ni,1)) = DN(E(ni,1)) + 1;
      %  DN(E(ni,2)) = DN(E(ni,2)) + 1; 
       
        TOC(E(ni,1)) = max([TOC(E(ni,1)), LEC(ni) .* OE(ni), LEOC(ni)]); 
        TOC(E(ni,2)) = max([TOC(E(ni,2)), LEC(ni) .* OE(ni), LEOC(ni)]); 
        TIC(E(ni,1)) = max([TIC(E(ni,1)), LEC(ni) .* IE(ni), LEIC(ni)]); 
        TIC(E(ni,2)) = max([TIC(E(ni,2)), LEC(ni) .* IE(ni), LEIC(ni)]); 
        
        % with angle weighting ...
        TOCP(E(ni,1)) = max([TOCP(E(ni,1)),LECP(ni) .* OE(ni)]); 
        TOCP(E(ni,2)) = max([TOCP(E(ni,1)),LECP(ni) .* OE(ni)]); 
        TICP(E(ni,1)) = max([TICP(E(ni,1)),LECP(ni) .* IE(ni)]); 
        TICP(E(ni,2)) = max([TICP(E(ni,1)),LECP(ni) .* IE(ni)]); 
      end
    else
      for ni=1:size(E,1)
        TOC(E(ni,1)) = max([TOC(E(ni,1)), LEC(ni) .* OE(ni), LEOC(ni)]); 
        TOC(E(ni,2)) = max([TOC(E(ni,2)), LEC(ni) .* OE(ni), LEOC(ni)]); 
        TIC(E(ni,1)) = max([TIC(E(ni,1)), LEC(ni) .* IE(ni), LEIC(ni)]); 
        TIC(E(ni,2)) = max([TIC(E(ni,2)), LEC(ni) .* IE(ni), LEIC(ni)]); 
      end
    end
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      for ni=1:size(EUOS,1)
         TOC( mod( EUOS(ni,1)-1 , size(SN,1) )+1) = LEUOC(ni); 
      end
    end
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      for ni=1:size(EUOS,1)
         TOC( mod(EUIS(ni,1)-1, size(SN,1) )+1) = LEUIC(ni); 
      end
    end
    %%
    clear LEC LEOC LIOC;
    if opt.slowdown
      slowdown  = max(1,2/j); 
    else
      slowdown  = 1; 
    end
opt.model=0;      
    if opt.model == 0 || opt.model == 2
      TOC  = single( spm_mesh_smooth(M,double(TOC), sf ))*1.4;   TOC = TOC / (slowdown/2); 
      TIC  = single( spm_mesh_smooth(M,double(TIC), sf ))*1.2;   TIC = TIC / (slowdown/2); 
      if opt.verb, fprintf('\n  TIC: %0.2f%s%0.2f, TOC: %0.2f%s%0.2f',mean(TIC),native2unicode(177, 'latin1'),std(TIC),mean(TOC),native2unicode(177, 'latin1'),std(TOC)); end
    end
    %%
    if opt.model
      % filter limits
      TOCP  = single( spm_mesh_smooth(M,double(TOCP), 1 ))*0.2;%1.0;% 1.5  
      TICP  = single( spm_mesh_smooth(M,double(TICP), 1 ))*0.8;

      % correction for intensities ...
      YI    = cat_surf_isocolors2(Y,VI); 
      YO    = cat_surf_isocolors2(Y,VO);  
      YppO  = cat_surf_isocolors2(Ypp,VO);  
     
      if opt.model == 1, fprintf('\n'); end
      if opt.verb, fprintf('  YIC: %0.2f%s%0.2f, YOC: %0.2f%s%0.2f',mean(YI),native2unicode(177, 'latin1'),std(YI),mean(YO),native2unicode(177, 'latin1'),std(YO)); end 

      WMth  = 3; YI   = max( -TICP , max(-1, min(0.5, YI - ((WMth/2 + Yl4/2) )  ))  ) / (slowdown);
      CSFth = 1; YO   = max( -TOCP , max(-1, min(0.5, ((CSFth/2 + Yl4/2) ) - YO ))  ) / (slowdown);% + 2*C
      CSFth = 0; Yppc = max( -TOCP , max(-0.05, min(0.05, 0.01 - YppO ))  ) / (slowdown);% + 2*C
      
      if 0
        YC = isocolors2(Y,V ); 
        YC = max( -0.5, min( 0.5, YC - Yl4 )) / (slowdown);
        YI = YI * 0.8 + 0.2 * YC; 
        YO = (YO * 0.8 - 0.2 * YC) .* min(1,YppO*20); 
      end
      YO = YO * 0.8 - 0.2 * Yppc; 
        
      if opt.verb, fprintf(', YIC: %0.2f%s%0.2f, YOC: %0.2f%s%0.2f',mean(YI),native2unicode(177, 'latin1'),std(YI),mean(YO),native2unicode(177, 'latin1'),std(YO)); end

      VOC = V - N .* repmat( TN/2 - YO ,1,3);    % outer surface 
      VIC = V + N .* repmat( TN/2 - YI ,1,3);    % inner surface

      YIC   = cat_surf_isocolors2(Y,VIC); 
      YOC   = cat_surf_isocolors2(Y,VOC);  
      YppOC = cat_surf_isocolors2(Ypp,VOC);  

      YI = YI .* ( abs(YI - (WMth/2 + Yl4/2))  > abs(YIC - (WMth/2 + Yl4/2)));
      YO = YO .* ( abs((CSFth/2 + Yl4/2) - YO) > abs((CSFth/2 + Yl4/2) - YOC) & YppOC>0);

      % filter correction 
      YO = single(spm_mesh_smooth(M,double(YO), sf )); 
      YI = single(spm_mesh_smooth(M,double(YI), sf ));

% if the point is/was perfect then do not change      
% if the new point is perfect / better than the old then simply use it?
      
      % combine 
      if opt.model == 1 % only intensity 
        TIC  = YI; 
        TOC  = YO; 
      elseif opt.model % combine
        mixing = max(0,min(0.5,-0.02 + 1.0*min(1,j/maxiter2)));%0.8;
        TIC  = mean( cat( 4, TIC*(1-mixing) , YI*mixing), 4); 
        TOC  = mean( cat( 4, TOC*(1-mixing) , YO*mixing), 4); 
      end  
    end
    
    
    % estimate first corrected inner and outer thickness 
    % Different levels of smoothing were use to have more effect on neighbors. 
    TOC  = single( spm_mesh_smooth(M,double( TOC  ), sf*2 )); % + spm_mesh_smooth(M,double( TOC  ), sf*4 ));  
    TIC  = single( spm_mesh_smooth(M,double( TIC  ), sf*2 )); % + spm_mesh_smooth(M,double( TIC  ), sf*4 ));
    TC   = TOC + TIC; TCsum = rms(TC(TC>0));

%{
    % correction in specified areas that also include a general
    % smoothness constrains of the cortical thickness
    TNC = TN - TOC/2 - TIC/2; 
    TNC = single(spm_mesh_smooth(M,double(TNC), sf*max(0.5,2 - 2*(j/maxiter2))) ); 
    TN(TC>0) = TNC(TC>0);   
    clear TC TNC flim;
%}
    
    % estimate new inner and outer surfaces
    VOC = V - N .* repmat( TN/2 - TOC ,1,3);    % outer surface 
    VIC = V + N .* repmat( TN/2 - TIC ,1,3);    % inner surface
    clear TOC TIC;
    
    % update thickness and surface
    TN  = sum( (VIC - VOC).^2 , 2) .^ 0.5;
    SN.vertices = mean(cat(3,VIC,VOC),3); 
    
    
    %% this is just for display and loop settings
    %  SX.vertices = VOC; SX.faces = S.faces; SX.facevertexcdata = TC; cat_surf_render2(SX);
    stopiterth = 0.00005; 
    if opt.debug && ( j==1 || mod(j,1)==0 || abs(TCsum)<0.01 || abs(TCsumo - TCsum)<stopiterth ) 
      TNM = TN>(mean(TN(:)) - 2*std(TN(:))) & TN<(mean(TN(:)) + 2*std(TN(:)));
      if ~opt.verb, fprintf('\n'); end
      try
        cat_io_cprintf('g5',sprintf('    remaining overlap:      %8.4f mm (Tlink: %4.2f%s%4.2f mm) %9.0fs',...
          TCsum,mean(TN(TNM)),native2unicode(177, 'latin1'),std(TN(TNM)),etime(clock,stime) )); stime = clock;
      end
    end
    if ( TCsum<0.005 || abs(TCsumo - TCsum)<stopiterth) && abs( mean(TN(NM)) - TNold )<0.001, break; end
    TCsumo = TCsum; TNold = mean(TN(NM)); 
  end
  
  % export cortical surfaces
  if opt.write, cat_surf_saveICO(SN,TN,mat,Pcs,sprintf('post_collcorr_%0.0fk',round( size(S.faces,1)/1000 / 10) * 10 ),0); else fprintf('\n'); end
  
  
  %% flip back
  if flipped, SN.faces = [SN.faces(:,1) SN.faces(:,3) SN.faces(:,2)]; SN.mati(7) = - SN.mati(7); end
  

end

function cat_surf_show_orthview(Psurf,Pm,color,cnames)
  fg = spm_figure('GetWin','Graphics');
  %fg = spm_figure('Create','SurfaceOverlay');%,Psurf);
  spm_figure('clear')

  id = 1;
  global st

  [pp,ff,ee] = spm_fileparts(Pm);
  hhm = spm_orthviews('Image',spm_vol(Pm));
  spm_orthviews('Caption',hhm,{'m*.nii (Intensity Normalized T1)'},'FontWeight','Bold');
  if ff(1)=='m', spm_orthviews('window',hhm,[0.3 1.03]); caxis([0.3,1.03]); end
  spm_orthviews('AddContext'); % need the context menu for mesh handling

  ov_mesh = 1; 
  for ix=1:numel(Psurf) 
    
    if ov_mesh && exist(Psurf{ix},'file')
      try
        spm_ov_mesh('display',id,Psurf{ix});
      catch
        fprintf('Please update to a newer version of spm12 for using this contour overlay\n');
        ov_mesh = 0;
        continue;
      end
    end
  end

  %% change line style
  if ov_mesh
    styles = {'b-','g-','r-','c-','m-','y-','w-','b.-','g.-','r.-','c.-','m.-','y.-','w.-'}; % need more if meshes were added
    names  = {'central';'pial';'white';'';'';''};
    if exist('color','var')
      styles(1:numel(color)) = color; 
      if exist('cnames','var')
        names(1:numel(color)) = cnames; 
      end  
    end    
    hM = findobj(st.vols{1}.ax{1}.cm,'Label','Mesh');
    UD = get(hM,'UserData');
    UD.width = [repmat(0.75,1,numel(UD.width) - numel(Psurf))  repmat(0.5,1,numel(Psurf))]; 
    UD.style = styles; %(1:numel(Psurf)); % need more if meshes were added
    set(hM,'UserData',UD);
    spm_ov_mesh('redraw',id);

    % TPM legend
    cc = axes('Position',[0.55 0.4 0.02 0.01],'Parent',fg);
    text(cc,0,1,[spm_str_manip(pp,'t') '/' ff ':']); 
    axis(cc,'off')
    for ix=1:numel(Psurf) 
      cc = axes('Position',[0.55 0.4 - 0.02*ix 0.02 0.01],'Parent',fg);
      plot(cc,[0 1],[1 1],styles{ix}); 
      text(cc,1.2,1,names{ix}); 
      axis(cc,'off')
    end
  end
end

function C = cat_surf_centroid(V,F,n)
% _________________________________________________________________________
% calculates the centroid of a region
% _________________________________________________________________________

  if ~exist('n','var'), n=1; end
  
  ndim = size(F,2);
  
  switch ndim
    case 2, ET = [1,2];
    case 3, ET = [1,2;1,3];
    case 4, ET = [1,2;1,3;1,4]; 
  end
    
  C = repmat( V(F(:,1),:) , 1, 1, n); 
  for e = 1:size(ET,1)
    ed = diff( cat( 3 , V(F(:,ET(e,1)),:) , V(F(:,ET(e,2)), :) ) , 1 , 3 ); 
    for ni = 1:numel(n)
      C(:,:,ni) = C(:,:,ni) + ni/(n+1) * ed;
    end
  end
end  

function [Yd,Yv] = cat_surf_vdist(S,V,M,opt)
% CAT surface rendering with PVE by distance approximation.
%
% [Yd,Yv] = cat_surf_render(S,V,opt)
%
%   Yd    .. distance map
%   Yv    .. surface to PVE map rendering
%   S     .. surface with vertices and faces
%   V     .. given volume or SPM volume structure
%   opt   .. option structure
%    .res .. higher surface resolution (0-default, 1-interp)
% 

% Improve speed by voxel-based pp of distance parts, if only Yv is relevant? 
%

  if ~exist('opt','var'), opt = struct(); end
  def.res  = 0;  
  def.fast = 0; 
  opt = cat_io_checkinopt(opt,def);

  % improve surface resolution?
  if opt.res
    % ...
    %S = cat_surf_fun('interp',S);
    S = cat_surf_meshinterp(S,opt.res);  
  end
  
  %% setup volume and transform vertices 
  if ~exist('V','var')
  % if not given create any volume
    Y  = false( round(max(S.vertices,[],1) - min(S.vertices)) + 10 );     
    Sv = S.vertices - repmat( min(S.vertices,[],1) - 5 , size(S.vertices,1)  , 1 );
  elseif isstruct(V)
    % modify coordinates by orientation matrix
    Y  = false( V.dims );    
    Sv = [S.vertices' ones(1,size(S.vertices,1))] .* V.mat;  
  elseif ndims(V)==3
    Y  = V; 
    Sv = S.vertices; 
  else
    % simply center the surface in the given volume
    os = round( (size(V) - (max(S.vertices,[],1) - min(S.vertices,[],1))) / 2 ); 
    Sv = S.vertices - repmat( min(S.vertices,[],1) + os , size(S.vertices,1)  , 1 );
  end

  
  %% estimate surface normals to have negative distances inside the surface
  Sn = spm_mesh_normals(S); 
  
  %% distance estimation 
  [VB,Svia] = unique( Sv , 'rows' );      % required for delaunay
  VN = Sn(Svia,:);                          % clear Sn
  VB = double(VB);                        % needed for delaunayn
  T  = delaunayn( VB );                   % delaunayn graph for faster dsearchn processing
  if exist('M','var')
    [VR(:,1),VR(:,2),VR(:,3)] = ind2sub(size(Y),find(M)); % x-y may be flipped!
  else
    [VR(:,1),VR(:,2),VR(:,3)] = ind2sub(size(Y),1:numel(Y)); % x-y may be flipped!
  end
  [VID,VDD] = dsearchn(VB,T,VR); clear T; % search nearest point with its distance
 % VB = single(VB); 
 % VM = VR - VB( VID ,:);                  % vector from surface point to voxel
  
  %% estimate if voxel is inside S 
  %  ... this is much to slow ...
  %  ... and this is unused
  %{
  VRc   = mat2cell(VR,ones(size(VR,1),1));
  VNc   = mat2cell(VN(VID,:),ones(size(VR,1),1));
  VRNa  = cellfun( @(u,v) acosd( (u*v') / (sum(u'.^2)^.5 * sum(v'.^2)^.5)) ,VRc,VNc,'UniformOutput',false);
  VSD   = cell2mat(VRNa);
  %}
  if exist('M','var')
    Yd = Y; Yd(M) = VDD;
    Yv = Y; Yv(M) = VID;
  else
    Yd = Y; Yd(1:numel(Y)) = VDD;
    Yv = Y; Yv(1:numel(Y)) = VID;
  end
  %% estimate surface normals to use a weighted
  %VMVR = mat2cell( cat( VM,VN,ones(size(VM)) , 3 ) , ones(1,size(VM,1)) ); 
  %VMVR = cellfun(@shiftdim,VMVR); 
  %VSD  = cellfun(@det,VMVR);
  
  %Yd = reshape(VDD,size(Y)) .* sign(90-reshape(VSD,size(Y))); 
  %Yv = min(1,max(0,Yd + 0.5)); 
end

function S=cat_surf_meshinterp(S,interp,method,distth)  
  if ~exist('interp','var'), interp = 1; else interp=single(interp); end
  if interp==0, return, end
  if ~exist('method','var'), method = 'linear'; end

  if ~isfield(S,'vertices')        || size(S.vertices,1)==0,        warning('Meshinterp:NoVertices','S has no vertices'); return; end
  if ~isfield(S,'faces')           || size(S.faces,1)==0,           warning('Meshinterp:NoFaces','S has no faces');       return; end
  if ~isfield(S,'facevertexcdata') || size(S.facevertexcdata,1)==0; else C=S.facevertexcdata; end
  if exist('C','var'); CT=(size(C,1)==size(S.vertices,1))+1; else CT=0; end
  
  V=S.vertices; F=single(S.faces); clear S; 
  
  for i=1:interp
    nV=single(size(V,1)); nF=single(size(F,1)); 
    
    NF=(1:nF)';
    
    switch method
      case 'linear'
       
        % addition vertices (middle of the edge)
        V1 = V(F(:,1),:) + 0.5*diff(cat(3,V(F(:,1),:),V(F(:,2),:)),1,3);
        V2 = V(F(:,2),:) + 0.5*diff(cat(3,V(F(:,2),:),V(F(:,3),:)),1,3);
        V3 = V(F(:,3),:) + 0.5*diff(cat(3,V(F(:,3),:),V(F(:,1),:)),1,3);

        % new faces which replace the old one
        F1 = [F(:,1),  nV + 2*nF + NF, nV +        NF];
        F2 = [F(:,2),  nV +        NF, nV +   nF + NF];
        F3 = [F(:,3),  nV +   nF + NF, nV + 2*nF + NF];
        F4 = [nV + NF, nV + 2*nF + NF, nV +   nF + NF];

        % colors
        if     CT==2, C=[C;nanmean(C(F(:,1),:),C(F(:,2),:));nanmean(C(F(:,2),:),C(F(:,3),:));nanmean(C(F(:,3),:),C(F(:,1),:))]; %#ok<AGROW>
        elseif CT==1, C=repmat(C,4,1);
        end
        
        V = [V;V1;V2;V3];  clear V1 V2 V3;    %#ok<AGROW>
        F = [F1;F2;F3;F4]; clear F1 F2 F3 F4; 

        % remove double vertices
        if CT==0, [V,F]   = reduce_points(V,F);
        else      [V,F,C] = reduce_points(V,F,C);
        end
        
      case 'dist'
        if ~exist('distth','var'), distth=sqrt(2); end
        
        % addition vertices (middle of the edge)
        E1 = diff(cat(3,V(F(:,1),:),V(F(:,2),:)),1,3);
        E2 = diff(cat(3,V(F(:,2),:),V(F(:,3),:)),1,3);
        E3 = diff(cat(3,V(F(:,3),:),V(F(:,1),:)),1,3);
        
        V1 = V(F(:,1),:) + repmat((sum(E1.^2,2).^0.5)>=distth,1,3) .* (0.5*E1);
        V2 = V(F(:,2),:) + repmat((sum(E2.^2,2).^0.5)>=distth,1,3) .* (0.5*E2);
        V3 = V(F(:,3),:) + repmat((sum(E3.^2,2).^0.5)>=distth,1,3) .* (0.5*E3);

        % new faces which replace the old one
        F1 = [F(:,1),  nV + 2*nF + NF, nV +        NF];
        F2 = [F(:,2),  nV +        NF, nV +   nF + NF];
        F3 = [F(:,3),  nV +   nF + NF, nV + 2*nF + NF];
        F4 = [nV + NF, nV + 2*nF + NF, nV +   nF + NF];
        
        % colors
        if     CT==2, C=[C;nanmean(C(F(:,1),:),C(F(:,2),:));nanmean(C(F(:,2),:),C(F(:,3),:));nanmean(C(F(:,3),:),C(F(:,1),:))]; %#ok<AGROW>
        elseif CT==1, C=repmat(C,4,1);
        end
        
        V = [V;V1;V2;V3];  clear V1 V2 V3;    %#ok<AGROW>
        F = [F1;F2;F3;F4]; clear F1 F2 F3 F4; 

        
        % remove double vertices
        if CT==0, [V,F]   = reduce_points(V,F);
        else      [V,F,C] = reduce_points(V,F,C);
        end

        % remove degnerated faces
        F((F(:,1)==F(:,2)) | (F(:,1)==F(:,3)) | (F(:,2)==F(:,3)),:)=[];
        
      otherwise 
        error('ERROR: Unknown method "%s"',method); 
    end
  end
  S.vertices = V; S.faces = double(F); if exist('C','var'), S.facevertexcdata = C; end
end

function [V,F,C]=reduce_points(V,F,C)
  try
    [V,~,j]  = unique(V, 'rows'); 
  catch %#ok<CTCH>
    V=single(V);
    [V,~,j]  = unique(V, 'rows'); 
  end
  j(end+1) = nan;
  F(isnan(F)) = length(j);
  if size(F,1)==1, F = j(F)'; if exist('C','var'), C=j(C)'; end    
  else             F = j(F);  if exist('C','var'), C=j(C);  end    
  end
end

function cdata2 = cat_surf_surf2vol2surf(S,S2,cdata,res)
 % create volume 
 % render data
 % aprax
 % proejcet 
end

function [Yp,Yt,vmat1,vmat1i] = cat_surf_surf2vol(S,Y,T,type,opt)
% cat_surf_surf2vol
% _________________________________________________________________________
%
% Render a surface mesh and/or its data to a (given) volume. 
% The volume Y can be used as mask to reduce processing time.
%
% [Yp,Yv,vmat,vmati] = cat_surf_surf2vol(S,Y,T,type,opt)
%
%   S       .. surface 
%   Y       .. volume for rendering 
%   T       .. thickness
%   type    .. ['pve','seg','val']
%     'pve' .. default that simply render the surface into a (given) volume
%     'pp'  .. create percentage position map
%     'seg' .. create a segmentation label map with the thickness map T
%     'val' .. just render the value (fast)
%   opt     .. additional parameter
%   .refine .. refine mesh to improve accuracy (default 1.8) 
%   .acc    .. render more layer to improve accuracy (default 0)
%   .verb   .. display
% _________________________________________________________________________
% Robert Dahnke 201911


% TODO:
%  - T as [1x2] and [T,T] vector with lower and upper enlargement
%  - render of negative values T + minT ... Yt - minT
%  - triangle/edge rendering and no refinement

  global mati vmati

  ftime = clock;
  
  if ~exist('Y','var'),    Y = []; end
  if ~exist('T','var'),    T = []; end
  if ~exist('type','var'), type = 'pve'; end
  if ~exist('opt','var'),  opt = struct(); end
  
  def.debug      = 0;     % debugging output vs. memory optimization
  def.refine     = 1.8;   % mesh refinement - more points more stable but slower
  def.acc        = 1;     % larger value - more layer - more exact but slower
  def.bdist      = 5;     % default volume 
  def.fsize      = 10;
  def.verb       = 1; 
  def.testseggmt = 0;
  def.interpBB   = struct('BB',[],'interpV',1,'mati',mati,'vmati',vmati);
  def.mat        = [];
  opt = cat_io_checkinopt(opt,def);

  if ~isempty(T) && size(S.vertices,1) ~= numel(T)
    error('The number of surface vertices and datavalues has to be equal.\n');
  end
  if strcmpi(type,'pp') && isempty(T)
    error('Position map creation requires a thickness map T.\n');
  end

  if ~isempty(opt.mat)
    vx_vol = sqrt(sum(opt.mat(1:3,1:3).^2));
  else
    vx_vol = [1 1 1];
  end
  opt.interpBB.interpV = vx_vol(1);
  
  %% save a temporary version of S and refine it
  if ~strcmpi(type,'val')
    stime = cat_io_cmd('  Refine mesh','g5','',opt.verb);
    Praw = [tempname '.gii'];
    save(gifti(struct('vertices',S.vertices,'faces',S.faces)),Praw); 
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,opt.refine .* mean(vx_vol)); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    Sr = gifti(Praw);
    delete(Praw);
  end
  
  
  %% transformation for create CS
  if 0 ~isempty(opt.interpBB.mati) && ~isempty(opt.interpBB.BB)
    S.vertices = S.vertices - repmat( opt.interpBB.BB([3,1,5]) - 1,size(S.vertices,1),1);            % correction for boundary box 
    S.vertices = S.vertices ./ repmat(abs(opt.interpBB.interpV ./ opt.interpBB.mati([8,7,9])),size(S.vertices,1),1);   % resolution adaption

    vmat1 = [0 0 0];

    if ~strcmpi(type,'val')
      Sr.vertices = Sr.vertices - repmat( opt.interpBB.BB([3,1,5]) - 1,size(Sr.vertices,1),1);            % correction for boundary box 
      Sr.vertices = Sr.vertices ./ repmat(abs(opt.interpBB.interpV ./ opt.interpBB.mati([8,7,9])),size(Sr.vertices,1),1);   % resolution adaption
    end
  elseif ~isempty(opt.mat)
    S = cat_surf_mat(S,opt.mat,1);
  end
  %%
  if isempty(Y)
    Y      = ones( round(max(S.vertices,[],1) - min(S.vertices)) + opt.bdist*2 ,'single');     
    vmat1  = -[min(S.vertices(:,1)) min(S.vertices(:,2)) min(S.vertices(:,3))] + opt.bdist; 
    vmat1i = repmat(min(S.vertices),size(S.vertices,1),1) - opt.bdist; 
  else   
    if sum(Y(:)) == 0, Y = Y + 1; end
    
    %{
    vmat1i = [0 0 0]; 
    vmat1  = -opt.interpBB.vmati(10:12); 
    
    S.vertices  = (opt.interpBB.vmati * [S.vertices'  ; ones(1,size(S.vertices ,1))])';
    if ~strcmpi(type,'val')
      Sr.vertices = (opt.interpBB.vmati * [Sr.vertices' ; ones(1,size(Sr.vertices,1))])';
    end
    if opt.interpBB.vmati(7)<0, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; end
    %}
    vmat1  = [0 0 0];% -opt.interpBB.vmati(10:12); 
    S.vertices = ([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat1'] * [S.vertices';ones(1,size(S.vertices,1))] )';
    if ~strcmpi(type,'val')
      Sr.vertices = ([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat1'] * [Sr.vertices';ones(1,size(Sr.vertices,1))] )';
    end
  end
  
  %% transfer thickness information to refined surface by volume rendering
  if ~isempty(T)
    if exist('stime','var')
      stime = cat_io_cmd('  Render thickness/data','g5','',opt.verb,stime);
    else
      stime = cat_io_cmd('  Render thickness/data','g5','',opt.verb);
    end
    Yt = nan(size(Y),'single'); 
    % render points
    I  = sub2ind(size(Y),...
        max(1,min(size(Y,1),round(S.vertices(:,1) + vmat1(1)))),...
        max(1,min(size(Y,2),round(S.vertices(:,2) + vmat1(2)))),...
        max(1,min(size(Y,3),round(S.vertices(:,3) + vmat1(3)))));
    Yt(I) = T; 
    if 1 % fill volume
      if all(T == round(T))
        [D,I] = cat_vbdist( single( ~isnan(Yt) ) , Y>0 ); 
        Yt  = Yt(I); 
      else
        Yt  = cat_vol_approx(Yt,1);
      end
    end
    %T    = isocolors(Yt,([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat'] * [So.vertices';ones(1,size(So.vertices,1))] )' ); % self projection
  else
    Yt = nan(size(Y),'single');
  end
  if strcmpi(type,'val')
    if (opt.verb) > 0,  fprintf('%5.0fs\n',etime(clock,stime)); end
    Yt  = Yt .* (Y>0);
    Yp  = Yt;
    return;
  end

  
 
  %% smooth the normals to avoid problems with self-intersections
  stime = cat_io_cmd('  Smooth normals','g5','',opt.verb,stime);
  Mr = spm_mesh_smooth(Sr);   
  smoothsurf = @(Y,s) [ ...        
    spm_mesh_smooth(Mr,double(Y(:,1)),s) , ...
    spm_mesh_smooth(Mr,double(Y(:,2)),s) , ...
    spm_mesh_smooth(Mr,double(Y(:,3)),s) ];
  Srn = spm_mesh_normals(Sr); 
  Srn = smoothsurf(Srn, min(100,max(50,opt.fsize * (numel(Sr.vertices)./numel(S.vertices))^0.5 ) ) );
  Srn = Srn ./ repmat( sum(Srn.^2,2).^0.5 , 1 , 3);
  
  
  %% render surface
  stime = cat_io_cmd('  Render final map','g5','',opt.verb,stime);
  switch lower(type)
    case {'pve','seg','pp'}
      %% render surface with PVE 
      ss      = max(1,opt.acc + 1) / 12;
      offset  = -0.25 : ss : 0.75 * round(opt.refine/0.75);
      % the transverse offset should include only few elements for runtime
      toffset = -1/3 * round(opt.refine/1.5) : (2/3) / max(1,round(opt.refine/1.5)) :  1/3 * round(opt.refine/1.5);
      Yp = zeros(size(Y),'single');
      for oi = 1:numel(offset)
        % middle point
        I = sub2ind(size(Y),...
            max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi) + vmat1(1)))),...
            max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi) + vmat1(2)))),...
            max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi) + vmat1(3)))));
        % diagonal elments
        if 1
          for ti1 = 1:numel(toffset)
            for ti2 = 1:numel(toffset)
              for ti3 = 1:numel(toffset)
                I = [I; sub2ind(size(Y),...
                    max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi) + toffset(ti1)*1/sqrt(3) + vmat1(1)))),...
                    max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi) + toffset(ti2)*1/sqrt(3) + vmat1(2)))),...
                    max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi) + toffset(ti3)*1/sqrt(3) + vmat1(3)))))];
              end
            end
          end
        end
        % direct neigbors
        if 1
          for ti1 = 1:numel(toffset)
            I = [I; sub2ind(size(Y),...
                max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi) + toffset(ti1) + vmat1(1)))),...
                max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi) + vmat1(2)))),...
                max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi) + vmat1(3)))))];
          end
          for ti2 = 1:numel(toffset)
            I = [I; sub2ind(size(Y),...
                max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi)                + vmat1(1)))),...
                max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi) + toffset(ti2) + vmat1(2)))),...
                max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi)                + vmat1(3)))))];
          end
          for ti3 = 1:numel(toffset)
            I = [I; sub2ind(size(Y),...
                max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi)                + vmat1(1)))),...
                max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi)                + vmat1(2)))),...
                max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi) + toffset(ti3) + vmat1(3)))))];
          end
        end
        % final rendering
        Yp(I) =  min(1, oi ./ sum( offset < 1 ) );
      end
    case 'pps'
      %% render percentage position map
      % read thickness value from map

      Tr    = isocolors(Yt,([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat1'] * [Sr.vertices';ones(1,size(Sr.vertices,1))] )' );
      %Tr    = isocolors(Yt,Sr.vertices);

      % create Ypp map
      ss      = max(1,opt.acc - 1) / round(4 * round(mean(T)));
      offset  = -0.25 : ss : 0.75 ; 
      toffset = -1/3 * round(opt.refine/1.5) : (2/3) / max(1,round(opt.refine/1.5)) :  1/3 * round(opt.refine/1.5);
      Yp = zeros(size(Y),'single');
      for oi = 1:numel(offset)
           % middle point
        I = sub2ind(size(Y),...
            max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi).* max(offset(oi)>0.5,Tr)./vx_vol(1) + vmat1(1)))),...
            max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi).* max(offset(oi)>0.5,Tr)./vx_vol(2) + vmat1(2)))),...
            max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi).* max(offset(oi)>0.5,Tr)./vx_vol(3) + vmat1(3)))));
        % diagonal elments
        if 1
          for ti1 = 1:numel(toffset)
            for ti2 = 1:numel(toffset)
              for ti3 = 1:numel(toffset)
                I = [I; sub2ind(size(Y),...
                    max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV + toffset(ti1)*1/sqrt(3) + vmat1(1)))),...
                    max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV + toffset(ti2)*1/sqrt(3) + vmat1(2)))),...
                    max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV + toffset(ti3)*1/sqrt(3) + vmat1(3)))))];
              end
            end
          end
        end
        % direct neigbors
        if 1
          for ti1 = 1:numel(toffset)
            I = [I; sub2ind(size(Y),...
                max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi).* max(offset(oi)>0.5,Tr)./vx_vol(1) + toffset(ti1) + vmat1(1)))),...
                max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi).* max(offset(oi)>0.5,Tr)./vx_vol(1)                + vmat1(2)))),...
                max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV                + vmat1(3)))))];
          end
          for ti2 = 1:numel(toffset)
            I = [I; sub2ind(size(Y),...
                max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV                + vmat1(1)))),...
                max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV + toffset(ti2) + vmat1(2)))),...
                max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV                + vmat1(3)))))];
          end
          for ti3 = 1:numel(toffset)
            I = [I; sub2ind(size(Y),...
                max(1,min(size(Y,1),round(Sr.vertices(:,1) + Srn(:,1)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV                + vmat1(1)))),...
                max(1,min(size(Y,2),round(Sr.vertices(:,2) + Srn(:,2)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV                + vmat1(2)))),...
                max(1,min(size(Y,3),round(Sr.vertices(:,3) + Srn(:,3)*offset(oi).* max(offset(oi)>0.5,Tr)./opt.interpBB.interpV + toffset(ti3) + vmat1(3)))))];
          end
        end

        Yp(I) = min( 1, oi./ sum( offset <= 0.5 )  );


      end
      Yp  = cat_vol_median3(Yp,smooth3(Yp>0 & Yp<1)>0.2,Yp>0,0.05);
      Yps = cat_vol_localstat(Yp,Yp>0 & Yp<1,1,1); Ypsw = min(1,smooth3(cat_vol_morph(Yps>0,'e'))*2); 
      Yp  = Yp .* (1 - Ypsw) + Ypsw .* Yps; clear Yps Ypsw;
  end
  
  % filling
  Yp0 = Yp; 
  for fi=1:2
    if sum( abs( Yp0(:) - Yp(:) ) ) / sum( Yp0(:) ) < 0.1 
      Yp(~cat_vol_morph( cat_vol_morph( Yp<(1-0.2*fi) ,'o',fi), 'l',[0.2 2]) & Yp==0 ) = 1;
    end
  end
  clear Yp0;
  % final filtering
  Yp = cat_vol_median3(Yp,smooth3(Yp>0 & Yp<1)>0.2,Yp>0,0.05);
  
  
  
  %% create segment map
  switch lower(type)     
    case 'pp'
      %% distance map for PVE
      %  optimized for real Ypp map
      Yd  = cat_vbdist( 1 - Yp ) - cat_vbdist( Yp , cat_vol_morph( Y , 'd', 2) );
      Yd(Yd<0) = Yd(Yd<0) + 0.5; % use correction
      Yd(Yd>0) = Yd(Yd>0) - 0.5;
      Yd  = Yd * opt.interpBB.interpV; 

      % create final segmenation 
      Yp  = max(0,min(1,(Yd + Yt/2) ./ Yt));
      
      Yp  = cat_vol_median3(Yp,smooth3(Yp>0 & Yp<1)>0.2,Yp>0,0.05);
      Yps = cat_vol_localstat(Yp,Yp>0 & Yp<1,1,1); Ypsw = min(1,smooth3(cat_vol_morph(Yps>0,'e'))*2); 
      Yp  = Yp .* (1 - Ypsw/2) + Ypsw/2 .* Yps; clear Yps Ypsw;
    case 'seg' 
      %% distance map for PVE
      %  optimized for thickness error
      Yd = cat_vbdist( 1 - Yp ) - cat_vbdist( Yp , cat_vol_morph( Y , 'd', 2) );
      Yd(Yd<0) = Yd(Yd<0) + 0.0; % no correction
      Yd(Yd>0) = Yd(Yd>0) - 0.0;
      Yd = Yd * opt.interpBB.interpV; 

      % create CSF background 
      [Yb,red] = cat_vol_resize(Yp,'reduceV',opt.interpBB.interpV,2,8,'meanm');
      Yb = cat_vol_morph( Yb>0 ,'ldc',16); 
      Yb = cat_vol_resize(Yb,'dereduceV',red);
      Yb = Yb | cat_vol_morph( Yp , 'd' ,2 );
      
      % create final segmentation 
      Yp = Yb + max(0,min(1, (Yd + Yt / 2) / opt.interpBB.interpV ) ) + max(0,min(1,(Yd - Yt / 2) / opt.interpBB.interpV ) ); 
      
      if opt.testseggmt
        %% test rendering 
        stime = cat_io_cmd('  Test segment map','g5','',opt.verb,stime);
        Ycd  = cat_vbdist( 2 - Yp , Yp<2.5 ); 
        Ywd  = cat_vbdist( Yp - 2 , Yp>1.5 );
        Ygmt = (Ywd + Ycd) * opt.interpBB.interpV; Ygmt(Ygmt>1000) = 0;
        Ypbt = min(Ywd + Ycd,cat_vol_pbtp(Yp,Ywd,Ycd)) * opt.interpBB.interpV; Ypbt(Ypbt>1000) = 0;

        % values
        rms = @(x) mean( x.^2 ) .^ 0.5;   
        Tgmt = cat_surf_isocolors2(Ygmt, ([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat1'] * [S.vertices';ones(1,size(S.vertices,1))] )' );
        Tpbt = cat_surf_isocolors2(Ypbt, ([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat1'] * [S.vertices';ones(1,size(S.vertices,1))] )' );
        
        fprintf('\n  T_surf:    %0.2f%s%0.2f (md=%0.2f)\n',mean(T),native2unicode(177, 'latin1'),std(T),median(T));
        fprintf('  T_direct:  %0.2f%s%0.2f (md=%0.2f, RMSE=%0.2f)\n',...
          mean(Ygmt(Ygmt(:)>0)),native2unicode(177, 'latin1'),std(Ygmt(Ygmt(:)>0)),median(Ygmt(Ygmt(:)>0)),rms(T - Tgmt));
        fprintf('  T_pbtfast  %0.2f%s%0.2f (md=%0.2f, RMSE=%0.2f)\n',...
          mean(Ypbt(Ypbt(:)>0)),native2unicode(177, 'latin1'),std(Ypbt(Ypbt(:)>0)),median(Ypbt(Ypbt(:)>0)),rms(T - Tpbt));
        cat_io_cmd(' ','g5','',opt.verb);
      end
  end
  
  if exist('Yt','var')
    Yt  = Yt .* (Y>0);
  end
  if opt.verb
    fprintf('%5.0fs\n',etime(clock,stime)); 
    cat_io_cmd(' ','g5','',opt.verb);
    fprintf('%5.0fs\n',etime(clock,ftime)); 
  end
end

function [V,vmat,vmati] = cat_surf_surf2vol_old(S,opt)
%% render inner surface area 
%  Render the volume V with V==1 within the surface. 
%  Use type=1 to render also the surface area with 0.5.
%  The transformation imat to create 
%  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)];     % matlab flip
%  SH.vertices = SH.vertices + imat;

  if ~exist('opt','var'), opt = struct(); end
  def.debug  = 0;     % debugging output vs. memory optimization 
  def.pve    = 1;     % 0 - no PVE; 1 - PVE; 1 - PP with thickness (opt.T)
                      % 2 - fill with surface texture values without interpolation and masking (==4)    
                      % 3 - fill with surface texture values with    interpolation and masking (==5)
  def.refine = 0.6;
  def.bdist  = 5; 
  def.mat    = [];
  def.res    = 1; % not yet ...
  
  opt = cat_io_checkinopt(opt,def);
  
  % save a temporary version of S and refine it
  if 1
  
    Praw = [tempname '.gii'];
    save(gifti(struct('vertices',S.vertices,'faces',S.faces)),Praw); 

    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,opt.refine .* opt.interpBB.interpV); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

    So = S; S = gifti(Praw);
    delete(Praw);
  else
    So = S;
  end
  
  if ~isempty(opt.mat)
    S = cat_surf_mat(S,mat); 
  end
  
  
  %% render surface points
  if ( isfield(opt,'V') || isfield(opt,'sizeV') ) && isfield(opt,'vmati')
    if isfield(opt,'V')
      V     = single(opt.V); 
    else
      V     = zeros(opt.sizeV,'single');
    end
    
    vmati = [0 0 0]; 
    vmat  = -opt.vmati(10:12); 
    
    S.vertices  = (opt.vmati * [S.vertices'  ; ones(1,size(S.vertices ,1))])';
    So.vertices = (opt.vmati * [So.vertices' ; ones(1,size(So.vertices,1))])';
    if opt.vmati(7)<0, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; end
  else
    if opt.pve == 0 
      V   = false( round(max(S.vertices,[],1) - min(S.vertices)) + opt.bdist*2 );    
    else
      V   = zeros( round(max(S.vertices,[],1) - min(S.vertices)) + opt.bdist*2 ,'single');     
    end
    vmat  = -[min(S.vertices(:,1)) min(S.vertices(:,2)) min(S.vertices(:,3))] + opt.bdist; 
    vmati = repmat(min(S.vertices),size(S.vertices,1),1) - opt.bdist; 
  end
  if opt.pve > 1 && opt.pve < 6
    % get surface data or give error
    if isfield(So,'cdata')
      cdata = So.cdata; 
    elseif  isfield(So,'facevertexcdata')
      cdata = So.facevertexcdata; 
    else
      error('cat_surf_fun:cat_surf_surf2vol:No datafield for filling');
    end
    
    % render data 
    I    = sub2ind(size(V),...
          max(1,min(size(V,1),round(S.vertices(:,1) + vmat(1)))),...
          max(1,min(size(V,2),round(S.vertices(:,2) + vmat(2)))),...
          max(1,min(size(V,3),round(S.vertices(:,3) + vmat(3)))));
    V(I) = 1; 
    V    = cat_vol_morph(V,'lc',1);  % closeing 
    
    % data filling
    if exist('cdata','var')
      Vv   = zeros( size(V) ,'single');  % same size so S and not So   
      I    = sub2ind(size(V),...
            max(1,min(size(V,1),round(So.vertices(:,1) + vmat(1)))),...
            max(1,min(size(V,2),round(So.vertices(:,2) + vmat(2)))),...
            max(1,min(size(V,3),round(So.vertices(:,3) + vmat(3)))));
      Vv(I) = cdata;
    end
    if opt.pve == 2 || opt.pve == 4
      [D,I] = vbdist(single(V )); Vv = Vv(I); clear D; 
      if opt.pve<4
        [D,I] = vbdist(single(~V)); Vv = Vv(I); clear D,
      end
    else
      Vv    = cat_vol_approx(Vv); 
    end
    
    % final masking
    if opt.pve > 3
      V    = Vv .* V; 
    else
      V    = Vv; 
    end 
    clear Vv; 
    
  elseif opt.pve == 0
    %%    
    I    = sub2ind(size(V),...
          max(1,min(size(V,1),round(S.vertices(:,1) + vmat(1)))),...
          max(1,min(size(V,2),round(S.vertices(:,2) + vmat(2)))),...
          max(1,min(size(V,3),round(S.vertices(:,3) + vmat(3)))));
    V(I) = 1; 

    V    = cat_vol_morph(V,'lc',1);  % closeing 
    V(I) = 0;                        % remove points of the surface
    V    = cat_vol_morph(V,'lab');   % final closing
  elseif opt.pve == 6
    % render a nice volume
    S2 = struct('faces',So.faces,'vertices',...
         [max(1,min(size(V,1),round(So.vertices(:,1) + vmat(1)))),...
          max(1,min(size(V,2),round(So.vertices(:,2) + vmat(2)))),...
          max(1,min(size(V,3),round(So.vertices(:,3) + vmat(3))))]);

    [D,I] = cat_surf_vdist(S2,V,V>0);
    D = D * opt.interpBB.interpV;
    
    YT = zeros(size(V),'single'); 
    YT( I>0 ) = opt.T( I( I>0 )); 
    YT = cat_vol_localstat(YT,YT>0,2,1);
    
    %%
    % refinement did not work and is also computation stupid and slow
    %[D2,I2] = cat_surf_vdist(S2,V>0 & D<2,V>0 & D<2,struct('res',1)); 
    %I(V>0 & D<2) = I2(V>0 & D<2); D(V>0 & D<2) = D2(V>0 & D<2);

    % % das geht natuer nach dem refinement nicht ...
    
%    opt.res
    if isfield(opt,'T')
      V = max(0,min(1, (D + T(I)/2) / T(I) ));
    else
      V = max(0,min(1, D + 0.5 ));
    end
    
  else %if opt.pve == 1
    %% fast PVE estimation by rendering multiple layer 
    
    
    %% part 1 .. render one surface with simple PVE
    
    % smooth the normals to avoid problems with self-intersections
    M = spm_mesh_smooth(S);   
    smoothsurf = @(V,s) [ ...        
      spm_mesh_smooth(M,double(V(:,1)),s) , ...
      spm_mesh_smooth(M,double(V(:,2)),s) , ...
      spm_mesh_smooth(M,double(V(:,3)),s) ];
    Sn = spm_mesh_normals(S); 
    Sn = smoothsurf(Sn,100);
    
    % render surface 
    offset = -0.25:0.1:1; 
    V=0*V;
    for oi = 1:numel(offset)
      I = sub2ind(size(V),...
          max(1,min(size(V,1),round(S.vertices(:,1) + Sn(:,1)*offset(oi) + vmat(1)))),...
          max(1,min(size(V,2),round(S.vertices(:,2) + Sn(:,2)*offset(oi) + vmat(2)))),...
          max(1,min(size(V,3),round(S.vertices(:,3) + Sn(:,3)*offset(oi) + vmat(3)))));
      V(I) = min(1,max( V(I) , oi./sum(offset<1))); 
    end
    V(~cat_vol_morph( V<1, 'l',[0.2 2]) & V==0 ) = 1;  % closing 
    
    
    
    %% part 2 .. surf2vol2surf
    Yt = zeros(size(V),'single'); 
    I  = sub2ind(size(V),...
        max(1,min(size(V,1),round(So.vertices(:,1) + vmat(1)))),...
        max(1,min(size(V,2),round(So.vertices(:,2) + vmat(2)))),...
        max(1,min(size(V,3),round(So.vertices(:,3) + vmat(3)))));
    Yt(I) = opt.T; 
    Yt    = cat_vol_approx(Yt,1); 
    
    %T = isocolors(Yt,([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat'] * [So.vertices';ones(1,size(So.vertices,1))] )' ); % self projection
    T  = isocolors(Yt,([0 1 0; 1 0 0; 0 0 1] *  [eye(3) vmat'] * [S.vertices';ones(1,size(S.vertices,1))] )' );
    
    %% create Ypp map
    offset = [-0.5:0.05:0.5 0.75:0.25:1]; 
    Vp=0*V;
    for oi = 1:numel(offset)
      I = sub2ind(size(V),...
          max(1,min(size(V,1),round(S.vertices(:,1) + Sn(:,1).*offset(oi) .* T./opt.interpBB.interpV + vmat(1)))),...
          max(1,min(size(V,2),round(S.vertices(:,2) + Sn(:,2).*offset(oi) .* T./opt.interpBB.interpV + vmat(2)))),...
          max(1,min(size(V,3),round(S.vertices(:,3) + Sn(:,3).*offset(oi) .* T./opt.interpBB.interpV + vmat(3)))));
      Vp(I) = min( 1.5, max( Vp(I) , offset(oi) + 0.5 )); 
    end
    Vp(~cat_vol_morph( Vp<1, 'l',[0.2 2]) & Vp==0 ) = 1;
    Vp(cat_vol_morph( Vp>1.2, 'lc',1) & Vp==0 ) = 1;
    Vp(cat_vol_morph( Vp>=1, 'lc',0) & Vp>0) = 1;
    Vp = min(1,Vp);
  end
  
  %%
  %SH = isosurface(V,0.6);        % create hull 
  %SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  %SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end

function alpha = cat_surf_edgeangle(N1,N2)
%cat_surf_fun>cat_surf_edgeangle Estimate angle between two vectors. 
%
%  alpha = cat_surf_edgeangle(N1,N2)
%  ________________________________________________________________________
%  Robert Dahnke 201909

  if 1 % fast version
    alpha = acosd( dot(N1,N2,2) ./ ( sum(N1.^2,2).^.5 .* sum(N2.^2,2).^.5 ));
  else
    %%
    alpha = zeros(size(N1,1),1,'single');
    for i=1:size(N1,1)
      a = N1(i,:); b = N2(i,:); 
      alpha(i) = acosd( dot(a,b) / (norm(a) * norm(b)) ); 
    end
  end
end

function N = cat_surf_normals(S) 
%cat_surf_normals Surface normals. 
%  Vertex normals of a triangulated mesh, area weighted, left-hand-rule 
%  N = patchnormals(FV) - struct with fields, faces Nx3 and vertices Mx3 
%  N: vertex normals as Mx3
%
%  https://de.mathworks.com/matlabcentral/fileexchange/24330-patch-normals
%  by Dirk-Jan Kroon
%
%  See also spm_mesh_normals.
%  ________________________________________________________________________
%  Dirk-Jan Kroon, Robert Dahnke, 2019

%  Well, I now use spm_mesh_normals, but maybe we need it some day.

  % face corners index 
  A = S.faces(:,1); 
  B = S.faces(:,2); 
  C = S.faces(:,3);

  % face normals 
  n = cross(S.vertices(A,:)-S.vertices(B,:),S.vertices(C,:)-S.vertices(A,:)); %area weighted

  % vertice normals 
  N = zeros(size(S.vertices));      % init vertex normals 
  for i = 1:size(S.faces,1)         % step through faces (a vertex can be reference any number of times) 
    N(A(i),:) = N(A(i),:) + n(i,:); % sum face normals 
    N(B(i),:) = N(B(i),:) + n(i,:); 
    N(C(i),:) = N(C(i),:) + n(i,:); 
  end
  
  % normalize
  Ns = sum(N.^2,2).^.5;
  N  = N ./ repmat(Ns,1,3); 
end

function S = cat_surf_mat(S,mat,invers,nin)
%cat_surf_fun>cat_surf_mat 
%  Apply transformation matrix mat to a surface structure S. 
%
%  S = cat_surf_fun('smat',S,mat[,inv,nin])
%  S = cat_surf_mat(S,mat[,inv,nin])
% 
%  S        .. surface structure with vertices and faces
%  mat      .. tranformation matrix
%              if empty only check for vertex inversion
%  invers   .. apply invers transformation (default = 0)
%  nin      .. set normal direction pointing to the inside (default = 1)
%
%  Used in cat_surf_createCS.
%  
%  See also spm_mesh_transform.
%  ------------------------------------------------------------------------
%  Robert Dahnke 201911

  if ~exist('invers','var'), invers = 0; end
  if ~exist('nin','var'), nin = 1; end 
  if ~exist('mat','var'), error('cat_surf_fun:cat_surf_mat:incompleteInput','Need transformation matrix as input!'); end 

  if ~isempty(mat)
    if isnumeric(S)
      % only vertices
      if invers
        S          = ( inv( mat ) * [S' ; ones(1,size(S,1))])'; 
      else
        S          = ( mat  * [S' ; ones(1,size(S,1))])'; 
      end
      S(:,4) = [];  
    else %if isstruct(S) % surface
      if invers
        S.vertices = ( inv( mat ) * [S.vertices' ; ones(1,size(S.vertices,1))])'; 
      else
        S.vertices = ( mat * [S.vertices' ; ones(1,size(S.vertices,1))])'; 
      end
      S.vertices(:,4) = []; 
    end  
  end
  
  % estimate the orientation of the surface by adding the normals and
  % testing which directions makes the surface great again :D
  if isstruct(S)
    N      = spm_mesh_normals(S); 
    CSnin  = mean( abs(S.vertices(:) - N(:) ) ) > mean( abs(S.vertices(:) + N(:) ) ); 

    % flip faces if normals point into the wrong direction
    if ( nin && ~CSnin ) ||  ( ~nin && CSnin )
      S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; 
    end
  end
end

function I = cat_surf_isocolors2(V,Y,mat,interp)
%cat_surf_fun>cat_surf_isocolors2 Map volume data to surface.
%  Calculates an interpolated value of a vertex V in Y.  
%
%  I = cat_surf_isocolors2(V,Y,mat,interp)
% 
%  I      .. facevertexcdata
%  V      .. vertices or surface struct
%  mat    .. transformation matrix
%  interp .. interpolation type ('nearest','linear)
%
%  See also isocolors, spm_mesh_project.
%  _____________________________________________________________________
%  Robert Dahnke, 2009-2019
  
  if isempty(Y), return; end
  if ~isnumeric(V) && isfield(V,'vertices'), V = V.vertices; end 
  if ~isnumeric(Y) && isfield(Y,'vertices'), Y = Y.vertices; end 
  if ~ismatrix(V) && ndims(Y)~=3, VR = Y; Y = V; V = VR; end; clear VR % flip
  if ndims(Y)~=3,  error('MATLAB:isocolor2:dimsY','Only 2 or 3 dimensional input of Y.'); end
  if ~exist('interp','var'); interp = 'linear'; end 
  
  if isa(V,'double'), V = single(V); end
  if isnumeric(Y)
    if ~isa(Y,'double')
      Y  = double(Y);
      YD = 0;
    else
      YD = 1;
    end
  else
    Ygii = Y; clear V
    Y    = single(Ygii.vertices); 
    YD   = 0; 
  end

  
  % inverse transformation of the given mat
  if exist('mat','var') && ~isempty(mat)
    V = cat_surf_mat(V,mat,1); 
  end
  
  
  nV   = size(V,1);
  ndim = size(V,2);
 
  
  % We have to process everything with double, thus larger images will cause memory issues.
  switch interp
    case 'nearest'
      V = max(1,min(round(V),repmat(size(Y),nV,1))); 
      I = Y( max(1,min(numel(Y), sub2ind(size(Y),V(:,2),V(:,1),V(:,3)))));
    case 'linear'
      nb  = repmat(shiftdim(double([0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]'),-1),nV,1);  
      enb = repmat(shiftdim((ones(8,1,'double')*[size(Y,2),size(Y,1),size(Y,3)])',-1),nV,1);  

      % calculate the weight of a neigbor (volume of the other corner) and
      w8b = reshape(repmat(V,1,2^ndim),[nV,ndim,2^ndim]); clear V;
      % if the streamline is near the boundary of the image you could be out of range if you add 1 
      n8b = min(floor(w8b) + nb,enb); clear enb
      n8b = max(n8b,1);
      w8b = flipdim(prod(abs(n8b - w8b),2),3);        

      % multiply this with the intensity value of R
      I = sum( Y( max(1,min(numel(Y), sub2ind(size(Y),n8b(:,2,:),n8b(:,1,:),n8b(:,3,:))))) .* w8b,3);
  end  
  if ~YD, I = single(I); end
end