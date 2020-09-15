function varargout = cat_surf_load(files,sides,data)
%cat_surf_load. Load gifti surface and combine them to one mesh.  
%
%  [S,sdata,fdata,rnames] = cat_surf_load(file,sides,loadmesh)
% 
%  S      .. surface patch structure with vertices, faces, and
%            facevertexcdata depending on the data variable
%  sdata  .. facevertexcdata of side
%  fdata  .. facevertexcdata of object
%  rnames .. region names in case of atlas files 
%
%  file   .. surface or texture 
%  sides  .. load specific sides 
%             'mesh' (default), 'lh', 'rh' ... 'meshcb', 'cb'
%  data   .. load surface mesh and/or surface data
%              1-mesh only, 2-texture only, 3-load both (default)
%  avg    .. use average mesh 
%              1-FreeSurfer, 2-Shooting 
%  type?  .. central, pial, white, hull

%#ok<*AGROW>

  % defaults
  if ~exist('sides','var'), data = 'mesh'; end
  if ~exist('data' ,'var'), data = 3; end
  if ~exist('avg'  ,'var'), avg  = 1; end  

  
  % check all input files
  files = cellstr(files);
%   for fi=1:numel(files)
%     if ~exist(files{fi},'file')
%       error('cat_surf_loadfiles:noFile','Input file %d "%s" does not exist.',fi,files{fi}); 
%     end
%   end


  % check side input string
  if iscell(sides) 
    if strfind(sides,'mesh')
      sides = unique( [sides, 'lh', 'rh'] ); 
    elseif strfind(sides,'meshcb')
      sides = unique( [sides, 'lh', 'rh', 'cb' ] ); 
    end  
  else
    if strcmp(sides,'mesh')
      sides = {'lh', 'rh'};
    elseif strcmp(sides,'meshcb')
      sides = {'lh', 'rh', 'cb'}; 
    end
  end
  
  
  %% check & prepare matlab patch data
  if data==1 || data==3
    S.vertices        = []; 
    S.faces           = []; 
  end
  if data==2 || data==3
    S.facevertexcdata = [];
  end
  sdata               = uint8([]); 
  odata               = uint8([]); 
  

  % check mesh type >> 1=individual, 2=resampled32k, 3=resampled164k
  %sinfo = cat_surf_info( files ); 
  %mtype = 1;
 
  for fi = 1:numel(files)
    for si = 1:numel(sides)      
      %%
      fname = char(cat_surf_rename(files{fi},'side',sides{si}));
      [pp,ff,ee] = spm_fileparts( fname ); 
      clear St; 
      if exist(fname, 'file')
        switch ee
          case '.gii'
            St = gifti(fname); 
            St = export(St,'patch'); 
            if (data==1 || data==3) && ~isfield(St,'vertices')
              % try to load individual surface
              fname2 = char(cat_surf_rename(fname,'dataname','central','ee','.gii',...
                'pp',strrep(pp,'atlases_surfaces','templates_surfaces')));
              % type
% >>>              
              if exist(fname2,'file')
                St = gifti(fname2);
              else % if normalized
                % load fiting template ...
                fname2 = char(cat_surf_rename(fname,'dataname','central','ee','.gii',...
                  'pp',strrep(pp,'atlases_surfaces','templates_surfaces')));
                St = gifti(fname2);
              end
            end
          case '.annot'
            % allways template data
            if data==1 || data==3
              fname2 = char(cat_surf_rename(fname,'dataname','central','ee','.gii',...
                'pp',strrep(pp,'atlases_surfaces','templates_surfaces')));
              St = gifti(fname2);
              St = export(St,'patch'); 
            end
            if data==2 || data==3
              [cx,cdata,ctab]  = cat_io_FreeSurfer('read_annotation',fname);
              St.facevertexcdata = cdata;
            end
          otherwise
% >>>                        
            %{
            snames = {'.central.','.inflated.','sphere','.hull.','.inner.','.outer.','.white.','.pial.'}; 
            if any( ~cellfun('isempty',strfind(snames,fname)))
            % FreeSurfer surface - no texture
              [St.vertices,St.faces]  = cat_io_FreeSurfer('read_surf',avg);
            
            else
            % FreeSurfer texture (individual data)
              St.
              
              if data==1 || data==3
                % individual surface
                St.vertices = [];
                St.faces    = [];
            end
            %}
            
        end
      else
        error('cat_surf_loadfiles:noFile2','The "%s" side of input file %d "%s" does not exist.',sides{si},fi,files{fi}); 
      end
     
      if data==1 || data==3
        S.faces     = [S.faces;    St.faces + size(S.vertices,1)]; 
        S.vertices  = [S.vertices; St.vertices]; 
      end
      if data==2 || data==3
        S.facevertexcdata = [S.facevertexcdata; St.facevertexcdata];
      end
      sdata = [sdata repmat( uint8(si) , numel(St.facevertexcdata) , 1) ]; 
      odata = [sdata repmat( uint8(fi) , numel(St.facevertexcdata) , 1) ];
    end
  end
  
  %%
  varargout{1} = S; 
  if nargout>1, varargout{2} = sdata; end
  if nargout>2, varargout{3} = odata; end
  if nargout>3 && exist('ctab','var'), varargout{4} = ctab; end
end


        