function varargout = cat_surf_parameters(job)
% ______________________________________________________________________
%
% cat_surf_parameters to extract surface parameters such as
% gyrification and cortical complexity.
%
%   varargout = cat_surf_parameters(job)
%
%   job.
%    .data_surf .. input data
%    .nproc     .. parallel processing (default 0)
%    .verb      .. verbose output (default cat_get_defaults('extopts.verb'))
%    .lazy      .. avoid reprocess of exist results (default 0)
%    .debug     .. (default cat_get_defaults('extopts.verb')>2)
%    = measures = 
%    .GI        .. estimate absolute mean curvature (default 0)
%    .FD        .. estimate fractal dimension (Yotter:2012; default 0)
%    .SD        .. estimate sulcal depth (default 0)
%    = experimental measures = (only cat_get_defaults('extopts.expertgui')>1)
%    .GIL       .. estimate Laplacian-based gyrification index 
%                  is a numeric in case of default users (default 0)
%                  is a structure in case of expert users 
%    .area      .. estimate area (not implemented; default 0)
%    .surfaces  .. further cortical surfaces
%     .IS       .. create inner surface (default 0)
%     .OS       .. create outer surface (default 0)
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id: cat_surf_parameters.m 1590 2020-03-23 07:38:58Z dahnke $

  SVNid = '$Rev: 1590 $';
 
  if nargin == 1
    P  = char(job.data_surf);
  else
    error('Not enough parameters.');
  end

  try
    if cat_io_matlabversion>20161, rng(0); else, randn('state',0); rand('state',0); end
  end
  
  % default structure
  def.fsavgDir    = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.trerr       = 0; % display errors
  def.nproc       = 0; % parallel processing
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 0; % do not reprocess exist results
  def.debug       = cat_get_defaults('extopts.verb')>2;
  % output parameter of validated measures 
  def.GI          = 0; % estimate absolute mean curvature
  def.FD          = 0; % estimate fractal dimension (Yotter:2012)
  def.SD          = 0; % estimate sulcal depth
  def.tGI         = 0; % Toro's GI
  % implemented but under test
  def.area        = 0; % estimate area
  def.gmv         = 0; % cortical volume
  % normalization of measures by scaling (in development)
  def.norm        = 31; % cat_surf_scaling/normalization [12-affine,1-radius,11-hullradius,2-area,21-hullarea,31-hullvolume];
  def.normprefix  = 'n';
  % further thickness measures by estimating the IS and OS by Tnormal that 
  % result in Tpbt = Tnormal = Tnear
  def.thickness.Tfs   = 0; % Freesurfer thickness metric = mean([ Tnear(IS) Tnear(OS) ],2) 
  def.thickness.Tmin  = 0; % mininmal   thickness metric = min([  Tnear(IS) Tnear(OS) ],2)
  def.thickness.Tmax  = 0; % maximal    thickness metric = max([  Tnear(IS) Tnear(OS) ],2) that is only 
  % experimental measures (cat_get_defaults('extopts.expertgui'))
  % def.GIL         = 0; % defined below due to numeric/structure definion
  % further surfaces
  def.surfaces.IS = 0; % create inner surface 
  def.surfaces.OS = 0; % create outer surface
  
  job = cat_io_checkinopt(job,def);
  if isfield(job,'Tfs'), job.thickness.Tfs = job.Tfs; end
  
  % estimate Laplacian-based gyrification index (including inward, outward, and generalized GI) 
  if ~isfield(job,'GIL'), job.GIL = 0; end
  
  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && numel(job.data_surf)>1 % no parallel processsing for just one file 
  
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    else
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    end
    return
  end

  % new banner
  if isfield(job,'process_index') && job.verb, spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Processed surfaces','Surfaces Completed');
  
  % just a counter for the progress bar
  sides     = {'l','r'};
  measuresn = job.GI + job.FD + job.SD + job.surfaces.IS + job.surfaces.OS + ...
    ( ( isnumeric(job.GIL) && job.GIL ) || ( isstruct(job.GIL) && job.GIL.GIL ) ); 
  measuresn = measuresn * numel(sides);
  
  
  % main loop
  for i=1:size(P,1)
    measuresi = 0; 
    pstr = sprintf(sprintf('%% %ds',max(10,round(log10(size(P,1))+3) * 2)),sprintf('%d/%d) ',i,size(P,1)));  
    nstr = repmat(' ',1,numel(pstr)); 
    
    % go through left and right hemisphere
    %try 
      for si=1:numel(sides)
        %% file names
        if si == 1 % lh
          Pname = deblank(P(i,:));
        else % rh
          Pname = cat_surf_rename(deblank(P(i,:)),'side','rh');
          Pname = Pname{1};
        end
        [pp,ff,ex] = spm_fileparts(Pname);
        name       = [ff ex];
        Psname     = fullfile( pp , strrep([ff ex],'central','scentral') ); 
        Pname2     = fullfile( pp , strrep([ff ex],'central','central2') );   % temporary surface with noise for failed save sulcal depth estimation 
        % dependencies
        PGI     = { 
          fullfile(pp,strrep(ff,'central','gyrification'));                    % unscaled absolute mean curvature   
          fullfile(pp,strrep(ff,'central',[job.normprefix 'gyrification']));   % scaled absolute mean curvature
          };
        PSD     = { 
          fullfile(pp,strrep(ff,'central','sqrtsulc'));
          fullfile(pp,strrep(ff,'central',[job.normprefix 'sqrtsulc']));
          };
        PFD     = fullfile(pp,strrep(ff,'central','fractaldimension'));
        % new affine normalized measures under test
        PtGI    = {
          fullfile(pp,strrep(ff,'central','toroGI'));
          fullfile(pp,strrep(ff,'central',[job.normprefix 'toroGI']));
          fullfile(pp,strrep(ff,'central',[job.normprefix 'xtoroGI']));
          };
        %
        % still under test
        Parea   = {
          fullfile(pp,strrep(ff,'central','area'));
          fullfile(pp,strrep(ff,'central',[job.normprefix 'area']));
          };
        Pgmv    = {
          fullfile(pp,strrep(ff,'central','gmv'));
          fullfile(pp,strrep(ff,'central',[job.normprefix 'gmv']));
          };
        % new experimental GIs
        PiGI    = fullfile(pp,strrep(ff,'central','inwardGI'));            
        PoGI    = fullfile(pp,strrep(ff,'central','outwardGI'));            
        PgGI    = fullfile(pp,strrep(ff,'central','generalizedGI'));        
        % thickness measures
        Ptfs    = fullfile(pp,strrep(ff,'central','thicknessfs'));      
        Ptmin   = fullfile(pp,strrep(ff,'central','thicknessmin'));      
        Ptmax   = fullfile(pp,strrep(ff,'central','thicknessmax'));      
        % other surfaces 
        PIS     = fullfile(pp,strrep([ff ex],'central','white'));         
        POS     = fullfile(pp,strrep([ff ex],'central','pial'));         
        Psphere = fullfile(pp,strrep(name,'central','sphere'));

        if job.verb && si==1, fprintf('\n%sExtract parameters for %s\n',pstr,Pname); end

        
        
        % normalization by affine information?
        if any( [ job.GI>1  && cat_io_rerun(PGI{2},Pname) ...
                  job.SD>1  && cat_io_rerun(PSD{2},Pname) ...
                  (job.tGI==2 || job.tGI==4) && cat_io_rerun(PtGI{2},Pname) ...
                  (job.tGI==3 || job.tGI==4) && cat_io_rerun(PtGI{3},Pname) ...
                  ] )
          cat_surf_scaling(struct('file',Pname,'norm',job.norm,'fname',Psname));
        end


        
        
        if job.area
          %% local surface area by nearest neighbor approach (see Winkler 2017)
          for areai = setdiff( (1:2) .* (job.area==[1 2] | job.area==[3 3]) ,0) 
            stime = clock; 
            if ~cat_io_rerun(Parea{areai},Pname) && job.lazy  
              if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Parea{areai},'link','cat_surf_display(''%s'')')); end
            else 
              Si   = gifti(Pname); 
              area = cat_surf_fun('area',Si); % in mm2
              cat_io_FreeSurfer('write_surf_data',Parea{areai},area); 
              if job.area>0
                cat_io_FreeSurfer('write_surf_data',Parea{areai},area / sum(area) * 1000 ); 
              end
              clear area Si;  
              if job.verb
                fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Parea{areai},'link','cat_surf_display(''%s'')')); 
              end
            end
            
            if nargout==1 && (areai==1 || areai==3), varargout{1}.([sides{si} 'Parea' ]){i} = Parea{areai}; end
            if nargout==1 && (areai==2 || areai==3), varargout{1}.([sides{si} 'Pareas']){i} = Parea{areai}; end
            measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
          end
        end

        
        
        
        if job.gmv
          %% local volume 
          for gmvi = setdiff( (1:2) .* (job.gmv==[1 2] | job.gmv==[3 3]) ,0)
            stime = clock; 
            if ~cat_io_rerun(Pgmv{gmvi},Pname) && job.lazy  
              if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Pgmv{gmvi},'link','cat_surf_display(''%s'')')); end
            else 
              Sw  = cat_surf_fun('whitevar',Pname);
              Sp  = cat_surf_fun('pialvar',Pname);
              gmv = cat_surf_fun('gmv',Sw,Sp); clear Sw Sp; 
              cat_io_FreeSurfer('write_surf_data',Pgmv{gmvi},gmv); clear gmv; 
              if job.verb
                fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Pgmv{gmvi},'link','cat_surf_display(''%s'')')); 
              end
            end
            if nargout==1 && (gmvi==1 || gmvi==3), varargout{1}.([sides{si} 'Pgmv' ]){i} = Pgmv{gmvi}; end  
            if nargout==1 && (gmvi==2 || gmvi==3), varargout{1}.([sides{si} 'Pgmvs']){i} = Pgmv{gmvi}; end  
            measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
          end
        end
  
        
        
        if job.GI
        %% gyrification index based on absolute mean curvature
          for GIi = setdiff( (1:2) .* (job.GI==[1 2] | job.GI==[3 3]) ,0) 
            if ~cat_io_rerun(PGI{GIi},Pname) && job.lazy  
              if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PGI{GIi},'link','cat_surf_display(''%s'')')); end
            else
              stime = clock; 
              cmd = sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',Pname,PGI{GIi});
              [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
              if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PGI{GIi},'link','cat_surf_display(''%s'')')); end
            end
            if nargout==1 && (GIi==1 || GIi==3), varargout{1}.([sides{si} 'PGI']){i}  = PGI{GIi}; end  
            if nargout==1 && (GIi==2 || GIi==3), varargout{1}.([sides{si} 'PGIs']){i} = PGI{GIi}; end  
          end
          
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        
        
        

        if job.SD
        %% sulcus depth
          for SDi = setdiff( (1:2) .* (job.SD==[1 2] | job.SD==[3 3]) , 0 ) 
            if ~cat_io_rerun(PSD{SDi},Pname) && job.lazy  
              if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PSD{SDi},'link','cat_surf_display(''%s'')')); end
            else
              stime = clock; 
              cmd = sprintf('CAT_SulcusDepth -sqrt "%s" "%s" "%s"',Pname,Psphere,PSD{SDi}); %-sqrt
              try
                [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug*0,job.trerr*0);
              catch
                % catch block that was required for some simulated datasets
                % and can probabely removed in future (RD 202002)
                S = gifti( Pname );

                S.vertices = S.vertices + 0.1 * (rand(size(S.vertices))-0.5); 
                MHS = spm_mesh_smooth(S); s=1;
                S.vertices = [ ...
                  spm_mesh_smooth(MHS,double(S.vertices(:,1)),s) , ...
                  spm_mesh_smooth(MHS,double(S.vertices(:,2)),s) , ...
                  spm_mesh_smooth(MHS,double(S.vertices(:,3)),s) ];
                clear MHS; 

                save( gifti(struct('faces',S.faces,'vertices',S.vertices)),Pname2,'Base64Binary'); clear S; 
                cmd = sprintf('CAT_SulcusDepth -sqrt "%s" "%s" "%s"',Pname2,Psphere,PSD{SDi}); %-sqrt
                delete(Pname2); 

                [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
              end
              if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PSD{SDi},'link','cat_surf_display(''%s'')')); end
            end
            if nargout==1 && (SDi==1 || SDi==3), varargout{1}.([sides{si} 'PSD' ]){i} = PSD{SDi}; end
            if nargout==1 && (SDi==2 || SDi==3), varargout{1}.([sides{si} 'PSDs']){i} = PSD{SDi}; end
          end
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end
        
        
        

        if job.FD
        %% fractal dimension using spherical harmonics
          if ~cat_io_rerun(PFD,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PFD,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cmd = sprintf('CAT_FractalDimension -sphere "%s" -nosmooth "%s" "%s" "%s"',Psphere,Pname,Psphere,PFD);
            [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PFD,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'PFD']){i} = PFD; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        
        
        
        if job.tGI
        %% Toro's gyrification index
          for tGIi = setdiff( (1:3) .* (job.tGI==[1 2 3] | job.tGI==[4 4 4]) , 0 ) 
            if ~cat_io_rerun(PtGI{tGIi},Pname) && job.lazy  
              if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PtGI{tGIi},'link','cat_surf_display(''%s'')')); end
            else
              stime = clock; 
              cmd = sprintf('CAT_DumpSurfaceRatio "%s" "%s" no_normalization ',Pname,PtGI{tGIi});
              [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
              if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PtGI{tGIi},'link','cat_surf_display(''%s'')')); end
            end
            if nargout==1 && (tGIi==1 || tGIi==4), varargout{1}.([sides{si} 'PtGI'  ]){i} = PtGI{tGIi}; end
            if nargout==1 && (tGIi==2 || tGIi==4), varargout{1}.([sides{si} 'PtGIs' ]){i} = PtGI{tGIi}; end
            if nargout==1 && (tGIi==3 || tGIi==4), varargout{1}.([sides{si} 'PtGIsx']){i} = PtGI{tGIi}; end
          end
          
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end
        
        
        
        
        
        %% Developer folding measures
        %  ------------------------------------------------------------------
        %  These approaches are still in development. 
        %  See cat_surf_gyrification for further information.
        %  ------------------------------------------------------------------
        if ( isnumeric(job.GIL) && job.GIL ) || ( isstruct(job.GIL) && job.GIL.GIL )
          %% gyrification index based on laplacian GI
          if isnumeric(job.GIL) % default user mode only support default values  
            GIL    = job.GIL;
            GILjob = struct('verb',job.verb); 
          else
            GIL    = job.GIL.GIL; 
            GILjob = job.GIL; 

            % new experimental GIs
            if isfield(job.GIL,'suffix')
              PiGI    = fullfile(pp,strrep(ff,'central',['inwardGI',job.GIL.suffix]));            
              PoGI    = fullfile(pp,strrep(ff,'central',['outwardGI',job.GIL.suffix]));            
              PgGI    = fullfile(pp,strrep(ff,'central',['generalizedGI',job.GIL.suffix])); 
              Phull   = fullfile(pp,strrep([ff '.gii'],'central',['hull',job.GIL.suffix])); 
              Pcore   = fullfile(pp,strrep([ff '.gii'],'central',['core',job.GIL.suffix])); 
            else
              Phull  = fullfile(pp,strrep([ff '.gii'],'central','hull')); 
              Pcore  = fullfile(pp,strrep([ff '.gii'],'central','core'));
            end
          end

          % run GI estimation if ~lazy or any output does not exist
          if job.lazy==0 || ...
            ( any(GIL==[1,4]) && ~exist(PiGI,'file') ) || ...
            ( any(GIL==[2,4]) && ~exist(PoGI,'file') ) || ...
            ( any(GIL==[3,4]) && ~exist(PgGI,'file') )


            % do not display GI processing details while you write into a file!
            if isfield(job,'process_index'), GILjob.verb = 0; end


            % process data
            stime = clock; 
            first = 1; 
            PGIL  = cat_surf_gyrification(Pname,GILjob);
          else
            first = 2;
            PGIL  = {PiGI,PoGI,PgGI};
            if ~exist(PiGI,'file'), PGIL(1) = ''; end
            if ~exist(PoGI,'file'), PGIL(2) = ''; end
            if ~exist(PgGI,'file'), PGIL(3) = ''; end
          end
          if nargout && exist('Phull','var')
            varargout{1}.([sides{si} 'Phull']){i} = Phull; 
            varargout{1}.([sides{si} 'Pcore']){i} = Pcore; 
          end
          
          %%
          type  = 'iog'; 
          for gi=1:numel(PGIL)
            if job.verb && ~isempty(PGIL{1})
              if first==1 
                fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PGIL{gi},'link','cat_surf_display(''%s'')')); 
                first = 0;
              elseif first==2
                fprintf('%sexist - Display %s\n',nstr,spm_file(PGIL{gi},'link','cat_surf_display(''%s'')')); 
              else
                fprintf('%s      - Display %s\n',nstr,spm_file(PGIL{gi},'link','cat_surf_display(''%s'')')); 
              end  
              if nargout==1,  varargout{1}.([sides{si} 'P' type(gi) 'GI']){i} = PGIL{gi}; end
            end
          end

          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end




        %% ----------------------------------------------------------------------
        %  Further thickness measures.
        %  ----------------------------------------------------------------------
        %existIOS = [exist(PIS,'file') exist(POS,'file')]; 

        if job.thickness.Tfs
          if ~cat_io_rerun(Ptfs,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Ptfs,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('Tfs',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Ptfs,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'Tfs']){i} = Ptfs; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        if job.thickness.Tmin
          if ~cat_io_rerun(Ptmin,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Ptmin,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('Tmin',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Ptmin,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'Tmin']){i} = Ptmin; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        if job.thickness.Tmax
          if ~cat_io_rerun(Ptmax,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Ptmax,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('Tmax',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Ptmax,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'Tmax']){i} = Ptmax; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        % delete temporary surface files
        %if existIOS
          % hier muesste ich noch die IS und OS ggf. aufraeumen
        %end

        %% ----------------------------------------------------------------------
        %  No measures, but I do not want another script. However, this leads
        %  to problems in batch processing, e.g. to resample and smooth the 
        %  results that are surfaces rather than textures (RD20190408). 
        %  ----------------------------------------------------------------------
        if job.surfaces.IS
          if ~cat_io_rerun(PIS,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PIS,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('white',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PIS,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'PIS']){i} = PIS; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        if job.surfaces.OS
          if ~cat_io_rerun(PIS,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(POS,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('pial',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(POS,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'POS']){i} = POS; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end


        if exist(Psname,'file'), delete(Psname); end
      end
      spm_progress_bar('Set',i);




      if isfield(job,'process_index') && job.verb
        fprintf('%sDone\n',nstr);
      end  
  %  catch
  %    if job.verb, cat_io_cprintf('err','%sERROR - Check data of %s\n',nstr,spm_file(P(i,:),'link','cat_surf_display(''%s'')')); end
  %  end 
  end
  if isfield(job,'process_index') && job.verb
    fprintf('\nDone\n');
  end  

  spm_progress_bar('Clear');  
  
  if nargout && ~exist('varargout','var'),  varargout{1} = struct(''); end
  
  % remove files that do not exist
  varargout{1} = cat_io_checkdepfiles( varargout{1} );
end
