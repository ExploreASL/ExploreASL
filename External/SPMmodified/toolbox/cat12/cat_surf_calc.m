function varargout = cat_surf_calc(job)
% ______________________________________________________________________
% Surface Calculation Tool - Only batch mode available. 
%
% [Psdata] = cat_surf_smooth(job)
%
% job.cdata      .. cellstr or cell of cellstr for multi-subject
%                   processing
% job.dataname   .. output name (def = 'ouput')
% job.outdir     .. output directory (if empty first subject directory) 
% job.expression .. texture calculation expression 
%                     's1 + s2' for dmtx==0
%                     'mean(S)' for dmtx==1 
% job.dmtx       .. use data matrix
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_surf_calc.m 1535 2019-12-13 15:40:18Z gaser $

  if strcmp(job,'selftest')
    cat_surf_calc_selftest;
  end
  
  if nargin == 1
    def.nproc           = 0; % multiple threads
    def.verb            = 0; % dispaly something
    def.lazy            = 0; % do not process anything if output exist (expert)
    def.datahandling    = 1; % 1-subjectwise, 2-datawise
    def.usefsaverage    = 1; % 
    def.assuregifti     = 0; % write gii output
    def.dataname        = 'output'; 
    def.usetexturefield = 0;
    job = cat_io_checkinopt(job,def);
  else
    error('Only batch mode'); 
  end
  
  % prepare output filename
  if iscellstr(job.cdata)
    sinfo = cat_surf_info(job.cdata{1});
  else
    sinfo = cat_surf_info(job.cdata{1}{1});
  end
  
  if strfind(sinfo.side,'mesh')
    if sinfo.resampled_32k
      job.fsaverage = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','mesh.central.freesurfer.gii');  
    else
      job.fsaverage = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','mesh.central.freesurfer.gii');  
    end
  else
    if sinfo.resampled_32k
      job.fsaverage = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','lh.central.freesurfer.gii');  
    else
      job.fsaverage = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.freesurfer.gii');  
    end
  end
  
  if ~isempty(job.outdir{1}), outdir = job.outdir{1}; else outdir=sinfo.pp; end  
  ee = sinfo.ee; if job.assuregifti, ee = '.gii'; end
  
  job.dataname = strrep(job.dataname,'.gii',''); % remove .gii extension

  % single or multi subject calculation
  if iscellstr(job.cdata)
    if isempty(outdir), outdir = fileparts(job.cdata{1}); end

    if job.usetexturefield 
      job.output = char(cat_surf_rename(job.cdata{1},...
        'preside','','pp',outdir,'name','','dataname',job.dataname,'ee',ee));  
    else
      job.output = fullfile(outdir,[job.dataname,ee]); 
    end
    
    % call surfcalc
    if strcmp(strrep(job.expression,' ',''),'s1') % this is just a copy
      copyfile(job.cdata{1},job.output);
    else
      job.verb = 1; 
      surfcalc(job);
    end
    fprintf('Output %s\n',spm_file(job.output,'link','cat_surf_display(''%s'')'));
  
  elseif job.datahandling==1
  % multisubject  
    for si = 1:numel(job.cdata{1}) 
      if ~isempty(outdir)
        soutdir = outdir;
      else
        soutdir = fileparts(job.cdata{1}{si});
      end

      job.output{si} = char(cat_surf_rename(job.cdata{1}{si},...
        'preside','','pp',soutdir,'dataname',job.dataname,'ee',ee));
    end

    % split job and data into separate processes to save computation time
    if job.nproc>0 && (~isfield(job,'process_index'))
      cat_parallelize(job,mfilename,'cdata');
      return
    elseif isfield(job,'printPID') && job.printPID 
      cat_display_matlab_PID
    end 
    
    if job.nproc==0 
      spm_progress_bar('Init',numel(job.cdata{1}),...
        sprintf('Surface Calculator\n%d',numel(job.cdata{1})),'Subjects Completed'); 
    end
    
    for si = 1:numel(job.cdata{1}) % for subjects
      sjob = rmfield(job,'cdata');
      sjob.verb = 0;
      
      % subject data 
      for ti = 1:numel(job.cdata) % for textures
        sjob.cdata{ti} = job.cdata{ti}{si};
      end
      
      sjob.output   = job.output{si};
      try
        if strcmp(strrep(job.expression,' ',''),'s1') % this is just a copy
          copyfile(sjob.cdata{1},job.output{si});
        else
          surfcalc(sjob);
        end
        fprintf('Output %s\n',spm_file(sjob.output,'link','cat_surf_display(''%s'')'));
      catch
        fprintf('Output %s failed\n',sjob.output);
      end
        
      if job.nproc==0 
        spm_progress_bar('Set',si);
      end
    end
    
    if job.nproc==0 
      spm_progress_bar('Clear');
    end
  else
    %%
    for di = 1:numel(job.cdata) 
      if ~isempty(outdir)
        soutdir = outdir;
      else
        soutdir = fileparts(job.cdata{di}{1});
      end
      
      sinfo = cat_surf_info(job.cdata{di});
      names = 'isequaln('; for si=1:numel(sinfo), names = [names '''' sinfo(si).name ''',']; end; names(end) = ')';  %#ok<AGROW>
      if numel(sinfo)>1 && eval(names)
        job.output{di,1} = char(cat_surf_rename(job.cdata{di}{1},...
          'preside','','pp',outdir,'name','','dataname',job.dataname,'ee',ee));  

      else
        job.output{di,1} = char(cat_surf_rename(job.cdata{di}{1},...
          'preside','','pp',soutdir,'name',job.dataname,'ee',ee));
      end
    end
    
    
    % split job and data into separate processes to save computation time
    if job.nproc>0 && (~isfield(job,'process_index'))
      cat_parallelize(job,mfilename,'cdata');
      return
    elseif isfield(job,'printPID') && job.printPID 
      cat_display_matlab_PID
    end 
    
    if job.nproc==0 
      spm_progress_bar('Init',numel(job.cdata{1}),...
        sprintf('Surface Calculator\n%d',numel(job.cdata{1})),'Subjects Completed'); 
    end

    for di = 1:numel(job.cdata) % for subjects
      sjob = rmfield(job,'cdata');
      sjob.verb = 0;
      
      % subject data 
      sjob.cdata  = job.cdata{di};
      sjob.output = job.output{di};
      try
        surfcalc(sjob);
        fprintf('Output %s\n',spm_file(sjob.output,'link','cat_surf_display(''%s'')'));
      catch
        fprintf('Output %s failed\n',sjob.output);
      end
        
      if job.nproc==0 
        spm_progress_bar('Set',di);
      end
    end
    
    if job.nproc==0 
      spm_progress_bar('Clear');
    end
    
  end
  
  
  
  if nargout
    if iscellstr(job.cdata)
      varargout{1}.output = {job.output};
    else % is cell of cellstrings 
      varargout{1}.output = job.output;
    end
  end
end
function surfcalc(job)
    
  def.debug     = cat_get_defaults('extopts.verb')>2;
  def.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.pbar      = 0;
  job = cat_io_checkinopt(job,def);

  %% calculation 
  [sinfo1,S1] = cat_surf_info(job.cdata{1},1);
  sinfo = cat_surf_info(job.cdata);
  
  cdata  = zeros([1,sinfo1(1).ncdata],'single');
  if sinfo1.datatype==3
    vdata = zeros([1,sinfo1(1).nvertices,3],'single'); 
  end

  % work on subsets ("slices") to save memory
  subsetsize = round(10e10 / numel(job.cdata));
  
  if job.verb
    spm_clf('Interactive'); 
    spm_progress_bar('Init',numel(job.cdata),...
      sprintf('Surface Calculator\n%s',job.output),'Input Surfaces Completed'); 
  end
  sdata = struct('dsize',[],'fsize',[],'vsize',[]); 
  for si = 1:ceil(sinfo1(1).nvertices/subsetsize)
    range = [ (si-1) * subsetsize + 1 , si * subsetsize ]; 
    range = min(range,sinfo1(1).nvertices);
    range = range(1):range(2);

    if job.dmtx
      S = zeros(numel(job.cdata),numel(range),1,'single');
    end
    if sinfo1.datatype==3
      V = zeros(numel(job.cdata),numel(range),3,'single');
    end


    %%
    for i=1:numel(job.cdata)
      if sinfo(i).ftype==1 
        GS = gifti(job.cdata{i});
        if isfield(GS,'cdata')
          try
            d  = reshape(GS.cdata,1,sinfo1(1).nvertices);
          catch
            if i>1 && numel(job.cdata{i-1})~=numel(job.cdata{i})
              fprintf('cat_surf_calc:gificdata',...
                'Number of elements must be equal:\n%12d entries in file %s\n%12d entries in file %s\n',...
                numel(job.cdata{i-1}), job.cdata{i-1}, numel(job.cdata{i}), job.cdata{i});
            end
          end
        else
          error('cat_surf_calc:gifticdata',...
            'No texture found in ''s''!',job.cdata{i}); 
        end
        sdata(i).dsize = size(GS.cdata); 
        if sinfo1.datatype==3
          V(i,:,:) = shiftdim(GS.vertices(range,:),-1);
          sdata(i).vsize = size(GS.vertices); 
          sdata(i).fsize = size(GS.faces); 
        end
      else
        d = cat_io_FreeSurfer('read_surf_data',job.cdata{i})';
        sdata(i).dsize = size(d); 
      end
      if i>1
        if any(sdata(i).dsize~=sdata(i-1).dsize)
          if sinfo(i).resampled==0
            error('cat_surf_calc:texturesize',...
              'Surface ''s%d'' (%s) does not match previous texture (non-resampled input)!%s',i,job.cdata{i}); 
          else            
            error('cat_surf_calc:texturesize',...
              'Surface ''s%d'' (%s) does not match previous texture!%s',i,job.cdata{i}); 
          end
        end
        if sinfo(i).datatype==3 && ...
          any(sdata(i).vsize~=sdata(i-1).vsize) || any(sdata(i).fsize~=sdata(i-1).fsize)
            error('cat_surf_calc:meshsize',...
              'Mesh ''s%d'' (%s) does not match to previous mesh!',i,job.cdata{i}); 
        end
      end
      d = d(1,range,1);


      if job.dmtx
        S(i,:) = d; 
      else
        eval(['s',num2str(i),'=d;']);
      end
      
      %% evaluate mesh 
      if sinfo1.datatype==3
        if job.usefsaverage
          CS = gifti(strrep(job.fsaverage,'lh.',sinfo1.side)); 
          vdata(1,range,:) = CS.vertices;  
        else
          vdata(1,range,:) = mean(V,1);
        end
      end
      
      if job.verb
        spm_progress_bar('Set',(si-1)*numel(job.cdata)/subsetsize + i);
      end      
    end

    %% evaluate texture
    try 
      eval(['cdata(range) = ' job.expression ';']);
    catch %#ok<CTCH>
      l = lasterror; %#ok<*LERR>
      error('%s\nCan''t evaluate "%s".',l.message,job.expression);
    end



    %spm_progress_bar('Set',si);
  end
  

  

  %% save texture
  ppn = fileparts(job.output); if ~exist(ppn,'dir'), mkdir(ppn); end
  if sinfo1.datatype==3 || strcmp(job.output(end-3:end),'.gii')
    if ~strcmp(job.output(end-3:end),'.gii'), job.output = [job.output '.gii']; end
    if sinfo1.datatype==3
      save(gifti(struct('vertices',shiftdim(vdata),'faces',S1{1}.faces,'cdata',cdata')),job.output,'Base64Binary');
    else
      save(gifti(struct('cdata',cdata)),job.output,'Base64Binary');
    end
  else
    cat_io_FreeSurfer('write_surf_data',job.output,cdata');
  end

  if job.verb
    spm_progress_bar('Clear');
  end
end

function cat_surf_calc_selftest
  % creation of a test directory with syntect (resampled) simple surface (cubes)
  %
  %   s15mm.[rl]h.thickness.resampled.C01.gii  % full gifti resampled
  %   s15mm.[rl]h.thickness.resampled.C02.gii  % full gifti resampled
  %   [rl]h.thickness.C01                      % FS texture > error 
  %   [rl]h.thickness.resampled.C01            % FS resampled texture 
  %
  %   [rl]h/beta0815.gii                       %  
  %   [rl]h/mask0815.gii                       %
  %
  % GUI - Parameter structure with(out) side handling 
  %
  % job cases:
  % - standard imcalc with a mix of GII and FS surfases 
  %   - different outputs (PATH,FS-FILE,GII-File)
  %   - only texture, both, use AVG, ...
  %   - different expressions with and without datamatrix  
  % 

end

























  