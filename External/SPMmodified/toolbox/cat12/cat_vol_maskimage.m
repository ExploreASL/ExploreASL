function Po = cat_vol_maskimage(job)
%
%
% Po = cat_vol_maskimage(job)
%
% job
%   .data    .. original images
%   .mask    .. lesion/brain mask images
%   .prefix  .. filename prefix 'masked_'
%

  def.returnOnlyFilename = 0; 
  def.prefix = 'masked_';
  def.mask   = {''};
  def.bmask  = {''};
  def.recalc = 0; 
  job = cat_io_checkinopt(job,def);

  job.data  = cellstr(job.data);
  job.mask  = cellstr(job.mask);
  job.bmask = cellstr(job.bmask);
  
  % convert to realy empty strings
  if numel(job.mask)==1  && isempty(job.mask{1}),  job.mask  = {}; end
  if numel(job.bmask)==1 && isempty(job.bmask{1}), job.bmask = {}; end
  
  % create output structure
  Po = cell(numel(job.data),1); 

  % SPM processbar
  if isfield(job,'process_index') && job.verb, spm('FnBanner',mfilename,SVNid); end
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data),'Image masking','Volumes Complete');
  
  % set filenames and remove old files
  for di=1:numel(job.data)
    % set filename
    [pp,ff,ee] = spm_fileparts(job.data{di}); 
    Po{di} = fullfile(pp,[job.prefix ff ee]);
     
    % for GUI output
    if job.returnOnlyFilename
      continue; 
    elseif exist(Po{di},'file') && job.recalc
      delete(Po{di}); 
    end
  end
  if job.returnOnlyFilename
    return
  end
  
  % error handling in case of missmatching number of files
  if numel(job.data)==1 && numel(job.mask)>1 && numel(job.bmask)==0
    % ok, multiple masks
  elseif( numel(job.data)~=numel(job.mask)  && numel(job.mask)>1  ) || ...
     ( numel(job.data)~=numel(job.bmask) && numel(job.bmask)>1 )
    error('SPM:CAT:cat_vo_maskimgage','Number of images and mask/brainmask images has to be equal or zero');
    
  elseif numel(job.mask)==0 && numel(job.bmask)==0
    for di=1:numel(job.data)
      fprintf('No masks. Just copy "%s"\n',job.data{di});
      [pp,ff,ee] = spm_fileparts(job.data{di}); 
      copyfile(fullfile(pp,[ff ee]),Po{di})
    end
    %fprintf('No mask images. Nothing useful to do!\n');
    return
  end
  
  % mask data
  for di=1:numel(job.data)
    % load image header to create output structure
    if exist(Po{di},'file') % if the output file already exist that use it to add a new mask!
      Vi = spm_vol(Po{di});
    else    
      Vi = spm_vol(job.data{di});
    end
    fprintf('Process "%s"\n',job.data{di});
    
    Vo = Vi; 
    Vo.fname = Po{di}; 

    % load other images and use imcalc to mask the images
    imcalcopt = struct('verb',0);
    if numel(job.data)==1 && numel(job.mask)>0 && numel(job.bmask)==0
      Vm = spm_vol(char(job.mask));
      cat_vol_imcalc([Vi;Vm],Vo,sprintf('i1 %s',...
        sprintf(' .* ( i%d<0.5 ) ',(1:numel(job.mask))+1)),imcalcopt); 
    elseif numel(job.data)==1 && numel(job.mask)>0 && numel(job.bmask)==0
      Vm = spm_vol(char(job.mask));
      Vb = spm_vol(char(job.bmask));
      cat_vol_imcalc([Vi;Vm;Vb],Vo,sprintf('i1 %s %s',...
        sprintf(' .* ( i%d<0.5 ) ',(1:numel(job.mask))+1),...
        sprintf(' .* ( i%d<0.5 ) ',(1:numel(job.bmask))+1)),imcalcopt); 
    elseif numel(job.mask)>0 && numel(job.bmask)==0
      Vm = spm_vol(job.mask{di});
      cat_vol_imcalc([Vi,Vm],Vo,'i1 .* ( i2<0.5 )',imcalcopt); 
    elseif numel(job.mask)==0 && numel(job.bmask)>0
      Vb = spm_vol(job.bmask{di});
      cat_vol_imcalc([Vi,Vb],Vo,'i1 .* ( i2>0.5 )',imcalcopt); 
    elseif numel(job.mask)>0 && numel(job.bmask)>0
      Vm = spm_vol(job.mask{di});
      Vb = spm_vol(job.bmask{di});
      cat_vol_imcalc([Vi,Vm,Vb],Vo,'i1 .* ( i2<0.5 | i3>0.5 )',imcalcopt); 
    end   
    
    spm_progress_bar('Set',di);
  end
  
  if isfield(job,'process_index') && job.verb, fprintf('Done\n'); end
  spm_progress_bar('Clear');
end