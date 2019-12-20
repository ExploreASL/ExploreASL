function varargout = cat_vol_headtrimming(job)
% ______________________________________________________________________
% Remove air around the head and convert the image data type to save disk-
% space but also to reduce memory-space and load/save times. Uses 99.99% of 
% the main intensity histogram to avoid problems due to outliers. Although 
% the internal scaling supports a relative high accuracy for the limited 
% number of bits, special values such as NAN and INF will be lost!
%
%   varargout = cat_vol_headtrimming(job)
%
%   job
%    .images .. filenames as cell of cellstr
%    .resdir .. result directory (char, default same as input files)
%    .pefix  .. filename prefix (char, default = 'trimmed_')
%    .addvox .. additional voxels around the box (default = 2); 
%    .pth    .. percentual threshold to estimate the box (0.1);
%    .avg    .. create the box on the the averag of ..
%                 avg = 0 .. use only first image (default)
%                 avg = 1 .. all
%                 avg > 1 .. use image 1 to avg
%    .verb   .. be verbose (default = 1) 
%    .ctype  .. 'uint16';
%    .range  .. 99.99;
%
% ______________________________________________________________________
%
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_vol_headtrimming.m 1339 2018-07-25 13:35:02Z gaser $

  SVNid = '$Rev: 1339 $';
  
  def.images  = {{}}; 
  def.resdir  = '';
  def.prefix  = 'trimmed_';
  def.postfix = '';
  def.addvox  = 2; 
  def.pth     = 0.4;
  def.avg     = 1; 
  def.verb    = 1; 
  def.open    = 2; 
  def.ctype   = 'uint16';
  def.range   = 99.99;
  def.range1  = 90;
  def.returnOnlyFilename  = 0;
  def.process_index = 1;
  job = cat_io_checkinopt(job,def);

  % choose prefix
  for di=1:numel(job.images)
    for si=1:numel(job.images{di})
      [pp,ff,ee,dd]        = spm_fileparts(job.images{di}{si});
      varargout{1}{di}{si} = fullfile(pp,[job.prefix ff job.postfix ee dd]);
    end
  end
  if job.returnOnlyFilename, return; end
  
  % check number of images
  if numel(job.images)>1
    nimgs = cellfun('length',job.images);
    if any(nimgs~=mean(nimgs))
      error('cat_vol_headtrimming:imgages',...
      ['The number of images of each set has to be equal where the i-th entry ' ...
       'of each set describes data of the same subject.']);
    end
  end
    
  if isfield(job,'process_index') && job.process_index && job.verb
    spm('FnBanner',mfilename,SVNid); 
  end
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.images{1}),'SANLM-Filtering','Volumes Complete');
  for di = 1:numel(job.images{1})
    if job.verb, fprintf('%58s: ',spm_str_manip(job.images{1}{di},'ra57')); end
    
    %% estimate trimming parameter
    for si = 1:numel(job.images), V(si) = spm_vol(job.images{si}{di}); end
    Y = zeros(V(1).dim,'single');
    for si = 1:max(1,min(numel(V,job.avg)))
      Y = Y + single(spm_read_vols(V(si))); 
    end
    Y = Y ./ max(1,min(numel(V,job.avg)));
    vx_vol  = sqrt(sum(V(1).mat(1:3,1:3).^2)); 
    [Y,hth] = cat_stat_histth(smooth3(Y),job.range1,0); 
    Y = (Y - hth(1)) ./ abs(diff(hth));
    
    % masking
    Yb = zeros(size(Y),'single'); Yb(2:end-1,2:end-1,2:end-1) = Y(2:end-1,2:end-1,2:end-1);
    Yb = smooth3(Yb)>job.pth; 
    Yb = cat_vol_morph(Yb,'do',job.open,vx_vol); 
    Yb = cat_vol_morph(Yb,'l',[10 0.1]); 
    [Yt,redB] = cat_vol_resize(Y,'reduceBrain',vx_vol,job.addvox,Yb); clear Yt; 
    
    % prepare update of AC orientation
    mati  = spm_imatrix(V(1).mat); mati(1:3) = mati(1:3) + mati(7:9).*redB.BB(1:2:end);
    
    %% trimming
    Vo = V; 
    for dxi = 1:numel(V)
      % optimize tissue intensity range
      if job.ctype
        job2 = struct('data',job.images{dxi}{di},'verb',0,'ctype',job.ctype,...
                      'range',job.range,'prefix',job.prefix,'postfix',job.postfix);
        P = cat_io_volctype(job2);
      else
        filecopy(job.images{dxi}{di},evarargout{1}{di}{si});
        P = varargout{1}{di}{si};
      end
      spm_progress_bar('Set',di + dxi/numel(V));
    
      Vo(dxi) = spm_vol(P{1}); 
      Yo = single(spm_read_vols(Vo(dxi)));
      Yo = cat_vol_resize(Yo,'reduceBrain',vx_vol,job.addvox,Yb); 
      Vo(dxi).mat = spm_matrix(mati);
      Vo(dxi).dim = redB.sizeTr;
      if exist(Vo(dxi).fname,'file'), delete(Vo(dxi).fname); end % delete required in case of smaller file size!
      spm_write_vol(Vo(dxi),Yo);
    end
    %%
    if job.verb, fprintf('saved %5.2f%%\n',100 - 100*((prod(redB.sizeTr) ./ prod(redB.sizeT)))); end
    spm_progress_bar('Set',di);
  end
  spm_progress_bar('Clear');
end