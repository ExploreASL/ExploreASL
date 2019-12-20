function varargout = cat_stat_check_cov2(job)
% cat_stat_check_cov to check covariance and image quality across sample.
% Use the CAT GUI or SPM Batch function for calls.
%
% Images have to be in the same orientation with same voxel size and
% dimension (e.g. spatially registered images), whereas surfaces have 
% to be same size (number of vertices).
%
% varargout = cat_stat_check_cov(job)
%
% job           .. spm job structure
%   .data_vol   .. volume input files
%   .data_surf  .. surface input files 
%   .c          ..
%  [.gap]       .. gap between slices (in case of volume input)
%
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id: cat_stat_check_cov2.m 1352 2018-08-10 16:11:27Z dahnke $

%#ok<*AGROW,*ASGLU,*TRYNC,*MINV,*INUSD,*INUSL>
  
  oldfig = findobj('type','figure','number',2); 
  if ~isempty(oldfig); delete(oldfig); end
  
   % create default SPM windows if required        
  if isempty(spm_figure('FindWin','Interactive'))
    spm('createintwin'); 
  end
  if isempty(spm_figure('FindWin','Graphics'))
    spm_figure('Create','Graphics',sprintf('%s: Graphics',spm('version'))); 
  %else
  %  spm_figure('Clear',spm_figure('FindWin','Graphics'))
  end   
  
  cat_io_cprintf('err','\nWARNING: cat_stat_check_cov2 is in an ealy development stage!\n\n')
  
%  cscc           .. cat_stat_check_cov data structure as unique global variable
%   .H            .. object/button handles
%   .pos          .. position values for GUI objects/buttons    
%   .YpY          .. covariance matrix 
% * .mean_cov     .. mean covariance matrix 
%   .files        .. CAT preprocessing files 
%    .data        .. normalized input files of cat_stat_check_cov
%                   (e.g. wmp1, thickness, curv, ...)
%    .org         .. original files used for preprocessing
%    .surf        .. surface files (thickness, mesh)
%    .surfr       .. resampled files (thickness, mesh)
%    .xml         .. XML data of CAT preprocessing
%    .log         .. log-file of CAT preprocessing
%    .pdf         .. pdf-file of CAT preprocessing (report figure) 
%    .jpg         .. pdf-file of CAT preprocessing (report figure) 
% 
%
  global filename ... structure from spm_str_manip with grouped filenames
         ... H        ... main structure with object/button handles
         pos      ... structure with position values for GUI objects/buttons
         YpY mean_cov         ... (sorted) covariance and mean covariance matrix
         mask1d mask2d                  ... masks for files on the trash list 
         trashlist trashhist trashhistf ... index list of subjects to remove and undo/redo list
         smask1d smask2d pmask1d pmask2d ...
         ind_sorted                     ... index lists 
         QM QM_names                    ... Quality measures from xml_files
         data_files org_files surf_files xml_files log_files pdf_files ...  different lists of filenames 
         data_array data_array_diff     ... slices of all subjects 
         img                            ... slice image(s) for GUI display
         FS FSi                         ... SPM fontsize 
         mesh_detected isscatter isxml sorted show_name useicons ... binary variables for GUI control 
         mn_data mx_data                ... minimum/maximum value of the mesh
         V Vo                           ... volume header structures
         Vchanged names_changed         ... modified volume header for 4D-structures
         sample protocol protocols      ... sample groups of each scan
         dataprefix                     ... prefix of the input files
         inorm                          ... normalize slice intensity even in normalized data
         n_samples  ...
         figcolor; % X X2 MD MD2; 
          % cbar img_alpha


  % positions & global font size option
  %  ------------------------------------------------------------------------
  trashlist     = []; % start with empty list 
  trashhist     = []; % history of trash operations (for undo)
  trashhistf    = []; % history of trash operations (for redo)
  sorted        = 0;  % show data by file order
  isscatter     = 0;  % active GUI surface plot
  show_name     = 0;  % show filenames in boxplot rather small dots
  inorm         = 1;  % normalize slice intensity even in normalized data
  ws            = spm('Winsize','Graphics');
  FS            = spm('FontSizes');
  FSi           = 8; 
  useicons      = 1; 
  figcolor      = [0.8 0.8 0.8]; % color of the figure

  gpos = get(spm_figure('FindWin','Graphics'),'Position'); 
  if gpos(2)>0 && gpos(2)<50, gpos = gpos(2); else gpos = 10; end
  popb = [0.038 0.035];                                          % size of the small buttons
  popm = 0.780;                                                  % x-position of the control elements
  posp = struct('naviui',0.835,'trashui',0.725,'checkui',0.780); % y-pos of major control elements
  pos = struct(...
    'fig',            [gpos(1) gpos(1) 1.4*ws(3) 1.2*ws(3)],... % figure
    'popup',          [10  10  200       100  ],... % popup in case of closing with non-empty trash list
    ...
    'corr',           [0.045 0.050 0.700 0.820],... % correlation matrix
    'slice',          [0.780 0.060 0.190 0.450],... % image plot
    'surfi',          [0.780 0.050 0.190 0.560],... % image plot
    'cbar',           [0.045 0.950 0.580 0.020],... % colorbar for correlation matrix
    'cbarfix',        [0.657 0.943 0.100 0.030],... % colorbar fix/auto option
    'showtrash',      [0.657 0.913 0.100 0.030],... % colorbar fix/auto option
    ... 
    'boxplot',        [0.100 0.055 0.880 0.915],... % boxplot axis
    'fnamesbox',      [0.830 0.003 0.160 0.038],... % show filenames in boxplot 
    ...
    'close',          [0.775 0.935 0.100 0.040],... % close button
    'help',           [0.875 0.935 0.100 0.040],... % help button
    ...
    'sort',           [0.772 0.880 0.110 0.050],... % list to use ordered matrix or Maha-distance 
    'boxp',           [0.872 0.880 0.110 0.050],... % list to display different variables as boxplot
    'samp',           [0.772 0.615 0.110 0.055],... % list to use ordered matrix or Maha-distance 
    'prot',           [0.872 0.615 0.110 0.055],... % list to display different variables as boxplot
    ...
    'alphabox',       [0.775 -0.001 0.200 0.030],... % show filenames in boxplot 
    'sslider',        [0.780 0.030 0.193 0.040],... % slider for z-slice  
    ...
    ... == navigation unit ==
    'scSelect',       [popm+popb(1)*0 posp.naviui popb],... % select (default) 
    'scZoomReset',    [popm+popb(1)*1 posp.naviui popb],... % standard zoom
    'scZoomIn',       [popm+popb(1)*2 posp.naviui popb],... % zoom in 
    'scZoomOut',      [popm+popb(1)*3 posp.naviui popb],... % zoom out
    'scPan',          [popm+popb(1)*4 posp.naviui popb],... % pan (moving hand)
    ...
    ... == tashlist unit ==
    'newtrash',       [popm+popb(1)*0 posp.trashui popb],... % new trash list
    'disptrash',      [popm+popb(1)*1 posp.trashui popb],... % print trash list
    'trash',          [popm+popb(1)*2 posp.trashui popb],... % add data to trash list
    'detrash',        [popm+popb(1)*3 posp.trashui popb],... % remove data from trash list
    'autotrash',      [popm+popb(1)*4 posp.trashui popb],... % button to mark data with low IQR
    ... second row?
    'undo',           [popm+popb(1)*0 posp.trashui-popb(2) popb],... % undo last trash list operation
    'redo',           [popm+popb(1)*1 posp.trashui-popb(2) popb],... % redo last trash list operation
    'trashrow',       [popm+popb(1)*2 posp.trashui-popb(2) popb],... % add data to trash list
    'detrashrow',     [popm+popb(1)*3 posp.trashui-popb(2) popb],... % button to mark data with low IQR
    'ziptrash',       [popm+popb(1)*4 posp.trashui-popb(2) popb],... % pack data from trash list
    ...'unziptrash',     [popm+popb(1)*4 posp.trashui-popb(2) popb],... % remove data from trash list
    ...
    ... == checklist unit ==
    'checkvol',       [popm+popb(1)*0 posp.checkui popb],... % open checkvol 
    'checksurf',      [popm+popb(1)*1 posp.checkui popb],... % open checksurf
    'checklog',       [popm+popb(1)*2 posp.checkui popb],... % open log-txt
    'checkxml',       [popm+popb(1)*3 posp.checkui popb],... % open xml-txt
    'checkpdf',       [popm+popb(1)*4 posp.checkui popb]);   % open pdf in external viewer
    ... 'checklow',     



  if nargin == 0, error('No argument given.'); end



  %% get all filenames from the data_vol/surf input
  %  ------------------------------------------------------------------------
  if isfield(job,'data_vol')
    datafield     = 'data_vol'; 
    datadir       = 'mri'; 
    mesh_detected = 0;
  elseif isfield(job,'data_surf')
    datafield     = 'data_surf'; 
    datadir       = 'surf'; 
    mesh_detected = 1;
  end
  % get all input scans/surfaces
  data_files = {}; 
  for i = 1:numel(job.(datafield))
    data_files = [data_files;job.(datafield){i}];
  end
  % number of samples and scans, trash mask arrays, sample array
  n_subjects  = numel(data_files);
  n_samples   = numel(job.(datafield));
  mask1d      = true(n_subjects,1);           % trash list mask 1D matrix
  mask2d      = true(n_subjects,n_subjects);  % trash list mask 2D matrix
  smask1d     = mask1d;
  smask2d     = mask2d;
  pmask1d     = mask1d;
  pmask2d     = mask2d;
  sample      = [];
  protocol    = []; 
  for i=1:n_samples
    sample = [sample, i*ones(1,size(job.(datafield){i},1))];
  end



  %% get the different files
  %  ------------------------------------------------------------------------
  spm_progress_bar('Init',n_samples,'Search files','subects completed')
  [filenames,fparts] = spm_str_manip(data_files,'trC'); 
  out_files  = data_files;
  org_files  = data_files;
  pdf_files  = data_files;
  xml_files  = data_files; 
  log_files  = data_files;
  surf_files = data_files;
  
  % get real prefix
  % - expect that all files have the same prefix
  % - fparts.s is not enough if all file start similar, eg. mwp1ADNI_*.nii
  [pp,ff,ee] = spm_fileparts(data_files{1});
  [pp1,pp2]  = spm_fileparts(pp); 
  orgfile    = cat_vol_findfiles(pp1,['*' cat_io_strrep(ff,fparts.s,'') '.nii'],struct('depth',1)); 
  if isempty(orgfile)
    orgfile  = cat_vol_findfiles(pp1,['*' cat_io_strrep(ff,fparts.s,'') '.img'],struct('depth',1)); 
  end
  [ppo,ffo,eeo] = spm_fileparts(orgfile{1});
  dataprefix = ff(1:strfind(ff,ffo)-1);
  
  % find files
  for i = 1:numel(data_files)
    [pp,ff,ee] = spm_fileparts(data_files{i});
    [pp1,pp2]  = spm_fileparts(pp); 

    % output files
    out_files{i} = fullfile(pp,ff,ee); 

    % set subdirectories
    if strcmp(pp2,datadir)
      reportdir = 'report';
      surfdir   = 'surf';
    else
      reportdir = ''; 
      surfdir   = '';
    end

    fname = cat_io_strrep(ff,dataprefix,'');
    % set original input files of the CAT preprocessing
    org_files{i} = fullfile(pp1,[fname '.nii']); 
    if ~exist(org_files{i},'file')
      org_files{i} = fullfile(pp1,[cat_io_strrep(ff,dataprefix,'') '.img']);
      if ~exist(org_files{i},'file')
        org_files{i} = ''; 
      end
    end

    % try to find the XML file if not given
    if isempty( char(job.data_xml) ) 
      xml_files{i} = fullfile(pp1,reportdir,...
        ['cat_' fname cat_io_strrep(ee,{'.nii','.img','.gii'},'.xml')]);
      if ~exist(xml_files{i},'file')
        xml_files{i} = ''; 
      end
    else
      xml_files = cellstr(job.data_xml);
    end

    % set report pdf
    pdf_files{i} = fullfile(pp1,reportdir,['catreport_' fname '.pdf']); 
    if ~exist(pdf_files{i},'file')
      pdf_files{i} = ''; 
    end

    % log files
    log_files{i} = fullfile(pp1,reportdir,['catlog_' fname '.txt']);
    if ~exist(log_files{i},'file')
      log_files{i} = ''; 
    end

    % surface files
    surf_files{i} = fullfile(pp1,surfdir,...
      ['lh.thickness.' fname ]);
    if ~exist(surf_files{i},'file')
      surf_files{i} = ''; 
    end

    spm_progress_bar('Set',i);  
  end
  spm_progress_bar('Clear');



  %% load header 
  %  ------------------------------------------------------------------------
  V  = spm_data_hdr_read(char(data_files));
  Vo = spm_data_hdr_read(char(org_files));



  %% load XML data
  %  ------------------------------------------------------------------------
  if isempty( xml_files )
    isxml     = 0;
    QM_names  = '';
  else
    isxml     = 1;

    if size(xml_files,1) ~= n_subjects
      error('XML-files must have the same number as sample size');
    end
   
    QM = nan(n_subjects,4 + (cat_get_defaults('extopts.expertgui')>1) + 2*mesh_detected);
    QM_names = {...
      'Noise rating (NCR)';...
      'Bias Rating (ICR)';...
      'Resoution Rating (RES)';...
      'Weighted overall image quality rating (IQR)';...
      'Protocol IQR difference (IQRp)'; ...
      'Euler number';...
      'Size of topology defects'};
    QM_names = QM_names(1:size(QM,2)); % remove Euler
    
    spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
    for i=1:n_subjects
      % get basename for xml- and data files
      [pth, xml_name]  = fileparts(deblank(xml_files{i}));
      [pth, data_name] = fileparts(V(i).fname);

      % remove leading 'cat_'
      xml_name = xml_name(5:end);

      % check for filenames
      if isempty(strfind(data_name,xml_name)) && ~isempty(xml_name)
        fprintf('Please check file names because of deviating subject names:\n %s vs. %s\n',...
          V(i).fname,xml_files{i});
      end

      xml = cat_io_xml(deblank(xml_files{i}));
      if isfield(xml,'qualityratings')
        QM(i,1:4)  = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.res_RMS xml.qualityratings.IQR];
        RMS(i,1) = xml.qualityratings.res_RMS;
      elseif isfield(xml,'QAM') % also try to use old version
        QM(i,1:4)  = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.qualityratings.res_RMS xml.QAM.QM.res_RMS xml.QAM.QM.IQR];
        RMS(i,1) = xml.QAM.res_RMS;
      else
        RMS(i,1) = nan; 
      end
      if cat_get_defaults('extopts.expertgui')>1
        QM(i,5)  = nan;
      end
      if mesh_detected && isfield(xml.subjectmeasures,'EC_abs')
        QM(i,end-1:end) = [xml.subjectmeasures.EC_abs xml.subjectmeasures.defect_size];
      end
      spm_progress_bar('Set',i);  
    end
    spm_progress_bar('Clear');
    
    % detect protocols by resolution
    pacc = 2; % larger values to detect many protocols, small values to have less
    RMS(isnan(RMS)) = 11; % avoid multiple NaN center
    [protocols,tmp,protocol] = unique(round(RMS*10^pacc)/10^pacc); clear pid RMS
    protocols(protocols==11) = 21; 
    protocols = protocols/2; % average mm rather than rating

    % remove last two columns if EC_abs and defect_size are not defined
    if mesh_detected && all(all(isnan(QM(:,end-1:end))))
      QM = QM(:,end-1:end);
    end

    % added protocol depending QA parameter
    if cat_get_defaults('extopts.expertgui')>1
      [Pth,rth,sq,rths,rthsc,sqs] = cat_tst_qa_cleaner_intern(QM(:,4),struct('site',{protocol},'figure',0));
      QM_names = char([QM_names;{'Protocol IQR difference (PIQR)'}]);
      QM(:,5) = rth(:,3) - QM(:,4);
    end

    % convert marks into rps rating
    mark2rps   = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
    markd2rpsd = @(mark) ( mark*10) + isnan(mark).*mark;
    QM(:,1:4)  = mark2rps(QM(:,1:4));
    if cat_get_defaults('exptops.expertgui')>1
      QM(:,5)    = markd2rpsd(QM(:,5));
    end
  end



  %% add constant to nuisance parameter
  %  ------------------------------------------------------------------------
  G = [];
  if ~isempty(job.c)
    for i=1:numel(job.c)
      G = [G job.c{i}];
    end
    if size(G,1) ~= n_subjects
      G = G';
    end
    G = [ones(n_subjects,1) G];
  end



  %% load surface data, prepare volume data loading
  %  ------------------------------------------------------------------------
  if mesh_detected
    % load surface texture data
    Y = spm_data_read(V)';

    % optional global scaling
    if isfield(job,'gSF')
      for i=1:numel(V)
        Y(:,2) = Y(:,2)*job.gSF(i);
      end
    end

    Y(isnan(Y)) = 0;

    % rescue unscaled data min/max
    mn_data = min(Y(:));
    mx_data = max(Y(:));
    Y = Y - repmat(mean(Y,2), [1 size(Y,2)]);

    % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
    if ~isempty(G) 
      Ymean = repmat(mean(Y), [n_subjects 1]);
      Y = Y - G*(pinv(G)*Y) + Ymean;
    end

    data_array = Y';
    YpY = (Y*Y')/n_subjects;

    % calculate residual mean square of mean adjusted Y
    Y = Y - repmat(mean(Y,1), [n_subjects 1]);
    data_array_diff = Y';

    %MSE = sum(Y.*Y,2);
  else
    if length(V)>1 && any(any(diff(cat(1,V.dim),1,1),1))
      error('images don''t all have same dimensions')
    end
    if max(max(max(abs(diff(cat(3,V.mat),1,3))))) > 1e-8
      error('images don''t all have same orientation & voxel size')
    end

    % consider image aspect ratio
    pos.slice(4) = pos.slice(4) * V(1).dim(2)/V(1).dim(1);

    slices = 1:job.gap:V(1).dim(3);

    dimx = length(1:job.gap:V(1).dim(1));
    dimy = length(1:job.gap:V(1).dim(2));
    Y    = zeros(n_subjects, prod(dimx*dimy));
    YpY  = zeros(n_subjects);
    %MSE  = zeros(n_subjects,1);
    data_array = zeros([V(1).dim(1:2) n_subjects]);



    %-Start progress plot
    %-----------------------------------------------------------------------
    spm_progress_bar('Init',V(1).dim(3),'Check correlation','planes completed')

    for j=slices

      M  = spm_matrix([0 0 j 0 0 0 job.gap job.gap job.gap]);

      for i = 1:n_subjects
        img = spm_slice_vol(V(i),M,[dimx dimy],[1 0]);
        img(isnan(img)) = 0;
        Y(i,:) = img(:);
        if isfield(job,'gSF')
          Y(i,:) = Y(i,:)*job.gSF(i);
        end
      end

      % make sure data is zero mean
      Y = Y - repmat(mean(Y,2), [1 prod(dimx*dimy)]);

      % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
      if ~isempty(G) 
        Ymean = repmat(mean(Y), [n_subjects 1]);
        Y = Y - G*(pinv(G)*Y) + Ymean;
      end

      YpY = YpY + (Y*Y')/n_subjects;

      % calculate residual mean square of mean adjusted Y
      Y = Y - repmat(mean(Y,1), [n_subjects 1]);

      %MSE = MSE + sum(Y.*Y,2);

      spm_progress_bar('Set',j);  

    end

    % correct filenames for 4D data
    if strcmp(V(1).fname, V(2).fname)
      names_changed = 1;
      Vchanged      = V;
      for i=1:n_subjects
        [pth,nam,ext] = spm_fileparts(V(i).fname);
        V(i).fname    = fullfile(pth, [nam sprintf('%04d',i) ext]);
      end
    else
      names_changed = 0;
    end

    spm_progress_bar('Clear');
  end
  clear Y



  %% normalize YpY and estimate mean_cov
  %  ------------------------------------------------------------------------
  d      = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
  dd     = d*d';
  YpY    = YpY./(dd+eps);
  t      = find(abs(YpY) > 1); 
  YpY(t) = YpY(t)./abs(YpY(t));
  YpY(1:n_subjects+1:end) = sign(diag(YpY));
  clear t d dd;

  % extract mean correlation for each data set
  mean_cov = zeros(n_subjects,1);
  for i=1:n_subjects
    cov0        = YpY(i,:);     % extract row for each subject
    cov0(i)     = [];           % remove cov with its own
    mean_cov(i) = mean(cov0);
  end
  clear cov0;



  %% output compressed filenames structure
  %  ------------------------------------------------------------------------
  fprintf('\n');
  fname_m = [];
  fname_tmp = cell(n_samples,1);
  fname_s   = cell(n_samples,1);
  fname_e   = cell(n_samples,1);
  for i=1:n_samples
    [tmp, fname_tmp{i}] = spm_str_manip(char(V(sample == i).fname),'C');
    fname_m = [fname_m; fname_tmp{i}.m];
    fname_s{i} = fname_tmp{i}.s;
    cat_io_cprintf('n','Compressed filenames sample %d: ',i);
    cat_io_cprintf('b',sprintf('%s %s \n',spm_str_manip(tmp,'f120'),...
      repmat('.',1,3*(numel(tmp)>120))));
  end
  filename = struct('s',{fname_s},'e',{fname_e},'m',{fname_m});
  clear fname_e fname_m fname_s fname_tmp tmp



  %% print suspecious files with cov>0.925
  %  ------------------------------------------------------------------------
  YpY_tmp = YpY - tril(YpY);
  [indx, indy] = find(YpY_tmp>0.925);
  [siv,si] = sort(YpY(sub2ind(size(YpY),indx,indy)),'descend');
  % if more than 25% of the data this points to longitudinal data of one subject and no warning will appear
  if ~isempty(indx) && (sqrt(length(indx)) < 0.25*n_subjects)
    fprintf('\nUnusual large correlation (check that subjects are not identical):\n');
    for i=si'
      % exclude diagonal
      if indx(i) ~= indy(i)
        % report file with lower mean correlation first
        if mean_cov(indx(i)) < mean_cov(indy(i))
          cat_io_cprintf('w',sprintf('  %0.4f',YpY(indx(i),indy(i)))); 
          cat_io_cprintf('n',' between ');
          cat_io_cprintf('b',filename.m{indx(i)}); cat_io_cprintf('n',' and ');
          cat_io_cprintf('b',filename.m{indy(i)}); fprintf('\n');
        else
          cat_io_cprintf('w',sprintf('  %0.4f',YpY(indy(i),indx(i)))); 
          cat_io_cprintf('n',' between ');
          cat_io_cprintf('b',filename.m{indy(i)}); cat_io_cprintf('n',' and ');
          cat_io_cprintf('b',filename.m{indx(i)}); fprintf('\n');
        end
      end
    end
  end



  %% sort data and estimate critical files
  %  ------------------------------------------------------------------------
  [mean_cov_sorted, ind_sorted] = sort(mean_cov,'descend');
  threshold_cov      = mean(mean_cov) - 2*std(mean_cov);
  n_thresholded      = find(mean_cov_sorted < threshold_cov,1,'first');
  if ~isempty(n_thresholded)
    fprintf('\nThese data have a mean correlation below 2 standard deviations. \n');
    fprintf('This does not necessarily mean that you have to exclude these data. \n');
    fprintf('However, these data have to be carefully checked:\n');
    for i=n_thresholded:n_subjects
      cat_io_cprintf('r',sprintf('  %0.4f ',mean_cov_sorted(i)));
      cat_io_cprintf('b',V(ind_sorted(i)).fname); fprintf('\n');
    end
  end



  %% output structure
  %  ------------------------------------------------------------------------
  if nargout>0
    varargout{1} = struct('table',[out_files,num2cell(mean_cov)],...
                          'covmat',YpY,...
                          'sorttable',[cellstr(V(ind_sorted).fname),num2cell(mean_cov_sorted)],...
                          'sortcovmat',YpY(ind_sorted,ind_sorted), ...
                          'cov',mean_cov,...
                          'threshold_cov',threshold_cov);
  end
  clear mean_cov_sorted threshold_cov;

  
  
  % check for replicates
  for i=1:n_subjects
    for j=1:n_subjects
      if (i>j) && (mean_cov(i) == mean_cov(j))
        try
          nami = deblank(V(i).fname);
          namj = deblank(V(j).fname);
          if strcmp(nami(end-1:end),',1')
            nami = nami(1:end-2);
          end 
          if strcmp(namj(end-1:end),',1')
            namj = namj(1:end-2);
          end 
          s = unix(['diff ' nami ' ' namj]);
          if (s==0), fprintf(['\nWarning: ' nami ' and ' namj ' are same files?\n']); end
        end
      end
    end
  end
  
  %% create figure
  %  ------------------------------------------------------------------------
  if mesh_detected
    create_figures(job)
  else
    create_figures(job,slices)
  end

 
return
%-End
%-----------------------------------------------------------------------
function create_figures(job,slices)
% -------------------------------------------------------------------------
% create figure
% -------------------------------------------------------------------------

  global H X X2 S V pos isxml QM QM_names MD MD2 protocol protocols figcolor ...
    FS FSi mesh_detected YpY mean_cov ind_sorted sample n_samples sorted
  
  H.graphics = spm_figure('FindWin','Graphics');
  H.figure   = figure(2);
  set(H.figure,'MenuBar','none','Position',pos.fig,...
    'NumberTitle','off','Resize','off','Visible','on','color',figcolor);    

  % move SPM Graphics figure to have no overlap 
  fpos = get(H.figure,'Position');
  gpos = get(H.graphics,'Position');
  set(H.graphics,'Position',[sum(fpos(1:2:3)) + 20 gpos(2:4)]);
  clear fpos gpos;

  if mesh_detected 
    set(H.figure,'Name','CAT Check Covarance: Click in image to display surfaces');
  else
    set(H.figure,'Name','CAT Check Covarance: Click in image to display slices');
  end

  % cursormode and update function
  H.dcm = datacursormode(H.figure);
  set(H.dcm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');

  % create two-area colormap
  colormap([jet(64); gray(64)]);

  % add colorbar without ticks and colorbar image
  H.cbar = axes('Position',pos.cbar,'Parent',H.figure,'Visible','off');
  image(H.cbar,1:64); set(get(H.cbar,'children'),'HitTest','off','Interruptible','off');
  set(H.cbar,'Ytick','','YTickLabel',''); 

  % set correlation matrix image as image
  H.corr = axes('Position',pos.corr,'Parent',H.figure,'Color',figcolor,...
    'Ytick','','YTickLabel','','XTickLabel','','XTick','','XTickLabel','');
  image(H.corr,64 * tril(YpY)); axis image;

  % scatter plot 
  if isxml
    H.scat(1) = axes('Position',pos.corr,'Parent',H.figure, ...
      'visible','off','Box','on','Color',[0.85 0.85 0.85]);

    if cat_get_defaults('extopts.expertgui')>1
      H.scat(2) = axes('Position',pos.corr,'Parent',H.figure, ...
      'visible','off','Box','on','Color',[0.85 0.85 0.85]);
    end
  end
  
  H.slice = axes('Position',pos.slice,'Parent',H.figure, ...
    'visible','off','Box','off','Color',figcolor,...
    'Ytick','','YTickLabel','','XTickLabel','','XTick','','XTickLabel','');
  if ~mesh_detected, axis off; end
  
  % add button for closing all windows
  H.close = uicontrol(H.figure,...
    'Units','normalized','position',pos.close,'Style','Pushbutton','callback',@closeWindows,...
    'string','Close','ToolTipString','Close windows','FontSize',FS(FSi),'ForegroundColor',[0.8 0 0]);

  % add button to open the help HTML
  H.help = uicontrol(H.figure,...
    'Units','normalized','position',pos.help,'Style','Pushbutton',...
    'string','Help','ToolTipString','Open help window','ForegroundColor',[0 0 0.8],'FontSize',FS(FSi),...
    'callback',['spm_help(''!Disp'',''/Users/dahnke/Documents/MATLAB/' ...
      'spm12/toolbox/cat12/html/cat_methods_QA.html'','''',H.graphics);']);


  % create popoup menu for SPM grafix window
  if isxml

    % estimate Mahalanobis distance between mean corr. and weighted overall quality
    X  = [mean_cov, QM(:,4)]; X(isnan(X)) = 0;  % mean correlation and IQR
    S  = cov(X);
    mu = mean(X);
    MD = (X-repmat(mu,[length(X),1]))*inv(S)*(X-repmat(mu,[length(X),1]))'; 
    MD = diag(MD);

    if cat_get_defaults('extopts.expertgui')>1
      X2  = [mean_cov, QM(:,5)]; X2(isnan(X2)) = 0; 
      S2  = cov(X2);
      mu2 = mean(X2);
      MD2 = (X2-repmat(mu2,[length(X2),1]))*inv(S2)*(X2-repmat(mu2,[length(X2),1]))';
      MD2 = diag(MD2);
    end  

    str  = { 'Boxplot...','Mean correlation',QM_names,'Mahalanobis distance'};
    tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation  ', 1} }; 
    for qmi = 1:size(QM,2)
      tmp = [ tmp , { {@show_mean_boxplot, QM(:,qmi), QM_names(qmi,:), -1} }]; 
    end
    tmp = [ tmp , { {@show_mean_boxplot, MD, 'Mahalanobis distance  ', -1} }];
  else
    str  = { 'Boxplot...','Mean correlation'};
    tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation  ', 1} };
  end
  if isxml && cat_get_defaults('extopts.expertgui')>1
      str = char([cellstr(str),{'Mahalanobis distance IQRp'}]);
      tmp = [ tmp , ...
             {{@show_mean_boxplot, MD2, 'Mahalanobis distance (IQRp)  ',2}}];  
  end
  H.boxp = uicontrol(H.figure,...
    'Units','normalized','position',pos.boxp,'Style','PopUp','callback','spm(''PopUpCB'',gcbo)',...
    'string',str,'ToolTipString','Display boxplot','FontSize',FS(FSi),'UserData',tmp);


  % create popoup menu for main check_cov window
  if isxml
    str  = { 'Image...','Mean Correlation: Order by selected filenames', ...
             'Mean Correlation: Sorted by mean correlation','Mahalanobis distance'};
    tmp  = { {@show_matrix, YpY, 0},...
             {@show_matrix, YpY(ind_sorted,ind_sorted), 1},...
             {@show_mahalanobis, X, MD, 1}};
    if cat_get_defaults('extopts.expertgui')>1
      str = char([cellstr(str),{'Mahalanobis distance IQRp'}]);
      tmp = [ tmp , ...
             {{@show_mahalanobis, X2, MD2, 2}}]; 
    end
  else
    str  = { 'Correlation matrix...','Order by selected filename',...
             'Sorted by mean correlation'};
    tmp  = { {@show_matrix, YpY, 0},...
             {@show_matrix, YpY(ind_sorted,ind_sorted), 1} };
  end

  H.sort = uicontrol(H.figure,...
    'Units','normalized','position',pos.sort,'Style','PopUp','UserData',tmp,...
    'callback','spm(''PopUpCB'',gcbo)','string',str,'ToolTipString','Sort matrix','FontSize',FS(FSi));

  onoff = {'on','off'};
 
  % choose only one sample for display
  str  = { 'Sample..',sprintf('full (%d)',numel(sample))};  tmp  = { {@show_sample, 0} }; 
  for i=1:n_samples, 
    str = [str,sprintf('S%d (%d)',i,sum(sample==i))]; 
    tmp = [ tmp , {{@show_sample, i}} ]; 
  end
  H.samp = uicontrol(H.figure,...
    'Units','normalized','position',pos.samp,'Style','PopUp','UserData',tmp,'enable',onoff{1 + (n_samples==1)},...
    'callback','spm(''PopUpCB'',gcbo)','string',str,'ToolTipString','Sort matrix','FontSize',FS(FSi));

  % choose center 
  %protocols   = 1:max(protocol); % only protocol ids
  str  = { 'Protocol..',sprintf('all (%d)',numel(protocol))};  tmp  = { {@show_protocol, 0} }; 
  for i=1:numel(protocols), 
    %str = [str,sprintf('P%03d',protocols(i))]; % only protocol ids
    str = [str,sprintf('P %5.2f (%d)',protocols(i),sum(protocol==i))]; 
    tmp = [ tmp , {{@show_protocol, i}} ]; 
  end
  H.prot = uicontrol(H.figure,...
    'Units','normalized','position',pos.prot,'Style','PopUp','UserData',tmp,'enable',onoff{1 + (numel(protocols)==1)},...
    'callback','spm(''PopUpCB'',gcbo)','string',str,'ToolTipString','Sort matrix','FontSize',FS(FSi));
  
  
  H.alphabox = uicontrol(H.figure,...
    'Units','normalized','position',pos.alphabox,'Style','CheckBox','callback',@update_alpha,...
    'string','Colorize diff. to sample mean','Value',1,...
    'ToolTipString','Colorize difference to sample mean (pos=green;neg=red)',...
    'Visible','off','BackgroundColor',figcolor,'FontSize',FS(FSi));

  H.cbarfix = uicontrol(H.figure,...
    'Units','normalized','Style','CheckBox','position',pos.cbarfix,'callback',{@checkbox_cbarfix},...
    'string','Fixed range','ToolTipString','Switch between fixed and auto-scaled colorbar',...
    'Value',1,'BackgroundColor',figcolor,'FontSize',FS(FSi));

  H.showtrash = uicontrol(H.figure,...
    'Units','normalized','Style','CheckBox','position',pos.showtrash,'callback',{@checkbox_showtrash},...
    'string','Show Trash','ToolTipString','Show trashed records',...
    'Value',1,'BackgroundColor',figcolor,'FontSize',FS(FSi));

  % add slider only for volume data
  if ~mesh_detected
    % voxelsize and origin
    vx   = sqrt(sum(V(1).mat(1:3,1:3).^2));
    Orig = V(1).mat\[0 0 0 1]';

    H.mm = uicontrol(H.figure,...
      'Units','normalized','position',pos.sslider,...
      ...'Min',(1 - Orig(3))*vx(3) ,'Max',(V(1).dim(3) - Orig(3))*vx(3),...
      'Min', -sum(slices<Orig(3)) * job.gap * vx(3),...
      'Max',  sum(slices>Orig(3)) * job.gap * vx(3),...
      'Style','slider','HorizontalAlignment','center',...
      'callback',@update_slices_array,...
      'ToolTipString','Select slice for display',...
      'SliderStep',[1 job.gap] / (V(1).dim(3)-1),'Visible','off');

    H.mm_txt = uicontrol(H.figure,...
      'Units','normalized','HorizontalAlignment','center',...
      'Style','text','BackgroundColor',figcolor,...
      'Position',[pos.sslider(1) pos.sslider(2)-0.005 0.2 0.02],...
      'String','0 mm','Visible','off','FontSize',FS(FSi));

    update_slices_array;
  end
  
  
  
  % == navigation unit ==
  H.naviuitext = uicontrol(H.figure,...
    'Units','normalized','Style','text','BackgroundColor',figcolor,...
    'Position',[pos.scSelect(1) pos.scSelect(2)+0.035 0.2 0.02],...
    'String','Navigation options','FontSize',FS(FSi));

  H.naviui.select = uicontrol(H.figure,...
    'Units','normalized','position',pos.scSelect,'callback','datacursormode(''on'')',...
    'Style','Pushbutton','enable','on','ToolTipString','Data selection');

  H.naviui.zoomReset = uicontrol(H.figure,...
    'Units','normalized','position',pos.scZoomReset,'callback','zoom out; datacursormode(''on'')',...
    'Style','Pushbutton','enable','on','ToolTipString','Reset zoom'); 

  H.naviui.zoomIn = uicontrol(H.figure,...
    'Units','normalized','position',pos.scZoomIn,'callback',@scZoomIn,...
    'Style','Pushbutton','enable','on','ToolTipString','Zoom in');

  H.naviui.zoomOut = uicontrol(H.figure,...
    'Units','normalized','position',pos.scZoomOut,'callback',@scZoomOut,...
    'Style','Pushbutton','enable','on','ToolTipString','Zoom out');

  H.naviui.pan = uicontrol(H.figure,...
    'Units','normalized','position',pos.scPan,'Enable','off','callback','pan on',...
    'Style','Pushbutton','enable','on','ToolTipString','Hand');



  % == check unit ==
  H.checkuitext = uicontrol(H.figure,...
    'Units','normalized','HorizontalAlignment','center','Style','text',...
    'BackgroundColor',figcolor,...
    'Position',[pos.checkvol(1) pos.checkvol(2)+0.035 0.2 0.02],...
    'String','View selected data','FontSize',FS(FSi));

  % add button to open one image with SPM check_reg
  H.checkui.vol = uicontrol(H.figure,...
    'Units','normalized','position',pos.checkvol,'callback',@checkvol,...
    'string','VOL','ToolTipString','Display original volume in SPM Graphics',...
    'Style','Pushbutton','FontSize',FS(FSi),'Enable','off');  

  % add button to open one image with SPM check_reg
  H.checkui.surf = uicontrol(H.figure,...
    'Units','normalized','position',pos.checksurf,'callback',@checksurf,...
    'string','SURF','ToolTipString','Display processed surfaces in own figure',...
    'Style','Pushbutton','FontSize',FS(FSi),'Enable','off');  

  % add button to open one image with SPM check_reg
  H.checkui.log = uicontrol(H.figure,...
    'Units','normalized','position',pos.checklog,'callback',@checklog,...
    'string','LOG','ToolTipString','Display log-file in SPM Graphics',...
    'Style','Pushbutton','FontSize',FS(FSi),'Enable','off');  

  % add button to open one image with SPM check_reg
  H.checkui.xml = uicontrol(H.figure,...
    'Units','normalized','position',pos.checkxml,'callback',@checkxml,...
    'string','XML','FontSize',FS(FSi),'ToolTipString','Display xml-file in SPM Graphics',...
    'Style','Pushbutton','Enable','off'); 

  % add button to open the pdf in an external viewer 
  H.checkui.pdf = uicontrol(H.figure,...
    'Units','normalized','position',pos.checkpdf,'callback',@checkpdf,...
    'ToolTipString','Display PDF report in external viewer',...
    'Style','Pushbutton','Enable','off');


  % == trashlist unit ==
  H.trashuitext = uicontrol(H.figure,...
    'Units','normalized','Style','text','Position',[pos.newtrash(1) pos.newtrash(2)+0.035 0.2 0.02],...
    'BackgroundColor',figcolor,'String','Trashlist operations','FontSize',FS(FSi));

  % add button for new garbage mask
  H.trashui.new = uicontrol(H.figure,...
    'Units','normalized','position',pos.newtrash,'callback',@newtrash,...
    'string','NEW','ForegroundColor',[ 0 0 0.8],'FontSize',FS(FSi),...
    'ToolTipString','Reset trash list','Style','Pushbutton','Enable','off');  

  % add button to set the active image as garbage
  H.trashui.trash = uicontrol(H.figure,...
    'Units','normalized','position',pos.trash,'callback',@trash,...
    'string','COL+','ForegroundColor',[0.8 0 0],'FontSize',FS(FSi),...
    'ToolTipString','Trash selected subject','Style','Pushbutton','Enable','off');

  % add button to remove the active image from garbage
  H.trashui.detrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.detrash,'callback',@detrash,...
    'string','COL-','ForegroundColor',[0 0.8 0],'FontSize',FS(FSi),...
    'ToolTipString','Trash selected subject','Style','Pushbutton','Enable','off');

  % add button for mask below threshold as garbage
  H.trashui.disptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.disptrash,'callback',@disptrash,...
    'string','VIEW','FontSize',FS(FSi),'ToolTipString','Display trash',...
    'Style','Pushbutton','Enable','off'); 

  H.trashui.autotrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.autotrash,'callback',@autotrash,...
    'string','IQRP','FontSize',FS(FSi),'ForegroundColor',[0.8 0 0.5],...
    'ToolTipString','Automatic IQR thresholding',...
    'Style','Pushbutton','Enable',onoff{(size(QM,2)<4) + 1}); 

  % == second row ==
  H.trashui.undo = uicontrol(H.figure,...
    'Units','normalized','position',pos.undo,'callback',@trashundo,...
    'Style','Pushbutton','Enable','off','ToolTipString','Undo last trashlist operation'); 

  H.trashui.redo = uicontrol(H.figure,...
    'Units','normalized','position',pos.redo,'callback',@trashredo,...
    'Style','Pushbutton','Enable','off','ToolTipString','Redo last trashlist operation'); 

  H.trashui.trashrow = uicontrol(H.figure,...
    'Units','normalized','position',pos.trashrow,'callback',@trashrow,...
    'string','ROW-','ForegroundColor',[0.8 0 0],'FontSize',FS(FSi),...
    'ToolTipString','Tresh row','Style','Pushbutton','Enable','off'); 

  H.trashui.detrashrow = uicontrol(H.figure,...
    'Units','normalized','position',pos.detrashrow,'callback',@detrashrow,...
    'string','ROW+','ForegroundColor',[0 0.8 0],'FontSize',FS(FSi),...
    'ToolTipString','Tresh row','Style','Pushbutton','Enable','off'); 

  H.trashui.ziptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.ziptrash,'callback',@ziptrash,...
    'string','DEL','ForegroundColor',[0.8 0 0],'FontWeight','bold',...
    'FontSize',FS(FSi),'ToolTipString','ZIP selected subject',...
    'Style','Pushbutton','Enable','off');

  % print real data
  show_matrix(YpY, sorted);

  set(H.figure,'Visible','on');  
  
  buttonupdate
 
  % redraw buttons
  pause(0.2) % Wait for the figure construction complete.
  warning off;  %#ok<WNOFF>
  jFig = get(H.figure, 'JavaFrame'); % get JavaFrame. You might see some warnings.
  warning on; %#ok<WNON>
  jWindow = jFig.fHG2Client.getWindow; % before 2011a it could be `jFig.fFigureClient.getWindow`. Sorry I cannot test. 
  jbh = handle(jWindow,'CallbackProperties'); % Prevent memory leak
  set(jbh,'ComponentMovedCallback',{@(~,~)(buttonupdate)});
  
  % show plot
  show_mean_boxplot(mean_cov,'Mean correlation  ',1);

return

function buttonupdate
%-----------------------------------------------------------------------
% This function print icons on buttons
%-----------------------------------------------------------------------
  global H

  % Update figure icons
  % close
  % help
  % check worst
  
  
  % == navi buttons ==
  buttonicon(H.naviui.select    ,'DC'  ,fullfile(matlabroot,'toolbox','matlab','icons','tool_data_cursor.png'));
  buttonicon(H.naviui.zoomReset ,'Zo'  ,fullfile(matlabroot,'toolbox','shared','dastudio','resources','glue','zoom_fit_view_mo.png'));
  buttonicon(H.naviui.zoomIn    ,'Z+'  ,fullfile(matlabroot,'toolbox','matlab','icons','tool_zoom_in.png'));
  buttonicon(H.naviui.zoomOut   ,'Z-'  ,fullfile(matlabroot,'toolbox','matlab','icons','tool_zoom_out.png'));
  buttonicon(H.naviui.pan       ,'H'   ,fullfile(matlabroot,'toolbox','matlab','icons','tool_hand.png'));
  
  % == check buttons ==
  %buttonicon(H.checkui.vol      ,'VOL' ,fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
  %buttonicon(H.checkui.surf     ,'SURF',fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
  %buttonicon(H.checkui.log      ,'PDF' ,fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
  %buttonicon(H.checkui.xml      ,'PDF' ,fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
  buttonicon(H.checkui.pdf      ,'PDF' ,fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
  
  % == trash buttons ==
 %/Applications/MATLAB_R2016a.app/toolbox/shared/comparisons/+comparisons/+internal/+text/deleted.png
 %RestoreOrphanedSignals_16
 %
  buttonicon(H.trashui.undo,'UNDO',fullfile(spm('dir'),'toolbox','cat12','html','icons','restore_24.png'));
  buttonicon(H.trashui.redo,'REDO',fullfile(spm('dir'),'toolbox','cat12','html','icons','Refresh_16.png'));

return

%-----------------------------------------------------------------------
function buttonicon(h,str,Picon) 
%-----------------------------------------------------------------------
% Function to print an image file Picon on the button with handle h. Use
% the string str if the global variable useicons<1 or other errors.
%-----------------------------------------------------------------------
  global FS FSi useicons

  usethisicon = useicons;
  if ~exist(Picon,'file')
    warning('Button Icon "%s" does not exist!',Picon); 
  end
  if ~exist('findjobj','file')
    warning('JAVA Function for  "%s" does not exist!',Picon);
    useicons = 0; 
  end
  if usethisicon
    try
      jButton = findjobj(h);
      jButton.setIcon(javax.swing.ImageIcon(Picon));
      jButton.setHorizontalTextPosition(javax.swing.SwingConstants.LEFT);
      jButton.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
    catch 
      usethisicon = 0;
    end
  end
  if usethisicon<1
    set(h,'string',str,'FontSize',FS(FSi));
  end
  if usethisicon>=1
    set(h,'string','');
  end
return

%-----------------------------------------------------------------------
function autotrash(obj, event_obj) 
%-----------------------------------------------------------------------
% 
%-----------------------------------------------------------------------
  global H pos QM trashlist  fc

if 0
  %%
  global protocol rps2mark
  
  %H.atfigure = figure;
 % d = dialog('Position',pos.popup,'Name','Remove low IQR data');
  %mark2rps   = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;

  rps2mark   = @(rps)  min(10.5,max(0.5,10.5 - rps/10)) + isnan(rps).*rps;
  
  d = dialog('Position',pos.popup,'Name','Remove low IQR data');
  uicontrol('Parent',d,'Style','text','Position',[20 60 160 20],...
     'String','Set treshold (default = 0.72)');
 
  H.atslider = uicontrol(d,...
      'position',[20 20 160 20],...
      'Min', 0.5,'Max', 2.0,'value',0.72,...
      'Style','slider','HorizontalAlignment','center',...
      'ToolTipString','Select slice for display',...
      'callback','global fc QM rps2mark; fc = get(H.atslider,''value''); cat_tst_qa_cleaner( rps2mark(QM(:,3)),struct(''figure'',3,''cf'',fc));',...
      'SliderStep',[0.01 0.05],'Visible','on');

  [Pth,rth,sq,rths,rthsc,sqs] = cat_tst_qa_cleaner_intern( rps2mark(QM(:,3)) , ...
    struct('site',{protocol},'figure',0,'fc',fc));
end


  del  = setdiff(find(QM(:,4)<-0.04*fc)',trashlist);
  if ~isempty(del)
    if isfield(pos,'x'), oldx = pos.x; end
    if isfield(pos,'y'), oldy = pos.y; H.y = rmfield(pos,'y'); end

    for di=1:numel(del)
      pos.x = del(di);
      trash
    end
    
    trashlist = [trashlist, del];

    if exist('oldx','var'), pos.x = oldx; end
    if exist('oldy','var'), pos.y = oldy; end
    
    set([H.trashui.new,H.trashui.disptrash,H.trashui.autotrash,H.trashui.ziptrash],'enable','on');
  else
    fprintf('Nothing to delete.\n');
  end
return

%-----------------------------------------------------------------------
function trash(obj, event_obj) 
%-----------------------------------------------------------------------
% Put a record on the trash list and mark them with a red cross in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global H pos trashlist trashhist mask1d mask2d isscatter
  if isfield(pos,'x') && all( trashlist~=pos.x ) 
    if exist('obj','var')
      trashlist = [trashlist pos.x];
      trashhist = [trashhist pos.x];
    end
    
    showtrash = get(H.showtrash,'Value');
    if ~showtrash && ~isscatter
      H.dcm.removeAllDataCursors;
    end
    
    if exist('obj','var')
      set(H.trashui.trash  ,'Enable','off');
      set(H.trashui.detrash,'Enable','on' );
      set(H.trashui.undo   ,'Enable','on' );
      set(H.trashui.new    ,'Enable','on' );
      set(H.trashui.disptrash,'Enable','on' );
      set(H.trashui.undo   ,'Enable','on' );
      set(H.trashui.redo   ,'Enable','off' );
    end
    if ~isscatter
      set(H.showtrash      ,'Enable','on' );
    end
    
    mask1d(pos.x)    = 0;
    mask2d(pos.x,:)  = 0;
    mask2d(:,pos.x)  = 0;
    
    for scati=1:numel(H.scat)
      scposx  = findobj(H.scat(scati),'type','scatter'); 
      scposxv = cell2mat(get(scposx,'UserData'));
      scposxi = find(scposxv==pos.x,1,'first');

      set(scposx(scposxi),'sizedatasource',get(scposx(scposxi),'marker'),...
        'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
        'ZDataSource','trash','MarkerFaceAlpha',0);
    end
    if ~isscatter
      update_matrix
    end
    
  
  end
return

%-----------------------------------------------------------------------
function trashrow(obj, event_obj) 
%-----------------------------------------------------------------------
% Put a record on the trash list and mark them with a red cross in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global H pos trashlist trashhist mask1d mask2d isscatter
  if isfield(pos,'y') && all( trashlist~=pos.y ) 
    trashlist = [trashlist pos.y];
    trashhist = [trashhist pos.y];
  
    set(H.trashui.trashrow  ,'Enable','off');
    set(H.trashui.detrashrow,'Enable','on' );
    set(H.trashui.undo      ,'Enable','on' );
    set(H.showtrash         ,'Enable','on' );
   
    mask1d(pos.y)    = 0;
    mask2d(pos.y,:)  = 0;
    mask2d(:,pos.y)  = 0;
    
    for scati=1:numel(H.scat)
      scposx  = findobj(H.scat,'type','scatter'); 
      scposxv = cell2mat(get(scposx,'UserData'));
      scposxi = find(scposxv==pos.y,1,'first');

      set(scposx(scposxi),'sizedatasource',get(scposx(scposxi),'marker'),...
        'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
        'ZDataSource','trash','MarkerFaceAlpha',0);
    end
    
    if ~isscatter
      update_matrix;
    end
    
    if isempty(trashlist), set(H.showtrash,'Enable','off'); end
  end
return

%-----------------------------------------------------------------------
function detrash(obj, event_obj)
%-----------------------------------------------------------------------
% Remove a record from trash list and restore the old look like in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global H pos trashlist trashhist mask1d mask2d isscatter
 
  if isfield(pos,'x') && any( trashlist==pos.x ) 
    if exist('obj','var')
      trashlist = setdiff(trashlist,pos.x);
      trashhist = [trashhist -pos.x];

      set(H.trashui.trash  ,'Enable','on' );
      set(H.trashui.detrash,'Enable','off');
    end

    % update matrix
    mask1d(pos.x)    = 1;
    mask2d(pos.x,:)  = sum(mask2d,1)>0;
    mask2d(:,pos.x)  = sum(mask2d,2)>0;

    % update scatter
    for scati=1:numel(H.scat)
      scposx  = findobj(H.scat(scati),'type','scatter'); 
      scposxv = cell2mat(get(scposx,'UserData'));
      scposxi = find(scposxv==pos.x,1,'first');

      set(scposx(scposxi),'marker',get(scposx(scposxi),'sizedatasource'),... 
        'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
    end

    if ~isscatter 
      update_matrix;
      if isempty(trashlist), set(H.showtrash,'Enable','off'); end
    end

  end
return

%-----------------------------------------------------------------------
function detrashrow(obj, event_obj)
%-----------------------------------------------------------------------
% Remove a record from trash list and restore the old look like in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global H pos trashlist trashhist mask1d mask2d isscatter

  if isfield(pos,'y') && any( trashlist==pos.y ) 
    trashlist = setdiff(trashlist,pos.y);
    trashhist = [trashhist -pos.y];
    
    set(H.trashui.trashrow  ,'Enable','on' );
    set(H.trashui.detrashrow,'Enable','off');
    
    % update matrix
    mask1d(pos.y)    = 1;
    mask2d(pos.y,:)  = sum(mask2d,1)>0;
    mask2d(:,pos.y)  = sum(mask2d,2)>0;
    
    % update scatter
    for scati=1:numel(H.scat)
      scposx  = findobj(H.scat(scati),'type','scatter'); 
      scposxv = cell2mat(get(scposx,'UserData'));
      scposxi = find(scposxv==pos.x,1,'first');

      set(scposx(scposxi),'marker',get(scposx(scposxi),'sizedatasource'),... 
        'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
    end
    
    if ~isscatter 
      update_matrix;
      if isempty(trashlist), set(H.showtrash,'Enable','off'); end
    end
    
  end
return

%-----------------------------------------------------------------------
function newtrash(obj, event_obj)
%-----------------------------------------------------------------------
% Create an empty trash list. 
%-----------------------------------------------------------------------
  global trashlist trashhist trashhistf H pos isscatter mask2d mask1d
  trashlist  = [];
  trashhist  = []; 
  trashhistf = []; 

  mask2d = true(size(mask2d));
  mask1d = true(size(mask1d));
   
  % find scatter objects 
  sc = findobj('ZDataSource','trash');
  for sci=1:numel(sc)
   set(sc(sci),'marker',get(sc(sci),'sizedatasource'),... 
        'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
  end
  
  if ~isscatter 
    update_matrix;
  end
  
  %trashhist = [trashhist -trashlist];
  %H.trashui.undo,H.trashui.redo,
  %H.trashui.trash,H.trashui.trashrow,,H.showtrasht
  unit = struct2cell(H.checkui); 
  set([unit{:}],'Enable','off');
  set([H.trashui.trash,H.trashui.trashrow,H.trashui.detrash,H.trashui.detrashrow,H.trashui.undo,H.trashui.redo],'Enable','off');
  if ~isempty(H.dcm.getCursorInfo)
    if isfield(pos,'x'), set(H.trashui.trash   ,'Enable','on'); end
    if isfield(pos,'y'), set(H.trashui.trashrow,'Enable','on'); end
  end
return

%-----------------------------------------------------------------------
function disptrash(obj, event_obj)
%-----------------------------------------------------------------------
% List all records of the trash list in the command window.
%-----------------------------------------------------------------------
  global trashlist data_files 

  fprintf('\nTrashlist:\n')
  for fi=1:numel(trashlist)
    fprintf('  %s\n',data_files{trashlist(fi)});
  end
return

%-----------------------------------------------------------------------
function trashundo(obj, event_obj)
%-----------------------------------------------------------------------
% List all records of the trash list in the command window.
%-----------------------------------------------------------------------
  global H pos trashlist trashhist trashhistf 
  
  if isfield(pos,'x'), oldx = pos.x; end
  if isfield(pos,'y'), oldy = pos.y; end
  if isfield(pos,'tar_mouseget'), oldtar_mouseget = pos.tar_mouseget; end
  pos.x = abs(trashhist(end)); 
  
  if ~isempty(trashlist) && trashhist(end)>0
    % last element was added to the trashlist 
    detrash
    trashlist(end) = [];
  else
    trash
    trashlist(end+1) = trashhist(end);
  end
  trashhistf(end+1) = trashhist(end); 
  trashhist(end)    = []; 
  
  set(H.trashui.redo,'enable','on');
  if isempty(trashhist), set(H.trashui.undo,'enable','off'); end
  
  if exist('oldx','var')
    pos.x = oldx; 
  else
    pos = rmfield(pos,'x'); 
  end
  if exist('oldy','var')
    pos.y = oldy; 
  elseif isfield(pos,'y') 
    pos = rmfield(pos,'y'); 
  end
  if exist('oldtar_mouseget','var')
    pos.tar_mouseget = oldtar_mouseget; 
  elseif isfield(pos,'tar_mouseget')
    pos = rmfield(pos,'tar_mouseget'); 
  end
return
%-----------------------------------------------------------------------
function trashredo(obj, event_obj)
%-----------------------------------------------------------------------
% List all records of the trash list in the command window.
%-----------------------------------------------------------------------
  global H pos trashlist trashhist trashhistf 
  
  if isfield(pos,'x'), oldx = pos.x; end
  if isfield(pos,'y'), oldy = pos.y; end
  if isfield(pos,'tar_mouseget'), oldtar_mouseget = pos.tar_mouseget; end
  pos.x = trashhistf(end); 
  
  if isempty(trashlist) && trashhistf(end)<0
    detrash
    trashlist(end) = [];
  else
    trash
    trashlist(end+1)  = trashhistf(end);
  end
  trashhist(end+1) = trashhistf(end); 
  trashhistf(end)  = [];
  
  set(H.trashui.undo,'enable','on');
  if isempty(trashhistf), set(H.trashui.redo,'enable','off'); end
  
  if exist('oldx','var')
    pos.x = oldx; 
  else
    pos = rmfield(pos,'x'); 
  end
  if exist('oldy','var')
    pos.y = oldy; 
  elseif isfield(pos,'y') 
    pos = rmfield(pos,'y'); 
  end
  if exist('oldtar_mouseget','var')
    pos.tar_mouseget = oldtar_mouseget; 
  elseif isfield(pos,'tar_mouseget')
    pos = rmfield(pos,'tar_mouseget'); 
  end
return

%-----------------------------------------------------------------------
function ziptrash(obj, event_obj)
%-----------------------------------------------------------------------
% Remove records and related files from the file system by zipping or
% storing in a separate directory (NOT READY). 
%-----------------------------------------------------------------------
  global X QM trashlist org_files 
  
  trashtime = datestr(clock,'yyyymmdd-HHMMSS');
  fprintf('Zip preprocessed data (trashtime = %s) of:\n',trashtime);
  for fi=numel(trashlist):-1:1
    fprintf('  Zip %s\n',org_files{trashlist(fi)});
    
    %% find preprocessed images 
    [pp,ff,ee] = spm_fileparts(org_files{trashlist(fi)});
    sim_files  = cat_vol_findfiles(pp,['*' ff '*.nii'],struct('maxdepth',1)); 
    sim_files  = [fullfile(pp,[ff ee]);setdiff(sim_files,fullfile(pp,[ff ee]))]; 
    for si=1:numel(sim_files); 
      [pps,ffs]    = spm_fileparts(sim_files{si});
      pp_files{si} = cat_vol_findfiles(pps,['*' ffs '*']);
      pp_files{si} = setdiff(pp_files{si},sim_files);
      for fsi = numel(pp_files{si}):-1:1
        ppfs = spm_str_manip(pp_files{si}{fsi},'hht');
        switch ppfs
          case 'err', pp_files{si}(fsi) = [];
        end
      end
      pp_files{si} = setdiff(pp_files{si},sim_files{si}); 
      if si==1
        pp_filescor = pp_files{si};
      else
        pp_filescor  = setdiff(pp_filescor,pp_files{si}); 
      end
    end
    
    
    %% zip the list and remove the files
    %  we need to go into the directory and use the short filenames 
    %  to opbtain save the relative path in the zip file!
    %odir = dir; cd(pp); 
    %pp_filescor1 = cat_io_strrep(pp_filescor,[pp filesep],''); 
    ffzip = fullfile(pp,sprintf('%s_cor%2.2f_IQR%2.2f_trashed%s',...
      ff,X(trashlist(fi),1),QM(trashlist(fi),3),trashtime));
    zip(ffzip,pp_filescor1,pp); 
    
    %for fii=1:numel(pp_filescor1), delete(pp_filescor1{fii}); end
    cd(odir); 
    
    
  end
  % update all variables :-/ 
  % or recreate job?

return

%-----------------------------------------------------------------------
function scZoomIn(obj, event_obj)
%-----------------------------------------------------------------------
  global H isscatter
  if isscatter
    hz = zoom(H.scat(H.scata));
  else
    hz = zoom(H.corr);
  end
  set(hz,'enable','on','direction','in');
return

%-----------------------------------------------------------------------
function scZoomOut(obj, event_obj)
%-----------------------------------------------------------------------
  global H isscatter
  if isscatter
    hz = zoom(H.scat(H.scata));
  else
    hz = zoom(H.corr);
  end
  set(hz,'enable','on','direction','out');
return 

%-----------------------------------------------------------------------
function closeWindows(obj, event_obj)
%-----------------------------------------------------------------------
% Close all windows. 
% Remove variables (NOT DONE).
%-----------------------------------------------------------------------
  global trashlist pos H
  %%
  posx = get(get(event_obj.Source,'Parent'),'Position');
  pos.popup(1:2) = [posx(1) + posx(3)*0.8, posx(1) + posx(4)*0.9];  
  if ~isempty(trashlist)
    d = dialog('Position',pos.popup,'Name','Close Check Sample');
    uicontrol('Parent',d,'Style','text','Position',[20 60 160 20],...
       'String','Trashlist not empty!');
    uicontrol('Parent',d,'TooltipString','Sopt closing',...
       'Position',[25 20 70 25],'String','Cancel','Callback','delete(gcf)');    
    uicontrol('Parent',d,'TooltipString','Close windows without last changes',...
       'Position',[100 20 70 25],'String','Close','ForegroundColor',[0.8 0 0],...
       'Callback','for i=2:26, try close(i); end; end; delete(gcf)');   
  else
    for i=2:26, try close(i); end; end; 

    spm_figure('Clear',H.graphics);
    
    clearvars -GLOBAL ...
       filename ... structure from spm_str_manip with grouped filenames
       H        ... main structure with object/button handles
       pos      ... structure with position values for GUI objects/buttons
       YpY mean_cov         ... (sorted) covariance and mean covariance matrix
       mask1d mask2d                  ... masks for files on the trash list 
       trashlist                      ... index list of subjects to remove
       ind_sorted   ... index lists 
       QM                             ... Quality measures from xml_files
       data_files org_files surf_files xml_files log_files pdf_files ...  different lists of filenames 
       data_array data_array_diff     ... slices of all subjects 
       img                            ... slice image(s) for GUI display
       FS FSi                         ... SPM fontsize 
       mesh_detected isscatter isxml sorted show_name useicons ... binary variables for GUI control 
       mn_data mx_data                ... minimum/maximum value of the mesh
       V Vo                           ... volume header structures
       Vchanged names_changed         ... modified volume header for 4D-structures
       sample                         ... sample groups of each scan
       dataprefix                     ... prefix of the input files
       inorm                          ... normalize slice intensity even in normalized data
       X X2 MD MD2; 
  end
return

%-----------------------------------------------------------------------
function id = mygetCursorInfo
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
  global H pos isscatter

  curs = H.dcm.getCursorInfo;
  
  if isscatter
    sc = unique([curs(:).Target]);
    id = get(sc,'UserData')'; 
    if iscell(id), id = cell2mat(id); end
  else
    id = unique([pos.x, pos.y]);
    %{
    pos  = reshape([curs(:).Position],numel(curs),2);   
    posx = unique(pos(:,1));
    posy = unique(pos(:,2));
  
    if sorted
      
    else
      
    end
    %}
  end
  
return

%-----------------------------------------------------------------------
function checkpdf(obj, event_obj)
%-----------------------------------------------------------------------
% Open PDF report of selected subjects. 
% This is only possible for using an extern viewer. 
% Hence, it would be useful to save a JPG or HTML file in cat_main.
%-----------------------------------------------------------------------
  global pdf_files 
  
  id = mygetCursorInfo;
  for i=1:numel(id), open(pdf_files{id(i)}); end
return

%-----------------------------------------------------------------------
function checksurf(obj, event_obj)
%-----------------------------------------------------------------------
% Open surface files of selected subjects. 
% This is very slow and some feedback would be useful. 
%-----------------------------------------------------------------------
  global pos surf_files isscatter
  
  if isscatter
    cat_surf_display(struct('data',surf_files{pos.x},'multisurf',1));
  else
    cat_surf_display(struct('data',char(surf_files(unique([pos.x,pos.y]))),'multisurf',1));
  end
  
  % give some feedback
return

%-----------------------------------------------------------------------
function checkvol(obj, event_obj)
%-----------------------------------------------------------------------
% Load the original image of selected files in SPM graphics window.
% Some further information or legend would be helpful.
%-----------------------------------------------------------------------
  global H pos org_files isscatter FS FSi
  
  spm_figure('Clear',H.graphics);
  spm_figure('Focus',H.graphics);
  spm_orthviews('Reset')
  gax = gca; set(gax,'Position',[0 0 1 1]); axis off;
  
  xeqy = pos.x == pos.y; 
  
  multi = 1; 
  if isscatter 
    if multi
      id = mygetCursorInfo';
      
      spm_check_registration(char(unique(org_files(id))));
      
      spm_orthviews('MaxBB')
    else
      ppos = [0.02 0.01 0.96 0.98];
      tpos = [0.50 0.98 0.96 0.02];
    end
  else
    if xeqy
      ppos = [0.02 0.01 0.96 0.98];
      tpos = [0.50 0.98 0.96 0.02];
    else
      ppos = [0.02 0.545 0.96 0.48 ; 0.02 0.010 0.96 0.48];
      tpos = [0.50 0.980 0.96 0.02 ; 0.50 0.475 0.96 0.02];
    end
  end
  
  if ~isscatter || ~multi || xeqy
    text(gax,tpos(1,1),tpos(1,2),spm_str_manip(org_files{pos.x},'k100'),...
      'FontSize',FS(FSi+1),'Color',[0 0 0.8],'LineStyle','none',...
      'HorizontalAlignment','center');

    hi1 = spm_orthviews('Image',spm_vol(org_files{pos.x}),ppos(1,:)); 
    spm_orthviews('AddContext',hi1);
  end
  
  if ~isscatter && ~xeqy
    text(gax,tpos(1,1),tpos(2,2),spm_str_manip(org_files{pos.y},'k100'),...
      'FontSize',FS(FSi+1),'Color',[0 0 0.8],'LineStyle','none',...
      'HorizontalAlignment','center');
    hi2 = spm_orthviews('Image',spm_vol(org_files{pos.y}),ppos(2,:)); 
    spm_orthviews('AddContext',hi2);
    spm_orthviews('MaxBB')
  end
  
return

%-----------------------------------------------------------------------
function checkxml(obj, event_obj)
%-----------------------------------------------------------------------
% Load XML report in SPM graphics window (see also checklog).
% This is just the first fast version of this function. 
% Finally, I want to use the xml structure from the file to print some
% specific informations similar to the CAT report in cat_main. 
%-----------------------------------------------------------------------
  global H pos xml_files isscatter

  % visdiff(xml_files{pos.x}, xml_files{pos.y},'text')  
  
  spm_figure('Clear',H.graphics); 
  spm_figure('Focus',H.graphics);
  axis off;
  
  if isscatter || (pos.x == pos.y)
    textbox = [0 0 1 1];
    files   = xml_files(pos.x); 
  else
    textbox = [0 0.5 1 0.5; 0 0 1 0.5];
    files   = xml_files([pos.y,pos.x]); 
  end
  
  % avoid some long useless text passages
  badtacks = {'software>','catlog>','atlas>','LAB>'};
  badmode = 0; bdid = badtacks;
  for fi=1:numel(files);
    fid = fopen(files{fi});
    ph  = uipanel(H.graphics,'Units','normalized','position',textbox(fi,:), ...
      'BorderWidth',0,'title',[spm_str_manip(files{fi},'k100') ' (extract)'],'ForegroundColor',[0 0 0.8]);
    lbh = uicontrol(ph,'style','listbox','Units','normalized',...
      'fontname','Fixedwidth','position',[ 0 0 1 1 ],'FontSize',9);
    indic = 1;
    indit = 1; 
    while 1
     tline = fgetl(fid);
     if ~ischar(tline), 
       break
     end
     for bi = 1:numel(badtacks), bdid{bi} = strfind(tline,badtacks{bi}); end
     if any(~cellfun('isempty',bdid))
       badmode = ~badmode; 
     end
     if ~badmode
       strings{indit}=tline; 
       indit = indit + 1;
     end
     indic = indic + 1;
    end
    fclose(fid);
    
    set(lbh,'string',strings);
    set(lbh,'Value',1);
    set(lbh,'Selected','on');
  end
return

%-----------------------------------------------------------------------
function show_sample(obj,event_obj)
%-----------------------------------------------------------------------
  global H sample pmask1d smask1d smask2d isscatter
  
  % set entry of the GUI element
  if obj>0
    set(H.samp,'Value',obj+2);
    smask1d(:) = sample==obj; 
  else
    set(H.samp,'Value',1); 
    smask1d(:) = true(size(smask1d)); 
  end
  
  smask2d = (single(smask1d) * single(smask1d'))>0; 
  
  groups        = unique(sample);
  symbols       = repmat('.',1:numel(groups));  % default symbol
  symbols(1:11) = 'o+^v<>ph*sd'; 
  
  % update scatter
  if isscatter
    if obj>0
      for scati=1:numel(H.scat)
        scpos  = [
          findobj(H.scat(scati),'type','scatter','marker',symbols(obj)); 
          findobj(H.scat(scati),'type','scatter','sizedatasource',symbols(obj))];
        scpos  = setdiff(scpos,H.sclegend);
        set(scpos,'Visible','on');
        scposn = setdiff(findobj(H.scat(scati),'type','scatter'),scpos);
        set(scposn,'Visible','off');
      end
    else
      for scati=1:numel(H.scat)
        set( setdiff( findobj(H.scat(scati),'type','scatter') , H.sclegend ),'Visible','on');
      end
    end
    for linei=1:numel(H.corrline)
      try
        indxy = get(H.corrline(linei),'UserData');
        if any(find(smask1d) == indxy(1)) && any(find(pmask1d) == indxy(1)) && ...
           any(find(smask1d) == indxy(2)) && any(find(pmask1d) == indxy(2))
          set(H.corrline(linei),'Visible','on');
        else
          set(H.corrline(linei),'Visible','off');
        end
      end
    end
  else
    update_matrix
  end
  
  H.dcm.removeAllDataCursors
return

%-----------------------------------------------------------------------
function show_protocol(obj,event_obj)
%-----------------------------------------------------------------------
  global H protocol pmask1d smask1d pmask2d isscatter protocols
  
  % set GUI element entry
  if obj>0
    set(H.prot,'Value',obj+2);
    pmask1d(:) = protocol==obj; 
  else
    set(H.prot,'Value',1);
    pmask1d(:) = true(size(pmask1d)); 
  end
  
  pmask2d = (single(pmask1d) * single(pmask1d'))>0; 
  
  % update scatter
  if isscatter
    if obj>0
      for scati=1:numel(H.scat)
        scpos  = findobj(H.scat(scati),'type','scatter','ZDataSource',num2str(obj,'%d')); 
        scpos  = setdiff(scpos,H.sclegend);
        set(scpos,'Visible','on');
        scposn = setdiff(findobj(H.scat(scati),'type','scatter'),scpos);
        set(scposn,'Visible','off');
      end
    else
      for scati=1:numel(H.scat)
        set( setdiff( findobj(H.scat(scati),'type','scatter') , H.sclegend ),'Visible','on');
      end
    end
    for linei=1:numel(H.corrline)
      try % deletet object > cleanup H?
        indxy = get(H.corrline(linei),'UserData');
        if any(find(pmask1d) == indxy(1)) && any(find(smask1d) == indxy(1)) && ...
           any(find(pmask1d) == indxy(2)) && any(find(smask1d) == indxy(2))
          set(H.corrline(linei),'Visible','on');
        else
          set(H.corrline(linei),'Visible','off');
        end
      end
    end
  else
    update_matrix
  end
  
  H.dcm.removeAllDataCursors
return

%-----------------------------------------------------------------------
function checklog(obj, event_obj)
%-----------------------------------------------------------------------
% Load the log-file from cat_main of the selected subjects into the SPM
% graphics window.
%-----------------------------------------------------------------------
  global H pos log_files isscatter
  
  spm_figure('Clear',H.graphics); 
  spm_figure('Focus',H.graphics);
  axis off;
  
  if isscatter || (pos.x == pos.y)
    textbox = [0 0 1 1];
    files   = log_files(pos.x); 
  else
    textbox = [0 0.5 1 0.5; 0 0 1 0.5];
    files   = log_files([pos.x,pos.y]); 
  end
  
  for fi=1:numel(files); 
    fid = fopen(files{fi});
    ph  = uipanel(H.graphics,'Units','normalized','position',textbox(fi,:), ...
      'BorderWidth',0,'title',spm_str_manip(files{fi},'k100'),'ForegroundColor',[0 0 0.8]);
    lbh = uicontrol(ph,'style','listbox','Units','normalized',...
      'fontname','Fixedwidth','position',[ 0 0 1 1 ],'FontSize',9);
    indic = 1;
    while 1
     tline = fgetl(fid);
     if ~ischar(tline), 
       break
     end
     strings{indic}=tline; 
     indic = indic + 1;
    end
    fclose(fid);
    set(lbh,'string',strings);
    set(lbh,'Value',1);
    set(lbh,'Selected','on');
  end
return

%-----------------------------------------------------------------------
function checkbox_names(obj, event_obj)
%-----------------------------------------------------------------------
  global H show_name
  show_name = get(H.chbox,'Value');
  show_mean_boxplot;
return
      
%-----------------------------------------------------------------------
function checkbox_showtrash(obj, event_obj)
%-----------------------------------------------------------------------
  global oldx oldy 
  global H pos isscatter %trashlist YpY

  showtrash = get(H.showtrash,'Value');
  
  if ~showtrash
    if isfield(pos,'x'), oldx = pos.x; else oldx = 0; end
    if isfield(pos,'y'), oldy = pos.y; else oldy = 0; end
  else
    if oldx>0, pos.x = oldx; end
    if oldy>0, pos.y = oldy; end
  end
  H.dcm.removeAllDataCursors;
  set(H.slice,'visible','off');  
 
  if ~showtrash 
    if ~isscatter
      if isfield(pos,'x'), pos = rmfield(pos,'x'); end
      if isfield(pos,'y'), pos = rmfield(pos,'y'); end
      H.dcm.removeAllDataCursors;
    end
  end
  
  update_matrix;
  
  if showtrash
    %createDatatip(H.dcm,get(H.corr,'children'),[oldx oldy]);
  else
    %createDatatip(H.dcm,get(H.corr,'children'),[oldx oldy]);
  end
return  

%-----------------------------------------------------------------------
function checkbox_cbarfix(obj, event_obj)
%-----------------------------------------------------------------------
  global isscatter
  
  if isscatter
    % do something
  else
    update_matrix;
  end
  
return  

%-----------------------------------------------------------------------
function show_mahalanobis(X,MD,scata)
%-----------------------------------------------------------------------
  global H FS FSi pos isscatter sample trashlist YpY mesh_detected  protocol
  
  if ~exist('scata','var')
    if isfield(H,'scata')
      scata = H.scata; 
    else
      scata = 1; H.scata = scata;
    end
  else
    H.scata = scata;
  end
  axis(H.scat(H.scata));
  %set(H.scat(H.scata),'

  % clear larger area and set background color to update labels and title
  if ~isscatter
    set([H.corr,get(H.corr,'children')],'visible','off','HitTest','off','Interruptible','off')
    if isfield(H,'cbarfix')   && ishandle(H.cbarfix),   set(H.cbarfix  ,'enable' ,'off'); end
    if isfield(H,'showtrash') && ishandle(H.showtrash), set(H.showtrash,'enable' ,'off'); end
    
    % remove sliders and text
    if isfield(pos,'tar_mouse')
      delete(pos.tar_mouse); 
      pos = rmfield(pos,'tar_mouse');
      set(H.alphabox,'Visible','off');

      if isfield(pos,'x'), pos = rmfield(pos,'x'); end
      if isfield(pos,'y'), pos = rmfield(pos,'y'); end

      if isfield(H,'slice')    && ishandle(H.slice),    set(H.slice        ,'Visible','off'); cla(H.slice); end
      if isfield(H,'sslider')  && ishandle(H.sslider),  set(H.sslider      ,'Visible','off'); end
      if isfield(H,'alphabox') && ishandle(H.alphabox), set(H.alphabox     ,'Visible','off'); end
      if isfield(H,'mm')       && ishandle(H.mm),       set([H.mm,H.mm_txt],'Visible','off'); end
      
      if ~mesh_detected
        set([H.mm,H.mm_txt],'Visible','off');
      end
      set([H.trashui.trash,H.trashui.trashrow,H.trashui.detrash,H.trashui.detrashrow],'Enable','off');
      unit = struct2cell(H.checkui); set([unit{:}],'Enable','off');  
    end
  end
  for i=setdiff(1:numel(H.scat),H.scata)
    set([H.scat(i);get(H.scat(i),'Children')],'visible','off','HitTest','off','Interruptible','off');
  end
  set([H.scat(H.scata);get(H.scat(H.scata),'children')],'visible','on','HitTest','on','Interruptible','on');
  try set([H.sclegend],'visible','off'); end
  set(findobj('type','Legend'),'visible','on');

  if isempty(get(H.scat(H.scata),'children'))
    % get very similar scans
    YpY_tmp = YpY - tril(YpY);
    [indx, indy] = find(YpY_tmp>0.925);

    groups        = unique(sample);
    symbols       = repmat('.',1:numel(groups));  % default symbol
    symbols(1:11) = 'o+^v<>ph*sd';                % need x for unset

    %% Create legend by creation of hidden dummy objects 
    %  display first object for the legend
    hold(H.scat(H.scata),'on'); 
    for gi=1:numel(groups)
      txt{gi} = sprintf('sample %d \n',gi);
      Xt = X(sample==groups(gi),:); 
      H.sclegend(gi) = scatter(H.scat(H.scata),Xt(1,1),Xt(1,2),30,[0 0 0],symbols(gi),'Linewidth',2);
      if ~isempty( strfind('osd^v<>ph', symbols(gi) ) )
        set(H.sclegend(gi),'MarkerFaceColor','flat','markerFaceAlpha',1/3);
      end
    end
    txt{end+1} = 'trashlist';
    H.sclegend(gi+1) = scatter(H.scat(H.scata),Xt(1,1),Xt(1,2),30,[1 0 0],'x','Linewidth',2,'Visible','off');
    if numel(indx)/size(YpY,1)<0.5 && numel(indx)>0
      txt{end+1}   = 'highly corr. scans'; 
      H.sclegend(gi+2) = plot(H.scat(H.scata),[X(indx(1),1);X(indy(1),1)],[X(indx(1),2);X(indy(1),2)],'Color',[0 0 0],'Linewidth',2);
    end
    % create legend
    hl = legend(H.scat(H.scata),txt,'location','southwest');
    set(get(hl,'title'),'string','Legend');
    hold(H.scat(H.scata),'off'); 
    set([H.sclegend],'visible','off');

    
    %%
    %X(isnan(X)) = 0;  
    %S  = cov(X); 
   % mu = mean(X);
   % MD = (X-repmat(mu,[size(X,1),1]))*inv(S)*(X-repmat(mu,[numel(X),1]))';
   % MD = diag(MD);
   
    % because we use a splitted colormap we have to set the color values explicitely
    MDs  = min(63,max(0,64*MD/10)); %(round(max(MD)/6)*6);
    C    = zeros(numel(MD),3);
    cmap = [jet(64); gray(64)];
    for i=1:numel(MD)
      C(i,:) = cmap(round( MDs(i) )+1,:);
    end

    hold(H.scat(H.scata),'on')

    % plot lines between similar objects
    if numel(indx)/size(YpY,1)<0.5 && numel(indx)
      for i=1:numel(indx)
        H.corrline(i) = plot(H.scat(H.scata),[X(indx(i),1);X(indy(i),1)],[X(indx(i),2);X(indy(i),2)],...
          '-','Color',repmat(0.9 - 0.9*((YpY_tmp(indx(i),indy(i))-0.925)/0.0725),1,3),...
          'LineWidth',2,'HitTest','off','Interruptible','off');
        set(H.corrline(i),'UserData',[indx(i) indy(i)]);
      end
    else
      H.corrline = struct(); 
    end

    % plot data entries
    I   = 1:size(X,1);
    for gi=1:numel(groups)
      It = I(sample==groups(gi)); 
      Xt = X(sample==groups(gi),:); 
      Ct = C(sample==groups(gi),:);
      Pt = protocol(sample==groups(gi),:);
      H.sc{scata} = cell(size(Xt,1),1);
      for sci=1:size(Xt,1)
        H.sc{scata}{gi}{sci} = scatter( H.scat(H.scata), ...
          Xt(sci,1), ...
          Xt(sci,2), ...
          30,...
          Ct(sci,:),...
          symbols(gi), ...
          'ZDataSource',num2str(Pt(sci),'%d'),...
          'UserData',It(sci),...
          'Linewidth',2);
        if ~isempty( strfind('osd^v<>ph', symbols(gi) ) )
          set(H.sc{scata}{gi}{sci},'MarkerFaceColor','flat','markerFaceAlpha',1/3);
        end
        if any(It(sci)==trashlist)
          set(H.sc{scata}{gi}{sci},'sizedatasource',get(H.sc{scata}{gi}{sci},'marker'),...
          'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
          'ZDataSource','trash','MarkerFaceAlpha',0);
        end
      end
    end

    hold(H.scat(H.scata),'off')
  end

  xlabel(H.scat(H.scata),'<----- Worst ---      Mean correlation       --- Best ------>  ','FontSize',FS(FSi),'FontWeight','Bold');
  ylabel(H.scat(H.scata),'<----- Worst ---      Weighted overall image quality rating     --- Best ------>  ','FontSize',FS(FSi),'FontWeight','Bold');
  title(H.scat(H.scata),'<--- Best -- Mahalanobis distance (Color) -- Worst ---->  ','FontSize',FS(FSi+2),'FontWeight','Bold');


  % update colorbar 
  cticks = 7;
  mn     = 0; 
  mx     = 10; %round(max(MD)/6)*6;
  ticks  = linspace(mn,mx,cticks);
  set(H.cbar,'XTick',1:63/(numel(ticks)-1):64,'XTickLabel',cellstr(num2str(ticks','%0.2f'))); %round(100*linspace(min(YpYt(:)),max(YpYt(:)),5))/100);
  %caxis(H.scat(H.scata),[mn mx])
  
  
  
  isscatter = 1;
  zoom reset
return

%-----------------------------------------------------------------------
function update_matrix(data, order)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
  global H YpY sorted mask1d mask2d ind_sorted mean_cov isscatter trashlist ...
    smask1d smask2d pmask1d pmask2d

  showtrash = get(H.showtrash,'Value'); 
 
  if ~exist('order','var')
    order = sorted; 
  else
    sorted = order; 
  end
  
  if ~exist('data' ,'var')
    if order
      data = YpY(ind_sorted,ind_sorted);
    else
      data = YpY;
    end
  end

  if order 
    mask1dt  = mask1d(ind_sorted); 
    mask2dt  = mask2d(ind_sorted,ind_sorted);
    
    smask1dt = smask1d(ind_sorted);
    smask2dt = smask2d(ind_sorted,ind_sorted);
    pmask1dt = pmask1d(ind_sorted); 
    pmask2dt = pmask2d(ind_sorted,ind_sorted);    
  else
    mask1dt = mask1d;
    mask2dt = mask2d;
    
    smask1dt = smask1d;
    smask2dt = smask2d;
    pmask1dt = pmask1d; 
    pmask2dt = pmask2d;    
  end
  
  
  if ~showtrash
    data = reshape(data(mask2dt & smask2dt & pmask2dt),...
      sum(mask1dt .* smask1dt .* pmask1dt),sum(mask1dt .* smask1dt .* pmask1dt)); 
    mask2dt = reshape(mask2dt(mask2dt & smask2dt & pmask2dt),...
      sum(mask1dt .* smask1dt .* pmask1dt),sum(mask1dt .* smask1dt .* pmask1dt)); 
  else
    data = reshape(data(smask2dt & pmask2dt),...
      sum(smask1dt .* pmask1dt),sum(smask1dt .* pmask1dt)); 
    mask2dt = reshape(mask2dt(smask2dt & pmask2dt),...
      sum(smask1dt .* pmask1dt),sum(smask1dt .* pmask1dt)); 
  end
  
   
  %% scale data 
  cticks  = 7;
  ltick   = cticks - 1; 
  htick   = ltick/2;
  try cbarfix = get(H.cbarfix,'Value'); catch, cbarfix = 1; end
  if cbarfix
    mx = 1.0;
    mn = 0.7; 
  elseif cbarfix==2
    mx = 1.0;
    md = median(mean_cov); 
    mn = mx - ltick * ceil((mx - md)/htick * 100)/100; 
  else
    md = round(median(mean_cov)*100)/100; 
    mx = md + htick * round(2*std(mean_cov) * 100)/100;
    mn = md - htick * round(2*std(mean_cov) * 100)/100; 
  end
  ticks = linspace(mn,mx,cticks);
 
  data_scaled = min(1,max(0,(data - mn)/(mx - mn)));

  % create image if not exist
  if isempty(get(H.corr,'children'))  
    image(H.corr,data_scaled);
  else
    if showtrash
      mylim = 0.5 + [0 size(YpY,1)] - [0 numel(mean_cov) - sum(smask1dt .* pmask1dt)]; 
    else
      mylim = (0.5 + [0 size(YpY,1)]) - [0 numel(mean_cov) - sum(mask1dt .* smask1dt .* pmask1dt)];
    end
    xlim(H.corr,mylim);
    ylim(H.corr,mylim);
  end
  
  % show only lower left triangle
  if showtrash
    set(get(H.corr,'children'),'Cdata',64 * (data_scaled + (1-mask2dt)) .* tril(data>0) );
  else
    set(get(H.corr,'children'),'Cdata',64 * (data_scaled) .* tril(data>0) );
  end
  
  % update colorbar 
  set(H.cbar,'XTick',1:63/(numel(ticks)-1):64,'XTickLabel',cellstr(num2str(ticks','%0.2f'))); %round(100*linspace(min(YpYt(:)),max(YpYt(:)),5))/100);

  % update axis limits
  if 0 && ~isscatter
    if showtrash
      mylim = 0.5 + [0 size(YpY,1)]; 
    else
      mylim = (0.5 + [0 size(YpY,1)]) - [0 numel(trashlist)];
    end
    xlim(H.corr,mylim);
    ylim(H.corr,mylim);
  end
  
  if isempty(H.dcm.getCursorInfo)
    unit = struct2cell(H.checkui); set([unit{cellfun(@ishandle,unit)}],'Enable','off');  
    set([H.trashui.trash;H.trashui.detrash;H.trashui.trashrow;H.trashui.detrashrow],'Enable','off');
  end
return

%-----------------------------------------------------------------------
function show_matrix(data, order)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
  global H FS FSi pos sorted isscatter %mesh_detected 
 
  set([H.corr,get(H.corr,'children')],'visible','on','HitTest','on','Interruptible','on');
  set(findobj('type','Legend'),'visible','off','HitTest','off','Interruptible','off');
  try
    for i=1:numel(H.scat), set([H.scat(i);get(H.scat(i),'Children')],...
        'visible','off','HitTest','off','Interruptible','off'); end; 
  end
  axis(H.corr);
  if isfield(H,'cbarfix')   && ishandle(H.cbarfix),   set(H.cbarfix  ,'enable' ,'on'); end
  if isfield(H,'showtrash') && ishandle(H.showtrash), set(H.showtrash,'enable' ,'on'); end
  try set([H.sclegend],'visible','off'); end
 
  % == set GUI fields ==
  
  % update buttons and remove cursor, slices/surfaces? and text
  if isfield(pos,'tar_mouse')
    delete(pos.tar_mouse); 
    pos = rmfield(pos,'tar_mouse'); 
  
    if isfield(pos,'x'), pos = rmfield(pos,'x'); end
    if isfield(pos,'y'), pos = rmfield(pos,'y'); end

    if isfield(H,'slice')    && ishandle(H.slice),    set(H.slice        ,'Visible','off'); cla(H.slice); end
    if isfield(H,'sslider')  && ishandle(H.sslider),  set(H.sslider      ,'Visible','off'); end
    if isfield(H,'alphabox') && ishandle(H.alphabox), set(H.alphabox     ,'Visible','off'); end
    if isfield(H,'mm')       && ishandle(H.mm),       set([H.mm,H.mm_txt],'Visible','off'); end

    set([H.trashui.trash,H.trashui.trashrow,H.trashui.detrash,H.trashui.detrashrow],'Enable','off');
    unit = struct2cell(H.checkui); set([unit{cellfun(@ishandle,unit)}],'Enable','off');  

    set(H.alphabox,'Visible','off');
  end
  
  isscatter = 0;

  % get sorting order
  sorted = order;
  
  % update image
  update_matrix(data,order)

  % update title and label elements
  if sorted
    xlabel(H.corr,'<----- Best ---      File Order      --- Worst ------>  ','FontSize',FS(FSi),'FontWeight','Bold');
    ylabel(H.corr,'<----- Worst ---      File Order      --- Best ------>  ','FontSize',FS(FSi),'FontWeight','Bold');
    title(H.corr,'Sorted Sample Correlation Matrix  ','FontSize',FS(FSi+2),'FontWeight','Bold');
  else
    xlabel(H.corr,'<----- First ---      File Order      --- Last ------>  ','FontSize',FS(FSi),'FontWeight','Bold');
    ylabel(H.corr,'<----- Last ---      File Order      --- First ------>  ','FontSize',FS(FSi),'FontWeight','Bold');
    title(H.corr,'Sample Correlation Matrix  ','FontSize',FS(FSi+2),'FontWeight','Bold');
  end

  zoom reset
return

%-----------------------------------------------------------------------
function show_mean_boxplot(data_boxp, name_boxp, quality_order, obj)
%-----------------------------------------------------------------------
  global H pos filename FS FSi sample bp mask1d show_name

  if 0
    % set GUI element
    mid = strfind( spm_str_manip(cellstr(get(H.boxp,'string')),'d'),cellstr(name_boxp));
    mid = find(cellfun('isempty',mid)==0,1,'first');
    set(H.boxp,'Value',mid);
  end
  
  if nargin == 0
    data_boxp     = bp.data;
    name_boxp     = bp.name;
    quality_order = bp.order;
  end

  spm_figure('Clear',H.graphics); 
  spm_figure('Focus',H.graphics); 

  H.boxplot = axes('Position',pos.boxplot,'Parent',H.graphics);
  set(H.graphics,'Renderer','OpenGL','color',[0.95 0.95 0.95]);

  %if isfield(H,'chbox'), nval = ~get(H.chbox,'value'), else nval = 1; end
  H.chbox = uicontrol(H.graphics,...
    'string','Show filenames','Units','normalized',...
    'position',pos.fnamesbox,'callback',@checkbox_names,...
    'Style','CheckBox','HorizontalAlignment','center',...
    'ToolTipString','Show filenames in boxplot','value',show_name,...
    'Interruptible','on','Visible','on','FontSize',FS(FSi));

  n_samples = max(sample);

  xpos = cell(1,n_samples);
  data = cell(1,n_samples);

  allow_violin = 2;

  %% create filenames
  hold on
  for i=1:n_samples
    indtype   = { mask1d' 'k.' [0 0 0] 10; ~mask1d' 'rx' [1 0 0] 3}; 
    gnames{i} = sprintf('S%d',i);
    for ii=1:size(indtype,1)
      ind  = find(sample == i & indtype{ii,1});
      if numel(ind>0)
        datap{i} = data_boxp(ind);
        if ii==1, data{i} = datap{i}; end

        if length(ind) < 8
          allow_violin = 0;
        end

        if n_samples == 1
          xpos{i} = (i-1)+2*(0:length(ind)-1)/(length(ind)-1);
        else
          xpos{i} = 0.5/length(ind) + 0.5+(i-1)+1*(0:length(ind)-1)/(length(ind));
        end

        if get(H.chbox,'value')
          for j=1:length(ind)
            H.fnames{j,i} = text(xpos{i}(j),datap{i}(j),...
              filename.m{ind(j)},'Color',indtype{ii,3},...
              'FontSize',FS(FSi-1),'HorizontalAlignment','center');
          end
        else
          for j=1:length(ind)
            H.fnames{j,i} = plot(xpos{i}(j),datap{i}(j),indtype{ii,2},'MarkerSize',indtype{ii,4});
          end
        end
      end
    end 
  end

  %% create boxplot
  opt = struct('groupnum',0,'ygrid',1,'box',1,'violin',allow_violin,'median',2,...
               'groupcolor',jet(n_samples),'names',{gnames},'xlim',[-.25 n_samples+1.25]); 
  if max(data_boxp) > min(data_boxp)
    ylim_add = 0.075;
    yamp = max(data_boxp) - min(data_boxp);
    ylim_min = min(data_boxp) - ylim_add*yamp;
    ylim_max = max(data_boxp) + ylim_add*yamp; 
    opt.ylim = [ylim_min ylim_max];
  end    
  cat_plot_boxplot(data,opt); box on;


  %% add colored labels and title
  if n_samples > 1
    [tmp,  tmp2] = spm_str_manip(char(filename.s),'C');
    title_str = sprintf('%s ',strrep(tmp,tmp2.s,''));
    %fprintf('\nCommon filename: %s\n',tmp);
  else
    title_str = sprintf('Common filename: %s*',spm_file(char(filename.s),'short25'));
  end
  title({['Boxplot: ' name_boxp],title_str},'FontSize',FS(FSi+1),'FontWeight','Bold');
  xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(10),...
      'FontWeight','Bold');

  xpos = -0.35 - n_samples*0.1;

  if quality_order == -2 
    % reverse order to have the good things allways on the top
    set(gca, 'YDir','reverse');
    quality_order = 1; 
    t = ylim_min; ylim_min = ylim_max; ylim_max = t; 
  end
  if (length(data_boxp) > 2)
    if quality_order > 0 
      text(xpos, ylim_min,'<----- Low rating (poor quality)  ','Color','red','Rotation',...
          90,'HorizontalAlignment','left','FontSize',FS(9),'FontWeight','Bold')
      text(xpos, ylim_max,'High rating (good quality) ------>  ','Color',[0 0.8 0],'Rotation',...
          90,'HorizontalAlignment','right','FontSize',FS(9),'FontWeight','Bold')
    else
      text(xpos, ylim_max,'Low rating (poor quality) ------>  ','Color','red','Rotation',...
          90,'HorizontalAlignment','right','FontSize',FS(9),'FontWeight','Bold')
      text(xpos, ylim_min,'<----- High rating (good quality)  ','Color',[0 0.8 0],'Rotation',...
          90,'HorizontalAlignment','left','FontSize',FS(9),'FontWeight','Bold')
    end
    text(xpos, (ylim_max+ylim_min)/2,sprintf('%s',name_boxp),'Color','black','Rotation',...
          90,'HorizontalAlignment','center','FontSize',FS(9),'FontWeight','Bold')
  end

  hold off

 
  bp = struct('data',data_boxp,'name',name_boxp,'order',quality_order);

return

%-----------------------------------------------------------------------
function update_alpha(obj, event_obj)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
  global H pos img img_alpha isscatter %mesh_detected

  alphaval = get(H.alphabox,'Value');

  % display image with 2nd colorbar (gray)
  image(H.slice,65 + img);
  set(H.slice,'Position',pos.slice  .* [1 1 1 1 - 0.5*isscatter],...
    'HitTest','off','Interruptible','off','visible','off');
  
  % prepare alpha overlays for red and green colors
  if alphaval > 0
    % get 2%/98% ranges of difference image
    range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

    hold(H.slice,'on')
    alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
    alpha_r = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
    hg = image(H.slice,alpha_g); set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),...
      'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
    %if ~mesh_detected, axis image; end
    hr = image(H.slice,alpha_r); set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),...
      'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
    %if ~mesh_detected, axis image; end
    hold(H.slice,'off')
  end

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
  global V Vchanged data_array data_array_diff H pos dataprefix inorm ...
    sorted ind_sorted isscatter names_changed img_alpha img

  if isfield(H,'mm')
    slice_mm = get(H.mm,'Value');
  else
    slice_mm = 0;
  end

  if names_changed
    P = Vchanged;
  else
    P = V;
  end

  vx   =  sqrt(sum(P(1).mat(1:3,1:3).^2));
  Orig = P(1).mat\[0 0 0 1]';
  sl   = round(slice_mm/vx(3)+Orig(3));

  % if slice is outside of image use middle slice
  if (sl>P(1).dim(3)) || (sl<1)
    sl = round(P(1).dim(3)/2);
  end
  set(H.mm,'Value',(sl-Orig(3))*vx(3));

  M  = spm_matrix([0 0 sl]);
  data_array_diff = data_array;

  %%
  for i = 1:length(V)
    img = spm_slice_vol(P(i),M,P(1).dim(1:2),[1 0]);
    img(isnan(img)) = 0;

    % rescue unscaled data
    data_array_diff(:,:,i) = img;
  
    if inorm %&& ~isempty(strfind(dataprefix,'wp'))
      imgscale = median(img(img ~= 0));
      data_array(:,:,i) = img/imgscale;
    else
      data_array(:,:,i) = img;
    end
  end
  %%
  imgscale = median(data_array(data_array ~= 0)); % scale image according to mean
  data_array = data_array / imgscale * 0.3;

  % calculate individual difference to mean image
  for i=1:size(data_array_diff,3)
    data_array_diff(:,:,i) = data_array_diff(:,:,i) - mean(data_array_diff,3);
  end

  % enhance contrast and scale image to 0..64
  %mn = min(data_array(:));
  %mx = max(data_array(:));
  data_array = min(64,max(1,63 * data_array + 1));

  if sorted
    if isfield(pos,'x')
      x = ind_sorted(pos.x);
      if ~isscatter
        y = ind_sorted(pos.y);
      end
    end
  else
    if isfield(pos,'x')
      x = pos.x;
      if ~isscatter
        y = pos.y;
      end
    end
  end

  % check whether mouse position is defined
  if isfield(pos,'x')
    if isscatter
      img       = data_array(:,:,x)';
      img_alpha = data_array_diff(:,:,x)';
    else
      img       = [data_array(:,:,y) data_array(:,:,x)]';
      img_alpha = [data_array_diff(:,:,y) data_array_diff(:,:,x)]';
    end

    % correct orientation
    img = rot90(img,2);
    img_alpha = rot90(img_alpha,2);

    % use gray scale colormap for values > 64
    update_alpha
  
    set(H.mm_txt,'String',sprintf('%0.1f mm',get(H.mm,'Value')));
  end

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj) 
%-----------------------------------------------------------------------
  global trashlist filename sample H X YpY ...
    data_array data_array_diff pos mesh_detected ind_sorted sorted ...
    isscatter img img_alpha mask1d ...
    org_files pdf_files log_files surf_files xml_files

  showtrash = get(H.showtrash,'Value');
  
  alphaval = get(H.alphabox,'Value');

  H.datatip = event_obj; 
  
  pos.pos_mouse = get(event_obj, 'Position');
  pos.tar_mouse = get(event_obj, 'Target');
  
  % Limit the number of datatips:
  % Although it is posible or maybe use to mark multiple objects for
  % comparision of volumes, surfaces and PDFs, the other check cases 
  % XML and LOG-files are worse to handle. 
  % Furthermore, it is elaboritiv to suport all trashlist operations.
  dcmlim = 0; % + 6*isscatter; % up to 6 plot objects?
  dcm = findall(gcf,'Type','hggroup','selected','off'); 
  if numel(dcm)>dcmlim, delete(dcm(dcmlim+1:end)); end
  set(findall(gcf,'Type','hggroup'),'Visible','on')
  
  if isscatter
    pos.x = find(X(:,1) == pos.pos_mouse(1));
    if isempty(pos.x)
      pos.x = find(X(:,2) == pos.pos_mouse(2));
    end

    % text info for data cursor window
    txt = {sprintf('S%d:%s',sample(pos.x),filename.m{pos.x})};

    set(H.slice,'Position',pos.slice .* [1 1 1 0.5],'Visible','on');

    x = pos.x;

  else % covariance matrix
    % check for valid mouse position
    if pos.pos_mouse(1) > pos.pos_mouse(2) || pos.pos_mouse(1)>length(sample) || pos.pos_mouse(2)>length(sample)
      %txt = {''}; 
      set([H.slice,H.mm,H.mm_txt,H.alphabox],'Visible','off');
      set(get(H.slice,'children'),'Visible','off');
      unit = struct2cell(H.checkui); set([unit{cellfun(@ishandle,unit)}],'Enable','off'); 
      set([H.trashui.trash,H.trashui.detrash,H.trashui.trashrow,H.trashui.detrashrow],'Enable','off');
      set(findall(gcf,'Type','hggroup'),'Visible','off')
      return
    end

    % save position of mouse
    if ~showtrash
      pos.x = find(cumsum(mask1d)==pos.pos_mouse(1),1,'first');
      pos.y = find(cumsum(mask1d)==pos.pos_mouse(2),1,'first');
    else
      pos.x = pos.pos_mouse(1);
      pos.y = pos.pos_mouse(2);
    end

    if sorted
      if isfield(pos,'x')
        x = ind_sorted(pos.x);
        y = ind_sorted(pos.y);
      end
    else
      if isfield(pos,'x')
        x = pos.x;
        y = pos.y;
      end
    end

    % text info for data cursor window
    if mesh_detected
      txt = {
        sprintf('Correlation: %3.3f',YpY(x,y)),...
        sprintf('Left:  S%d:%s',sample(x),filename.m{x}),...
        sprintf('Right: S%d:%s',sample(y),filename.m{y})};
    else
      txt = {
        sprintf('Correlation:  %3.3f',YpY(x,y)), ...
        sprintf('Column (Top): S%d:%s',sample(x),filename.m{x}), ...
        sprintf('Row (Bottom): S%d:%s',sample(y),filename.m{y})};
    end

    set(H.slice,'Position',pos.slice,'Visible','on');
  end

  % == check unit ==
  onoff = {'on','off'};
  if isscatter
    set(H.checkui.vol ,'Enable',onoff{ isempty(org_files{pos.x})+1  });
    set(H.checkui.surf,'Enable',onoff{ isempty(surf_files{pos.x})+1 });
    set(H.checkui.xml ,'Enable',onoff{ isempty(xml_files{pos.x})+1  });
    set(H.checkui.log ,'Enable',onoff{ isempty(log_files{pos.x})+1  });
    set(H.checkui.pdf ,'Enable',onoff{ isempty(pdf_files{pos.x})+1  });
  else
    set(H.checkui.vol ,'Enable',onoff{ (isempty(org_files{pos.x})  | isempty(org_files{pos.y}))  + 1 });
    set(H.checkui.surf,'Enable',onoff{ (isempty(surf_files{pos.x}) | isempty(surf_files{pos.y})) + 1 });
    set(H.checkui.xml ,'Enable',onoff{ (isempty(xml_files{pos.x})  | isempty(xml_files{pos.y}))  + 1 });
    set(H.checkui.log ,'Enable',onoff{ (isempty(log_files{pos.x})  | isempty(log_files{pos.y}))  + 1 });
    set(H.checkui.pdf ,'Enable',onoff{ (isempty(pdf_files{pos.x})  | isempty(pdf_files{pos.y}))  + 1 });
  end

  % == trash list unit ==

  if ~isempty(pos.x) 
    if all( trashlist~=pos.x ) 
      set(H.trashui.trash  ,'Enable','on' );
      set(H.trashui.detrash,'Enable','off');
    else
      set(H.trashui.trash  ,'Enable','off');
      set(H.trashui.detrash,'Enable','on' );
    end
  %else
  %  set([H.trashui.trash,H.trashui.detrash],'Enable','off');
  end
  if ~isscatter && isfield(pos,'y') && ~isempty(pos.y)
    if all( trashlist~=pos.y ) 
      set(H.trashui.trashrow  ,'Enable','on');
      set(H.trashui.detrashrow,'Enable','off');
    else
      set(H.trashui.trashrow  ,'Enable','off');
      set(H.trashui.detrashrow,'Enable','on');
    end
  end
  if ~isempty(trashlist)
    set([H.trashui.new,H.trashui.disptrash,H.trashui.ziptrash],'Enable','on');
  else
    set([H.trashui.new,H.trashui.disptrash,H.trashui.ziptrash],'Enable','off');
  end  

  set(H.alphabox,'Visible','on');



  if mesh_detected 
    % use indexed 2D-sheet to display surface data as image
    % check surface size to use indexed 2D map
    if (length(data_array(:,x)) == 163842) || (length(data_array(:,x)) == 32492)
      % combined surface 
      
      if (length(data_array(:,x)) == 163842) 
        ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','fsavg.index2D_256x128.txt'));
      else
        ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','fsavg.index2D_256x128.txt'));
      end
      
      if isscatter
        img = reshape(data_array(ind,x),[256,128]);
      else
        img = [reshape(data_array(ind,x),[256,128]) reshape(data_array(ind,y),[256,128])];
      end
      img = circshift(img,128);
      
      % alpha overlay
      if isscatter
        img_alpha = reshape(data_array_diff(ind,x),[256,128]);
      else
        img_alpha = [reshape(data_array_diff(ind,x),[256,128]) reshape(data_array_diff(ind,y),[256,128])];
      end
      img_alpha = circshift(img_alpha,128);
  
        
    elseif (length(data_array(:,x)) == 327684) || (length(data_array(:,x)) == 64984)
      is32k = length(data_array(:,x)) == 64984;
      if is32k
        atlasdir = 'atlases_surfaces_32k';
        tmpdir   = 'templates_surfaces_32k';
      else
        atlasdir = 'atlases_surfaces';
        tmpdir   = 'templates_surfaces'; 
      end
      lrb = [size(data_array,1)/2,size(data_array,1)/2+1];
      lab = fullfile(spm('dir'),'toolbox','cat12',atlasdir,'lh.aparc_DK40.freesurfer.annot'); 
      ind = spm_load(fullfile(spm('dir'),'toolbox','cat12',tmpdir,'fsavg.index2D_256x128.txt'));
      
      %% average both atlas maps to keep it simpler
      [vr{1},lb{1},tb{1}] = cat_io_FreeSurfer('read_annotation',lab);
      [vr{2},lb{2},tb{2}] = cat_io_FreeSurfer('read_annotation',strrep(lab,'lh.','rh.'));
      atlas_lh = circshift(reshape(lb{1}(ind),[256,128]),64)';
      atlas_rh = circshift(reshape(lb{2}(ind),[256,128]),64)';
      atlas    = cat_vol_median3(single(cat(3,atlas_lh,atlas_rh))); 
      [gx,gy]  = gradient(atlas(:,:,1)); 
      atlasmsk = single(((abs(gx) + abs(gy))./atlas(:,:,1))<0.05); 
      atlasmsk(atlasmsk==0) = nan; 
      atlasmsk = flipud(atlasmsk); 
      
      %% load data array entry
      data_array_x_lh = data_array(1:lrb(2),x);
      data_array_x_rh = data_array(lrb(2):end,x);
      if ~isscatter
        data_array_y_lh = data_array(1:lrb(1),y);
        data_array_y_rh = data_array(lrb(2):end,y);
      end
         
      % rotate the image (') and shift it (64) to have a posterior cutting edge 
      img_x_lh = flipud(circshift(reshape(data_array_x_lh(ind),[256,128]),64)');
      img_x_rh = flipud(circshift(reshape(data_array_x_rh(ind),[256,128]),64)');
      if ~isscatter
        img_y_lh = flipud(circshift(reshape(data_array_y_lh(ind),[256,128]),64)');
        img_y_rh = flipud(circshift(reshape(data_array_y_rh(ind),[256,128]),64)');
      end
      if isscatter
        img = [img_x_lh .* atlasmsk; img_x_rh .* atlasmsk];
      else
        img = [img_x_lh .* atlasmsk; img_y_lh .* atlasmsk; ...
               img_x_rh .* atlasmsk; img_y_rh .* atlasmsk];
      end
     
      %% alpha overlay
      data_array_x_lh = data_array_diff(1:lrb(2),x);
      data_array_x_rh = data_array_diff(lrb(2):end,x);
      if ~isscatter
        data_array_y_lh = data_array_diff(1:lrb(1),y);
        data_array_y_rh = data_array_diff(lrb(2):end,y);
      end
      img_x_lh = flipud(circshift(reshape(data_array_x_lh(ind),[256,128]),64)');
      img_x_rh = flipud(circshift(reshape(data_array_x_rh(ind),[256,128]),64)');
      if ~isscatter
        img_y_lh = flipud(circshift(reshape(data_array_y_lh(ind),[256,128]),64)');
        img_y_rh = flipud(circshift(reshape(data_array_y_rh(ind),[256,128]),64)');
      end
      if isscatter
        img_alpha = [img_x_lh .* atlasmsk; img_x_rh .* fliplr(atlasmsk)];
      else
        img_alpha = [img_x_lh .* atlasmsk; img_y_lh .* fliplr(atlasmsk); ...
                     img_x_rh .* atlasmsk; img_y_rh .* fliplr(atlasmsk)];
      end
      
    else
      if isscatter
        img = data_array(:,x)';
        % alpha overlay
        img_alpha = data_array_diff(:,x)';
      else
        img = [data_array(:,y) data_array(:,x)]';
        % alpha overlay
        img_alpha = [data_array_diff(:,y) data_array_diff(:,x)]';
      end
    end

    % scale img to 0..64
    sd = std(data_array(:));
    mn = -sd*2;
    mx = sd*2;
    img = 64*((img - mn)/(mx-mn));
  else
    % add slider for colume data
    set(H.mm,'Visible','on');
    set(H.mm_txt,'Visible','on');
    if isscatter
      img = data_array(:,:,x)';
      % alpha overlay
      img_alpha = data_array_diff(:,:,x)';
    else
      img = [data_array(:,:,y) data_array(:,:,x)]';
      % alpha overlay
      img_alpha = [data_array_diff(:,:,y) data_array_diff(:,:,x)]';
    end
  end

%  set(H.cbar,'TickLength',[0 0],'XTickLabel',linspace(0.7,1.0,7)); %round(100*linspace(min(YpY(mask2d(:))),max(YpY(mask2d(:))),5))/100);


  % correct orientation
  img = rot90(img,2);
  img_alpha = rot90(img_alpha,2);

  % display image with 2nd colorbar (gray)
  image(H.slice,65 + img); 
  set(H.slice,'visible','off');
  
  % prepare alpha overlays for red and green colors
  if alphaval > 0
    %% get 2%/98% ranges of difference image
    range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

    hold(H.slice,'on');
    alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
    alpha_r   = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
    hg = image(H.slice,alpha_g); 
    set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),...
      'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
    %if ~mesh_detected, axis image; end
    hr = image(H.slice,alpha_r); 
    set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),...
      'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
    %if ~mesh_detected, axis image; end
    hold(H.slice,'off');
  end

  if mesh_detected
    xlabel('2D surface maps');
  end

return

%-----------------------------------------------------------------------
function varargout = cat_tst_qa_cleaner_intern(data,opt)
%% THIS FUNCTION IS TO FAT - SEPARATE AND CLEAN IT! 
%  Do not forget to remove old external version from SVN if this is done.
%  _____________________________________________________________________
%  Estimate quality grades of given rating of one (or more) protocols
%  with 2 to 6 grads to separate passed, (unassignable) and failed 
%  images, by finding the first peak in the image quality histogram  
%  and using its width (standard deviation) in a limited range. 
%  If multiple protocols are used, than use the site variable opt.site 
%  and use the site depending output rths.
%
%  The passed range can be variated by opt.cf with lower values for harder 
%  and higher values for softer thresholds (more passed images), where 
%  opt.cf=1, describes a range that is similar to about 1% BWP noise that 
%  is equal to 5 rps.
%  ROC evaluation showed that opt.cf=0.72 allows the best separation of 
%  images without and with artifacts, but if the majority of your data 
%  include light artifacts (e.g. by movements in young children) that 
%  a softer weighing, e.g. opt.cf=2, is preferable (maximum is 4). 
%
%  Use the selftest with randomly generated data to get a first impression:
%    cat_tst_qa_cleaner('test')
%  _____________________________________________________________________
%
%  This tool is still in development / undert test:
%   * the combination of different sites is not finished
%   * multiside output required a 'stacked' output
%
%  [Pth,rth,sq,rths,rthsc,sqs] = cat_tst_qa_remover(data[,opt])
%
%    Pth      .. global threshold for passed images 
%                (for odd grades this is in the middle of the unassignable)
%    rth      .. all global threshold(s) between grads
%    sq       .. estimated first peak and its std, where the std depend on
%                the number of grades!
%    rths     .. site depending thresholds between grads of each input 
%    rthsc    .. site depending thresholds between grads of each input 
%                (global corrected, removed further low quality data)
%    sqs      .. site depending first peaks and stds of passed data 
%
%    data     .. array of quality ratings or xml-files
%    opt      .. option structure
%     .grads  .. number of grads (2:6, default=6, see below)
%     .cf     .. factor for harder/softer thresholds (defaults=0.72)
%     .figure .. display histogramm with colored ranges of grads
%                 1 - use current figure
%                 2 - create new figure (default)
%                 3 - use one test figure (default in the selftest)
%  _____________________________________________________________________
%
%  Grades:
%    2 grads:
%      P   passed
%      F   failed
%    3 grads:
%      P   passed
%      U   unassignable
%      F   failed
%    4 grads:
%      P+  clear passed 
%      P-  just passed
%      F+  just failed
%      F-  clear failed
%    5 grads:
%      P+  clear passed 
%      P-  just passed
%      U   unassignable
%      F+  just failed
%      F-  clear failed
%    6 grads (default):
%      P+  clear passed 
%      P   passed 
%      P-  just passed
%      F+  just failed
%      F   failed
%      F-  clear failed
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_stat_check_cov2.m 1352 2018-08-10 16:11:27Z dahnke $ 

  global figcolor

  clear th; 
  if ~exist('opt','var'), opt = struct(); end
  def.cf        = 0.72;                 % normalization factor for rating 
  def.grads     = 6;                    % number of grads (default = 6)
  def.model     = 1;                    % model used for rating
  def.figure    = 2;                    % figure=2 for new/own figure
  def.smooth    = 0;                    % smoothing of output data
  def.siterf    = 1000000;              % round factor to identify similar resolution level 
  def.siteavgperc = [0.10 0.90];        % ?
  opt = cat_io_checkinopt(opt,def); 
  opt.cf = max( 0 , min( 4 , opt.cf )); % limit of cf
  
  % test options
  %opt.model = 2;
  %opt.grads = 6;
  
  % if no intput is given use SPM select to get some xml-files
  if ~exist('data','var') || isempty(data)
    data = cellstr(spm_select(inf,'XML','select qa XML-files',{},pwd,'^cat_.*')); 
  elseif ischar(data)
    data = cellstr(data);
  end
  if isempty(data) || (iscell(data) && all(cellfun('isempty',data)))
    if nargout>=1, varargout{1} = 3; end
    if nargout>=2, varargout{2} = 3; end
    if nargout>=3, varargout{3} = [2.5 0.5]; end
    if nargout>=4, varargout{4} = 3*ones(size(data)); end
    if nargout>=5, varargout{5} = 3*ones(size(data)); end
    if nargout>=6, varargout{6} = repmat([2.5 0.5],numel(data),1); end

    return;
  end
  if iscell(data) && ~strcmp(data{1},'test')
    fprintf('Load XML data');
    P = data; 
    xml = cat_io_xml(data,struct(),'read',1); clear data; 
    for di=1:numel(xml)
      opt.site(di,1) = xml(di).qualityratings.res_RMS; 
      data(di,1)     = xml(di).qualityratings.NCR; 
    end,
  end
  

  % --------------------------------------------------------------------
  % If a site variable is given (e.g. by the RMS resolution) then call
  % the cleanup for each subset. The threshold will be collected in a 
  % vector [markthss x opt.grads] with the same length as data. 
  % Nevertheless an average threshold will is estimated as average of 
  % the percentual range give by opt.siteavgperc with e.g. [0.1 0.9] to
  % concider 80% of the data.
  %  -------------------------------------------------------------------
  if isfield(opt,'site')
    if numel(opt.site)~=numel(data),
      error('cat_tst_qa_cleaner_intern:numelsitedata','Numer of elements in data and opt.site have to be equal.\n');
    end
    opt.site = round(opt.site*opt.siterf)/opt.siterf; 
    sites    = unique(opt.site); 
    markth   = zeros(numel(sites),opt.grads-1); 
    markths  = zeros(numel(data),opt.grads-1); 
    siteth   = zeros(numel(data),2); 
    for si=1:numel(sites)
      sdatai = find(opt.site==sites(si));
      opts = opt; 
      opts = rmfield(opts,'site');
      opts.figure = 0; 
      [Sth,markth(si,:),out{1:4}] = cat_tst_qa_cleaner_intern(data(sdatai),opts); 
      markths(sdatai,:) = repmat(markth(si,:),numel(sdatai),1); 
      siteth(sdatai,:)  = out{4}; 
    end
    % estimate global threshold
    markthss = sortrows(markth);
    th = cat_stat_nanmean(markthss(max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(1)))):...
                                   max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(2)))),:),1); 
    sd  = out{3}; 
    thx = out{4}; 
    % modify local rating based on the global one                 
    markths2 = markths;
    markths2 = min(markths2,1.2*repmat(th,size(markths2,1),1)); % higher thresholds even for sides with low rating 
    markths2 = max(markths2,0.8*repmat(th,size(markths2,1),1)); % lower  thresholds even for sides with high rating 
    d  = data; 

  else
    %  -----------------------------------------------------------------
    %  Simulate data, if no data is given by several normal distributed
    %  random numbers.
    %  -----------------------------------------------------------------
    if exist('data','var') && ~(iscell(data) && strcmp(data{1},'test'))
      d = data; 
      if numel(d)==0, 
        if nargout>=1, varargout{1} = nan; end
        if nargout>=2, varargout{2} = nan(1,opt.grads); end
        if nargout>=3, varargout{3} = nan(1,2); end
        if nargout>=4, varargout{4} = nan(size(data)); end
        if nargout>=5, varargout{5} = nan(size(data)); end
        if nargout>=6, varargout{6} = nan(size(data)); end
        return;
      end
    elseif iscell(data) && strcmp(data{1},'test')
      % Testcases with different quality ratings
      scans      = 100; % number of scans (per site) for simulation
      testcase   = round(rand(1)*10);
      randoffset = 0.5*randn(1,4);

      switch testcase
        case 0 % good quality, no outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.15)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.03)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.02))];
         case 1 % good quality, with average outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.40)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.15)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.05))];
        case 2 % good-average quality, with outlier group 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];
        case 3 % good-average quality, without outlier group 
          d = [2.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               3.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))]; 
        case 4 % average to low quality, with light falloff  
          d = [3.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               3.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 5 % high to good quality, with light falloff  
          d = [1.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               2.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 6 % high quality, no outlier
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.13)), ...
               3.0 + randoffset(3) + 0.3*randn(1,round(scans*0.05)), ...
               5.0 + randoffset(4) + 0.3*randn(1,round(scans*0.02))];
        case 7 % good quality with second average peak 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];   
        case 8 % good quality with second low quality peak 
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(2) + 0.2*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];    
        case 9 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.2*randn(1,round(scans*0.20)), ...
               3.0 + randoffset(2) + 0.3*randn(1,round(scans*0.20)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.0 + randoffset(4) + 0.8*randn(1,round(scans*0.50))];          
        case 10 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.10)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(4) + 1.0*randn(1,round(scans*0.60))];           
      end

      % remove high quality outlier and set them to normal
      cor = max(1,median(d)-std(d)/2);
      md= d<(cor); d(md) = cor + 0.05*randn(1,sum(md));
      
      % set selftest figure
      opt.figure = 3; 
    end

    
  %% Models
  %  -------------------------------------------------------------------
  %  I start with several ideas that all based on a similar idea: to 
  %  find the first peak that is given by the subset of images without
  %  inferences and to use the variance of this peak for further scaling
  %  of subsets for other grads. As far as IQR is already scaled, we 
  %  can limit the variance value ... e.g. the rating has an error of 
  %  0-2 rps (0.0-0.2 mark points) that is very low for high-quality data
  %  and higher for low-quality data. Due to our the general subdivion 
  %  of the rating scale in +,o, and - (e.g. B+,B,B-) we got a subrange 
  %  of 3.33 rps (1/3 mark points) that gives some kind of upper limit.
  %  -------------------------------------------------------------------
    thx = nan; sd = nan; th = zeros(1,opt.grads-1); 
    switch opt.model
      case 0
        % only global thresholding ... 
        % this is just to use the color bar output 
        thx = 3; 
        sd  = 1; 
        th  = 1.5:1:100;
        th(6:end) = []; 
      case 1 
        % kmeans model:
        % * estimate peaks based on the histogram
        % * mix the first and second peak until it fits to 30% of the data 
        %   or until the number of loops is similar the number of peaks 
        % * use the std give by one BWP noise level (0.5) to describe the 
        %   variance the passed interval.
        
        hx = hist(d,0.5:1:5.5);
        peaks = sum(hx>(max(hx)/5))*3;
        [thx,sdx] = kmeans3D(d,peaks); sdx = sdx./thx;
        for i=1:peaks
          if sum(d<thx(i))/numel(d) < 0.3
            thx(1) = cat_stat_nanmean(thx(1:2));
            sdx(1) = cat_stat_nanstd(d(d<thx(1)));
          end
        end
        sd    = 0.25 / (opt.grads/2) * opt.cf; % 0.5 = 1% BWP noise
        th(1) = thx(1) - sdx(1) + 2*sd(1); %- mean(sdx(1:min(3,numel(sdx)))) 
        for i = 2:opt.grads-1
          th(i) = th(i-1) + 2*sd(1); % 
        end
      case 2
        % similar to case 1, but with std optimization based on the data 
        % ... surprisingly the simple model 1 works better
        
        hx = hist(d,0.5:1:5.5); 
        %for i=1:1, hx(2:end-1) = cat_stat_nanmean(cat(1,hx(1:end-2),hx(2:end-1),hx(3:end)),1); end
        peaks = sum(hx>(max(hx)/5))*3;
        [thx,sdx] = kmeans3D(d,peaks); sdx = sdx./thx;
        for i=1:peaks
          %if numel(thx)>i && sum(d<thx(i))/numel(d) < 0.05
          %  thx(1) = []; sdx(1) = [];
          if sum(d<thx(i))/numel(d) < 0.3 %numel(thx)>i && 
            thx(1) = cat_stat_nanmean(thx(1:2)); 
            sdx(1) = cat_stat_nanstd(d(d<thx(1)));
          end
        end
        sdx(1) = cat_stat_nanstd(d(d<thx(1)));
        [thx,sdx] = kmeans3D(d(d<=(max([min(d),thx(1)+sdx(1)]))),3); thx=thx(2); sdx=sdx(2);  %sdx = sdx./thx;
        sd    = min(1/3,max(1/6,sdx(1))) / (opt.grads/2) * opt.cf; % 0.5 = 1% BWP noise*16
        th(1) = thx(1) - sdx(1) + 2*sd(1);
        for i = 2:opt.grads-1
          th(i) = th(i-1) + 2*sd(1); % 2*2/3*
        end
      
    end
    
    markths  = repmat(mean(th(floor(opt.grads/2):ceil(opt.grads/2))),size(data));
    markths2 = markths;
    siteth   = repmat([thx(1) sd],numel(data),1); 
  end
  
  
%% Print
%  ---------------------------------------------------------------------
%  This part is just for to plot a colorated histogram and the percents
%  of images in each group.
%  ---------------------------------------------------------------------
  if opt.figure
    if opt.figure==2
      f = figure;
      set(f,'color','w')
    elseif opt.figure==3
      f = findobj('type','figure','name','qa_cleaner_test');
      if isempty(f), figure('name','qa_cleaner_test'); else figure(f(1)); clf(f(1)); end
    end
    box on;
    
    %figure
    ss = 0.05; 
    [h,r]  = hist(d,0.5:ss:10.5); 
    for i=1:opt.smooth, h(2:end-1) = cat_stat_nanmean(cat(1,h(1:end-2),h(2:end-1),h(3:end)),1); end
    sh = 1; %sum(h);
    
    % background histogram (all data)
    %bar(r,h/sh,'facecolor',figcolor,'edgecolor','none');
    %fill(r,h/sh,figcolor,'edgecolor','none');
    hold on
    
    yl = [0 max(h)+1]; ylim(yl);
    % main grid
    for i=1.5:6,       plot([i i],ylim,'color',figcolor); end
    switch numel(th)
      case 1
        hx = h; hx(r> th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(1))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 2
        hx = h; hx(r>=th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1) | r>th(2)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% unassignable' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 3
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(2))/numel(d)*100),'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% failed+',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed-',sum(d>=th(3))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 4
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(4)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(2))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% check' ,sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100)          ,'color',[0.7  0.0  0.0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.71,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.67,sprintf('%5.2f%% unassignable',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.63,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.59,sprintf('%5.2f%% failed-',sum(d>=th(4))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 5
        % plot
        testbar=0; % it would be cool to use bars but they failed at least in MATLAB R2013 and killed the axis positions...
        if testbar==1
          hx = h; hx(r>=th(1)+ss) = 0;           bar(r,hx/sh,'facecolor',[0.0  0.5  0.2],'edgecolor','none','barwidth',1);
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; bar(r,hx/sh,'facecolor',[0.4  0.7  0.1],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; bar(r,hx/sh,'facecolor',[0.7  0.8  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; bar(r,hx/sh,'facecolor',[0.9  0.6  0.4],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; bar(r,hx/sh,'facecolor',[0.75 0.3  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(5)-ss) = 0;           bar(r,hx/sh,'facecolor',[0.6  0.15 0.1],'edgecolor','none','barwidth',1);      
        else
          hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(5)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none'); 
        end
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(3))/numel(d)*100) ,'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% passed-',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.55,sprintf('%5.2f%% failed' ,sum(d>=th(4) & d<th(5))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.50,sprintf('%5.2f%% failed-',sum(d>=th(5))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
    end
    xlim([min(r),6.5]); 
    
    % subgrid
    for i=5/6:1/3:6.4, plot([i i],[0 0.03]*max(ylim),'color',[0.2 0.2 0.2]); end
        
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    
    
    % colored main grads
    FS = get(gca,'Fontsize')*1.3;
    set(gca,'XTick',0.5:1:6.5,'XTickLabel',{'100','90','80','70','60','50','40'},'TickLength',[0.02 0.02]);
    % further color axis objects...
    axA = copyobj(gca,gcf); axB = copyobj(axA,gcf); axC = copyobj(gca,gcf); 
    axD = copyobj(gca,gcf); axE = copyobj(gca,gcf); axF = copyobj(gca,gcf);
    % set colors...
    set(axA,'YTick',[],'XTickLabel',{},'XTick',1,'XColor',color(QMC,1),'Color','none','XTicklabel','A','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axB,'YTick',[],'XTickLabel',{},'XTick',2,'XColor',color(QMC,2),'Color','none','XTicklabel','B','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axC,'YTick',[],'XTickLabel',{},'XTick',3,'XColor',color(QMC,3),'Color','none','XTicklabel','C','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axD,'YTick',[],'XTickLabel',{},'XTick',4,'XColor',color(QMC,4),'Color','none','XTicklabel','D','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axE,'YTick',[],'XTickLabel',{},'XTick',5,'XColor',color(QMC,5),'Color','none','XTicklabel','E','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axF,'YTick',[],'XTickLabel',{},'XTick',6,'XColor',color(QMC,6),'Color','none','XTicklabel','F','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    hold off; 
    
    if isfield(opt,'site') && numel(sites>1);
      title(sprintf('Histogram (cf=%0.2f) - global treshold for multisite output (n=%d)',opt.cf,numel(sites)),'Fontsize',FS);
    else
      title(sprintf('Histogram (cf=%0.2f)',opt.cf),'Fontsize',FS);
    end
    xlabel('IQR (rps)','Fontsize',FS); 
    ylabel('number of scans','Fontsize',FS); 
  end
  %%
  MarkColor = cat_io_colormaps('marks+',40); 
  if isfield(opt,'site') && numel(sites)>1, globcorr = ' (global corrected)'; else globcorr = ''; end
  if exist('P','var')
    files = P(data<=markths2(:,3)); 
    fprintf('PASSED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 0
      iqrs = [xml(data<=markths2(:,3)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    else
      
    end
    
    % bad files ...
    files = P(data>markths2(:,3) & data<=markths2(:,4)); 
    fprintf('FAILED+%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,3) & data<=markths2(:,4)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,4) & data<=markths2(:,5)); 
    iqrs  = [xml(data>markths2(:,4) & data<=markths2(:,5)).qualityratings];
    if 1
      fprintf('FAILED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,5)); 
    fprintf('FAILED-%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,5)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
  end
  
  
  %% create output
  if nargout>=1, varargout{1} = mean(th(floor(opt.grads/2):ceil(opt.grads/2))); end
  if nargout>=2, varargout{2} = th; end
  if nargout>=3, varargout{3} = [thx(1) sd(1)]; end
  if nargout>=4, varargout{4} = markths;  end
  if nargout>=5, varargout{5} = markths2; end
  if nargout>=6, varargout{6} = siteth; end
  if nargout>=7, varargout{7} = sites; end
  
  if 0
    %%
    b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
                'value',zeta, 'min',0, 'max',1);
    bgcolor = f.Color;
    bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                    'String','0','BackgroundColor',bgcolor);
    bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                    'String','1','BackgroundColor',bgcolor);
    bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                    'String','Damping Ratio','BackgroundColor',bgcolor);
  end
return
