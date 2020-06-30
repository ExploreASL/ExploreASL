function varargout = cat_stat_check_cov2(job)
% cat_stat_check_cov 
% _________________________________________________________________________
% Function to check covariance and image quality across samples.
% Use the CAT GUI or SPM batch mod for initialization.
%
% Selected volumes have to be in the same orientation with same voxel 
% size and dimension (e.g., spatially registered images), whereas surfaces  
% have to be same size (number of vertices; i.e., resampled surfaces).
%
% varargout = cat_stat_check_cov(job)
%
%  job           .. spm job structure
%    .data_vol   .. volume input files
%    .data_surf  .. surface input files 
%    .c          .. nuisance variables (cell)
%   [.gap]       .. gap between slices (in case of volume input)
%  varargout     ..
%
% _________________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id: cat_stat_check_cov2.m 1612 2020-05-03 12:18:52Z gaser $

%#ok<*AGROW,*ASGLU,*TRYNC,*MINcscc.data.V,*INUSD,*INUSL,*MINV>


% ------------------------------------------------------------------------- 
% Extra development
%
%	- (X) CAT version & parameter control
%       - Show cat-version in datatip 
%         > required version management update in cat_tst_qa
%       - Use CAT version number as nuisance variable 
%         > ask Christian .. maybe as flag
%	- (X) Create CAT segmentation parameter as nuisance variable:
%       n-dimensional distance of biasstr, APP, GCUT+cleanup, LASstr, regstr?
%
% ------------------------------------------------------------------------- 
%
% Development: 
%	- (5) Load surfaces in graphic window or close to it (sorted). 
%	- (5) Data trashing 
%       - message window if no images were selected
% - (5) Colorbar: 
%       - Define Colorbar (min-max, auto=sd-factor), (+-)-buttons?
%       - Auto/fixed colorbar in scatterplot > redraw function?
%       - sample / protocol / global scaling ...
%       > scatter plot update function
%	- (4) Use resampled surface if possible (faster) 
%       - if more than one > which one? > fastest & less smoothed
% - (1) Autotrash > groups with percentual number of files ?	
% - (1) add LEFT/RIGHT, COLUMN/ROW legend to surface slice print
%
% ------------------------------------------------------------------------- 
%
% Possible extensions:
% - (2) Options button with a submenu to set default display options and behavior 
%   		- open original vs. normalized volumes / surfaces
%     	- affine registered volumes
%       - use global scale for processing / visualization (initial parameter)
%       - display datatip variables 
%       - choose colormap
%       - choose window colors
% - (0) Icons for (sorted) Matrixplot, Maha Dist, Worst MNC Cases, 
%       (Sub)Boxplots (MNC, QA, nuisance) ?
%	- (1) overlay with half transparency for 1-2 SD
% - (1) Data viewing GUI	
%       - open all similar objects in scatter plot? - Use colored edges?
%       - (5) show load progress in case of surface display 
%             > not so easy > cat_surf_render!
%       - (4) Button to enlarge the slice-figure in the SPM graphics window 
%             (close to slider)
%       - (2) multiview of highly correlated datasets in scatterplot with 
%             modified check icons (Button for on/off?)
% - (1) try to get all figure handles even in case of additional surface 
%       figures (and not try to close figure 2:26)
% - (1) error message if all scans of a group are missing 
%	- (1) mark critical selections, e.g of protokol/sample (red background) with redraw
%	- (0) more space for the title of the boxplot (in case of complex group paths)
%
% ------------------------------------------------------------------------- 
% 
% Rejected ideas:
% - Use ROI-files to check regions? 
% 	> This would be a separate call of Checkcov that is not required yet.
%
% - Longitudinal mode?
% 	> Too complex right now due to manifold data structure and varying
%    	number of scans.
%
% - Use multiple datatips? 
% 	> Multiple datatips would be nice to view and (de)exclude multiple  
%     objects but will have unclear behavior for the slice window. 
% 	> Deselection is not intuitive! The required number variates (many for
%     delete, less for view).
% 	> This would be elaborate and confusing and it is easier to use only 
%     one at all!
%
% - Save/load functions (buttons) for the selection list?
% 	> This requires  a strong interaction between the disk data and virtual
%     data structure. 
% 	> Elaborative and not really important.  
%
% - Boxplot dependency for ?Data selection?, i.e. show only data visible
%   in matrix/scatter plot?
% 	> No, because this generally removes the samples and this is not required. 
%
% - Use of a table with parameters (MNC, IQR, PIQR, MD, nuisance, autotrash, 
%   tissue volume, cat_pp_version, cat_pp_para, Euler, ...) in an extra 
%   figure or the SPM Graphics window?
% 	> No, this is not required and elaborative, because Checkcov is a graphic
%     tool that already combines such information in abstract figures 
%     (covar matrix, Mahanalobis distance).
%
% - Show autotrash in datatip
% 	> No, because it depend on further autotrash options
%
% - Display additional data (IQR,MNC, ?) in/close to sliceplot
% 	> No, the data was added to the datatip. Overall there is not enough 
%     space for this, requiring an additional figure ...
%
% ------------------------------------------------------------------------- 

  if nargin == 0, error('No argument given.'); end
  
  % remove old figure
  oldfig = findobj('type','figure','number',2); 
  if ~isempty(oldfig); delete(oldfig); end
  
  % create default SPM windows if required        
  if isempty(spm_figure('FindWin','Interactive'))
    spm('createintwin'); 
  end
  if isempty(spm_figure('FindWin','Graphics'))
    spm_figure('Create','Graphics',sprintf('%s: Graphics',spm('version'))); 
  end   
  
  cat_io_cprintf('err','\nWARNING: cat_stat_check_cov2 is in development!\n\n')
  
  
% ------------------------------------------------------------------------- 
% Overview of the main global datastructure of cscc of cat_stat_check_cov 
% ------------------------------------------------------------------------- 
%  cscc .. cat_stat_check_cov (cscc) data structure as short unique global 
%          variable that nobody else use
%    .H             	.. object/button/image/axis handles
%      .mesh_detected .. volume vs. surface mode (=1)
%      .isscatter     .. Covariance matrix vs. Mahanalobis plot (=1)
%      .isxml         .. CAT XML are availabl
%      .sorted        .. sorted Covariance matrix?
%      .show_name     .. show names in boxplot
%      .inorm         .. normalize image intensity (for display) ???
%      .cbar        .. Colorbar
%      ...          .. see create_figure subfunction
%    .pos           .. position values for GUI objects/buttons   
%      ...          .. see initialization below
%    .data
%      .YpY         .. covariance matrix 
%      .QM          .. Quality measures from xml files
%      .QMnames     .. name of the cscc.data.QM rows
%      .mean_cov    .. mean covariance matrix 
%      .img         .. slice image(s) for GUI display
%      .img_alpha   .. slice image(s) for GUI display
%      .cata        .. minimum/maximum value of the mesh
%      .V           .. header of the input files (volume/surface)
%      .Vo          .. header of the original volume files (used for
%                      preprocessing
%      .img_alpa    ..
%      *Vchanged       .. modified volume header for 4D-structures
%      *Vchanged_names .. modified volume header for 4D-structures
%      .X           .. data structure for Mahalobis distance
%      .X2          .. data structure for Mahalobis distance
%      .MD          .. data structure for Mahalobis distance
%      .MD2         .. data structure for Mahalobis distance
%    .files         .. CAT preprocessing files 
%      .data        .. normalized input files of cat_stat_check_cov
%                      (e.g. wmp1, thickness, curv, ...)
%      .org         .. original files used for preprocessing
%      .surf        .. surface files (thickness, mesh)
%      .surfr       .. resampled files (thickness, mesh)
%      .xml         .. cscc.data.XML data of CAT preprocessing
%      .log         .. log-file of CAT preprocessing
%      .pdf         .. pdf-file of CAT preprocessing (report figure) 
%      .jpg         .. pdf-file of CAT preprocessing (report figure)
%      .fname       .. structure from spm_str_manip with grouped filenames
%      .dataprefix  .. only the prefix of the data input
%    .select
%      .trashlist   .. index list of subjects to remove 
%      .trashhist   .. undo trash list 
%      .trashhistf  .. redo trash list
%      .trash1d     .. 1D mask for files on the trash list
%      .trash2d     .. 2D mask for files on the trash list 
%      .samp1d      .. 1D mask for cscc.datagroups.sample view 
%      .samp2d      .. 2D mask for cscc.datagroups.sample view
%      .prot1d      .. 1D mask for cscc.datagroups.protocol view 
%      .prot2d      .. 2D mask for cscc.datagroups.protocol view 
%    .datagroups
%      .n_samples
%      .sample
%      .protocol
%      .protocols 
%    .display       .. display and print variables
%      .WS          .. SPM window size
%      .FS          .. list of SPM font size
%      .FSi         .. main selection of SPM font size
%      .figcolor    .. background color of main figure
%      .useicons    .. use button icons (use JAcscc.data.VA workaround)
% -------------------------------------------------------------------------  
  
  clearvars -GLOBAL cscc;       % clear old
  global cscc;                  % create new
  cscc.job = job;
  cscc.job.expertgui   = cat_get_defaults('extopts.expertgui');
  cscc.H = struct();
  cscc.H.cbarfix.Value = 0;     % inactive
  
  
  cscc.display = struct( ...             display and print variables
    'WS', spm('Winsize','Graphics'), ... SPM window size
    'FS', spm('FontSizes'), ...          list of SPM font size
    'FSi', 8, ...                        main selection of SPM font size
    'figcolor',[0.8 0.8 0.8], ...        background color of main figure
    'useicons',1); %                     use button icons (use JAcscc.data.VA workaround!)

  
  cscc.select.trashlist     = []; % start with empty list 
  cscc.select.trashhist     = []; % history of trash operations (for undo)
  cscc.select.trashhistf    = []; % history of trash operations (for redo)
  
  
  cscc.H.sorted             = 0;  % show data by file order
  cscc.H.isscatter          = 0;  % active GUI surface plot
  cscc.H.show_name          = 0;  % show filenames in boxplot rather small dots
  cscc.H.inorm              = 1;  % normalize slice intensity even in normalized data
 
  
  % volume or surfaces
  if isfield(job,'data_vol')
    datafield     = 'data_vol'; 
    datadir       = 'mri'; 
    cscc.H.mesh_detected = 0;
  elseif isfield(job,'data_surf')
    datafield     = 'data_surf'; 
    datadir       = 'surf'; 
    cscc.H.mesh_detected = 1;
  end
  cscc.files.datafield = datafield; 
  
  
  % positions & global font size option
  cscc.display.WP = get(spm_figure('FindWin','Graphics'),'Position'); 
  if cscc.display.WP(2)>0 && cscc.display.WP(2)<50
    cscc.display.WP = cscc.display.WP(2);
  else
    cscc.display.WP = 10;
  end
  
  cscc.H.show_violin = 1;  
  
  popb = [0.038 0.035];  % size of the small buttons
  popm = 0.780;          % x-position of the control elements
  
  
  cscc.posp = struct('naviui',0.835,'trashui',0.725,'checkui',0.780); % y-pos of major control elements
  cscc.pos  = struct(...
    'fig',            [cscc.display.WP(1) cscc.display.WP(1) ...
                       1.4*cscc.display.WS(3) 1.2*cscc.display.WS(3)],... % figure
    'popup',          [10  10  200       100  ],... % popup in case of closing with non-empty trash list
    'trashpopup',     [0.7*cscc.display.WS(3) 0.4*cscc.display.WS(3) ...
                       1.2*cscc.display.WS(3) 0.9*cscc.display.WS(3)],... % figure
    ...
    'corr',           [0.045 0.050 0.700 0.820],... % correlation matrix
    'slice',          [0.780 0.060 0.190 0.450] - cscc.H.mesh_detected*[0 0.01 0 0],... % image plot
    ...'surfi',          [0.780 0.050 0.190 0.560],... % image plot
    'cbar',           [0.045 0.950 0.700 0.020],... % colorbar for correlation matrix
    ...'cbar',           [0.045 0.950 0.580 0.020],... % colorbar for correlation matrix (active cbarfix)
    ...'cbarfix',        [0.657 0.943 0.100 0.030],... % colorbar fix/auto option (active cbarfix)  
    ... see also  cscc.H.cbarfix.Value  and  checkbox_cbarfix  function!
    ... 
    'boxplot',        [0.100 0.055 0.880 0.915],... % boxplot axis
    'refresh',        [0.100 0.003 0.055 0.032],... % refresh boxplot
    'worst',          [0.150 0.003 0.200 0.032],... % show worst (in boxplot) 
    'fnamesbox',      [0.830 0.001 0.160 0.032],... % show filenames in boxplot 
    'plotbox',        [0.830 0.022 0.160 0.032],... % switch between boxplot and violin plot 
    ...
    'close',          [0.775 0.935 0.100 0.040],... % close button
    'help',           [0.875 0.935 0.100 0.040],... % help button
    ...
    'sort',           [0.772 0.880 0.110 0.050],... % list to use ordered matrix or Maha-distance 
    'boxp',           [0.872 0.880 0.110 0.050],... % list to display different variables as boxplot
    'samp',           [0.772 0.615 0.110 0.055],... % list to use ordered matrix or Maha-distance 
    'prot',           [0.872 0.615 0.110 0.055],... % list to display different variables as boxplot
    'showtrash',      [0.776 0.610 0.110 0.025],... % colorbar fix/auto option
    ...
    'alphabox',       [0.775 -0.001 0.200 0.030] + cscc.H.mesh_detected*[0 0.016 0 0],... % show filenames in boxplot 
    'sslider',        [0.780 0.030 0.193 0.040],... % slider for z-slice  
    ...
    ... == navigation unit ==
    'scSelect',       [popm+popb(1)*0 cscc.posp.naviui popb],... % select (default) 
    'scZoomReset',    [popm+popb(1)*1 cscc.posp.naviui popb],... % standard zoom
    'scZoomIn',       [popm+popb(1)*2 cscc.posp.naviui popb],... % zoom in 
    'scZoomOut',      [popm+popb(1)*3 cscc.posp.naviui popb],... % zoom out
    'scPan',          [popm+popb(1)*4 cscc.posp.naviui popb],... % pan (moving hand)
    ...
    ... == tashlist unit ==
    'newtrash',       [popm+popb(1)*0 cscc.posp.trashui popb],... % new trash list
    'disptrash',      [popm+popb(1)*1 cscc.posp.trashui popb],... % print trash list
    'trash',          [popm+popb(1)*2 cscc.posp.trashui popb],... % add data to trash list
    'detrash',        [popm+popb(1)*3 cscc.posp.trashui popb],... % remove data from trash list
    'autotrash',      [popm+popb(1)*4 cscc.posp.trashui popb],... % button to mark data with low IQR
    ... second row
    'undo',           [popm+popb(1)*0 cscc.posp.trashui-popb(2) popb],... % undo last trash list operation
    'redo',           [popm+popb(1)*1 cscc.posp.trashui-popb(2) popb],... % redo last trash list operation
    'trashrow',       [popm+popb(1)*2 cscc.posp.trashui-popb(2) popb],... % add data to trash list
    'detrashrow',     [popm+popb(1)*3 cscc.posp.trashui-popb(2) popb],... % button to mark data with low IQR
    'ziptrash',       [popm+popb(1)*4 cscc.posp.trashui-popb(2) popb],... % pack data from trash list
    ...
    ... == checklist unit ==
    'checkvol',       [popm+popb(1)*0 cscc.posp.checkui popb],... % open checkvol 
    'checksurf',      [popm+popb(1)*1 cscc.posp.checkui popb],... % open checksurf
    'checklog',       [popm+popb(1)*2 cscc.posp.checkui popb],... % open log-txt
    'checkxml',       [popm+popb(1)*3 cscc.posp.checkui popb],... % open xml-txt
    'checkpdf',       [popm+popb(1)*4 cscc.posp.checkui popb]);   % open pdf in external viewer
    % checksurfn?

  


  %% get all filenames from the data_vol/surf input
  %  ----------------------------------------------------------------------
  
  % get all input scans/surfaces
  cscc.files.data = {}; cscc.files.dataid = zeros(0,3); 
  for i = 1:numel(job.(datafield))
    cscc.files.data   = [cscc.files.data;cellstr(char(job.(datafield){i}))];
    cscc.files.dataid = [cscc.files.dataid; ...
     [(numel(cscc.files.dataid) + (1:numel(job.(datafield){i})) - 1)', ...
      repmat(i,numel(job.(datafield){i}),1), ...
      (1:numel(job.(datafield){i}))']];
  end
  
  
  
  %  define trash directory 
  %  ----------------------------------------------------------------------
  trashdirname = '+cat_checkcov_excluded'; 
  if ~isfield(job,'trashdirhome')
    [dirnames,dparts] = spm_str_manip(char(cscc.files.data(:)),'hC'); 
    [pp,dd,ee] = spm_fileparts(dparts.s); dd = [dd ee]; 
  else
    pp = job.trashdirhome;
  end
  cscc.trashdir = fullfile(pp,trashdirname);  
  try 
    if ~exist(cscc.trashdir,'dir'), mkdir(cscc.trashdir); end
  catch
    cat_io_cprintf('warn','Was not able to create default exclusion directory:\n  %s',cscc.trashdir); 
    cscc.trashdir = fullfile(pp,dd,trashdirname);
    try 
      if ~exist(cscc.trashdir,'dir'), mkdir(cscc.trashdir); end
    catch
      cat_io_cprintf('warn','Was not able to create alternative exclusion directory:\n  %s',cscc.trashdir); 
      cscc.trashdir = fullfile(spm_select(1,'dir','Select exclusion home directory'),trashdirname);
      try 
        if ~exist(cscc.trashdir,'dir'), mkdir(cscc.trashdir); end
      catch
        error('cat_stat_check_cov2:mktrashdir',...
          'Was not able to create exclusion directory (check writing rights):\n  %s',cscc.trashdir); 
      end
    end
  end  
  fprintf('Exclusion directory: \n  %s\n\n',cscc.trashdir)

  
  
  % remove files that do not exist
  cscc.oldtrash = cat_vol_findfiles(cscc.trashdir,'*-*-*',struct('depth',1)); 
  isfirst = 1; 
  for fi=numel(cscc.files.data):-1:1
    [pp,ff,ee] = spm_fileparts(cscc.files.data{fi}); 
    
    % remove MAC OS hidden files
    if strcmp(ff(1:2),'._')
      cscc.files.data(fi) = []; 
      job.(datafield){cscc.files.dataid(fi,2)}(cscc.files.dataid(fi,3)) = []; 
      if isfield(job,'data_xml') && ~isempty(job.data_xml) && numel(job.data_xml)>fi
        job.data_xml(fi) = []; 
      end
      if ~isempty(job.c)
        for i=1:numel(job.c)
          job.c{i}(fi) = [];
        end
      end
      continue
    end
    
    
    % remove non existing OS hidden files
    if ~exist(fullfile(pp,[ff ee]),'file')
      if isfirst
        fprintf('Miss input data:\n')
        isfirst=0;
      end
      fprintf('  %s\n',cscc.files.data{fi}); 
      cscc.files.data(fi) = []; 
      job.(datafield){cscc.files.dataid(fi,2)}(cscc.files.dataid(fi,3)) = []; 
      if ~isempty(job.data_xml) && numel(job.data_xml)>fi
        job.data_xml(fi) = []; 
      end
      if ~isempty(job.c)
        for i=1:numel(job.c)
          job.c{fi}(i) = [];
        end
      end
    end
  end
  
  
  
  % number of scans, samples and prtocols; trash mask arrays, sample array
  n_subjects  = numel(cscc.files.data);
  cscc.datagroups.n_samples   = numel(job.(datafield));
  cscc.datagroups.n_subjects  = numel(cscc.files.data);
  cscc.select.trash1d          = true(n_subjects,1);           % trash list mask 1D matrix
  cscc.select.trash2d         = true(n_subjects,n_subjects);  % trash list mask 2D matrix
  cscc.select.samp1d          = cscc.select.trash1d;
  cscc.select.samp2d          = cscc.select.trash2d;
  cscc.select.prot1d          = cscc.select.trash1d;
  cscc.select.prot2d          = cscc.select.trash2d;
  cscc.datagroups.sample      = [];
  cscc.datagroups.protocol    = []; 
  for i = 1:cscc.datagroups.n_samples
    cscc.datagroups.sample = ...
      [cscc.datagroups.sample, i * ones(1,size(job.(datafield){i},1))];
  end

  


  %% get the different files
  %  ----------------------------------------------------------------------
  spm_progress_bar('Init',numel(cscc.files.data),'Search files','subjects completed')
  spm_figure('GetWin','Interactive');
  [filenames,fparts] = spm_str_manip(cscc.files.data,'trC'); 
  cscc.files.out   = cscc.files.data;
  cscc.files.org   = cscc.files.data;
  cscc.files.pdf   = cscc.files.data;
  cscc.files.jpg   = cscc.files.data;
  cscc.files.xml   = cscc.files.data; 
  cscc.files.log   = cscc.files.data;
  cscc.files.surf  = cscc.files.data;
  cscc.files.surfr = cscc.files.data;
  
  % get real prefix
  % - expect that all files have the same prefix
  % - fparts.s is not enough if all file start similar, eg. mwp1ADNI_*.nii
  [tmp,lfname] = max(cellfun('length',fparts.m)); % longest filename
  [pp,ff,ee] = spm_fileparts(cscc.files.data{lfname});
  [pp1,pp2]  = spm_fileparts(pp); 
  orgfile    = cat_vol_findfiles(pp1,...
    ['*' cat_io_strrep(ff,fparts.s,'') '.nii'],struct('depth',1)); 
  if isempty(orgfile)
    orgfile  = cat_vol_findfiles(pp1,...
      ['*' cat_io_strrep(ff,fparts.s,'') '.img'],struct('depth',1)); 
  end
  [tmp,sfname]  = min(cellfun('length',orgfile)); % shortest of the longest
  [ppo,ffo,eeo] = spm_fileparts(orgfile{sfname});
  cscc.files.dataprefix = ff(1:strfind(ff,ffo)-1);
  clear orgfile;
  
  % find files
  for i = 1:numel(cscc.files.data)
    [pp,ff,ee] = spm_fileparts(cscc.files.data{i});
    [pp1,pp2]  = spm_fileparts(pp); 

    % output files
    cscc.files.out{i} = fullfile(pp,ff,ee); 

    % set subdirectories
    if strcmp(pp2,datadir)
      reportdir = 'report';
      surfdir   = 'surf';
    else
      reportdir = ''; 
      surfdir   = '';
    end

    fname = cat_io_strrep(ff,cscc.files.dataprefix,'');
    % set original input files of the CAT preprocessing
    cscc.files.org{i} = fullfile(pp1,[fname '.nii']); 
    if ~exist(cscc.files.org{i},'file')
      cscc.files.org{i} = fullfile(pp1,...
        [cat_io_strrep(ff,cscc.files.dataprefix,'') '.img']);
      if ~exist(cscc.files.org{i},'file')
         f1 = cat_vol_findfiles(pp1,[fname '.nii']);
         if ~isempty(f1) && exist(f1{1},'file')
           cscc.files.org{i} = f1{1}; 
         else %if ~exist(cscc.files.org{i},'file')
           cscc.files.org{i} = ''; 
         end
      end
    end
    
    cscc.files.p0{i} = fullfile(pp,['p0' fname '.nii']);
    if ~exist(cscc.files.p0{i},'file')
      cscc.files.p0{i} = ''; 
    end
    
    % try to find the cscc.data.XML file if not given
    if ~isfield(job,'data_xml') || isempty( char(job.data_xml) ) 
      cscc.files.xml{i} = fullfile(pp1,reportdir,['cat_' fname ...
        cat_io_strrep(ee,{'.nii','.img','.gii'},'.xml')]);
      if ~exist(cscc.files.xml{i},'file')
        cscc.files.xml{i} = ''; 
      end
    else
      cscc.files.xml = cellstr(job.data_xml);
    end

    % set report pdf
    cscc.files.pdf{i} = fullfile(pp1,reportdir,['catreport_' fname '.pdf']); 
    if ~exist(cscc.files.pdf{i},'file')
      cscc.files.pdf{i} = ''; 
    end

    % set report jpg
    cscc.files.jpg{i} = fullfile(pp1,reportdir,['catreportj_' fname '.jpg']); 
    if ~exist(cscc.files.jpg{i},'file')
      cscc.files.jpg{i} = ''; 
    end
    
    % log files
    cscc.files.log{i} = fullfile(pp1,reportdir,['catlog_' fname '.txt']);
    if ~exist(cscc.files.log{i},'file')
      cscc.files.log{i} = ''; 
    end

    % surface files
    cscc.files.surf{i} = fullfile(pp1,surfdir,...
      ['lh.thickness.' fname ]);
    if ~exist(cscc.files.surf{i},'file')
      cscc.files.surf{i} = ''; 
    end

    % resampled surface files
    if exist(fullfile(pp1,surfdir),'dir')
      try
        surf = cat_vol_findfiles(fullfile(pp1,surfdir),...
          ['s*.thickness.resampled.' fname '.gii'],struct);
      catch
        surf = {}; 
      end 
      if ~isempty(surf)
        cscc.files.surfr{i} = surf{1};
      else
        cscc.files.surfr{i} = '';
      end
    end
    spm_progress_bar('Set',i);  
  end
  spm_progress_bar('Clear');
  
  if all(cellfun('isempty',cscc.files.org))
    cat_io_cprintf('warn','Failed to find the original files!\n');
  end
  if all(cellfun('isempty',cscc.files.pdf))
    cat_io_cprintf('warn','Failed to find the report files!\n');
  end
  if all(cellfun('isempty',cscc.files.xml))
    cat_io_cprintf('warn','Failed to find the XML files!\n');
  end

  
  
  %% load header 
  %  ----------------------------------------------------------------------
  % Load data scan by scan because the surfaces V-structure already include
  % the texture data that slow down loading. Moreover, this is more save if 
  % some org-files are missing.
  spm_progress_bar('Init',numel(cscc.files.data),'Load data','subjects completed')
  spm_figure('GetWin','Interactive');
  for i=1:numel(cscc.files.data)
    cscc.data.V(i)  = spm_data_hdr_read(cscc.files.data{i});
    cscc.data.Vo(i) = spm_data_hdr_read(cscc.files.org{i});
    spm_progress_bar('Set',i);
  end
  spm_progress_bar('Clear');
    

  %% load cscc.data.XML data
  %  ----------------------------------------------------------------------
  if all(cellfun('isempty',cscc.files.xml))
    cscc.H.isxml        = 0;
    cscc.data.QM_names  = '';
  else
    cscc.H.isxml        = 1;

    if size(cscc.files.xml,1) ~= n_subjects
      error('XML-files must have the same number as cscc.datagroups.sample size');
    end
   
    cscc.data.QM = nan(n_subjects,4 + ...
      (cscc.job.expertgui) + 2*cscc.H.mesh_detected);
    cscc.data.QM_names = {...
      'Noise rating (NCR)';...
      'Bias rating (ICR)';...
      'Resoution rating (RES)';...
      'Weighted overall image quality rating (IQR)';...
      'Protocol-based IQR (PIQR)'; ...
      'Euler number';...
      'Size of topology defects (TDS)'};
    cscc.data.QM_names = cscc.data.QM_names(1:min( numel(cscc.data.QM_names) , size(cscc.data.QM,2) )); % remove Euler
    
    spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
    spm_figure('GetWin','Interactive');
    for i=1:n_subjects
      % get basename for xml- and data files
      [pth, xml_name]  = fileparts(deblank(cscc.files.xml{i}));
      [pth, data_name] = fileparts(cscc.data.V(i).fname);

      % remove leading 'cat_'
      xml_name = xml_name(5:end);

      % check for filenames
      if isempty(strfind(data_name,xml_name)) && ~isempty(xml_name)
        fprintf('Please check file names because of deviating subject names:\n %s vs. %s\n',...
          cscc.data.V(i).fname,cscc.files.xml{i});
      end

      xml = cat_io_xml(deblank(cscc.files.xml{i}));
      if isempty(cscc.files.org{i}) 
        cscc.files.org{i} = xml.filedata.fname;
      end
      
      cscc.data.catrev = xml.software.version_cat;
      if isfield(xml,'qualityratings')
        cscc.data.QM(i,1:4)  = [xml.qualityratings.NCR xml.qualityratings.ICR ...
          xml.qualityratings.res_RMS xml.qualityratings.IQR];
        RMS(i,1) = xml.qualityratings.res_RMS;
      elseif isfield(xml,'QAM') % also try to use old version
        cscc.data.QM(i,1:4)  = [xml.QAM.cscc.data.QM.NCR xml.QAM.cscc.data.QM.ICR ...
          xml.qualityratings.res_RMS xml.QAM.cscc.data.QM.res_RMS xml.QAM.cscc.data.QM.IQR];
        RMS(i,1) = xml.QAM.res_RMS;
      else
        RMS(i,1) = nan; 
      end
      if cscc.job.expertgui
        cscc.data.QM(i,5)  = nan;
      end
      if cscc.H.mesh_detected && isfield(xml.subjectmeasures,'EC_abs')
        cscc.data.QM(i,end-1:end) = ...
          [xml.subjectmeasures.EC_abs xml.subjectmeasures.defect_size];
      end
      spm_progress_bar('Set',i);  
    end
    spm_progress_bar('Clear');
    
    % detect cscc.datagroups.protocols by resolution
    pacc = 2; % larger values to detect many cscc.datagroups.protocols, small values to have less
    RMS(isnan(RMS)) = 11; % avoid cscc.H.multiple NaN center
    [cscc.datagroups.protocols,tmp,cscc.datagroups.protocol] = ...
      unique(round(RMS*10^pacc)/10^pacc); clear pid RMS
    cscc.datagroups.protocols(cscc.datagroups.protocols==11) = 21; 
    cscc.datagroups.protocols = cscc.datagroups.protocols/2; % average mm rather than rating

    % remove last two columns if EC_abs and defect_size are not defined
    if cscc.H.mesh_detected && all(all(isnan(cscc.data.QM(:,end-1:end))))
      cscc.data.QM = cscc.data.QM(:,end-1:end);
    end

    % added cscc.datagroups.protocol depending QA parameter
    
    [Pth,rth,sq,rths,rthsc,sqs] = cat_tst_qa_cleaner_intern(...
      cscc.data.QM(:,4),struct('site',{cscc.datagroups.protocol},'figure',0));
    cscc.data.QM_names = [cscc.data.QM_names;{'Protocol IQR difference (PIQR)'}];
    cscc.data.QM(:,5) = rth(:,1) - cscc.data.QM(:,4);
  
    % convert marks into rps rating
    mark2rps   = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
    markd2rpsd = @(mark) ( mark*10) + isnan(mark).*mark;
    cscc.data.QM(:,1:4)  = mark2rps(cscc.data.QM(:,1:4));
    if cat_get_defaults('exptops.expertgui')>1
      cscc.data.QM(:,5)    = markd2rpsd(cscc.data.QM(:,5));
    end
  end



  %% add constant to nuisance parameter
  %  ----------------------------------------------------------------------
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
  %  ----------------------------------------------------------------------
  if cscc.H.mesh_detected
    % load surface texture data
    Y = spm_data_read(cscc.data.V)';
    cscc.data.data_array_org = Y'; 
    
    % optional global scaling
    if isfield(job,'gSF')
      for i=1:numel(cscc.data.V)
        Y(:,2) = Y(:,2)*job.gSF(i);
      end
    end

    Y(isnan(Y)) = 0;

    % rescue unscaled data min/max
    cscc.data.cdata(1) = min(Y(:));
    cscc.data.cdata(2) = max(Y(:));
    Y = Y - repmat(mean(Y,2), [1 size(Y,2)]);

    % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
    if ~isempty(G) 
      [indinf,tmp] = find(isinf(G) | isnan(G));
      if ~isempty(indinf)
        fprintf('Nuisance parameter for %s is Inf or NaN.\n',V(indinf).fname);
        return
      end
      Ymean = repmat(mean(Y), [n_subjects 1]);
      Y = Y - G*(pinv(G)*Y) + Ymean;
    end

    cscc.data.data_array = Y';
    cscc.data.YpY = (Y*Y')/n_subjects;
    
    % calculate residual mean square of mean adjusted Y
    Y = Y - repmat(mean(Y,1), [n_subjects 1]);
    cscc.data.data_array_diff = Y';

    %MSE = sum(Y.*Y,2);
  else
    if length(cscc.data.V)>1 && any(any(diff(cat(1,cscc.data.V.dim),1,1),1))
      error('images don''t all have same dimensions')
    end
    if max(max(max(abs(diff(cat(3,cscc.data.V.mat),1,3))))) > 1e-8
      error('images don''t all have same orientation & voxel size')
    end

    % consider image aspect ratio
    cscc.pos.slice(4) = cscc.pos.slice(4) * cscc.data.V(1).dim(2) / cscc.data.V(1).dim(1);

    slices = 1:job.gap:cscc.data.V(1).dim(3);

    dimx = length(1:job.gap:cscc.data.V(1).dim(1));
    dimy = length(1:job.gap:cscc.data.V(1).dim(2));
    Y    = zeros(n_subjects, prod(dimx*dimy));
    cscc.data.YpY  = zeros(n_subjects);
    %MSE  = zeros(n_subjects,1);
    cscc.data.data_array = zeros([cscc.data.V(1).dim(1:2) n_subjects]);



    %-Start progress plot
    %-----------------------------------------------------------------------
    spm_progress_bar('Init',cscc.data.V(1).dim(3),'Check correlation','planes completed')
    spm_figure('GetWin','Interactive');

    for j=slices

      M  = spm_matrix([0 0 j 0 0 0 job.gap job.gap job.gap]);

      for i = 1:n_subjects
        cscc.data.img = spm_slice_vol(cscc.data.V(i),M,[dimx dimy],[1 0]);
        cscc.data.img(isnan(cscc.data.img)) = 0;
        Y(i,:) = cscc.data.img(:);
        if isfield(job,'gSF')
          Y(i,:) = Y(i,:)*job.gSF(i);
        end
      end

      % make sure data is zero mean
      Y = Y - repmat(mean(Y,2), [1 prod(dimx*dimy)]);

      % remove nuisance and add mean again 
      % (otherwise correlations are quite small and misleading)
      if ~isempty(G) 
        Ymean = repmat(mean(Y), [n_subjects 1]);
        Y = Y - G*(pinv(G)*Y) + Ymean;
      end

      cscc.data.YpY = cscc.data.YpY + (Y*Y')/n_subjects;

      % calculate residual mean square of mean adjusted Y
      Y = Y - repmat(mean(Y,1), [n_subjects 1]);

      %MSE = MSE + sum(Y.*Y,2);

      spm_progress_bar('Set',j);  

    end

    % correct filenames for 4D data
    if strcmp(cscc.data.V(1).fname, cscc.data.V(2).fname)
      cscc.data.Vchanged_names = 1;
      cscc.data.Vchanged      = cscc.data.V;
      for i=1:n_subjects
        [pth,nam,ext] = spm_fileparts(cscc.data.V(i).fname);
        cscc.data.V(i).fname    = fullfile(pth, [nam sprintf('%04d',i) ext]);
      end
    else
      cscc.data.Vchanged_names = 0;
    end

    spm_progress_bar('Clear');
  end
  clear Y



  %% normalize cscc.data.YpY and estimate cscc.data.mean_cov
  %  ----------------------------------------------------------------------
  d      = sqrt(diag(cscc.data.YpY)); % sqrt first to avoid under/overflow
  dd     = d*d';
  cscc.data.YpY    = cscc.data.YpY./(dd+eps);
  t      = find(abs(cscc.data.YpY) > 1); 
  cscc.data.YpY(t) = cscc.data.YpY(t)./abs(cscc.data.YpY(t));
  cscc.data.YpY(1:n_subjects+1:end) = sign(diag(cscc.data.YpY));
  clear t d dd;

  % extract mean correlation for each data set
  cscc.data.mean_cov = zeros(n_subjects,1);
  for i=1:n_subjects
    cov0        = cscc.data.YpY(i,:);     % extract row for each subject
    cov0(i)     = [];           % remove cov with its own
    cscc.data.mean_cov(i) = mean(cov0);
  end
  clear cov0;



  %% output compressed filenames structure
  %  ----------------------------------------------------------------------
  fprintf('\n');
  fname_m = [];
  fname_tmp = cell(cscc.datagroups.n_samples,1);
  fname_s   = cell(cscc.datagroups.n_samples,1);
  fname_e   = cell(cscc.datagroups.n_samples,1);
  for i=1:cscc.datagroups.n_samples
    [tmp, fname_tmp{i}] = ...
      spm_str_manip(char(cscc.data.V(cscc.datagroups.sample == i).fname),'C');
    fname_m = [fname_m; fname_tmp{i}.m];
    fname_s{i} = fname_tmp{i}.s;
    cat_io_cprintf('n','Compressed filenames sample %d: ',i);
    cat_io_cprintf('b',sprintf('%s %s \n',spm_str_manip(tmp,'f120'),...
      repmat('.',1,3*(numel(tmp)>120))));
  end
  cscc.files.fname = struct('s',{fname_s},'e',{fname_e},'m',{fname_m});
  clear fname_e fname_m fname_s fname_tmp tmp



  %% print suspecious files with high cov
  % use slightly higher threshold for (smoothed) mesh data
  %  ------------------------------------------------------------------------
  cscc.data.YpY_tmp = cscc.data.YpY - tril(cscc.data.YpY);
  if cscc.H.mesh_detected
    [indx, indy] = find(cscc.data.YpY_tmp>0.950 & cscc.data.YpY_tmp < (1-eps));
  else
    [indx, indy] = find(cscc.data.YpY_tmp>0.925);
  end
  [siv,si] = sort(cscc.data.YpY(sub2ind(size(cscc.data.YpY),indx,indy)),'descend');
  % if more than 25% of the data this points to longitudinal data of one 
  % subject and no warning will appear
  cscc.data.islongitudinal = (length(indx) > 0.25*n_subjects);
  if ~isempty(indx)
    if ~cscc.data.islongitudinal
      fprintf('\nUnusual large correlation (check that subjects are not identical):\n');
      for i=si'
        % exclude diagonal
        if indx(i) ~= indy(i)
          cat_io_cprintf('w',sprintf('  %0.4f',cscc.data.YpY(indx(i),indy(i)))); 
          cat_io_cprintf('n',' between ');
          cat_io_cprintf('b',cscc.files.fname.m{indx(i)}); cat_io_cprintf('n',' and ');
          cat_io_cprintf('b',cscc.files.fname.m{indy(i)}); fprintf('\n');
        end
      end
    else
      fprintf('\nMany unusual large correlations were found (e.g. common in longitudinal data).\n');
    end
  end

  [indx, indy] = find(cscc.data.YpY_tmp == 1);
  % give warning that data are identical
  if ~isempty(indx)
    fprintf('\nWARNING: Data of these subjects are identical!\n');
    for i=1:length(indx)
      if cscc.datagroups.n_samples > 1
        fprintf('%s (sample %d) and %s (sample %d)\n',filename.m{indx(i)},sample(indx(i)),filename.m{indy(i)},sample(indy(i)));
        cat_io_cprintf('w',sprintf('  %0.4f',cscc.data.YpY(indx(i),indy(i)))); 
        cat_io_cprintf('n',' between ');
        cat_io_cprintf('b',sprintf('%s (sample %d)',cscc.files.fname.m{indx(i)},cscc.datagroups.sample(indx(i)))); cat_io_cprintf('n',' and ');
        cat_io_cprintf('b',sprintf('%s (sample %d)',cscc.files.fname.m{indy(i)},cscc.datagroups.sample(indy(i)))); fprintf('\n');
      else
        cat_io_cprintf('w',sprintf('  %0.4f',cscc.data.YpY(indx(i),indy(i)))); 
        cat_io_cprintf('n',' between ');
        cat_io_cprintf('b',cscc.files.fname.m{indx(i)}); cat_io_cprintf('n',' and ');
        cat_io_cprintf('b',cscc.files.fname.m{indy(i)}); fprintf('\n');
      end
    end
  end



  %% sort data and estimate critical files
  %  ------------------------------------------------------------------------
  [cscc.data.mean_cov_cscc.H.sorted, cscc.data.ind_sorted] = sort(cscc.data.mean_cov,'descend');
  threshold_cov      = mean(cscc.data.mean_cov) - 2*std(cscc.data.mean_cov);
  n_thresholded      = find(cscc.data.mean_cov_cscc.H.sorted < threshold_cov,1,'first');
  if ~isempty(n_thresholded)
    fprintf('\nThese data have a mean correlation below 2 standard deviations. \n');
    fprintf('This does not necessarily mean that you have to exclude these data. \n');
    fprintf('However, these data have to be carefully checked:\n');
    for i=n_thresholded:n_subjects
      cat_io_cprintf('r',sprintf('  %0.4f ',cscc.data.mean_cov_cscc.H.sorted(i)));
      cat_io_cprintf('b',cscc.data.V(cscc.data.ind_sorted(i)).fname); fprintf('\n');
    end
  end



  %% output structure
  %  ------------------------------------------------------------------------
  if nargout>0
    varargout{1} = struct(...
      'table',[cscc.files.out,num2cell(cscc.data.mean_cov)],...
      'covmat',cscc.data.YpY,...
      'sorttable',[cellstr(cscc.data.V(cscc.data.ind_sorted).fname),...
        num2cell(cscc.data.mean_cov_cscc.H.sorted)],...
      'sortcovmat',cscc.data.YpY(cscc.data.ind_sorted,cscc.data.ind_sorted), ...
      'cov',cscc.data.mean_cov,...
      'threshold_cov',threshold_cov);
  end
  clear cscc.data.mean_cov_cscc.H.sorted threshold_cov;

  
  
  % check for replicates
  for i=1:n_subjects
    for j=1:n_subjects
      if (i>j) && (cscc.data.mean_cov(i) == cscc.data.mean_cov(j))
        try
          nami = deblank(cscc.data.V(i).fname);
          namj = deblank(cscc.data.V(j).fname);
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
  %  ----------------------------------------------------------------------
  if cscc.H.mesh_detected
    create_figures(job)
  else
    create_figures(job,slices)
  end
  show_boxplot;

 
return
%-End
%--------------------------------------------------------------------------
function create_figures(job,slices)
% -------------------------------------------------------------------------
% Create the main GUI figures of cat_stat_check_cov2.
% -------------------------------------------------------------------------

  global cscc 
  
  
  cscc.H.graphics = spm_figure('FindWin','Graphics');
  cscc.H.figure   = figure(2);
  set(cscc.H.figure,'MenuBar','none','Position',cscc.pos.fig,...
    'NumberTitle','off','Resize','off','Visible','on',...
    'HitTest','off','Interruptible','off',...
    'color',cscc.display.figcolor,'CloseRequestFcn',@closeWindows); %,'cat_stat_check_cov2(''closeWindows'')'); 


  % move SPM Graphics figure to have no overlap 
  fcscc.pos = get(cscc.H.figure,'Position');
  cscc.display.WP = get(cscc.H.graphics,'Position');
  set(cscc.H.graphics,'Position',[sum(fcscc.pos(1:2:3)) + 20 cscc.display.WP(2:4)]);
  clear fcscc.pos cscc.display.WP;

  if cscc.H.mesh_detected 
    set(cscc.H.figure,'Name','CAT Check Covariance: Click in image to display surfaces');
  else
    set(cscc.H.figure,'Name','CAT Check Covariance: Click in image to display slices');
  end

  
  % cursormode and update function
  cscc.H.dcm = datacursormode(cscc.H.figure);
  set(cscc.H.dcm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');

  
  % create two-area colormap
  colormap([jet(64); gray(64)]);

  
  % add colorbar without ticks and colorbar image
  cscc.H.cbar = axes('Position',cscc.pos.cbar,'Parent',...
    cscc.H.figure,'Visible','off');
  image(cscc.H.cbar,1:64); set(get(cscc.H.cbar,'children'),...
    'HitTest','off','Interruptible','off');
  set(cscc.H.cbar,'Ytick','','YTickLabel',''); 

  
  % set correlation matrix image as image
  cscc.H.corr = axes('Position',cscc.pos.corr,'Parent',cscc.H.figure,...
    'Color',cscc.display.figcolor,'Ytick','','YTickLabel','',...
    'XTickLabel','','XTick','','XTickLabel','');
  image(cscc.H.corr,64 * tril(cscc.data.YpY)); axis image;

  
  % scatter plot 
  if cscc.H.isxml
    cscc.H.scat(1) = axes('Position',cscc.pos.corr,'Parent',cscc.H.figure, ...
      'visible','off','Box','on','Color',[0.85 0.85 0.85]);

    if cscc.job.expertgui
      cscc.H.scat(2) = axes('Position',cscc.pos.corr,'Parent',cscc.H.figure, ...
      'visible','off','Box','on','Color',[0.85 0.85 0.85]);
    end
  end
  
  
  % slice axis
  cscc.H.slice = axes('Position',cscc.pos.slice,'Parent',cscc.H.figure, ...
    'visible','off','Box','off','Color',cscc.display.figcolor,...
    'Ytick','','YTickLabel','','XTickLabel','',...
    'XTick','','XTickLabel','');
  if ~cscc.H.mesh_detected, axis off; end
  
  
  % add button for closing all windows
  cscc.H.close = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.close,'Style','Pushbutton',...
    'callback',@closeWindows,'string','Close','ToolTipString','Close windows',...
    'FontSize',cscc.display.FS(cscc.display.FSi),'ForegroundColor',[0.8 0 0],...
    'Visible','off');
  
  
  % check worst scans >> maybe on the SPM graphics figure?!  
  % this would allow to show the worst data specific for the boxplot data!
  cscc.H.worst = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.close,'Style','Pushbutton',...
    'HorizontalAlignment','center','callback',@check_worst_data,...
    'string','Check worst','ToolTipString','Display most deviating files (MNC)',...
    'FontSize',cscc.display.FS(cscc.display.FSi),'ForegroundColor',[0.8 0 0]);

  
  % add button to open the HTML help
  cscc.H.help = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.help,'Style','Pushbutton',...
    'string','Help','ToolTipString','Open help window',...
    'ForegroundColor',[0 0 0.8],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'callback',['global cscc; spm_help(''!Disp'',fullfile(spm(''Dir''),''toolbox'...
    ',''cat12'',''html'',''cat_tools_checkcov.html'','''',cscc.H.graphics);']);

    
    
    

  %% create popoup menu for SPM grafix window
  
  % MNC
  str  = { 'Boxplot','Mean correlation'};
  tmp  = { {@show_boxplot, cscc.data.mean_cov, 'Mean correlation (MNC) ', 1} };
  
  if cscc.H.isxml
    
    % manual control of output
    showEQM      = cscc.job.expertgui; % other QM (noise,bias,res)
    showPIQR     = cscc.job.expertgui; % PIQR and MD2 plot
    
    % add QM header and IQR
    if showEQM || ( size(cscc.data.QM,2)>6 && cscc.H.mesh_detected )
      str = [str,{'Quality measures:'}]; 
      tmp = [ tmp , { { @show_boxplot, cscc.data.QM(:,4), cscc.data.QM_names{4}, 1 } } ]; 
    
      % IQR
      str = [ str , { '    Overall image quality rating (IQR)'  } ];
      tmp = [ tmp , { { @show_boxplot, cscc.data.QM(:,4), cscc.data.QM_names{4}, 1 } } ]; 
    else
      % only IQR
      str = [ str , { 'Overall image quality rating'  } ];
      tmp = [ tmp , { { @show_boxplot, cscc.data.QM(:,4), cscc.data.QM_names(4,:), 1 } } ]; 
    end  
    
    % add other QM
    % PIQR
    if showPIQR
      str = [ str , { '    Protocol-based image quality rating (PIQR)'  } ];
      tmp = [ tmp , { { @show_boxplot, cscc.data.QM(:,5), cscc.data.QM_names{5}, 1 } } ]; 
    end
    % NCR, ICR, RES
    if showEQM   
      for qmi = 1:3
        str = [ str , { ['    ' cscc.data.QM_names{qmi}] } ];
        tmp = [ tmp , { {@show_boxplot, cscc.data.QM(:,qmi), cscc.data.QM_names{qmi}, 1} }]; 
      end
    end
    % surfaces QM measures
    if size(cscc.data.QM,2)>6 && cscc.H.mesh_detected
      for qmi = 6:7
        str = [ str , { ['    ' cscc.data.QM_names{qmi}] } ];
        tmp = [ tmp , { {@show_boxplot, cscc.data.QM(:,qmi), cscc.data.QM_names{qmi}, 2} }]; 
      end
    end
    
    
    % estimate Mahalanobis distance between mean corr. and weighted overall quality
    cscc.data.X = [cscc.data.mean_cov, cscc.data.QM(:,4)];
    cscc.data.X(isnan(cscc.data.X)) = 0;  % mean correlation and IQR
    S   = cov(cscc.data.X);
    mu  = mean(cscc.data.X);
    Xmu = cscc.data.X - repmat(mu,[length(cscc.data.X),1]); 
    cscc.data.MD = diag(Xmu * inv(S) * Xmu');

    % if PIQR is used than we need further variables (similar to IQR)
    if showPIQR
      cscc.data.X2  = [cscc.data.mean_cov, cscc.data.QM(:,5)]; 
      cscc.data.X2(isnan(cscc.data.X2)) = 0; 
      S2   = cov(cscc.data.X2);
      mu2  = mean(cscc.data.X2);
      Xmu2 = cscc.data.X2 - repmat(mu2,[length(cscc.data.X2),1]); 
      cscc.data.MD2 = diag(Xmu2 * inv(S2) * Xmu2');
    end  

    % add PIQR
    if showPIQR 
      str = [str,{'Mahalanobis distances:'}]; %tmp = [ tmp , {{@sprintf,''}}]; 
      tmp = [ tmp , { {@show_boxplot, cscc.data.MD , 'Mahalanobis distance with IQR   ', 2 } } ]; 
      str = [ str , { '    with IQR' } ];
      tmp = [ tmp , { {@show_boxplot, cscc.data.MD , 'Mahalanobis distance with IQR   ', 2 } } ];  
      str = [ str , { '    with PIQR' } ];
      tmp = [ tmp , { {@show_boxplot, cscc.data.MD2, 'Mahalanobis distance with PIQR  ', 2 } } ];  
    else
      str = [ str , { 'Mahalanobis distance' } ];
      tmp = [ tmp , { {@show_boxplot, cscc.data.MD, 'Mahalanobis distance  ', 2 } } ];  
    end
    
  end
  
  
  % add nuisance variable(s) 
  % maybe as separate button?
  if isfield(job,'c') & ~isempty(job.c); 
    if numel(job.c)>1
      str = [str,{'Nuisance variables:'}]; %tmp = [ tmp , {{@sprintf,''}}]; 
    else
      str = [str,{'Nuisance variable'}];
    end
    tmp = [ tmp , {{@show_boxplot, job.c{1} ,sprintf('Nuisance variable %d  ',1), 0}}]; 
    if numel(job.c)>1
      for ci = 1:numel(1)
        str = [str,{sprintf('    Variable %d',ci)}]; 
        tmp = [ tmp , {{@show_boxplot, job.c{ci} ,sprintf('Nuisance variable %d  ',ci), 0}}]; 
      end
    end
  end
  
  
  % final create
  cscc.H.boxp = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.boxp,'Style','PopUp',...
    'callback','spm(''PopUpCB'',gcbo)','string',str,'UserData',tmp,...
    'ToolTipString','Display boxplot','FontSize',cscc.display.FS(cscc.display.FSi));


  
  
  
  
  %% create popoup menu for main check_cov window
  if cscc.H.isxml
    str  = { 'Plot',...
             'Corr. matrix order by selected files', ...
             'Corr. matrix sorted by mean corr.', ...
             'Mahalanobis distance'};
    tmp  = { {@show_matrix, cscc.data.YpY, 0},...
             {@show_matrix, cscc.data.YpY(cscc.data.ind_sorted,cscc.data.ind_sorted), 1},...
             {@show_mahalanobis, cscc.data.X, cscc.data.MD, 1}};
    if showPIQR
      str{4} = [ str{4} ' (IQR)'];
      str = char([cellstr(str),{'Mahalanobis distance PIQR'}]);
      tmp = [ tmp , ...
             {{@show_mahalanobis, cscc.data.X2, cscc.data.MD2, 2}}]; 
    end
  else
    str  = { 'Plot',...
             'Order by data selection',...
             'Sorted by mean correlation'};
    tmp  = { {@show_matrix, cscc.data.YpY, 0},...
             {@show_matrix, cscc.data.YpY(cscc.data.ind_sorted,cscc.data.ind_sorted), 1} };
  end

  cscc.H.sort = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.sort,'Style','PopUp',...
    'UserData',tmp,'callback','spm(''PopUpCB'',gcbo)','string',str,...
    'ToolTipString','Sort matrix','FontSize',cscc.display.FS(cscc.display.FSi));

  onoff = {'on','off'};
 
  
  
  
  %% == navigation unit ==
  cscc.H.selectuitext = uicontrol(cscc.H.figure,...
    'Units','normalized','Style','text','BackgroundColor',cscc.display.figcolor,...
    'Position',[cscc.pos.samp(1) cscc.pos.samp(2)+0.055 0.2 0.02],...
    'String','Data selection','FontSize',cscc.display.FS(cscc.display.FSi));

  % choose only one cscc.datagroups.sample for display
  str  = { 'Sample',sprintf('full (%d)',numel(cscc.datagroups.sample))}; 
  tmp  = { {@show_sample, 0, 1} }; 
  for i=1:cscc.datagroups.n_samples, 
    str = [str,sprintf('S%d (%d)',i,sum(cscc.datagroups.sample==i))]; 
    tmp = [ tmp , {{@show_sample, i, 1}} ]; 
  end
  cscc.H.samp = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.samp,'Style','PopUp',...
    'UserData',tmp,'enable',onoff{1 + (cscc.datagroups.n_samples==1)},...
    'callback','spm(''PopUpCB'',gcbo)','string',str,...
    'ToolTipString','Select sample to display','FontSize',cscc.display.FS(cscc.display.FSi));

  % choose center 
  if ~isfield(cscc.datagroups,'protocols')
    cscc.datagroups.protocols   = 1:max(cscc.datagroups.protocol); % only protocol ids
  end 
  str  = { 'Protocol',sprintf('all (%d)',numel(cscc.datagroups.protocol))};
  tmp  = { {@show_protocol, 0} }; 
  for i=1:numel(cscc.datagroups.protocols), 
    %str = [str,sprintf('P%03d',cscc.datagroups.protocols(i))]; % only protocol ids
    str = [str,sprintf('%5.2f (%d)',cscc.datagroups.protocols(i),...
      sum(cscc.datagroups.protocol==i))]; 
    tmp = [ tmp , {{@show_protocol, i}} ]; 
  end
  cscc.H.prot = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.prot,'Style','PopUp',...
    'UserData',tmp,'enable',onoff{ min( 2 , 1 + (numel(cscc.datagroups.protocols)==1) + ...
       1 - cscc.job.expertgui ) },...
    'callback','spm(''PopUpCB'',gcbo)','string',str,...
    'ToolTipString','Select protocol by its RMS resolution','FontSize',cscc.display.FS(cscc.display.FSi));
  
  
  %%
  cscc.H.alphabox = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.alphabox,'Style','CheckBox',...
    'callback',@update_alpha,...
    'string','Colorize diff. to sample mean','Value',1,...
    'ToolTipString','Colorize difference to sample mean (pos=green;neg=red)',...
    'Visible','off','BackgroundColor',cscc.display.figcolor,...
    'FontSize',cscc.display.FS(cscc.display.FSi));

  %{ 
  cscc.H.cbarfix = uicontrol(cscc.H.figure,...
    'Units','normalized','Style','CheckBox','position',cscc.pos.cbarfix,...
    'callback',{@checkbox_cbarfix},...
    'string','Fixed range','ToolTipString','Switch between fixed and auto-scaled colorbar',...
    'Value',0,'BackgroundColor',cscc.display.figcolor,...
    'FontSize',cscc.display.FS(cscc.display.FSi));
  %}
  
  cscc.H.showtrash = uicontrol(cscc.H.figure,...
    'Units','normalized','Style','CheckBox','position',cscc.pos.showtrash,...
    'callback',{@checkbox_showtrash},'enable','off',...
    'string','Show excluded','ToolTipString','Show excluded records',...
    'Value',0,'BackgroundColor',cscc.display.figcolor,'FontSize',...
    cscc.display.FS(cscc.display.FSi));

  
  
  
  %% add slider only for volume data
  if ~cscc.H.mesh_detected
    % voxelsize and origin
    vx   = sqrt(sum(cscc.data.V(1).mat(1:3,1:3).^2));
    Orig = cscc.data.V(1).mat\[0 0 0 1]';

    cscc.H.mm = uicontrol(cscc.H.figure,...
      'Units','normalized','position',cscc.pos.sslider,...
      ...'Min',(1 - Orig(3))*vx(3) ,'Max',(cscc.data.V(1).dim(3) - Orig(3))*vx(3),...
      'Min', -sum(slices<Orig(3)) * job.gap * vx(3),...
      'Max',  sum(slices>Orig(3)) * job.gap * vx(3),...
      'Style','slider','HorizontalAlignment','center',...
      'callback',@update_slices_array,...
      'ToolTipString','Select slice for display',...
      'SliderStep',[1 job.gap] / (cscc.data.V(1).dim(3)-1),'Visible','off');

    cscc.H.mm_txt = uicontrol(cscc.H.figure,...
      'Units','normalized','HorizontalAlignment','center',...
      'Style','text','BackgroundColor',cscc.display.figcolor,...
      'Position',[cscc.pos.sslider(1) cscc.pos.sslider(2)-0.005 0.2 0.02],...
      'String','0 mm','Visible','off','FontSize',cscc.display.FS(cscc.display.FSi));

    update_slices_array;
  end
  
  
  
  
  %% == navigation unit ==
  cscc.H.naviuitext = uicontrol(cscc.H.figure,...
    'Units','normalized','Style','text','BackgroundColor',cscc.display.figcolor,...
    'Position',[cscc.pos.scSelect(1) cscc.pos.scSelect(2)+0.035 0.2 0.02],...
    'String','Navigation options','FontSize',cscc.display.FS(cscc.display.FSi));

  cscc.H.naviui.select = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.scSelect,'callback','datacursormode(''on'')',...
    'Style','Pushbutton','enable','on','ToolTipString','Data selection');

  cscc.H.naviui.zoomReset = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.scZoomReset,...
    'callback','zoom out; datacursormode(''on'')',...
    'Style','Pushbutton','enable','on','ToolTipString','Reset zoom'); 

  cscc.H.naviui.zoomIn = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.scZoomIn,'callback',...
    ['global cscc; ' ...
     'if cscc.H.isscatter, ' ...
     '   hz = zoom(cscc.H.scat(cscc.H.scata)); ' ...
     ' else ' ...
     '   hz = zoom(cscc.H.corr);' ...
     ' end;' ...
     'set(hz,''enable'',''on'',''direction'',''in'')'], ... 
    'Style','Pushbutton','enable','on','ToolTipString','Zoom in');

  cscc.H.naviui.zoomOut = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.scZoomOut,'callback',...
    ['global cscc; ' ...
     'if cscc.H.isscatter, ' ...
     '   hz = zoom(cscc.H.scat(cscc.H.scata)); ' ...
     ' else ' ...
     '   hz = zoom(cscc.H.corr);' ...
     ' end;' ...
     'set(hz,''enable'',''on'',''direction'',''out'')'], ...
    'Style','Pushbutton','enable','on','ToolTipString','Zoom out');

  cscc.H.naviui.pan = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.scPan,'Enable','off','callback','pan on',...
    'Style','Pushbutton','enable','on','ToolTipString','Hand');


  

  %% == check unit ==
  cscc.H.checkuitext = uicontrol(cscc.H.figure,...
    'Units','normalized','HorizontalAlignment','center','Style','text',...
    'BackgroundColor',cscc.display.figcolor,...
    'Position',[cscc.pos.checkvol(1) cscc.pos.checkvol(2)+0.035 0.2 0.02],...
    'String','View selected data','FontSize',cscc.display.FS(cscc.display.FSi));

  % add button to open one image with SPM check_reg
  cscc.H.checkui.vol = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.checkvol,'callback',@checkvol,...
    'string','VOL','ToolTipString','Display original volume(s) in SPM graphics window',...
    'Style','Pushbutton','FontSize',cscc.display.FS(cscc.display.FSi),'Enable','off');  

  % add button to open one image with SPM check_reg
  cscc.H.checkui.surf = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.checksurf,'callback',@checksurf,...
    'string','SURF','ToolTipString','Display surface(s) in own figures',...
    'Style','Pushbutton','FontSize',cscc.display.FS(cscc.display.FSi),'Enable','off');  

  % add button to open one image with SPM check_reg
  cscc.H.checkui.log = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.checklog,'callback',@checklog,...
    'string','LOG','ToolTipString','Display log-file in SPM Graphics',...
    'Style','Pushbutton','FontSize',cscc.display.FS(cscc.display.FSi),'Enable','off');  

  % add button to open one image with SPM check_reg
  cscc.H.checkui.xml = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.checkxml,'callback',@checkxml,...
    'string','XML','FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Display xml-file in SPM Graphics',...
    'Style','Pushbutton','Enable','off'); 

  % add button to open the pdf in an external viewer 
  cscc.H.checkui.pdf = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.checkpdf,'callback',@checkpdf,...
    'ToolTipString','Display PDF report in external viewer',...
    'Style','Pushbutton','Enable','off');


  
  
  %% == trashlist unit ==
  cscc.H.trashuitext = uicontrol(cscc.H.figure,...
    'Units','normalized','Style','text',...
    'Position',[cscc.pos.newtrash(1) cscc.pos.newtrash(2)+0.035 0.2 0.02],...
    'BackgroundColor',cscc.display.figcolor,'String','Exclusion operations',...
    'FontSize',cscc.display.FS(cscc.display.FSi));

  % add button for new garbage mask
  cscc.H.trashui.new = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.newtrash,'callback',@newtrash,...
    'string','NEW','ForegroundColor',[ 0 0 0.8],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Reset exclusion list','Style','Pushbutton','Enable','off');  

  % add button to set the active image as garbage
  cscc.H.trashui.trash = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.trash,'callback',@trash,...
    'string','TIP+','ForegroundColor',[0.8 0 0],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Exclude record','Style','Pushbutton','Enable','off');

  cscc.H.trashui.trashcol = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.trash,'callback',@trash,...
    'string','COL-','ForegroundColor',[0.8 0 0],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Exclude column','Style','Pushbutton','Enable','off'); 

   % add button to remove the active image from garbage
  cscc.H.trashui.detrash = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.detrash,'callback',@detrash,...
    'string','TIP-','ForegroundColor',[0 0.8 0],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Include record','Style','Pushbutton','Enable','off');
  
  cscc.H.trashui.detrashcol = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.detrash,'callback',@detrash,...
    'string','COL+','ForegroundColor',[0 0.8 0],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Include column','Style','Pushbutton','Enable','off'); 
  
  % add button for mask below threshold as garbage
  cscc.H.trashui.disptrash = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.disptrash,'callback',@disptrash,...
    'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Print exclusion list in command window',...
    'string','VIEW','Style','Pushbutton','Enable','off'); 

  if isfield(cscc.data,'QM')
		cscc.H.trashui.autotrash = uicontrol(cscc.H.figure,...
			'Units','normalized','position',cscc.pos.autotrash,'callback',@autotrash,...
			'string','AUTO','FontSize',cscc.display.FS(cscc.display.FSi),...
			'ToolTipString','Automatic exclusion','ForegroundColor',[0 0.8 0],...
			'Style','Pushbutton','Enable',onoff{(size(cscc.data.QM,2)<4) + 1}); 
  end

  % == second row ==
  cscc.H.trashui.undo = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.undo,'callback',@trashundo,...
    'Style','Pushbutton','Enable','off','ToolTipString','Undo last exclusion operation'); 

  cscc.H.trashui.redo = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.redo,'callback',@trashredo,...
    'Style','Pushbutton','Enable','off','ToolTipString','Redo last exclusion operation'); 

  cscc.H.trashui.trashrow = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.trashrow,'callback',@trashrow,...
    'string','ROW-','ForegroundColor',[0.8 0 0],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Exclude row','Style','Pushbutton','Enable','off'); 

  cscc.H.trashui.detrashrow = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.detrashrow,'callback',@detrashrow,...
    'string','ROW+','ForegroundColor',[0 0.8 0],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Include row','Style','Pushbutton','Enable','off'); 

  cscc.H.trashui.ziptrash = uicontrol(cscc.H.figure,...
    'Units','normalized','position',cscc.pos.ziptrash,'callback',@ziptrash,...
    'string','DEL','ForegroundColor',[0.8 0 0],'FontWeight','bold',...
    'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Move excluded records to exclusion directory',...
    'Style','Pushbutton','Enable','off');

  
  
  
  %% print real data
  show_matrix(cscc.data.YpY, cscc.H.sorted);

  set(cscc.H.figure,'Visible','on');  
  try % catch main figure closing
    buttonupdate

    % redraw buttons
    pause(0.2) % Wait for the figure construction complete.
    warning off;  %#ok<WNOFF>
    jFig = get(cscc.H.figure, 'JavaFrame'); % get JavaFrame. You might see some warnings.
    warning on; %#ok<WNON>
    jWindow = jFig.fHG2Client.getWindow; % before 2011a it could be `jFig.fFigureClient.getWindow`. 
    jbh = handle(jWindow,'CallbackProperties'); % Prevent memory leak
    set(jbh,'ComponentMovedCallback',{@(~,~)(buttonupdate)});
  
    % show plot
    show_boxplot(cscc.data.mean_cov,'Mean correlation  ',1);
    set(cscc.H.figure,'HitTest','on','Interruptible','on');
  end
  
return

function buttonupdate
%-----------------------------------------------------------------------
% Put icons on buttons of the main GUI. 
%-----------------------------------------------------------------------
  global cscc

  % Update figure icons
  % close
  % help
  
  
  %% == navi buttons ==
  buttonicon(cscc.H.naviui.select,'DC'  ,...
    fullfile(matlabroot,'toolbox','matlab','icons','tool_data_cursor.png'));
  buttonicon(cscc.H.naviui.zoomReset,'Zo'  ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','tool_fit.png'));
  buttonicon(cscc.H.naviui.zoomIn,'Z+'  ,...
    fullfile(matlabroot,'toolbox','matlab','icons','tool_zoom_in.png'));
  buttonicon(cscc.H.naviui.zoomOut,'Z-'  ,...
    fullfile(matlabroot,'toolbox','matlab','icons','tool_zoom_out.png'));
  buttonicon(cscc.H.naviui.pan,'H'   ,...
    fullfile(matlabroot,'toolbox','matlab','icons','tool_hand.png'));

  
  
  %% == check buttons ==
  buttonicon(cscc.H.checkui.vol,'VOL' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_spm_view.png'));
  buttonicon(cscc.H.checkui.surf,'SURF',...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_surfc.png'));
  buttonicon(cscc.H.checkui.log ,'LOG' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_cat_log.png'));
  buttonicon(cscc.H.checkui.xml,'XML' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_cat_xml.png'));
  buttonicon(cscc.H.checkui.pdf,'PDF' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_cat_report_surf.png'));

  
  
  %% == trash buttons ==
  buttonicon(cscc.H.trashui.new,'NEW'  ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_new.png'));
  buttonicon(cscc.H.trashui.disptrash,'VIEW' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_export.png'));
  buttonicon(cscc.H.trashui.undo,'UNDO' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','undo.png'));
  buttonicon(cscc.H.trashui.redo,'REDO' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','redo.png'));
  buttonicon(cscc.H.trashui.trash,'TIP-' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_tip_rm.png'));
  buttonicon(cscc.H.trashui.detrash,'TIP+' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_tip_add.png'));
  set([cscc.H.trashui.trash,cscc.H.trashui.detrash],'visible','off'); 
  buttonicon(cscc.H.trashui.trashcol,'COL-' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_col_rm.png'));
  buttonicon(cscc.H.trashui.detrashcol,'COL+' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_col_add.png'));
  buttonicon(cscc.H.trashui.trashrow,'ROW-' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_row_rm.png'));
  buttonicon(cscc.H.trashui.detrashrow,'ROW+' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_row_add.png'));
  buttonicon(cscc.H.trashui.autotrash,'AUTO' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_auto.png'));
  buttonicon(cscc.H.trashui.ziptrash,'DEL' ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','file_delete.png'));

return
%-----------------------------------------------------------------------

function buttonicon(h,str,Picon,useboth) 
%-----------------------------------------------------------------------
% Function to print an image file Picon on the button with handle h. Use
% the string str if the global variable cscc.display.useicons<1 or other 
% errors.
%-----------------------------------------------------------------------
  global cscc


  usethisicon = cscc.display.useicons;

  set(h,'Visible','on');  

  if ~exist(Picon,'file')
    warning('Button Icon "%s" does not exist!',Picon); 
  end

  if ~exist('useboth','var') 
    useboth = 0; 
  end

  if ~exist('cat_io_findjobj','file')
    warning('JAVA Function for  "%s" does not exist!',Picon);
    cscc.display.useicons = 0; 
  end

  set(h,'string','');

  if usethisicon
    try
      jButton = cat_io_findjobj(h);
      jButton.setIcon(javax.swing.ImageIcon(Picon));
      jButton.setHorizontalTextPosition(javax.swing.SwingConstants.RIGHT);
      jButton.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
    catch 
      usethisicon = 0;
    end
  end

  try
    if usethisicon<1 || useboth 
      set(h,'string',str,'FontSize',cscc.display.FS(cscc.display.FSi));
    else
      set(h,'string','');
    end
  end
return
%-----------------------------------------------------------------------

function estimateExclusion(varargin)
%-----------------------------------------------------------------------
% Estimate critical scans that should be excluded and list them on the
% global variable "at".  
%-----------------------------------------------------------------------

  global cscc at

  useIQR = isfield(cscc.data,'QM') & cscc.H.ETB.PIQRrslider.isEnabled; 
  
  % if a gui is available use it
  if isfield(cscc.H,'ETB') && nargin>0;
    %MNCath  = get(cscc.H.ETB.MNCaslider,'Value')/10;
    MNCrth  = get(cscc.H.ETB.MNCrslider,'Value')/10;
    %PIQRath = get(cscc.H.ETB.PIQRaslider,'Value');
    PIQRrth = get(cscc.H.ETB.PIQRrslider,'Value')/10;
  else
    MNCrth  = at.MNCrth;
    PIQRrth = at.PIQRrth;
  end
  MNCath  = at.MNCath*2;
  PIQRath = at.PIQRath*2;
  
  %  image quality 
  %  ----------------------------------------------------------------------
  datat    = false(size(cscc.data.QM,1),1); 
  datat(cscc.select.trashlist) = 1; 
  
  if useIQR
    delIQRa  = find( cscc.data.QM(:,5) < -at.PIQRath );
    mdPIQR   = cat_stat_nanmedian(cscc.data.QM(:,5));
    sdPIQR   = cat_stat_nanstd(cscc.data.QM(:,5));
    dataf    = (cscc.data.QM(:,5) < (mdPIQR - PIQRath) ) | ...        
               (cscc.data.QM(:,5) < (mdPIQR - PIQRrth * sdPIQR)); 
    dataf(datat) = 0; 
    delIQR   = unique( [ delIQRa ; find( dataf ) ] );
    dataf(delIQRa)=1;
    if isfield(cscc.H,'ETB');
      cscc.H.ETB.histPIQRp.Data = cscc.data.QM(~dataf & ~datat,5);
      cscc.H.ETB.histPIQRf.Data = cscc.data.QM( dataf,5);
      cscc.H.ETB.histPIQRt.Data = cscc.data.QM( datat,5);
    end
    datafiqr=dataf;  
  end
                   
  %  removal of data based on the mean covarance
  %  ----------------------------------------------------------------------
  %  negative outlier that does not fit well to this dataset
  %  - using relative weighting based on the median and the standard 
  %    deviation (eg. scans lower than 2*sdMNC from mdMNC) 
  %  - using a absolute value (eg. < mdMNC - 0.05 ) 
  mdMNC   = cat_stat_nanmedian(cscc.data.mean_cov);
  sdMNC   = cat_stat_nanstd(cscc.data.mean_cov);
  dataf   = (cscc.data.mean_cov < (mdMNC - MNCath) ) | ...         % records below absolute lower threshold
            (cscc.data.mean_cov < (mdMNC - MNCrth * sdMNC) );
  delMNCl = find( dataf );    % records below relative lower threshold
  
   
  % create final exclusion list for MNC 
  datad = false(size(dataf)); 
  if ~cscc.data.islongitudinal
    % positive outlier that fit well to another dataset (dublicate entries / rescans)
    delMNCh = find((tril(cscc.data.YpY,-1) > (mdMNC + MNCath) ) & ...      % records above absolute higher threshold
                   (tril(cscc.data.YpY,-1) > (mdMNC + MNCrth * sdMNC) ) );
     
    % get x ad y record
    [delMNCh(:,2),delMNCh(:,1)] = ind2sub(size(cscc.data.YpY),delMNCh);
    
    % remove both entries if they are in different groups
    groups     = any(cscc.datagroups.sample(delMNCh) ~= ...
                   repmat(cscc.datagroups.sample(delMNCh(:,1))',1,2),2); 
    multigroup = delMNCh([groups;groups]);
    delMNCh    = delMNCh(:,2); 
    
    datad([delMNCh; multigroup(:)]) = 1;  
    delMNC     = [delMNCl;delMNCh;multigroup(:)];
  else
    delMNC     = delMNCl; 
  end

  if isfield(cscc.H,'ETB')
    cscc.H.ETB.histMNCp.Data = cscc.data.mean_cov(~dataf & ~datad & ~datat);
    cscc.H.ETB.histMNCf.Data = cscc.data.mean_cov( dataf);
    cscc.H.ETB.histMNCd.Data = cscc.data.mean_cov( datad);
    cscc.H.ETB.histMNCt.Data = cscc.data.mean_cov( datat);
  end
  
  % final combination
  if useIQR
    del = unique([delIQR',delMNC']); 
  else
    del = unique(delMNC'); 
  end
  del = setdiff( del , cscc.select.trashlist);
  
  if isfield(cscc.H,'ETB') && isfield(cscc.H.ETB,'tab') && any(cscc.H.ETB.tab(:).isvalid)
    %%
    isdel{4}=false(size(dataf)); isdel{4}(unique([cscc.select.trashlist,del]))=1; isdel{3}=~isdel{4};
    for i=3:size(cscc.H.ETB.tab,1)
      for j=2:size(cscc.H.ETB.tab,2)
        if cscc.H.ETB.tab(1,j).isvalid
          switch cscc.H.ETB.tab(1,j).String(1)
            case 'R' 
              set(cscc.H.ETB.tab(i,j),'String',sprintf('%0.2f%% (%d)',...
                100 * (sum(isdel{i})) ./ numel(isdel{i}),sum(isdel{i})));
            case 'G'
              k = str2double(cscc.H.ETB.tab(1,j).String(2:end));
              set(cscc.H.ETB.tab(i,j),'String',sprintf('%d',sum(cscc.datagroups.sample'==k & isdel{i})));
            case 'N'
              k = str2double(cscc.H.ETB.tab(1,j).String(2:end));
              set(cscc.H.ETB.tab(i,j),'String',sprintf('%0.2f',mean(cscc.job.c{k}(isdel{i}))));
            case 'M'
              set(cscc.H.ETB.tab(i,j),'String',sprintf('%0.3f',mean(cscc.data.mean_cov(isdel{i}))));
            case 'I'
              set(cscc.H.ETB.tab(i,j),'String',sprintf('%0.2f',mean(cscc.data.QM(isdel{i},4))));
            case 'P'
              set(cscc.H.ETB.tab(i,j),'String',sprintf('%0.2f',mean(cscc.data.QM(isdel{i},5))));
            case 'E'
              set(cscc.H.ETB.tab(i,j),'String',sprintf('%0.2f',mean(cscc.data.QM(isdel{i},6))));
          end
        end
      end  
    end
  end
  if isfield(cscc.H,'ETB') && isfield(cscc.H.ETB,'tab') 
    %%
    legends = {'histMNClegend','histPIQRlegend'};
    for i=1:numel(legends)
      if isfield(cscc.H.ETB,legends{i}) && any(cscc.H.ETB.(legends{i})(:).isvalid)
        for j=1:numel(cscc.H.ETB.(legends{i}).String)
          switch cscc.H.ETB.(legends{i}).String{j}(1:2)
            case 'pa', 
              if i==1
                cscc.H.ETB.(legends{i}).String{j} = sprintf('passed (%d)',sum(~dataf & ~datad));
              else
                cscc.H.ETB.(legends{i}).String{j} = sprintf('passed (%d)',sum(~datafiqr));
              end
            case 'ex',
              if i==1
                cscc.H.ETB.(legends{i}).String{j} = sprintf('excluded (%d)',sum( dataf & ~datad & ~datat));
              else
                cscc.H.ETB.(legends{i}).String{j} = sprintf('excluded (%d)',sum( datafiqr & ~datat));
              end
            case 're', cscc.H.ETB.(legends{i}).String{j} = sprintf('rescans (%d)',sum(datad & ~datat));
            case 'al', cscc.H.ETB.(legends{i}).String{j} = sprintf('alr. exc. (%d)',sum(datat));
          end
        end
      end
    end
  end
  
  at.del = del; 
return
%-----------------------------------------------------------------------

function autotrashGUI
%-----------------------------------------------------------------------
% GUI to control the exlusion of scans by low mean correlation or low 
% image quality.
%-----------------------------------------------------------------------
    global cscc 
    
    useIQR = isfield(cscc.data,'QM') && 1; 
  
    
    %% move to begin if finished
cscc.pos.trashpopup = [cscc.H.figure.Position(1:2) 1.2*cscc.display.WS(3) 0.8*cscc.display.WS(3)];
cscc.pos.trashpopup(1) = cscc.H.figure.Position(1) + (cscc.H.figure.Position(3) - cscc.pos.trashpopup(3))/2;
cscc.pos.trashpopup(2) = cscc.H.figure.Position(2) + (cscc.H.figure.Position(4) - cscc.pos.trashpopup(4)) - 24;        
                      
    % use the dialog if finished
    if isfield( cscc.H ,'ETB' ) && isfield( cscc.H.ETB ,'main' ) && ...
      ~isempty( cscc.H.ETB.main ) && isvalid(cscc.H.ETB.main )
      fpos = get(cscc.H.ETB.main,'Position'); 
      delete(cscc.H.ETB.main);
    else
      fpos = cscc.pos.trashpopup; 
    end
      
    %cscc.H.ETB.main = dialog('Position',cscc.pos.trashpopup,'Name','Remove problematic data');
    cscc.H.ETB.main  = figure('Position',fpos,'Name','Remove problematic data','MenuBar','none'); 
    
     % == MNC histogram ==
    mncrvals = [5 2 100 max(std(cscc.data.mean_cov)*400,150*(median(cscc.data.mean_cov) - min(cscc.data.mean_cov)))];
    mncbvals = [1 0 7 4];
    %cscc.H.ETB.mncbvalues = 2.^( (mncbvals(1):1:mncbvals(3)) - (mncbvals(3) - 4) )/100; % automatic
    cscc.H.ETB.mncbvalues = [0.001 0.002 1/300 0.005 0.01 0.02 1/30];
    cscc.H.ETB.histMNC      = subplot('position',[0.05 0.5 0.425 0.450]);
    cscc.H.ETB.histMNCrange = 0 : cscc.H.ETB.mncbvalues(mncbvals(4)) : 1;  
      
    if 0 %useIQR
      cscc.H.ETB.histMNCcb = uicontrol(cscc.H.ETB.main ,...
        'Units','normalized','Style','CheckBox','position',[0.05 0.96 0.03 0.03],'callback','',...
        'string','Fixed range','ToolTipString','Use MNC Thresholding',...
        'Value',1,'FontSize',cscc.display.FS(cscc.display.FSi));
      cscc.H.ETB.histPIQRcb = uicontrol(cscc.H.ETB.main ,...
        'Units','normalized','Style','CheckBox','position',[0.52 0.96 0.03 0.03],'callback','',...
        'string','Fixed range','ToolTipString','Use PIQR Thresholding',...
        'Value',1,'FontSize',cscc.display.FS(cscc.display.FSi));
    end
   
    % plot histrogam
    hold on; 
    cscc.H.ETB.histMNCp = histogram(cscc.H.ETB.histMNC,...
      cscc.data.mean_cov(cscc.data.mean_cov(:) > 0),cscc.H.ETB.histMNCrange); 
    cscc.H.ETB.histMNCf = histogram(cscc.H.ETB.histMNC,...
      cscc.data.mean_cov(cscc.data.mean_cov(:) < 0),cscc.H.ETB.histMNCrange); 
    cscc.H.ETB.histMNCd = histogram(cscc.H.ETB.histMNC,...
      cscc.data.mean_cov(cscc.data.mean_cov(:) < 0),cscc.H.ETB.histMNCrange); 
    cscc.H.ETB.histMNCt = histogram(cscc.H.ETB.histMNC,...
      cscc.data.mean_cov(cscc.data.mean_cov(:) < 0),cscc.H.ETB.histMNCrange); 
    set(cscc.H.ETB.histMNCp,'FaceColor',[0   0.7 0  ],'EdgeColor',[0 0.7   0  ],'LineStyle','none','FaceAlpha',1); 
    set(cscc.H.ETB.histMNCf,'FaceColor',[1.0 0   0  ],'EdgeColor',[0.9 0   0  ],'LineStyle','none','FaceAlpha',1); 
    set(cscc.H.ETB.histMNCd,'FaceColor',[0.0 0.5 1.0],'EdgeColor',[1.0 0.4 0.2],'LineStyle','none','FaceAlpha',1); 
    set(cscc.H.ETB.histMNCt,'FaceColor',[0.5 0   0  ],'EdgeColor',[0.5 0   0  ],'LineStyle','none','FaceAlpha',1); 
    cscc.H.ETB.histMNClegend = legend(...
      [cscc.H.ETB.histMNCp,cscc.H.ETB.histMNCf,cscc.H.ETB.histMNCd,cscc.H.ETB.histMNCt],...
      {'passed','excluded','rescans','alr. exc.'});
    hold off;

    % add title and set limits
    title('mean correlation (MNC)'); box on; grid on
    ylim( cscc.H.ETB.histMNC , ylim .* [1 1.2]); 
    xlim( cscc.H.ETB.histMNC , max(0,min(1,median(cscc.data.mean_cov)*1.01 + ...
      mncrvals(4)/10 * [-1 1/3] * median(cscc.data.mean_cov)/10 )));
    set(gca,'Xdir','reverse'); 
  
    
    
  
    % IQR histogram
    % print inactive looking elements if not available
    if useIQR, Color = zeros(1,3); else Color = ones(1,3) * 0.6; end
    if useIQR
      piqrvals = [5 1 10*max([20;abs(cscc.data.QM(:,5))])^(1/1.2) 10*max([3;1.8*abs(cscc.data.QM(:,5))])^(1/1.2)];
      piqbvals = [1 0 7 4];
      cscc.H.ETB.piqrvals = piqrvals;
      %cscc.H.ETB.piqbvalues = 2.^( (piqbvals(1):1:piqbvals(3)) - (piqbvals(3) - 2) );
      cscc.H.ETB.piqbvalues = [0.01 0.02 1/30 0.05 0.1 0.2 1/3 0.5];

      cscc.H.ETB.histPIQR = subplot('position',[0.525 0.5 0.425 0.450]);
      cscc.H.ETB.histPIQRrange = ...
        [ fliplr( 0 : -cscc.H.ETB.piqbvalues(piqbvals(4)) : min(-8,floor(min(cscc.data.QM(:,5))) ))  ...
          cscc.H.ETB.piqbvalues(piqbvals(4)) : cscc.H.ETB.piqbvalues(piqbvals(4)) : max(32,ceil(max(cscc.data.QM(:,5)))) ];  
      
      hold on
      cscc.H.ETB.histPIQRp = histogram(cscc.H.ETB.histPIQR,...
        cscc.data.QM(cscc.data.QM(:,5) < 3,5),cscc.H.ETB.histPIQRrange); 
      cscc.H.ETB.histPIQRf = histogram(cscc.H.ETB.histPIQR,...
        cscc.data.QM(cscc.data.QM(:,5) > 3,5),cscc.H.ETB.histPIQRrange); 
      cscc.H.ETB.histPIQRt = histogram(cscc.H.ETB.histPIQR,...
        cscc.data.QM(cscc.data.QM(:,5) > 3,5),cscc.H.ETB.histPIQRrange); 
      hold off
      set(cscc.H.ETB.histPIQRp,'FaceColor',[0   0.7 0],'EdgeColor',[0 0.4 0],'LineStyle','none','FaceAlpha',1); 
      set(cscc.H.ETB.histPIQRf,'FaceColor',[1   0   0],'EdgeColor',[0.4 0 0],'LineStyle','none','FaceAlpha',1); 
      set(cscc.H.ETB.histPIQRt,'FaceColor',[0.5 0   0],'EdgeColor',[0.5 0 0],'LineStyle','none','FaceAlpha',1); 
    
      title('protocol-base image quality rating (PIQR)'); box on; grid on
      ylim( ylim .* [1 1.2]); 
      xlim( cscc.H.ETB.histPIQR , max(min([-piqrvals(3)*2;cscc.data.QM(:,5)]),...
         min(max([piqrvals(3);cscc.data.QM(:,5)]),piqrvals(3)/100)^1.2 * [-1 1/2] ) )
      set(gca,'Xdir','reverse')
      cscc.H.ETB.histPIQRlegend = legend(...
        [cscc.H.ETB.histPIQRp,cscc.H.ETB.histPIQRf,cscc.H.ETB.histPIQRt],...
        {'passed','excluded','alr. exc.'});
    else 
      cscc.H.ETB.histPIQR = subplot('position',[0.525 0.5 0.425 0.450],...
        'Color',ones(1,3) * 0.97,'XColor',Color,'YColor',Color);
      hold on
      plot(cscc.H.ETB.histPIQR,[0,1],[1,0],'Color',Color); 
      plot(cscc.H.ETB.histPIQR,[0,1],[0,1],'Color',Color); 
      hold off
      box on;
      title('protocol-base image quality rating (PIQR)','Color',Color)
    end
    
    
    
    
    % == View sliders ==
    % MNC view slider
    posm = [0.03 0.37 0.235 0.045]; 
    post = round((posm + [0 0.045 0 -0.015]) .* cscc.pos.trashpopup([3:4,3:4])); 
    poss = round(posm .* cscc.pos.trashpopup([3:4,3:4]));
    cscc.H.ETB.MNCvslider = javax.swing.JSlider; 
    javacomponent(cscc.H.ETB.MNCvslider,poss);
    uicontrol('Parent',cscc.H.ETB.main,'Style','text','position',post,'String','window aperture');
    set(cscc.H.ETB.MNCvslider,'Minimum',mncrvals(1),'Maximum',mncrvals(3),'Value',mncrvals(4),'snapToTicks',0,...
      'minorTickSpacing',mncrvals(2),'MajorTickSpacing',mncrvals(1),'PaintTicks',false, 'PaintLabels',false); 
    set(cscc.H.ETB.MNCvslider,...
      'StateChangedCallback',...
      ['global cscc;' ...
       'xlim( cscc.H.ETB.histMNC , max(0,min(1,median(cscc.data.mean_cov)*1.01+ '...
       '  (double(get(cscc.H.ETB.MNCvslider,''Value''))/10)^1.2 * [-1 1/3] * median(cscc.data.mean_cov)/10 ))); ' ...
       'set(cscc.H.ETB.histMNC,''Xdir'',''reverse'');']);

    % PIQR view slider
    posm = [0.51 0.37 0.235 0.045];
    post = round((posm + [0 0.045 0 -0.015]) .* cscc.pos.trashpopup([3:4,3:4])); 
    poss = round(posm .* cscc.pos.trashpopup([3:4,3:4]));
    uicontrol('Parent',cscc.H.ETB.main,'Style','text','position',post,'String','window aperture','ForegroundColor',Color);
    cscc.H.ETB.PIQRvslider = javax.swing.JSlider; 
    javacomponent(cscc.H.ETB.PIQRvslider,poss);
    if useIQR
      set(cscc.H.ETB.PIQRvslider,'Minimum',piqrvals(1),'Maximum',piqrvals(3),'Value',piqrvals(4),'snapToTicks',0,...
        'minorTickSpacing',piqrvals(2),'MajorTickSpacing',piqrvals(1),'PaintTicks',false, 'PaintLabels',false); 
      set(cscc.H.ETB.PIQRvslider,...
        'StateChangedCallback',...
        ['global cscc;' ...
         'xlim( cscc.H.ETB.histPIQR , max(min([-cscc.H.ETB.piqrvals(3)*2;cscc.data.QM(:,5)]),' ...
          ' min(max([cscc.H.ETB.piqrvals(3);cscc.data.QM(:,5)]), ' ...
          ' (double(get(cscc.H.ETB.PIQRvslider,''Value''))/10)^1.2 * [-1 1/2] ) ) ); '...
         'set(cscc.H.ETB.histPIQR,''Xdir'',''reverse'');']);
    else
      cscc.H.ETB.PIQRvslider.disable;
    end
    % change slider font size ???
    % fs = cscc.H.ETB.MNCvslider.getFont; 
    % fs.setSize = cscc.display.FS(cscc.display.FSi-1); 
    % cscc.H.ETB.MNCvslider.setFont(fs); % nothing to set the font-size
    % create new font? but how?
    
    

    % == Bin sliders ==
    % MNC bin slider
    posm = [0.26 0.37 0.235 0.045]; 
    post = round((posm + [0 0.045 0 -0.015]) .* cscc.pos.trashpopup([3:4,3:4])); 
    poss = round(posm .* cscc.pos.trashpopup([3:4,3:4]));
    cscc.H.ETB.MNCbslider = javax.swing.JSlider; 
    javacomponent(cscc.H.ETB.MNCbslider,poss);
    uicontrol('Parent',cscc.H.ETB.main,'Style','text','position',post,'String','bin width');
    set(cscc.H.ETB.MNCbslider,'Minimum',mncbvals(1),'Maximum',mncbvals(3),'Value',mncbvals(4),'snapToTicks',1,...
      'minorTickSpacing',mncbvals(2),'MajorTickSpacing',mncbvals(1),'PaintTicks',0, 'PaintLabels',false); 
    set(cscc.H.ETB.MNCbslider,'StateChangedCallback',['global cscc; ' ...
      'val = cscc.H.ETB.mncbvalues(get(cscc.H.ETB.MNCbslider,''Value'')); ' ...
      'cscc.H.ETB.histMNCp.BinWidth = val;' ...
      'cscc.H.ETB.histMNCf.BinWidth = val;' ...
      'cscc.H.ETB.histMNCt.BinWidth = val;' ...
      'cscc.H.ETB.histMNCd.BinWidth = val;' ...
      'ylim( cscc.H.ETB.histMNC , ''auto'' ); ylim( cscc.H.ETB.histMNC , ylim( cscc.H.ETB.histMNC ) .* [1 1.2]); ']);

    % PIQR bin slider
    posm = [0.74 0.37 0.235 0.045];
    post = round((posm + [0 0.045 0 -0.015]) .* cscc.pos.trashpopup([3:4,3:4])); 
    poss = round(posm .* cscc.pos.trashpopup([3:4,3:4]));
    uicontrol('Parent',cscc.H.ETB.main,'Style','text','position',post,'String','bin width','ForegroundColor',Color);
    cscc.H.ETB.PIQRbslider = javax.swing.JSlider; 
    javacomponent(cscc.H.ETB.PIQRbslider,poss);
    if useIQR
      set(cscc.H.ETB.PIQRbslider,'Minimum',piqbvals(1),'Maximum',piqbvals(3),'Value',piqbvals(4),'snapToTicks',1,...
        'minorTickSpacing',piqbvals(2),'MajorTickSpacing',piqbvals(1),'PaintTicks',0, 'PaintLabels',false); 
      set(cscc.H.ETB.PIQRbslider,'StateChangedCallback',['global cscc;' ...
        'val = cscc.H.ETB.piqbvalues(get(cscc.H.ETB.PIQRbslider,''Value'')); ' ...
        'cscc.H.ETB.histPIQRp.BinWidth = val;' ...
        'cscc.H.ETB.histPIQRf.BinWidth = val;' ...
        'cscc.H.ETB.histPIQRt.BinWidth = val;' ...
        'ylim( cscc.H.ETB.histPIQR , ''auto'' ); ylim( cscc.H.ETB.histPIQR , ylim( cscc.H.ETB.histPIQR ) .* [1 1.2]); ']);
    else
      cscc.H.ETB.PIQRvslider.disable;
    end
    
    
    
    % == Exclusion sliders ==
  
    % MNC exclusion slider
    posm = [0.03 0.26 0.465 0.06]; %0.465
    post = round((posm + [0 0.068 0 -0.03]) .* cscc.pos.trashpopup([3:4,3:4])); 
    poss = round(posm .* cscc.pos.trashpopup([3:4,3:4]));
    cscc.H.ETB.MNCrslider = javax.swing.JSlider; mncrvals = [10 1 50 30];
    javacomponent(cscc.H.ETB.MNCrslider,poss);
    set(cscc.H.ETB.MNCrslider,'Minimum',mncrvals(1),'Maximum',mncrvals(3),'Value',mncrvals(4),'snapToTicks',1,...
      'minorTickSpacing',mncrvals(2),'MajorTickSpacing',mncrvals(1),'PaintTicks',true, 'PaintLabels',false); 
    set(cscc.H.ETB.MNCrslider,'StateChangedCallback',@estimateExclusion);
    uicontrol('Parent',cscc.H.ETB.main,'Style','text','position',post,'String','standard deviation limit');
    for i=1:5
      uicontrol('Parent',cscc.H.ETB.main,'Style','text','String',sprintf('%d',i),'Fontsize',9,...
            'Units','normalized','position', [ posm(1:2) + [0.007 + (i-1)*0.1056 -0.03] 0.03 0.03]);
    end
    
    % PIQR exclusion slider
    posm = [0.51 0.26 0.465 0.06]; 
    post = round((posm +[0 0.068 0 -0.03]) .* cscc.pos.trashpopup([3:4,3:4])); 
    poss = round(posm .* cscc.pos.trashpopup([3:4,3:4]));
    uicontrol('Parent',cscc.H.ETB.main,'Style','text','position',post,'String','standard deviation limit','ForegroundColor',Color);
    cscc.H.ETB.PIQRrslider = javax.swing.JSlider; mncrvals = [10 1 50 30];
    javacomponent(cscc.H.ETB.PIQRrslider,poss);
    set(cscc.H.ETB.PIQRrslider,'Minimum',mncrvals(1),'Maximum',mncrvals(3),'Value',mncrvals(4),'snapToTicks',1,...
      'minorTickSpacing',mncrvals(2),'MajorTickSpacing',mncrvals(1),'PaintTicks',true, 'PaintLabels',false); %,...
    set(cscc.H.ETB.PIQRrslider,'StateChangedCallback',@estimateExclusion);
    if ~useIQR, cscc.H.ETB.PIQRrslider.disable; end
    for i=1:5
      uicontrol('Parent',cscc.H.ETB.main,'Style','text','String',sprintf('%d',i),'Fontsize',9,...
            'Units','normalized','position', [ posm(1:2) + [0.007 + (i-1)*0.1056 -0.03] 0.03 0.03]);
    end
  
    % === table and buttongroup ===
    cscc.H.ETB.bg = uibuttongroup(cscc.H.ETB.main,...
      'Position',[-0.01 -0.01 1.02 0.22]);

    % table with rows for total, passed and failed scans
    % with MNC, IQR PIQR, Euler, Nuisance1-3 as for loop table
    posm = [0.03 0.03 0 0.04];
    clear tab; 
    tab(1,1) = struct('str','Group'    ,'col',[0   0   0  ],'width',0.09,'ad',0); 
    tab(2,1) = struct('str','Total:'   ,'col',[0   0   0  ],'width',0.09,'ad',0); 
    tab(3,1) = struct('str','Accepted:','col',[0   0.5 0  ],'width',0.09,'ad',0); 
    tab(4,1) = struct('str','Excluded:','col',[0.9 0   0  ],'width',0.09,'ad',0); 
    
    tab(1,2) = struct('str','Ratio'   ,'col',[0   0   0  ],'width',...
                 0.05 + 0.02*cscc.datagroups.n_samples,'ad',0.01); 
    tab(2,2) = struct('str',sprintf('100%% (%d)',cscc.datagroups.n_subjects),...
      'col',[0 0 0],'width',0,'ad',0); 
    for i=1:min(10,cscc.datagroups.n_samples)
      tab(1,end+1) = struct('str',sprintf('G%d',i),'col',[0 0 0],'width',...
        0.015*ceil(cscc.datagroups.n_subjects^(1/10)),'ad',0.02*(i==1)); 
      tab(2,end) = struct('str',sprintf('%d',sum(cscc.datagroups.sample==i)),...
        'col',[0 0 0],'width',0,'ad',0); 
    end
    
    % mean correlation
    tab(1,end+1) = struct('str','MNC'     ,'col',[0   0   0  ],'width',0.06,'ad',0.02); 
    tab(2,end)   = struct('str',sprintf('%0.3f',mean( cscc.data.mean_cov)),...
      'col',[0   0   0  ],'width',0.06,'ad',0.02); 
    
    % xml variables
    if useIQR
      tab(1,end+1) = struct('str','IQR','col',[0   0   0  ],'width',0.06,'ad',0); 
      tab(2,end)   = struct('str',sprintf('%0.2f',mean( cscc.data.QM(:,4))),...
        'col',[0   0   0  ],'width',0.06,'ad',0.02); 
      tab(1,end+1) = struct('str','PIQR','col',[0   0   0  ],'width',0.06,'ad',0); 
      tab(2,end)   = struct('str',sprintf('%0.2f',mean( cscc.data.QM(:,5))),...
        'col',[0   0   0  ],'width',0.06,'ad',0.02); 
      if size(cscc.data.QM,2)>5
        tab(1,end+1) = struct('str','Euler','col',[0   0   0  ],'width',0.06,'ad',0); 
        tab(2,end)   = struct('str',sprintf('%0.2f',mean( cscc.data.QM(:,6))),...
          'col',[0   0   0  ],'width',0.06,'ad',0.02); 
      end
    end
    
    % nuisance variables
    for i=1:min(3,numel(cscc.job.c))
      tab(1,end+1) = struct('str',sprintf('N%d',i),'col',[0 0 0],'width',0.07,'ad',0.02*(i==1)); 
      tab(2,end)   = struct('str',sprintf('%0.2f',mean(cscc.job.c{i})),'col',[0 0 0],'width',0,'ad',0); 
    end
    
    % scale table 
    adx = (0.8 - (sum([tab(1,:).width]) + sum([tab(1,:).ad])) ) / size(tab,2);
    for j=2:size(tab,2), tab(1,j).ad = tab(1,j).ad + adx; end
    
    % print table
    for i=1:size(tab,1)
      for j=1:size(tab,2)
        cscc.H.ETB.tab(i,j) = uicontrol('Parent',cscc.H.ETB.main,'Style','text',...
          'position',round((posm + [ sum([tab(1,1:j-1).width]) + sum([tab(1,1:j).ad]) ...
            (4-i)*0.035 tab(1,j).width 0]) .* cscc.pos.trashpopup([3:4,3:4])),...
          'HorizontalAlignment','right','Fontsize',cscc.display.FS(cscc.display.FSi+4),...
          'ForegroundColor',tab(i,1).col); 
        if ~isempty(tab(i,j).str), set(cscc.H.ETB.tab(i,j),'String',tab(i,j).str); end
        if ~isempty(tab(i,j).col), set(cscc.H.ETB.tab(i,j),'ForegroundColor',tab(i,j).col); end
        if i==1, set(cscc.H.ETB.tab(i,j),'FontWeight','bold','ForegroundColor',[0.3 0.3 0.3]); end
        
      end
    end
    %
    estimateExclusion;
   

    posm = [0.91 0.035 0.06 0.07];
    % add button to open the help HTML

    cscc.H.ETB.auto = uicontrol(cscc.H.ETB.main,...
      'Units','normalized','position',posm + [-posm(3)+0.005 posm(4)-0.01 0 0],'Style','Pushbutton',...
      'string','Auto','ToolTipString','Set defaults','ForegroundColor',[0 0.8 0],...
      'FontSize',cscc.display.FS(cscc.display.FSi+4),'enable','off',...
      'callback','');

    % add button to open the help HTML
    cscc.H.ETB.help = uicontrol(cscc.H.ETB.main,...
      'Units','normalized','position',posm + [0 posm(4)-0.01 0 0],'Style','Pushbutton',...
      'string','Help','ToolTipString','Open help of this window','ForegroundColor',[0 0 0.8],...
      'FontSize',cscc.display.FS(cscc.display.FSi+4),'callback',...
      ['global cscc; ' ...
       'if ~isfield(cscc.H,''helpfig'') || ~isvalid(cscc.H.helpfig) ' ...
       '  cscc.H.helpfig = spm_figure(''Create'',''CAThelp'',''CAT help''); ' ...
       '  cscc.H.helpfig.Position = cscc.H.graphics.Position;' ...
       'end; ' ...
       'spm_help(''!Disp'',fullfile(spm(''Dir''),''toolbox'',''cat12'',''html'...
       ',''cat_tools_checkcov.html'','''',cscc.H.helpfig);']);

    % reject button
    closefunction = ['global cscc; ' ...
       'delete(cscc.H.ETB.main); ' ... 
       'cscc.H = rmfield(cscc.H,''ETB'');' ...
       'clearvars -global at;'];
    cscc.H.ETB.reject = uicontrol(cscc.H.ETB.main,...
      'Units','normalized','position',posm + [-posm(3)+0.005 0 0 0],'Style','Pushbutton',...
      'string','Reject','ToolTipString','Close window without settings','ForegroundColor',[0.8 0 0],...
      'FontSize',cscc.display.FS(cscc.display.FSi+4),'callback', closefunction );
    cscc.H.ETB.main.CloseRequestFcn = closefunction;

    % apply button
    cscc.H.ETB.apply = uicontrol(cscc.H.ETB.main,...
      'Units','normalized','position',posm,'Style','Pushbutton',...
      'string','Apply','ToolTipString','Apply settings and close window','ForegroundColor',[0 0.6 0 ],...
      'FontSize',cscc.display.FS(cscc.display.FSi+4),'callback',@autotrash);

    %
    buttonicon(cscc.H.ETB.auto   ,'Auto'   ,fullfile(spm('dir'),'toolbox','cat12','html','icons','trash_auto.png'));
    buttonicon(cscc.H.ETB.help   ,'Help'   ,fullfile(spm('dir'),'toolbox','cat12','html','icons','status_help.png'));
    buttonicon(cscc.H.ETB.reject ,'Reject' ,fullfile(spm('dir'),'toolbox','cat12','html','icons','status_failed.png'));
    buttonicon(cscc.H.ETB.apply  ,'Apply'  ,fullfile(spm('dir'),'toolbox','cat12','html','icons','status_passed.png'));
  
  
  

return
%-----------------------------------------------------------------------

function autotrash(obj, event_obj) 
%-----------------------------------------------------------------------
% Subfunction of cat_check_cov2 that calls another GUI (autotrashGUI) or
% directly estimates the worst datasets that should be excluded. Uses an
% additional global variable "at" for default parameter and the new 
% exclusion list that is finally added to the trashlist.
%-----------------------------------------------------------------------
  global cscc at

  % set default values
  if ~exist('at','var') || isempty(at)
    at.MNCath  = 0.05;
    at.MNCrth  = 2;
    at.PIQRath = 4;
    at.PIQRrth = 2;

    if ~isfield(at,'del') && 1 %cscc.job.expertgui 
      autotrashGUI; 
      return
    else
      % estimate list
      estimateExclusion
      if isempty(at.del)
        fprintf('Nothing to delete.\n');
      end
    end
  end
  
  
  if isfield(at,'del') && ~isempty(at.del)
    if isfield(cscc.pos,'x'), oldx = cscc.pos.x; end
    if isfield(cscc.pos,'y'), oldy = cscc.pos.y; cscc.H.y = rmfield(cscc.pos,'y'); end

    % remove the different datasets
    for di=1:numel(at.del)
      cscc.pos.x = at.del(di);
      trash
    end
    
    % sort by MD
    [ss,si] = sort(cscc.data.MD(at.del),'descend'); 
    at.del = at.del(si);
    
    cscc.select.trashlist  = [cscc.select.trashlist, at.del];
    cscc.select.trashhist  = [cscc.select.trashhist, at.del];
    cscc.select.trashhistf = []; 
    
    if exist('oldx','var'), cscc.pos.x = oldx; end
    if exist('oldy','var'), cscc.pos.y = oldy; end
    
    set([cscc.H.trashui.new,cscc.H.trashui.disptrash,cscc.H.trashui.undo...
      cscc.H.trashui.autotrash,cscc.H.trashui.ziptrash],'enable','on');
    set(cscc.H.trashui.redo,'enable','off');
  
    % update boxplot
    show_boxplot
    
    % set fields
    if isfield(cscc.H,'showtrash') && ishandle(cscc.H.showtrash), set(cscc.H.showtrash,'Enable','on'); end
  end
  
  % cleanup
  if isfield(cscc.H,'ETB')
    if isfield(cscc.H.ETB,'main') && isvalid(cscc.H.ETB.main )
      delete(cscc.H.ETB.main)
    end
    cscc.H = rmfield(cscc.H,'ETB'); 
  end
  
  clearvars -global at;
return

%-----------------------------------------------------------------------
function trash(obj, event_obj) 
%-----------------------------------------------------------------------
% Puts a record on the trash list and marks it with a red cross in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global cscc 
  
  if isfield(cscc.pos,'x') && all( cscc.select.trashlist~=cscc.pos.x ) 
    if exist('obj','var')
      cscc.select.trashlist = [cscc.select.trashlist cscc.pos.x];
      cscc.select.trashhist = [cscc.select.trashhist cscc.pos.x];
    end
    
    showtrash = get(cscc.H.showtrash,'Value');
    if ~showtrash && ~cscc.H.isscatter
      cscc.H.dcm.removeAllDataCursors;
    end
    
    if exist('obj','var')
      set(cscc.H.trashui.undo     ,'Enable','on' );
      set(cscc.H.trashui.new      ,'Enable','on' );
      set(cscc.H.trashui.disptrash,'Enable','on' );
      set(cscc.H.trashui.ziptrash ,'Enable','on' );
      set(cscc.H.trashui.undo     ,'Enable','on' );
      set(cscc.H.trashui.redo     ,'Enable','off' );
    end
    if ~cscc.H.isscatter
      set(cscc.H.trashui.trashcol  ,'Enable','off');
      set(cscc.H.trashui.detrashcol,'Enable','on' );
      set(cscc.H.showtrash         ,'Enable','on' );
    else
      set(cscc.H.trashui.trash    ,'Enable','off');
      set(cscc.H.trashui.detrash  ,'Enable','on' );
    end
    
    cscc.select.trash1d(cscc.pos.x)    = 0;
    cscc.select.trash2d(cscc.pos.x,:)  = 0;
    cscc.select.trash2d(:,cscc.pos.x)  = 0;
    
    for scati=1:numel(cscc.H.scat)
      sc.posx  = findobj(cscc.H.scat(scati),'type','scatter'); 
      sc.posxv = cell2mat(get(sc.posx,'UserData'));
      sc.posxi = find(sc.posxv==cscc.pos.x,1,'first');

      set(sc.posx(sc.posxi),'sizedatasource',...
        get(sc.posx(sc.posxi),'marker'),...
        'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
        'ZDataSource','trash','MarkerFaceAlpha',0);
    end
    if ~cscc.H.isscatter
      update_matrix
    end
    
  
  end
return

%-----------------------------------------------------------------------
function trashrow(obj, event_obj) 
%-----------------------------------------------------------------------
% Puts a record on the trash list and marks it with a red cross in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global cscc
  
  if isfield(cscc.pos,'y') && all( cscc.select.trashlist~=cscc.pos.y ) 
    cscc.select.trashlist = [cscc.select.trashlist cscc.pos.y];
    cscc.select.trashhist = [cscc.select.trashhist cscc.pos.y];
  
    set(cscc.H.trashui.trashrow  ,'Enable','off');
    set(cscc.H.trashui.detrashrow,'Enable','on' );
    set(cscc.H.trashui.ziptrash  ,'Enable','on' );
    set(cscc.H.trashui.undo      ,'Enable','on' );
    set(cscc.H.showtrash         ,'Enable','on' );
   
    cscc.select.trash1d(cscc.pos.y)    = 0;
    cscc.select.trash2d(cscc.pos.y,:)  = 0;
    cscc.select.trash2d(:,cscc.pos.y)  = 0;
    
    for scati=1:numel(cscc.H.scat)
      sc.posx  = findobj(cscc.H.scat,'type','scatter'); 
      sc.posxv = cell2mat(get(sc.posx,'UserData'));
      sc.posxi = find(sc.posxv==cscc.pos.y,1,'first');

      set(sc.posx(sc.posxi),'sizedatasource',...
        get(sc.posx(sc.posxi),'marker'),...
        'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
        'ZDataSource','trash','MarkerFaceAlpha',0);
    end
    
    if ~cscc.H.isscatter
      update_matrix;
    end
   
  end
return

%-----------------------------------------------------------------------
function detrash(obj, event_obj)
%-----------------------------------------------------------------------
% Removes a record from trash list and restores the old look like in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global cscc 
  
  if isfield(cscc.pos,'x') && any( cscc.select.trashlist==cscc.pos.x ) 
    if exist('obj','var')
      cscc.select.trashlist = setdiff(cscc.select.trashlist,cscc.pos.x);
      cscc.select.trashhist = [cscc.select.trashhist -cscc.pos.x];

      if cscc.H.isscatter
        set(cscc.H.trashui.trash  ,'Enable','on' );
        set(cscc.H.trashui.detrash,'Enable','off');
      else
        set(cscc.H.trashui.trashcol  ,'Enable','on' );
        set(cscc.H.trashui.detrashcol,'Enable','off');
      end
    end

    % update matrix
    cscc.select.trash1d(cscc.pos.x)    = 1;
    cscc.select.trash2d(cscc.pos.x,:)  = sum(cscc.select.trash2d,1)>0;
    cscc.select.trash2d(:,cscc.pos.x)  = sum(cscc.select.trash2d,2)>0;

    % update scatter
    for scati=1:numel(cscc.H.scat)
      sccscc.posx  = findobj(cscc.H.scat(scati),'type','scatter'); 
      sccscc.posxv = cell2mat(get(sccscc.posx,'UserData'));
      sccscc.posxi = find(sccscc.posxv==cscc.pos.x,1,'first');

      set(sccscc.posx(sccscc.posxi),'marker',...
        get(sccscc.posx(sccscc.posxi),'sizedatasource'),... 
        'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
    end

    if ~cscc.H.isscatter 
      update_matrix;
      if isempty(cscc.select.trashlist), set(cscc.H.showtrash,'Enable','off'); end
    end

    if isempty(cscc.select.trashlist)
      set([cscc.H.showtrash,cscc.H.trashui.ziptrash],'Enable','off');
    end
  end
return

%-----------------------------------------------------------------------
function detrashrow(obj, event_obj)
%-----------------------------------------------------------------------
% Removes a record from trash list and restores the old look like in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global cscc 

  if isfield(cscc.pos,'y') && any( cscc.select.trashlist==cscc.pos.y ) 
    cscc.select.trashlist = setdiff(cscc.select.trashlist,cscc.pos.y);
    cscc.select.trashhist = [cscc.select.trashhist -cscc.pos.y];
    
    set(cscc.H.trashui.trashrow  ,'Enable','on' );
    set(cscc.H.trashui.detrashrow,'Enable','off');
    
    % update matrix
    cscc.select.trash1d(cscc.pos.y)    = 1;
    cscc.select.trash2d(cscc.pos.y,:)  = sum(cscc.select.trash2d,1)>0;
    cscc.select.trash2d(:,cscc.pos.y)  = sum(cscc.select.trash2d,2)>0;
    
    % update scatter
    for scati=1:numel(cscc.H.scat)
      sccscc.posx  = findobj(cscc.H.scat(scati),'type','scatter'); 
      sccscc.posxv = cell2mat(get(sccscc.posx,'UserData'));
      sccscc.posxi = find(sccscc.posxv==cscc.pos.y,1,'first');

      marker = get(sccscc.posx(sccscc.posxi),'sizedatasource');
      if isempty(strfind(marker,'+o*.xsdv^><ph')), marker = '.'; end
      set(sccscc.posx(sccscc.posxi),'marker',marker,... 
        'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
    end
    
    if ~cscc.H.isscatter 
      update_matrix;
      if isempty(cscc.select.trashlist), set(cscc.H.showtrash,'Enable','off'); end
    end
    
    if isempty(cscc.select.trashlist)
      set([cscc.H.showtrash,cscc.H.trashui.ziptrash],'Enable','off');
    end
  end
return

%-----------------------------------------------------------------------
function newtrash(obj, event_obj)
%-----------------------------------------------------------------------
% Create an empty trash list. 
%-----------------------------------------------------------------------
  global cscc 
  cscc.select.trashlist  = [];
  cscc.select.trashhist  = []; 
  cscc.select.trashhistf = []; 

  cscc.select.trash2d = true(size(cscc.select.trash2d));
  cscc.select.trash1d = true(size(cscc.select.trash1d));
   
  % find scatter objects 
  sc = findobj('ZDataSource','trash');
  for sci=1:numel(sc)
   set(sc(sci),'marker',get(sc(sci),'sizedatasource'),... 
        'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
  end
  
  if ~cscc.H.isscatter 
    update_matrix;
  end
  
  %cscc.select.trashhist = [cscc.select.trashhist -cscc.select.trashlist];
  %cscc.H.trashui.undo,cscc.H.trashui.redo,
  %cscc.H.trashui.trash,cscc.H.trashui.trashrow,,cscc.H.showtrasht
  unit = struct2cell(cscc.H.checkui); 
  set([unit{:}],'Enable','off');
  set([cscc.H.trashui.trash,cscc.H.trashui.detrash,...
       cscc.H.trashui.trashcol,cscc.H.trashui.detrashcol,...
       cscc.H.trashui.trashrow,cscc.H.trashui.detrashrow,...
       cscc.H.trashui.undo,cscc.H.trashui.redo,cscc.H.showtrash,...
       cscc.H.trashui.disptrash,cscc.H.trashui.ziptrash],'Enable','off');
  if ~isempty(cscc.H.dcm.getCursorInfo)
    if isfield(cscc.pos,'x'), set([cscc.H.trashui.trash,cscc.H.trashui.trashcol],'Enable','on'); end
    %if isfield(cscc.pos,'y'), set(cscc.H.trashui.trashrow,'Enable','on'); end
  end
return

%-----------------------------------------------------------------------
function disptrash(obj, event_obj)
%-----------------------------------------------------------------------
% List all records of the trash list in the command window.
%-----------------------------------------------------------------------
  global cscc

  fprintf('\nExcluded files:\n');
  for fi=1:numel(cscc.select.trashlist)
    [pp,ff,ee,dd] = spm_fileparts(cscc.files.data{cscc.select.trashlist(fi)}); 
    fprintf('  %s\n',fullfile(pp,[ff ee]));
  end
return

%-----------------------------------------------------------------------
function trashundo(obj, event_obj)
%-----------------------------------------------------------------------
% List all records of the trash list in the command window.
%-----------------------------------------------------------------------
  global cscc
  
  if isfield(cscc.pos,'x'), oldx = cscc.pos.x; end
  if isfield(cscc.pos,'y'), oldy = cscc.pos.y; end
  if isfield(cscc.pos,'tar_mouseget'), oldtar_mouseget = cscc.pos.tar_mouseget; end
  cscc.pos.x = abs(cscc.select.trashhist(end)); 
  
  if ~isempty(cscc.select.trashlist) && cscc.select.trashhist(end)>0
    % last element was added to the cscc.select.trashlist 
    detrash
    cscc.select.trashlist(end) = [];
  else
    trash
    cscc.select.trashlist(end+1) = cscc.select.trashhist(end);
  end
  cscc.select.trashhistf(end+1) = cscc.select.trashhist(end); 
  cscc.select.trashhist(end)    = []; 
  
  set(cscc.H.trashui.redo,'enable','on');
  if isempty(cscc.select.trashhist), set(cscc.H.trashui.undo,'enable','off'); end
  
  if exist('oldx','var')
    cscc.pos.x = oldx; 
  else
    cscc.pos = rmfield(cscc.pos,'x'); 
  end
  if exist('oldy','var')
    cscc.pos.y = oldy; 
  elseif isfield(cscc.pos,'y') 
    cscc.pos = rmfield(cscc.pos,'y'); 
  end
  if exist('oldtar_mouseget','var')
    cscc.pos.tar_mouseget = oldtar_mouseget; 
  elseif isfield(cscc.pos,'tar_mouseget')
    cscc.pos = rmfield(cscc.pos,'tar_mouseget'); 
  end
  if isempty(cscc.select.trashlist)
    set([cscc.H.showtrash,cscc.H.trashui.ziptrash],'Enable','off');
  else
    set([cscc.H.showtrash,cscc.H.trashui.ziptrash],'Enable','on');
  end
return

%-----------------------------------------------------------------------
function trashredo(obj, event_obj)
%-----------------------------------------------------------------------
% List all records of the trash list in the command window.
%-----------------------------------------------------------------------
  global cscc
  
  if isfield(cscc.pos,'x'), oldx = cscc.pos.x; end
  if isfield(cscc.pos,'y'), oldy = cscc.pos.y; end
  if isfield(cscc.pos,'tar_mouseget'), oldtar_mouseget = cscc.pos.tar_mouseget; end
  cscc.pos.x = cscc.select.trashhistf(end); 
  
  if isempty(cscc.select.trashlist) && cscc.select.trashhistf(end)<0
    detrash
    cscc.select.trashlist(end) = [];
  else
    trash
    cscc.select.trashlist(end+1)  = cscc.select.trashhistf(end);
  end
  cscc.select.trashhist(end+1) = cscc.select.trashhistf(end); 
  cscc.select.trashhistf(end)  = [];
  
  set(cscc.H.trashui.undo,'enable','on');
  if isempty(cscc.select.trashhistf), set(cscc.H.trashui.redo,'enable','off'); end
  
  if exist('oldx','var')
    cscc.pos.x = oldx; 
  else
    cscc.pos = rmfield(cscc.pos,'x'); 
  end
  if exist('oldy','var')
    cscc.pos.y = oldy; 
  elseif isfield(cscc.pos,'y') 
    cscc.pos = rmfield(cscc.pos,'y'); 
  end
  if exist('oldtar_mouseget','var')
    cscc.pos.tar_mouseget = oldtar_mouseget; 
  elseif isfield(cscc.pos,'tar_mouseget')
    cscc.pos = rmfield(cscc.pos,'tar_mouseget'); 
  end
  if isempty(cscc.select.trashlist)
    set([cscc.H.showtrash,cscc.H.trashui.ziptrash],'Enable','off');
  else
    set([cscc.H.showtrash,cscc.H.trashui.ziptrash],'Enable','on');
  end
return

%-----------------------------------------------------------------------
function ziptrash(obj, event_obj)
%-----------------------------------------------------------------------
% Remove records and related files from the file system by zipping or
% storing in a separate directory (NOT READY). 
%-----------------------------------------------------------------------
  global cscc 
  
  % dialog box with question to proceed
  d = dialog('Position',cscc.pos.popup,'Name','Remove unfitting data',...
    'Position',[cscc.H.figure.Position(1:2) + cscc.H.figure.Position(3:4).*[1 0.7] - ...
    cscc.pos.popup(3:4)/2 cscc.pos.popup(3:4)]);
  uicontrol('Parent',d,'Style','text','Position',[20 60 160 20],...
     'String',{sprintf('Copy %d scans into trash directory!',numel(cscc.select.trashlist)),...
     });
  uicontrol('Parent',d,'TooltipString','Stop operation!','selected','on',...
     'Position',[20 20 50 25],'String','No','Foregroundcolor',[ 0.8 0 0 ],...
     'Callback','delete(gcf);');    
  uicontrol('Parent',d,'TooltipString','Go on!','selected','on',...
     'Position',[75 20 50 25],'String','Yes','Foregroundcolor',[ 0 0.6 0 ],...
     'Callback',@emptytrash); 
  uicontrol('Parent',d,'TooltipString','More information to this operation',...
     'Position',[130 20 50 25],'String','Help','Foregroundcolor',[ 0 0 0.8 ],...
     'Callback',['delete(gcf); spm_help(''!Disp'',fullfile(spm(''Dir''),''toolbox'...
     ',''cat12'',''html'',''cat_tools_checkcov.html'','''',spm_figure(''GetWin'',''Graphics'')); ']); %    
return

%-----------------------------------------------------------------------
function emptytrash(obj, event_obj)
%-----------------------------------------------------------------------
% Remove records and related files from the file system by zipping or
% storing in a separate directory (NOT READY). 
%-----------------------------------------------------------------------
  global cscc 
 
  close(get(obj,'Parent'))
  set(cscc.H.figure,'Visible','off')
  
  testtrash = 0;
  
  % unique main trashsubdir
  tid = char([48:57 65:90 97:122]); 
  newtrashdir = [datestr(clock,'yyyymmdd_HHMMSS') '_' ...
    tid(floor( rand(1,6)*numel(tid-1) + 1))];
  newtrashfile = fullfile(cscc.trashdir,...
    ['restore_' newtrashdir '.m']);
  
  % otherwise we go and by creating the new specific trash directory
  try mkdir(fullfile(cscc.trashdir,newtrashdir)); end
  
    
  %% trash
  
    
  
%-----------------------------------------------------------------------
% open popup to avoid interaction
%-----------------------------------------------------------------------
%%
  spm_figure('GetWin','Interactive');
  
  if 0
    popup = dialog('Position',cscc.pos.popup,'Name','Remove unfitting data',...
        'Position',[cscc.H.figure.Position(1:2) + cscc.H.figure.Position(3:4).*[1 0.7] - ...
        cscc.pos.popup(3:4)/2 cscc.pos.popup(3:4)]); 
    uicontrol('Parent',popup,'Style','text','Position',[20 30 160 40],...
       'String',{sprintf('Copy scans into trash directory!',numel(cscc.select.trashlist)),...
       'Do not interrupt!','This window close when finished.'});
  end
  
  %%
  spm_progress_bar('Init',numel(cscc.select.trashlist),'Search related files','subjects completed')
  spm_figure('GetWin','Interactive');
  for fi=1:numel(cscc.select.trashlist)
    %% find preprocessed data and create data structure for store thrash 
    [pp,ff,ee]           = spm_fileparts(cscc.files.org{cscc.select.trashlist(fi)});
    trash(fi).dir        = pp;
    trash(fi).fname      = cscc.files.org{cscc.select.trashlist(fi)}; 
    trash(fi).MNC        = cscc.data.mean_cov(cscc.select.trashlist(fi));
    trash(fi).MD         = cscc.data.MD(cscc.select.trashlist(fi));
    try
      trash(fi).IQR      = cscc.data.QM(cscc.select.trashlist(fi),4);
      trash(fi).PIQR     = cscc.data.QM(cscc.select.trashlist(fi),5);
    end
    
    
    %% find related files
    %  --------------------------------------------------------------------
    %  This is the most difficult part that has no perfect solution due to 
    %  the different types how data could be stored (see also finding the 
    %  original image). The are two major ways to store data (i) with no 
    %  subdirectories or (ii) with many. In case of many directories, the 
    %  filesname are possible identical and the subject specific
    %  information is stored in the directory names. Names could further
    %  include a pre and suffix with e.g. protocol information etc. 
    %    1) ../[pre][name][post].* 
    %    2) ../dir_1/../dir_n/t1.*
    %  
    %  Steps: 
    %   1) find all files that have include the original file name
    %   2) remove files that include this name more than once, ie. in set 
    %      of scans {"1.nii", "2.nii", ... , "11.nii", ...} the file "11"
    %      is not a processed version of "1".
    %   3) ... 
    %      
    %  --------------------------------------------------------------------
    % 1) find all files with the main filename
    trash(fi).sim_files  = cat_vol_findfiles(pp,['*' ff '*.nii'],struct('maxdepth',1)); 
    % 2) find files that inlude the main filename more than once and remove them
    trash(fi).sim_files0 = cat_vol_findfiles(pp,['*' ff '*' ff '*.nii'],struct('maxdepth',1)); 
    trash(fi).sim_files  = setdiff(trash(fi).sim_files,trash(fi).sim_files0); 
    % 3) ...
    trash(fi).sim_files  = [fullfile(pp,[ff ee]); setdiff(trash(fi).sim_files, fullfile(pp,[ff ee]))]; 
    for si=1:numel(trash(fi).sim_files); 
      [pps,ffs]    = spm_fileparts(trash(fi).sim_files{si});
      pp_files{si} = cat_vol_findfiles(pps,['*' ffs '*']);
      pp_files0    = cat_vol_findfiles(pps,['*' ffs '*' ffs '*']);
      pp_files{si} = setdiff(pp_files{si},pp_files0);
      pp_files{si} = setdiff(pp_files{si},trash(fi).sim_files);
      for fsi = numel(pp_files{si}):-1:1
        ppfs = spm_str_manip(pp_files{si}{fsi},'hht');
        switch ppfs
          case 'err', pp_files{si}(fsi) = [];
        end
      end
      pp_files{si} = setdiff(pp_files{si},trash(fi).sim_files{si}); 
      if si==1
        pp_filescor = pp_files{si};
      else
        pp_filescor = setdiff(pp_filescor,pp_files{si}); 
      end
    end
    trash(fi).sim_files = [trash(fi).sim_files; pp_filescor];
    spm_progress_bar('Set',fi);
  end
  spm_progress_bar('Clear');
  
  % replace leading path by the new trashdirectory
  [tmp,grouphome] = spm_str_manip(cscc.files.fname.s,'C');
  for fi=1:numel(cscc.select.trashlist)
    trash(fi).sim_filest = cat_io_strrep(trash(fi).sim_files,...
      grouphome.s,fullfile(cscc.trashdir,[newtrashdir filesep]));
  end
  
  
  %% create restore script
   
  % start with some help 
  script = {
   ['%% restore_' newtrashdir '.m'] 
    '%  -------------------------------------------------------------------'
    '%  This is an automaticly generated script to restore files that were '
    '%  previously removed in a quality control process based on the mean'
    '%  covariance (MNC) and the protocol image qualtiy rating (PIQR). The'
    '%  files are stored in the equaly named subdirectory and listed below.'
    '%  ------------------------------------------------------------------'
    '%  '
    ' '
    '%% List of scans: '
    '%  ------------------------------------------------------------------'
    '%  Structure with the name of the removed scan, its mean covariance '
    '%  (MNC), image quality rating (IQR), the protocol image qualtiy '
    '%  rating (PIQR), Mahanalobis distance and file to file lists one ' 
    '%  of the orinal files and one of the trashed files.'
    '%  ------------------------------------------------------------------'
  };

  % add a data structur that include major information especial the filelists
  for fi=1:numel(cscc.select.trashlist)
    script{end+1} = sprintf('scan(%d).fname = ''%s'';',fi,trash(fi).fname);
    script{end+1} = sprintf('scan(%d).MNC   = %0.4f;' ,fi,trash(fi).MNC  );
    script{end+1} = sprintf('scan(%d).IQR   = %0.2f;' ,fi,trash(fi).IQR  );
    script{end+1} = sprintf('scan(%d).PIQR  = %0.2f;' ,fi,trash(fi).PIQR );
    script{end+1} = sprintf('scan(%d).MD    = %0.2f;' ,fi,trash(fi).MD   );
    script{end+1} = sprintf('scan(%d).files = {',fi);
    for fii=1:numel(trash(fi).sim_files)
      if testtrash
        script{end+1} = sprintf('  ''%s'';',cat_io_strrep(trash(fi).sim_files{fii},...
          grouphome.s,fullfile(cscc.trashdir,['restore' filesep])));
      else
        script{end+1} = sprintf('  ''%s'';',trash(fi).sim_files{fii});
      end
    end
    script{end+1} = '  };';
    script{end+1} = sprintf('scan(%d).tfiles = {',fi);
    for fii=1:numel(trash(fi).sim_filest)
      script{end+1} = sprintf('  ''%s'';',trash(fi).sim_filest{fii});
    end
    script{end+1} = '  };';
  end
  
  
%% FUTURE DEVELOPMENT
%  ------------------------------------------------------------------------
% the idea is to add a GUI with a menu that allows do select different
% scans for restoring
%{
script(end+1) = {
''
'' % ask for type (1) 'restore all' (2) 'deside each by question' (3) 'deside each by id-area'
'' % sort by quality?
''
%}
%  ------------------------------------------------------------------------



  % this part now describes the major operation for restore the previously
  % defined filelists
  script = [script;{
    '%% restore files'
    '%  ------------------------------------------------------------------'
   ['fprintf(''Restore %d files of ' newtrashdir ':\\n'',numel(scan));']
    'for si = 1:numel(scan)'
    '  fprintf(''  Restore "%s":'',scan(si).fname ); si_err=0;'
    '  for i = 1:numel(scan(si).files)'
    '    if exist(scan(si).tfiles{i},''file'') && ... '
    '      ~exist(scan(si).files{i},''file'')'
    '      pp = spm_fileparts(scan(si).files{i});'
    '      if ~exist(pp,''dir''), mkdir(pp); end'
    '      movefile(scan(si).tfiles{i},scan(si).files{i});'
    '    else'
    '      % error message? '
    '      si_err = 1;'
    '    end'
    '  end'
    '  if si_err'
    '    cat_io_cprintf(''err'','' ERROR.\\n'');'
    '  else'
    '    cat_io_cprintf([0 0.8 0],'' OK.\\n'');'
    '  end'
    'end'
    ''
   ['cat_io_rmdir(''' fullfile(cscc.trashdir,newtrashdir) ''')']
    ''
    '%  ------------------------------------------------------------------'
    '%  end of script'
  }];
   
    
  % save the final script
  fid = fopen(newtrashfile,'w');
  for li=1:numel(script), fprintf(fid,[strrep(script{li},'%','%%') '\n']); end
  fclose(fid); 
  
  
  %%  
  spm_progress_bar('Init',numel(cscc.select.trashlist),'Trash related files','subjects completed')
  spm_figure('GetWin','Interactive');
  for fi=1:numel(cscc.select.trashlist)
    for fii=1:numel(trash(fi).sim_files)
      fprintf('  Trash: %s\n',trash(fi).fname)
      pp = spm_fileparts(trash(fi).sim_filest{fii});
      if ~exist(pp,'dir'), mkdir(pp); end
      if testtrash
        copyfile(trash(fi).sim_files{fii},trash(fi).sim_filest{fii});
      else
        movefile(trash(fi).sim_files{fii},trash(fi).sim_filest{fii});
      end
    end
    
    %% zip the list and remove the files
    %  we need to go into the directory and use the short filenames 
    %  to opbtain save the relative path in the zip file!
    if 0
      pp_filescor1 = cat_io_strrep(pp_filescor,[pp filesep],''); 
      ffzip = fullfile(pp,sprintf('%s_cor%2.2f_IQR%2.2f_trashed%s',...
        ff,cscc.data.X(cscc.select.trashlist(fi),1),...
        cscc.data.QM(cscc.select.trashlist(fi),3),trashtime));
      zip(ffzip,pp_filescor1,pp); 
      for fii=1:numel(pp_filescor1), delete(pp_filescor1{fii}); end
    end
    
    spm_progress_bar('Set',fi);
  end
  spm_progress_bar('Clear');
  
% close popup
%delete(popup);
  

  %%
  for fi=1:numel(trash)
    cat_io_rmdir(trash(fi).dir)
  end
    
  %% restore cat_stat_check_cov with updated filelist
  job2 = cscc.job;
  for fi = sort(cscc.select.trashlist,'descend')
    job2.(cscc.files.datafield){cscc.files.dataid(fi,2)}(cscc.files.dataid(fi,3)) = []; 
    if ~isempty(job2.data_xml) && numel(job2.data_xml)>fi
      job2.data_xml(fi) = []; 
    end
    if ~isempty(job2.c) 
      for ci=1:numel(job2.c)
        job2.c{ci}(fi) = []; 
      end
    end
  end
  
  cat_stat_check_cov2(job2);
return

%-----------------------------------------------------------------------
function closeWindows(obj, event_obj)
%-----------------------------------------------------------------------
% Close all windows and remove variables.
%-----------------------------------------------------------------------
  global cscc
  
  if strcmp( event_obj.Source.Type , 'figure')
    cscc.posx = get(event_obj.Source,'Position');
  else
    cscc.posx = get(get(event_obj.Source,'Parent'),'Position');
  end
  cscc.pos.popup(1:2) = [cscc.posx(1) + cscc.posx(3)*0.8, cscc.posx(1) + cscc.posx(4)*0.9]; 
  
  if isfield(cscc,'select') && isfield(cscc.select,'trashlist') && ~isempty(cscc.select.trashlist)
    d = dialog('Position',cscc.pos.popup,'Name','Close Checkcov');
    uicontrol('Parent',d,'Style','text','Position',[20 60 160 20],...
       'String','Trashlist not empty!');
    uicontrol('Parent',d,'TooltipString','Sopt closing',...
       'Position',[25 20 70 25],'String','Cancel','Callback','delete(gcf)');    
    uicontrol('Parent',d,'TooltipString','Close windows without last changes',...
       'Position',[100 20 70 25],'String','Close','ForegroundColor',[0.8 0 0],...
       'Callback',[ % see standard closing in else case
          'global scss;'...
          'set(cscc.H.graphics,''Color'',[1 1 1]);' ... 
          'spm_figure(''Clear'',cscc.H.graphics); ' ...
          'set(cscc.H.graphics,''Position'',cscc.display.WP); ' ...
          'for i=3:26, try close(i); end; end;' ...
          'clearvars -GLOBAL cscc;']);   
  else
    % clear SPM graphics window
    try set(cscc.H.graphics,'Color',[1 1 1]); end
    try   
      spm_figure('Clear',cscc.H.graphics); 
      set(cscc.H.graphics,'Position',cscc.display.WP); 
    catch
      try 
        spm_figure('Clear',spm_figure('FindWin','Graphics'));
        set(spm_figure('FindWin','Graphics'),'Position',cscc.display.WP); 
      end
    end
    
    % clear other surface window figures
    for i=3:26, try close(i); end; end; 
    delete(gcbf); % remove figure with closerequest

    % clear vars
    clearvars -GLOBAL cscc; 
  end
return

%-----------------------------------------------------------------------
function id = mygetCursorInfo
%-----------------------------------------------------------------------
% 
%-----------------------------------------------------------------------
  global cscc

  curs = cscc.H.dcm.getCursorInfo;
  
  if cscc.H.isscatter
    sc = unique([curs(:).Target]);
    id = get(sc,'UserData')'; 
    if iscell(id), id = cell2mat(id); end
  else
    id = unique([cscc.pos.x, cscc.pos.y]);
    %{
    cscc.pos  = reshape([curs(:).Position],numel(curs),2);   
    cscc.posx = unique(cscc.pos(:,1));
    cscc.posy = unique(cscc.pos(:,2));
  
    if cscc.H.sorted
      
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
  global cscc 
  
  id = mygetCursorInfo;
  if all(~cellfun('isempty',cscc.files.jpg(id))) || numel(id)>2
    %%
    spm_figure('Clear',cscc.H.graphics);
    spm_figure('Focus',cscc.H.graphics);
    
    if cscc.H.isscatter || cscc.pos.x == cscc.pos.y
      ppos = [0 0 1 1];
      jpg  = imread(cscc.files.jpg{cscc.pos.x}); 
      set(gca,'position',ppos(1,:));
      gpos = cscc.H.graphics.Position; 
      [Xq,Yq] = meshgrid(1:size(jpg,2)/gpos(3)/2:size(jpg,2),...
                         1:size(jpg,1)/gpos(4)/2:size(jpg,1));
      jpgi = zeros([size(Xq,1) size(Xq,2) 3],'uint8');
      for i=1:3, jpgi(:,:,i) = uint8(interp2(single(jpg(:,:,i)),Xq,Yq,'linear')); end
      image(cscc.H.graphics.CurrentAxes,jpgi);
      set(gca,'Visible','off'); 
    else
      ppos = [0.0 0.502 1.0 0.498 ; 0.0 0.0 1.0 0.498];
      jpg  = {imread(cscc.files.jpg{cscc.pos.x});imread(cscc.files.jpg{cscc.pos.y})};
      for fi=1:2
        ax(fi) = subplot(2,1,fi); 
        set(ax(fi),'position',ppos(fi,:));
        gpos = cscc.H.graphics.Position; 
        [Xq,Yq] = meshgrid(1:size(jpg{fi},2)/gpos(3)/2:size(jpg{fi},2),...
                           1:size(jpg{fi},1)/gpos(4)/2:size(jpg{fi},1));
        jpgi = zeros([size(Xq,1) size(Xq,2) 3],'uint8');
        for i=1:3, jpgi(:,:,i) = uint8(interp2(single(jpg{fi}(:,:,i)),Xq,Yq,'linear')); end
        image(cscc.H.graphics.CurrentAxes,jpgi);
        set(gca,'Visible','off'); axis equal tight; zoom(2); axis fill
      end
      pan yon
      linkaxes(ax)
    end
  else
    for i=1:numel(id), open(cscc.files.pdf{id(i)}); end
  end
return

%-----------------------------------------------------------------------
function checksurf(obj, event_obj)
%-----------------------------------------------------------------------
% Open surface files of selected subjects. 
% This is very slow and some feedback would be useful. 
%-----------------------------------------------------------------------
  global cscc
  
  spm_progress_bar('Init',2 - cscc.H.isscatter,'Load surfaces','subjects completed')
  if cscc.H.isscatter
    h = cat_surf_display(struct('data',cscc.files.surf{cscc.pos.x},'multisurf',1));
    h.Position(1:2) = cscc.H.graphics.Position .* [1.05 0.95];
  else
    h = cat_surf_display(struct('data',char(cscc.files.surf(unique([cscc.pos.x,cscc.pos.y]))),'multisurf',1));
  end
  spm_progress_bar('Clear');
  
  % give some feedback because loading take so long
return

%-----------------------------------------------------------------------
function checkvol(obj, event_obj)
%-----------------------------------------------------------------------
% Load the original image of selected files in SPM graphics window.
% Some further information or legend would be helpful.
%-----------------------------------------------------------------------
  global cscc st 
  
  spm_figure('Clear',cscc.H.graphics);
  spm_figure('Focus',cscc.H.graphics);
  spm_orthviews('Reset')
  [zl,rl] = spm_orthviews('ZoomMenu');
  if size(zl,1) > 1
    zl = zl';
    rl = rl';
  end
  if numel(zl)==8
    zl = [zl(1:end-2) 60 zl(end-1:end)];
    rl = [rl(1:end-2)  1 rl(end-1:end)];
  end
  spm_orthviews('ZoomMenu',zl,rl); 
  job.colormapc = flipud(cat_io_colormaps('BCGWHcheckcov'));
  job.prop  = 0.2; 
 
  %%
  cscc.H.multi = 1; 
  if cscc.H.isscatter 
    xeqy = 0;
    
    if cscc.H.multi
      id = mygetCursorInfo';
      
      spm_check_registration(char(unique(cscc.files.org(id))));
      
      
      if exist(cscc.files.p0{cscc.pos.x},'file')
        spm_orthviews('addtruecolourimage',1,cscc.files.p0{cscc.pos.x},...
          job.colormapc,job.prop,0,5);
      end
      
      vx_vol    = sqrt(sum(st.vols{1}.mat(1:3,1:3).^2));
      Ysrc      = cat_vol_resize(st.vols{1}.private.dat,'reduceV',vx_vol,vx_vol * 2,2,'meanm');
      [Ysrc,th] = cat_stat_histth(Ysrc,99);
      spm_orthviews('window',1,th + [ 0.1*-diff(th) 0.3*diff(th)] );   
      
    else
      ppos = [0.02 0.01 0.96 0.98];
      tpos = [0.50 0.98 0.96 0.02];
    end
  else
    xeqy = cscc.pos.x == cscc.pos.y; 

    if xeqy
      ppos = [0.02 0.01 0.96 0.98];
      tpos = [0.50 0.98 0.96 0.02];
    else
      ppos = [0.02 0.545 0.96 0.48 ; 0.02 0.010 0.96 0.48];
      tpos = [0.50 0.980 0.96 0.02 ; 0.50 0.475 0.96 0.02];
    end
  end
     
      
  if ~cscc.H.isscatter || ~cscc.H.multi || xeqy
    
    hi1 = spm_orthviews('Image',spm_vol(cscc.files.org{cscc.pos.x}),ppos(1,:)); 
    spm_orthviews('AddContext',hi1);
    
    gax = axes('Visible','off','Position',[0 0 1 1]);
    text(gax,tpos(1,1),tpos(1,2),spm_str_manip(cscc.files.org{cscc.pos.x},'k100'),...
      'FontSize',cscc.display.FS(cscc.display.FSi+1),'Color',[0 0 0.8],'LineStyle','none',...
      'HorizontalAlignment','center');
    
    if exist(cscc.files.p0{cscc.pos.x},'file')
      spm_orthviews('addtruecolourimage',1,cscc.files.p0{cscc.pos.x},...
        job.colormapc,job.prop,0,5);
    end
  end
  
  if ~cscc.H.isscatter && ~xeqy
    
    hi2 = spm_orthviews('Image',spm_vol(cscc.files.org{cscc.pos.y}),ppos(2,:));
    spm_orthviews('AddContext',hi2);
  
    gax = axes('Visible','off','Position',[0 0 1 1]);
    text(gax,tpos(1,1),tpos(2,2),spm_str_manip(cscc.files.org{cscc.pos.y},'k100'),...
      'FontSize',cscc.display.FS(cscc.display.FSi+1),'Color',[0 0 0.8],'LineStyle','none',...
      'HorizontalAlignment','center');
    
    if exist(cscc.files.p0{cscc.pos.y},'file')
      spm_orthviews('addtruecolourimage',2,cscc.files.p0{cscc.pos.y},...
        job.colormapc,job.prop,0,5);
    end

  end
  spm_orthviews('Reposition',[0 0 0]); 
  spm_orthviews('Zoom',120)
  
  spm_orthviews('redraw');  
return

%-----------------------------------------------------------------------
function check_worst_data(obj, event_obj) 
%-----------------------------------------------------------------------
% Old check worst function. The spm_input could be replaced by an popup 
% window. A specification of the data range would rather than the x worst 
% images would be useful. 
%-----------------------------------------------------------------------

  global cscc bp 

  if isempty(spm_figure('FindWin','Interactive')), spm('createintwin'); end

  if isempty(bp)
    data  = cscc.data.V;
    name  = 'Mean correlation';
  else
    data  = bp.data;
    name  = bp.name; 
  end
  switch deblank(name)
    case 'Noise rating (NCR)',             name = 'NCR'; 
    case 'Bias rating (ICR)',              name = 'ICR'; 
    case 'Resoution rating (RES)',         name = 'RES'; 
    case 'Weighted overall image quality rating (IQR)', name = 'IQR';
    case 'Euler number',                   name = 'Euler';
    case 'Size of topology defects (TDS)', name = 'TDS';
    case 'Mean correlation',               name = 'MNC'; 
    case 'Protocol-based IQR (PIQR)',      name = 'PIQR'; 
    case 'Mahalanobis distance',           name = 'MD';
    case 'Mahalanobis distance (IQR)',     name = 'MD IQR';
    case 'Mahalanobis distance (PIQR)',    name = 'MD PIQR';
  end
  if ~cscc.H.showtrash.Value
    data = setdiff(data,cscc.select.trashlist);
  end
  [sdata,order] = sort(data,'descend');
 
  
  n = numel(data);
  number = min([n 24]);
  number = spm_input('How many files?',1,'e',number);
  number = max([number 1]);
  number = min([number 24]);
  number = min([number n]);
  
  sdata = sdata(n:-1:1);
  list  = cellstr(char(cscc.data.V(order(n:-1:1)).fname));
  list2 = list(1:number);

  if cscc.H.mesh_detected
    % display single meshes and correct colorscale of colorbar
    for i=1:number
      h = cat_surf_render('Disp',deblank(list2(i,:)));

      % shift each figure slightly
      if i==1
          pos = get(h.figure,'Position');
      else
          pos = pos - [20 20 0 0];
      end

      % remove menubar and toolbar, use filename as title
      set(h.figure,'MenuBar','none','Toolbar','none',...
        'Name',spm_file(list2{i},'short50'),...
        'NumberTitle','off','Position',pos);
      cat_surf_render('ColourMap',h.axis,jet);
      cat_surf_render('ColourBar',h.axis,'on');
      cat_surf_render('CLim',h,[mn_data mx_data]);
    end
  else
    % break the filename into smaller pieces to display it close to the figure 
    [fnames0,fnames] = spm_str_manip(list2,'C'); 
    if isempty(fnames)
      fnames.m = list2;
      fnames.s = '';
    end
    for fi=1:numel(fnames.m)
      if ~isempty(fnames.s) && numel(spm_str_manip(fnames.m{fi},'H'))>1
        fnames.m{fi} = ['.' filesep fnames.m{fi}]; 
      end
      pp = spm_str_manip(fnames.m{fi},'Hl199'); pp(end+1) = filesep;
      ff = spm_str_manip(fnames.m{fi},'t');

      ssep = round(70 ./ sqrt(numel(fnames.m))); 
      lpp  = round( ceil( numel(pp)/(ssep/2) ) * (ssep/2) ); 
      lff  = round( ceil( numel(ff)/(ssep/2) ) * (ssep/2) );
      lpp  = lpp:-ssep:1; lpp(lpp>numel(pp)) = []; 
      lff  = lff:-ssep:1; lff(lff>numel(ff)) = []; 
      for pi = lpp, pp = sprintf('%s\n%s',pp(1:pi),pp(pi+1:end)); end
      for pi = lff, ff = sprintf('%s\n%s',ff(1:pi),ff(pi+1:end)); end
      
      % some values 
      val = sprintf('%s=%0.2f',name,sdata(fi)); 
      
      if numel(pp)>7
        list3{fi} = sprintf('%s\n%s\n%s',val,pp,ff);
      else
        list3{fi} = sprintf('%s\n%s%s',val,pp,ff);
      end
    end
    
    %
    spm_check_registration(char(list2));
    
    % add caption
    for fi=1:numel(fnames.m);
      spm_orthviews('Caption',fi,list3{fi},...
        'FontSize',cscc.display.FS(cscc.display.FSi-1));
    end
    spm_orthviews('Resolution',0.2);
    set(cscc.H.boxp,'Visible','on');
  end
  
  figure(cscc.H.figure);
return

%-----------------------------------------------------------------------
function checkxml(obj, event_obj)
%-----------------------------------------------------------------------
% Load XML report in SPM graphics window (see also checklog).
% This is just the first fast version of this function. 
% Finally, I want to use the xml structure from the file to print some
% specific information similar to the CAT report in cat_main. 
%-----------------------------------------------------------------------
  global cscc

  % visdiff(cscc.files.xml{cscc.pos.x}, cscc.files.xml{cscc.pos.y},'text')  
  
  spm_figure('Clear',cscc.H.graphics); 
  spm_figure('Focus',cscc.H.graphics);
  axis off;
  
  if cscc.H.isscatter || (cscc.pos.x == cscc.pos.y)
    textbox = [0 0 1 1];
    files   = cscc.files.xml(cscc.pos.x); 
  else
    textbox = [0 0.5 1 0.5; 0 0 1 0.5];
    files   = cscc.files.xml([cscc.pos.y,cscc.pos.x]); 
  end
  
  % avoid some long useless text passages
  badtacks = {'software>','catlog>','H.data>','LAB>'};
  badmode = 0; bdid = badtacks;
  for fi=1:numel(files);
    fid = fopen(files{fi});
    ph  = uipanel(cscc.H.graphics,'Units','normalized','position',textbox(fi,:), ...
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
% Function to control the visibility of the samples in the main plots.
% See also show_protocol.
%-----------------------------------------------------------------------

  global cscc
  
  % set entry of the GUI element
  if exist('obj','var')
    if obj>0
      set(cscc.H.samp,'Value',obj + 2);
      cscc.select.samp1d(:) = cscc.datagroups.sample==obj; 
    else
      set(cscc.H.samp,'Value',1); 
      cscc.select.samp1d(:) = true(size(cscc.select.samp1d)); 
    end
  else 
    obj = get(cscc.H.samp,'Value') - 2;
    if obj>0
      cscc.select.samp1d(:) = cscc.datagroups.sample==obj; 
    else
      cscc.select.samp1d(:) = true(size(cscc.select.samp1d)); 
    end
  end
  
  cscc.select.samp2d = (single(cscc.select.samp1d) * single(cscc.select.samp1d'))>0; 
  
  groups        = unique(cscc.datagroups.sample);
  symbols       = repmat('.',1,numel(groups));  % default symbol
  symbols(1:11) = 'o+^v<>ph*sd'; 
  
  % update scatter
  if cscc.H.isscatter
    if 0 %obj>0
      % find all (also deleted) scatter points of the active scatter plot
      sc.pos = [
        findobj(cscc.H.scat(cscc.H.scata),'type','scatter','marker',symbols(obj)); 
        findobj(cscc.H.scat(cscc.H.scata),'type','scatter','sizedatasource',symbols(obj))];
    else
      sc.pos = findobj(cscc.H.scat(cscc.H.scata),'type','scatter');
    end
    sc.posn = findobj(cscc.H.scat(cscc.H.scata),'type','scatter'); 
    
    % remove legend objects
    sc.pos  = setdiff(sc.pos ,[cscc.H.sclegend{cscc.H.scata}]); 
    sc.posn = setdiff(sc.posn,[cscc.H.sclegend{cscc.H.scata}]); 
    
    % remove objects without ID
    sc.pos(  cellfun('isempty',{sc.pos.UserData})  ) = [];      
    sc.posn( cellfun('isempty',{sc.posn.UserData}) ) = [];  
    
    %% remove inactive protocols / samples
    [pset,pseti] = setdiff( [sc.pos.UserData] , find( cscc.select.prot1d==0 ) ); sc.pos = sc.pos( pseti ); 
    [pset,pseti] = setdiff( [sc.pos.UserData] , find( cscc.select.samp1d==0 ) ); sc.pos = sc.pos( pseti ); 
    
    % 
    sc.posn = setdiff(sc.posn,sc.pos); 

    % remove trashed objects
    if ~cscc.H.showtrash.Value && ~isempty(sc.pos)
      sc.pos(strcmp({sc.pos.ZDataSource},'trash'))=[];
    end
    
    set(sc.pos ,'Visible','on');
    set(sc.posn,'Visible','off');
    
    
    for linei=1:numel(cscc.H.corrline{cscc.H.scata})
      indxy = get(cscc.H.corrline{cscc.H.scata}(linei),'UserData');
      if any(indxy(1) == [sc.pos.UserData]) && any(indxy(2) == [sc.pos.UserData])
        set(cscc.H.corrline{cscc.H.scata}(linei),'Visible','on');
      else
        set(cscc.H.corrline{cscc.H.scata}(linei),'Visible','off');
      end
    end
    
    xlim(cscc.H.scat(cscc.H.scata),'auto'); 
    xticks = get(cscc.H.scat(cscc.H.scata),'xtick'); 
    xlim([xticks(1) - diff(xticks(1:2)/2) xticks(end) + diff(xticks(1:2)/2)]);  

    ylim(cscc.H.scat(cscc.H.scata),'auto');
    yticks = get(cscc.H.scat(cscc.H.scata),'ytick'); 
    ylim([yticks(1) - diff(yticks(1:2)/2) yticks(end) + diff(yticks(1:2)/2)]);  

    zoom(cscc.H.scat(cscc.H.scata),'reset');
  else
    update_matrix
  end
  
  cscc.H.dcm.removeAllDataCursors
return

%-----------------------------------------------------------------------
function show_protocol(obj,event_obj)
%-----------------------------------------------------------------------
% Function to control the visibility of protocols in the main plots.
% See also show_sample.
%-----------------------------------------------------------------------

  global cscc
  
  % set GUI element entry
  if exist('obj','var')
    if obj>0
      set(cscc.H.prot,'Value',obj + 2);
      cscc.select.prot1d(:) = cscc.datagroups.protocol==obj; 
    else
      set(cscc.H.prot,'Value',1);
      cscc.select.prot1d(:) = true(size(cscc.select.prot1d)); 
    end
  else 
    obj = get(cscc.H.prot,'Value') - 2;
    if obj>0
      cscc.select.prot1d(:) = cscc.datagroups.protocol==obj; 
    else
      cscc.select.prot1d(:) = true(size(cscc.select.prot1d)); 
    end
  end
  
  cscc.select.prot2d = (single(cscc.select.prot1d) * single(cscc.select.prot1d'))>0; 
  
  % update scatter
  if cscc.H.isscatter
    if obj>0
      % find all (also deleted) scatter points of the active scatter plot
      sc.pos = findobj(cscc.H.scat(cscc.H.scata),'type','scatter','ZDataSource',num2str(obj,'%d')); 
    else
      sc.pos = findobj(cscc.H.scat(cscc.H.scata),'type','scatter');
    end
    sc.posn = findobj(cscc.H.scat(cscc.H.scata),'type','scatter');
    
    % remove legend objects
    sc.pos  = setdiff(sc.pos ,[cscc.H.sclegend{cscc.H.scata}]); 
    sc.posn = setdiff(sc.posn,[cscc.H.sclegend{cscc.H.scata}]); 

    sc.posn( cellfun('isempty',{sc.posn.UserData}) ) = [];   
    if ~isempty(sc.pos)
      % remove objects without ID
      sc.pos(  cellfun('isempty',{sc.pos.UserData}) )  = [];      
      
      % remove inactive protocols / samples
      [pset,pseti] = setdiff( [sc.pos.UserData] , find( cscc.select.prot1d==0 )); sc.pos = sc.pos( pseti ); 
    
      [pset,pseti] = setdiff( [sc.pos.UserData] , find( cscc.select.samp1d==0 )); sc.pos = sc.pos( pseti ); 
    end    
    
    sc.posn = setdiff(sc.posn,sc.pos); 

    if ~cscc.H.showtrash.Value && ~isempty(sc.pos)
      sc.pos(strcmp({sc.pos.ZDataSource},'trash'))=[];
    end
    
    set(sc.pos ,'Visible','on');
    set(sc.posn,'Visible','off');

    % remove lines if not both points are visible
    for linei=1:numel(cscc.H.corrline{cscc.H.scata})
      indxy = get(cscc.H.corrline{cscc.H.scata}(linei),'UserData');
      
      if ~isempty(sc.pos) && ...
        (any(indxy(1) == [sc.pos.UserData]) && any(indxy(2) == [sc.pos.UserData]))
        set(cscc.H.corrline{cscc.H.scata}(linei),'Visible','on');
      else
        set(cscc.H.corrline{cscc.H.scata}(linei),'Visible','off');
      end
    end

    
    xlim(cscc.H.scat(cscc.H.scata),'auto'); 
    xticks = get(cscc.H.scat(cscc.H.scata),'xtick'); 
    xlim([xticks(1) - diff(xticks(1:2)/2) xticks(end) + diff(xticks(1:2)/2)]);  

    ylim(cscc.H.scat(cscc.H.scata),'auto');
    yticks = get(cscc.H.scat(cscc.H.scata),'ytick'); 
    ylim([yticks(1) - diff(yticks(1:2)/2) yticks(end) + diff(yticks(1:2)/2)]);  

    zoom(cscc.H.scat(cscc.H.scata),'reset');
  else
    update_matrix
  end
  
  cscc.H.dcm.removeAllDataCursors
return

%-----------------------------------------------------------------------
function checklog(obj, event_obj)
%-----------------------------------------------------------------------
% Load the log-file from cat_main of the selected subjects into the SPM
% graphics window.
%-----------------------------------------------------------------------
  global cscc
  
  spm_figure('Clear',cscc.H.graphics); 
  spm_figure('Focus',cscc.H.graphics);
  axis off;
  
  if cscc.H.isscatter || (cscc.pos.x == cscc.pos.y)
    textbox = [0 0 1 1];
    files   = cscc.files.log(cscc.pos.x); 
  else
    textbox = [0 0.5 1 0.5; 0 0 1 0.5];
    files   = cscc.files.log([cscc.pos.x,cscc.pos.y]); 
  end
  
  for fi=1:numel(files); 
    fid = fopen(files{fi});
    ph  = uipanel(cscc.H.graphics,'Units','normalized','position',textbox(fi,:), ...
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
function checkbox_showtrash(obj, event_obj)
%-----------------------------------------------------------------------
  global oldx oldy 
  global cscc

  showtrash = get(cscc.H.showtrash,'Value');
  
  if ~showtrash
    if isfield(cscc.pos,'x'), oldx = cscc.pos.x; else oldx = 0; end
    if isfield(cscc.pos,'y'), oldy = cscc.pos.y; else oldy = 0; end
  else
    if oldx>0, cscc.pos.x = oldx; end
    if oldy>0, cscc.pos.y = oldy; end
  end
  cscc.H.dcm.removeAllDataCursors;
  set(cscc.H.slice,'visible','off');  
 
  if ~showtrash 
    if ~cscc.H.isscatter
      if isfield(cscc.pos,'x'), cscc.pos = rmfield(cscc.pos,'x'); end
      if isfield(cscc.pos,'y'), cscc.pos = rmfield(cscc.pos,'y'); end
      cscc.H.dcm.removeAllDataCursors;
    end
  end
  
  if cscc.H.isscatter
    sc.pos = findobj(cscc.H.scat(cscc.H.scata),'type','scatter','ZDataSource','trash'); 

    if ~showtrash
      set(sc.pos,'visible','off');
    else
      set(sc.pos,'visible','on');
    end
    
    % remove lines if not both points are visible
    for linei=1:numel(cscc.H.corrline{cscc.H.scata})
      indxy = get(cscc.H.corrline{cscc.H.scata}(linei),'UserData');
      
      if ~isempty(sc.pos) && ~showtrash && ...
        (any(indxy(1) == [sc.pos.UserData]) || any(indxy(2) == [sc.pos.UserData]))
        set(cscc.H.corrline{cscc.H.scata}(linei),'Visible','off');
      else
        set(cscc.H.corrline{cscc.H.scata}(linei),'Visible','on');
      end
    end

    
    xlim(cscc.H.scat(cscc.H.scata),'auto'); 
    xticks = get(cscc.H.scat(cscc.H.scata),'xtick'); 
    xlim([xticks(1) - diff(xticks(1:2)/2) xticks(end) + diff(xticks(1:2)/2)]);  

    ylim(cscc.H.scat(cscc.H.scata),'auto');
    yticks = get(cscc.H.scat(cscc.H.scata),'ytick'); 
    ylim([yticks(1) - diff(yticks(1:2)/2) yticks(end) + diff(yticks(1:2)/2)]);  

    zoom(cscc.H.scat(cscc.H.scata),'reset');
  else
    update_matrix;
  end
  
  %{
  if showtrash
    %createDatatip(cscc.H.dcm,get(cscc.H.corr,'children'),[oldx oldy]);
  else
    %createDatatip(cscc.H.dcm,get(cscc.H.corr,'children'),[oldx oldy]);
  end
  %}
return  

%{
%-----------------------------------------------------------------------
function checkbox_cbarfix(obj, event_obj)
%-----------------------------------------------------------------------
  global cscc
  
  if cscc.H.isscatter
    % do something
  else
    update_matrix;
  end
  
return  
%}

%-----------------------------------------------------------------------
function show_mahalanobis(X,MD,scata)
%-----------------------------------------------------------------------
  global cscc
    
  if ~exist('scata','var')
    if isfield(cscc.H,'scata')
      scata = cscc.H.scata; 
    else
      scata = 1; cscc.H.scata = scata;
    end
  else
    cscc.H.scata = scata;
  end
  axis(cscc.H.scat(cscc.H.scata));
  

  
  if ~cscc.H.isscatter
    set([cscc.H.corr,get(cscc.H.corr,'children')],'visible','off','HitTest','off','Interruptible','off')
    if isfield(cscc.H,'cbarfix')   && ishandle(cscc.H.cbarfix),   set(cscc.H.cbarfix  ,'enable' ,'off'); end
    if isfield(cscc.H,'showtrash') && ishandle(cscc.H.showtrash), set(cscc.H.showtrash,'Value' ,1); end
    
    % remove sliders and text
    if isfield(cscc.pos,'tar_mouse')
      delete(cscc.pos.tar_mouse); 
      cscc.pos = rmfield(cscc.pos,'tar_mouse');
      set(cscc.H.alphabox,'Visible','off');

      if isfield(cscc.pos,'x'), cscc.pos = rmfield(cscc.pos,'x'); end
      if isfield(cscc.pos,'y'), cscc.pos = rmfield(cscc.pos,'y'); end

      if isfield(cscc.H,'slice')    && ishandle(cscc.H.slice)
        set(cscc.H.slice,'Visible','off'); cla(cscc.H.slice); 
      end
      if isfield(cscc.H,'sslider')  && ishandle(cscc.H.sslider)
        set(cscc.H.sslider,'Visible','off');
      end
      if isfield(cscc.H,'alphabox') && ishandle(cscc.H.alphabox)
        set(cscc.H.alphabox,'Visible','off'); 
      end
      if isfield(cscc.H,'mm') && ishandle(cscc.H.mm),
        set([cscc.H.mm,cscc.H.mm_txt],'Visible','off');
      end
      
      if ~cscc.H.mesh_detected
        set([cscc.H.mm,cscc.H.mm_txt],'Visible','off');
      end
    end
      
    set([cscc.H.trashui.trash,cscc.H.trashui.detrash,...
         cscc.H.trashui.trashcol,cscc.H.trashui.detrashcol,...
         cscc.H.trashui.trashrow,cscc.H.trashui.detrashrow],'Visible','off','Enable','off');
    set([cscc.H.trashui.trash,cscc.H.trashui.detrash],'Visible','on');

    unit = struct2cell(cscc.H.checkui); set([unit{:}],'Enable','off');  

  end
  
  for i=setdiff(1:numel(cscc.H.scat),cscc.H.scata)
    set([cscc.H.scat(i);get(cscc.H.scat(i),'Children')],...
      'visible','off','HitTest','off','Interruptible','off');
  end
  set([cscc.H.scat(cscc.H.scata);get(cscc.H.scat(cscc.H.scata),'children')],...
    'visible','on','HitTest','on','Interruptible','on');
  try set([cscc.H.sclegend{:}],'visible','off'); end
  set(findobj('type','Legend'),'visible','on');

  if isempty(get(cscc.H.scat(cscc.H.scata),'children'))
    % get very similar scans
    cscc.data.YpY_tmp = cscc.data.YpY - tril(cscc.data.YpY);
    [indx, indy] = find(cscc.data.YpY_tmp>0.925);

    groups        = unique(cscc.datagroups.sample);
    symbols       = repmat('.',1,numel(groups));  % default symbol
    symbols(1:11) = 'o+^v<>ph*sd';                % need x for unset

    %% Create legend by creation of hidden dummy objects 
    %  display first object for the legend
    hold(cscc.H.scat(cscc.H.scata),'on'); 
    for gi=1:numel(groups)
      txt{gi} = sprintf('sample %d \n',gi);
      Xt = cscc.data.X(cscc.datagroups.sample==groups(gi),:); 
      cscc.H.sclegend{scata}(gi) = scatter(cscc.H.scat(cscc.H.scata),...
        Xt(1,1),Xt(1,2),30,[0 0 0],symbols(gi),'Linewidth',2);
      if ~isempty( strfind('osd^v<>ph', symbols(gi) ) )
        set(cscc.H.sclegend{scata}(gi),'MarkerFaceColor','flat','markerFaceAlpha',1/3);
      end
    end
    txt{end+1} = 'excluded';
    cscc.H.sclegend{scata}(gi+1) = scatter(cscc.H.scat(cscc.H.scata),...
      Xt(1,1),Xt(1,2),30,[1 0 0],'x','Linewidth',2,'Visible','off');
    if numel(indx)/size(cscc.data.YpY,1)<0.5 && numel(indx)>0
      txt{end+1}   = 'highly corr. scans'; 
      cscc.H.sclegend{scata}(gi+2) = plot(cscc.H.scat(cscc.H.scata),...
        [cscc.data.X(indx(1),1);cscc.data.X(indy(1),1)],...
        [cscc.data.X(indx(1),2);cscc.data.X(indy(1),2)],'Color',[0 0 0],'Linewidth',2);
    end
    % create legend
    hl = legend(cscc.H.scat(cscc.H.scata),txt,'location','southwest');
    set(get(hl,'title'),'string','Legend');
    hold(cscc.H.scat(cscc.H.scata),'off'); 
    set([cscc.H.sclegend{scata}],'visible','off');

    
    %%
    %cscc.data.X(isnan(cscc.data.X)) = 0;  
    %S  = cov(cscc.data.X); 
   % mu = mean(cscc.data.X);
   % cscc.data.MD = (cscc.data.X-repmat(mu,[size(cscc.data.X,1),1]))*inv(S)*(cscc.data.X-repmat(mu,[numel(cscc.data.X),1]))';
   % cscc.data.MD = diag(cscc.data.MD);
   
    % because we use a splitted colormap we have to set the color values explicitely
    MDs  = min(63,max(0,64*MD/10)); %(round(max(cscc.data.MD)/6)*6);
    C    = zeros(numel(MD),3);
    cmap = [jet(64); gray(64)];
    for i=1:numel(MD)
      C(i,:) = cmap(round(MDs(i) )+1,:);
    end

    hold(cscc.H.scat(cscc.H.scata),'on')

    % plot lines between similar objects
    if numel(indx)/size(cscc.data.YpY,1)<0.5 && numel(indx)
      for i=1:numel(indx)
        cscc.H.corrline{scata}(i) = plot( cscc.H.scat(cscc.H.scata),...
          [X(indx(i),1);X(indy(i),1)],[X(indx(i),2);X(indy(i),2)],'-','Color',...
          repmat(cscc.display.figcolor(1)*0.9 - cscc.display.figcolor(1)*0.9 * ...
            ((cscc.data.YpY_tmp(indx(i),indy(i)) - 0.925)/0.0725),1,3),...
          'LineWidth',2,'HitTest','off','Interruptible','off');
        set(cscc.H.corrline{scata}(i),'UserData',[indx(i) indy(i)]);
      end
    else
      cscc.H.corrline{scata} = plot([]); 
    end

    % plot data entries
    I   = 1:size(X,1);
    for gi=1:numel(groups)
      It = I(cscc.datagroups.sample==groups(gi)); 
      Xt = X(cscc.datagroups.sample==groups(gi),:); 
      Ct = C(cscc.datagroups.sample==groups(gi),:);
      Pt = cscc.datagroups.protocol(cscc.datagroups.sample==groups(gi),:);
      cscc.H.sc{scata} = cell(size(Xt,1),1);
      for sci=1:size(Xt,1)
        cscc.H.sc{scata}{gi}{sci} = scatter( cscc.H.scat(cscc.H.scata), ...
          Xt(sci,1), ...
          Xt(sci,2), ...
          30,...
          Ct(sci,:),...
          symbols(gi), ...
          'ZDataSource',num2str(Pt(sci),'%d'),...
          'UserData',It(sci),...
          'Linewidth',2);
        if ~isempty( strfind('osd^v<>ph', symbols(gi) ) )
          set(cscc.H.sc{scata}{gi}{sci},'MarkerFaceColor','flat','markerFaceAlpha',1/3);
        end
        if any(It(sci)==cscc.select.trashlist)
          set(cscc.H.sc{scata}{gi}{sci},'sizedatasource',get(cscc.H.sc{scata}{gi}{sci},'marker'),...
          'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
          'ZDataSource','trash','MarkerFaceAlpha',0);
        end
      end
    end

    hold(cscc.H.scat(cscc.H.scata),'off')
  end

  xlabel(cscc.H.scat(cscc.H.scata),...
    '<----- Worst ---      Mean correlation       --- Best ------>  ',...
      'FontSize',cscc.display.FS(cscc.display.FSi),'FontWeight','Bold');
  ylabel(cscc.H.scat(cscc.H.scata),...
    '<----- Worst ---      Weighted overall image quality rating     --- Best ------>  ',...
    'FontSize',cscc.display.FS(cscc.display.FSi),'FontWeight','Bold');
  title(cscc.H.scat(cscc.H.scata),...
    '<--- Smallest -- Mahalanobis distance (Color) -- Largest ---->  ','FontSize',...
    cscc.display.FS(cscc.display.FSi+2),'FontWeight','Bold');


  % update colorbar 
  cticks = 7;
  mn     = 0; 
  mx     = 10; %round(max(cscc.data.MD)/6)*6;
  ticks  = linspace(mn,mx,cticks);
  set(cscc.H.cbar,'XTick',1:63/(numel(ticks)-1):64,...
    'XTickLabel',cellstr(num2str(ticks','%0.2f'))); 
  %round(100*linspace(min(cscc.data.YpYt(:)),max(cscc.data.YpYt(:)),5))/100);
  %caxis(cscc.H.scat(cscc.H.scata),[mn mx])
  
  cscc.H.isscatter = 1;
  
  show_sample
  %show_protocol

  xlim(cscc.H.scat(cscc.H.scata),'auto'); 
  xticks = get(cscc.H.scat(cscc.H.scata),'xtick'); 
  xlim([xticks(1) - diff(xticks(1:2)/2) xticks(end) + diff(xticks(1:2)/2)]);  

  ylim(cscc.H.scat(cscc.H.scata),'auto');
  yticks = get(cscc.H.scat(cscc.H.scata),'ytick'); 
  ylim([yticks(1) - diff(yticks(1:2)/2) yticks(end) + diff(yticks(1:2)/2)]);  
  
  zoom reset
return

%-----------------------------------------------------------------------
function update_matrix(data, order)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
  global cscc

  showtrash = get(cscc.H.showtrash,'Value'); 
 
  if ~exist('order','var')
    order = cscc.H.sorted; 
  else
    cscc.H.sorted = order; 
  end
  
  if ~exist('data' ,'var')
    if order
      data = cscc.data.YpY(cscc.data.ind_sorted,cscc.data.ind_sorted);
    else
      data = cscc.data.YpY;
    end
  end

  if order 
    cscc.select.trash1dt = cscc.select.trash1d(cscc.data.ind_sorted); 
    cscc.select.trash2dt = cscc.select.trash2d(cscc.data.ind_sorted,cscc.data.ind_sorted);
    
    cscc.select.samp1dt  = cscc.select.samp1d(cscc.data.ind_sorted);
    cscc.select.samp2dt  = cscc.select.samp2d(cscc.data.ind_sorted,cscc.data.ind_sorted);
    cscc.select.prot1dt  = cscc.select.prot1d(cscc.data.ind_sorted); 
    cscc.select.prot2dt  = cscc.select.prot2d(cscc.data.ind_sorted,cscc.data.ind_sorted);    
  else
    cscc.select.trash1dt = cscc.select.trash1d;
    cscc.select.trash2dt = cscc.select.trash2d;
    
    cscc.select.samp1dt  = cscc.select.samp1d;
    cscc.select.samp2dt  = cscc.select.samp2d;
    cscc.select.prot1dt  = cscc.select.prot1d; 
    cscc.select.prot2dt  = cscc.select.prot2d;    
  end
  
  
  if ~showtrash
    data = reshape(data(cscc.select.trash2dt & cscc.select.samp2dt & cscc.select.prot2dt),...
      sum(cscc.select.trash1dt .* cscc.select.samp1dt .* cscc.select.prot1dt),...
      sum(cscc.select.trash1dt .* cscc.select.samp1dt .* cscc.select.prot1dt)); 
    cscc.select.trash2dt = reshape(cscc.select.trash2dt( ...
      cscc.select.trash2dt & cscc.select.samp2dt & cscc.select.prot2dt),...
      sum(cscc.select.trash1dt .* cscc.select.samp1dt .* cscc.select.prot1dt),...
      sum(cscc.select.trash1dt .* cscc.select.samp1dt .* cscc.select.prot1dt)); 
  else
    data = reshape(data(cscc.select.samp2dt & cscc.select.prot2dt),...
      sum(cscc.select.samp1dt .* cscc.select.prot1dt),sum(cscc.select.samp1dt .* cscc.select.prot1dt)); 
    cscc.select.trash2dt = reshape(cscc.select.trash2dt(cscc.select.samp2dt & cscc.select.prot2dt),...
      sum(cscc.select.samp1dt .* cscc.select.prot1dt),sum(cscc.select.samp1dt .* cscc.select.prot1dt)); 
  end
  
   
  %% scale data 
  cticks  = 7;
  ltick   = cticks - 1; 
  htick   = ltick/2;
  if cscc.H.cbarfix.Value==1
    mx = 1.0;
    if cscc.H.mesh_detected
      mn = 0.8; 
    else
      mn = 0.7; 
    end
  elseif cscc.H.cbarfix.Value==2
    mx = 1.0;
    md = median(cscc.data.mean_cov); 
    mn = mx - ltick * ceil((mx - md)/htick * 100)/100; 
  else
    md = round(median(cscc.data.mean_cov)*100)/100; 
    dd = min( round( (1 - md) / htick * 100) /100 * htick, ...
              htick * round(2 * std(cscc.data.mean_cov) * 100 ) / 100 ); 
    dd = round(dd * 100 ) / 100; 
    mx = md + dd;
    mn = md - dd; 
  end
  ticks = linspace(mn,mx,cticks);
 
  data_scaled = min(1,max(0,(data - mn)/(mx - mn)));

  % create image if not exist
  if isempty(get(cscc.H.corr,'children'))  
    image(cscc.H.corr,data_scaled);
  else
    if showtrash
      mylim = 0.5 + [0 size(cscc.data.YpY,1)] - [0 numel(cscc.data.mean_cov) - ...
        sum(cscc.select.samp1dt .* cscc.select.prot1dt)]; 
    else
      mylim = (0.5 + [0 size(cscc.data.YpY,1)]) - [0 numel(cscc.data.mean_cov) - ...
        sum(cscc.select.trash1dt .* cscc.select.samp1dt .* cscc.select.prot1dt)];
    end
    if diff(mylim)==0
      % deselect all images
      set(get(cscc.H.corr,'children'),'Cdata',65);
      xlim(cscc.H.corr,[0.5 1.5]);
      ylim(cscc.H.corr,[0.5 1.5]);
      axis off
      return
    end
    xlim(cscc.H.corr,mylim);
    ylim(cscc.H.corr,mylim);
  end
  
  % show only lower left triangle
  bg = 119;
  if showtrash
    set(get(cscc.H.corr,'children'),'Cdata',64 * (data_scaled + 65/64*(1-cscc.select.trash2dt)) .* tril(data>0) + bg*(~tril(data>0)) );
  else
    set(get(cscc.H.corr,'children'),'Cdata',64 * (data_scaled) .* tril(data>0) + bg*(~tril(data>0)) );
  end
  
  % update colorbar 
  set(cscc.H.cbar,'XTick',1:63/(numel(ticks)-1):64,...
    'XTickLabel',cellstr(num2str(ticks','%0.2f'))); 
  %round(100*linspace(min(cscc.data.YpYt(:)),max(cscc.data.YpYt(:)),5))/100);

  % update axis limits
  if 0 && ~cscc.H.isscatter
    if showtrash
      mylim = 0.5 + [0 size(cscc.data.YpY,1)]; 
    else
      mylim = (0.5 + [0 size(cscc.data.YpY,1)]) - [0 numel(cscc.select.trashlist)];
    end
    xlim(cscc.H.corr,mylim);
    ylim(cscc.H.corr,mylim);
  end
  
  if isempty(cscc.H.dcm.getCursorInfo)
    unit = struct2cell(cscc.H.checkui); 
    set([unit{cellfun(@ishandle,unit)}],'Enable','off');  
    set([cscc.H.trashui.trash;cscc.H.trashui.detrash; ...
         cscc.H.trashui.trashcol;cscc.H.trashui.detrashcol; ...
         cscc.H.trashui.trashrow;cscc.H.trashui.detrashrow],'Enable','off');
  end
return

%-----------------------------------------------------------------------
function show_matrix(data, order)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
  global cscc
 
  set([cscc.H.corr,get(cscc.H.corr,'children')],'visible','on','HitTest','on','Interruptible','on');
  set(findobj('type','Legend'),'visible','off','HitTest','off','Interruptible','off');
  try
    for i=1:numel(cscc.H.scat), set([cscc.H.scat(i);get(cscc.H.scat(i),'Children')],...
        'visible','off','HitTest','off','Interruptible','off'); end; 
  end
  axis(cscc.H.corr);
  if isfield(cscc.H,'cbarfix')   && ishandle(cscc.H.cbarfix)
    set(cscc.H.cbarfix  ,'enable' ,'on'); 
  end
  if isfield(cscc.H,'showtrash') && ishandle(cscc.H.showtrash) 
    if isempty(cscc.select.trashlist)
      set(cscc.H.showtrash,'enable' ,'off'); 
    else
      set(cscc.H.showtrash,'enable' ,'on'); 
    end
  end
  try set([cscc.H.sclegend{:}],'visible','off'); end
 
  % == set GUI fields ==
  
  % update buttons and remove cursor, slices/surfaces? and text
  if isfield(cscc.pos,'tar_mouse')
    delete(cscc.pos.tar_mouse); 
    cscc.pos = rmfield(cscc.pos,'tar_mouse'); 
  
    if isfield(cscc.pos,'x'), cscc.pos = rmfield(cscc.pos,'x'); end
    if isfield(cscc.pos,'y'), cscc.pos = rmfield(cscc.pos,'y'); end

    if isfield(cscc.H,'slice')    && ishandle(cscc.H.slice),    set(cscc.H.slice        ,'Visible','off'); cla(cscc.H.slice); end
    if isfield(cscc.H,'sslider')  && ishandle(cscc.H.sslider),  set(cscc.H.sslider      ,'Visible','off'); end
    if isfield(cscc.H,'alphabox') && ishandle(cscc.H.alphabox), set(cscc.H.alphabox     ,'Visible','off'); end
    if isfield(cscc.H,'mm')       && ishandle(cscc.H.mm),       set([cscc.H.mm,cscc.H.mm_txt],'Visible','off'); end

    set(cscc.H.alphabox,'Visible','off');
  end
  
  set([cscc.H.trashui.trash,cscc.H.trashui.detrash,...
       cscc.H.trashui.trashcol,cscc.H.trashui.detrashcol,...
       cscc.H.trashui.trashrow,cscc.H.trashui.detrashrow],'Enable','off','visible','on');
  set([cscc.H.trashui.trash,cscc.H.trashui.detrash],'Visible','off'); 
  unit = struct2cell(cscc.H.checkui); set([unit{cellfun(@ishandle,unit)}],'Enable','off');  

  
  cscc.H.isscatter = 0;

  % get sorting order
  cscc.H.sorted = order;
  
  % update image
  update_matrix(data,order)

  % update title and label elements
  if cscc.H.sorted
    xlabel(cscc.H.corr,'<----- Best ---      File Order      --- Worst ------>  ',...
      'FontSize',cscc.display.FS(cscc.display.FSi),'FontWeight','Bold');
    ylabel(cscc.H.corr,'<----- Worst ---      File Order      --- Best ------>  ',...
      'FontSize',cscc.display.FS(cscc.display.FSi),'FontWeight','Bold');
    title(cscc.H.corr,'Sorted Sample Correlation Matrix  ','FontSize',...
      cscc.display.FS(cscc.display.FSi+2),'FontWeight','Bold');
  else
    xlabel(cscc.H.corr,'<----- First ---      File Order      --- Last ------>  ',...
      'FontSize',cscc.display.FS(cscc.display.FSi),'FontWeight','Bold');
    ylabel(cscc.H.corr,'<----- Last ---      File Order      --- First ------>  ',...
    'FontSize',cscc.display.FS(cscc.display.FSi),'FontWeight','Bold');
    title(cscc.H.corr,'Sample Correlation Matrix  ','FontSize',...
    cscc.display.FS(cscc.display.FSi+2),'FontWeight','Bold');
  end

  zoom reset
  return
  
%-----------------------------------------------------------------------
function checkbox_names(obj, event_obj)
%-----------------------------------------------------------------------
  global cscc
  cscc.H.show_name = get(cscc.H.chbox,'Value');
  show_boxplot;
return

%-----------------------------------------------------------------------
function checkbox_plot(obj, event_obj)
%-----------------------------------------------------------------------
  global cscc
  cscc.H.show_violin = get(cscc.H.plotbox,'Value');
  show_boxplot;
return

%-----------------------------------------------------------------------
function show_boxplot(data_boxp, name_boxp, quality_order, obj)
%-----------------------------------------------------------------------
  global cscc bp

  if nargin <3
    data_boxp     = bp.data;
    name_boxp     = bp.name;
    quality_order = bp.order;
  else
    bp.data       = data_boxp;
    bp.name       = name_boxp;
    bp.order      = quality_order;
  end

  if iscell(name_boxp), name_boxp = name_boxp{1}; end
  
  spm_figure('Clear',cscc.H.graphics); 
  spm_figure('Focus',cscc.H.graphics); 

  cscc.H.boxplot = axes('Position',cscc.pos.boxplot,'Parent',cscc.H.graphics);
  set(cscc.H.graphics,'Renderer','OpenGL','color',[0.95 0.95 0.95]);

  %if isfield(cscc.H,'chbox'), nval = ~get(cscc.H.chbox,'value'), else nval = 1; end
 
  cscc.H.refresh = uicontrol(cscc.H.graphics,...
    'Units','normalized','position',cscc.pos.refresh,'callback',@show_boxplot,...
    'string','REF','ForegroundColor',[0 0.8 0],'FontSize',cscc.display.FS(cscc.display.FSi),...
    'ToolTipString','Include row','Style','Pushbutton','Enable','on'); 
  buttonicon(cscc.H.refresh,'REF'  ,...
    fullfile(spm('dir'),'toolbox','cat12','html','icons','refresh.png'));
  
  switch deblank(name_boxp)
    case 'Noise rating (NCR)',             name = 'NCR'; 
    case 'Bias rating (ICR)',              name = 'ICR'; 
    case 'Resoution rating (RES)',         name = 'RES'; 
    case 'Weighted overall image quality rating (IQR)', name = 'IQR';
    case 'Euler number',                   name = 'Euler';
    case 'Size of topology defects (TDS)', name = 'TDS';
    case 'Mean correlation',               name = 'MNC'; 
    case 'Protocol-based IQR (PIQR)',      name = 'PIQR'; 
    case 'Mahalanobis distance',           name = 'MD';
    case 'Mahalanobis distance (IQR)',     name = 'MD-IQR';
    case 'Mahalanobis distance (PIQR)',    name = 'MD-PIQR';
    otherwise,                             name = ''; 
  end
  
  cscc.H.worst = uicontrol(cscc.H.graphics,...
    'Units','normalized','position',cscc.pos.worst,'Style','Pushbutton',...
    'HorizontalAlignment','center','callback',@check_worst_data, ..., bp.data, bp.name},...
    'string',['Check worst ' name],'ToolTipString','Display most deviating files',...
    'FontSize',cscc.display.FS(cscc.display.FSi),'ForegroundColor',[0.8 0 0]);

  cscc.H.chbox = uicontrol(cscc.H.graphics,...
    'string','Show filenames','Units','normalized',...
    'position',cscc.pos.fnamesbox,'callback',@checkbox_names,...
    'Style','CheckBox','HorizontalAlignment','center',...
    'ToolTipString','Show filenames in boxplot','value',cscc.H.show_name,...
    'Interruptible','on','Visible','on','FontSize',cscc.display.FS(cscc.display.FSi));

  cscc.datagroups.n_samples = max(cscc.datagroups.sample);

  xcscc.pos = cell(1,cscc.datagroups.n_samples);
  data = cell(1,cscc.datagroups.n_samples);

  %% create filenames
  hold on
  allow_violin = 1;
  for i=1:cscc.datagroups.n_samples
    indtype   = { cscc.select.trash1d' 'k.' [0 0 0] 10; ~cscc.select.trash1d' 'rx' [1 0 0] 3}; 
    gnames{i} = sprintf('S%d',i);
    for ii=1:size(indtype,1)
      ind  = find(cscc.datagroups.sample == i & indtype{ii,1});
      if numel(ind>0)
        datap{i} = data_boxp(ind);
        if ii==1, data{i} = datap{i}; end

				if length(ind)<10
					allow_violin = 0; 
          cscc.H.show_violin = 0;
				end
				
        if cscc.datagroups.n_samples == 1
          xcscc.pos{i} = (i-1)+2*(0:length(ind)-1)/(length(ind)-1);
        else
          xcscc.pos{i} = 0.5/length(ind) + 0.5+(i-1)+1*(0:length(ind)-1)/(length(ind));
        end

        if get(cscc.H.chbox,'value')
          for j=1:length(ind)
            cscc.H.fnames{j,i} = text(xcscc.pos{i}(j),datap{i}(j),...
              cscc.files.fname.m{ind(j)},'Color',indtype{ii,3},...
              'FontSize',cscc.display.FS(cscc.display.FSi-1),'HorizontalAlignment','center');
          end
        else
          for j=1:length(ind)
            cscc.H.fnames{j,i} = plot(xcscc.pos{i}(j),datap{i}(j),indtype{ii,2},'MarkerSize',indtype{ii,4});
          end
        end
      end
    end 
  end

% allow violin plot onl if samples are all large enough
if allow_violin
	cscc.H.plotbox = uicontrol(cscc.H.graphics,...
			'string','Violinplot','Units','normalized',...
			'position',cscc.pos.plotbox,'callback',@checkbox_plot,...
			'Style','CheckBox','HorizontalAlignment','center',...
			'ToolTipString','Switch to Violinplot','value',cscc.H.show_violin,...
			'Interruptible','on','Visible','on','FontSize',cscc.display.FS(cscc.display.FSi));
end

  %% create boxplot
  opt = struct('groupnum',0,'ygrid',1,'box',1,'violin',2*cscc.H.show_violin,'median',2,...
               'groupcolor',jet(cscc.datagroups.n_samples),'names',{gnames},...
               'xlim',[-.25 cscc.datagroups.n_samples+1.25]); 
  if max(data_boxp) > min(data_boxp)
    ylim_add = 0.075;
    yamp     = max(data_boxp) - min(data_boxp);
    ylim_min = min(data_boxp) - ylim_add * yamp;
    ylim_max = max(data_boxp) + ylim_add * yamp; 
    norm     = round(20/yamp);
    if norm==0
      norm = round(yamp/20); 
    end
    opt.ylim = round( [ylim_min ylim_max] * norm ) / norm;
  end    
  cat_plot_boxplot(data,opt); box on;


  %% add colored labels and title
  if cscc.datagroups.n_samples > 1
    [tmp,  tmp2] = spm_str_manip(char(cscc.files.fname.s),'C'); %spm_file(char(tmp2.s),'short25')
    if length(tmp)>80,tmp = sprintf('../%s',strrep(tmp,tmp2.s,'')); end
    title_str = sprintf('Common filename: %s*',tmp);
  else
    title_str = sprintf('Common filename: %s*',spm_file(char(cscc.files.fname.s),'short25'));
  end
  title({['Boxplot: ' name_boxp],title_str},'FontSize',cscc.display.FS(cscc.display.FSi+1),'FontWeight','Bold');
  xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',cscc.display.FS(10),...
      'FontWeight','Bold');

  xcscc.pos = -0.35 - cscc.datagroups.n_samples*0.1;

  if quality_order >= 2 
    % reverse order to have the good things allways on the top
    set(gca, 'YDir','reverse');
    quality_order = 1; 
    t = ylim_min; ylim_min = ylim_max; ylim_max = t; 
  end
  if (length(data_boxp) > 2)
    if quality_order > 0 
      text(xcscc.pos, ylim_min,'<----- Low rating (poor quality)  ','Color','red','Rotation',...
          90,'HorizontalAlignment','left','FontSize',cscc.display.FS(9),'FontWeight','Bold')
      text(xcscc.pos, ylim_max,'High rating (good quality) ------>  ','Color',[0 0.8 0],'Rotation',...
          90,'HorizontalAlignment','right','FontSize',cscc.display.FS(9),'FontWeight','Bold')
    elseif quality_order < 0 
			if strfind(name_boxp,'Mahalanobis')
				text(xpos, ylim_max,'Largest distance to sample ------>  ','Color','red','Rotation',...
						90,'HorizontalAlignment','right','FontSize',FS,'FontWeight','Bold')
				text(xpos, ylim_min,'<----- Smallest distance to sample  ','Color','green','Rotation',...
						90,'HorizontalAlignment','left','FontSize',FS,'FontWeight','Bold')
			else
				text(xpos, ylim_max,'Low rating (poor quality) ------>  ','Color','red','Rotation',...
						90,'HorizontalAlignment','right','FontSize',FS,'FontWeight','Bold')
				text(xpos, ylim_min,'<----- High rating (good quality)  ','Color','green','Rotation',...
						90,'HorizontalAlignment','left','FontSize',FS,'FontWeight','Bold')
			end
    end
    text(xcscc.pos, (ylim_max+ylim_min)/2,sprintf('%s',name_boxp),'Color','black','Rotation',...
          90,'HorizontalAlignment','center','FontSize',cscc.display.FS(9),'FontWeight','Bold')
  end

  hold off


return

%-----------------------------------------------------------------------
function update_alpha(obj, event_obj)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
  global cscc

  alphaval = get(cscc.H.alphabox,'Value');

  % display image with 2nd colorbar (gray)
  image(cscc.H.slice,65 + cscc.data.img);
  set(cscc.H.slice,'Position',cscc.pos.slice  .* [1 1 1 1 - 0.5*cscc.H.isscatter],'visible','off'); 
  
  if ~cscc.H.mesh_detected
    set(cscc.H.slice,'HitTest','off','Interruptible','off');
  end
  
  % prepare alpha overlays for red and green colors
  if alphaval > 0
    % get 2%/98% ranges of difference image
    range = cat_vol_iscaling(cscc.data.img_alpha(:),[0.02 0.98]);

    hold(cscc.H.slice,'on')
    if ~cscc.H.mesh_detected
      alpha_g = cat(3, zeros(size(cscc.data.img_alpha)), ...
        alphaval*ones(size(cscc.data.img_alpha)), zeros(size(cscc.data.img_alpha)));
      alpha_r = cat(3, alphaval*ones(size(cscc.data.img_alpha)), ...
        zeros(size(cscc.data.img_alpha)), zeros(size(cscc.data.img_alpha)));
    else
      alpha_r = cat(3, zeros(size(cscc.data.img_alpha)), ...
        alphaval*ones(size(cscc.data.img_alpha)), zeros(size(cscc.data.img_alpha)));
      alpha_g = cat(3, alphaval*ones(size(cscc.data.img_alpha)), ...
        zeros(size(cscc.data.img_alpha)), zeros(size(cscc.data.img_alpha)));
    end
    hg = image(cscc.H.slice,alpha_g); 
    set(hg, 'AlphaData', cscc.data.img_alpha .* (cscc.data.img_alpha>range(2))); 
    if ~cscc.H.mesh_detected, set(hg,'HitTest','off','Interruptible','off','AlphaDataMapping','scaled'); end
    hr = image(cscc.H.slice,alpha_r); 
    set(hr, 'AlphaData',-cscc.data.img_alpha .* (cscc.data.img_alpha<range(1))); 
    if ~cscc.H.mesh_detected, set(hr,'HitTest','off','Interruptible','off','AlphaDataMapping','scaled'); end
    hold(cscc.H.slice,'off')
  end

  if cscc.H.mesh_detected
    if alphaval > 0
      colormap(cscc.H.slice,[gray(64);gray(64)]); 
    else
      colormap(cscc.H.slice,[gray(64);jet(64)]);
    end  
  end
return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
  global cscc

  if isfield(cscc.H,'mm')
    slice_mm = get(cscc.H.mm,'Value');
  else
    slice_mm = 0;
  end

  if cscc.data.Vchanged_names
    P = cscc.data.Vchanged;
  else
    P = cscc.data.V;
  end

  vx   =  sqrt(sum(P(1).mat(1:3,1:3).^2));
  Orig = P(1).mat\[0 0 0 1]';
  sl   = round(slice_mm/vx(3)+Orig(3));

  % if slice is outside of image use middle slice
  if (sl>P(1).dim(3)) || (sl<1)
    sl = round(P(1).dim(3)/2);
  end
  set(cscc.H.mm,'Value',(sl-Orig(3))*vx(3));

  M  = spm_matrix([0 0 sl]);
  cscc.data.data_array_diff = cscc.data.data_array;

  %%
  for i = 1:length(cscc.data.V)
    cscc.data.img = spm_slice_vol(P(i),M,P(1).dim(1:2),[1 0]);
    cscc.data.img(isnan(cscc.data.img)) = 0;

    % rescue unscaled data
    cscc.data.data_array_diff(:,:,i) = cscc.data.img;
  
    if cscc.H.inorm %&& ~isempty(strfind(cscc.files.dataprefix,'wp'))
      cscc.data.imgscale = median(cscc.data.img(cscc.data.img ~= 0));
      cscc.data.data_array(:,:,i) = cscc.data.img/cscc.data.imgscale;
    else
      cscc.data.data_array(:,:,i) = cscc.data.img;
    end
  end
  %%
  cscc.data.imgscale = median(cscc.data.data_array(cscc.data.data_array ~= 0)); % scale image according to mean
  cscc.data.data_array = cscc.data.data_array / cscc.data.imgscale * 0.3;

  % calculate individual difference to mean image
  for i=1:size(cscc.data.data_array_diff,3)
    cscc.data.data_array_diff(:,:,i) = cscc.data.data_array_diff(:,:,i) - mean(cscc.data.data_array_diff,3);
  end

  % enhance contrast and scale image to 0..64
  %mn = min(cscc.data.data_array(:));
  %mx = max(cscc.data.data_array(:));
  cscc.data.data_array = min(64,max(1,63 * cscc.data.data_array + 1));

  if cscc.H.sorted
    if isfield(cscc.pos,'x')
      x = cscc.data.ind_sorted(cscc.pos.x);
      if ~cscc.H.isscatter
        y = cscc.data.ind_sorted(cscc.pos.y);
      end
    end
  else
    if isfield(cscc.pos,'x')
      x = cscc.pos.x;
      if ~cscc.H.isscatter
        y = cscc.pos.y;
      end
    end
  end

  % check whether mouse position is defined
  if isfield(cscc.pos,'x')
    if cscc.H.isscatter
      cscc.data.img       = cscc.data.data_array(:,:,x)';
      cscc.data.img_alpha = cscc.data.data_array_diff(:,:,x)';
    else
      cscc.data.img       = [cscc.data.data_array(:,:,y) cscc.data.data_array(:,:,x)]';
      cscc.data.img_alpha = [cscc.data.data_array_diff(:,:,y) cscc.data.data_array_diff(:,:,x)]';
    end

    % correct orientation
    cscc.data.img = rot90(cscc.data.img,2);
    cscc.data.img_alpha = rot90(cscc.data.img_alpha,2);

    % use gray scale colormap for values > 64
    update_alpha
  
    set(cscc.H.mm_txt,'String',sprintf('%0.1f mm',get(cscc.H.mm,'Value')));
  end

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj) 
%-----------------------------------------------------------------------
  global cscc

  alphaval = get(cscc.H.alphabox,'Value');
  cscc.H.surfori  = 0; 
  
  if gca ~= cscc.H.slice 
    showtrash = get(cscc.H.showtrash,'Value');


    cscc.H.datatip = event_obj; 

    cscc.pos.pos_mouse = get(event_obj, 'Position');
    cscc.pos.tar_mouse = get(event_obj, 'Target');

    % Limit the number of datatips:
    % Although it is cscc.posible or maybe use to mark cscc.H.multiple objects for
    % comparision of volumes, surfaces and PDFs, the other check cases 
    % cscc.data.XML and LOG-files are worse to handle. 
    % Furthermore, it is elaboritiv to suport all cscc.select.trashlist operations.
    dcmlim = 0; % + 6*cscc.H.isscatter; % up to 6 plot objects?
    dcm = findall(gcf,'Type','hggroup','selected','off'); 
    if numel(dcm)>dcmlim, delete(dcm(dcmlim+1:end)); end
    set(findall(gcf,'Type','hggroup'),'Visible','on')



    % check for valid mouse position
    if cscc.H.isscatter
      cscc.pos.x = find(cscc.data.X(:,1) == cscc.pos.pos_mouse(1));
      if isempty(cscc.pos.x)
        cscc.pos.pos_mouse = get(event_obj, 'Position');
        cscc.pos.x = find(cscc.data.X(:,2) == cscc.pos.pos_mouse(2));
      end

      %{
      if isfield(cscc.job,'c')
        for ci=1:size(size(cscc.job.c,2))
          if round(cscc.job.c(x,ci)) == cscc.job.c(cscc.pos.x,ci)
            nuisx = sprintf('; N%d=%d',ci,cscc.job.c(cscc.pos.x,ci)); 
          else
            nuisx = sprintf('; N%d=%d',ci,cscc.job.c(cscc.pos.x,ci)); 
          end
        end
        if size(cscc.job.c,2)>3
          nuisx2 = ['\n              ' nuisx(3:end)]; nuisx = []; 
        end
      end
      %}
      
      % text info for data cursor window
      txt = {sprintf('S%d:%s',cscc.datagroups.sample(cscc.pos.x),cscc.files.fname.m{cscc.pos.x});
             sprintf('MNC=%5.3f',cscc.data.mean_cov(cscc.pos.x))};
      txt{2} = sprintf('%s; IQR=%5.2f',txt{2},cscc.data.QM(cscc.pos.x,4));
      if numel(cscc.datagroups.protocols)>1
        txt{2} = sprintf('%s (P=%0.2f)',txt{2},...
          cscc.datagroups.protocols(cscc.datagroups.protocol(cscc.pos.x)));
      end
      txt{2} = sprintf('%s; MD=%0.2f' ,txt{2},cscc.data.MD(cscc.pos.x));
     
      set(cscc.H.slice,'Position',cscc.pos.slice .* [1 1 1 0.5],'Visible','on');

      x = cscc.pos.x;
    else % covariance matrix
     
      if cscc.pos.pos_mouse(1) > cscc.pos.pos_mouse(2) || ...
          cscc.pos.pos_mouse(1)>length(cscc.datagroups.sample) || cscc.pos.pos_mouse(2)>length(cscc.datagroups.sample)
        %txt = {''}; 
        set([cscc.H.slice,cscc.H.alphabox],'Visible','off');
        if ~cscc.H.mesh_detected
          set([cscc.H.mm,cscc.H.mm_txt],'Visible','off'); 
        end
        set(get(cscc.H.slice,'children'),'Visible','off');
        unit = struct2cell(cscc.H.checkui); set([unit{cellfun(@ishandle,unit)}],'Enable','off'); 
        set([cscc.H.trashui.trash,cscc.H.trashui.detrash,...
             cscc.H.trashui.trashcol,cscc.H.trashui.detrashcol,...
             cscc.H.trashui.trashrow,cscc.H.trashui.detrashrow],'Enable','off');
        set(findall(gcf,'Type','hggroup'),'Visible','off')
        return
      elseif cscc.pos.pos_mouse(1)==1 && cscc.pos.pos_mouse(2)==1 && ...
        get(get(cscc.H.corr,'children'),'CData')==65
        txt = 'Error: No datapoints available.\nModify data selection!';
        return
      end

      % save position of mouse
      cscc.pos.x = find(cumsum(max(showtrash,cscc.select.trash1d) & ...
        cscc.select.samp1d & cscc.select.prot1d)==cscc.pos.pos_mouse(1),1,'first');
      cscc.pos.y = find(cumsum(max(showtrash,cscc.select.trash1d) & ...
        cscc.select.samp1d & cscc.select.prot1d)==cscc.pos.pos_mouse(2),1,'first');
    
      
      if cscc.H.sorted
        if isfield(cscc.pos,'x')
          x = cscc.data.ind_sorted(cscc.pos.x); cscc.pos.x = x; 
          y = cscc.data.ind_sorted(cscc.pos.y); cscc.pos.y = y; 
        end
      else
        if isfield(cscc.pos,'x')
          x = cscc.pos.x;
          y = cscc.pos.y;
        end
      end

      % text info for data cursor window
      if isfield(cscc.data,'QM') & size(cscc.data.QM,2)>=4
        QMtxtx = sprintf('; IQR=%5.2f; MD=%5.2f',cscc.data.QM(x,4),cscc.data.MD(x));
        QMtxty = sprintf('; IQR=%5.2f; MD=%5.2f',cscc.data.QM(y,4),cscc.data.MD(y));
      end
      nuisx = ''; nuisy = ''; 

      if isfield(cscc.job,'c') & ~isempty(cscc.job.c)
        for ci=1:size(cscc.job.c,2)
          if round(cscc.job.c{ci}(x)) == cscc.job.c{ci}(x)
            nuisx = sprintf('; N_%d=%d',ci,cscc.job.c{ci}(x)); 
            nuisy = sprintf('; N%d=%d',ci,cscc.job.c{ci}(y)); 
          else
            nuisx = sprintf('; N%d=%0.2f',ci,cscc.job.c{ci}(x)); 
            nuisy = sprintf('; N%d=%0.2f',ci,cscc.job.c{ci}(y)); 
          end
        end
        if size(cscc.job.c,2)>2 % new separate line
          nuisx = sprintf('Column (Top):  %s',nuisx(3:end)); 
          nuisy = sprintf('Row (Bottom):  %s',nuisy(3:end)); 
        end
      end

      if isfield(cscc.data,'QM') & size(cscc.data.QM,2)>=4
        txt = {
          sprintf('Correlation:      %3.3f',cscc.data.YpY(x,y)),...
          sprintf('Column (Top):  S%d:%s',cscc.datagroups.sample(x),cscc.files.fname.m{x}), ...
          sprintf('Row (Bottom):  S%d:%s',cscc.datagroups.sample(y),cscc.files.fname.m{y}), ...
          sprintf('Column (Top):  MNC=%5.3f%s%s',cscc.data.mean_cov(x),QMtxtx,nuisx),...
          sprintf('Row (Bottom):  MNC=%5.3f%s%s',cscc.data.mean_cov(y),QMtxty,nuisy)...
        };
      else
        txt = {
          sprintf('Correlation:      %3.3f',cscc.data.YpY(x,y)),...
          sprintf('Column (Top):  S%d:%s',cscc.datagroups.sample(x),cscc.files.fname.m{x}), ...
          sprintf('Row (Bottom):  S%d:%s',cscc.datagroups.sample(y),cscc.files.fname.m{y}), ...
          sprintf('Column (Top):  MNC=%5.3f%s',cscc.data.mean_cov(x),nuisx),...
          sprintf('Row (Bottom):  MNC=%5.3f%s',cscc.data.mean_cov(y),nuisy)...
        };
      end
      set(cscc.H.slice,'Position',cscc.pos.slice,'Visible','on');
    
    end

    % == check unit ==
    onoff = {'on','off'};
    if cscc.H.isscatter
      set(cscc.H.checkui.vol ,'Enable',onoff{ isempty(cscc.files.org{cscc.pos.x})+1  });
      set(cscc.H.checkui.surf,'Enable',onoff{ isempty(cscc.files.surf{cscc.pos.x})+1 });
      set(cscc.H.checkui.xml ,'Enable',onoff{ isempty(cscc.files.xml{cscc.pos.x})+1  });
      set(cscc.H.checkui.log ,'Enable',onoff{ isempty(cscc.files.log{cscc.pos.x})+1  });
      set(cscc.H.checkui.pdf ,'Enable',onoff{ isempty(cscc.files.pdf{cscc.pos.x})+1  });
    else
      set(cscc.H.checkui.vol ,'Enable',onoff{ (isempty(cscc.files.org{cscc.pos.x})  | isempty(cscc.files.org{cscc.pos.y}))  + 1 });
      set(cscc.H.checkui.surf,'Enable',onoff{ (isempty(cscc.files.surf{cscc.pos.x}) | isempty(cscc.files.surf{cscc.pos.y})) + 1 });
      set(cscc.H.checkui.xml ,'Enable',onoff{ (isempty(cscc.files.xml{cscc.pos.x})  | isempty(cscc.files.xml{cscc.pos.y}))  + 1 });
      set(cscc.H.checkui.log ,'Enable',onoff{ (isempty(cscc.files.log{cscc.pos.x})  | isempty(cscc.files.log{cscc.pos.y}))  + 1 });
      set(cscc.H.checkui.pdf ,'Enable',onoff{ (isempty(cscc.files.pdf{cscc.pos.x})  | isempty(cscc.files.pdf{cscc.pos.y}))  + 1 });
    end

    % == trash list unit ==
    if ~isempty(cscc.pos.x) 
      if all( cscc.select.trashlist~=cscc.pos.x ) 
        set([cscc.H.trashui.trash  ,cscc.H.trashui.trashcol  ],'Enable','on' );
        set([cscc.H.trashui.detrash,cscc.H.trashui.detrashcol],'Enable','off');
      else
        set([cscc.H.trashui.trash  ,cscc.H.trashui.trashcol  ],'Enable','off' );
        set([cscc.H.trashui.detrash,cscc.H.trashui.detrashcol],'Enable','on');
      end
    %else
    %  set([cscc.H.trashui.trash,cscc.H.trashui.detrash],'Enable','off');
    end
    if ~cscc.H.isscatter && isfield(cscc.pos,'y') && ~isempty(cscc.pos.y)
      if all( cscc.select.trashlist~=cscc.pos.y ) 
        set(cscc.H.trashui.trashrow  ,'Enable','on');
        set(cscc.H.trashui.detrashrow,'Enable','off');
      else
        set(cscc.H.trashui.trashrow  ,'Enable','off');
        set(cscc.H.trashui.detrashrow,'Enable','on');
      end
    end
    if ~isempty(cscc.select.trashlist)
      set([cscc.H.trashui.new,cscc.H.trashui.disptrash,cscc.H.trashui.ziptrash],'Enable','on');
    else
      set([cscc.H.trashui.new,cscc.H.trashui.disptrash,cscc.H.trashui.ziptrash],'Enable','off');
    end  

    set(cscc.H.alphabox,'Visible','on');



    if cscc.H.mesh_detected 
      % use indexed 2D-sheet to display surface data as image
      % check surface size to use indexed 2D map
      if (length(cscc.data.data_array(:,x)) == 163842) || (length(cscc.data.data_array(:,x)) == 32492)
        % combined surface 

        if (length(cscc.data.data_array(:,x)) == 163842) 
          ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','fsavg.index2D_256x128.txt'));
        else
          ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','fsavg.index2D_256x128.txt'));
        end

        if cscc.H.isscatter
          cscc.data.img = reshape(cscc.data.data_array(ind,x),[256,128]);
        else
          cscc.data.img = [reshape(cscc.data.data_array(ind,x),[256,128]) reshape(cscc.data.data_array(ind,y),[256,128])];
        end
        cscc.data.img = circshift(cscc.data.img,128);

        % alpha overlay
        if cscc.H.isscatter
          cscc.data.img_alpha = reshape(cscc.data.data_array_diff(ind,x),[256,128]);
        else
          cscc.data.img_alpha = [reshape(cscc.data.data_array_diff(ind,x),[256,128]) ...
            reshape(cscc.data.data_array_diff(ind,y),[256,128])];
        end
        cscc.data.img_alpha = circshift(cscc.data.img_alpha,128);


      elseif (length(cscc.data.data_array(:,x)) == 327684) || (length(cscc.data.data_array(:,x)) == 64984)
        is32k = length(cscc.data.data_array(:,x)) == 64984;
        if is32k
          atlasdir = 'atlases_surfaces_32k';
          tmpdir   = 'templates_surfaces_32k';
        else
          atlasdir = 'atlase_surfaces';
          tmpdir   = 'templates_surfaces'; 
        end
        lrb = [size(cscc.data.data_array,1)/2,size(cscc.data.data_array,1)/2+1];
        lab = fullfile(spm('dir'),'toolbox','cat12',atlasdir,'lh.aparc_DK40.freesurfer.annot'); 
        ind = spm_load(fullfile(spm('dir'),'toolbox','cat12',tmpdir,'fsavg.index2D_256x128.txt'));

        %% average both atlas maps to keep it simpler
        [vr{1},lb{1},tb{1}] = cat_io_FreeSurfer('read_annotation',lab);
        [vr{2},lb{2},tb{2}] = cat_io_FreeSurfer('read_annotation',strrep(lab,'lh.','rh.'));
        cscc.data.atlastable = {tb{1}.table(:,5), tb{1}.struct_names}; 
        cscc.data.atlas_lh = circshift(reshape(lb{1}(ind),[256,128]),64)';
        cscc.data.atlas_rh = circshift(reshape(lb{2}(ind),[256,128]),64)';
        cscc.data.atlas1   = cat_vol_median3(single(cat(3,cscc.data.atlas_lh,cscc.data.atlas_rh))); 
        [gx,gy]  = gradient(cscc.data.atlas1(:,:,1)); 
        cscc.data.atlasmsk = single(((abs(gy) + abs(gx))./cscc.data.atlas1(:,:,1))<0.05); 
        cscc.data.atlasmsk(cscc.data.atlasmsk==0) = nan; 
        cscc.data.atlasmsk = flipud(cscc.data.atlasmsk); 
        %cscc.data.atlas1   = flipud(circshift(reshape(cscc.data.atlas1(ind),[256,128]),64)');
        cscc.data.atlas1   = flipud(cscc.data.atlas1(:,:,1)); 
        
        %% load data array entry
        cscc.data.data_array_x_lh = cscc.data.data_array(1:lrb(2),x);
        cscc.data.data_array_x_rh = cscc.data.data_array(lrb(2):end,x);
        if ~cscc.H.isscatter
          cscc.data.data_array_y_lh = cscc.data.data_array(1:lrb(1),y);
          cscc.data.data_array_y_rh = cscc.data.data_array(lrb(2):end,y);
        end

        % rotate the image (') and shift it (64) to have a cscc.posterior cutting edge 
        cscc.data.img_x_lh = flipud(circshift(reshape(cscc.data.data_array_x_lh(ind),[256,128]),64)');
        cscc.data.img_x_rh = flipud(circshift(reshape(cscc.data.data_array_x_rh(ind),[256,128]),64)');
        if ~cscc.H.isscatter
          cscc.data.img_y_lh = flipud(circshift(reshape(cscc.data.data_array_y_lh(ind),[256,128]),64)');
          cscc.data.img_y_rh = flipud(circshift(reshape(cscc.data.data_array_y_rh(ind),[256,128]),64)');
        end
        %%
        if cscc.H.surfori
          if cscc.H.isscatter
            cscc.data.img   = [cscc.data.img_x_lh .* cscc.data.atlasmsk', cscc.data.img_x_rh .* cscc.data.atlasmsk'];
            cscc.data.atlas = [cscc.data.atlas1'; cscc.data.atlas1'];
          else
            cscc.data.img   = [cscc.data.img_x_lh' .* cscc.data.atlasmsk', cscc.data.img_y_lh' .* cscc.data.atlasmsk'; ...
                               cscc.data.img_x_rh' .* cscc.data.atlasmsk', cscc.data.img_y_rh' .* cscc.data.atlasmsk'];
            cscc.data.atlas = [cscc.data.atlas1', cscc.data.atlas1'; ...
                               cscc.data.atlas1', cscc.data.atlas1'];
          end
        else
          if cscc.H.isscatter
            cscc.data.img   = [cscc.data.img_x_lh' .* cscc.data.atlasmsk', fliplr(cscc.data.img_x_rh' .* cscc.data.atlasmsk')];
            cscc.data.imgo  = [cscc.data.img_x_lh', fliplr(cscc.data.img_x_rh')];
            cscc.data.atlas = [cscc.data.atlas1', fliplr(cscc.data.atlas1')];
          else
            cscc.data.img   = [cscc.data.img_x_lh' .* cscc.data.atlasmsk', fliplr(cscc.data.img_x_rh' .* cscc.data.atlasmsk'); ...
                               cscc.data.img_y_lh' .* cscc.data.atlasmsk', fliplr(cscc.data.img_y_rh' .* cscc.data.atlasmsk')];
            cscc.data.imgo  = [cscc.data.img_x_lh', fliplr(cscc.data.img_x_rh'); ...
                               cscc.data.img_y_lh', fliplr(cscc.data.img_y_rh')];
            cscc.data.atlas = [cscc.data.atlas1' fliplr(cscc.data.atlas1'); ...
                               cscc.data.atlas1' fliplr(cscc.data.atlas1')];
          end
        end
        
        %% alpha overlay
        cscc.data.data_array_x_lh = cscc.data.data_array_diff(1:lrb(2),x);
        cscc.data.data_array_x_rh = cscc.data.data_array_diff(lrb(2):end,x);
        if ~cscc.H.isscatter
          cscc.data.data_array_y_lh = cscc.data.data_array_diff(1:lrb(1),y);
          cscc.data.data_array_y_rh = cscc.data.data_array_diff(lrb(2):end,y);
        end
        cscc.data.img_x_lh = flipud(circshift(reshape(cscc.data.data_array_x_lh(ind),[256,128]),64)');
        cscc.data.img_x_rh = flipud(circshift(reshape(cscc.data.data_array_x_rh(ind),[256,128]),64)');
        if ~cscc.H.isscatter
          cscc.data.img_y_lh = flipud(circshift(reshape(cscc.data.data_array_y_lh(ind),[256,128]),64)');
          cscc.data.img_y_rh = flipud(circshift(reshape(cscc.data.data_array_y_rh(ind),[256,128]),64)');
        end
        
        if cscc.H.surfori  
          if cscc.H.isscatter
            cscc.data.img_alpha = [cscc.data.img_x_lh .* cscc.data.atlasmsk; cscc.data.img_x_rh .* fliplr(cscc.data.atlasmsk)];
          else
            cscc.data.img_alpha = [cscc.data.img_x_lh .* cscc.data.atlasmsk; cscc.data.img_y_lh .* fliplr(cscc.data.atlasmsk); ...
                         cscc.data.img_x_rh .* cscc.data.atlasmsk; cscc.data.img_y_rh .* fliplr(cscc.data.atlasmsk)];
          end
        else
          if cscc.H.isscatter
            cscc.data.img_alpha   = [cscc.data.img_x_lh' .* cscc.data.atlasmsk', fliplr(cscc.data.img_x_rh' .* cscc.data.atlasmsk')];
          else
            cscc.data.img_alpha   = [cscc.data.img_x_lh' .* cscc.data.atlasmsk', fliplr(cscc.data.img_x_rh' .* cscc.data.atlasmsk'); ...
                           cscc.data.img_y_lh' .* cscc.data.atlasmsk', fliplr(cscc.data.img_y_rh' .* cscc.data.atlasmsk')];
          end
        end

      else
        if cscc.H.isscatter
          cscc.data.img = cscc.data.data_array(:,x)';
          % alpha overlay
          cscc.data.img_alpha = cscc.data.data_array_diff(:,x)';
        else
          cscc.data.img = [cscc.data.data_array(:,y) cscc.data.data_array(:,x)]';
          % alpha overlay
          cscc.data.img_alpha = [cscc.data.data_array_diff(:,y) cscc.data.data_array_diff(:,x)]';
        end
      end

      % scale cscc.data.img to 0..64
      if 0 %strfind(cscc.files.dataprefix,'thickness')
        mn = 0.0;
        mx = 5.0; 
      else
        sd = std(cscc.data.data_array(:));
        mn = -sd*2;
        mx =  sd*2;
      end
      cscc.data.img  = 64*((cscc.data.img - mn)/(mx-mn));
    else
      % add slider for colume data
      set(cscc.H.mm,'Visible','on');
      set(cscc.H.mm_txt,'Visible','on');
      if cscc.H.isscatter
        cscc.data.img = cscc.data.data_array(:,:,x)';
        % alpha overlay
        cscc.data.img_alpha = cscc.data.data_array_diff(:,:,x)';
      else
        cscc.data.img = [cscc.data.data_array(:,:,y) cscc.data.data_array(:,:,x)]';
        % alpha overlay
        cscc.data.img_alpha = [cscc.data.data_array_diff(:,:,y) cscc.data.data_array_diff(:,:,x)]';
      end
    end

  %  set(cscc.H.cbar,'TickLength',[0 0],'XTickLabel',linspace(0.7,1.0,7)); 
  %round(100*linspace(min(cscc.data.YpY(cscc.select.trash2d(:))),max(cscc.data.YpY(cscc.select.trash2d(:))),5))/100);


    % correct orientation
    cscc.data.img       = rot90(cscc.data.img,2);
    cscc.data.img_alpha = rot90(cscc.data.img_alpha,2);
    if cscc.H.mesh_detected
      cscc.data.imgo   = rot90(cscc.data.imgo,2);
      cscc.data.atlas   = rot90(cscc.data.atlas,2);
     end    

    % display image with 2nd colorbar (gray)
    image(cscc.H.slice,65 + cscc.data.img); 
    set(cscc.H.slice,'visible','off');

    % prepare alpha overlays for red and green colors
    if alphaval > 0
      %% get 2%/98% ranges of difference image
      range = cat_vol_iscaling(cscc.data.img_alpha(:),[0.02 0.98]);

      hold(cscc.H.slice,'on');
      alpha_g = cat(3, zeros(size(cscc.data.img_alpha)), ...
        alphaval*ones(size(cscc.data.img_alpha)), zeros(size(cscc.data.img_alpha)));
      alpha_r = cat(3, alphaval*ones(size(cscc.data.img_alpha)), ...
        zeros(size(cscc.data.img_alpha)), zeros(size(cscc.data.img_alpha)));
      hg = image(cscc.H.slice,alpha_g); 
      set(hg, 'AlphaData', cscc.data.img_alpha.*(cscc.data.img_alpha>range(2)));
      if ~cscc.H.mesh_detected, set(hg,'HitTest','off','Interruptible','off','AlphaDataMapping','scaled'); end
      hr = image(cscc.H.slice,alpha_r); 
      set(hr, 'AlphaData',-cscc.data.img_alpha.*(cscc.data.img_alpha<range(1)));
      if ~cscc.H.mesh_detected, set(hg,'HitTest','off','Interruptible','off','AlphaDataMapping','scaled'); end
      hold(cscc.H.slice,'off');
    end

    if cscc.H.mesh_detected
      xlabel('2D surface maps');
    end
  else
    cscc.pos.pos_mouse = get(event_obj, 'Position');
    if cscc.H.surfori
      ax   = mod(cscc.pos.pos_mouse(1) - 1,256) + 1;
      ay   = mod(cscc.pos.pos_mouse(2) - 1,128) + 1;
    else
      if cscc.pos.pos_mouse(1)<128
        ay = cscc.pos.pos_mouse(1);
      else
        ay = 257 - cscc.pos.pos_mouse(1);
      end
      ax   = mod(cscc.pos.pos_mouse(2) - 1,256);
    end
    rid  = cscc.data.atlas(ax,ay);
    ridi = find(rid==cscc.data.atlastable{1},1); 
    
    if cscc.pos.pos_mouse(1)<128, side = 'lh'; else side = 'rh'; end
    if ~isempty(ridi)
      roi = [side '.' cscc.data.atlastable{2}{ridi}];
    else
      roi = 'unknown';
    end
    
    if alphaval > 0  
      dt = cscc.data.imgo(ax,ay);
    else
      dt = cscc.data.img_alpha(ax,ay);
    end
    
    txt = {sprintf('%s (%d,%d)',roi,ax,ay) sprintf('Value: %0.2f',dt) }; 
  end
return

%-----------------------------------------------------------------------
function varargout = cat_tst_qa_cleaner_intern(data,opt)
%% Tcscc.HIS FUNCTION IS TO FAT - SEPARATE AND CLEAN IT! 
%  Do not forget to remove old external version from Scscc.data.VN if this is done.
%  _____________________________________________________________________
%  Estimate quality grades of given rating of one (or more) cscc.datagroups.protocols
%  with 2 to 6 grads to separate passed, (unassignable) and failed 
%  images, by finding the first peak in the image quality histogram  
%  and using its width (standard deviation) in a limited range. 
%  If cscc.H.multiple cscc.datagroups.protocols are used, than use the site variable opt.site 
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
%   * cscc.H.multiside output required a 'stacked' output
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
% $Id: cat_stat_check_cov2.m 1612 2020-05-03 12:18:52Z gaser $ 

  global cscc

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
    data = cellstr(spm_select(inf,'XML','select qa cscc.data.XML-files',{},pwd,'^cat_.*')); 
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
    fprintf('Load cscc.data.XML data');
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
    %bar(r,h/sh,'facecolor',cscc.display.figcolor,'edgecolor','none');
    %fill(r,h/sh,cscc.display.figcolor,'edgecolor','none');
    hold on
    
    yl = [0 max(h)+1]; ylim(yl);
    % main grid
    for i=1.5:6,       plot([i i],ylim,'color',cscc.display.figcolor); end
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
    cscc.display.FS = get(gca,'Fontsize')*1.3;
    set(gca,'XTick',0.5:1:6.5,'XTickLabel',{'100','90','80','70','60','50','40'},'TickLength',[0.02 0.02]);
    % further color axis objects...
    axA = copyobj(gca,gcf); axB = copyobj(axA,gcf); axC = copyobj(gca,gcf); 
    axD = copyobj(gca,gcf); axE = copyobj(gca,gcf); axF = copyobj(gca,gcf);
    % set colors...
    set(axA,'YTick',[],'XTickLabel',{},'XTick',1,...
      'XColor',color(QMC,1),'Color','none',...
      'XTicklabel','A','TickLength',[0 0],...
      'Fontsize',cscc.display.FS,'Fontweight','bold');
    set(axB,'YTick',[],'XTickLabel',{},'XTick',2,...
      'XColor',color(QMC,2),'Color','none',...
      'XTicklabel','B','TickLength',[0 0],...
      'Fontsize',cscc.display.FS,'Fontweight','bold');
    set(axC,'YTick',[],'XTickLabel',{},'XTick',3,...
      'XColor',color(QMC,3),'Color','none',...
      'XTicklabel','C','TickLength',[0 0],...
      'Fontsize',cscc.display.FS,'Fontweight','bold');
    set(axD,'YTick',[],'XTickLabel',{},'XTick',4,...
      'XColor',color(QMC,4),'Color','none',...
      'XTicklabel','D','TickLength',[0 0],...
      'Fontsize',cscc.display.FS,'Fontweight','bold');
    set(axE,'YTick',[],'XTickLabel',{},'XTick',5,...
      'XColor',color(QMC,5),'Color','none',...
      'XTicklabel','E','TickLength',[0 0],...
      'Fontsize',cscc.display.FS,'Fontweight','bold');
    set(axF,'YTick',[],'XTickLabel',{},'XTick',6,...
      'XColor',color(QMC,6),'Color','none',...
      'XTicklabel','F','TickLength',[0 0],...
      'Fontsize',cscc.display.FS,'Fontweight','bold');
    hold off; 
    
    if isfield(opt,'site') && numel(sites>1);
      title(sprintf('Histogram (cf=%0.2f) - global treshold for cscc.H.multisite output (n=%d)',...
        opt.cf,numel(sites)),'Fontsize',cscc.display.FS);
    else
      title(sprintf('Histogram (cf=%0.2f)',opt.cf),...
        'Fontsize',cscc.display.FS);
    end
    xlabel('IQR (rps)','Fontsize',cscc.display.FS); 
    ylabel('number of scans','Fontsize',cscc.display.FS); 
  end
  %%
  MarkColor = cat_io_colormaps('marks+',40); 
  if isfield(opt,'site') && numel(sites)>1, globcorr = ' (global corrected)'; else globcorr = ''; end
  if exist('P','var')
    files = P(data<=markths2(:,3)); 
    fprintf('PASSED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    
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
  
 
return
