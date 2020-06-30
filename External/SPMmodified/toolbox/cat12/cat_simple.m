function output = cat_simple(job)
% Configuration file for simplyfied preprocessing jobs
% _________________________________________________________________________
% Robert Dahnke
% $Id: cat_simple.m 1582 2020-03-11 16:55:11Z gaser $


  % defaults settings
  def.catversion    = 'estwrite';
  def.tpm           = fullfile(spm('dir'),'TPM','TPM.nii');
  def.nproc         = cat_get_defaults('extopts.nproc');
  def.debug         = 0;
  def.ignoreErrors  = 1;
  job = cat_io_checkinopt(job,def); 
  
  if isdeployed, job.nproc = 0; end
  
  expert    = cat_get_defaults('extopts.expertgui'); 
  proc_surf = isfield(job.surface,'yes');

  % ROI atlas setup
  % -----------------------------------------------------------------------
  % volume default ROIs
  if isfield(job,'ROImenu') % expert/developer GUI that allows control each atlas map 
    if isfield(job.ROImenu,'atlases')
      % image output
      def.atlases = job.ROImenu.atlases;
    end
    job = cat_io_checkinopt(job,def);
  end
  
  % surface default ROIs
  if ~proc_surf && isfield(job,'sROImenu') && ~isfield(job.sROImenu,'noROI')
    cat_io_cprintf('err', ...
      ['\nYou cannot select surface ROI analysis without surface processing.\n' ...
       'Select surface processing or turn of surface ROI analysis.\n\n']);
    output = {};
    return
  end
  
  if isfield(job,'sROImenu') % expert/developer GUI that allows control for each atlas map 
    if isfield(job.sROImenu,'satlases')
      % image output
      def.satlases = job.sROImenu.satlases;
    end
    job = cat_io_checkinopt(job,def);
  end
  
  
  %% longitudinal data handling
  %  ----------------------------------------------------------------------
  if isfield(job,'subj') && isfield(job,'data')
    long      = 1; 
    pfield    = 'data'; 
  elseif isfield(job,'datalong') && isfield(job.datalong,'timepoints')
    % translate estimate into long structure
    job.catversion = strrep(job.catversion,'estwrite','long');
    long      = 1; 
    pfield    = 'subj'; 
    job.data  = {};

    for ti = 1:numel(job.datalong.timepoints)
      for si = 1:numel(job.datalong.timepoints{ti})
        job.subj(si).mov{ti,1} = job.datalong.timepoints{ti}{si}; 
      end
      job.data = [job.data job.datalong.timepoints{ti}];
    end

% ######   
% check input for same number of subjects for timepoint definition
% ######
  elseif isfield(job,'datalong') && isfield(job.datalong,'subjects')
    % translate estimate into long structure
    job.catversion = strrep(job.catversion,'estwrite','long');
    long      = 1; 
    pfield    = 'subj'; 
    job.data  = {};

    for si = 1:numel(job.datalong.subjects)
      for ti = 1:numel(job.datalong.subjects{si})
        job.subj(si).mov{ti,1} = job.datalong.subjects{si}{ti}; 
      end
      job.data = [job.data; job.datalong.subjects{si}];
    end
  else
    long      = 0; 
    pfield    = 'data'; 
  end
  
  
  
  %% split job and data into separate processes to save computation time
  %  ----------------------------------------------------------------------
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && (numel(job.data)>1)
    if nargout==1
      output = cat_parallelize(job,mfilename,pfield);
    else
      cat_parallelize(job,mfilename,pfield);
    end
    return
  end
 
  
  
  % preprocessing / long
  % -----------------------------------------------------------------------
  % specification of catversion and their batch-tag and name setting
  switch job.catversion
    case 'long',              estwrite = 'Segment longitudinal data:';
    case 'long1173',          estwrite = 'CAT12.1: Segment longitudinal data R1173 (2017/09)';
    case 'long1173plus',      estwrite = 'CAT12.3: Segment longitudinal data R1392 (2018/12)';
    case 'long1445',          estwrite = 'CAT12.6: Segment longitudinal data R1445 (2019/03)';
    case 'long1585',          estwrite = 'CAT12.7: Segment longitudinal data R1585 (2020/03)';
    case 'estwrite',          estwrite = 'CAT12: Segmentation:';
    case 'estwrite1173',      estwrite = 'CAT12.1: Segmentation R1173 (2017/09):'; %.1
    case 'estwrite1173plus',  estwrite = 'CAT12.3: Segmentation R1173 (2018/12):'; %.3
    case 'estwrite1445',      estwrite = 'CAT12.6: Segmentation R1445 (2019/03)';
    case 'estwrite1585',      estwrite = 'CAT12.7: Segmentation R1585 (2020/03)';
  end
  
  
  
  mbi = 1;

  % set input files (see handling for longitudinal data above)
  if long
    matlabbatch{mbi}.spm.tools.cat.(job.catversion).subj = job.subj; 
    temp   = [job.subj(:).mov]; 
    [Cs,C] = spm_str_manip(temp(:),'C'); 
    nsub   = numel(job.data); 
  else
    matlabbatch{mbi}.spm.tools.cat.(job.catversion).data = job.data;
    [Cs,C] = spm_str_manip(job.data,'C'); 
    nsub   = numel(job.data); 
  end
  if isempty(C), clear C; C.s = Cs; end
  mdir = spm_fileparts(C.s); % define main data directory for later ...
  
  % set templates:
  % - keep in mind all TPMs should be in the same MNI space to support the 
  %   standard atlases and that also these special templates has to be
  %   recognised by CAT (currently is deactives ROI processing! RD201904)
  switch job.tpm
    case 'adults'
      matlabbatch{mbi}.spm.tools.cat.(job.catversion).opts.tpm =  ...
        {fullfile(spm('dir'),'tpm','TPM.nii')};
    case 'childen'
      matlabbatch{mbi}.spm.tools.cat.(job.catversion).opts.tpm = ...
        {fullfile(spm('dir'),'toolbox','cat12','templates_volumes','TPM_Age11.5.nii')};
  end
  
  % setting of standard fields
  matlabbatch{mbi}.spm.tools.cat.(job.catversion).nproc          = 0;
  matlabbatch{mbi}.spm.tools.cat.(job.catversion).output.GM.mod  = 1; 
  matlabbatch{mbi}.spm.tools.cat.(job.catversion).output.WM.mod  = 1; 
  matlabbatch{mbi}.spm.tools.cat.(job.catversion).output.surface = proc_surf;
  
  % set registration
  matlabbatch{mbi}.spm.tools.cat.(job.catversion).extopts.registration = job.registration; 

  if expert
    matlabbatch{mbi}.spm.tools.cat.(job.catversion).extopts.admin.ignoreErrors = job.ignoreErrors; 
  end
  
  if isfield(job,'atlases')
    if long
      matlabbatch{mbi}.spm.tools.cat.(job.catversion).ROImenu.atlases = job.atlases; % only volume ROIs
    else
      matlabbatch{mbi}.spm.tools.cat.(job.catversion).output.ROImenu.atlases = job.atlases; % only volume ROIs
    end
  end
  
  % developer and new version only! 
  % for fast tests of the whole pipeline of the developer mode
  if job.debug 
    switch job.catversion
      case {'estwrite1445','long1445','estwrite1585','long1585'}
        matlabbatch{mbi}.spm.tools.cat.(job.catversion).extopts.registration.regstr       = eps;                 % fast shooting only in new versions! 
        matlabbatch{mbi}.spm.tools.cat.(job.catversion).extopts.admin.lazy                = 0;                   % use lazy .. did not work yet - missing DEPs
        matlabbatch{mbi}.spm.tools.cat.(job.catversion).output.surface                    = double(proc_surf)*7; % use 0.8 mm for pbt ans fast registration (=7)
    end
    matlabbatch{mbi}.spm.tools.cat.(job.catversion).extopts.segmentation.restypes.fixed   = [2 0];               % use only 2 mm for VBM preprocessing
  end
  
  % here we have to (re)move some fields!
  if long
    matlabbatch{mbi}.spm.tools.cat.tools = matlabbatch{1}.spm.tools.cat;
    matlabbatch{mbi}.spm.tools.cat = rmfield(matlabbatch{1}.spm.tools.cat,(job.catversion)); 
    surf = matlabbatch{mbi}.spm.tools.cat.tools.(job.catversion).output.surface;
    matlabbatch{mbi}.spm.tools.cat.tools.(job.catversion) = rmfield(matlabbatch{mbi}.spm.tools.cat.tools.(job.catversion),'output'); 
    matlabbatch{mbi}.spm.tools.cat.tools.(job.catversion).output.surface = surf; clear surf; 
    matlabbatch{mbi}.spm.tools.cat.tools.(job.catversion).modulate = 1;
  end
    
  % for dependency objects
  if isfield(job,'process_index')
    addpath(fullfile(spm('dir'),'matlabbatch'));
  end
 
  
  
  % further modalities
  % -----------------------------------------------------------------------
  % bad topic this is for the next decade ...
  if isfield(job,'mod')
    cat_io_cprintf('err','Additional modalities are not supported right now!\n'); 
  end
  
  
  
  % volumetric smoothing
  % -----------------------------------------------------------------------
  vsmooth = job.fwhm_vol;
  for si = 1:numel(vsmooth)
		if long
		  for di = 1:numel(job.datalong.subjects)
				mbi = mbi + 1;  
				matlabbatch{mbi}.spm.spatial.smooth.data = cfg_dep(...
				sprintf('%s Segmented longitudinal data (Subj %d)',estwrite, di), ...
				substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
				substruct('.','sess', '()',{1}, '.','files'));
				matlabbatch{mbi}.spm.spatial.smooth.fwhm    = repmat( vsmooth(si) ,1,3);
				matlabbatch{mbi}.spm.spatial.smooth.prefix  = sprintf('s%d',vsmooth(si));
			end
		else
			for ti = 1:2
				mbi = mbi + 1;
				matlabbatch{mbi}.spm.spatial.smooth.data = cfg_dep( ...
					sprintf('%s mwp%d Image',estwrite,ti), ...
					substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
					substruct('.','tiss', '()',{ti}, '.','mwp', '()',{':'}));
				matlabbatch{mbi}.spm.spatial.smooth.fwhm    = repmat( vsmooth(si) ,1,3);
				matlabbatch{mbi}.spm.spatial.smooth.prefix  = sprintf('s%d',vsmooth(si));
			end
		end
  end

  
  
  
  
  
  
  % surfaces 
  % -----------------------------------------------------------------------
  if proc_surf
    % estimate surface measures 
    % ---------------------------------------------------------------------
    mbi = mbi + 1; 
    surf_mbi = mbi;
    if long
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.data_surf   = ...
         cfg_dep(sprintf('%s Left Central Surfaces',estwrite), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','surf', '()',{':'}));
    else
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.data_surf   = ...
        cfg_dep(sprintf('%s Left Central Surface',estwrite), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('()',{1}, '.','lhcentral', '()',{':'})); 
    end
    if expert>1 %developer
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.GIL.GIL     = 0;
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.surfaces.IS = 0;
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.surfaces.OS = 0;
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.area        = 0;
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.gmv         = 0;
    elseif expert % expert
      % nothing right now
    end
    matlabbatch{mbi}.spm.tools.cat.stools.surfextract.GI            = 1;
    matlabbatch{mbi}.spm.tools.cat.stools.surfextract.FD            = 0;
    matlabbatch{mbi}.spm.tools.cat.stools.surfextract.SD            = 0;
    matlabbatch{mbi}.spm.tools.cat.stools.surfextract.nproc         = 0;
    if job.debug
      matlabbatch{mbi}.spm.tools.cat.stools.surfextract.lazy        = 1; 
    end
  


    % surface smoothing
    % -----------------------------------------------------------------------
		if proc_surf
			ssmooth1  = job.surface.yes.fwhm_surf1;
			ssmooth2  = job.surface.yes.fwhm_surf2;
		end
    measures = {};
    
    % prepare the datafields depending on the internal selection above
    if isfield(matlabbatch{mbi}.spm.tools.cat.stools.surfextract,'GI') && matlabbatch{mbi}.spm.tools.cat.stools.surfextract.GI
      measures = [measures; {'gyrification', 'lPGI'}]; 
    end
    
    if isfield(matlabbatch{mbi}.spm.tools.cat.stools.surfextract,'FD') && matlabbatch{mbi}.spm.tools.cat.stools.surfextract.FD
      measures = [measures; {'fractal dimension', 'lPFD'}]; 
    end
    
    if isfield(matlabbatch{mbi}.spm.tools.cat.stools.surfextract,'SD') && matlabbatch{mbi}.spm.tools.cat.stools.surfextract.SD
      measures = [measures; {'sulcal depth', 'lPSD'}]; 
    end
    
    if isfield(matlabbatch{mbi}.spm.tools.cat.stools.surfextract,'GIL') && ...
     ((isstruct(matlabbatch{mbi}.spm.tools.cat.stools.surfextract.GIL) && ...
       matlabbatch{mbi}.spm.tools.cat.stools.surfextract.GIL.GIL) || ...
      (isnumeric(matlabbatch{mbi}.spm.tools.cat.stools.surfextract.GIL) && ...
       matlabbatch{mbi}.spm.tools.cat.stools.surfextract.GIL))
      measures = [measures; {'generalizedGI', 'lPgGI'}]; 
      if expert
        measures = [measures; {'inwardGI', 'lPiGI'}];
        measures = [measures; {'outwarGI', 'lPoGI'}]; 
      end
    end
    
    if isfield(matlabbatch{mbi}.spm.tools.cat.stools.surfextract,'area') && matlabbatch{mbi}.spm.tools.cat.stools.surfextract.area
      measures = [measures; {'surface area', 'lPara'}]; 
    end
    
    if isfield(matlabbatch{mbi}.spm.tools.cat.stools.surfextract,'gmv') && matlabbatch{mbi}.spm.tools.cat.stools.surfextract.gmv
      measures = [measures; {'surface GM volume', 'lPGIL'}]; 
    end
    
    % create surface smoothing batch for thickness
    for si = 1:numel(ssmooth1)
      mbi = mbi + 1;  
      % thickness
      if long
        matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.data_surf(1) = ...
         cfg_dep(sprintf('%s Left Thickness',estwrite), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','thick', '()',{':'}));
      else
        matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.data_surf(1) = ...
          cfg_dep(sprintf('%s Left Thickness',estwrite), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('()',{1}, '.','lhthickness', '()',{':'}));
      end
            
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.merge_hemi   = 1;
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.mesh32k      = 1;
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.fwhm_surf    = ssmooth1(si);
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.nproc        = 0;
    end

    for si = 1:numel(ssmooth2)
      mbi = mbi + 1;        
      % further parameters
      for mi = 1:size(measures,1)
        matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.data_surf(mi) = ...
          cfg_dep(sprintf('Extract additional surface parameters: Left %s',measures{mi,1}), ...
          substruct('.','val', '{}',{surf_mbi}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('()',{1}, '.',measures{mi,2}, '()',{':'}));
      end
      
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.merge_hemi   = 1;
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.mesh32k      = 1;
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.fwhm_surf    = ssmooth2(si);
      matlabbatch{mbi}.spm.tools.cat.stools.surfresamp.nproc        = 0;
    end
    
    
    
    % surface ROI data mapping
    if 1
      rdata  = cat_get_defaults('extopts.satlas');
      rfiles = cell(size(rdata,1),1);
      for ri = 1:size(rdata,1) 
        if exist(rdata{ri,2},'file')
          rfiles{ri} = rdata{ri,2};
        else
          rfiles{ri} = cat_vol_findfiles(fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces_32k'),['lh.' rdata{ri,2} '*annot']); 
        end
      end
      for mi = 1:size(measures,1)
        mbi = mbi + 1;  
        matlabbatch{mbi}.spm.tools.cat.stools.surf2roi.cdata{mi} = ...
          cfg_dep(sprintf('Extract additional surface parameters: Left %s',measures{mi,1}), ...
          substruct('.','val', '{}',{surf_mbi}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('()',{1}, '.',measures{mi,2}, '()',{':'}));
        if expert>1
          matlabbatch{mbi}.spm.tools.cat.stools.surf2roi.rdata = rfiles;
          matlabbatch{mbi}.spm.tools.cat.stools.surf2roi.nproc = 0;
        end
      end
    end
  end
  
  
  
  % estimate TIV & CGW
  % -----------------------------------------------------------------------
  % * I am not sure if this is useful here because it is befor any removal
  %   of files
  % * this estimation should be replaced in future by a more general script
  % -----------------------------------------------------------------------
  mydata = datestr(clock,'YYYYddmm-HHMM'); 
  mbi = mbi + 1;  
  matlabbatch{mbi}.spm.tools.cat.tools.calcvol.data_xml(1)    = cfg_dep(...
    sprintf('%s CAT Report',estwrite),...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','catreport', '()',{':'})); 
  matlabbatch{mbi}.spm.tools.cat.tools.calcvol.calcvol_TIV    = 0;
  matlabbatch{mbi}.spm.tools.cat.tools.calcvol.calcvol_name   = fullfile(mdir,sprintf('TIV_%s_N%d.txt',mydata,nsub)); 

  
  % estimate IQR 
  mbi = mbi + 1;  
  matlabbatch{mbi}.spm.tools.cat.tools.iqr.data_xml           = cfg_dep(...
    sprintf('%s CAT Report',estwrite),...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','catreport', '()',{':'})); 
  matlabbatch{mbi}.spm.tools.cat.tools.iqr.iqr_name           = fullfile(mdir,sprintf('IQR_%s_N%d.txt',mydata,nsub)); 
  
  
  % Delete temporary files? 
  % - No, the preprocessing results are required or can be used later 
  
  
  %% create output
  %  ----------------------------------------------------------------------
  if cat_get_defaults('extopts.subfolders')
    roifolder  = 'label';
    mrifolder  = 'mri';
    surffolder = 'surf';
    repfolder  = 'report';
  else
    roifolder  = '';
    mrifolder  = '';
    surffolder = '';
    repfolder  = 'report';
  end
  
  % unsmoothed segmentations
  for fi = 1:numel(job.data)
    output.mwp1{fi} = spm_file(job.data,'prefix',fullfile(mrifolder,'mwp1'));
    output.mwp2{fi} = spm_file(job.data,'prefix',fullfile(mrifolder,'mwp2'));
  end
  
  if proc_surf && exist('measures','var')
    for mi = 1:size(measures,1)
      for fi = 1:numel(job.data)
        output.measures{mi,1}{fi} = spm_file(job.data,'prefix',fullfile(surffolder,sprintf('mesh.%s.',measures{mi,1})),'ext','.gii');
      end
    end
  end
  
  % smoothed data
  for si = 1:numel(vsmooth)
    for fi = 1:numel(job.data)
      output.(sprintf('s%dmwp1',vsmooth(si))){fi} = spm_file(job.data,'prefix',fullfile(mrifolder,sprintf('s%dmwp1',vsmooth(si))));
      output.(sprintf('s%dmwp2',vsmooth(si))){fi} = spm_file(job.data,'prefix',fullfile(mrifolder,sprintf('s%dmwp2',vsmooth(si))));
    end
  end
  
  if proc_surf && exist('measures','var')
    for mi = 1:size(measures,1)
      if strcmp(measures{mi,1},'thickness')
        for si = 1:numel(ssmooth1)
          for fi = 1:numel(job.data)
            output.(sprintf('s%d%s',ssmooth1(si),measures{mi,1})){fi} = ...
              spm_file(job.data,'prefix',fullfile(surffolder,sprintf('s%dmm.mesh.%s.',ssmooth1(si),measures{mi,1})),'ext','.gii');
          end
        end
      else
        for si = 1:numel(ssmooth2)
          for fi = 1:numel(job.data)
            output.(sprintf('s%d%s',ssmooth2(si),measures{mi,1})){fi} = ...
              spm_file(job.data,'prefix',fullfile(surffolder,sprintf('s%dmm.mesh.%s.',ssmooth2(si),measures{mi,1})),'ext','.gii');
          end
        end
      end
    end
  end
  
  % xml data
  for fi = 1:numel(job.data)
    output.catroi{fi} = spm_file(job.data,'prefix',fullfile(roifolder,'catROI_'),'ext','.xml');
  end
  
  for fi = 1:numel(job.data)
    output.catxml{fi} = spm_file(job.data,'prefix',fullfile(repfolder,'cat_'),'ext','.xml');
  end
    
  spm_jobman('initcfg'); % reinitialize
  spm_jobman('run',matlabbatch); %,inputs{:});
end
