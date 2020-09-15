function cat_main_roi(job,trans,Ycls,Yp0) 
% ______________________________________________________________________
%  ROI Partitioning:
%  This part estimates individual measurements for different ROIs.
%  The ROIs are defined in the CAT normalized space and there are three 
%  ways to estimate them: (1) in (internal) subject space, (2) in 
%  normalized space (where also the VBM is done and that is defined 
%  by extopts.vox, and (3) in the atlas space.
%  Estimation in normalized space is more direct and avoids further
%  transformations or individual adaptions. The way over the subject space
%  has the advantage that individual anatomical refinements are possible, 
%  but this has to be done and evaluated for each atlas and it is not so 
%  simple at the end. Another thing (that came up later) was the evaluation
%  in atlas space, that most relevant because some atlas maps use now
%  higher resolutions to describe fine structures. 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_roi.m 1589 2020-03-19 14:41:22Z dahnke $
  
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  [pth,nam] = spm_fileparts( trans.native.Vo.fname ); 
    
  % in case of SPM input segmentation we have to add the name here to have a clearly different naming of the CAT output 
  if isfield(job,'spmpp'), nam = ['c1' nam]; end
  
  % definition of subfolders
  if job.extopts.subfolders
    labelfolder   = 'label';
  else
    labelfolder   = '';
  end

  vx_vol  = sqrt( sum( trans.native.Vi.mat(1:3,1:3).^2 ) ); % voxel size of the processed image
  vx_voln = sqrt( sum( trans.affine.mat(1:3,1:3).^2 ) );    % voxel size of the registration space

  % print progress
  stime = cat_io_cmd('ROI estimation'); if job.extopts.verb, fprintf('\n'); end 

  
  % get atlases maps that should be evaluated
  FAF = job.extopts.atlas; 
  FA  = {}; fai = 1;
  AN  = fieldnames(job.output.atlases);
  for ai = 1:numel(AN)
    fafi = find(cellfun('isempty',strfind(FAF(:,1),[AN{ai} '.']))==0,1); 
    if ~isempty(fafi) && (isempty(job.output.atlases.(AN{ai})) || job.output.atlases.(AN{ai})) 
      FA(fai,:) = FAF(fafi,:);  %#ok<AGROW>
      fai = fai+1; 
    end
  end
  if isempty(FA)
    % deactivate output
    FN = job.output.atlas; 
    for ai = 1:numel(AN)
      job.output.atlas.(FN{ai}) = 0; 
    end
  else
    % get atlas resolution 
    % we sort the atlases to reduce data resampling
    VA = spm_vol(char(FA(:,1))); 
    for ai=1:numel(VA), VAvx_vol(ai,:) = sqrt(sum(VA(ai).mat(1:3,1:3).^2)); end  %#ok<AGROW>
    [VAs,VAi] = sortrows(VAvx_vol);  %#ok<ASGLU>
    FA = FA(VAi,:); VA = VA(VAi,:); VAvx_vol = VAvx_vol(VAi,:); 
  end
  
  %stime2= clock;
  for ai=1:size(FA,1)
    %%
    if ai==1 || any(VAvx_vol(ai,:)~=VAvx_vol(ai-1,:)) || any(VA(ai).dim~=VA(ai-1).dim)
      % resample data in atlas resolution for the first time or if the atlas resolution changes
      
      % map data to actual template space
      if ai==1
        stime2  = cat_io_cmd('  Data mapping to normalized atlas space','g5','', job.extopts.verb-1); 
      else
        stime2  = cat_io_cmd('  Data mapping to normalized atlas space','g5','', job.extopts.verb-1,stime2); 
      end
      
      %wVv = spm_vol(FA{ai,1});   % atlas volume information 
      transw = rmfield(trans.warped,'yx');     % dartel/shooting deformation data
      transw.odim = VA(ai).dim;
      transw.M1   = VA(ai).mat;
      mati        = spm_imatrix( trans.affine.mat - VA(ai).mat) ; 
      vdim        = spm_imatrix( VA(ai).mat ); % trans.affine.mat ); 
      matit       = mati(1:3) ./ vdim(7:9); %./vx_vol; 
      for i=1:3, transw.y(:,:,:,i) = transw.y(:,:,:,i) * job.extopts.vox ./ VAvx_vol(ai,i);  end
      transw.y    = cat(4,transw.y(:,:,:,1) + matit(1), transw.y(:,:,:,2) + matit(2), transw.y(:,:,:,3) + matit(3) );
    end
    
    %% haveing the new map you can project the segmentation
    wYp0     = cat_vol_ROInorm(Yp0 ,transw,1,0,FA);
    wYcls    = cat_vol_ROInorm(Ycls,transw,1,1,FA);
    wYa      = cat_vol_ROInorm([],[],ai,0,FA);

    if debug % this does not work in case of trimmed atlases set do not include the full brain 
      fprintf('\n%8s %8s %8s\n','GM','WM','CSF') 
      fprintf('%8.2f %8.2f %8.2f\n', [cat_stat_nansum(wYcls{1}(:)),cat_stat_nansum(wYcls{2}(:)),cat_stat_nansum(wYcls{3}(:))] / 1000); 
      fprintf('%8.2f %8.2f %8.2f\n', [cat_stat_nansum(Ycls{1}(:)) ,cat_stat_nansum(Ycls{2}(:)) ,cat_stat_nansum(Ycls{3}(:)) ] * prod(vx_vol) / 1000 / 255); 
      fprintf('\n');
    end      
 
    [px,atlas] = fileparts(FA{ai,1}); clear px; %#ok<ASGLU>
    stime2 = cat_io_cmd(sprintf('  ROI estimation of ''%s'' atlas',atlas),'g5','', job.extopts.verb-1,stime2);

    % map atlas to actual template space 
    %wYa   = cat_vol_ROInorm([],[],ai,0,FA);

    if job.extopts.WMHC == 3   
      FA{ai,3} = unique( [FA{ai,3} {'wmh'}] ); %#ok<AGROW>
    end
    
    %% extract ROI data
    csv   = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],FA{ai,3},FA);  % volume
    
    % thickness
    if exist('Yth1','var') 
    % For thickness we want to avoid values in non-cortical regions such as 
    % the ventricles or regions with relative low GM volume or high CSF volume. 
      csv  = cat_vol_ROIestimate(wYp0,wYa,wYth1,ai,'T',csv,{'gm'},job.extopts.atlas); %.*wYmim
      % correct for ventricular regions that use the 'Ven' keyword.
      ven  = find(cellfun('isempty',strfind( csv(:,2) , 'Ven'))==0); 
      csv(ven,end) = {nan};  %#ok<FNDSB>
      % correct for regions with relative low GM (<10%) or high CSF volume (>50%). 
      csvf = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],{'csf','gm','wm'},FA);
      vola = [nan,nan,nan;cell2mat(csvf(2:end,3:end))]; 
      volr = vola ./ repmat(sum(vola,2),1,3); 
      csv(volr(:,2)<0.1 | volr(:,2)>0.5,end) = {nan}; 
    end

    % xml-export one file for all (this is a structure)
    ROI.(atlas) = csv;

  end 
  if exist('ROI','var')
	  % write results
	  catROI = cat_roi_fun('csvtab2xmlroi',ROI);
	  cat_io_xml(fullfile(pth,labelfolder,['catROI_' nam '.xml']),catROI,'write');
	  
	  cat_io_cmd(' ','g5','',job.extopts.verb,stime2);
	  cat_io_cmd('','n','',1,stime);
  end
  
return

%=======================================================================
function wYv = cat_vol_ROInorm(Yv,warped,ai,mod,FA)
% ----------------------------------------------------------------------
% normalized space:  
% ----------------------------------------------------------------------
% for normalized space no further adaptions are available, but 
% a masking based on the tissue map can be used
% ----------------------------------------------------------------------
 
  % load mask (and complete undefined parts)
  if isempty(Yv)
    % no input - load atlas
    [pp,ff,ee] = spm_fileparts(FA{ai,1});
    FAai1 = fullfile(pp,[ff ee]); 
    if ~exist(FAai1,'file')
      error('cat:cat_main:missAtlas','Miss cat atlas-file ''%s''!',FA{ai,1});
    end
    % try multiple times, because of read error in parallel processing
    for i=1:5
      try
        wVv = spm_vol(FA{ai,1});
        wYv = spm_read_vols(wVv);
        break
      catch 
        pause(0.5)
      end
    end

    if 0
      % resample atlas, if the atlas resolution differs from the actual template resolution
      for i=1:5
        try
        %if any(wVv.mat(:) ~= warped.nmat(:)) || any( wVv.dim ~= warped.ndim )
        % wVv2 = wVv; wVv2.mat = warped.nmat; wVv2.dim = warped.ndim;  
        %  [t,wYv] = cat_vol_imcalc([wVv2,wVv],wVv2,'i2',struct('interp',0,'verb',0));
        if any(wVv.mat(:) ~= warped.omat(:)) || any( wVv.dim ~= warped.odim )
          wVv2 = wVv; wVv2.mat = warped.omat; wVv2.dim = warped.odim;  
          [t,wYv] = cat_vol_imcalc([wVv2,wVv],wVv2,'i2',struct('interp',0,'verb',0));
        else
          wYv = spm_read_vols(wVv);
        end
        catch
          pause(0.5)
        end
      end
    end
    wYv = cat_vol_ctype(wYv,wVv(1).private.dat.dtype);
  else
    % map image to atlas space
    %for yi=1:numel(warped.ress), warped.y(:,:,:,yi) = warped.y(:,:,:,yi) * warped.ress(yi); end
    if mod==0
      old = 0;
      
      %Vlai = spm_vol(FA{ai,1});
      %vx_vol_Vlai   = sqrt(sum(Vlai.mat(1:3,1:3).^2));
      %vx_vol_Vdef   = sqrt(sum(warped.M1(1:3,1:3).^2));
      
      if old % this did not work for the thickness map?  
        [wYv,w] = spm_diffeo('push',Yv,warped.y,warped.odim(1:3)); spm_field('boundary',1);
        wYv = spm_field(w,wYv,[sqrt(sum(warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]);
      elseif 0 %any( vx_vol_Vdef > vx_vol_Vlai*2)
        try
          %% increase resolution - this caused many memory error messages (RD201911) 
          fc   = ceil(vx_vol_Vdef / vx_vol_Vlai);
          ddim = (size(warped.y)-1)*fc+1; ddim(4)=[]; 
          eyev = eye(4); eyev(1:end-1) = eyev(1:end-1) * 1/fc;
          Yy   = zeros([ddim 3],'single');                        
          for k1=1:3
            for i=1:ddim(3)
              Yy(:,:,i,k1) = single(spm_slice_vol(warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]),ddim(1:2),[1,NaN])); % adapt for res
            end
          end
          Yvi  = zeros(ddim,'single'); 
          for i=1:ddim(3)
            Yvi(:,:,i) = single(spm_slice_vol(Yv(:,:,:),eyev*spm_matrix([0 0 i]),ddim(1:2),[1,NaN])); % adapt for res
          end

          [wYv,w]  = spm_diffeo('push',Yvi,Yy,warped.odim(1:3));
          % divide by jacdet to get unmodulated data
          wYv = wYv./(w+0.001); 
        catch
          cat_io_cprintf('warn','\n  Possible memory problems use default push operation.'); 
          [wYv,w]  = spm_diffeo('push',Yv,warped.y,warped.odim(1:3));
          % divide by jacdet to get unmodulated data
          wYv = wYv./(w+0.001); 
        end
      else      
        [wYv,w]  = spm_diffeo('push',Yv,warped.y,warped.odim(1:3));
        % divide by jacdet to get unmodulated data
        wYv = wYv./(w+0.001); 
      end
    elseif mod==1 && iscell(Yv) % tissue case
      nicemapping = 1;
      if nicemapping
        % Modulation using spm_diffeo and push introduces aliasing artifacts,
        % thus we use the def2det function of the inverted deformations to obtain the old and 
        % in my view a more appropriate jacobian determinant 
        % The 2nd reason to use the old modulation is compatibility with cat_vol_defs.m
        Yy = spm_diffeo('invdef',warped.y,warped.odim,eye(4),warped.M0);
        w  = spm_diffeo('def2det',Yy)/det(warped.M0(1:3,1:3)); clear Yy;
       % w  = w * mean(warped.ress(yi)); 
        % ensure that jacobian det is positive (no clue why some times the sign is switched)
        if mean(w(~isnan(w))) < 0, w = -w; end 
        w(:,:,[1 end]) = NaN; w(:,[1 end],:) = NaN; w([1 end],:,:) = NaN;
      end
      
      wYv = cell(1,numel(Yv));
      for i=1:numel(Yv)
        if nicemapping 
          [wYv{i},w2] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3)); 
          % divide by jacdet to get unmodulated data
          wYv{i} = wYv{i}./(w2+0.001); 
          wYv{i} = wYv{i} .* w;
        else
          % simple push
          [wYv{i},w] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3)); 
        end
      end
    else
      error('unknown case');
    end
    %spm_smooth(wYv,wYv,[1 1 1]); 
  end
  
  
return
%=======================================================================

%=======================================================================
function csv = cat_vol_ROIestimate(Yp0,Ya,Yv,ai,name,csv,tissue,FA)
% ----------------------------------------------------------------------
% estimate values
% ----------------------------------------------------------------------


% load atlas-csv-file

  [pp,ff] = fileparts(FA{ai,1});
  csvf = fullfile(pp,[ff '.csv']);

  if isempty(csv) 
    if exist(csvf,'file')
      csv = cat_io_csv(csvf,'','',struct('delimiter',';')); 
    else
      IDs = unique(Ya); 
      cat_io_cprintf('warn',sprintf('\n    Cannot find ''%s'' csv-file with region names! ',ff)); 
      %ROIid;ROIabbr;ROIname;ROIbaseid;ROIbasename;Voxel;Volume;XYZ
      csv = [{'ROIid','ROIabbr','ROIname'}; ...
        num2cell(IDs) ...
        cellstr([repmat('ROI',numel(IDs)) num2str(IDs,'%03d')]) ...
        cellstr([repmat('ROI',numel(IDs)) num2str(IDs,'%03d')])];
    end
    
    % remove empty rows and prepare structure names
    if size(csv,2)>2, csv(:,3:end)=[]; end
    for ri=size(csv,1):-1:1
      if isempty(csv{ri,1}) || isempty(csv{ri,2}) 
        csv(ri,:)=[];
      elseif csv{ri,1}==0
        csv(ri,:)=[];
      end       
    end
  end
  name = genvarname(strrep(strrep(name,'-','_'),' ','_'));
  
  
  
  %% volume case
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));   
  % other maps with masks
  for ti=1:numel(tissue)
    switch name(1)
      case 'V' % volume
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
        for ri=2:size(csv,1)
          switch lower(tissue{ti})
            case 'csf',    Ymm = single(Yv{3}) .* single(Ya==csv{ri,1});
            case 'gm',     Ymm = single(Yv{1}) .* single(Ya==csv{ri,1});
            case 'wm',     Ymm = single(Yv{2}) .* single(Ya==csv{ri,1});
            case 'wmh',    Ymm = single(Yv{7}) .* single(Ya==csv{ri,1}); 
            case 'brain',  Ymm = single(Yv{1} + Yv{2} + Yv{3} + Yv{7}) .* single(Ya==csv{ri,1});
            case 'tissue', Ymm = single(        Yv{2} + Yv{3} + Yv{7}) .* single(Ya==csv{ri,1});
            case '',       Ymm = single(Ya==csv{ri,1});
          end
          csv{ri,end} = 1/1000 * cat_stat_nansum(Ymm(:));
        end
      case 'c'
        return
      otherwise % 
        csv{1,end+1} = strrep([name tissue{ti}],'Tgm','ct');  %#ok<AGROW>
        switch lower(tissue{ti})
          case 'csf',    Ymm = Yp0toC(Yp0,1); 
          case 'gm',     Ymm = Yp0toC(Yp0,2); 
          case 'wm',     Ymm = Yp0toC(Yp0,3); 
          case 'wmh',    Ymm = Yp0toC(Yp0,4); 
          case 'brain',  Ymm = Yp0>0.5;
          case 'tissue', Ymm = Yp0>1.5;
          case '',       Ymm = true(size(Yp0));
        end
        for ri=2:size(csv,1)
          csv{ri,end} = cat_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
        end
    end
  end
  
return