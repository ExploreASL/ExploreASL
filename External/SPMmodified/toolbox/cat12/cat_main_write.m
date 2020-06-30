function cat_warnings = cat_main_write(Ym,Ymi,Ycls,Yp0,Yl1,job,res,trans,cat_warnings)
% ______________________________________________________________________
% Write volumetric preprocessing results.
%  
%   cat_warnings = cat_main_write(Ym,Ycls,Yp0b,job,trans,cat_warnings)
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_write.m 1577 2020-03-09 17:36:03Z dahnke $


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  %VT  = res.image(1);  % denoised/interpolated n*.nii % RD 20190328
  VT0 = res.image0(1); % original 
  [pth,nam,ee,ie] = spm_fileparts(VT0.fname); 
  
  % in case of SPM input segmentation we have to add the name here to have a clearly different naming of the CAT output 
  if isfield(res,'spmpp'), nam = ['c1' nam]; VT0.fname = fullfile(pth,[nam,ee,ie]); end

  %d = VT.dim(1:3); % RD 20190328
  tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)]; 

  % definition of subfolders
  if job.extopts.subfolders
    mrifolder     = 'mri';
  else
    mrifolder     = '';
  end
  
  stime = cat_io_cmd('Write result maps');

  % bias, noise and global corrected without masking for subject space and with masking for other spaces 
  cat_io_writenii(VT0,Ym,mrifolder,'m', ...Dartel
    'bias and noise corrected, global intensity normalized','uint16',[0,0.0001], ... 
    min([1 0 2],[job.output.bias.native job.output.bias.warped job.output.bias.dartel]),trans);
  cat_io_writenii(VT0,Ym.*(Yp0>0.1),mrifolder,'m', ... 
    'bias and noise corrected, global intensity normalized (masked due to normalization)','uint16',[0,0.0001], ...
    min([0 1 0],[job.output.bias.native job.output.bias.warped job.output.bias.dartel]),trans);

  % bias, noise and local intensity corrected without masking for subject space and with masking for other spaces 
  cat_io_writenii(VT0,Ymi,mrifolder,'mi', ...
    'bias and noise corrected, local intensity normalized','uint16',[0,0.0001], ... 
    min([1 0 2],[job.output.las.native job.output.las.warped job.output.las.dartel]),trans);
  cat_io_writenii(VT0,Ymi.*(Yp0>0.1),mrifolder,'mi', ... 
    'bias and noise corrected, local intensity normalized (masked due to normalization)','uint16',[0,0.0001], ...
    min([0 1 0],[job.output.las.native job.output.las.warped job.output.las.dartel]),trans);

  % Yp0
  cat_io_writenii(VT0,Yp0,mrifolder,'p0','label map','uint8',[0,5/255],job.output.label,trans);
  
  % partitioning
  cat_io_writenii(VT0,Yl1,mrifolder,'a0','brain atlas map for major structures and sides',...
    'uint8',[0,1],job.output.atlas,trans);

  % class maps 1-3
  fn = {'GM','WM','CSF','head','head','background'};
  for clsi=1:3
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
      min([1 1 0 3],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
      job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'single',[0,1],...
      min([0 0 3 0],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
      job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
  end

  % write WMH class maps
  if job.extopts.WMHC>=3 && job.extopts.WMHCstr>0 && ~job.inv_weighting;
    cat_io_writenii(VT0,single(Ycls{7})/255,mrifolder,'p7','WMH tissue map','uint8',[0,1/255],...
      min([1 1 0 3],[job.output.WMH.native job.output.WMH.warped ...
      job.output.WMH.mod job.output.WMH.dartel]),trans); % 1 0 0 0
    cat_io_writenii(VT0,single(Ycls{7})/255,mrifolder,'p7','WMH tissue map','single',[0,1],...
      min([0 0 3 0],[job.output.WMH.native job.output.WMH.warped ...
      job.output.WMH.mod job.output.WMH.dartel]),trans); % 0 1 2 2


  elseif any([job.output.WMH.native job.output.WMH.warped ...
      job.output.WMH.mod job.output.WMH.dartel])
    % should write, but can not ...
    if job.extopts.WMHC<3
      cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_WMHC:output','Cannot write WMH images because no seperate class is used (set WMHC>=3).');
    elseif job.extopts.WMHCstr==0
      cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_WMHC:output','Cannot write WMH images because WMHCstr is 0.');
    elseif job.inv_weighting
      cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_WMHC:output','Cannot write WMH images because inverse contrast was detected.');
    end
    % 
  end 
  % export lesions
  if job.extopts.SLC>0 && numel(Ycls)>7
    cat_io_writenii(VT0,single(Ycls{8}),mrifolder,'p8','stroke lesion map','uint8',[0,1/255],...
      min([1 1 0 3],[job.output.SL.native job.output.SL.warped ...
      job.output.SL.mod job.output.SL.dartel]),trans); % 1 0 0 0
    cat_io_writenii(VT0,single(Ycls{8}),mrifolder,'p8','stroke lesion map','single',[0,1],...
      min([0 0 3 0],[job.output.WMH.native job.output.SL.warped ...
      job.output.SL.mod job.output.SL.dartel]),trans); % 1 0 0 0
  end

  % developer - intensity scaled tissue classe maps
  % ----------------------------------------------------------------------
  % The strong normalization of the T1 data can directly be used as tissue
  % segmentation. The Ymi images is scaled to get similar maps for each 
  % tissue class, with good visible differences in the sulci.
  job.output.intsegments = job.extopts.experimental;
  if job.output.intsegments
    if (any(tc(:)) || job.extopts.WMHC==3 && job.extopts.WMHCstr>0 && ~job.inv_weighting); 

      % intensity scaled tissue maps
      Yclsi = cell(1,3);
      for clsi=1:3
        clsid = [2 3 1];
        %Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*5/255;  % RD 20190328
        Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
        Yclsi{clsi} = Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),clsid(clsi));
        switch clsi
          case 1
            Yclsi{clsi} = Yclsi{clsi} .* (Ycls{clsi}>0) + ...
              Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),3) .* (Ycls{2}==0 & (Ycls{1}>0));
          case 2 
            Yclsi{clsi} = Yclsi{clsi} .* (Ycls{clsi}>0) +  ...
              Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2) .* (Ycls{2}==255 & ~cat_vol_morph(Ycls{2}>192,'e'));
            Ywmhp = cat_vol_ctype( Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2) .* cat_vol_morph(Ycls{2}>192,'e') * 255);
          case 3 
            Yclsi{clsi} = Yclsi{clsi} + ...
              (Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2)) .* (Ycls{1}==0 & (Ycls{3}>0)) + ...
              (Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),3)) .* (Ycls{2}==0 & (Ycls{3}>0));
        end
        Yclsi{clsi} = cat_vol_ctype(Yclsi{clsi} * 255);
      end
      clear Yp0; 

      % class maps 1-3
      % Yclss = single(Yclsi{1} + Yclsi{2} + Yclsi{3} + Ywmhp + Ywmh)/255; % + single(Ywmh)/255;
      for clsi=1:3
        cat_io_writenii(VT0,single(Yclsi{clsi})/255,mrifolder,sprintf('pi%d',clsi),...
          sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
          min([1 1 0 3],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
          job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
        cat_io_writenii(VT0,single(Yclsi{clsi})/255,mrifolder,sprintf('pi%d',clsi),...
          sprintf('%s tissue map',fn{clsi}),'single',[0,1],...
          min([0 0 2 0],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
          job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
      end
      clear Yclsi; 

      % write WMH class maps
      if job.extopts.WMHC>=3 && job.extopts.WMHCstr>0 && ~job.inv_weighting;
        cat_io_writenii(VT0,(single(Ywmhp) + single(Ycls{7}))/255,mrifolder,...
          'pi7','WMH tissue map','uint8',[0,1/255],...
          min([1 1 0 2],[job.output.WMH.native job.output.WMH.warped ...
          job.output.WMH.mod job.output.WMH.dartel]),trans); % 1 0 0 0
        cat_io_writenii(VT0,(single(Ywmhp) + single(Ycls{7}))/255,mrifolder,...
          'pi7','WMH tissue map','single',[0,1],...
          min([0 0 2 0],[job.output.WMH.native job.output.WMH.warped ...
          job.output.WMH.mod job.output.WMH.dartel]),trans); % 0 1 2 2
      end 
      clear Ywmhp;
    end
  else
    clear Yp0; 
  end
  % ----------------------------------------------------------------------


  % classe maps 4-6 (for full TPM/template creation, e.g. for apes)
  if any(cell2mat(struct2cell(job.output.TPMC)'))
    for clsi=4:6
      cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
        min([1 1 0 3],[job.output.TPMC.native job.output.TPMC.warped ...
        job.output.TPMC.mod job.output.TPMC.dartel]),trans);
      cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'single',[0,1],...
        min([0 0 3 0],[job.output.TPMC.native job.output.TPMC.warped ...
        job.output.TPMC.mod job.output.TPMC.dartel]),trans);
    end
  end
  %clear cls clsi fn Ycls; % we need this maps later for the ROIs

  % write jacobian determinant
  if job.output.jacobian.warped
    %%
    if res.do_dartel==2 % shooting
      dt = trans.jc.dt2; 
      dx = 10; % smaller values are more accurate, but large look better; 
      [D,I] = cat_vbdist(single(~(isnan(dt) | dt<0 | dt>100) )); D=min(1,D/min(dx,max(D(:)))); 
      dt = dt(I); dt = dt .* ((1-D) + D); dt(isnan(dt))=1; 
      %dt = 1/max(eps,dt); 
      clear y0 D I
    else %dartel
      [y0, dt] = spm_dartel_integrate(reshape(trans.jc.u,[trans.warped.odim(1:3) 1 3]),[1 0], 6);  %#ok<ASGLU>
      clear y0
    end
    N      = nifti;
    N.dat  = file_array(fullfile(pth,mrifolder,['wj_', nam, '.nii']),trans.warped.odim(1:3),...
               [spm_type('float32') spm_platform('bigend')],0,1,0);
    N.mat  = trans.warped.M1;
    N.mat0 = trans.warped.M1;
    N.descrip = ['Jacobian' VT0.descrip];
    create(N);
    N.dat(:,:,:) = dt;
  end


  %% write atlas output
  if any(cell2mat(struct2cell(job.output.atlas)'))
    % get atlases
    FAF = job.extopts.atlas; 
    FA  = {}; fai = 1;
    AN  = fieldnames(job.output.atlases);
    for ai = 1:numel(AN)
      fafi = find(cellfun('isempty',strfind(FAF(:,1),[AN{ai} '.']))==0);
      if ~isempty(fafi) && job.output.atlases.(AN{ai}), FA(fai,:) = FAF(fafi,:); fai = fai+1; end %#ok<AGROW>
    end
    
    for ai=1:size(FA,1)
      [px,atlas] = fileparts(FA{ai,1});  %#ok<ASGLU>

      % map atlas in native space
      Vlai = spm_vol(FA{ai,1});
      if any( Vlai.dim ~= trans.warped.odim )
        % interpolation
        Vlai = spm_vol(FA{ai,1});
        yn = numel(trans.warped.y); 
        p  = ones([4,yn/3],'single'); 
        p(1,:) = trans.warped.y(1:yn/3);
        p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
        p(3,:) = trans.warped.y(yn/3*2+1:yn);
        amat   = Vlai.mat \ trans.warped.M1; 
        p      = amat(1:3,:) * p;

        Yy = zeros([res.image(1).dim(1:3),3],'single'); 
        Yy(1:yn/3)        = p(1,:);
        Yy(yn/3+1:yn/3*2) = p(2,:);
        Yy(yn/3*2+1:yn)   = p(3,:);

        Yy = double(Yy); 
      else
        Yy = double(trans.warped.y);
      end
      Ylai = cat_vol_ctype(spm_sample_vol(Vlai,Yy(:,:,:,1),Yy(:,:,:,2),Yy(:,:,:,3),0));
      Ylai = reshape(Ylai(:),trans.native.Vi.dim); 
      if ~debug, clear Yy; end

      % write map (mri as tissue subforder and mri_atals as ROI subfolder)
      if isempty(mrifolder), amrifolder = ''; else amrifolder = 'mri_atlas'; end
      cat_io_writenii(VT0,Ylai,amrifolder,[atlas '_'],[atlas ' original'],...
        'uint8',[0,1],job.output.atlas,trans);
      if ~debug, clear Vlai Ylai; end
    end
  end


  %% deformations y - dartel > subject
  if job.output.warps(1)
      Yy        = spm_diffeo('invdef',trans.warped.y,trans.warped.odim,eye(4),trans.warped.M0);
      N         = nifti;
      N.dat     = file_array(fullfile(pth,mrifolder,['y_', nam, '.nii']),[trans.warped.odim(1:3),1,3],'float32',0,1,0);
      N.mat     = trans.warped.M1;
      N.mat0    = trans.warped.M1;
      N.descrip = 'Deformation';
      create(N);
      N.dat(:,:,:,:,:) = reshape(Yy,[trans.warped.odim,1,3]);
      if ~debug, clear Yy; end
  end

  %% deformation iy - normalized > subject 
  if job.output.warps(2) 
    if any(trans.native.Vo.dim~=trans.native.Vi.dim)
      %%
      vx_voli  = sqrt(sum(trans.native.Vi.mat(1:3,1:3).^2));  
      vx_volo  = sqrt(sum(trans.native.Vo.mat(1:3,1:3).^2));
      eyev = eye(4); eyev([1 6 11]) = eyev([1 6 11]) .* vx_volo./vx_voli; 
      Yy2  = zeros([trans.native.Vo.dim 1 3],'single');                        
      for k1=1:3
        for i=1:trans.native.Vo.dim(3),
          Yy2(:,:,i,:,k1) = trans.warped.M1(k1,4) + trans.warped.M1(k1,k1) * ...
            single(spm_slice_vol(trans.warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]), ...
            trans.native.Vo.dim(1:2),[1,NaN])); % adapt for res
        end
      end
    else 
      yn = numel(trans.warped.y); 
      p  = ones([4,yn/3],'single'); 
      p(1,:) = trans.warped.y(1:yn/3);
      p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
      p(3,:) = trans.warped.y(yn/3*2+1:yn);
      p      = trans.warped.M1(1:3,:) * p;

      Yy2 = zeros([res.image(1).dim(1:3),1,3],'single'); 
      Yy2(1:yn/3)        = p(1,:);
      Yy2(yn/3+1:yn/3*2) = p(2,:);
      Yy2(yn/3*2+1:yn)   = p(3,:);
    end
    clear p; 

    % f2 = spm_diffeo('resize', f1, dim)
    % write new output
    Ndef      = nifti;
    Ndef.dat  = file_array(fullfile(pth,mrifolder,['iy_', nam, '.nii']),[res.image0(1).dim(1:3),1,3],...
                [spm_type('float32') spm_platform('bigend')],0,1,0);
    Ndef.mat  = res.image0(1).mat;
    Ndef.mat0 = res.image0(1).mat;
    Ndef.descrip = 'Inverse Deformation';
    create(Ndef);
    Ndef.dat(:,:,:,:,:) = Yy2;

    if ~debug, clear Yy2; end
  end
  fprintf('%5.0fs\n',etime(clock,stime));




end