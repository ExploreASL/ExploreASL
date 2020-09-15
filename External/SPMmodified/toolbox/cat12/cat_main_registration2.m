function [trans,reg] = cat_main_registration2(job,res,Ycls,Yy,tpmM,Ylesion)
% ______________________________________________________________________
%  Spatial registration function of cat_main preprocessing that include
%  the SPM DARTEL and (optimized) SHOOTING registration approaches. 
%
%  There are 4 image spaces with image different dimensions and resolution: 
%    (1) individual/original  (IR)  - properties of the (interpolated) anatomical image 
%    (2) template             (TR)  - properties of the registration template image
%    (3) registration         (RR)  - properties for the registration 
%    (4) output/analyse       (AR)  - final resolution of the normalized images
%
%  DARTEL runs depending on the output resolution and RR = AR, whereas 
%  the original DARTEL runs always on the TR. It takes typically about  
%  3 Minutes/Subject for RR = 1.5 mm. 
%  SHOOTING needs about 10 Minutes/Subject for RR = 1.5 mm and much more 
%  for higher RR and especially the low smooth deformation are important. 
%  The optimized version used therefore reduced resolution level to speedup
%  the iteration and allow more smooth deformations and reduce the final 
%  costs. The changes between the deformation and the matching compared to 
%  the templates are used as iteration criteria.
%
%  Main control parameter:
%    job.extopts.regstr ..
%      * Main cases:
%        0  .. DARTEL, 
%        eps - 1  .. Optimized Shooting with low (eps; fast) to high
%                    to high quality (1; slow) with 0.5 as default
%        2  .. Optimized Shooting with fixed resolutions (3:(3-TR)/4:TR)
%        3  .. Optimized Shooting with fixed resolutions (TR/2:TR/4:TR)
%        4  .. Default Shooting
% 
%      * Fixed resolution level with interpolation:
%        11 .. hard deformations
%        12 .. medium deformations
%        13 .. soft deformations
%
%      * Fixed resolution level without interpolation (max res = TR):
%        21 .. hard deformations
%        22 .. medium deformations
%        23 .. soft deformations
%
%  Structure:
%    [trans,reg] = cat_main_registration(job,res,Ycls,Yy,tpmM)
%   
%    trans .. output variable that is used in cat_main and cat_io_writenii
%    reg   .. output variable that include the summarize of the deformation
% 
%    job   .. SPM job structure with cat_main parameter
%    res   .. SPM parameter structure with further fields from CAT processing
%    Ycls  .. tissue classification that is used for deformation
%    Yy    .. old initial SPM deformation 
%    tpmM  .. template files
% ______________________________________________________________________
%
%  The TR is expected to be between 1.0 and 1.5 mm. Higher resolution are 
%  possible, but not will not be become a standard CAT templates.
%
%  RRs higher than the TR allow small improvements, but with improper costs.
%  Tests with 0.5 mm RR lead to an shooting error (to high determinant).
%
%  Lower RR lead to smoother deformations and were the focus of the
%  optimization and allow further deformation cases.
%
%  There are different ways to use lower resolutions: 
%    (LR = lowest resolution - should be better than 3 mm!):
%    1) flexible step size with fixed lower limit: 
%       a) LR : (LR - TR)/5 : TR 
%    2) upper limit
%       a) fixed levels without/without limit by TR oder by user (deformation limit)
%          max( RR , 3.0 : -0.5 : 1.0 )
%          max( RR , 2.5 : -0.5 : 0.5 )
%       b) dynamic levels
%    3) TR limit
%       a)  TR+SS*5 : -SS : TR    with SS = 0.25 - 0.5
%    4) LR limit
%
%  Although, it is possible to use the deformation to estimate finer
%  rigid/affine transformation, we do not do this because rigid/affine 
%  output should run without spatial registration.
%
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_registration2.m 1574 2020-03-02 23:55:15Z gaser $

  if ~nargin, help cat_main_registration; return; end

  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_main_registration'); debug = 1; break; end; end

  
  % output
  trans = struct();
  reg   = struct();
  def.nits          = 64; 
  def.opt.ll1th     = 0.001;    % smaller better/slower
  def.opt.ll3th     = 0.002; 
  def.fast          = 0;        % no report, no output
  def.affreg        = 1;        % additional affine registration
  def.iterlim       = inf;      % limit inner iterations
  def.rigidShooting = 0;        % rigid vs. affine Shooting 
  if ~isfield(job.extopts,'reg'), job.extopts.reg = struct(); end
  dreg = cat_io_checkinopt(job.extopts.reg,def); 
  
  
  % error in parallel processing - can't find shooting files (2016/12)
  % switch to shooting directory
  if ~exist('spm_shoot_defaults','file') || ~exist('spm_shoot_update','file') 
    olddir = cd; cd(fullfile(spm('dir'),'toolbox','Shoot'));
  end


  % this is just for me to create different templates and can be removed in a final version
  fast   = dreg.iterlim;                   % limit iterations per template level to test if processing work in principle 
  if numel(job.extopts.vox)>1 || numel(job.extopts.regstr)>1 || (isfield(job,'export') && job.export)
    job.extopts.multigreg = 1;
    export = 1;  % write files in sub-directories
  else
    job.extopts.multigreg = 0;
    export = 0; 
  end

  % additional affine registration of the GM segment
  if dreg.affreg
    %%
    obj.image         = res.image; 
    obj.image.pinfo   = repmat([255;0],1,size(Ycls{1},3));
    obj.image.dt      = [spm_type('UINT8') spm_platform('bigend')];
    obj.image.dat     = cat_vol_ctype(single(Ycls{1})/2 + single(Ycls{2}));  
    if isfield(obj.image,'private'), obj.image = rmfield(obj.image,'private'); end
    obj.samp = 1.5; 
    obj.fwhm = 1;
    if isfield(obj,'tpm'), obj = rmfield(obj,'tpm'); end
    obj.tpm = spm_load_priors8(res.tpm(1:2));
    Affine  = spm_maff8(obj.image,obj.samp,obj.fwhm,obj.tpm,res.Affine,job.opts.affreg,80);
    res.Affine = Affine;
    spm_progress_bar('Clear');
  end

  do_req = res.do_dartel; 
  % this is the main loop for different parameter
  for regstri = numel(job.extopts.regstr):-1:1;
    for voxi = numel(job.extopts.vox):-1:1; 
      if (numel(job.extopts.regstr) || numel(job.extopts.vox)) && job.extopts.verb, fprintf('\n\n'); end 
      
      
      %% set dartel/shooting templates
      if job.extopts.regstr(regstri)==0
        job.extopts.templates = job.extopts.darteltpms;
      else
        job.extopts.templates = job.extopts.shootingtpms;
      end
      res.tpm2 = cell(1,numel(job.extopts.templates)); 
      Vtmp = spm_vol(job.extopts.templates{1}); tmpM = Vtmp(1).mat; 
      n1 = max(min(2,numel(Vtmp)),numel(Vtmp)-1);  % use GM and WM for shooting
     
      %% registration main parameter
      lowres                     = 2.5;                   % lowest resolution .. best between 2 and 3 mm 
      tpmres                     = abs(tpmM(1));          % TPM resolution 
      tempres                    = abs(tmpM(1));          % template resolution 
      reg(regstri).opt.nits      = dreg.nits;             % registration interation (shooting default = 24)
      reg(regstri).opt.vxreg     = tpmres;                % regularisation parameter that original depend on the template resolution
      reg(regstri).opt.rres      = tempres;               % final registration resolution 
      reg(regstri).opt.stepsize  = (lowres - reg(regstri).opt.rres)/4;  % stepsize of reduction 
      reg(regstri).opt.resfac    = (lowres : -reg(regstri).opt.stepsize : reg(regstri).opt.rres) / reg(regstri).opt.rres; % reduction factor 
      reg(regstri).opt.ll1th     = dreg.opt.ll1th;                 % smaller better/slower
      reg(regstri).opt.ll3th     = dreg.opt.ll3th;                 % smaller better/slower 
      reg(regstri).opt.regstr    = job.extopts.regstr;   

    
      if job.extopts.regstr(regstri)==0
      % Dartel
        res.do_dartel            = 1; 
        reg(regstri).opt.rres    = tempres; %job.extopts.vox(voxi);
      elseif job.extopts.regstr(regstri)>0 && job.extopts.regstr(regstri)<=1
      % Optimized Shooting - manual limit 
        reg(regstri).opt.stepsize  = (lowres - tempres)/4;  % stepsize of reduction 
        reg(regstri).opt.resfac    = (lowres : -reg(regstri).opt.stepsize : tempres) / tempres; % reduction factor 

        reg(regstri).opt.ll1th     = 0.0010 + 0.10*(1-job.extopts.regstr(regstri));   % smaller better/slower
        reg(regstri).opt.ll3th     = 0.0001 + 0.10*(1-job.extopts.regstr(regstri));   % smaller better/slower 
      elseif job.extopts.regstr(regstri)==4 
      % Default Shooting  
        reg(regstri).opt.rres        = tempres;           % registration resolution depending on template resolution 
        reg(regstri).opt.stepsize    = 0;                 % stepsize of reduction 
        reg(regstri).opt.nits        = 24;                % Dartel default interation number
        reg(regstri).opt.resfac      = ones(1,5);         % reduction factor 
        reg(regstri).opt.ll1th       = 0;                 % smaller better/slower
        reg(regstri).opt.ll3th       = 0;                 % smaller better/slower 
      elseif job.extopts.regstr(regstri)==5 % this may not work because you finally need another interpolation! 
      % based on vox  
        tempres                      = job.extopts.vox(voxi);   % template resolution 
        reg(regstri).opt.rres        = job.extopts.vox(voxi);   % registration resolution depending on template resolution 
        reg(regstri).opt.stepsize    = (tempres*2 - tempres)/4;  % stepsize of reduction 
        reg(regstri).opt.resfac      = (tempres*2 : -reg(regstri).opt.stepsize : job.extopts.vox(voxi)) / job.extopts.vox(voxi); % reduction factor 
      elseif job.extopts.regstr(regstri)==2 
      % Optimized Shooting - manual limit 
        reg(regstri).opt.stepsize    = (lowres - tempres)/4;  % stepsize of reduction 
        reg(regstri).opt.resfac      = (lowres : -reg(regstri).opt.stepsize : tempres) / tempres; % reduction factor 
      elseif job.extopts.regstr(regstri)==3
      % Optimized Shooting - dynamic limit (depending on template resolution)  
        reg(regstri).opt.stepsize    = (tempres/2 - tempres)/4;  % stepsize of reduction 
        reg(regstri).opt.resfac      = (tempres/2 : -reg(regstri).opt.stepsize : tempres) / tempres; % reduction factor 
      elseif job.extopts.regstr(regstri)==6
        %%
        reg(regstri).opt.nits      = 64;;                   % registration interation (shooting default = 24)
        reg(regstri).opt.ll1th     = 0.0001;                 % smaller better/slower
        reg(regstri).opt.ll3th     = 0.0005;                 % smaller better/slower 
        %reg(regstri).opt.resfac    =  [4 3 2 1.5 1] ; 
        %reg(regstri).opt.resfac    =  [1 1 1 1 1]*1.5; 
        %reg(regstri).opt.resfac    =  [3 2 2.5/1.5 2/1.5 1]; 
        reg(regstri).opt.resfac    =  [2 2 2 2 2]; 
        %reg(regstri).opt.resfac    =  [4 4 4 3 2]; % n?
        %reg(regstri).opt.resfac    =  [4 4 3 2 1]; % 
        %reg(regstri).opt.resfac    =  [4 4 4 4 2]; % 
        %reg(regstri).opt.resfac    =  [4 4 4 4 4];
        %reg(regstri).opt.resfac    =  [5 5 5 5 5]; % not good
        %reg(regstri).opt.resfac    =  [6 6 6 6 6]; % not good
        %reg(regstri).opt.resfac    =  [3 2.5 2 1.5 1]; 
        %reg(regstri).opt.resfac    =  [2 1.75 1.5 1.25 1]; 
        %reg(regstri).opt.stepsize  = (lowres - reg(regstri).opt.rres)/4;  % stepsize of reduction 
        reg(regstri).opt.rres      =  reg(regstri).opt.resfac(end) * tempres;
        %reg(regstri).opt.resfac    =  max(reg(regstri).opt.resfac,reg(regstri).opt.rres)  / reg(regstri).opt.rres; 
        reg(regstri).opt.stepsize  =  max( abs( diff( reg(regstri).opt.resfac ) ) );  
        reg(regstri).opt.resfac    =  reg(regstri).opt.resfac  / reg(regstri).opt.resfac(end);  
      else
        % futher test cases
        highres = 1.0; % 0.5
        switch job.extopts.regstr(regstri)
        % -----------------------------------------------------------------
        % absolute fixed resolutions and reduction
        % -----------------------------------------------------------------
        % This allows identical smooth iterations for all levels and some
        % kind of frequency/deformation limit. 
        % The resolution levels are 1.0:0.5:3.0 mm, but can be changed to 
        % to 0.5:0.5:2.5 mm.
        % default = 12 | 22, expert = 11:13 | 21:23
        % -----------------------------------------------------------------
          case {11,12,13,14,15,16,17} 
            % independent of the TR and therefore can include interpolation
            reg(regstri).opt.rres     = 1.0 + 0.5 * (job.extopts.regstr(regstri) - 11); 
            reg(regstri).opt.stepsize = 0.5; 
            reg(regstri).opt.ll1th    = 0.005 * reg(regstri).opt.rres;                 % smaller better/slower
            reg(regstri).opt.ll3th    = 0.010 * reg(regstri).opt.rres;                 % smaller better/slower 
            reg(regstri).opt.resfac   = max( highres + 4*reg(regstri).opt.stepsize : -reg(regstri).opt.stepsize : highres , reg(regstri).opt.rres ) / reg(regstri).opt.rres;   
          case {21,22,23,24,25}
            % dependent on the TR without interpolation interpolation
            reg(regstri).opt.rres     = max(tempres,1.0 + 0.5 * (job.extopts.regstr(regstri) - 21)); 
            reg(regstri).opt.stepsize = 0.5; 
            reg(regstri).opt.ll1th    = 0.005 * reg(regstri).opt.rres;                 % smaller better/slower
            reg(regstri).opt.ll3th    = 0.010 * reg(regstri).opt.rres;                 % smaller better/slower 
            reg(regstri).opt.resfac   = max( highres + 4*reg(regstri).opt.stepsize : -reg(regstri).opt.stepsize : highres , reg(regstri).opt.rres ) / reg(regstri).opt.rres;   
        % -----------------------------------------------------------------
        % There are further cases, but they all have some drawbacks:
        % * Using a fixed stepsize with different final resolutions run 
        %   into the problematic low resolution (<3 mm).
        % * Using a additive or multiplicative template depending reduction
        %   will be bad for the GUI and it is better to have some fixed
        %   levels.
        % -----------------------------------------------------------------
          otherwise
             error('cat_main_registration:incorrectparameter','Incorrect value of "regres".\n');
        end
        
        % not required from a theoretic point of view but ...
        reg(regstri).opt.ll1th = reg(regstri).opt.ll1th * tempres/tpmres; 
        reg(regstri).opt.ll3th = reg(regstri).opt.ll3th * tempres/tpmres; 
      end
      
      
      %% manual setting of shooting parameter  
      if 0
        job.extopts.vox(voxi)          = 1.5;                   %#ok<UNRCH> % output resolution ...
        reg(regstri).opt.stepsize      = max(eps,0.5);          % stepsize of reduction
        reg(regstri).opt.rres          = 1.0;                   % registration resolution expert parameter ...
        reg(regstri).opt.vxreg         = 1.5;                   % regularisation parameter that original depend on the template resolution
        reg(regstri).opt.nits          = 5;                     % registration interation (shooting default = 24)
        reg(regstri).opt.ll1th         = 0.01;                  % smaller better/slower
        reg(regstri).opt.ll3th         = 0.04;                  % smaller better/slower 
        res.do_dartel                  = 2;                     % method, 1 - Dartel, 2 - Shooting
        fast                           = 10;                    % inner iterations limit
      end
     

      %% set default dartel/shooting templates in debug mode 
      if 0 %debug
        % only in case of the default templates
        job.extopts.templates = templates; 
        if job.extopts.expertgui==2 && ...
           (res.do_dartel==1 && job.extopts.regstr(regstri)>0) || ...
           (res.do_dartel==2 && job.extopts.regstr(regstri)==0)
          if res.do_dartel==2 && job.extopts.regstr(regstri)==0
            cat_io_cprintf('warn','Switch to default Dartel Template.\n');
            job.extopts.templates      = cat_vol_findfiles(fullfile(spm('dir'),'toolbox','cat12','templates_volumes'),'Template_*_IXI555_MNI152.nii'); 
            job.extopts.templates(end) = []; 
            reg(regstri).opt.rres = job.extopts.vox(voxi); 
          elseif res.do_dartel==1 && job.extopts.regstr(regstri)>0
            cat_io_cprintf('warn','Switch to default Shooting Template.\n');
            job.extopts.templates = cat_vol_findfiles(fullfile(spm('dir'),'toolbox','cat12','templates_volumes'),'Template_*_IXI555_MNI152_GS.nii',struct('depth',1)); 
          end
        end
        res.tpm2 = cell(1,numel(job.extopts.templates)); 
      end
      run2 = struct(); 
      for j=1:numel(res.tpm2)
        for i=1:n1, run2(i).tpm = sprintf('%s,%d',job.extopts.templates{j},i);end
        res.tpm2{j} = spm_vol(char(cat(1,run2(:).tpm)));
      end


      
      
      %% resolutions:
      %tpmres = abs(tpmM(1));                                                    % template resolution 
      tmpres = abs(tmpM(1));                                                    % template resolution 
      regres = reg(regstri).opt.rres; if isinf(regres), regres = tmpres; end    % registration resolution
      newres = job.extopts.vox(voxi); if isinf(newres), newres = tmpres; end    % output resolution


      % mat matrices for different spaces
      M0     = res.image.mat;                                                   % for (interpolated) individual volume
      M1d    = tmpM;                                                            % for template 
      imat=spm_imatrix(tmpM); imat(7:9)=imat(7:9) * newres/tmpres; M1=spm_matrix(imat); 
      imat=spm_imatrix(tmpM); imat(7:9)=imat(7:9) * regres/tmpres; M1r=spm_matrix(imat); 
      imat=spm_imatrix(tpmM); imat(7:9)=imat(7:9) * newres/tpmres; M1t=spm_matrix(imat); 
      imat=spm_imatrix(tpmM); imat(7:9)=imat(7:9) * regres/tpmres; M1tr=spm_matrix(imat); 
%M1t = M1; M1tr = M1r;  % tpm origin is diffent in SPM and review
      
     % if job.extopts.regstr(regstri)==0, M1r=M1; end % Dartel only!
    
      % image dimension 
      VT  = res.image(1);
      if isfield(res,'imagesc'); VT0 = res.imagec(1); else VT0 = res.image0(1); end
      idim = VT.dim(1:3);                                                     % (interpolated) input image resolution
      %idim = res.image(1).dim(1:3);                                          % input image resolution
      odim = floor(res.tpm2{1}(1).dim * tmpres/newres/2)*2+1;                     % output image size
      rdim = floor(res.tpm2{1}(1).dim * tmpres/regres/2)*2+1;                     % registration image size
      if job.extopts.regstr(regstri)==0, rdim = odim; end % Dartel only!


      % to write a correct x=-1 output image, we have to be sure that the x value of the bb is negative
      if res.bb(1)<res.bb(2), bbt=res.bb(1); res.bb(1)=res.bb(2); res.bb(2)=bbt; clear bbt; end

      % matrix to create the affine/rigide transformation matrices
      mm   = [[res.bb(1,1) res.bb(2,1) res.bb(1,1) res.bb(2,1) res.bb(1,1) res.bb(2,1) res.bb(1,1) res.bb(2,1); 
               res.bb(1,2) res.bb(1,2) res.bb(2,2) res.bb(2,2) res.bb(1,2) res.bb(1,2) res.bb(2,2) res.bb(2,2);
               res.bb(1,3) res.bb(1,3) res.bb(1,3) res.bb(1,3) res.bb(2,3) res.bb(2,3) res.bb(2,3) res.bb(2,3)]; 
               ones(1,8)];
      vx2  = M1\mm;
      vx2d = M1d\mm;
      vx3  = M1t\mm; %ones(4,8); vx3([5,13,21,29])=odim(1); vx3([10,14,26,30])=odim(2); vx3([19,23,27,31])=odim(3); % output image
      mat  = mm/vx3; 

      % individual parameters
      trans.native.Vo = VT0;
      trans.native.Vi = res.image(1);

      % affine parameters
      imat=spm_imatrix(vx2d/vx3); imat(7:9)=imat(7:9) * regres/tmpres; Mad=spm_matrix(imat); 
      % ExploreASL fix - rounding the results of inversion to increase reproducibility between different OS and Matlab versions
      % Ma            = M0\inv(res.Affine)*M1t*vx2/vx3;                            % individual to registration space
      Ma            = xASL_round(M0\inv(res.Affine)*M1t*vx2/vx3,12);                             % individual to registration space
      mat0a         = res.Affine\M1t*vx2/vx3;                                    % mat0 for affine output
      mata          = mm/vx3;                                                    % mat  for affine output
      trans.affine  = struct('odim',odim,'mat',mata,'mat0',mat0a,'M',Ma);        % structure for cat_io_writenii
      clear mata; 

      % rigid parameters
      [M3,R]        = spm_get_closest_affine( affind(rgrid( idim ) ,M0) , ...
                        affind(Yy,tpmM) , single(Ycls{1})/255); clear M3;        %#ok<ASGLU>
      % ExploreASL fix - rounding the results of inversion to increase reproducibility between different OS and Matlab versions                  
      % Mr            = M0\inv(R)*M1t*vx2/vx3;                                     % transformation from subject to registration space
      Mr            = xASL_round(M0\inv(R)*M1t*vx2/vx3,12);                                      % transformation from subject to registration space
      mat0reg       = R\M1tr*vx2/vx3;                                             
      mat0r         = R\M1t*vx2/vx3;                                             % mat0 for rigid ouput
      matr          = mm/vx3;                                                    % mat  for rigid ouput
      trans.rigid   = struct('odim',odim,'mat',matr,'mat0',mat0r,'M',Mr);        % structure for cat_io_writenii
      % ExploreASL fix - rounding the results of inversion to increase reproducibility between different OS and Matlab versions
      %res.rigid     = M0\inv(R); 
      res.rigid     = xASL_round(M0\inv(R),12); 
      clear matr; 

      % save old spm normalization used for atlas map
      trans.atlas.Yy = Yy; 

      % rigid vs. affine input in registration: 0 - rigid, 1 - affine (var only used by Dartel) 
      job.extopts.regra = 1; %job.extopts.regstr(regstri)==0; 
      if job.extopts.regra==0
        % rigid
        Mar = Mr; 
        TAR = R;
      else
        % affine
        Mar = Ma;
        TAR = res.Affine;
      end

      %n1    = 2;    % use GM/WM for dartel                
      if numel(Ycls)>n1, Ycls(n1+1:end) = []; end
      if exist('Ylesion','var') && sum(Ylesion(:))>0
        Yclso = Ycls; 
        %Ycls{1}(Ylesion) = 128; 
        %Ycls{2}(Ylesion) = 128; 
        %Ycls{1} = cat_vol_ctype( single(Ycls{1}) + 50.*Ylesion ); 
        %Ycls{2} = cat_vol_ctype( single(Ycls{1}) + 200.*Ylesion ); 
       % clear Ylesion;
      end

      


      if res.do_dartel && do_req
      %  -------------------------------------------------------------------
      %  do Dartel / Shooting registration
      %  -------------------------------------------------------------------
        if job.extopts.regstr(regstri)==0
        %  -----------------------------------------------------------------
        %  Dartel spatial normalization to given template
        %  -----------------------------------------------------------------

          stime = cat_io_cmd(sprintf('Dartel registration with %0.2f mm on a %0.2f mm Template',newres,tempres)); 
          fprintf('\n  Template: "%s"\n',job.extopts.templates{1})
      
          reg(regstri).opt.rres = newres; % 

          % dartel parameter 1
          rform = 0;    % regularization form: 0 - Linear Elastic Energy
          code  = 2;    % multinomial
          lmreg = 0.01; % LM regularization
          cyc   = 3;    % cycles
          its   = 3;    % relaxation iterations (inner iteration)
          n1    = 2;    % use GM/WM for dartel
          if fast, its = min(3,max(1,min(its,fast))); end % subiteration

          % rparam .. regularization parameters: mu, lambda, id
          % K      .. time steps
          param = struct('K',{0 0 1 2 4 6},'its',its, ...
            'rparam',{[4 2 1e-6],[2 1 1e-6],[1 0.5 1e-6],[0.5 0.25 1e-6],[0.25 0.125 1e-6],[0.25 0.125 1e-6]});

          % initialize varibles and load anatomical image in registration space 
          f = zeros([rdim(1:3) 2],'single');
          g = zeros([rdim(1:3) 2],'single');
          u = zeros([rdim(1:3) 3],'single');
          for k1=1:n1
            for i=1:rdim(3),
              f(:,:,i,k1) = single(spm_slice_vol(single(Ycls{k1}),Mar*spm_matrix([0 0 i]),rdim(1:2),[1,NaN])/255);
            end
          end
          
          if exist('Ylesion','var') && sum(Ylesion(:))>0
            Ylesion = single(Ylesion); 
            ls = zeros(rdim(1:3),'single');
            for i=1:rdim(3),
              ls(:,:,i) = single(spm_slice_vol(single(Ylesion),Mar*spm_matrix([0 0 i]),rdim(1:2),[1,NaN]));
            end
            ls(isnan(ls)) = 0;
          end
          
          %% iterative processing
          % ---------------------------------------------------------------------
          it0 = 0;  % main iteration number for output
          reg(regstri).dtc = zeros(1,6);
          for it = 1:numel(param)
            prm   = [rform, param(it).rparam, lmreg, cyc, its, param(it).K, code];
            % load new template for this iteration
            for k1=1:n1
              for i=1:rdim(3),
                g(:,:,i,k1) = single(spm_slice_vol(res.tpm2{it}(k1),Mad*spm_matrix([0 0 i]),rdim(1:2),[1,NaN]));
              end
              
              if exist('Ylesion','var') && sum(Ylesion(:))>0
                f(:,:,:,k1) = f(:,:,:,k1) .* (1-ls) + g(:,:,:,k1) .* ls;
              end
            end
            
            for j = 1:param(it).its,
              it0 = it0 + 1;
              [u,ll] = dartel3(u,f,g,prm);
              reg(regstri).lld(it0,:)  = ll ./ [prod(rdim) prod(rdim) newres^3]; 
              reg(regstri).lldf(it0,:) = [reg(regstri).lld(it0,1) / prod(regres),ll(1),ll(2),ll(1)+ll(2),ll(3)]; 
              cat_io_cprintf(sprintf('g%d',5+2*(mod(it0,its)==1)),...
                sprintf('% 5d | %6.4f | %8.0f %8.0f %8.0f %8.3f \n',it0,reg(regstri).lldf(it0,:)));
              
              if it0==1 % simplified! use values after first iteration rather than before
                [y0, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6); clear y0; %#ok<ASGLU>
                reg(regstri).ll(1,:)    =  reg(regstri).lldf(it0,:); 
                dtx = dt; 
                dtx(dtx>eps & dtx<1)    = 1./dtx(dtx>eps & dtx<1); 
                reg(regstri).dtc(1)     =  mean(abs(dtx(:)-1)); 
                reg(regstri).rmsdtc(1)  =  mean((dtx(:)-1).^2).^0.5;
                dtg = cat_vol_grad(single(dtx)); 
                reg(regstri).rmsgdt(1)  = mean((dtg(:)).^2).^0.5;
                clear dtg dt dtx; 
              end
            end 
           
            [y0, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6); clear y0;  %#ok<ASGLU>
            reg(regstri).ll(it+1,:)    =  reg(regstri).lldf(it0,:); 
            reg(regstri).dtc(it+1)     =  mean(abs(dt(:)-1)); 
            reg(regstri).rmsdtc(it+1)  =  mean((dt(:)-1).^2).^0.5;
            dtg = cat_vol_grad(single(dt)); 
            reg(regstri).rmsgdt(it+1)  = mean((dtg(:)).^2).^0.5;
            clear dtg; 
            % xASL hack - saving the intermediate steps of the DARTEL registration
			if job.extopts.xasl_savesteps
				xASL_wrp_DARTELSaveIntermedTrans(Yy,u,odim,rdim,idim,Mar,mat,M0,M1,VT0.fname,it);
			end
          end
          reg(regstri).rmsdt         = mean((dt(:)-1).^2).^0.5; 
          reg(regstri).dt            = mean(abs(dt(:)-1));
          if ~debug; clear f g; end


          %% jacobian 
          if job.output.jacobian.warped || debug || (isfield(job.extopts,'multigreg') && job.extopts.multigreg)
            trans.jc = struct('u',u,'odim',idim); 
          end

          % deformation
          y0      = spm_dartel_integrate(reshape(u,[rdim(1:3) 1 3]),[0 1], 6); clear u;
          prm     = [3 3 3 0 0 0];
          Coef    = cell(1,3);
          Coef{1} = spm_bsplinc(y0(:,:,:,1),prm);
          Coef{2} = spm_bsplinc(y0(:,:,:,2),prm);
          Coef{3} = spm_bsplinc(y0(:,:,:,3),prm);
          clear y0;
          [t1,t2] = ndgrid(1:idim(1),1:idim(2),1); t3 = 1:idim(3);

          Yyd = Yy; 
          for z=1:idim(3)
            [t11,t22,t33] = defs2(Coef,z,Mar,prm,t1,t2,t3);
            Yyd(:,:,z,1) = t11;
            Yyd(:,:,z,2) = t22;
            Yyd(:,:,z,3) = t33;
          end
          clear Coef t1 t2 t3 t11 t22 t33 z
          M = mat\M1;
          for i=1:size(Yyd,3),
            t1          = Yyd(:,:,i,1);
            t2          = Yyd(:,:,i,2);
            t3          = Yyd(:,:,i,3);
            Yyd(:,:,i,1) = M(1,1)*t1 + M(1,2)*t2 + M(1,3)*t3 + M(1,4);
            Yyd(:,:,i,2) = M(2,1)*t1 + M(2,2)*t2 + M(2,3)*t3 + M(2,4);
            Yyd(:,:,i,3) = M(3,1)*t1 + M(3,2)*t2 + M(3,3)*t3 + M(3,4);
          end
          clear t1 t2 t3 M; 

          vx_vols  = sqrt(sum(M0(1:3,1:3).^2));  
          vx_volt  = sqrt(sum(M1(1:3,1:3).^2));  
          interpol = any(vx_vols>vx_volt) + any(vx_vols/2>vx_volt);
          if interpol
            Yx = zeros(min([inf inf inf 3],size(Yyd)*2^interpol - 2^interpol + 1),'single');
            if interpol
              for i=1:3, Yx(:,:,:,i) = single(interp3(Yyd(:,:,:,i),interpol,'cubic')); end % linear
            end
          else
            Yx = Yyd;
          end

          trans.warped = struct('y',Yyd,'yx',Yx,'odim',odim,'M0',M0,'M1',M1,'M2',M1\TAR*M0,'dartel',1);
          clear Yyd; 


          if ~debug
            if job.extopts.verb<1, fprintf(sprintf('%s',repmat('\b',1,it0*47 + 2))); fprintf('\n'); end
            fprintf('%s %4.0fs\n',repmat(' ',1,66),etime(clock,stime)); 
          else
            cat_io_cmd(' ','',''); cat_io_cmd('','','',job.extopts.verb,stime); 
          end

        elseif job.extopts.regstr(regstri)>0
        %% -----------------------------------------------------------------
        %  Geodesic Shooting with Template 0 to 4 
        %  -----------------------------------------------------------------
          dreg.rigidShooting = 0; % rigid vs. affine Shooting 
 
          % multiresolution parameter
          % this part may require further work 
          tempres2 = reg(regstri).opt.resfac * regres;  % registration resolution
          if numel(reg(regstri).opt.ll3th)~=numel(tempres2), reg(regstri).opt.ll3th = repmat(reg(regstri).opt.ll3th(1),numel(tempres2)); end
          rdims   = zeros([numel(reg(regstri).opt.resfac),3]); 
          for ri=1:numel(reg(regstri).opt.resfac), rdims(ri,:) = floor(res.tpm2{1}(1).dim * tmpres/regres / reg(regstri).opt.resfac(ri)); end
          
          Mrregs  = cell(size(reg(regstri).opt.resfac)); 
          imat=spm_imatrix(tmpM); imat(7:9)=imat(7:9) * regres/tmpres; M1rr=spm_matrix(imat); 
          vx2rr = M1rr\mm; 
   
          for ri=1:numel(reg(regstri).opt.resfac)
            vx3rr = ones(4,8); vx3rr([5,13,21,29])=rdims(ri,1); vx3rr([10,14,26,30])=rdims(ri,2); vx3rr([19,23,27,31])=rdims(ri,3); % registration image
            vxtpm = tpmM\mm; % registration image
            Mads{ri} = (tmpM\mm)/vx3rr;
            if dreg.rigidShooting
              if ri==1, mat0reg = R\M1rr * vx2rr/vxtpm; end 
              % ExploreASL fix - rounding the results of inversion to increase reproducibility between different OS and Matlab versions
              %Mrregs{ri} = M0\inv(R)*M1rr*vx2rr/vx3rr; %;    
              Mrregs{ri} = xASL_round(M0\inv(R)*M1rr*vx2rr/vx3rr,12);
            else
              if ri==1, mat0reg = res.Affine\M1rr * vx2rr/vxtpm; end 
              % ExploreASL fix - rounding the results of inversion to increase reproducibility between different OS and Matlab versions
              %Mrregs{ri} = M0\inv(res.Affine)*M1rr*vx2rr/vx3rr;  
              Mrregs{ri} = xASL_round(M0\inv(res.Affine)*M1rr*vx2rr/vx3rr,12);
            end
            if ri==1, Mys{ri} = eye(4); else Mys{ri}= eye(4); Mys{ri}(1:12) = Mys{ri}(1:12) * reg(regstri).opt.resfac(ri)/reg(regstri).opt.resfac(ri-1); end;
          end
       
          % if saffine, Mrregs{1} = M0\inv(res.Affine)*M1rr*vx2rr/vx3rr; end       
          if job.extopts.verb
            if reg(regstri).opt.stepsize>10^-3 || regres~=tmpres
              if job.extopts.regstr(regstri)>0 && job.extopts.regstr(regstri)<=1
                stime   = cat_io_cmd(sprintf('Optimized Shooting registration with %0.2f:%0.2f:%0.2f mm (regstr=%0.2f)',...
                  tempres2(1),reg(regstri).opt.stepsize,tempres2(end),job.extopts.regstr(regstri))); 
              else
                stime   = cat_io_cmd(sprintf('Optimized Shooting registration with %0.2f:%0.2f:%0.2f mm',...
                  tempres2(1),reg(regstri).opt.stepsize,tempres2(end))); 
              end
            else
              stime   = cat_io_cmd(sprintf('Default Shooting registration with %0.2f mm',reg(regstri).opt.rres)); 
            end
            fprintf('\n  Template: "%s"\n',job.extopts.templates{1})
          end    

          % shooting parameter
          sd = spm_shoot_defaults;          % load shooting defaults
          if fast, reg(regstri).opt.nits = min(reg(regstri).opt.nits,5*fast); end  % at least 5 iterations to use each tempalte
          if  job.extopts.regstr(regstri)==6
            % need finer schedule for coarse to fine for ll3 adaptive threshold
            nits       = reg(regstri).opt.nits;;        % default was 24  
            lam        = 1/(nits*0.2);                  % Decay of coarse to fine schedule (default 0.5)
            inter      = 24;                            % Scaling of parameters at first iteration
            sd.sched   = (inter-1) * exp(-lam * ((1:(nits+1))-1)) + 1;
% die variante bis ...4 ist sch?n weich und kommt an die 0.07 ran ... 
% das reicht aber wohl nicht aus um 0.06 zu erreichen!
            sd.sched   = ((sd.sched - min(sd.sched)) / (max(sd.sched) - min(sd.sched)) * (inter - 1) + 4) ;
            maxoil     = 8;                            % Maximum number of time steps for integration
            sd.eul_its = round((0:(nits-1))*(maxoil-0.5001)/(nits-1)+1); % Start with fewer steps
            nits       = numel(sd.sched)-1;            % Shooting iterations 
            tmpl_no    = floor(((1:nits)-1)/(nits-1)*(numel(res.tpm2)-0.51))+1; % Sort out which template for each iteration (default = round with more hr-iter)
            %tmpl_no    = floor( (((1:nits)-1)/(nits-1)).^(1/2) * (numel(res.tpm2)-0.51)) + 1; 
% lief mit 512 ganz gut - ca. 1/5 1.5 mmm Operationen ... vieleicht shift bissel vergr??ern          
            %sd.rparam  = [1e-2 0.01 0.2 0.05 0.2];  % Regularisation parameters for deformation
            % abs. disp., laplacian, bending engery, linear elasticity mu, le lambda [1e-4 0.001 0.2 0.05 0.2]
            sd.rparam  = [1e-2 0.01 0.2 0.05 0.2];
            
          elseif (job.extopts.regstr(regstri)>0 || regres~=tmpres) && job.extopts.regstr(regstri)~=4 
            % need finer schedule for coarse to fine for ll3 adaptive threshold
            nits       = reg(regstri).opt.nits;        % default was 24  
            lam        = 1/(nits*0.2); %1/4;                          % Decay of coarse to fine schedule (default 0.5)
            inter      = 32;                           % Scaling of parameters at first iteration
            sd.sched   = (inter-1)*exp(-lam*((1:(nits+1))-1))+1; %sd.sched = sd.sched/sd.sched(end);
            sd.sched   = ((sd.sched - min(sd.sched)) / (max(sd.sched) - min(sd.sched)) * (inter-1) + 1) ;
            maxoil     = 8;                            % Maximum number of time steps for integration
            sd.eul_its = round((0:(nits-1))*(maxoil-0.5001)/(nits-1)+1); % Start with fewer steps
            nits       = numel(sd.sched)-1;            % Shooting iterations 
            tmpl_no    = floor(((1:nits)-1)/(nits-1)*(numel(res.tpm2)-0.51))+1; % Sort out which template for each iteration (default = round with more hr-iter)
            %tmpl_no    = floor( (((1:nits)-1)/(nits-1)).^(nits.^(1/10)) * (numel(res.tpm2)-0.5)) + 1; 
          else
            nits       = reg(regstri).opt.nits;         
            tmpl_no    = round(((1:nits)-1)/(nits-1)*(numel(res.tpm2)-0.51))+1; 
            lam        = nan; 
            inter      = nan;
          end
          reg(regstri).shooting = struct('nits',nits,'lam',lam,'inter',inter,'sched',sd.sched,'eul_its',sd.eul_its,'tmpl_no',tmpl_no');

          
          %% The actual work
          % ---------------------------------------------------------------------
 %R=res.Affine;
          spm_progress_bar('Init',nits,'Shooting progress','iterations')
    
          prm2 = 1; it = 1; itd = 1; reg(regstri).dtc = zeros(1,5); ll  = zeros(1,2);
          while it<=nits || (prm2>1.1 && (it<=nits*1.1 || prm2>(3+it/nits) ) )
  
            itime = clock;
            if it==nits+1, cat_io_cprintf([.2 .2 1],'Stoppage time!\n'); end
            if it==1 || (tmpl_no(min(nits,it))~=tmpl_no(min(nits,it)-1)) 
              ti  = tmpl_no(min(nits,it)); %ittime(it) = clock;

              if debug && it>1, fo=f{1}; end %#ok<NASGU> % just for debugging    

 %itu = 0;
 %itx = 0; 
              %% load rigide/affine data
              % this has to be done by linear interpolation
              % 
              f = cell(1,n1+1); f{n1+1} = ones(rdims(ti,1:3),'single'); 
              for k1=1:n1
                Yclsk1 = single(Ycls{k1}); 
                f{k1}  = zeros(rdims(ti,1:3),'single');
                fnl    = zeros(rdims(ti,1:3),'single');   
                if reg(regstri).opt.resfac(ti)>1, spm_smooth(Yclsk1,Yclsk1,repmat((reg(regstri).opt.resfac(ti)-1) * 2,1,3)); end
                for i=1:rdims(ti,3),
                  f{k1}(:,:,i) = single(spm_slice_vol(Yclsk1,Mrregs{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN])/255);
                  fnl(:,:,i)   = single(spm_slice_vol(Yclsk1,Mrregs{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[5,NaN])/255);
                end
                f           = combine_linear_spline(f,fnl,Mrregs{ti},k1); clear fnl;
                msk         = ~isfinite(f{k1});
                f{k1}(msk)  = 0;
                clear msk; 
              end
              if debug, fx = f{1}; end %#ok<NASGU> % just for debugging

              %% template
              g = cell(1,n1+1); g{n1+1} = ones(rdims(ti,1:3),'single');
              for k1=1:n1
                g{k1} = zeros(rdims(ti,1:3),'single');
                tpm2k1 = res.tpm2{ti}(k1).private.dat(:,:,:,k1); 
                if reg(regstri).opt.resfac(ti)>1, spm_smooth(tpm2k1,tpm2k1,repmat((reg(regstri).opt.resfac(ti)-1) * 2,1,3)); end
                for i=1:rdims(ti,3),
                  g{k1}(:,:,i) = single(spm_slice_vol(tpm2k1,Mads{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN]));
                end
                g{k1}(isnan(g{k1}(:))) = min(g{k1}(:)); % remove boundary interpolation artefact
                g{n1+1} = g{n1+1} - g{k1};
                if debug && k1==1, gx = g{1}; end %#ok<NASGU> % just for debugging
                g{k1} = spm_bsplinc(log(g{k1}), sd.bs_args);
              end
              g{n1+1} = log(max(g{n1+1},eps)); 
%%
              if exist('Ylesion','var') && sum(Ylesion(:))>0
                Yclsk1 = single(Ylesion); 
                if reg(regstri).opt.resfac(ti)>1, spm_smooth(Yclsk1,Yclsk1,repmat((reg(regstri).opt.resfac(ti)-1) * 2,1,3)); end
                ls = zeros(rdims(ti,1:3),'single'); lsnl = ls; 
                for i=1:rdims(ti,3),
                  ls(:,:,i)   = single(spm_slice_vol(Yclsk1,Mrregs{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN])); 
                  lsnl(:,:,i) = single(spm_slice_vol(Yclsk1,Mrregs{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[5,NaN])); 
                end
                ls(isnan(ls)) = 0; % ExploreASL fix - to avoid having NaNs
                lsnl(isnan(lsnl)) = 0; % ExploreASL fix - to avoid having NaNs
                ls = combine_linear_spline(ls,lsnl,Mrregs{ti}); clear lsnl; 
                for k1=1:numel(g)-1
                  tpm2k1 = res.tpm2{ti}(k1).private.dat(:,:,:,k1); 
                  t = zeros(rdims(ti,1:3),'single'); tnl = t; 
                  for i=1:rdims(ti,3),
                    t(:,:,i)   = single(spm_slice_vol(tpm2k1,Mads{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN]));
                    tnl(:,:,i) = single(spm_slice_vol(tpm2k1,Mads{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[5,NaN]));
                  end
                  t = combine_linear_spline(t,tnl,Mads{ti}); clear tnl; 
                  
                  f{k1} = f{k1} .* (1-ls) + t .* ls;
                  
                  msk          = ~isfinite(f{k1});
                  f{k1}(msk)   = 0;
                  f{n1+1}      = f{n1+1} - f{k1}; 
                  f{n1+1}(msk) = 0.00001;
                end
                
                %Ycls{1} = cat_vol_ctype( single(Ycls{1}) + 50.*Ylesion ); 
                %Ycls{2} = cat_vol_ctype( single(Ycls{1}) + 200.*Ylesion ); 
               % clear Ylesion;
              else
                for k1=1:numel(f)-1
                   f{end}     = f{end} - f{k1}; 
                end
                msk          = ~isfinite(f{k1});
                f{n1+1}(msk) = 0.00001;
              end
              
              %% loading segmentation and creating of images vs. updating these maps
              ll  = zeros(1,3);
              if it==1
                % create shooting maps
                %y   = affind( squeeze( reshape( affind( spm_diffeo('Exp',zeros([rdims(ti,:),3],'single'),[0 1]), ...
                %      mat0reg), [rdims(ti,:),1,3] ) ) , inv(mat0reg)); clear def;                          % deformation field
                y   = affind( squeeze( reshape( affind( spm_diffeo('Exp',zeros([rdims(ti,:),3],'single'),[0 1]), ...
                      mat0reg), [rdims(ti,:),1,3] ) ) , xASL_round(inv(mat0reg),12)); clear def;     % deformation field                 
                u   = zeros([rdims(ti,:) 3],'single');                                                     % flow field
                dt  = ones(rdims(ti,:),'single');                                                          % jacobian
              elseif any(rdims(ti,:)~=rdims(ti-1,:))
                % updates only for changed resolutions

                % update resolution of shooting maps
    %             if debug
    %               for i=1:rdims(ti,3),
    %                  fox(:,:,i) = single(spm_slice_vol(fo,Mys{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN])) / Mys{ti}(1); % adapt for res
    %               end
    %             end

                %% size update y - deformation field
                yo = y;  
                y  = zeros([rdims(ti,:) 3],'single');            
                for k1=1:3
                  ynl = zeros(rdims(ti,:),'single');        
                  for i=1:rdims(ti,3),
                    y(:,:,i,k1) = single( spm_slice_vol(yo(:,:,:,k1),Mys{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN]) ) / Mys{ti}(1); % adapt for res
                    ynl(:,:,i)  = single( spm_slice_vol(yo(:,:,:,k1),Mys{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[5,NaN]) ) / Mys{ti}(1); % adapt for res
                  end
                  y = combine_linear_spline(y,ynl,Mys{ti},k1); clear ynl; 
                end
                y(~isfinite(y))=y(find(~isfinite(y))+1); % use neighbor value in case of nan
                if ~debug, clear yo; end

                %% size update u - flow field
                uo = u; 
                u  = zeros([rdims(ti,:) 3],'single');
                for k1=1:3
                  unl = zeros(rdims(ti,:),'single');        
                  for i=1:rdims(ti,3),
                    u(:,:,i,k1) = single(spm_slice_vol(uo(:,:,:,k1),Mys{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN])) / Mys{ti}(1); % (tempres(ti) / tempres(ti-1))^2; % adapt for res 
                    unl(:,:,i)  = single(spm_slice_vol(uo(:,:,:,k1),Mys{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[5,NaN])) / Mys{ti}(1); 
                  end
                  u = combine_linear_spline(u,unl,Mys{ti},k1); clear unl; 
                end
                u(~isfinite(u)) = eps;
                if ~debug, clear uo; end

                % size update dt
                dto = dt;
                dt  = zeros(rdims(ti,:),'single');
                for i=1:rdims(ti,3),
                  dt(:,:,i) = single(spm_slice_vol(dto,Mys{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN]));
                  dtnl(:,:,i) = single(spm_slice_vol(dto,Mys{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[5,NaN]));
                end
                dt = combine_linear_spline(dt,dtnl,Mys{ti}); clear dtnl; 
                dt(~isfinite(dt)) = 1;
                if ~debug, clear dto; end

                % [ux,ll(1),ll(2),ll(3)] = spm_shoot_update(g,f,u,y,dt,prm,sd.bs_args,sd.scale);
                %if ti>4
                %  ito= it; it = inf; continue;
                %end
              end
            end





            % More regularisation in the early iterations, as well as a less accurate approximation in the integration.
            % No, similar regularisation works in our case better and avoid to trap into local maxima.  
           % vxreg    = repmat(reg(regstri).opt.vxreg,1,3);  %
            vxreg    = repmat(reg(regstri).opt.vxreg .* reg(regstri).opt.resfac(tmpl_no(min(nits,it))),1,3);
            %prm      = [vxreg, sd.rparam * sd.sched(it+1)  * prm2  * prod(vxreg)]; 
%            prm      = [vxreg, min(sd.rparam * sd.sched(1)    *    1 / max(1,prod(vxreg(1:2)/4)) , ...
%                                   sd.rparam * sd.sched(it+1) * prm2 / max(1,prod(vxreg(1:2)/4)) ) ]; 
            prm      = [vxreg, min(sd.rparam * sd.sched(1)              *    1 / vxreg(1)/1.5 , ... / max(1,prod(vxreg/1.5^3))
                                   sd.rparam * sd.sched(min(nits,it)+1) * prm2 / vxreg(1)/1.5 ) ]; % / max(1,prod(vxreg/1.5^3))
            %prm(4:8) =  prm(4:8) ./ [4^tmpl_no(it) 2^(tmpl_no(it)-1) 1 1 1];;
            prm(4:8) =  prm(4:8) ./ [100*min(nits,it)/nits 1 1 1 1];;
            int_args = [sd.eul_its(min(nits,it)), sd.cyc_its]; 
            prm2 = prm2*0.9 + 0.1*1; % average with one to go back to the oringal value

            
            % Gauss-Newton iteration to re-estimate deformations for this subject
            llo=ll; 
            if job.extopts.verb
              if it<=nits
                pcol = sprintf('g%d',5+2*(it==1 || (tmpl_no(min(nits,it))~=tmpl_no(min(nits,it)-1))));
              else
                pcol = [0.3 0.3 0.7];
              end
              if reg(regstri).opt.stepsize<=10^-3
                cat_io_cprintf(pcol,sprintf('% 5d |',it));
              else
                cat_io_cprintf(pcol,sprintf('% 5d | %0.2f |',it,tempres2(ti)));
              end
            end

            if 1 %ti<0
              if job.extopts.regstr(regstri)==4 
                [txt,u,ll(1),ll(2)] = evalc('spm_shoot_update(g,f,u,y,dt,prm,sd.bs_args,sd.scale)');  %#ok<ASGLU>
                [y,J] = spm_shoot3d(u,prm,int_args); 
                dt = spm_diffeo('det',J); clear J
              else
                for i = 1%:max(1,round( (6 - tmpl_no(it))^2 / (1 * 5) ))
                  [txt,u,ll(1),ll(2)] = evalc('spm_shoot_update(g,f,u,y,dt,prm,sd.bs_args,sd.scale)');  %#ok<ASGLU>
                  [y,J] = spm_shoot3d(u,prm,int_args); 
                  dt = spm_diffeo('det',J); clear J                 
                end
              end
              if job.extopts.verb
                cat_io_cprintf(pcol,sprintf('%7.4f%8.4f%8.4f', ...
                  ll(1)/numel(u), ll(2)/numel(u), (ll(1)+ll(2))/numel(u))); % * prod(vxreg)
                if isfield(job,'export') && job.export
                  cat_io_cprintf(pcol,...
                    sprintf(' |%6.2f%7.2f | %7.4f%8.4f%8.4f%8.4f%8.4f%8.4f', prm2, sd.sched(min(nits,it)) * prm2 , prm(4:end) ));
                end
              end
               
              %%
              wc  = 0; maxsubiter = 5;
              dtm = ( 10 + max(0,min(20,20 * it./nits)) ) * tempres2(ti); 
              while job.extopts.regstr(regstri)~=4 && wc<=maxsubiter && ... %(it + min(nits/10,20))<nits && ...
                  ( any(~isfinite(dt(:))) || any(dt(:)>dtm) || any(dt(:)<1/dtm) )
                fprintf('.'); 
                wc  = wc + 1; 
               
                if wc==1 && debug
                  uo=u; yo=y; dto=dt; 
                elseif wc>maxsubiter && debug
                  fprintf('>'); u=uo; y=yo; dt=dto;  %#ok<NASGU>
                end
                
                
                %% local Laplace filter for continuous field smoothness
% war ok ... nochmal eine stufe schw?cher?   
if job.extopts.regstr(regstri)==6
                mxu = max(1,min(dtm,2.^( (3 - wc) ))); %mxu = (1.3).^(11 - wc); % dieser Wert war zu stark bei hohen Iterationen!
                for i=1:3  
                  if 1%tmpl_no(min(it,nits))<tmpl_no(end)
                    u(:,:,:,i) = cat_vol_smooth3X(u(:,:,:,i),2 - ti./nits);% 0.5 + ti./nits);
                  end
                 % u(:,:,:,i) = ( cat_vol_laplace3R(u(:,:,:,i) / mxu - min(min(min(u(:,:,:,i)))) / mxu,...
                 %   true(size(dt)) , 0.1) + min(min(min(u(:,:,:,i)))) / mxu ) * mxu;
                end
end

                %% update parameter
                prm2 = min(48,max(1 , prm2 * (4 - 2*( ti./nits) ) ) ); 
                prm(4:end) = min(sd.rparam * sd.sched(1),prm(4:end) * prm2); 
               
                [y,J] = spm_shoot3d(u,prm,int_args);  %#ok<ASGLU>
                dt = spm_diffeo('det',J); clear J     %#ok<NASGU>
                
                % reprocess this iteration
                [txt,u,ll(1),ll(2)] = evalc('spm_shoot_update(g,f,u,y,dt,prm,sd.bs_args,sd.scale)');  %#ok<ASGLU>
                [y,J] = spm_shoot3d(u,prm,int_args); 
                dt = spm_diffeo('det',J); clear J txt

                %
                if wc==0
                  if ~exist('lllo','var')
                    lll = ll; 
                  else
                    if  (lll(1) - ll(1))/numel(u)<0.01 && ...
                        (lll(2) - ll(2))/numel(u)<0.01 
                      it = max(it+1,find([tmpl_no,nits]>tmpl_no(it),1,'first')); 
                      reg(regstri).ll(ti,1:4) = [ll(1)/numel(dt) ll(2)/numel(dt) (ll(1)+ll(2))/numel(dt) ll(2)]; 
                      reg(regstri).dtc(ti) = mean(abs(dt(:)-1)); 
                      fprintf('\b>'); 
                    end
                  end
                end
                
                % repeat some interations
                % it = floor( ( it + find( tmpl_no == tmpl_no(it) , 1 ,'first'))/2); 
              end
              
              
              % accelelerate processing     
              if job.extopts.regstr(regstri)~=4  && ...
                 (llo(1) - ll(1))/numel(u)<0.005 && ...
                 (llo(2) - ll(2))/numel(u)<0.005 
                 %tmpl_no(it) < max(tmpl_no)
                fprintf('+');
                if prm2<1
                  prm2 = prm2 / 1.01; 
                else
                  prm2 = prm2 / 1.25;
                end
              elseif wc==0 && ...
                ~(it==1 || (tmpl_no(min(nits,it))~=tmpl_no(min(nits,it)-1))) && ...
                ( any(dt(:)>dtm*0.66) || any(dt(:)<1/(dtm*0.66)) ) 
                  fprintf('-');
                if 0 && (any(dt(:)>dtm*0.9) || any(dt(:)<1/(dtm*0.9)) ) 
                  for i=1:3  
                    spm_smooth(u(:,:,:,i),u(:,:,:,i),2*[1 1 1]); %./vxreg)
                  end
                  prm2 = min(120,max(1 , prm2 * 1.5 ) ); 
                else                  
                  prm2 = min(120,max(1 , prm2 * 1.25 ) ); 
                end
              end
              
              
              %{
              wc = 0; dtm = 200; 
              while job.extopts.regstr(regstri)~=4 && any(~isfinite(dt(:)) | dt(:)>dtm | dt(:)<1/dtm) && wc<11
                fprintf('!'); wc=wc+1;
              
                if wc>4
                  fprintf('x');
                end

                usemask = 0; 
                if usemask
                  Yg   = cat_vol_grad(dt); 
                  Ygu = zeros(size(Yg),'single'); 
                  for i=1:3, Ygu = max(Ygu,cat_vol_grad(dt)); end
                  dtm2 = dtm*2/3;
                  mask = cat_vol_smooth3X(~isfinite(dt) | dt>dtm*dtm2 | dt<1/dtm2 | Yg>5 | Ygu>5 , 4 / reg(regstri).opt.vxreg ); % reg(regstri).opt.vxreg^0.5); 
                  mask = mask ./ max(mask(:));
                end
                
                for i=1:3, 
                  mxu = 2;
                  u(:,:,:,i) = cat_vol_laplace3R(u(:,:,:,i) / mxu,abs(u(:,:,:,1))>0.05,0.1) * mxu;
                  if tmpl_no(min(nits,it))<=max(tmpl_no)
                    if usemask
                      us = cat_vol_smooth3X(u(:,:,:,i), 4 ); % / reg(regstri).opt.vxreg ); % ); 
                      u(:,:,:,i) = us.*(mask) + (1-mask).*u(:,:,:,i); 
                    else
                      u(:,:,:,i) = cat_vol_smooth3X(u(:,:,:,i), 4 ); % /reg(regstri).opt.vxreg 
                    end
                  else
                   % ti = inf;
                  end
                end
                if tmpl_no(min(nits,it))<=max(tmpl_no)
                  [txt,u,ll(1),ll(2),ll(3)] = evalc('spm_shoot_update(g,f,u,y,dt,prm,sd.bs_args,sd.scale)');  %#ok<ASGLU>
                  [y,J] = spm_shoot3d(u,prm,int_args); 
                  dt = spm_diffeo('det',J); clear J

                  prm2 = prm2 * 1.2; prm(4:end) = prm(4:end) * prm2; 
                  itu = itu + 0.5;  
                  if itu < 6 - tmpl_no(min(nits,it))
                    %it = find( tmpl_no == tmpl_no(it) , 1 ,'first'); 
                  else
                    if itx<40
                      it = it - 1; itx = itx + 1; 
                    end
                  end
                end
              end
              %}
              
              cat_io_cprintf(sprintf('g%d',5+2*(it==1 || (tmpl_no(min(nits,it))~=tmpl_no(min(nits,it)-1)))),'\n')
              
            else
              % debugging
              fprintf('%3.0f | %s \n',etime(clock,itime),sprintf('%0.8f ',cat_stat_nanmean(g{1}(:)),cat_stat_nanmean(f{1}(:)),...
                cat_stat_nanmean(u(:)),cat_stat_nanmean(y(:)),cat_stat_nanmean(dt(:)))); %#ok<UNRCH>
              ll = zeros(1,3);
            end

            % save iteration parameter for later analysis
            if it==1 || (tmpl_no(min(nits,it))~=tmpl_no(min(nits,it)-1)) 
              reg(regstri).ll(ti,1:4)  = [ll(1)/numel(dt) ll(2)/numel(dt) (ll(1)+ll(2))/numel(dt) ll(2)]; 
              dtx = dt; 
              dtx(dtx>eps & dtx<1)     = 1./dtx(dtx>eps & dtx<1); 
              reg(regstri).dtc(ti)     = mean(abs(dtx(:)-1)); 
              reg(regstri).rmsdtc(ti)  = mean((dtx(:)-1).^2).^0.5;
              dtg = cat_vol_grad(single(dtx)); 
              reg(regstri).rmsgdt(ti)  = mean((dtg(:)).^2).^0.5;
              clear dtx;
              clear dtg; 
            end

            % default Shooting error detection
            if any(~isfinite(dt(:)) | dt(:)>10000 | dt(:)<1/10000)
              cat_io_cprintf('err',sprintf('Problem with Shooting (dets: %g .. %g)\n', min(dt(:)), max(dt(:)))); it=inf;
            end

            % xASL hack - saving the intermediate steps of the GS registration
			if rdim(1) > 80
				genIt = find(ceil(nits*[1 2 3 4 5 6]/6) == it);
				if ~isempty(genIt)
					xASL_wrp_GSSaveIntermedTrans(y,idim,odim,rdim,M0,M1,R,M1t,M1r,VT0.fname,genIt(1));
				end
			end
			
            % avoid unneccessary iteration
            if job.extopts.regstr(regstri)>0 && job.extopts.regstr(regstri)~=4 && job.extopts.regstr(regstri)~=6 && ...
               ll(1)<0.1 && ...
                (abs(llo(1)-ll(1))/numel(u)<0.0001 && abs(llo(2)-ll(2))/numel(u)<0.0001 && ...
                ( ti>1 || (ti==1 && ll(1)/numel(u)<1 && ll(1)/max(eps,llo(1))<1 && ll(1)/max(eps,llo(1))>(1-0.01) )) && ...
                ( ll(1)/numel(u)<1 && ll(1)/max(eps,llo(1))<1 && ll(1)/max(eps,llo(1))>(1-reg(regstri).opt.ll1th) ))
              it = max(it+1,find([tmpl_no,nits]>tmpl_no(it),1,'first')); 
              reg(regstri).ll(ti,1:4) = [ll(1)/numel(dt) ll(2)/numel(dt) (ll(1)+ll(2))/numel(dt) ll(2)]; 
              reg(regstri).dtc(ti) = mean(abs(dt(:)-1)); 

            else
              it = it+1; itd = itd+1; 
            end
            spm_progress_bar('Set',it);
          end
          spm_progress_bar('Clear');

          % some parameter for later ..
          dtx = dt; 
          dtx(dtx>eps & dtx<1)       = 1./dtx(dtx>eps & dtx<1); 
          reg(regstri).rmsdt         = mean((dtx(:)-1).^2).^0.5; 
          reg(regstri).dt            = mean(abs(dtx(:)-1));
          reg(regstri).dtc(ti+1)     = mean(abs(dtx(:)-1)); 
          reg(regstri).rmsdtc(ti+1)  = mean((dtx(:)-1).^2).^0.5; 
          dtg = cat_vol_grad(single(dtx)); 
          reg(regstri).rmsgdt(ti+1)  = mean((dtg(:)).^2).^0.5;
          clear dtg; 
          reg(regstri).ll(ti+1,1:4)  = [ll(1)/numel(dt) ll(2)/numel(dt) (ll(1)+ll(2))/numel(dt) ll(2)]; 
          clear dt1; 




          %% preparte output
          if job.extopts.verb
            if job.extopts.regstr(regstri)==0
              cat_io_cmd(sprintf('Dartel registration with %0.2f mm takes',tempres(1)));
            elseif reg(regstri).opt.stepsize>10^-3  
              cat_io_cmd(sprintf('Shooting registration with %0.2f:%0.2f:%0.2f mm takes',tempres2(1),reg(regstri).opt.stepsize,tempres2(end))); 
            else
              cat_io_cmd(sprintf('Shooting registration with %0.2f mm takes',tempres2(1))); 
            end
            itime = cat_io_cmd(sprintf('  Prepare output'),'','',job.extopts.verb,stime);
          end
          
          
          %% update for output resolution 
          if any(odim ~= rdim)
            eyev = eye(4); eyev(1:end-1) = eyev(1:end-1) * M1t(1)./M1r(1);
            yid  = zeros([odim 3],'single'); 
            %[D,I] = cat_vbdist(single(~(isnan(y)))); y = y(I); 
            for k1=1:3
              yidnl = zeros(odim,'single');        
              for i=1:odim(3),
                yid(:,:,i,k1) = single(spm_slice_vol(y(:,:,:,k1),eyev*spm_matrix([0 0 i]),odim(1:2),[1,NaN])) / eyev(1); 
                yidnl(:,:,i)  = single(spm_slice_vol(y(:,:,:,k1),eyev*spm_matrix([0 0 i]),odim(1:2),[5,NaN])) / eyev(1); % adapt for res
              end
              yid = combine_linear_spline(yid,yidnl,eyev,k1);   
            end
          else
            yid = y; 
          end
         
          if dreg.rigidShooting
              % ExploreASL fix - rounding the results of inversion to increase reproducibility between different OS and Matlab versions
			  %yi  = spm_diffeo('invdef',yid,idim,inv(M1t\R*M0),eye(4));  % output yi in anatomical resolution
			  yi  = spm_diffeo('invdef',yid,idim,xASL_round(inv(M1t\R*M0),12),eye(4));  % output yi in anatomical resolution
		  else
			  % ExploreASL fix - rounding the results of inversion to increase reproducibility between different OS and Matlab versions
			  %yi  = spm_diffeo('invdef',yid,idim,inv(M1t\res.Affine*M0),eye(4));           % output yi in anatomical resolution
			  yi  = spm_diffeo('invdef',yid,idim,xASL_round(inv(M1t\res.Affine*M0),12),eye(4));           % output yi in anatomical resolution
          end
          dt2 = spm_diffeo('def2det',yid); if ~debug, clear yid; end  
          

          % interpolation for improved output ... need update 2012/12
          if dreg.fast
            yx=yi; 
          else
            if 1
              vx_vols  = sqrt(sum(M0(1:3,1:3).^2));  
              vx_volt  = sqrt(sum(M1(1:3,1:3).^2));  
              interpol = any(vx_vols>vx_volt) + any(vx_vols/2>vx_volt);
              if interpol
                yx = zeros(min([inf inf inf 3],size(yi)*2^interpol - 2^interpol + 1),'single');
                if interpol
                  for i=1:3, yx(:,:,:,i) = single(interp3(yi(:,:,:,i),interpol,'cubic')); end % linear
                end
                %yx = yx; % * regres/newres;
              else
                yx = yi; % * regres/newres;
              end
            else
              % size update y - deformation field
              yx = zeros([odims 3],'single'); 
              My = eye(4); My(1:12) = My(1:12) * regres/newres;  Mys{ri}    = eye(4);
              Mys{ri}(1:12)    = Mys{ri}(1:12) * reg(regstri).opt.resfac(ri)/reg(regstri).opt.resfac(ri-1);
              for k1=1:3
                yxnl = zeros(odim,'single');
                for i=1:odims(ti,3),
                  yx(:,:,i,k1) = single(spm_slice_vol(yo(:,:,:,k1),My*spm_matrix([0 0 i]),odims(ti,1:2),[1,NaN])) * regres/newres; % adapt for res
                  yxnl(:,:,i)  = single(spm_slice_vol(yo(:,:,:,k1),My*spm_matrix([0 0 i]),odims(ti,1:2),[5,NaN])) * regres/newres; % adapt for res
                end
                yx = combine_linear_spline(yx,yxnl,My,k1); clear yxnl; 
                yx(~isfinite(yx) | yx<0) = eps; drawnow
              end
            end
          end
          if dreg.rigidShooting
            trans.warped = struct('y',yi,'yx',yx,'odim',odim,'M0',M0,'M1',M1,'M2',M1\R*M0,'dartel',2);
          else
            trans.warped = struct('y',yi,'yx',yx,'odim',odim,'M0',M0,'M1',M1,'M2',M1\res.Affine*M0,'dartel',2);
          end
          if job.output.jacobian.warped, trans.jc = struct('u',u,'odim',odim,'dt2',dt2); end % u ist nicht auf vox angepasst!

          if job.extopts.verb
            cat_io_cmd('','','',job.extopts.verb,itime); 
            cat_io_cmd(' ','',''); cat_io_cmd('','','',job.extopts.verb,stime); 
          end
        end


        % Report
        if (job.extopts.verb>1  && ~dreg.fast) || export

          %% preparte output directory
          if job.extopts.subfolders, mrifolder = 'mri'; else mrifolder = ''; end
          [pth,nam] = spm_fileparts(VT0.fname); 
          [temppp,tempff] = spm_fileparts(job.extopts.templates{1});  %#ok<ASGLU>
          if job.extopts.regstr(regstri)==0
            testfolder = sprintf('Dartel_%s_rr%0.1f_default',tempff,newres);
          elseif job.extopts.regstr(regstri)==4
            testfolder = sprintf('Shooting_%s_rr%0.1f_or%0.1f_default',tempff,regres,newres);
          else 
            if reg(regstri).opt.stepsize>10^-3 
              testfolder = sprintf('Shooting_%s_tr%0.1f_rr%0.1f-%0.1f_or%0.1f_regstr%0.1f_it%0.0f',tempff,tempres,tempres2([1,5]),newres,job.extopts.regstr(regstri),reg(regstri).opt.nits);
            else
              testfolder = sprintf('Shooting_%s_tr%0.1f_rr%0.1f_or%0.1f_regstr%0.1f_it%0.0f',tempff,tempres,regres,newres,job.extopts.regstr(regstri),reg(regstri).opt.nits);
            end
          end
          
          % display registration power profil 
          if job.extopts.expertgui==2
            mdisplay = [2 2 2 2 2]; 
          elseif job.extopts.expertgui==1
            mdisplay = [0 2 0 2 0]; 
          else
            mdisplay = [0 1 0 1 0]; 
          end  
          if job.extopts.regstr(regstri)==0, dartelfac = 1.5; else dartelfac = 1.0; end
          if mdisplay(1)
            fprintf('Registration power: \n'); 
            fprintf('%30s','Jacobian determinant: '); 
            QMC   = cat_io_colormaps('marks+',17);
            reg(regstri).reldtc = reg(regstri).dtc / max(reg(regstri).dtc); 
            color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
            if mdisplay(1)>1
              for dti=1:numel(reg(regstri).dtc)-1
                cat_io_cprintf( color(QMC,( 1 - reg(regstri).reldtc(dti)/dartelfac/0.25 ) *6),sprintf('%0.3f ',reg(regstri).reldtc(dti))); 
              end
              fprintf('| ');
            end
            cat_io_cprintf( color(QMC,(reg(regstri).dt - 0.05)/dartelfac/0.25 * 6), sprintf(' %0.6f ',reg(regstri).dt));
            fprintf('\n'); 
          end
          
          if mdisplay(2)
            fprintf('%30s','Jacobian determinant (RMS): '); 
            QMC   = cat_io_colormaps('marks+',17);
            reg(regstri).relrmsdtc = reg(regstri).rmsdtc; %/max(reg(regstri).rmsdtc); 
            color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
            if mdisplay(2)>1
              for dti=1:numel(reg(regstri).relrmsdtc)-1
                cat_io_cprintf( color(QMC,reg(regstri).relrmsdtc(dti)/dartelfac/0.5*6),sprintf('%0.3f ',reg(regstri).relrmsdtc(dti))); 
              end
              fprintf('| ');
            end
            cat_io_cprintf( color(QMC,(reg(regstri).rmsdt)/dartelfac/0.5 * 6), sprintf(' %0.6f ',reg(regstri).rmsdt));
            fprintf('\n'); 
          end
          
          if mdisplay(3)
            fprintf('%30s','Jacobian determinant'' (RMS): '); 
            QMC   = cat_io_colormaps('marks+',17);
            color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
            if mdisplay(3)>1
              for dti=1:numel(reg(regstri).rmsgdt)-1
                cat_io_cprintf( color(QMC,reg(regstri).rmsgdt(dti)/dartelfac/0.5*6),sprintf('%0.3f ',reg(regstri).rmsgdt(dti))); 
              end
              fprintf('| ');
            end
            cat_io_cprintf( color(QMC,(reg(regstri).rmsgdt(end))/dartelfac/0.5 * 6), sprintf(' %0.6f ',reg(regstri).rmsgdt(end)));
            fprintf('\n'); 
          end
          
          if mdisplay(4)
            % this work very well
            fprintf('%30s','Template Matching: '); 
            QMC   = cat_io_colormaps('marks+',17);
            color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
            if mdisplay(4)>1
              for dti=1:size(reg(regstri).ll,1)-1, 
                cat_io_cprintf( color(QMC, (reg(regstri).ll(dti,1) - 0.05) / 0.15 * 6),sprintf('%0.3f ',reg(regstri).ll(dti,1))); 
              end
              fprintf('| ');
            end
            cat_io_cprintf( color(QMC,(reg(regstri).ll(end,1) - 0.05)/0.15 * 6), sprintf(' %0.6f ',reg(regstri).ll(end,1)));
            fprintf('\n'); 
          end
          
          reg(regstri).cbr  = diff(reg(regstri).relrmsdtc(1:end-1)) ./ -diff(reg(regstri).ll(1:end-1,1)');
          reg(regstri).scbr = sum(diff(reg(regstri).relrmsdtc(1:end-1)) ./ -diff(reg(regstri).ll(1:end-1,1)')); 
          reg(regstri).mcbr = mean(diff(reg(regstri).relrmsdtc(1:end-1)) ./ -diff(reg(regstri).ll(1:end-1,1)')); 
          if mdisplay(5)
            % this work very well
            fprintf('%30s','Cost Benefit Ratio (CBR): ');
            QMC   = cat_io_colormaps('marks+',17);
            color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
            if mdisplay(5)>1
              for dti=1:size(reg(regstri).cbr,2), 
                cat_io_cprintf( color(QMC, (reg(regstri).cbr(dti))/2),sprintf('%0.3f ',reg(regstri).cbr(dti))); 
              end
              fprintf('| ');
            end
            cat_io_cprintf( color(QMC,( reg(regstri).mcbr)/2), sprintf(' %0.6f ', reg(regstri).mcbr));
            fprintf('\n'); 
          end
          
          
          % write xml
          if export
            cat_io_xml(fullfile(pth,mrifolder,testfolder,['reg_', nam, '.xml']),reg(regstri))
          end  
        end
          
          
        if export && ~dreg.fast
          % write output
          stime = cat_io_cmd(sprintf('Write Output with %0.2f mm',job.extopts.vox(voxi)));

                   
          
          % tissue ouptut
          fn = {'GM','WM','CSF'};
          if exist('Yclso','var'), Ycls = Yclso; end
          for clsi=1%max(1,2*(export-1))
            if export>1
              %% template creation
              cat_io_writenii(VT0,single(Ycls{clsi})/255,fullfile(mrifolder,testfolder),sprintf('p%d',clsi),...
                sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],[0 0 0 2],trans);
            end
          end
          %%
         for clsi=1:2
            cat_io_writenii(VT0,single(Ycls{clsi})/255,fullfile(mrifolder,testfolder),sprintf('p%d',clsi),...
              sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],[0 0 0 0],trans);
            cat_io_writenii(VT0,single(Ycls{clsi})/255,fullfile(mrifolder,testfolder),sprintf('p%d',clsi),...
              sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],[0 0 1 0],trans);
         end
           %%
          
            % hier ist kein unterschied per definition ...  
            %cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
            %  sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],[0 0 2 0],trans);  


          
          %% write jacobian determinant
          if job.extopts.regstr(regstri)>0 % shooting
            if debug, dt2o=dt2; end %#ok<NASGU>
            dx = 10; % smaller values are more accurate, but large look better; 
            [D,I] = cat_vbdist(single(~(isnan(dt2) | dt2<0 | dt2>100) )); D=min(1,D/min(dx,max(D(:)))); 
            dt2 = dt2(I); dt2 = dt2 .* ((1-D) + D); 
            dt2(isnan(dt2))=1; 
          else %dartel
            [y0, dt2] = spm_dartel_integrate(reshape(trans.jc.u,[trans.warped.odim(1:3) 1 3]),[1 0], 6); %#ok<ASGLU>
            clear y0
          end
          
          % create nifti
          N         = nifti;
          N.dat     = file_array(fullfile(pth,mrifolder,testfolder,['wj_', nam, '.nii']),trans.warped.odim(1:3),...
                      [spm_type('float32') spm_platform('bigend')],0,10/256^2,0);
          N.mat     = M1;
          N.mat0    = M1;
          N.descrip = ['Jacobian' VT0.descrip];
          create(N);
          N.dat(:,:,:) = dt2;

          
          cat_io_cmd('','',''); cat_io_cmd('','','',job.extopts.verb,stime); 


          %% deformations y - dartel > subject
          if job.output.warps(1)
              Yy2       = spm_diffeo('invdef',trans.warped.yx,trans.warped.odim,eye(4),trans.warped.M0);
              N         = nifti;
              N.dat     = file_array(fullfile(pth,mrifolder,testfolder,['y_', nam, '.nii']),[trans.warped.odim(1:3),1,3],'float32',0,1,0);
              N.mat     = trans.warped.M1;
              N.mat0    = trans.warped.M1;
              N.descrip = 'Deformation';
              create(N);
              N.dat(:,:,:,:,:) = reshape(Yy2,[trans.warped.odim,1,3]);
              clear Yy2; 
          end

         
          %% deformation iy - subject > dartel
          if job.output.warps(2)
            if any(trans.native.Vo.dim~=trans.native.Vi.dim)
              %%
              vx_voli  = sqrt(sum(trans.native.Vi.mat(1:3,1:3).^2));  
              vx_volo  = sqrt(sum(trans.native.Vo.mat(1:3,1:3).^2));
              eyev = eye(4); eyev([1 6 11]) = eyev([1 6 11]) .* vx_volo./vx_voli; 
              Yy2  = zeros([trans.native.Vo.dim 1 3],'single');                        
              for k1=1:3
                Yy2nl = zeros(trans.native.Vo.dim,'single');    
                for i=1:trans.native.Vo.dim(3),
                  Yy2(:,:,i,:,k1) = trans.warped.M1(k1,4) + trans.warped.M1(k1,k1) * ...
                    single(spm_slice_vol(trans.warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]), ...
                    trans.native.Vo.dim(1:2),[1,NaN])); % adapt for res
                  Yy2nl(:,:,i) = trans.warped.M1(k1,4) + trans.warped.M1(k1,k1) * ...
                    single(spm_slice_vol(trans.warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]), ...
                    trans.native.Vo.dim(1:2),[5,NaN])); % adapt for res
                end
                Yy2 = combine_linear_spline(Yy2,Yy2nl,eyev,k1); clear Yy2nl; 
              end
            else 
              yn     = numel(trans.warped.y); 
              p      = ones([4,yn/3],'single'); 
              p(1,:) = trans.warped.y(1:yn/3);
              p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
              p(3,:) = trans.warped.y(yn/3*2+1:yn);
              p      = M1(1:3,:) * p;

              Yy2 = zeros([trans.native.Vo.dim(1:3),1,3],'single'); 
              Yy2(1:yn/3)        = p(1,:);
              Yy2(yn/3+1:yn/3*2) = p(2,:);
              Yy2(yn/3*2+1:yn)   = p(3,:);
            end
            clear p; 

            % f2 = spm_diffeo('resize', f1, dim)
            % write new output
            Ndef      = nifti;
            Ndef.dat  = file_array(fullfile(pth,mrifolder,testfolder,['iy_', nam, '.nii']),[trans.native.Vo.dim,1,3],...
                        [spm_type('float32') spm_platform('bigend')],0,1,0);
            Ndef.mat  = res.image0(1).mat;
            Ndef.mat0 = res.image0(1).mat;
            Ndef.descrip = 'Inverse Deformation';
            create(Ndef);
            Ndef.dat(:,:,:,:,:) = Yy2;
            clear Yy2;
          end
        end
      end
    end
  end

  % back to old directory 
  if exist('olddir','var'), cd(olddir); end
end
%=======================================================================
function lin = combine_linear_spline(lin,spl,mat,k1)  
% Function to combine a linear and spline interpolated version of an image
% that replace the bad values of the spline image close to the image 
% boundary by the simple linear version. 

% Moreover it replace NaNs by close image values.

  repnan = 1; % this has to be more adaptiv to avoid lines in the Determinant
  
  % helping function to find the image boundary and add some space
  mskid = @(Y,r) max( ones(1,6) + round(1/r),min( ...
    [ size(Y,1) size(Y,1) size(Y,2) size(Y,2) size(Y,3) size(Y,3)] - round(1/r),[ 
    find( ~isnan(Y(:,round(size(Y,1)/2),round(size(Y,3)/2))),1,'first') + round(1/r)  ...
    find( ~isnan(Y(:,round(size(Y,2)/2),round(size(Y,3)/2))),1,'last')  - round(1/r)  ...
    find( ~isnan(Y(round(size(Y,1)/2),:,round(size(Y,3)/2))),1,'first') + round(1/r)  ...
    find( ~isnan(Y(round(size(Y,1)/2),:,round(size(Y,3)/2))),1,'last')  - round(1/r)  ... 
    find( ~isnan(Y(round(size(Y,1)/2),round(size(Y,2)/2),:)),1,'first') + round(1/r)  ...
    find( ~isnan(Y(round(size(Y,1)/2),round(size(Y,2)/2),:)),1,'last')  - round(1/r)  ]));
  
  
  % create mask that include to good non-linear interpolated values and
  % avoid the values close to the image edge.
  vx_vol = sqrt(sum(mat(1:3,1:3).^2));    
  try
    mskidv = mskid(spl,abs(vx_vol(1)/2)); % amout of interpolation * 2 (inverse becuase of reduction)
  catch
    mskidv = [1 size(spl,1) 1 size(spl,2) 1 size(spl,3)]; 
  end
  msk = false(size(spl)); 
  msk(mskidv(1):mskidv(2),mskidv(3):mskidv(4),mskidv(5):mskidv(6)) = true; 
  
  % combin both images
  if exist('k1','var')
    if isnumeric(lin)
      if ndims(lin)==4
         % set NaN to closest value
         if repnan
          [D,I] = cat_vbdist(single(~(isnan(lin(:,:,:,k1))))); clear D;  %#ok<ASGLU>
          lint = lin(:,:,:,k1); 
          lint = lint(I); spl = spl(I); clear I; 
          lin(:,:,:,k1) = lint; clear lint;
         end
         
        lin(:,:,:,k1)   = lin(:,:,:,k1) .* (1-msk) + spl .* msk; 
        
        lini = lin(:,:,:,k1); 
        lini(~isfinite(lini))=lini(mod(find(~isfinite(lini))+1,numel(lini)));
        lin(:,:,:,k1) = lini;  
      elseif ndims(lin)==5
        if repnan
          [D,I] = cat_vbdist(single(~(isnan(lin(:,:,:,:,k1))))); clear D; %#ok<ASGLU> 
          lint = lin(:,:,:,:,k1); 
          lint = lint(I); spl = spl(I); clear I; 
          lin(:,:,:,:,k1) = lint; clear lint;
        end
        
        lin(:,:,:,:,k1) = lin(:,:,:,:,k1) .* (1-msk) + spl .* msk; 
        
        lini = lin(:,:,:,:,k1); 
        lini(~isfinite(lini))=lini(mod(find(~isfinite(lini))+1,numel(lini)));
        lin(:,:,:,:,k1) = lini;  
      end
     
    else
      if repnan
        [D,I] = cat_vbdist(single(~(isnan(lin{k1})))); %#ok<ASGLU> 
        lin{k1} = lin{k1}(I); spl = spl(I); clear D I; 
      end
      
      lin{k1}(:,:,:) = lin{k1}(:,:,:) .* (1-msk) + spl .* msk; 
      
      lin{k1}(~isfinite(lin{k1}))=lin{k1}(find(~isfinite(lin{k1}))+1);
    end
  else
     % set NaN to closest value
    if repnan
      [D,I] = cat_vbdist(single(~(isnan(lin))));  %#ok<ASGLU>
      lin = lin(I); spl = spl(I); clear D I; 
    end
    
    lin(:,:,:) = lin(:,:,:) .* (1-msk) + spl .* msk; 
    
    lin(~isfinite(lin))=lin(mod(find(~isfinite(lin))+1,numel(lin)));
  end
  
end
%=======================================================================
function x = rgrid(d)
  x = zeros([d(1:3) 3],'single');
  [x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
  for i=1:d(3),
      x(:,:,i,1) = x1;
      x(:,:,i,2) = x2;
      x(:,:,i,3) = single(i);
  end
end
%=======================================================================
function y1 = affind(y0,M)
  y1 = zeros(size(y0),'single');
  for d=1:3,
      y1(:,:,:,d) = y0(:,:,:,1)*M(d,1);
      y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,2)*M(d,2);
      y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,3)*M(d,3) + M(d,4);
  end
end
%=======================================================================
function [x1,y1,z1] = defs2(sol,z,M,prm,x0,y0,z0)
  iM = inv(M);
  z01 = z0(z)*ones(size(x0));

  x1a  = iM(1,1)*x0 + iM(1,2)*y0 + iM(1,3)*z01 + iM(1,4);
  y1a  = iM(2,1)*x0 + iM(2,2)*y0 + iM(2,3)*z01 + iM(2,4);
  z1a  = iM(3,1)*x0 + iM(3,2)*y0 + iM(3,3)*z01 + iM(3,4);

  x1 = spm_bsplins(sol{1},x1a,y1a,z1a,prm);
  y1 = spm_bsplins(sol{2},x1a,y1a,z1a,prm);
  z1 = spm_bsplins(sol{3},x1a,y1a,z1a,prm);
end
%=======================================================================
