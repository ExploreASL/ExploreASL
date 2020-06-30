function varargout = cat_vol_sanlm(varargin)
% Spatial Adaptive Non Local Means (SANLM) Denoising Filter
%_______________________________________________________________________
% Filter a set of images and add the prefix 'sanlm_'.
% Missing input (data) will call GUI or/and use defaults. 
%
% Examples:
%   cat_vol_sanlm(struct('data','','prefix','n','rician',0));
%   cat_vol_sanlm(struct('data','','NCstr',[-1:0.5:1,inf,-inf));
%
% Input:
%  job - harvested job data structure (see matlabbatch help)
% 
% Output:
%  out - computation results, usually a struct variable.
%
%  cat_vol_sanlm(job)
%
%  job
%   .data             .. set of images 
%   .prefix           .. prefix for filtered images (default = 'sanlm_') 
%   .suffix           .. suffix for filtered images (default = '')
%                        'NCstr' will add the used parameter 
%   .verb             .. verbose processing (default = 1)
%   .spm_type         .. file datatype (default single = 16);
%   .replaceNANandINF .. replace NAN by 0, -INF by minimum and INF by maximum
%   .rician           .. noise distribution
%   .intlim           .. intensity limitation (default = 0.9999)
%   .addnoise         .. add noise to noiseless regions
%     Add minimal amount of noise in regions without any noise to avoid
%     problems of image segmentation routines. The value defines the 
%     strength of the noise by the percentage of the mean signal intensity. 
%   .NCstr            .. strength of noise correction (default = -inf) 
%     A value of 1 used the full correction by the SANLM filter. Values 
%     between 0 and 1 mix the original and the filtered images, whereas
%     INF estimated a global value depending on the changes of the SANLM
%     filter that reduce to strong filtering in high quality data.
%     Negative values work on a local rather than a global level. 
%     A value of -INF is recommend, but you define a set of values (see
%     example) for further changes depending on your data. 
%     
%      0                  .. no denoising 
%      1                  .. full denoising (original sanlm)
%      2                  .. "light":   NCstr=-0.5, red=0, fred=0; miter=0 
%      3 | -inf           .. "medium":  NCstr=-1.0, red=1, fred=0; miter=0
%      4                  .. "strong":  NCstr= 1.0, red=1, fred=1; miter=1
%      0 < job.NCstr < 1  .. automatic global correction with user weighting
%     -9 < job.NCstr < 0  .. automatic local correction with user weighting
%      inf                .. global automatic correction
%     -inf                .. local automatic correction
%
%   .returnOnlyFilename   .. just to get the resulting filenames for SPM 
%                            batch mode (default = 0)
%   .resolutionDependency .. resolution depending filter strength 
%     Use the .resolutionDependencyRange Parameter (default = 1)
%   .resolutionDependencyRange .. [full-correction no-correction]
%     Limit the filter size depending on the general brain size, where 
%     filtering of images with 2.5 mm voxel size and higher will remove 
%     important anatomical information (default = [1 2.5]).    
%   .red                  .. low resolution filtering (if high-res data)
%   .fred                 .. force resolution reduction 
%   .iter                 .. additional iterations on the reduced resolution
%                            (default = 0)
%   .miter                .. additional main iterations of the full filter
%                            (default = 0)
%
% Some MR images were interpolated or use a limited frequency spectrum to 
% support higher spatial resolution with acceptable scan-times 
% (eg. 0.5x0.5x1.5 mm on a 1.5 Tesla scanner). However, this can result in
% "low-frequency" noise that can not be handled by the standard NLM filter.
% Hence, an additional filtering step is used on a reduces resolution that
% uses an internal call of this routine with direct image in- an output.
% 
%   src = cat_vol_sanlm(job,V,i,src)
%
% As far as filtering of low resolution data will also remove anatomical
% information the filter uses by default maximal one reduction with a  
% resolution limit of 1.6 mm. I.e. a 0.5x0.5x1.5 mm image is reduced 
% to 1.0x1.0x1.5 mm, whereas a 0.8x0.8x0.4 mm images is reduced to 
% 0.8x0.8x0.8 mm and a 1x1x1 mm dataset is not reduced at all. 
%
%_______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id: cat_vol_sanlm1445.m 1470 2019-05-29 11:49:16Z gaser $

  
    
    if nargin == 0 
        varargin{1}.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
        if isempty(char(varargin{1}.data)); return; end
    end
    
    % default optoins
    def.verb                        = 2;         % be verbose
    def.prefix                      = 'sanlm_';  % prefix
    def.suffix                      = '';        % suffix
    def.replaceNANandINF            = 1;         % replace NAN and INF
    def.spm_type                    = 16;        % file datatype (default single)
    def.NCstr                       = -Inf;      % 0 - no denoising, eps - light denoising, 
                                                 % 1 - maximum denoising, inf = auto; 
    def.rician                      = 0;         % use inf for GUI
    def.intlim                      = 0.9999;    % general intensity limitation to remove strong outlier
    def.resolutionDependency        = 1;         % resolution depending filter strength
    def.resolutionDependencyRange   = [1 2.5];   % [full-correction no-correction]
    def.relativeIntensityAdaption   = 1;         % use intensity to limit relative corrections (0 - none, 1 - full)
    def.relativeIntensityAdaptionTH = 1;         % larger values for continuous filter strength
    def.relativeFilterStengthLimit  = 1;         % limit the noise correction by the relative changes 
                                                 % to avoid over-filtering in low intensity regions
    def.outlier                     = 1;         % threshold to define outlier voxel to filter them with full strength
    def.addnoise                    = 1;         % option to add a minimal amount of noise in regions without noise
    def.returnOnlyFilename          = 0;         % just to get the resulting filenames for SPM batch mode
    def.red                         = 1;         % number of reductions (be careful using values greater 1!)
    def.fred                        = 0;         % force reduce
    def.iter                        = 0;         % additional inner iterations on the reduced resolution
    def.miter                       = 0;         % additional main iterations of the full filter
    varargin{1} = cat_io_checkinopt(varargin{1},def);
    
    % special cases of the CAT GUI
    if isfield(varargin{1},'nlmfilter') 
      if isfield(varargin{1}.nlmfilter,'optimized')
        varargin{1} = cat_io_checkinopt(varargin{1}.nlmfilter.optimized,varargin{1});
      elseif isfield(varargin{1}.nlmfilter,'expert')
        varargin{1} = cat_io_checkinopt(varargin{1}.nlmfilter.expert,varargin{1});
      end
    end
    switch varargin{1}.NCstr
      case 2,        varargin{1}.NCstr =  -0.5; varargin{1}.red = 0; varargin{1}.fred = 0; varargin{1}.miter = 0; % light
      case {3,-inf}, varargin{1}.NCstr =  -1.0; varargin{1}.red = 1; varargin{1}.fred = 0; varargin{1}.miter = 0; % medium
      case 4,        varargin{1}.NCstr =   1.0; varargin{1}.red = 1; varargin{1}.fred = 1; varargin{1}.miter = 1; varargin{1}.iter = 1;% strong
    end 
    
    if nargin <= 1 && isstruct(varargin{1}) % job structure input
        if nargout>0
          varargout = cat_vol_sanlm_file(varargin{1});
        else
          cat_vol_sanlm_file(varargin{1});
        end
    else % image input
        varargout{1} = cat_vol_sanlm_filter(varargin{:});
    end

end

%_______________________________________________________________________
function varargout = cat_vol_sanlm_file(job)
    SVNid = '$Rev: 1470 $';
    
    if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
    else
       job.data = cellstr(job.data);
    end
    if isempty(char(job.data)); if nargout>0, varargout{1} = {{''}}; end; return; end


    % map GUI data
    if isfield(job,'nlmfilter') 
        if isfield(job.nlmfilter,'classic') 
            job.NCstr = 1; 
        elseif isfield(job.nlmfilter,'optimized') 
            FN = fieldnames(job.nlmfilter.optimized); 
            for fni=1:numel(FN)
                if isfield(job.nlmfilter.optimized,FN{fni})
                    job.(FN{fni}) = job.nlmfilter.optimized.(FN{fni}); 
                end
            end
            job.NCstr = -abs(job.nlmfilter.optimized.NCstr); 
        end
    end

    % parameter limitations
    if ~isinf(job.NCstr), job.NCstr = max(-9.99,min(1,job.NCstr)); end        % guarantee values from -9.99 to 1 or inf
    job.resolutionDependency        = max(0,min(9.99,job.resolutionDependency));
    job.relativeIntensityAdaption   = max(0,min(9.99,job.relativeIntensityAdaption));
    job.relativeIntensityAdaptionTH = max(0,min(9.99,job.relativeIntensityAdaptionTH));
    job.relativeFilterStengthLimit  = max(0,min(9.99,job.relativeFilterStengthLimit));
    job.outlier                     = max(0,min(9.99,job.outlier));
    job.addnoise                    = max(0,min(9.99,job.addnoise));

    % create automatic filenames by the parameter
    if ~isempty(strfind(job.prefix,'PARA'))
        nprefix = strrep(job.prefix,'PARA',''); 
        if numel(nprefix)>0 && ~(strcmp(nprefix(end),'_') || strcmp(nprefix(end),'-')), nprefix = [nprefix '_']; end
        if job.NCstr>=0
            job.prefix = sprintf('%sNC%0.2f_',nprefix,job.NCstr); 
        elseif isinf(job.NCstr) && sign(job.NCstr)==-1
            job.prefix = sprintf('%sNC%0.2f_',nprefix,3); 
        else
            job.prefix = sprintf('%sNC%-0.2f_RN%d_RD%d_RIA%0.2f_RR%d_FR%d_RNI%d_OL%0.2f_PN%0.1f_iterm%d_iter%d',...
                nprefix, job.NCstr , job.rician , job.resolutionDependency , job.red , job.fred , job.relativeIntensityAdaption , ...
                job.replaceNANandINF , job.outlier , job.addnoise , job.iterm , job.iter); 
        end
    elseif ~isempty(strfind(job.suffix,'PARA'))
        if numel(job.NCstr)==1 && job.NCstr>=0 
            job.suffix = sprintf('_NC%0.2f',job.NCstr); 
        elseif isinf(job.NCstr) && sign(job.NCstr)==-1
            job.prefix = sprintf('_NC%0.2f',3); 
        else
            job.suffix = sprintf('_NC%-0.2f_RN%d_RD%d_RIA%0.2f_RR%d_FR%d_RNI%d_OL%0.2f_PN%0.1f_iterm%d_iter%d',...
                job.NCstr , job.rician , job.resolutionDependency , job.red , job.fred , job.relativeIntensityAdaption , ...
                job.replaceNANandINF , job.outlier , job.addnoise , job.iterm , job.iter); 
        end
    end
    % just to get the resulting filenames for SPM batch mode
    if job.returnOnlyFilename
        for i = 1:numel(job.data)
            [pth,nm,xt,vr]  = spm_fileparts(deblank(job.data{i})); 
            varargout{1}{i} = fullfile(pth,[job.prefix nm job.suffix xt vr]);
        end
        return
    end

    V  = spm_vol(char(job.data));

    % new banner
    if isfield(job,'process_index') && job.verb, spm('FnBanner',mfilename,SVNid); end
    spm_clf('Interactive'); 
    spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');

    for i = 1:numel(job.data)
        cat_vol_sanlm_filter(job,V,i);
    end

    if isfield(job,'process_index') && job.verb, fprintf('Done\n'); end
    spm_progress_bar('Clear');
end

%_______________________________________________________________________
function src2 = cat_vol_sanlm_filter(job,V,i,src)
    Vo = V;

    QMC   = cat_io_colormaps('marks+',17);
    color = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);

    % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
    dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end


    [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname)); %#ok<ASGLU>

    stime = clock; stimef = clock;
    vx_vol  = sqrt(sum(V(i).mat(1:3,1:3).^2));

    if nargin<4
      src = single(spm_read_vols(V(i)));
    else
      src = single(src);
    end
        
		for im=1:1+job.miter;   
			% prevent NaN and INF
			if job.replaceNANandINF
				src(isnan(src)) = 0;
				src(isinf(src) & src<0) = min(src(:));
				src(isinf(src) & src>0) = max(src(:));
			else
				if sum(isnan(src(:))>0), nanmsk = isnan(src); end
				if sum(isinf(src(:))>0), infmsk = int8( isinf(src) .* sign(src) ); end
			end
			
			% histogram limit
			src = cat_stat_histth(src,job.intlim); 
			
			% use intensity normalisation because cat_sanlm did not filter values below ~0.01 
			th  = max( cat_stat_nanmean( src(src(:)>cat_stat_nanmean(src(src(:)>0))) ) , ...
								 abs(cat_stat_nanmean( src(src(:)<abs(cat_stat_nanmean(src(src(:)<0)))))) );
							 
			% low resolution filtering
			if job.red>0 && (any(vx_vol<0.8) || job.fred )
				[srcr,resr]   = cat_vol_resize( src ,'reduceV',vx_vol,min(2.2*(job.fred+1),min(vx_vol)*2.3),32,'median');
				
				jobr            = job;
				jobr.red        = job.red - 1;     
				jobr.miter      = 0; 
				jobr.addnoise   = 0;                        % no additional noise on lower resolution 
				jobr.resolutionDependency      = 1;         % resolution depending filter strength
				jobr.resolutionDependencyRange = [1 1.6];   % [full-correction no-correction]
				jobr.outlier    = jobr.outlier*2; 
				jobr.NCstr      = -prod(3-resr.vx_red)*2 .* ...                 
													min(1,1 - (mean(resr.vx_volr) - jobr.resolutionDependencyRange(1) )) / ...
													diff(jobr.resolutionDependencyRange); 
	
				if 1
					% larger var => more information
					Ygr   = cat_vol_grad(srcr/th,resr.vx_volr,0); 
					grsd  = std(Ygr(Ygr(:)>0)); 
					grrel = numel(Ygr(Ygr(:)>grsd))/numel(Ygr(Ygr(:)>0)); 
					jobr.NCstr = jobr.NCstr * grrel*10; 
				end
	
				srco = src;
				if any(resr.vx_red>1) && any( resr.vx_volr < jobr.resolutionDependencyRange(2)*(job.fred+1) ) && jobr.NCstr~=0
					% first block
					Vr = V(i); Vmat = spm_imatrix(Vr.mat); Vmat(7:9) = Vmat(7:9).*resr.vx_red; Vr.mat = spm_matrix(Vmat);
					srcR  = cat_vol_resize(srcr,'dereduceV',resr,'cubic');
					for iter=1:1+job.iter
						srcr  = cat_vol_sanlm_filter(jobr,Vr,1,srcr);
					end
					srcRS = cat_vol_resize(srcr,'dereduceV',resr,'cubic');
					src   = (src - srcR) + srcRS; 
					clear srcRS srcr srcR;
	
					% second block
					if 1 % the second displaced filtering helps to deduce low frequency noise a little bit more 
						[srcr,resr] = cat_vol_resize( smooth3(src(2:end,2:end,2:end)) ,'reduceV',...
							vx_vol,min(2.2*(job.fred+1),min(vx_vol)*2.3),32,'median');
						srcRr = cat_vol_resize(srcr,'dereduceV',resr,'cubic');
						srcR  = src; srcR(2:end,2:end,2:end) = srcRr; 
						for iter=1:1+job.iter
							srcr  = cat_vol_sanlm_filter(jobr,Vr,1,srcr);
						end
						srcRr = cat_vol_resize(srcr,'dereduceV',resr,'cubic');
						srcRS = src; srcRS(2:end,2:end,2:end) = srcRr; 
						src   = src + (src - srcR) + srcRS; 
						clear srcRS srcr srcR;
						src = src / 2;
					end
				end
				
				NCstrr = 15 * abs(cat_stat_nanmean(abs(src(:)/th - srco(:)/th))); 
			else
				NCstrr = 0; 
			end
			if job.verb>1 && (nargin>3 || NCstrr>0) 
				cat_io_cprintf('g5',sprintf('R%1d) %0.2fx%0.2fx%0.2f mm:  ',job.red,vx_vol)); stime = clock; 
			end
			
			
			% the real noise filter
			if im==1, srco = src; end
			src = (src / th) * 100; 
			cat_sanlm(src,3,1,job.rician);
			src = (src / 100) * th; 
			if job.verb>1 && (nargin>3 || NCstrr>0) && im<1+job.miter
				if job.verb>1 && (nargin>3 || NCstrr>0), cat_io_cprintf('g5',sprintf('  %5.0fs\n',etime(clock,stime))); end
			end
		end
			
    % measures changes
    NCrate = cat_stat_nanmean(abs(src(:)/th - srco(:)/th)); 
    
    % set actual filter rate - limit later!
    % the factor 15 was estimated on the BWP 
    NCstr                         = job.NCstr; 
    NCstr(isinf(NCstr) & NCstr>0) = 15 * NCrate;   
    NCstr(isinf(NCstr) & NCstr<0) = -1;   
   
   
    for NCstri = 1:numel(NCstr)
       
        if NCstr(NCstri)<0
        % adaptive local denoising 

            %% prepare local map
            % use less filtering for low-res data to avoid anatomical blurring ???
            NCs   = max(eps,abs(src - srco)/th); 
            
            % preserve anatomical details by describing the average changes
            % and not the strongest - this reduce the ability of artifact
            % correction!
            stdNC = std(NCs(NCs(:)~=0)); 
            NCsm  = cat_vol_median3(NCs,NCs>stdNC,true(size(NCs))); % replace outlier
            [NCsr,resT2] = cat_vol_resize(NCsm,'reduceV',vx_vol,2,32,'meanm'); clear NCsm; 
            NCsr  = cat_vol_localstat(NCsr,true(size(NCs)),1,1);
            NCsr  = cat_vol_smooth3X(NCsr,1/mean(resT2.vx_volr));
            NCsr  = cat_vol_resize(NCsr,'dereduceV',resT2); 
            NCso  = NCs;  
            
            % no correction of local abnormal high values (anatomy)
            NCs = NCsr + (NCso>stdNC & NCso<=stdNC*4 & NCso>NCsr*2 & NCso<NCsr*16) .* (-NCsr); 
            NCs = cat_vol_smooth3X(NCs,2/mean(resT2.vx_vol));  
            % weighting
            NCs = NCs ./ mean(NCs(:)); NCs = max(0,min(1,NCs * (15*NCrate) * abs(NCstr(NCstri)))); 
            
            
            
            % Volume dependency with lower filtering in images with lower resolution. 
            % Even low noise correction in structural images with low resolution
            % can remove important anatomical details. Hence, filter strength is 
            % here linear adapted by the average voxel resolution.
            % For volumeDependencyRange = [1 2.5]  images with 1 mm or better are fully filtered, 
            % whereas images with 2.5 mm or lower resolution are not filtered
            if job.resolutionDependency
              NCs = NCs .* max( 0 , min( 1 , 1 - ...
                ( mean(vx_vol) - job.resolutionDependencyRange(2) ) / ...
                diff(job.resolutionDependencyRange) ) ); 
            end

            
                    
            %% relative changes
            %  The SANLM filter is often very successful in the background
            %  and removed nearly all noise. However, routines such as the
            %  SPM Unified Segmentation expect Gaussian distribution in all
            %  regions and is troubled by regions with too low variance. 
            %  Hence, a relative limitation of SANLM correction is added 
            %  here that is based on the bias reduced image intensity. 
            %
            %  NCi defines the normalized intensity to estimate the relative change rate.
            %  The log10 function is used to reduce bias effects.
            %  A small resolution dependent (anatomical) blurring is used to avoid artifacts
            %  but large values would be counterproductive (job.relativeIntensityAdaption). 
            %  The average change rate in signal region is estimated (mNCs) and used as limit
            %  for corrections (job.relativeFilterStengthLimit), where higher values allow 
            %  stronger filtering. 
            [NCi,range]  = cat_stat_histth(src,0.99); % lower values > more similar filtering
            NCi  = max(eps,log10( 1 + (NCi + range(1)) / diff(range) * 7 + 3 )); % bias corr + intensity normalization 
            NCi  = cat_vol_smooth3X( NCi , job.relativeIntensityAdaption / mean(vx_vol)); % smoothing
            if job.relativeIntensityAdaption>0 && ...
               job.relativeFilterStengthLimit && ~isinf(job.relativeFilterStengthLimit)
                NCsi = NCs ./ max(eps,NCi); 
                mNCs = cat_stat_nanmean( NCsi(src(:)>th/2 & NCsi(:)>0 )) * ...
                          job.relativeFilterStengthLimit * ...
                          max(1,min(4,4 - job.relativeIntensityAdaption*2)); % lower boundary for strong adaption
                NCsi = min( NCsi , mNCs ) .* NCi; 
                
                % Finally, both images were mixed
                NCs  = NCs  * (1 - job.relativeIntensityAdaption) + ... % linear average model to contoll filter strength
                       NCsi *      job.relativeIntensityAdaption; 
                if ~debug, clear NCsi; end
            end  
            

            
            % heavy outlier / artifacts
            if job.outlier>0
                NCi = min(1,max(0,NCi .* (NCso - ( (stdNC*2) / job.outlier ) ) ./ ((stdNC*2) / job.outlier ))); 
                NCs = max(NCs, NCi); 
                if ~debug, clear NCi; end
                if ~debug, clear NCso; end
            end  
            
            
            
            if debug, src2 = srco.*(1-NCs) + src.*NCs; end
            % ds('d2','',vx_vol,src/th,srco/th,srco2/th, NCs,160)

            
            
            % mix original and noise corrected image
            src2 = srco.*(1-NCs) + src.*NCs; 
            NCstr(NCstri) = -cat_stat_nanmean(NCs(:)); 
            if ~debug, clear NCs; end
            
        elseif NCstr(NCstri)==1
        % no adaption (original filter)
            src2   = src; 
        
        elseif NCstr(NCstri)>0
        % simple global denoising

            NCstr(NCstri) = min(1,max(0,NCstr(NCstri)));
            
            % mix original and noise corrected image
            src2   = srco*(1-NCstr(NCstri)) + src*NCstr(NCstri); 
            
        else
        % no denoising ... nothing to do
            src2 = src; 
            
        end

        
        
        %% add noise
        if job.addnoise
          % Small adaption for inhomogeneity to avoid too much noise in
          % regions with low signal intensity.
          sth  = cat_vol_smooth3X(log10(2 + 8*src/th),4/mean(vx_vol)) * th; 
          
          % Correction only of regions with less noise and with (src~=0) to 
          % avoid adding of noise in skull-stripped data. This may lead to
          % problems with the skull-stripping detection in cat_run_job!
          % Also important in case of ADNI.
          src2 = src2 + max( 0 , min(1 , cat_vol_smooth3X( ...
                 ( job.addnoise.*sth/100 ) - abs(srco - src) , 4/mean(vx_vol) ) ./ ( job.addnoise.*sth/100 ) )) .* ...
                 ( src~=0 ) .* ... save skull-stripping / defacing regions
                 (randn(size(src)) * job.addnoise.*sth/100);  
          if ~debug, clear sth; end
        end
        if numel(NCstr)==1 && ~debug, clear src srco; end
        
        
        
        %% restore NAN and INF
        if exist('nanmsk','var'), src2(nanmsk) = nan; end
        if exist('infmsk','var'), src2(infmsk==-1) = -inf; src2(infmsk==1) = inf; end
        if job.verb>1 && (nargin>3 || NCstrr>0), cat_io_cprintf('g5',sprintf('  %5.0fs\n',etime(clock,stime))); end
          
        if nargin==4
          return; 
        end
        
        % use only float precision
        Vo(i).fname   = fullfile(pth,[job.prefix nm job.suffix '.nii' vr]);
        Vo(i).descrip = sprintf('%s SANLM filtered (NCstr=%-4.2f > %0.2f)',...
          V(i).descrip,job.NCstr(NCstri),abs(NCstr(NCstri)) + NCstrr);
        Vo(i).dt(1)   = 16; % default - changes later if required 
        if exist(Vo(i).fname,'file'); delete(Vo(i).fname); end
        spm_write_vol(Vo(i), src2);
        spm_progress_bar('Set',i);
        
        
        
        %% if single should be not used, the image has to be converted ... 
        if job.spm_type~=16
          ctype.data   = Vo(i).fname; 
          
          if job.spm_type
            ctype.ctype  = job.spm_type; 
          else
            ctype.ctype  = V(i).dt(1); 
          end
          ctype.range  = 99.99;
          ctype.prefix = ''; 
          cat_io_volctype(ctype);
        end
        
        
        
        %% display result and link images for comparison
        if job.verb
            % I am not sure if this is intuitive. Maybe someone will think
            % that red means failed filtering ...
            %   green > low filtering 
            %   red   > strong filtering
            if NCstr(NCstri)>0 || isinf(NCstr(NCstri))
                fprintf('NCstr = '); 
                cat_io_cprintf( color( ( ( abs(NCstr(NCstri)) ) * 6 )) , ...
                    sprintf('%- 5.2f > %4.2f', job.NCstr(NCstri) , abs(NCstr(NCstri)) ));
                cat_io_cprintf( [0 0 0] , ', '); % restore default color!
            end
            
            % this is a long string but it loads the original and the filtered
            % image for comparison
            fprintf('%5.0fs, Output %s\n',etime(clock,stimef),...
              spm_file(Vo(i).fname,'link',sprintf(...
              ['spm_figure(''Clear'',spm_figure(''GetWin'',''Graphics'')); ' ...
               'spm_orthviews(''Reset''); ' ... remove old settings
               'ho = spm_orthviews(''Image'',''%s'' ,[0 0.51 1 0.49]); ',... top image
               'hf = spm_orthviews(''Image'',''%%s'',[0 0.01 1 0.49]);', ... bottom image
               'spm_orthviews(''Caption'', ho, ''original''); ', ... caption top image
               'spm_orthviews(''Caption'', hf, ''filtered''); ', ... caption bottom image
               'spm_orthviews(''AddContext'',ho); spm_orthviews(''AddContext'',hf); ', ... % add menu
               'spm_orthviews(''Zoom'',40);', ... % zoom in
              ],V(i).fname)));
        end
        
        
    end
end

