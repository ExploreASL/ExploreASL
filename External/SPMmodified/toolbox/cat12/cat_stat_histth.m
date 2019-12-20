function varargout = cat_stat_histth(src,percent,opt)
% ______________________________________________________________________
% Remove outliers based on the histogram and replace them by the new
% limits. E.g. some MRI images have some extremely high or low values 
% that can trouble other functions that try to work on the full given 
% input range of the data. Removing only 0.2% of the data often often
% helps to avoid problems without removing important informations. 
% 
% The function can also print a histogram, box- or violin plot of 
% the given data (using cat_plot_boxplot) and give some basic values.
%
% [res,ths] = cat_stat_histth(src,percent,verb)
%  
%   src     .. input data
%   res     .. limited input data 
%   ths     .. estimated thresholds
%   percent .. included values (default) = 0.998; 
%   verb    .. 0 - none 
%              1 - histogram
%              2 - histogram without boundary values
%              3 - box-plot
%              4 - violin-plot
%              5 - violin- and box-plot 
%
% Expert options: 
%   
%  [res,ths] = cat_stat_histth(src,percent, opt )
%
%   opt      .. structur with further fields
%    .verb   .. see above 
%    .fs     .. font size 
%    .hbins  .. bins for histogram estimation 
%    .vacc   .. limit number of elements used in the violin plot
%               e.g. vacc=100 means src(1:100:end)
%
% Examples: 
%   s=10; b = randn(s,s,s); cat_stat_histth(b,0.9,4);
%   s=10; b = rand(s,s,s);  cat_stat_histth(b,0.9,5);
% ______________________________________________________________________
% $Id: cat_stat_histth.m 1309 2018-04-23 14:19:28Z dahnke $


  %% check input
  if nargin==0, help cat_stat_histth; return; end
  if ~exist('src','var') || isempty(src); 
    varargout{1} = src;
    varargout{2} = nan(1,2); 
    return; 
  end
  
  if ~exist('opt','var'), opt = struct(); end
  if ~isstruct(opt)
    verb = opt; clear opt; opt.verb = verb; 
  end
  def.verb        = 0;
  def.fs          = 16; 
  def.hbins       = 10000; 
  def.vacc        = max(1,min(100000,round(numel(src)/1000))); % reduce elements in violin plot
  opt = cat_io_checkinopt(opt,def); 

  if nargin==0, help cat_stat_histth; return; end
  if ~exist('percent','var') || isempty(percent)
    tol = 0.002; 
  else
    if percent<=1
      tol = 1 - percent; 
    elseif percent<=100
      tol = 1 - percent/100; 
    else 
      error('cat_stat_histth:percent','Percent has to be in the range of 1 to 100');  
    end
  end
  
  
  % histogram
  [hsrc,hval] = hist(src(~isinf(src(:)) & ~isnan(src(:)) & src(:)<3.4027e+38),opt.hbins);
  hp          = cumsum(hsrc)./sum(hsrc); 

  
  % lower limit
  if opt.verb, srco=src; end
  if min(src(:))~=0
    ths(1) = hval(max([1,find(hp>tol,1,'first')])); 
    src(src<ths(1)) = ths(1); 
  else
    ths(1) = 0; 
  end
  % upper limit
  ths(2) = hval( min( [numel(hval) , find(hp<(1-tol),1,'last') ])); 
  src(src>ths(2)) = ths(2); 
  
  %% display
  if opt.verb
    % get figure
    fh = findobj('name','cat_stat_histth');
    if isempty(fh), figure('name','cat_stat_histth'); else figure(fh); clf; end
    
    % create main plot
    subplot('Position',[0.10 0.06 0.64 0.86]); 
    if opt.verb == 1
      hist(src(:),100); 
    elseif opt.verb == 2
      hist(src(src(:)>ths(1) & src(:)<ths(2)),100); 
    else
      src2 = src(src(:)>ths(1) & src(:)<ths(2));
      if opt.verb == 3
        src2     = src2(1:opt.vacc:end); 
        boxwidth = 0.8;
      else
        src2     = src2(:); 
        boxwidth = 1.0;
      end
      cat_plot_boxplot( { src2 } , struct('violin',opt.verb - 3,'boxwidth',boxwidth)); 
    end
    if opt.verb <= 2
      title('Histogram','fontsize',opt.fs); 
    elseif opt.verb == 3
      title('Boxplot','fontsize',opt.fs); 
    elseif opt.verb == 4
      title('Violinplot without boundary values','fontsize',opt.fs); 
    elseif opt.verb == 5
      title('Violin/Boxplot without boundary values','fontsize',opt.fs); 
    end
    set(gca,'fontsize',opt.fs);
    
    
    % print some values
    annotation('textbox',[.75 0.06 0.23 0.86],'String',sprintf([...
      'max:  %10.4f \nmin:  %10.4f \nmedian:%9.4f \nmean: %10.4f \nstd:  %10.4f \n\n' ...
      'lth:  %10.4f \nhth:  %10.4f \n\n' ...
      'lth_{99}:%10.4f \nhth_{99}:%10.4f \n\n' ...
      'lth_{95}:%10.4f \nhth_{95}:%10.4f \n\n' ...
      ],[max(srco(:)),min(srco(:)),cat_stat_nanmedian(srco(:)),...
      cat_stat_nanmean(srco(:)),cat_stat_nanstd(srco(:))],...
      ths, ...
      [hval(find(hp>0.005,1,'first')),hval(find(hp<0.995,1,'last'))], ...
      [hval(find(hp>0.025,1,'first')),hval(find(hp<0.975,1,'last'))]), ...
      'FitBoxToText','On','fontsize',opt.fs*0.85,'linestyle','none','fontname','fixedwidth'); 
  end
  
  if nargout>=1, varargout{1} = src; end
  if nargout>=2, varargout{2} = ths; end
end