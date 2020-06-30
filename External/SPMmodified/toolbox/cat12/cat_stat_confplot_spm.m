% [signal_change, xyz] = cat_stat_confplot_spm(SPM,xSPM,hReg,names,Ic)
%
% SPM, xSPM, hReg - parameters saved in workspace
% name     - optional names of columns given as {'name1','name2'...}
% Ic       - number of contrast (usually 1 for effects of interest)   
%
% signal_change  - signal change
% xyz            - coordinates of local cluster maximum     
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_stat_confplot_spm.m 1563 2020-02-11 16:00:46Z gaser $

global xY Hc x

try
  [xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
  spm_XYZreg('SetCoords',xyz,hReg);
catch
  [hReg xSPM SPM] = spm_results_ui('Setup');
  [xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
  spm_XYZreg('SetCoords',xyz,hReg);
end

CI = 1.6449;         % = spm_invNcdf(1 - 0.05);

%-Colour specifications
%-----------------------------------------------------------------------
Col   = [0 0 0; .8 .8 .8; 1 0 0];

Cplot = 'Parameter estimates';

%-Specify VOI
%-----------------------------------------------------------------------
xY.def = spm_input('VOI definition...',1,'b',...
      {'sphere','box','cluster','voxel'},[],4);
Q = ones(1,size(xSPM.XYZmm,2));

% default font size and boxplot
if ~exist('Hc','var') | (exist('Hc','var') & isempty(Hc))
  Hc.FS = 18;          % font size
  Hc.LW = 2;           % line width
  Hc.MS = 10;          % marker size
  Hc.boxplot = 1;      % show boxplot
  Hc.medianplot = 0;   % show median line
  Hc.rawdata = 0;      % show raw data
  Hc.adjust = 1;       % adjust raw data
  Hc.connected = 0;    % show connected lines for repeated anova
  Hc.legend = 1;       % show legend
  Hc.line = cell(1,1); % trend line property 
end

if ~exist('colored','var')
    colored = spm_input('Boxplot','+1', 'm',['Colored|Define colors|Grey'],[2 1 0],1);
  switch colored
  case 0
    groupcolor = [0.7 0.7 0.7];
  case 1
    groupcolor = spm_input('Colors','!+0','r',[],[n_effects 3]);
  case 2
    groupcolor = [];
  end
else
  groupcolor = [];
end

switch xY.def

  case 'sphere'
  %---------------------------------------------------------------
  xY.spec = spm_input('VOI radius (mm)','!+0','r',0,1,[0,Inf]);
  d     = [xSPM.XYZmm(1,:) - xyz(1);
  xSPM.XYZmm(2,:) - xyz(2);
  xSPM.XYZmm(3,:) - xyz(3)];
  Q     = find(sum(d.^2) <= xY.spec^2);
  XYZstr = sprintf(' averaged in sphere (radius %d mm)', xY.spec);
  xY.string = sprintf('sphere_%dmm_at_%g_%g_%gmm',xY.spec,xyz);

  case 'box'
  %---------------------------------------------------------------
  xY.spec = spm_input('box dimensions [x y z] {mm}',...
      '!+0','r','0 0 0',3);
  Q     = find(all(abs(xSPM.XYZmm - xyz*Q) <= xY.spec(:)*Q/2));
  XYZstr = sprintf(' averaged in box dimensions (%3.2f %3.2f %3.2f)', xY.spec);
  xY.string = sprintf('box_%g_%g_%gmm_at_%g_%g_%gmm',xY.spec,xyz);

  case 'cluster'
  %---------------------------------------------------------------
  [x0 i] = spm_XYZreg('NearestXYZ',xyz,xSPM.XYZmm);
  A     = spm_clusters(xSPM.XYZ);
  Q     = find(A == A(i));
  XYZstr = sprintf(' averaged in cluster');
  xY.string = sprintf('cluster_at_%g_%g_%gmm',x0);

  case 'voxel'
  %---------------------------------------------------------------
  d     = [xSPM.XYZmm(1,:) - xyz(1);
  xSPM.XYZmm(2,:) - xyz(2);
  xSPM.XYZmm(3,:) - xyz(3)];
  d2 = sum(d.^2);
  Q = find(d2==min(d2));
  XYZstr = sprintf(' in voxel');
  xY.string = sprintf('voxel_at_%g_%g_%gmm',xyz);
end

XYZ     = xSPM.XYZ(:,Q);    % coordinates

%-Parameter estimates:   beta = xX.pKX*xX.K*y;
%-Residual mean square: ResMS = sum(R.^2)/xX.trRV
%---------------------------------------------------------------

beta0  = spm_get_data(SPM.Vbeta, XYZ);
beta   = mean(beta0,2);

try
  fprintf('Read raw data...');
  y = spm_get_data(SPM.xY.VY, XYZ);
  fprintf(sprintf('%s',repmat('\b',1,150)));
  fprintf(sprintf('%s',repmat(' ',1,150)));
  Hc.y_found = 1;
catch
  fprintf('No raw data found! Please check that you have not moved your data.\n');
  Hc.y_found = 0;
  try, close(Hc.h12); end
end

ResMS  = spm_get_data(SPM.VResMS,XYZ);
ResMS  = mean(ResMS,2);
Bcov   = ResMS*SPM.xX.Bcov;
Bcov   = Bcov;

% determine which contrast
%---------------------------------------------------------------
if ~exist('Ic','var')
  Ic = spm_input('Which contrast?','!+1','m',{SPM.xCon.name});
end

xCon = SPM.xCon(Ic);

TITLE = {Cplot XYZstr};

% find contrast and related columns in design matrix
%-------------------------------------------------------------- 
c0 = xCon.c;

[indi, indj] = find(c0~=0);
ind_X = unique(indi)';
X = SPM.xX.X;
X = X(:,ind_X);
n_effects = size(X,2);

covariate = 0;

% check for covariates
if ~isempty(SPM.xX.iC) & n_effects <= 2
  for i=1:n_effects
    % contrast is defined at entries of iC
    if ~isempty(find(ind_X(i) == SPM.xX.iC))
      covariate = 1;
    else
      covariate = 0;
    end
  end
end

c0 = c0(ind_X,:);

if ~exist('names','var')
  define_names = spm_input('Define names?',1,'yes|use numbers',[1 0],1);
  if define_names
    names = [];
    for i=1:n_effects
      new_name = spm_input(['Name for parameter ' num2str(i)],1,'s');
      names = strvcat(names,new_name);
    end
  else
    names = num2str((1:n_effects)');
  end
end

% compute contrast of parameter estimates and 90% C.I.
%-------------------------------------------------------------- 
signal_change0 = xCon.c'*beta0;
signal_change  = xCon.c'*beta;
CI = CI*sqrt(diag(xCon.c'*Bcov*xCon.c));

if ~exist('repeated_anova','var')
  repeated_anova = ~isempty(SPM.xX.iB);
  
  if repeated_anova
    [rw,cl] = find(SPM.xX.I == length(SPM.xX.iB)); % find column which codes subject factor (length(xX.iB) -> n_subj)
    
    % expect that subject factors are 2nd colum, group 3rd column, time 4th column
    if cl(1) == 2
      n_groups = max(SPM.xX.I(:,3));
      count = 0;
      for i=1:n_groups
        ind_times{i} = count + (1:max(SPM.xX.I(find(SPM.xX.I(:,3)==i),4)));
        count = count + max(SPM.xX.I(find(SPM.xX.I(:,3)==i),4));
      end
      n_time = max(SPM.xX.I(:,4));
      n_groupsxtime = n_groups*n_time;
      
      if n_groupsxtime ~= n_effects
        repeated_anova = 0;
      end
    else
      repeated_anova = 0;
    end
  end
end

% GUI figure
%--------------------------------------------------------------

% estimate maximum available window width for the two windows
Ms = spm('WinSize','0',1);
menu_width = 150;
max_width = floor((Ms(3) - menu_width)/2);
max_width = min(max_width, 800);

Hc.h10 = figure(10);
clf

set(Hc.h10,'Position',[0 800 menu_width 550],'MenuBar','none','NumberTitle','off');
hNewButton = uicontrol(Hc.h10,...
    'Position',[20 500 110 20],...
    'Callback','cat_stat_confplot_spm',...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Plot');
hClearButton = uicontrol(Hc.h10,...
    'position',[20 460 110 20],...
    'Callback','clear -globalvar Hc;clear names Ic colored groupcolor repeated_anova',...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Reset variables');
hSaveButton = uicontrol(Hc.h10,...
    'position',[20 420 110 20],...
    'Callback',{@save_image},...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Save images');
hCloseButton = uicontrol(Hc.h10,...
    'position',[20 380 110 20],...
    'Callback','clear -globalvar Hc;clear names Ic colored groupcolor repeated_anova,close(10,11,12)',...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Close windows');
htext1 = uicontrol(Hc.h10,...
    'position',[20 340 60 20],...
    'Style','Text',...
    'String','Font Size');
hedit1 = uicontrol(Hc.h10,...
    'position',[80 340 50 20],...
    'Callback',{@set_font_size},...
    'Interruptible','on',...
    'Style','Edit',...
    'String',num2str(Hc.FS));
htext2 = uicontrol(Hc.h10,...
    'position',[20 300 60 20],...
    'Visible','off',...
    'Style','Text',...
    'String','Line width');
hedit2 = uicontrol(Hc.h10,...
    'position',[80 300 50 20],...
    'Callback',{@set_line_width},...
    'Interruptible','on',...
    'Style','Edit',...
    'Visible','off',...
    'String',num2str(Hc.LW));
htext3 = uicontrol(Hc.h10,...
    'position',[20 270 60 20],...
    'Visible','off',...
    'Style','Text',...
    'String','Marker Size');
hedit3 = uicontrol(Hc.h10,...
    'position',[80 270 50 20],...
    'Callback',{@set_marker_size},...
    'Interruptible','on',...
    'Style','Edit',...
    'Visible','off',...
    'String',num2str(Hc.MS));
hAdjustData = uicontrol(Hc.h10,...
    'position',[20 240 110 20],...
    'Callback',(@adjust_data),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',Hc.adjust,...
    'ToolTip','Adjust raw data for potential covariates and subject effects',...
    'String','Adjust Raw Data');
hShowBoxplot = uicontrol(Hc.h10,...
    'position',[20 220 110 20],...
    'Callback',(@show_boxplot),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',Hc.boxplot,...
    'String','Show Boxplot');
hShowRawdata = uicontrol(Hc.h10,...
    'position',[20 200 110 20],...
    'Callback',(@show_rawdata),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',Hc.rawdata,...
    'String','Show Raw Data');
hShowMedianplot = uicontrol(Hc.h10,...
    'position',[20 180 110 20],...
    'Callback',(@show_medianplot),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',Hc.medianplot,...
    'String','Show Medianplot');
hShowConnected = uicontrol(Hc.h10,...
    'position',[20 160 110 20],...
    'Callback',(@show_connected),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',Hc.connected,...
    'String','Connect Lines');
hShowLegend = uicontrol(Hc.h10,...
    'position',[20 140 110 20],...
    'Callback',(@show_legend),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',Hc.legend,...
    'String','Show Legend');

if Hc.y_found
  set(hAdjustData,'Visible','on');
  
  if covariate
    set(htext2,'Visible','on');
    set(hedit2,'Visible','on');
    set(htext3,'Visible','on');
    set(hedit3,'Visible','on');

    if repeated_anova
      set(hShowConnected,'Visible','on');
      set(hShowLegend,'Visible','on');
    end
  else
    set(hShowBoxplot,'Visible','on');
    set(hShowRawdata,'Visible','on');

    if repeated_anova
      set(hShowMedianplot,'Visible','on');
      set(htext2,'Visible','on');
      set(hedit2,'Visible','on');
    end
  end
end

% % signal change plot
%--------------------------------------------------------------

if ~exist('Hc','var') | (exist('Hc','var') & ~isfield(Hc,'h11'))
  Hc.h11 = figure(11);
  set(Hc.h11,'Position',[menu_width 800 max_width 550],'NumberTitle','off','MenuBar','none');
else
  Hc.h11 = figure(11);
end

cla
hold on

% estimates
%--------------------------------------------------------------
h = bar(signal_change');
set(h,'FaceColor',Col(2,:));

% standard error
%--------------------------------------------------------------
for j = 1:length(signal_change)
  line([j j],([CI(j) 0 - CI(j)] + signal_change(j)),...
        'LineWidth',2,'Color',Col(3,:))
end

title(TITLE,'FontSize',14,'FontWeight','bold')
ylabel('parameter estimate','FontSize',12)
set(gca,'XLim',[0.4 (length(signal_change) + 0.6)],'XTick',1:length(signal_change));

if exist('names','var')
  if size(names,1) == length(signal_change)
    set(gca,'XTickLabel',names,'TickLabelInterpreter','none');
  end
end

hold off

if repeated_anova
  if n_groups > 3
    Hc.col = jet(n_groups);
  else
    Hc.col = [1 0 0;0 0 1; 0 1 0];
  end
else
  Hc.col = [1 0 0;0 0 1; 0 1 0];
end

% prepare raw values for boxplot
%--------------------------------------------------------------
if Hc.y_found

  % adjust raw data
  if Hc.adjust

		% define subject effects and potential covariates
		G_columns = [SPM.xX.iB SPM.xX.iC];
		
		% only consider nuisance parameters and parameters where
		% contrast is defined
		for i=1:n_effects
			G_columns(find(G_columns==ind_X(i))) = [];
		end
		
		% remove nuisance effects from data
		if ~isempty(G_columns)
			G = SPM.xX.X(:,G_columns);
			G = G - mean(G);
			y = y - G*(pinv(G)*y);
		end
  end
  
  % use mean inside cluster
  y = mean(y,2);
  
  y_label = 'raw signal';
  
  % estimate group means for correction for repeated anovas or interaction designs
  % expect that subject factors are 2nd colum, group 3rd column, time 4th column
  if Hc.adjust & (repeated_anova | ((n_groups > 1) & covariate))
		mean_group = zeros(n_groups,1);
		count_times = 1;
		for i=1:n_groups
			ind_group = find(SPM.xX.I(:,3) == i);
			if repeated_anova
			  % find subjects effects in that group
				ind_subj = unique(SPM.xX.I(ind_group,2));
				n_subj_group = numel(ind_subj);
				n_times = max(SPM.xX.I(ind_group,4));
				mean_group(i) = sum(beta(SPM.xX.iH(count_times:(count_times+n_times-1))))/n_times + ...
				  sum(beta(SPM.xX.iB(ind_subj)))/n_subj_group;
				count_times = count_times + n_times;
			else
				mean_group(i) = beta(SPM.xX.iH(i));
			end
			y(ind_group,:) = y(ind_group,:) - mean(y(ind_group,:)) + mean_group(i);
		end
  end
  
  yy = cell(n_effects,1);
  for i=1:n_effects
    yy{i} = y(find(X(:,i)~=0),:);
  end
  
  if ~exist('Hc','var') | (exist('Hc','var') & ~isfield(Hc,'h12'))
    Hc.h12 = figure(12);
    set(Hc.h12,'Position',[max_width+menu_width 800 max_width 550],'NumberTitle','off','MenuBar','none');
  else
    Hc.h12 = figure(12);
  end
  
  cla
  
  if Hc.rawdata
    vshowdata = 1;
  else
    vshowdata = 0;
  end
  
  if Hc.boxplot
    vbox = 1;
    voutliers = 1;
  else
    vbox = 0;
    voutliers = 0;
  end

  vstruct = struct('showdata',vshowdata,'box',vbox,'outliers',voutliers);
  if ~isempty(groupcolor)
    vstruct = setfield(vstruct,'groupcolor',groupcolor);
  end


  if isempty(yy{1}), return, end
  
  title_name = 'raw data ';
  if Hc.adjust, title_name = ['adjusted ' title_name]; end

  if covariate

    % previous plot must be deleted
    clf
    
    xx = cell(n_effects,1);
    % use existing x-variable if available
    if exist('x','var') & numel(x)==size(X,1)
      xx_array = [min(x) max(x)]; 
      for i=1:n_effects
        xx{i} = X(H.SPM{1}.xX.I(:,3)==i);
      end
      x0 = x;
    else
      xx_array = [min(X(X~=0)) max(X(X~=0))]; 
      for i=1:n_effects
        xx{i} = X(H.SPM{1}.xX.I(:,3)==i,i);
      end
      x0 = sum(X,2);
    end
        
    hold on
    for i=1:n_effects
      x2 = xx{i};
      y2 = mean(yy{i},2);
      plot(x2,y2,'.','MarkerSize',Hc.MS,'Color',Hc.col(i,:));

      P = polyfit(x2,y2,1);       
      % plot trend line
      Hc.line{i} = plot(xx_array,polyval(P,xx_array),'Color',Hc.col(i,:),'LineWidth',Hc.LW);
    end

    % for repeated anovas also plot connected lines if defined
    if repeated_anova & Hc.connected
      % coding of subject factor should be hopefully always 2nd column of xX.I
      n_subjects = max(SPM.xX.I(:,2));

      y0 = mean(y,2);
      for i=1:n_subjects

        ind = find(SPM.xX.I(:,2) == i);
        x_tmp = x0(ind);
        y_tmp = y0(ind);
        
        if ~isempty(ind)
          line(x_tmp,y_tmp,'Color',Hc.col(SPM.xX.I(ind(1),3),:));
        end
      end
    end
    
    hold off

    if Hc.legend
      if n_effects == 2
        legend(names(1,:),['Fit ' names(1,:)],names(2,:),['Fit ' names(2,:)])
      else
        legend(names(1,:),['Fit ' names(1,:)])
      end
    end

    TITLE = {['Scatterplot of ' title_name] XYZstr};
  else
    cat_plot_boxplot(yy,vstruct);
    TITLE = {['Boxplot of ' title_name] XYZstr};
    set(gca,'XLim',[0.4 (length(signal_change) + 0.6)],'XTick',1:length(signal_change));
    
    if exist('names','var')
      if size(names,1) == length(signal_change)
        set(gca,'XTickLabel',names,'TickLabelInterpreter','none');
      end
    end  
    
    if repeated_anova & Hc.medianplot
      hold on
  
      plot_data = zeros(n_effects,1);
      count = 1;
      for i=1:n_groups
        for j=1:length(ind_times{i})
          plot_data(count) = median(yy{ind_times{i}(j)});
          count = count + 1;
        end
        Hc.line{i} = plot(ind_times{i},plot_data(ind_times{i}),'Color',Hc.col(i,:),'LineWidth',Hc.LW);
      end
      hold off
    end
  end
  
  title(TITLE,'FontSize',14,'FontWeight','bold')
  ylabel(y_label,'FontSize',12)

    
  set(gca(Hc.h12),'FontSize',Hc.FS);
end

set(gca(Hc.h11),'FontSize',Hc.FS);

%==========================================================================
function set_font_size(obj, event_obj)

global Hc

Hc.FS = str2num(get(obj,'String'));

if isempty(Hc.FS) | numel(Hc.FS)>1
  fprintf('Error: Please enter a single number for defining font size\n');
else
  set(gca(Hc.h11),'FontSize',Hc.FS);
  if Hc.y_found
    set(gca(Hc.h12),'FontSize',Hc.FS);
  end
end

end

%==========================================================================
function set_line_width(obj, event_obj)

global Hc

Hc.LW = str2num(get(obj,'String'));

if isempty(Hc.LW) | numel(Hc.LW)>1
  fprintf('Error: Please enter a single number for defining line width\n');
else
  if Hc.y_found
    if ~isempty(Hc.line{1})
      for i=1:numel(Hc.line)
        set(Hc.line{i},'LineWidth',Hc.LW);        
      end
    end
  end
end

end

%==========================================================================
function set_marker_size(obj, event_obj)

global Hc

Hc.MS = str2num(get(obj,'String'));

if isempty(Hc.MS) | numel(Hc.MS)>1
  fprintf('Error: Please enter a single marker for defining plot symbols\n');
else
  if Hc.y_found
    ic = findobj(Hc.h12, 'Type', 'line');
    for i=1:numel(ic)
      % Marker has more than two points
      if numel(get(ic(i),'XData')) > 2
        set(ic(i),'MarkerSize',Hc.MS);
      end
    end
  end
end

end

%==========================================================================
function adjust_data(obj, event_obj, filename)

global Hc

Hc.adjust = get(obj, 'Value');

end

%==========================================================================
function show_connected(obj, event_obj, filename)

global Hc

Hc.connected = get(obj, 'Value');

end

%==========================================================================
function show_rawdata(obj, event_obj, filename)

global Hc

if Hc.boxplot | Hc.medianplot
  Hc.rawdata = get(obj, 'Value');
end

end

%==========================================================================
function show_medianplot(obj, event_obj, filename)

global Hc

if Hc.rawdata | Hc.boxplot
  Hc.medianplot = get(obj, 'Value');
end

end

%==========================================================================
function show_boxplot(obj, event_obj, filename)

global Hc

if Hc.rawdata | Hc.medianplot
  Hc.boxplot = get(obj, 'Value');
end

end

%==========================================================================
function show_legend(obj, event_obj, filename)

global Hc

Hc.legend = get(obj, 'Value');

end

%==========================================================================
function save_image(obj, event_obj, filename)

global xY Hc

if ~exist('filename', 'var')
    
  filename = xY.string;
    
  [filename, newpth] = uiputfile({ ...
      '*.png' 'PNG files (*.png)'}, 'Save as', filename);
else
    [pth, nam, ext] = fileparts(filename);
    if isempty(pth), pth = cd; end
    if isempty(nam)
        [filename, newpth] = uiputfile({ ...
            '*.png' 'PNG files (*.png)'}, 'Save as', nam);
    else
        filename = fullfile(pth, nam);
        newpth = pth;
    end
end

% remove potential .png
filename = regexprep(filename,'.png','');

try
  % keep background color
  set(Hc.h10, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
  hh = getframe(Hc.h11);
  img = hh.cdata;
  col = colormap;
  saved_file = fullfile(newpth,['estimates_' filename '.png']);
  imwrite(img,col,saved_file);
  fprintf('File %s saved.\n',saved_file);
catch
  fprintf('File %s could not be saved.\n',saved_file);
end

if Hc.y_found
  try
    % keep background color
    set(Hc.h12, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
    hh = getframe(Hc.h12);
    img = hh.cdata;
    col = colormap;
    saved_file = fullfile(newpth,['boxplot_' filename '.png']);
    imwrite(img,col,saved_file);
    fprintf('File %s saved.\n',saved_file);
  catch
    fprintf('File %s could not be saved.\n',saved_file);
  end
end

end