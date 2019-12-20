function spm_LST
% LST Toolbox wrapper to call lst functions
%_______________________________________________________________________
% Paul Schmidt, 2015/08/04

addpath(fullfile(spm('dir'),'toolbox','LST'));
rev = 'ersion 2.0.15';
ps_LST_update(2)

SPMid = spm('FnBanner',mfilename,rev);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','LST');
url = fullfile(spm('Dir'),'toolbox','LST','lst.html');
spm_help('!Disp',url,'',Fgraph,'LST: Lesion segmentation tool for SPM');


fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
'Label',	'LST',...
'Separator',	'on',...
'Tag',		'LST',...
'HandleVisibility','on');
h1  = uimenu(h0,...
'Label',	'Lesion segmentation (LGA)',...
'Separator',	'off',...
'Tag',		'Lesion segmentation (LGA)',...
'HandleVisibility','on');
h11  = uimenu(h1,...
'Label',	'Lesion segmentation',...
'Separator',	'off',...
'Tag',		'Estimate PVE label and segment lesions',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.lga'');',...
'HandleVisibility','on');
h12  = uimenu(h1,...
'Label',	'Determination of optimal threshold',...
'Separator',	'off',...
'Tag',		'Determination of optimal initial threshold',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.doit'');',...
'HandleVisibility','on');
h2  = uimenu(h0,...
'Label',	'Lesion segmentation (LPA)',...
'Separator',	'off',...
'Tag',		'Prediction model based on a binary classifier',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.lpa'');',...
'HandleVisibility','on');
h3  = uimenu(h0,...
'Label',	'Longitudinal pipeline',...
'Separator',	'on',...
'Tag',		'Longitudinal lesion segmentation and quantification of change',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.long'');',...
'HandleVisibility','off');
h4  = uimenu(h0,...
'Label',	'Lesion filling',...
'Separator',	'off',...
'Tag',		'Lesion filling',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.filling'');',...
'HandleVisibility','on');
h5  = uimenu(h0,...
'Label',	'Extract values of interest',...
'Separator',	'off',...
'Tag',		'Total lesion volume',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.tlv'');',...
'HandleVisibility','on');
h6  = uimenu(h0,...
'Label',	'Create binary lesion maps',...
'Separator',	'off',...
'Tag',		'Thresholding of lesion probability maps',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.thresholding'');',...
'HandleVisibility','on');
h7  = uimenu(h0,...
'Label',	'Merge HTML reports',...
'Separator',	'off',...
'Tag',		'Merge existing HTML reports from multiple subjects',...
'CallBack','spm_jobman(''interactive'','''',''spm.tools.LST.merge_reports'');',...
'HandleVisibility','off');
%{
h8  = uimenu(h0,...
'Label',	'Create overlay',...
'Separator',	'off',...
'Tag',		'Thresholding of lesion probability maps',...
'CallBack','try,ps_LST_create_gif;end',...
'HandleVisibility','on');
%}
h8  = uimenu(h0,...
'Label',	'Check for updates',...
'Separator',	'on',...
'Tag',		'Check for updates',...
'CallBack','spm(''alert'',evalc(''ps_LST_update(1)''),''LST Update'');',...
'HandleVisibility','off');
h9  = uimenu(h0,...
'Label',	'Open manual',...
'Separator',	'off',...
'Tag',		'Manual',...
'CallBack','try,open(fullfile(spm(''dir''),''toolbox'',''LST'',''doc'',''LST_documentation.pdf''));end',...
'HandleVisibility','off');
h10  = uimenu(h0,...
'Label',	'LST website',...
'Separator',	'off',...
'Tag',		'LST website',...
'CallBack',['set(gcbf,''Pointer'',''Watch''),',...
			'web(''http://www.statistical-modeling.de/lst.html'',''-browser'');',...
			'set(gcbf,''Pointer'',''Arrow'')'],...
'HandleVisibility','off');