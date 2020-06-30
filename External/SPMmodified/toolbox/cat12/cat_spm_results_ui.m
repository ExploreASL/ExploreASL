function varargout = cat_spm_results_ui(varargin)
% User interface for SPM/PPM results: Display and analysis of regional effects
% FORMAT [hReg,xSPM,SPM] = cat_spm_results_ui('Setup',[xSPM])
%
% hReg   - handle of MIP XYZ registry object
%          (see spm_XYZreg.m for details)
% xSPM   - structure containing specific SPM, distribution & filtering details
%          (see spm_getSPM.m for contents)
% SPM    - SPM structure containing generic parameters
%          (see spm_spm.m for contents)
%
% NB: Results section GUI CallBacks use these data structures by name,
%     which therefore *must* be assigned to the correctly named variables.
%__________________________________________________________________________
%
% The SPM results section is for the interactive exploration and
% characterisation of the results of a statistical analysis.
%
% The user is prompted to select a SPM{T} or SPM{F}, that is thresholded at
% user specified levels. The specification of the contrasts to use and the
% height and size thresholds are described in spm_getSPM.m. The resulting
% SPM is then displayed in the Graphics window as a maximum intensity
% projection, alongside the design matrix and contrasts employed.
%
% The cursors in the MIP can be moved (dragged) to select a particular
% voxel. The three mouse buttons give different drag and drop behaviour:
% Button 1 - point & drop; Button 2 - "dynamic" drag & drop with
% co-ordinate & SPM value updating; Button 3 - "magnetic" drag & drop,
% where the cursor jumps to the nearest suprathreshold voxel in the MIP,
% and shows the value there.
% See spm_mip_ui.m, the MIP GUI handling function for further details.
%
% The design matrix and contrast pictures are "surfable": Click and drag
% over the images to report associated data. Clicking with different
% buttons produces different results. Double-clicking extracts the
% underlying data into the base workspace.
% See spm_DesRep.m for further details.
%
% The current voxel specifies the voxel, suprathreshold cluster, or
% orthogonal planes (planes passing through that voxel) for subsequent
% localised utilities.
%
% A control panel in the Interactive window enables interactive exploration
% of the results.
%
% p-values buttons:
%   (i) volume   - Tabulates p-values and statistics for entire volume.
%                                            - see spm_list.m
%  (ii) cluster  - Tabulates p-values and statistics for nearest cluster.
%                - Note that the cursor will jump to the nearest
%                  suprathreshold voxel, if it is not already at a
%                  location with suprathreshold statistic.
%                                            - see spm_list.m
% (iii)  S.V.C   - Small Volume Correction:
%                  Tabulates p-values corrected for a small specified
%                  volume of interest. (Tabulation by spm_list.m)
%                                            - see spm_VOI.m
%
% Data extraction buttons:
% Eigenvariate/CVA
%                - Extracts the principal eigenvariate for small volumes 
%                  of interest; or CVA of data within a specified volume
%                - Data can be adjusted or not for eigenvariate summaries
%                - If temporal filtering was specified (fMRI), then it is
%                  the filtered data that is returned.
%                - Choose a VOI of radius 0 to extract the (filtered &)
%                  adjusted data for a single voxel. Note that this vector
%                  will be scaled to have a 2-norm of 1. (See spm_regions.m
%                  for further details.)
%                - The plot button also returns fitted and adjusted
%                  (after any filtering) data for the voxel being plotted.)
%                - Note that the cursor will jump to the nearest voxel for
%                  which raw data was saved.
%                                            - see spm_regions.m
%
% Visualisation buttons:
%   (i) plot     - Graphs of adjusted and fitted activity against
%                  various ordinates.
%                - Note that the cursor will jump to the nearest
%                  suprathreshold voxel, if it is not already at a
%                  location with suprathreshold statistic.
%                - Additionally, returns fitted and adjusted data to the
%                  MATLAB base workspace.
%                                               - see spm_graph.m
%  (ii) overlays - Popup menu: Overlays of filtered SPM on a structural image
%     -   slices - Slices of the thresholded statistic image overlaid
%                  on a secondary image chosen by the user. Three
%                  transverse slices are shown, being those at the
%                  level of the cursor in the z-axis and the two
%                  adjacent to it.           - see spm_transverse.m
%     - sections - Orthogonal sections of the thresholded statistic
%                  image overlaid on a secondary image chosen by the user.
%                  The sections are through the cursor position.
%                                            - see spm_sections.m
%     -   render - Render blobs on previously extracted cortical surface
%                                            - see spm_render.m
% (iii) save     - Write out thresholded SPM as image
%                                            - see spm_write_filtered.m
%
% The current cursor location can be set by editing the co-ordinate widgets
% at the bottom of the Interactive window. (Note that many of the results
% section facilities are "linked" and can update co-ordinates. E.g.
% clicking on the co-ordinates in a p-value listing jumps to that location.)
%
% Graphics appear in the bottom half of the Graphics window, additional
% controls and questions appearing in the Interactive window.
%
%                           ----------------
%
% The MIP uses a template outline in MNI space. Consequently for the
% results section to display properly the input images to the statistics
% section should be in MNI space.
%
% Similarly, secondary images should be aligned with the input images used
% for the statistical analysis.
%
%                           ----------------
%
% In addition to setting up the results section, cat_spm_results_ui.m sets
% up the results section GUI and services the CallBacks. FORMAT
% specifications for embedded CallBack functions are given in the main
% body of the code.
%__________________________________________________________________________
% Copyright (C) 1996-2018 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston & Andrew Holmes
% $Id: cat_spm_results_ui.m 1610 2020-04-29 16:49:12Z gaser $
 
 
%==========================================================================
% - FORMAT specifications for embedded CallBack functions
%==========================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take.          )
%
% cat_spm_results_ui sets up and handles the SPM results graphical user
% interface, initialising an XYZ registry (see spm_XYZreg.m) to co-ordinate
% locations between various location controls.
%
%__________________________________________________________________________
%
% FORMAT [hreg,xSPM,SPM] = cat_spm_results_ui('Setup')
% Query SPM and setup GUI. 
%
% FORMAT [hreg,xSPM,SPM] = cat_spm_results_ui('Setup',xSPM)
% Query SPM and setup GUI using a xSPM input structure. This allows to run
% results setup without user interaction. See spm_getSPM for details of
% allowed fields.
%
% FORMAT hReg = cat_spm_results_ui('SetupGUI',M,DIM,xSPM,Finter)
% Setup results GUI in Interactive window
% M       - 4x4 transformation matrix relating voxel to "real" co-ordinates
% DIM     - 3 vector of image X, Y & Z dimensions
% xSPM    - structure containing xSPM. Required fields are:
% .Z      - minimum of n Statistics {filtered on u and k}
% .XYZmm  - location of voxels {mm}
% Finter  - handle (or 'Tag') of Interactive window (default 'Interactive')
% hReg    - handle of XYZ registry object
%
% FORMAT cat_spm_results_ui('DrawButts',hReg,DIM,xSPM,Finter,WS,FS)
% Draw GUI buttons
% hReg    - handle of XYZ registry object
% DIM     - 3 vector of image X, Y & Z dimensions
% xSPM    - structure containing xSPM. Required fields are:
% Finter  - handle of Interactive window
% WS      - WinScale  [Default spm('WinScale') ]
% FS      - FontSizes [Default spm('FontSizes')]
%
% FORMAT hFxyz = cat_spm_results_ui('DrawXYZgui',M,DIM,xSPM,xyz,hReg)
% Setup editable XYZ control widgets at foot of Interactive window
% M      - 4x4 transformation matrix relating voxel to "real" co-ordinates
% DIM    - 3 vector of image X, Y & Z dimensions
% xSPM   - structure containing SPM; Required fields are:
% .Z     - minimum of n Statistics {filtered on u and k}
% .XYZmm - location of voxels {mm}
% xyz    - Initial xyz location {mm}
% hReg   - handle of XYZ registry object
% hFxyz  - handle of XYZ control - the frame containing the edit widgets
%
% FORMAT cat_spm_results_ui('EdWidCB')
% Callback for editable XYZ control widgets
%
% FORMAT cat_spm_results_ui('UpdateSPMval',hFxyz)
% FORMAT cat_spm_results_ui('UpdateSPMval',UD)
% Updates SPM value string in Results GUI (using data from UserData of hFxyz)
% hFxyz - handle of frame enclosing widgets - the Tag object for this control
% UD    - XYZ data structure (UserData of hFxyz).
%
% FORMAT xyz = cat_spm_results_ui('GetCoords',hFxyz)
% Get current co-ordinates from editable XYZ control
% hFxyz - handle of frame enclosing widgets - the Tag object for this control
% xyz   - current co-ordinates {mm}
% NB: When using the results section, should use XYZregistry to get/set location
%
% FORMAT [xyz,d] = cat_spm_results_ui('SetCoords',xyz,hFxyz,hC)
% Set co-ordinates to XYZ widget
% xyz   - (Input) desired co-ordinates {mm}
% hFxyz - handle of XYZ control - the frame containing the edit widgets
% hC    - handle of calling object, if used as a callback. [Default 0]
% xyz   - (Output) Desired co-ordinates are rounded to nearest voxel if hC
%         is not specified, or is zero. Otherwise, caller is assumed to
%         have checked verity of desired xyz co-ordinates. Output xyz returns
%         co-ordinates actually set {mm}.
% d     - Euclidean distance between desired and set co-ordinates.
% NB: When using the results section, should use XYZregistry to get/set location
%
% FORMAT hFxyz = cat_spm_results_ui('FindXYZframe',h)
% Find/check XYZ edit widgets frame handle, 'Tag'ged 'hFxyz'
% h     - handle of frame enclosing widgets, or containing figure [default gcf]
%         If ischar(h), then uses spm_figure('FindWin',h) to locate named figures
% hFxyz - handle of confirmed XYZ editable widgets control
%         Errors if hFxyz is not an XYZ widget control, or a figure containing
%         a unique such control
%
% FORMAT cat_spm_results_ui('PlotUi',hAx)
% GUI for adjusting plot attributes - Sets up controls just above results GUI
% hAx - handle of axes to work with
%
% FORMAT cat_spm_results_ui('PlotUiCB')
% CallBack handler for Plot attribute GUI
%
% FORMAT Fgraph = cat_spm_results_ui('Clear',F,mode)
% Clears results subpane of Graphics window, deleting all but semi-permanent
% results section stuff
% F      - handle of Graphics window [Default spm_figure('FindWin','Graphics')]
% mode   - 1 [default] - clear results subpane
%        - 0           - clear results subpane and hide results stuff
%        - 2           - clear, but respect 'NextPlot' 'add' axes
%                        (which is set by `hold on`)
% Fgraph - handle of Graphics window
%
% FORMAT hMP = cat_spm_results_ui('LaunchMP',M,DIM,hReg,hBmp)
% Prototype callback handler for integrating MultiPlanar toolbox
%
% FORMAT cat_spm_results_ui('Delete',h)
% deletes HandleGraphics objects, but only if they're valid, thus avoiding
% warning statements from MATLAB.
%
% modified version of
% spm_results_ui.m r7388
%__________________________________________________________________________
 
SVNid = '$Rev: 1610 $'; 

%-Condition arguments
%--------------------------------------------------------------------------
if nargin == 0, Action='Setup'; else Action=varargin{1}; end
useCAT = 2; % 0-like SPM, 1-surface handling, 2-cat_surf_renderer  
 
global result_ui_varargout use_tfce mesh_detect

% prevent that TFCE is called if not yet installed
if ~exist(fullfile(spm('dir'),'toolbox','TFCE'),'dir')
  use_tfce = 0;
end
 
%==========================================================================
switch lower(Action), case 'setup'                         %-Set up results
%==========================================================================
 
    %-Initialise
    %----------------------------------------------------------------------
    spm('FnBanner',mfilename,SVNid);
    try
      dcm = datacursormode(spm_figure('FindWin','Graphics')); 
      set(dcm,'Enable','off','UpdateFcn',[]); 
      
      spm_figure('Clear',spm_figure('FindWin','Graphics')); 
    end
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
    spm_clf('Satellite');
 
    %-Get thresholded xSPM data and parameters of design
    %======================================================================
    if nargin > 1
        if use_tfce
            [SPM,xSPM] = tfce_getSPM(varargin{2});
        else
            [SPM,xSPM] = spm_getSPM(varargin{2});
        end
    else
        [spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
        swd = spm_file(spmmatfile,'fpath');
        load(fullfile(swd,'SPM.mat'));
        [Ic,xCon] = spm_conman(SPM,'T&F',Inf,'    Select contrast(s)...');
        SPM.Ic = Ic; SPM.xCon = xCon;
        SPM.swd = swd;
        
        % check for existing TFCE results for this contrast
        if numel(Ic)==1 & exist(fullfile(swd,sprintf('%s_log_p_%04d.nii',xCon(Ic).STAT,Ic))) || ...
                          exist(fullfile(swd,sprintf('%s_log_p_%04d.gii',xCon(Ic).STAT,Ic)))
						stat_str = {'TFCE',xCon(Ic).STAT};
						statType = spm_input('Type of statistic',1,'m',...
								sprintf('TFCE (non-parametric)|%s (non-parametric)|%s (SPM parametric)',...
								xCon(Ic).STAT,xCon(Ic).STAT),[],1);
						if statType < 3
								use_tfce = 1;
								SPM.statType = stat_str{statType};
								[SPM,xSPM] = tfce_getSPM(SPM);
								xSPM.statType = stat_str{statType};
						else
								use_tfce = 0;
								[SPM,xSPM] = spm_getSPM(SPM);
								xSPM.statType = xCon(Ic).STAT;
						end
        else
						use_tfce = 0;
						[SPM,xSPM] = spm_getSPM(SPM);
						xSPM.statType = xCon(Ic).STAT;
        end
    end
 
    if isempty(xSPM) 
        varargout = {[],[],[]};
        return;
    end
 
    if spm_mesh_detect(xSPM.Vspm)
      mesh_detect = 1;
    else
      mesh_detect = 0;
    end

    %-Ensure pwd = swd so that relative filenames are valid
    %----------------------------------------------------------------------
    cd(SPM.swd)
    
    %-Get space information
    %======================================================================
    M   = SPM.xVol.M;
    DIM = SPM.xVol.DIM;
    
    %-Space units
    %----------------------------------------------------------------------
    try
        try
            units = SPM.xVol.units;
        catch
            units = xSPM.units;
        end
    catch
        try
            Modality = spm('CheckModality');
        catch
            Modality = {'PET','FMRI','EEG'};
            selected = spm_input('Modality: ','+1','m',Modality);
            Modality = Modality{selected};
            spm('ChMod',Modality);
        end
        if strcmp(Modality,'EEG')
            datatype = {...
                'Volumetric (2D/3D)',...
                'Scalp-Time',...
                'Scalp-Frequency',...
                'Time-Frequency',...
                'Frequency-Frequency'};
            selected = spm_input('Data Type: ','+1','m',datatype);
            datatype = datatype{selected};
        else
            datatype = 'Volumetric (2D/3D)';
        end
        
        switch datatype
            case 'Volumetric (2D/3D)'
                units    = {'mm' 'mm' 'mm'};
            case 'Scalp-Time'
                units    = {'mm' 'mm' 'ms'};
            case 'Scalp-Frequency'
                units    = {'mm' 'mm' 'Hz'};
            case 'Time-Frequency'
                units    = {'Hz' 'ms' ''};
            case 'Frequency-Frequency'
                units    = {'Hz' 'Hz' ''};
            otherwise
                error('Unknown data type.');
        end
    end

    
    % CAT.begin    
    % -------------------------------------------------------------------------
    % ToDo:
    % * add batch call with render settings > cat_conf_stools
    % * add atlas data coursor > cat_surf_render
    % * add colorbar > cat_surf_render
    % * fix contrast box boundaries
    % * full atlas integration > spm_atlas, spm_XYZreg, spm_list (lot of work) 
    % -------------------------------------------------------------------------
    if useCAT
        [pp,ff,ee] = spm_fileparts(xSPM.Vspm.fname);
        if strcmp(ee,'.gii')
          datatype = 'Surface (3D)'; 
         % units    = {'mm' 'mm'}; 
        end
        
        % change surface ?
        % - first we have to use the Shooting Surface for the statistic to
        %   obtain more meaningful MNI coordinates 
        % - here we have to change to the used surface (FSaverage)
        if spm_mesh_detect(xSPM.Vspm) & exist('spm_cat12')
          FSavg = '.freesurfer.gii'; 
          GSavg = '.Template_T1_IXI555_MNI152_GS.gii';
          if ischar(SPM.xVol.G)
            SPM.xVol.G = strrep(SPM.xVol.G,GSavg,FSavg);
          end
        end
        
        % change coordinates to spherial ??
        %XYZmmFS = xSPM.
        
    end
    % -------------------------------------------------------------------------
    % CAT.end

    if spm_mesh_detect(xSPM.Vspm)
        DIM(3) = Inf; % force 3D coordinates
    elseif DIM(3) == 1
        units{3} = '';
        if DIM(2) == 1
            units{2} = '';
        end
    end
    xSPM.units      = units;
    SPM.xVol.units  = units;
    
     
    %-Setup Results User Interface; Display MIP, design matrix & parameters
    %======================================================================
 
    %-Setup results GUI
    %----------------------------------------------------------------------
    spm_clf(Finter);
    spm('FigName',['SPM{',xSPM.STAT,'}: Results'],Finter,CmdLine);
    hReg      = cat_spm_results_ui('SetupGUI',M,DIM,xSPM,Finter);
 
    %-Setup design interrogation menu
    %----------------------------------------------------------------------
    hDesRepUI = spm_DesRep('DesRepUI',SPM);
 
    %-Setup contrast menu
    %----------------------------------------------------------------------
    hConUI    = cat_spm_results_ui('SetupConMenu',xSPM,SPM,Finter);
    
    %-Atlas menu
    %----------------------------------------------------------------------
    if ~spm_mesh_detect(xSPM.Vspm)
        hAtlasUI = cat_spm_results_ui('SetupAtlasMenu',Finter);
    end
    
    %-Setup Maximum intensity projection (MIP) & register
    %----------------------------------------------------------------------
    FS     = spm('FontSizes');
    hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
    if spm_mesh_detect(xSPM.Vspm)
       % block to call cat_surf_render rather than spm_mesh_render
        if useCAT>1 & exist('cat_surf_render')
            hMax = cat_surf_render('Disp',SPM.xVol.G,'Parent',hMIPax,'Results',1);
            tmp  = zeros(1,prod(xSPM.DIM));
            tmp(xSPM.XYZ(1,:)) = xSPM.Z;
            hMax = cat_surf_render('Overlay',hMax,tmp);
            hMax = cat_surf_render('ColourMap',hMax,jet);
            hMax = cat_surf_render('Register',hMax,hReg);
        else
            hMax = spm_mesh_render('Disp',SPM.xVol.G,'Parent',hMIPax);
            tmp  = zeros(1,prod(xSPM.DIM));
            tmp(xSPM.XYZ(1,:)) = xSPM.Z;
            hMax = spm_mesh_render('Overlay',hMax,tmp);
            hMax = spm_mesh_render('Register',hMax,hReg);
            hMax = spm_mesh_render('ColourMap',hMax,jet);
            hMax = spm_mesh_render('View',hMax,'top');
        end
    elseif isequal(units(2:3),{'' ''})
        set(hMIPax, 'Position',[0.05 0.65 0.55 0.25]);
        [allS,allXYZmm] = spm_read_vols(xSPM.Vspm);
        plot(hMIPax,allXYZmm(1,:),allS,'Color',[0.6 0.6 0.6]);
        set(hMIPax,'NextPlot','add');
        MIP = NaN(1,xSPM.DIM(1));
        MIP(xSPM.XYZ(1,:)) = xSPM.Z;
        XYZmm = xSPM.M(1,:)*[1:xSPM.DIM(1);zeros(2,xSPM.DIM(1));ones(1,xSPM.DIM(1))];
        plot(hMIPax,XYZmm,MIP,'b-+','LineWidth',2);
        plot(hMIPax,[XYZmm(1) XYZmm(end)],[xSPM.u xSPM.u],'r');
        clim = get(hMIPax,'YLim');
        axis(hMIPax,[sort([XYZmm(1) XYZmm(end)]) 0 clim(2)]);
        %set(hMIPax,'XTick',[],'YTick',[]);
    else
        hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,M,DIM,hMIPax,units);
        spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
    end
    
    if xSPM.STAT == 'P'
        str = xSPM.STATstr;
    else
        str = ['SPM\{',xSPM.STATstr,'\}'];
    end
    text(240,260,str,...
        'Interpreter','TeX',...
        'FontSize',FS(14),'Fontweight','Bold',...
        'Parent',hMIPax)
 
 
    %-Print comparison title
    %----------------------------------------------------------------------
    hTitAx = axes('Parent',Fgraph,...
        'Position',[0.02 0.96 0.96 0.04],...
        'Visible','off');
 
    text(0.5,0.5,xSPM.title,'Parent',hTitAx,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'FontWeight','Bold','FontSize',FS(14))
 
 
    %-Print SPMresults: Results directory & thresholding info
    %----------------------------------------------------------------------
    hResAx = axes('Parent',Fgraph,...
        'Position',[0.05 0.55 0.45 0.05],...
        'DefaultTextVerticalAlignment','baseline',...
        'DefaultTextFontSize',FS(9),...
        'DefaultTextColor',[1,1,1]*.7,...
        'Units','points',...
        'Visible','off');
    AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
    h     = text(0,24,'SPMresults:','Parent',hResAx,...
        'FontWeight','Bold','FontSize',FS(14));
    text(get(h,'Extent')*[0;0;1;0],24,spm_file(SPM.swd,'short30'),'Parent',hResAx)

    try
        thresDesc = xSPM.thresDesc;
        if use_tfce
            text(0,12,sprintf('Height threshold %s',thresDesc),'Parent',hResAx)
        else
            if strcmp(xSPM.STAT,'P')
                text(0,12,sprintf('Height threshold %s',thresDesc),'Parent',hResAx)
            else
                text(0,12,sprintf('Height threshold %c = %0.6f  {%s}',xSPM.STAT,xSPM.u,thresDesc),'Parent',hResAx)
            end
        end
    catch
        text(0,12,sprintf('Height threshold %c = %0.6f',xSPM.STAT,xSPM.u),'Parent',hResAx)
    end
    if spm_mesh_detect(xSPM.Vspm), str = 'vertices'; else str = 'voxels'; end
    text(0,00,sprintf('Extent threshold k = %0.0f %s',xSPM.k,str), 'Parent',hResAx)
 
 
    %-Plot design matrix
    %----------------------------------------------------------------------
    hDesMtx   = axes('Parent',Fgraph,'Position',[0.65 0.55 0.25 0.25]);
    hDesMtxIm = image((SPM.xX.nKX + 1)*32,'Parent',hDesMtx);
    xlabel(hDesMtx,'Design matrix','FontSize',FS(10))
    set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')',...
        'UserData',struct(...
        'X',        SPM.xX.xKXs.X,...
        'fnames',   {reshape({SPM.xY.VY.fname},size(SPM.xY.VY))},...
        'Xnames',   {SPM.xX.name}))
 
    %-Plot contrasts
    %----------------------------------------------------------------------
    nPar   = size(SPM.xX.X,2);
    xx     = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
    nCon   = length(xSPM.Ic);
    xCon   = SPM.xCon;
    if nCon
        dy     = 0.15/max(nCon,2);
        hConAx = axes('Parent',Fgraph, 'Position',[0.65 (0.80 + dy*.1) 0.25 dy*(nCon-.1)],...
            'Tag','ConGrphAx','Visible','off');
        if use_tfce & xSPM.invResult
          str = 'inverse contrast';
        else
          str = 'contrast';
        end
        if nCon > 1, str = [str 's']; end
        title(hConAx,str)
        htxt   = get(hConAx,'title');
        set(htxt,'FontSize',FS(10),'FontWeight','normal','Visible','on','HandleVisibility','on')
    end
 
    for ii = nCon:-1:1
        hCon = axes('Parent',Fgraph, 'Position',[0.65 (0.80 + dy*(nCon - ii +.1)) 0.25 dy*.9]);
        if xCon(xSPM.Ic(ii)).STAT == 'T' && size(xCon(xSPM.Ic(ii)).c,2) == 1
 
            %-Single vector contrast for SPM{t} - bar
            %--------------------------------------------------------------
            yy = [zeros(1,nPar);repmat(xCon(xSPM.Ic(ii)).c',2,1);zeros(1,nPar)];
            if use_tfce & xSPM.invResult, yy = -yy; end
            h  = patch(xx,yy,[1,1,1]*.5,'Parent',hCon);
            set(hCon,'Tag','ConGrphAx',...
                'Box','off','TickDir','out',...
                'XTick',spm_DesRep('ScanTick',nPar,10) - 0.5,'XTickLabel','',...
                'XLim', [0,nPar],...
                'YTick',[-1,0,+1],'YTickLabel','',...
                'YLim',[min(xCon(xSPM.Ic(ii)).c),max(xCon(xSPM.Ic(ii)).c)] +...
                [-1 +1] * max(abs(xCon(xSPM.Ic(ii)).c))/10  )
 
        else
 
            %-F-contrast - image
            %--------------------------------------------------------------
            h = image((xCon(xSPM.Ic(ii)).c'/max(abs(xCon(xSPM.Ic(ii)).c(:)))+1)*32,...
                'Parent',hCon);
            set(hCon,'Tag','ConGrphAx',...
                'Box','on','TickDir','out',...
                'XTick',spm_DesRep('ScanTick',nPar,10),'XTickLabel','',...
                'XLim', [0,nPar]+0.5,...
                'YTick',[0:size(SPM.xCon(xSPM.Ic(ii)).c,2)]+0.5,...
                'YTickLabel','',...
                'YLim', [0,size(xCon(xSPM.Ic(ii)).c,2)]+0.5 )
 
        end
        ylabel(hCon,num2str(xSPM.Ic(ii)),'FontSize',FS(10),'FontWeight','normal')
        set(h,'ButtonDownFcn','spm_DesRep(''SurfCon_CB'')',...
            'UserData', struct( 'i',    xSPM.Ic(ii),...
                                'h',    htxt,...
                                'xCon', xCon(xSPM.Ic(ii))))
    end
 
 
    %-Store handles of results section Graphics window objects
    %----------------------------------------------------------------------
    H  = get(Fgraph,'Children');
    H  = findobj(H,'flat','HandleVisibility','on');
    H  = findobj(H);
    Hv = get(H,'Visible');
    set(hResAx,'Tag','PermRes','UserData',struct('H',H,'Hv',{Hv}))
 

    TabDat = call_list('List',xSPM,hReg);

    %-Finished results setup
    %----------------------------------------------------------------------
    varargout = {hReg,xSPM,SPM,TabDat};
    spm('Pointer','Arrow')
 
 
    
  
% -------------------------------------------------------------------------
% This block is required for postprocessing the SPM result figure. 
% It include fixes to avoid rotation of non 3D elements and some other tiny 
% changes to the contrast & design matrix.

    if isfield(xSPM,'G')
        % create SPM result table and fix elements to avoid rotation of tables
        spm_list_cleanup; 

        % corrections for top elements
        hRes.Fgraph       = spm_figure('GetWin','Graphics');
        hRes.FgraphC      = get( hRes.Fgraph ,'children'); 
        hRes.FgraphAx     = findobj( hRes.FgraphC,'Type','Axes');
        hRes.FgraphAxPos  = cell2mat(get( hRes.FgraphAx , 'Position'));
        hRes.Ftext        = findobj(hRes.Fgraph,'Type','Text');
        % fine red lines of the SPM result table
        hRes.Fline        = findobj(hRes.Fgraph,'Type','Line','Tag','');% ,'UIcontextMenu',[]);
        hRes.FlineAx      = get(hRes.Fline,'parent');
    
        % find the SPM result texts to fix them against rotation
        hres.Ftext      = findobj(hRes.Fgraph,'Type','Text','Color',[0.7 0.7 0.7]); 
        hRes.Ftext3dres = get(  hres.Ftext ,'parent'); 
       % for axi=1:numel(hRes.Ftext),      set( hRes.Ftext(axi),'Color',[0.2 0.2 0.2]); end
        for axi = 1:numel(hRes.Ftext3dres), set( hRes.Ftext3dres{axi},'HitTest','off'); end
        for axi = 1:numel(hRes.FlineAx ),   set( hRes.FlineAx{axi},'visible','off'); end
        
        % make nice contrast box that is a bit larger than the original boxes
        hRes.Fcons        = hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.65 & hRes.FgraphAxPos(:,2) > 0.6 ) ;  
        for axi = 1:numel( hRes.Fcons ), l = get( hRes.Fcons(axi) , 'ylim'); set(hRes.Fcons(axi) , 'box','on','ylim', round(l) + [-0.015 0.015]); end
    
        % remove non integer values
        hRes.Fdesm        = hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.65 & hRes.FgraphAxPos(:,2) < 0.6 ) ;  
        xt = get(hRes.Fdesm,'xtick'); xt(round(xt)~=xt) = []; set(hRes.Fdesm,'xtick',xt); 
    
        %
        hRes.Fval         = hRes.FgraphAx( hRes.FgraphAxPos(:,1) > 10); 
        hRes.Fsurf        = hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.05);
        hRes.Flabels      = [ hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.65); hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.02)]; 
        for axi = 1:numel( hRes.Flabels ), set(hRes.Flabels(axi),'HitTest','off'); end
    
        if nargout==0, 
            fprintf( ...
             ['\n' ...
              '========================================================================\n' ...
              '  You have to call cat_spm_results_ui with all output parameters: \n' ...
              '    [hReg,xSPM,SPM] = cat_spm_results_ui; \n\n' ...
              '  Otherwise, the menu and tables will not work properly! Call now: \n    %s\n'...
              '========================================================================\n\n'], ...
              spm_file('[hReg,xSPM,SPM] = cat_spm_results_ui(''Output'');',...
                'link','[hReg,xSPM,SPM] = cat_spm_results_ui(''Output'');'));  
          
            result_ui_varargout = varargout; 
            clear varargout; 
        end
    end
% -------------------------------------------------------------------------
    
    %======================================================================
    case 'output'                            %-Set up results section GUI
    %======================================================================
        fprintf('Updated result_ui output.\n'); 
        varargout = result_ui_varargout; 
 
    %======================================================================
    case 'setupgui'                            %-Set up results section GUI
    %======================================================================
        % hReg = cat_spm_results_ui('SetupGUI',M,DIM,xSPM,Finter)
        if nargin < 5, Finter='Interactive'; else Finter = varargin{5}; end
        if nargin < 4, error('Insufficient arguments'), end
        M      = varargin{2};
        DIM    = varargin{3};
        xSPM   = varargin{4};
        Finter = spm_figure('GetWin',Finter);
        WS     = spm('WinScale');
        FS     = spm('FontSizes');
 
        %-Create frame for Results GUI objects
        %------------------------------------------------------------------
        hPan = uipanel('Parent',Finter,'Title','','Units','Pixels',...
                'Position',[001 001 400 190].*WS,...
                'BorderType','Line', 'HighlightColor',[0 0 0],...
                'BackgroundColor',spm('Colour'));
        hReg = uipanel('Parent',hPan,'Title','','Units','Pixels',...
                'BorderType','Etchedin', ...
                'Position',[005 005 390 180].*WS,...
                'BackgroundColor',[179 179 179]/255);
 
        %-Initialise registry in hReg frame object
        %------------------------------------------------------------------
        [hReg,xyz] = spm_XYZreg('InitReg',hReg,M,DIM,[0;0;0]);
 
        %-Setup editable XYZ widgets & cross register with registry
        %------------------------------------------------------------------
        hFxyz      = cat_spm_results_ui('DrawXYZgui',M,DIM,varargin{4},xyz,hReg);
        spm_XYZreg('XReg',hReg,hFxyz,'cat_spm_results_ui');
 
        %-Set up buttons for results functions
        %------------------------------------------------------------------
        cat_spm_results_ui('DrawButts',hReg,DIM,xSPM,Finter,WS,FS);

        if spm_check_version('matlab','7.11') ~= 0
            drawnow; % required to force "ratio locking"
            set(findobj(hPan),'Units','Normalized','FontUnits','Normalized');
        end
        
        varargout  = {hReg};
 
    %======================================================================
    case 'drawbutts'   %-Draw results section buttons in Interactive window
    %======================================================================
        % cat_spm_results_ui('DrawButts',hReg,DIM,xSPM,Finter,WS,FS)
        %
        if nargin < 4, error('Insufficient arguments'), end
        hReg = varargin{2};
        DIM  = varargin{3};
        xSPM = varargin{4};
        if nargin < 5, Finter = spm_figure('FindWin','Interactive');
        else Finter = varargin{5}; end
        if nargin < 6, WS = spm('WinScale');  else  WS = varargin{6}; end
        if nargin < 7, FS = spm('FontSizes'); else  FS = varargin{7}; end

        %-p-values
        %------------------------------------------------------------------
        hPan = uipanel('Parent',hReg,'Title','p-values','Units','Pixels',...
            'Position',[005 085 110 092].*WS,...
            'BorderType','Beveledout', ...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
				if use_tfce
					callback = 'TabDat = tfce_list(''List'',xSPM,hReg);cat_spm_results_ui(''spm_list_cleanup'',hReg);';
				else
					callback = 'TabDat = spm_list(''List'',xSPM,hReg);cat_spm_results_ui(''spm_list_cleanup'',hReg);';
				end
        uicontrol('Parent',hPan,'Style','PushButton','String','whole brain',...
            'Units','Pixels',...
            'FontSize',FS(10),...
            'ToolTipString',...
            'Tabulate summary of local maxima, p-values & statistics',...
            'Callback',callback,...
            'Interruptible','on','Enable','on',... 
            'Position',[005 055 100 020].*WS);
				if use_tfce
					callback = 'TabDat = tfce_list(''ListCluster'',xSPM,hReg);cat_spm_results_ui(''spm_list_cleanup'',hReg);';
				else
					callback = 'TabDat = spm_list(''ListCluster'',xSPM,hReg);cat_spm_results_ui(''spm_list_cleanup'',hReg);';
				end
        uicontrol('Parent',hPan,'Style','PushButton','String','current cluster',...
            'Units','Pixels',...
            'FontSize',FS(10),...
            'ToolTipString',...
            'Tabulate p-values & statistics for local maxima of nearest cluster',...
            'Callback',callback,...
            'Interruptible','off','Enable','on',...
            'Position',[005 030 100 020].*WS); 
 

        %-SPM area - used for Volume of Interest analyses
        %------------------------------------------------------------------
        if spm_mesh_detect(xSPM.Vspm) | use_tfce
            Enable = 'off';
        else
            Enable = 'on';
        end
        
        uicontrol('Parent',hPan,'Style','PushButton','String','small volume',...
            'Units','Pixels',...
            'FontSize',FS(10),...
            'ToolTipString',['Small Volume Correction - corrected p-values ',...
            'for a small search region'],...
            'Callback','TabDat = spm_VOI(SPM,xSPM,hReg);',...
            'Interruptible','off','Enable',Enable,...  
            'Position',[005 005 100 020].*WS);

        hPan = uipanel('Parent',hReg,'Title','Multivariate','Units','Pixels',...
            'Position',[120 085 150 092].*WS,...
            'BorderType','Beveledout', ...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
        uicontrol('Parent',hPan,'Style','PushButton','String','eigenvariate',...
            'Position',[005 055 069 020].*WS,...
            'ToolTipString',...
            'Responses (principal eigenvariate) in volume of interest',...
            'Callback','[Y,xY] = spm_regions(xSPM,SPM,hReg)',...
            'Interruptible','on','Enable',Enable,...
            'FontSize',FS(10));
        uicontrol('Parent',hPan,'Style','PushButton','String','CVA',...
            'Position',[076 055 069 020].*WS,...
            'ToolTipString',...
            'Canonical variates analysis for the current contrast and VOI',...
            'Callback','CVA = spm_cva_ui('''',xSPM,SPM)',...
            'Interruptible','off','Enable',Enable,...
            'FontSize',FS(10));            

  
        %-Visualisation
        %------------------------------------------------------------------
        hPan = uipanel('Parent',hReg,'Title','Display','Units','Pixels',...
            'Position',[275 085 110 092].*WS,...
            'BorderType','Beveledout',...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
        uicontrol('Parent',hPan,'Style','PushButton','String','plot',...
            'FontSize',FS(10),...
            'ToolTipString','plot data & contrasts at current voxel',...
            'Callback','[Y,y,beta,Bcov] = spm_graph_ui(xSPM,SPM,hReg);',...
            'Interruptible','on','Enable','on',...
            'Position',[005 055 100 020].*WS,...
            'Tag','plotButton');
 
        str0 = {'overlays...',...
            'slices', ...  
            'sections', ...
            'CAT T1 IXI555 GS', ...
            'montage',... 
            'render',...  
            'previous sections',...
            'previous render'};
        tstr0 = { 'overlay filtered SPM on another image: ',...
            '3 slices / ',...
            'slice overlay /', ...
            'T1 avg slice overlay /', ...
            'ortho sections / ', ...
            'render /', ...
            'previous ortho sections /', ...
            'previous surface rendering'};

        tmp0 = {'spm_transverse(''set'',xSPM,hReg)',...
            'spm_sections(xSPM,hReg);',...
            ['spm_sections(xSPM,hReg,fullfile(spm(''dir''),''toolbox'',''cat12'','...
             '''templates_volumes'',''Template_T1_IXI555_MNI152_GS.nii''));',...
             'cat_spm_results_ui(''spm_list_cleanup'');',...
            ],... 
             {@myslover},...
            ['spm_render(   struct( ''XYZ'',    xSPM.XYZ,',...
            '''t'',     xSPM.Z'',',...
            '''mat'',   xSPM.M,',...
            '''dim'',   xSPM.DIM))'],...
            ['global prevsect;','spm_sections(xSPM,hReg,prevsect);'],...
            ['global prevrend;','if ~isstruct(prevrend)',...
            'prevrend = struct(''rendfile'','''',''brt'',[],''col'',[]); end;',...            
            'spm_render(    struct( ''XYZ'',    xSPM.XYZ,',...
            '''t'',     xSPM.Z'',',...
            '''mat'',   xSPM.M,',...
            '''dim'',   xSPM.DIM),prevrend.brt,prevrend.rendfile)']};
            
        if spm_mesh_detect(xSPM.Vspm)
            ind = [1 3 7];
        else
            ind = [1:3 5:numel(str0)];
        end
        
        % add average image only if cat12 exists
        if exist('spm_cat12')
          ind = sort([ind 4]);
        end

        % select entries
        str  = str0(ind);
        tstr = tstr0(ind);
        tmp  = tmp0(ind(2:end)-1);
        
        uicontrol('Parent',hPan,'Style','popupmenu','String',str,...
            'FontSize',FS(10),...
            'ToolTipString',cat(2,tstr{:}),...
            'Callback','spm(''PopUpCB'',gcbo)',...
            'UserData',tmp,...
            'Interruptible','on','Enable','on',...
            'Position',[005 030 100 020].*WS);
 
        str = {'save...',...
               'thresholded SPM',...
               'all clusters (binary)',...
               'all clusters (n-ary)',...
               'current cluster'};
        tmp = {{@mysavespm, 'thresh' },...
               {@mysavespm, 'binary' },...
               {@mysavespm, 'n-ary'  },...
               {@mysavespm, 'current'}};
        
        uicontrol('Parent',hPan,'Style','popupmenu','String',str,...
            'FontSize',FS(10),...
            'ToolTipString','Save as image',...
            'Callback','spm(''PopUpCB'',gcbo)',...
            'UserData',tmp,...
            'Interruptible','on','Enable','on',...
            'Position',[005 005 100 020].*WS);
 
        %-ResultsUI controls
        %------------------------------------------------------------------
        uicontrol('Parent',hReg,'Style','PushButton','String','clear',...
            'ToolTipString','Clear results subpane',...
              'FontSize',FS(9),'ForegroundColor','r',...
            'Callback',['cat_spm_results_ui(''Clear''); ',...
              'spm_input(''!DeleteInputObj''),',...
              'spm_clf(''Satellite'')'],...
            'Interruptible','on','Enable','on',...
            'DeleteFcn','spm_clf(''Graphics'')',...
            'Position',[280 050 048 020].*WS);
 
        uicontrol('Parent',hReg,'Style','PushButton','String','exit',...
            'ToolTipString','Exit the results section',...
            'FontSize',FS(9),'ForegroundColor','r',...
            'Callback','cat_spm_results_ui(''close'')',...
            'Interruptible','on','Enable','on',...
            'Position',[332 050 048 020].*WS);
 
 
    %======================================================================
    case 'setupconmenu'                               %-Setup Contrast Menu
    %======================================================================
    % cat_spm_results_ui('SetupConMenu',xSPM,SPM,Finter)
    if nargin < 4, Finter = 'Interactive'; else Finter = varargin{4}; end
    if nargin < 3, error('Insufficient arguments'), end
    xSPM = varargin{2};
    SPM  = varargin{3};
    Finter = spm_figure('GetWin',Finter);
    hC   = uimenu(Finter,'Label','Contrasts', 'Tag','ContrastsUI');
    
    if use_tfce
        hC1  = uimenu(hC,'Label','Change Contrast');
        
        for i=1:numel(SPM.xCon)
            % check whether TFCE results were found
            if exist(fullfile(SPM.swd,sprintf('%s_log_p_%04d.nii',xSPM.STAT,i))) || ...
               exist(fullfile(SPM.swd,sprintf('%s_log_p_%04d.gii',xSPM.STAT,i)))
                xSPM2 = xSPM;
                xSPM2.invResult = 0;
                hC2 = uimenu(hC1,'Label',[SPM.xCon(i).STAT, ': ', SPM.xCon(i).name], ...
                    'UserData',struct('Ic',i),...
                    'Callback',{@mychgcon,xSPM2});
                if any(xSPM.Ic == i) & ~xSPM.invResult
                    set(hC2,'ForegroundColor',[0 0 1],'Checked','on');
                end
                xSPM2 = xSPM;
                xSPM2.invResult = 1;
                hC3 = uimenu(hC1,'Label',['(inverse contrast) ',SPM.xCon(i).STAT, ': ', SPM.xCon(i).name], ...
                    'UserData',struct('Ic',i),...
                    'Callback',{@mychgcon,xSPM2});
                if any(xSPM.Ic == i) & xSPM.invResult
                    set(hC3,'ForegroundColor',[0 0 1],'Checked','on');
                end
            end
        end
    else
        hC1  = uimenu(hC,'Label','New Contrast...',...
            'UserData',struct('Ic',0),...
            'Callback',{@mychgcon,xSPM});
        hC1 = uimenu(hC,'Label','Change Contrast');
        for i=1:numel(SPM.xCon)
            hC2 = uimenu(hC1,'Label',[SPM.xCon(i).STAT, ': ', SPM.xCon(i).name], ...
                'UserData',struct('Ic',i),...
                'Callback',{@mychgcon,xSPM});
            if any(xSPM.Ic == i)
                set(hC2,'ForegroundColor',[0 0 1],'Checked','on');
            end
        end
        hC1 = uimenu(hC,'Label','Previous Contrast',...
            'Accelerator','P',...
            'UserData',struct('Ic',xSPM.Ic-1),...
            'Callback',{@mychgcon,xSPM});
        if xSPM.Ic-1<1, set(hC1,'Enable','off'); end
        hC1 = uimenu(hC,'Label','Next Contrast',...
            'Accelerator','N',...
            'UserData',struct('Ic',xSPM.Ic+1),...
            'Callback',{@mychgcon,xSPM});
        if xSPM.Ic+1>numel(SPM.xCon), set(hC1,'Enable','off'); end
    end
    
    hC1 = uimenu(hC,'Label','Significance level','Separator','on');
    xSPMtmp = xSPM; xSPMtmp.thresDesc = '';
    uimenu(hC1,'Label','Change...','UserData',struct('Ic',xSPM.Ic),...
        'Callback',{@mychgcon,xSPMtmp}); 
    
    if strcmp(xSPM.STAT,'P')
        xSPMtmp = xSPM; xSPMtmp.thresDesc = 'LogBF';
        uimenu(hC1,'Label','Set LogBF','UserData',struct('Ic',xSPM.Ic),...
            'Callback',{@mychgcon,xSPMtmp});
    else
        xSPMtmp = xSPM; xSPMtmp.thresDesc = 'p<0.05 (FWE)';
        uimenu(hC1,'Label','Set to 0.05 (FWE)','UserData',struct('Ic',xSPM.Ic),...
            'Callback',{@mychgcon,xSPMtmp});
        xSPMtmp = xSPM; xSPMtmp.thresDesc = 'p<0.001 (unc.)';
        uimenu(hC1,'Label','Set to 0.001 (unc.)','UserData',struct('Ic',xSPM.Ic),...
            'Callback',{@mychgcon,xSPMtmp});
    end
    
    uimenu(hC1,'Label',[xSPM.thresDesc ', k=' num2str(xSPM.k)],...
        'Enable','off','Separator','on');
    
    hC1 = uimenu(hC,'Label','Default cluster listings in result table',...
        'Separator','on',...
            'Callback',{@mycluster,xSPM,8,3});

    hC1 = uimenu(hC,'Label','Detailed cluster listings in result table',...
        'Separator','off',...
            'Callback',{@mycluster,xSPM,4,16});

    hC1 = uimenu(hC,'Label','Multiple display...',...
        'Separator','on',...
        'Callback',{@mycheckres,xSPM});
        
    varargout = {hC};
    
        
    %======================================================================
    case 'setupatlasmenu'                                %-Setup Atlas Menu
    %======================================================================
        % cat_spm_results_ui('SetupAtlasMenu',Finter)
    
        Finter = varargin{2};
        
        %hC  = uicontextmenu;
        hC   = uimenu(Finter,'Label','Atlas', 'Tag','AtlasUI');
        
        hC1  = uimenu(hC,'Label','Label using');
        
        list = spm_atlas('List','installed');
        for i=1:numel(list)
            uimenu(hC1,'Label',list(i).name,...
                'Callback',sprintf('call_list(''label'',''%s''); cat_spm_results_ui(''spm_list_cleanup'');',list(i).name));
        end
        if isempty(list), set(hC1,'Enable','off'); end
        
        %hC2  = uimenu(hC,'Label','Download Atlas...',...
        %    'Separator','on',...
        %    'Callback','spm_atlas(''install'');');
        
        varargout = {hC};
    

        
    %======================================================================
    case 'spm_list_cleanup'
    %======================================================================
        % cat_spm_results_ui('spm_list_cleanup',hReg)
        
        if ~mesh_detect, return; end
        if nargin>1
          spm_list_cleanup(varargin{2}); 
        else
          spm_list_cleanup
        end

      
    %======================================================================
    case 'setupsurfatlasmenu'                            %-Setup Atlas Menu
    %======================================================================
        % spm_results_ui('SetupSurfAtlasMenu',Finter)
        % RD202003 NOT WORKING YET 
        % requires update of spm_list('label',atlas) and spm_atlas  
    
        Finter = varargin{2};
        Vspm   = varargin{3};
        
        %hC  = uicontextmenu;
        hC   = uimenu(Finter,'Label','Atlas', 'Tag','AtlasUI');
        
        hC1  = uimenu(hC,'Label','Label using');
        
        %%
        satlases  = cat_get_defaults('extopts.satlas');
        expert    = cat_get_defaults('extopts.expertgui');
        listi = 1; 
        for si = 1:size(satlases,1)
            [pp,ff,ee] = spm_fileparts(satlases{si,1}); 
            if expert >= satlases{si,3}
                list(listi)      = struct('file',satlases{si,2},'name',satlases{si,1}); 
                list(listi).file = strrep(list(listi).file,'lh.','mesh.');
                if Vspm.dim(1) == 64984
                   list(listi).file = strrep(list(listi).file,'atlases_surfaces','atlases_surfaces_32k'); 
                end
                listi = listi + 1;
            end
        end
        for i=1:numel(list)
            uimenu(hC1,'Label',list(i).name,...
                'Callback',sprintf('spm_list(''label'',''%s''); cat_spm_results_ui(''spm_list_cleanup''); ',list(i).name));
        end
        if isempty(list), set(hC1,'Enable','off'); end
        
        %hC2  = uimenu(hC,'Label','Download Atlas...',...
        %    'Separator','on',...
        %    'Callback','spm_atlas(''install'');');
        
        varargout = {hC};
    
        
    %======================================================================
    case 'drawxyzgui'                                   %-Draw XYZ GUI area
    %======================================================================
        % hFxyz = cat_spm_results_ui('DrawXYZgui',M,DIM,xSPM,xyz,hReg)
        if nargin<6,  hReg=spm_XYZreg('FindReg','Interactive');
        else hReg=varargin{6}; end
        if nargin < 5, xyz=[0;0;0]; else xyz=varargin{5}; end
        if nargin < 4, error('Insufficient arguments'), end
        DIM     = varargin{3};
        M       = varargin{2};
        xyz     = spm_XYZreg('RoundCoords',xyz,M,DIM);
 
        %-Font details
        %------------------------------------------------------------------
        WS      = spm('WinScale');
        FS      = spm('FontSizes');
        PF      = spm_platform('fonts');
 
        %-Create XYZ control objects
        %------------------------------------------------------------------
        hFxyz = uipanel('Parent',hReg,'Title','co-ordinates','Units','Pixels',...
            'Position',[005 005 265 040].*WS,...
            'BorderType','Beveledout',...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
 
        uicontrol('Parent',hReg,'Style','Text','String','x =',...
            'Position',[015 010 024 018].*WS,...
            'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
            'HorizontalAlignment','Center');
        hX = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(1)),...
            'ToolTipString','enter x-coordinate',...
            'Position',[039 010 056 020].*WS,...
            'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
            'HorizontalAlignment','Right',...
            'Tag','hX',...
            'Callback','cat_spm_results_ui(''EdWidCB'')');
 
        uicontrol('Parent',hReg,'Style','Text','String','y =',...
            'Position',[100 010 024 018].*WS,...
            'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
            'HorizontalAlignment','Center');
        hY = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(2)),...
            'ToolTipString','enter y-coordinate',...
            'Position',[124 010 056 020].*WS,...
            'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
            'HorizontalAlignment','Right',...
            'Tag','hY',...
            'Callback','cat_spm_results_ui(''EdWidCB'')');
 
        if DIM(3) ~= 1
        uicontrol('Parent',hReg,'Style','Text','String','z =',...
            'Position',[185 010 024 018].*WS,...
            'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
            'HorizontalAlignment','Center');
        hZ = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(3)),...
            'ToolTipString','enter z-coordinate',...
            'Position',[209 010 056 020].*WS,...
            'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
            'HorizontalAlignment','Right',...
            'Tag','hZ',...
            'Callback','cat_spm_results_ui(''EdWidCB'')');
        else
        hZ = [];
        end
        
        %-Statistic value reporting pane
        %------------------------------------------------------------------
        hPan = uipanel('Parent',hReg,'Title','statistic','Units','Pixels',...
            'Position',[275 005 110 040].*WS,...
            'BorderType','Beveledout', ...
            'ShadowColor',[0.5 0.5 0.5],...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
        hSPM = uicontrol('Parent',hPan,'Style','Text','String','',...
            'Position',[005 001 100 020].*WS,...
            'FontSize',FS(10),...
            'HorizontalAlignment','Center');
 
 
        %-Store data
        %------------------------------------------------------------------
        set(hFxyz,'Tag','hFxyz','UserData',struct(...
            'hReg', [],...
            'M',    M,...
            'DIM',  DIM,...
            'XYZ',  varargin{4}.XYZmm,...
            'Z',    varargin{4}.Z,...
            'hX',   hX,...
            'hY',   hY,...
            'hZ',   hZ,...
            'hSPM', hSPM,...
            'xyz',  xyz ));
 
        set([hX,hY,hZ],'UserData',hFxyz)
        varargout = {hFxyz};
 
 
    %======================================================================
    case 'edwidcb'                          %-Callback for editable widgets
    %======================================================================
        % cat_spm_results_ui('EdWidCB')
 
        hC    = gcbo;
        d     = find(strcmp(get(hC,'Tag'),{'hX','hY','hZ'}));
        hFxyz = get(hC,'UserData');
        UD    = get(hFxyz,'UserData');
        xyz   = UD.xyz;
        nxyz  = xyz;
 
        o = evalin('base',['[',get(hC,'String'),']'],'sprintf(''error'')');
        if ischar(o) || length(o)>1
            warning(sprintf('%s: Error evaluating ordinate:\n\t%s',...
                mfilename,lasterr))
        else
            nxyz(d) = o;
            nxyz = spm_XYZreg('RoundCoords',nxyz,UD.M,UD.DIM);
        end
 
        if abs(xyz(d)-nxyz(d))>0
            UD.xyz = nxyz; set(hFxyz,'UserData',UD)
            if ~isempty(UD.hReg), spm_XYZreg('SetCoords',nxyz,UD.hReg,hFxyz); end
            set(hC,'String',sprintf('%.3f',nxyz(d)))
            cat_spm_results_ui('UpdateSPMval',UD)
        end
 
        
    %======================================================================
    case 'updatespmval'                           %-Update SPM value in GUI
    %======================================================================
        % cat_spm_results_ui('UpdateSPMval',hFxyz)
        % cat_spm_results_ui('UpdateSPMval',UD)
        if nargin<2, error('insufficient arguments'), end
        if isstruct(varargin{2}), UD=varargin{2}; else UD = get(varargin{2},'UserData'); end
        i  = spm_XYZreg('FindXYZ',UD.xyz,UD.XYZ);
        if isempty(i), str = ''; else str = sprintf('%6.2f',UD.Z(i)); end
        set(UD.hSPM,'String',str);
 
 
    %======================================================================
    case 'getcoords'             % Get current co-ordinates from XYZ widget
    %======================================================================
        % xyz = cat_spm_results_ui('GetCoords',hFxyz)
        if nargin<2, hFxyz='Interactive'; else hFxyz=varargin{2}; end
        hFxyz     = cat_spm_results_ui('FindXYZframe',hFxyz);
        varargout = {getfield(get(hFxyz,'UserData'),'xyz')};
 
 
    %======================================================================
    case 'setcoords'                       % Set co-ordinates to XYZ widget
    %======================================================================
        % [xyz,d] = cat_spm_results_ui('SetCoords',xyz,hFxyz,hC)
        if nargin<4, hC=NaN; else hC=varargin{4}; end
        if nargin<3, hFxyz=cat_spm_results_ui('FindXYZframe'); else hFxyz=varargin{3}; end
        if nargin<2, error('Set co-ords to what!'); else xyz=varargin{2}; end
 
        %-If this is an internal call, then don't do anything
        if isequal(hFxyz,hC), return, end
 
        UD = get(hFxyz,'UserData');
 
        %-Check validity of coords only when called without a caller handle
        %------------------------------------------------------------------
        if ~ishandle(hC)
            [xyz,d] = spm_XYZreg('RoundCoords',xyz,UD.M,UD.DIM);
            if d>0 && nargout<2, warning(sprintf(...
                '%s: Co-ords rounded to nearest voxel centre: Discrepancy %.2f',...
                mfilename,d))
            end
        else
            d = [];
        end
 
        %-Update xyz information & widget strings
        %------------------------------------------------------------------
        UD.xyz = xyz; set(hFxyz,'UserData',UD)
        set(UD.hX,'String',sprintf('%.2f',xyz(1)))
        set(UD.hY,'String',sprintf('%.2f',xyz(2)))
        set(UD.hZ,'String',sprintf('%.2f',xyz(3)))
        cat_spm_results_ui('UpdateSPMval',UD);
 
        %-Tell the registry, if we've not been called by the registry...
        %------------------------------------------------------------------
        if (~isempty(UD.hReg) && ~isequal(UD.hReg,hC))
            spm_XYZreg('SetCoords',xyz,UD.hReg,hFxyz);
        end
 
        %-Return arguments
        %------------------------------------------------------------------
        varargout = {xyz,d};
 

    %======================================================================
    case 'findxyzframe'                                  % Find hFxyz frame
    %======================================================================
        % hFxyz = cat_spm_results_ui('FindXYZframe',h)
        % Sorts out hFxyz handles
        if nargin<2, h='Interactive'; else h=varargin{2}; end
        if ischar(h), h=spm_figure('FindWin',h); end
        if ~ishandle(h), error('invalid handle'), end
        if ~strcmp(get(h,'Tag'),'hFxyz'), h=findobj(h,'Tag','hFxyz'); end
        if isempty(h), error('XYZ frame not found'), end
        if length(h)>1, error('Multiple XYZ frames found'), end
        varargout = {h};


    %======================================================================
    case 'plotui'                               %-GUI for plot manipulation
    %======================================================================
        % cat_spm_results_ui('PlotUi',hAx)
        if nargin<2, hAx=gca; else hAx=varargin{2}; end
 
        WS = spm('WinScale');
        FS = spm('FontSizes');
        Finter=spm_figure('FindWin','Interactive');
        figure(Finter)
 
        %-Check there aren't already controls!
        %------------------------------------------------------------------
        hGraphUI = findobj(Finter,'Tag','hGraphUI');
        if ~isempty(hGraphUI)           %-Controls exist
            hBs = get(hGraphUI,'UserData');
            if hAx==get(hBs(1),'UserData')  %-Controls linked to these axes
                return
            else                %-Old controls remain
                delete(findobj(Finter,'Tag','hGraphUIbg'))
            end
        end
 
        %-Frames & text
        %------------------------------------------------------------------
        hGraphUIbg = uipanel('Parent',Finter,'Title','','Tag','hGraphUIbg',...
            'BackgroundColor',spm('Colour'),...
            'BorderType','Line', 'HighlightColor',[0 0 0],...
            'Units','Pixels','Position',[001 195 400 055].*WS);
        hGraphUI   = uipanel('Parent',hGraphUIbg,'Title','','Tag','hGraphUI',...
            'BorderType','Etchedin', ...
            'BackgroundColor',[179 179 179]/255,...
            'Units','Pixels','Position',[005 005 390 046].*WS);
        hGraphUIButtsF = uipanel('Parent',hGraphUI,'Title','plot controls',...
            'Units','Pixels','Position',[005 005 380 039].*WS,...
            'BorderType','Beveledout', ...
            'FontAngle','Italic',...
            'FontSize',FS(10),...
            'ForegroundColor',[1 1 1],...
            'BackgroundColor',[179 179 179]/255);
 
        %-Controls
        %------------------------------------------------------------------
        h1 = uicontrol('Parent',hGraphUIButtsF,'Style','CheckBox','String','hold',...
            'ToolTipString','toggle hold to overlay plots',...
            'FontSize',FS(10),...
            'Value',double(strcmp(get(hAx,'NextPlot'),'add')),...
            'Callback',[...
            'if get(gcbo,''Value''), ',...
            'set(get(gcbo,''UserData''),''NextPlot'',''add''), ',...
            'else, ',...
            'set(get(gcbo,''UserData''),''NextPlot'',''replace''), ',...
            'end'],...
            'Interruptible','on','Enable','on',...
            'Tag','holdButton',...
            'Position',[005 005 070 020].*WS);
        h2 = uicontrol('Parent',hGraphUIButtsF,'Style','CheckBox','String','grid',...
            'ToolTipString','toggle axes grid',...
            'FontSize',FS(10),...
            'Value',double(strcmp(get(hAx,'XGrid'),'on')),...
            'Callback',[...
            'if get(gcbo,''Value''), ',...
            'set(get(gcbo,''UserData''),''XGrid'',''on'','...
            '''YGrid'',''on'',''ZGrid'',''on''), ',...
            'else, ',...
            'set(get(gcbo,''UserData''),''XGrid'',''off'','...
            '''YGrid'',''off'',''ZGrid'',''off''), ',...
            'end'],...
            'Interruptible','on','Enable','on',...
            'Position',[080 005 070 020].*WS);
        h3 = uicontrol('Parent',hGraphUIButtsF,'Style','CheckBox','String','Box',...
            'ToolTipString','toggle axes box',...
            'FontSize',FS(10),...
            'Value',double(strcmp(get(hAx,'Box'),'on')),...
            'Callback',[...
            'if get(gcbo,''Value''), ',...
            'set(get(gcbo,''UserData''),''Box'',''on''), ',...
            'else, ',...
            'set(get(gcbo,''UserData''),''Box'',''off''), ',...
            'end'],...
            'Interruptible','on','Enable','on',...
            'Position',[155 005 070 020].*WS);
        h4 = uicontrol('Parent',hGraphUIButtsF,'Style','popupmenu',...
            'ToolTipString','edit axis text annotations',...
            'FontSize',FS(10),...
            'String',{'text','Title','Xlabel','Ylabel'},...
            'Callback','cat_spm_results_ui(''PlotUiCB'')',...
            'Interruptible','on','Enable','on',...
            'Position',[230 005 070 020].*WS);
        h5 = uicontrol('Parent',hGraphUIButtsF,'Style','popupmenu',...
            'ToolTipString','change various axes attributes',...
            'FontSize',FS(10),...
            'String',{'attrib','LineWidth','XLim','YLim','handle'},...
            'Callback','cat_spm_results_ui(''PlotUiCB'')',...
            'Interruptible','off','Enable','on',...
            'Position',[305 005 070 020].*WS);
 
        %-Handle storage for linking, and DeleteFcns for linked deletion
        %------------------------------------------------------------------
        set(hGraphUI,'UserData',[h1,h2,h3,h4,h5])
        set([h1,h2,h3,h4,h5],'UserData',hAx)
 
        set(hGraphUIbg,'UserData',...
            [hGraphUI,hGraphUIButtsF,h1,h2,h3,h4,h5],...
            'DeleteFcn','cat_spm_results_ui(''Delete'',get(gcbo,''UserData''))')
        set(hAx,'UserData',hGraphUIbg,...
            'DeleteFcn','cat_spm_results_ui(''Delete'',get(gcbo,''UserData''))')


    %======================================================================
    case 'plotuicb'
    %======================================================================
        % cat_spm_results_ui('PlotUiCB')
        hPM = gcbo;
        v   = get(hPM,'Value');
        if v==1, return, end
        str = cellstr(get(hPM,'String'));
        str = str{v};
 
        hAx = get(hPM,'UserData');
        switch str
            case 'Title'
                h = get(hAx,'Title');
                set(h,'String',spm_input('Enter title:',-1,'s+',get(h,'String')))
            case 'Xlabel'
                h = get(hAx,'Xlabel');
                set(h,'String',spm_input('Enter X axis label:',-1,'s+',get(h,'String')))
            case 'Ylabel'
                h = get(hAx,'Ylabel');
                set(h,'String',spm_input('Enter Y axis label:',-1,'s+',get(h,'String')))
            case 'LineWidth'
                lw = spm_input('Enter LineWidth',-1,'e',get(hAx,'LineWidth'),1);
                set(hAx,'LineWidth',lw)
            case 'XLim'
                XLim = spm_input('Enter XLim',-1,'e',get(hAx,'XLim'),[1,2]);
                set(hAx,'XLim',XLim)
            case 'YLim'
                YLim = spm_input('Enter YLim',-1,'e',get(hAx,'YLim'),[1,2]);
                set(hAx,'YLim',YLim)
            case 'handle'
                varargout={hAx};
            otherwise
                warning(['Unknown action: ',str])
        end
 
        set(hPM,'Value',1)
 
 
    %======================================================================
    case 'clear'                                    %-Clear results subpane
    %======================================================================
        % Fgraph = cat_spm_results_ui('Clear',F,mode)
        % mode 1 [default] usual, mode 0 - clear & hide Res stuff, 2 - RNP
 
        if nargin<3, mode=1; else mode=varargin{3}; end
        if nargin<2, F='Graphics'; else F=varargin{2}; end
        F = spm_figure('FindWin',F);
 
        %-Clear input objects from 'Interactive' window
        %------------------------------------------------------------------
        %spm_input('!DeleteInputObj')
 
 
        %-Get handles of objects in Graphics window & note permanent results objects
        %------------------------------------------------------------------
        H = get(F,'Children');                          %-Get contents of window
        H = findobj(H,'flat','HandleVisibility','on');  %-Drop GUI components
        h = findobj(H,'flat','Tag','PermRes');          %-Look for 'PermRes' object
 
        if ~isempty(h)
            %-Found 'PermRes' object
            % This has handles of permanent results objects in it's UserData
            tmp  = get(h,'UserData');
            HR   = tmp.H;
            HRv  = tmp.Hv;
        else
            %-No trace of permanent results objects
            HR   = [];
            HRv  = {};
        end
        H = setdiff(H,HR);              %-Drop permanent results obj
 
 
        %-Delete stuff as appropriate
        %------------------------------------------------------------------
        if mode==2  %-Don't delete axes with NextPlot 'add'
            H = setdiff(H,findobj(H,'flat','Type','axes','NextPlot','add'));
        end
 
        delete(H)
        %set(F,'resize','on');set(F,'resize','off')
 
        if mode==0  %-Hide the permanent results section stuff
            set(HR,'Visible','off')
        else
            set(HR,{'Visible'},HRv)
        end
 
        
    %======================================================================
    case 'close'                                            %-Close Results
    %======================================================================
        set(spm_figure('GetWin','Graphics'),'Color',[1 1 1]);
        spm_clf('Interactive');
        spm_clf('Graphics');
        close(spm_figure('FindWin','Satellite'));
        evalin('base','clear');
    
    
    %======================================================================
    case 'launchmp'                            %-Launch multiplanar toolbox
    %======================================================================
        % hMP = cat_spm_results_ui('LaunchMP',M,DIM,hReg,hBmp)
        if nargin<5, hBmp = gcbo; else hBmp = varargin{5}; end
        hReg = varargin{4};
        DIM  = varargin{3};
        M    = varargin{2};
 
        %-Check for existing MultiPlanar toolbox
        hMP  = get(hBmp,'UserData');
        if ishandle(hMP)
            figure(ancestor(hMP,'figure'));
            varargout = {hMP};
            return
        end
 
        %-Initialise and cross-register MultiPlanar toolbox
        hMP = spm_XYZreg_Ex2('Create',M,DIM);
        spm_XYZreg('Xreg',hReg,hMP,'spm_XYZreg_Ex2');
 
        %-Setup automatic deletion of MultiPlanar on deletion of results controls
        set(hBmp,'Enable','on','UserData',hMP)
        set(hBmp,'DeleteFcn','cat_spm_results_ui(''delete'',get(gcbo,''UserData''))')
 
        varargout = {hMP};
 
 
    %======================================================================
    case 'delete'                           %-Delete HandleGraphics objects
    %======================================================================
        % cat_spm_results_ui('Delete',h)
        h = varargin{2};
        delete(h(ishandle(h)));
 
 
    %======================================================================
    otherwise
    %======================================================================
        error('Unknown action string')
 
end

%==========================================================================
function varargout = call_list(varargin)
%==========================================================================
global use_tfce

if use_tfce
  varargout = {tfce_list(varargin{:})};
else
  varargout = {spm_list(varargin{:})};
end

%==========================================================================
function mycluster(obj,evt,xSPM,distmin,nbmax)
%==========================================================================

spm_get_defaults('stats.results.volume.distmin',distmin);
spm_get_defaults('stats.results.volume.nbmax',nbmax);

hReg = evalin('base','hReg;');
TabDat = call_list('List',xSPM,hReg);
figure(spm_figure('GetWin','Interactive'));

%==========================================================================
function mychgcon(obj,evt,xSPM)
%==========================================================================
global use_tfce

xSPM2.swd   = xSPM.swd;
try, xSPM2.units = xSPM.units; end
xSPM2.Ic    = getfield(get(obj,'UserData'),'Ic');
if isempty(xSPM2.Ic) || all(xSPM2.Ic == 0), xSPM2 = rmfield(xSPM2,'Ic'); end
xSPM2.Im    = xSPM.Im;
xSPM2.pm    = xSPM.pm;
xSPM2.Ex    = xSPM.Ex;
xSPM2.title = '';
if ~isempty(xSPM.thresDesc)
    if strcmp(xSPM.STAT,'P')
        % These are soon overwritten by spm_getSPM
        xSPM2.thresDesc = xSPM.thresDesc;
        xSPM2.u = xSPM.u;
        xSPM2.k = xSPM.k;
        % xSPM.STATstr contains Gamma
    else
        td = regexp(xSPM.thresDesc,'p\D?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
        if isempty(td)
            td = regexp(xSPM.thresDesc,'\w=(?<u>[\.\d]+)','names');
            td.thresDesc = 'none';
        end
        if strcmp(td.thresDesc,'unc.'), td.thresDesc = 'none'; end
        xSPM2.thresDesc = td.thresDesc;
        if use_tfce
            xSPM2.invResult = xSPM.invResult;
        end
        xSPM2.statType = xSPM.statType;
        xSPM2.u     = str2double(td.u);
        xSPM2.k     = xSPM.k;
    end
end
hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
xyz  = spm_XYZreg('GetCoords',hReg);
[hReg,xSPM,SPM, TabDat] = cat_spm_results_ui('setup',xSPM2);
spm_XYZreg('SetCoords',xyz,hReg);
spm_list_cleanup;
assignin('base','hReg',hReg);
assignin('base','xSPM',xSPM);
assignin('base','SPM',SPM);
assignin('base','TabDat',TabDat);
figure(spm_figure('GetWin','Interactive'));

%==========================================================================
function mycheckres(obj,evt,xSPM)
%==========================================================================
spm_check_results([],xSPM);

%==========================================================================
function mysavespm(action)
%==========================================================================
xSPM = evalin('base','xSPM;');
XYZ  = xSPM.XYZ;

switch lower(action)
    case 'thresh'
        Z = xSPM.Z;
        
    case 'binary'
        Z = ones(size(xSPM.Z));
        
    case 'n-ary'
        if ~isfield(xSPM,'G')
            Z       = spm_clusters(XYZ);
            num     = max(Z);
            [n, ni] = sort(histc(Z,1:num), 2, 'descend');
            n       = size(ni);
            n(ni)   = 1:num;
            Z       = n(Z);
        else
            C       = NaN(1,size(xSPM.G.vertices,1));
            C(xSPM.XYZ(1,:)) = ones(size(xSPM.Z));
            C       = spm_mesh_clusters(xSPM.G,C);
            Z       = C(xSPM.XYZ(1,:));
        end
        
    case 'current'
        [xyzmm,i] = spm_XYZreg('NearestXYZ',...
            cat_spm_results_ui('GetCoords'),xSPM.XYZmm);
        cat_spm_results_ui('SetCoords',xSPM.XYZmm(:,i));
        
        if ~isfield(xSPM,'G')
            A   = spm_clusters(XYZ);
            j   = find(A == A(i));
            Z   = ones(1,numel(j));
            XYZ = xSPM.XYZ(:,j);
        else
            C   = NaN(1,size(xSPM.G.vertices,1));
            C(xSPM.XYZ(1,:)) = ones(size(xSPM.Z));
            C   = spm_mesh_clusters(xSPM.G,C);
            C   = C==C(xSPM.XYZ(1,i));
            Z   = C(xSPM.XYZ(1,:));
        end
        
    otherwise
        error('Unknown action.');
end

if isfield(xSPM,'G')
    F     = spm_input('Output filename',1,'s');
    if isempty(spm_file(F,'ext'))
        F = spm_file(F,'ext','.gii');
    end
    F     = spm_file(F,'CPath');
    M     = gifti(xSPM.G);
    C     = zeros(1,size(xSPM.G.vertices,1));
    C(xSPM.XYZ(1,:)) = Z; % or use NODE_INDEX
    M.cdata = C;
    save(M,F);
    if exist('cat_surf_render','file')
        cmd   = 'cat_surf_render(''Disp'',''%s'')';
    else
        cmd   = 'spm_mesh_render(''Disp'',''%s'')';
    end
else
    V   = spm_write_filtered(Z, XYZ, xSPM.DIM, xSPM.M,...
    sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k));
    cmd = 'spm_image(''display'',''%s'')';
    F   = V.fname;
end

fprintf('Written %s\n',spm_file(F,'link',cmd));                         %-#

%==========================================================================
function myslover
%==========================================================================
spm_input('!DeleteInputObj');
xSPM = evalin('base','xSPM;');

so = slover;
[img,sts] = spm_select(1,'image','Select image for rendering on');
if ~sts, return; end
so.img.vol = spm_vol(img);
%obj.img.type = 'truecolour';
%obj.img.cmap = gray;
%[mx,mn] = slover('volmaxmin', obj.img.vol);
%obj.img.range = [mn mx];
so.img.prop = 1;

so = add_spm(so,xSPM);

so.transform = deblank(spm_input('Image orientation', '+1', ...
    'Axial|Coronal|Sagittal', char('axial','coronal','sagittal'), 1));
so = fill_defaults(so);
slices = so.slices;
so.slices = spm_input('Slices to display (mm)', '+1', 'e', ...
    sprintf('%0.0f:%0.0f:%0.0f',slices(1),mean(diff(slices)),slices(end)));

so.figure = spm_figure('GetWin', 'SliceOverlay');
so = paint(so);
assignin('base','so',so);

%==========================================================================
function spm_list_cleanup(hReg)
%==========================================================================
global mesh_detect

if ~mesh_detect, return; end

hRes.Fgraph       = [spm_figure('FindWin','Graphics'),spm_figure('FindWin','Satellite')];

% fine red lines of the SPM result table
hRes.Fline        = findobj(hRes.Fgraph,'Type','Line','Tag','');% ,'UIcontextMenu',[]);
hRes.FlineAx      = get(hRes.Fline,'parent');

set(hRes.Fline,'HitTest','off'); %
for axi = 1:numel( hRes.FlineAx ),  rotate3d(hRes.FlineAx{axi},'off'); end
for axi = 1:numel( hRes.FlineAx ),  set(hRes.FlineAx{axi},'visible','off'); end
 
%%
hRes.Img        = get(findobj(hRes.Fgraph,'Type','Image','Tag','Transverse'),'parent');
for axi = 1:numel( hRes.Img  ),  rotate3d(hRes.Img{axi},'off'); end

%% find the SPM string within the surface axis 
hRes.Ftext        = findobj(hRes.Fgraph,'Type','Text');
stext             = get(hRes.Ftext,'String'); 
hRes.Ftext3dspm   = findobj(hRes.Fgraph,'Type','Text','String', ...
	stext{ find(~cellfun('isempty',strfind(stext,'SPM\{'))) } );
set(hRes.Ftext3dspm,'visible','off','HitTest','off');

%% get backgroundcolor
bgc = get(spm_figure('FindWin','Graphics'),'Color'); 
% get low contrast texts 
Ftextcol = cell2mat(get(hRes.Ftext,'Color'));
%% invert text
rms = @(x) mean(x.^2,2).^0.5; 
bgcdist = rms(Ftextcol - repmat(bgc, numel(hRes.Ftext), 1)); 
bgcdist = abs(bgcdist)<0.3;
if ~isempty(bgcdist) 
	for fi = 1:numel(bgcdist)
		if bgcdist(fi)>0
			set( hRes.Ftext( fi ) , 'Color' , min(1,max(0, 1 - Ftextcol( fi ,:))) );  
		end
	end
end

%% update spm_XYZreg XYZ update function
if ~exist('hReg','var')
	ListXYZ=findobj('ButtonDownFcn','spm_XYZreg(''SetCoords'',get(gcbo,''UserData''),hReg,1);');
	for i=1:numel(ListXYZ)
		set(ListXYZ(i),'ButtonDownFcn',[get(ListXYZ(i),'ButtonDownFcn') ' cat_spm_results_ui(''spm_list_cleanup'');']); 
	end
end
%==========================================================================

