function str = cat_main_reportstr(job,res,qa,cat_warnings) 
% ______________________________________________________________________
% 
% Prepare text output for CAT report. This function may heavily change
% due to parameter changes. However, this only effects the output of the
% CAT report. Called from cat_main. 
%
%   str = cat_main_reportstr(job,res,qa,cat_warnings)
%
%   str          .. cellstrings with 3 elements with two strings
%   str{1}       .. full width table 
%   str{2:3}     .. half width table 
%
%   job          .. SPM/CAT parameter structure
%   res          .. SPM preprocessing structure
%   qa           .. CAT quality and subject information structure
%   cat_warnings .. CAT warning structure
%
%   See also cat_main_reportfig and cat_main_reportcmd.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_reportstr.m 1577 2020-03-09 17:36:03Z dahnke $

%#ok<*AGROW>

  QMC   = cat_io_colormaps('marks+',17);
  color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
   
  %mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
  marks2str   = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
  mark2rps    = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
  grades      = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  mark2grad   = @(mark) grades{min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3)))};


  % CAT GUI parameter:
  % --------------------------------------------------------------------
  str{1} = [];

  % use red output if a beta version is used
  catv = qa.software.revision_cat; isbeta = strfind(lower(catv),'beta'); 
  if ~isempty(isbeta), catv = [catv(1:isbeta-1) '\color[rgb]{0.8 0 0}' catv(isbeta:isbeta+3 ) '\color[rgb]{0.8 0 0}' catv(isbeta+4:end)]; end
  
  % ExploreASL Hack
  if ~isfield(qa.software,'version_matlab')
      dotIndex = strfind(version,'.');
      version_matlab = version;
      qa.software.version_matlab = version_matlab(1:dotIndex(1)+1);
  end
    
  % 1 line: Matlab, SPM12, CAT12 version number and GUI and experimental mode 
  str{1} = [str{1} struct('name', 'Version: OS / Matlab / SPM12 / CAT12:','value',...
    sprintf('%s / %s / %s / %s (%s)',qa.software.system,qa.software.version_matlab,...
    qa.software.version_spm,qa.software.version_cat,catv))];
  % add CAT segmentation version if not current
  if ~isempty(qa.software.version_segment)
    str{1}(end).name = [str{1}(end).name(1:end-1) ' / seg:']; 
    str{1}(end).value = [str{1}(end).value ' / \color[rgb]{0 0.2 1}' qa.software.version_segment '']; 
  end 
  % write GUI mode
  if     job.extopts.expertgui==1, str{1}(end).value = [str{1}(end).value '\bf\color[rgb]{0 0.2 1}e']; 
  elseif job.extopts.expertgui==2, str{1}(end).value = [str{1}(end).value '\bf\color[rgb]{0 0.2 1}d'];
  end  
  % write experimental flag
  if str{1}(end).value(end)==')'
    if job.extopts.experimental, str{1}(end).value = [str{1}(end).value(1:end-1) '\bf\color[rgb]{0 0.2 1}x\color[rgb]{0 0 0})']; end  
  else
    if job.extopts.experimental, str{1}(end).value = [str{1}(end).value '\bf\color[rgb]{0 0.2 1}x']; end  
  end

  
  % 2 lines: TPM, Template, Normalization method with voxel size
  str{1} = [str{1} struct('name', 'Tissue Probability Map:','value',strrep(spm_str_manip(res.tpm(1).fname,'k40d'),'_','\_'))];
  if res.do_dartel
    if job.extopts.regstr==0 % Dartel
      str{1} = [str{1} struct('name', 'Dartel Registration to: ',...
                        'value',strrep(spm_str_manip(job.extopts.darteltpm{1},'k40d'),'_','\_'))];
    elseif job.extopts.regstr==4 % Dartel
      str{1} = [str{1} struct('name', 'Shooting Registration to: ',...
                        'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
    else
      if job.extopts.expertgui==0
        str{1} = [str{1} struct('name','Optimized Shooting Registration to:',...
                          'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
      else
        str{1} = [str{1} struct('name', sprintf('Optimized Shooting Registration (regstr:%s) to :',sprintf('%g ',job.extopts.regstr)),...
                          'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
      end
    end
  end

  % 1 line 1: Affreg
  str{1} = [str{1} struct('name', 'affreg:','value',sprintf('%s',job.opts.affreg))];
  % 1 line 2: APP

  APPstr = {'none','light','full','','','animal'}; APPstr{1071} = 'rough'; APPstr{1145} = 'rough(new)'; 
  str{1}(end).name  = [str{1}(end).name(1:end-1) ' / APP ']; 
  str{1}(end).value = [str{1}(end).value sprintf(' / %s',APPstr{job.extopts.APP+1})];

  % 1 line 3: biasstr / biasreg+biasfwhm
  if job.opts.biasstr>0
    biasstr = {'ultralight','light','medium','strong','heavy'};
    str{1}(end).name  = [str{1}(end).name(1:end-1) ' / biasstr '];  
    str{1}(end).value = [str{1}(end).value sprintf(' / %s',biasstr{round(job.opts.biasstr*4)+1})];
    if job.extopts.expertgui % add the value
      str{1}(end).value = [str{1}(end).value sprintf('(%0.2f;breg:%0.2f;bfwhm:%0.2f)',job.opts.biasstr,job.opts.biasreg,job.opts.biasfwhm)]; 
    end
  else
    str{1}(end).name  = [str{1}(end).name(1:end-1) ' / biasreg / biasfwhm'];
    str{1}(end).value = [str{1}(end).value sprintf(' / %0.2f / %0.2f',job.opts.biasreg,job.opts.biasfwhm)]; 
  end
  if isfield(job.opts,'acc') && job.opts.acc>0
    str{1} = [str{1} struct('name', '','value','')];
    accstr = {'ultra low','low','std','high','ultra high'};
    str{1}(end).name  = [str{1}(end).name(1:end-1) 'SPM accuracy (samp/tol) '];  
    str{1}(end).value = [str{1}(end).value sprintf(' / %s',accstr{round(job.opts.acc*4)+1})];
    if job.extopts.expertgui % add the value
      str{1}(end).value = [str{1}(end).value sprintf('%0.2f (%0.2f/%0.2f)',job.opts.acc,job.opts.samp,job.opts.tol)]; 
    end
  else
    if job.extopts.expertgui
      str{1} = [str{1} struct('name', '','value','')];
      str{1}(end).name  = [str{1}(end).name(1:end-1) 'SPM accuracy (samp/tol)'];
      str{1}(end).value = [str{1}(end).value sprintf('%0.2f / %0.2f',job.opts.samp,job.opts.tol)]; 
    end
  end


  % 1 line: adaptive noise parameter ( MRFstr + SANLM + NCstr )
  NCstr.labels = {'none','full','light','medium','strong','heavy'};
  NCstr.values = {0 1 2 -inf 4 5}; 
  defstr  = {'none','ultralight','light','medium','strong','heavy',... sanlm vs. isarnlm
             'ultralight+','ultralight+','light+','medium+','strong+','heavy+'};
  defstrm = @(x) defstr{ round(max(0,min(2,x))*4) + 1 + (x>0) + (x>1)};
  str{1} = [str{1} struct('name', 'Noise reduction:','value','')]; 
  if job.extopts.NCstr==0 
    if job.extopts.mrf==0
      str{1}(end).value = 'no noise correction';
    else
      if job.extopts.expertgui==0
        str{1}(end).value = 'MRF'; 
      else
        str{1}(end).value = sprintf('MRF(%0.2f)',job.extopts.mrf); 
      end  
    end
  else
    str{1}(end).value = sprintf('SANLM(%s)',NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')});
  end

  if job.extopts.NCstr~=0 && job.extopts.mrf
    if job.extopts.expertgui==0
      str{1}(end).value = '+MRF'; 
    else
      str{1}(end).value = sprintf('+MRF(%0.2f)',job.extopts.mrf); 
    end 
  end


  % 1 line(s): LASstr / GCUTstr / CLEANUPstr
  str{1}(end).name  = 'LASstr / GCUTstr / CLEANUPstr:';
  gcutstr  = {'none','SPM','GCUT','APRG'};  
  if ~job.extopts.expertgui
    str{1}(end).value = sprintf('%s / %s / %s',defstrm(job.extopts.LASstr),...
      gcutstr{ceil(job.extopts.gcutstr+2)},defstrm(job.extopts.cleanupstr)); 
  else
    str{1}(end).value = sprintf('%s(%0.2f) / %s(%0.2f) / %s(%0.2f)',...
      defstrm(job.extopts.LASstr),job.extopts.LASstr,gcutstr{ceil(job.extopts.gcutstr+2)},...
      job.extopts.gcutstr,defstrm(job.extopts.cleanupstr),job.extopts.cleanupstr); 
  end
  restype = char(fieldnames(job.extopts.restypes));
  if job.extopts.expertgui
    str{1} = [str{1} struct('name', 'KAMAP / WMHC / SLC / collcorr / restype:','value',...
           sprintf('%d / %d / %d / %d / %s',...
          job.extopts.spm_kamap,job.extopts.WMHC,job.extopts.SLC,job.extopts.collcorr,restype))];
  else
    str{1} = [str{1} struct('name', 'restype:','value',sprintf('%s',restype))];
  end
  if ~strcmp('native',restype)
    str{1}(end).value = [str{1}(end).value sprintf('(%0.2f %0.2f)',job.extopts.restypes.(restype))];
  end 

  % line 8: surface parameter
  if job.output.surface
    str{1} = [str{1} struct('name', 'Voxel resolution (original > internal > PBT; vox):',...
           'value',sprintf('%4.2fx%4.2fx%4.2f > %4.2fx%4.2fx%4.2f > %4.2f mm%s; %4.2f mm ', ... 
           qa.qualitymeasures.res_vx_vol,qa.qualitymeasures.res_vx_voli,job.extopts.pbtres,native2unicode(179, 'latin1'),job.extopts.vox(1)))];
  else
    str{1} = [str{1} struct('name', 'Voxel resolution (original > intern; vox):',...
           'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s; %4.2f mm', ...
           qa.qualitymeasures.res_vx_vol,native2unicode(179, 'latin1'),qa.qualitymeasures.res_vx_voli,native2unicode(179, 'latin1'),job.extopts.vox(1)))];
  end       
  % str{1} = [str{1} struct('name', 'Norm. voxel size:','value',sprintf('%0.2f mm',job.extopts.vox))]; % does not work yet 


  % Image Quality measures:
  % --------------------------------------------------------------------
  str{2} =       struct('name', '\bfImage and Preprocessing Quality:','value',''); 
  str{2} = [str{2} struct('name',' Resolution:','value',marks2str(qa.qualityratings.res_RMS,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.res_RMS),mark2grad(qa.qualityratings.res_RMS))))];
  str{2} = [str{2} struct('name',' Noise:','value',marks2str(qa.qualityratings.NCR,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.NCR),mark2grad(qa.qualityratings.NCR))))];
  str{2} = [str{2} struct('name',' Bias:','value',marks2str(qa.qualityratings.ICR,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.ICR),mark2grad(qa.qualityratings.ICR))))]; % not important and more confusing 
  str{2} = [str{2} struct('name','\bf Weighted average (IQR):','value',marks2str(qa.qualityratings.IQR,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR))))];
  if isfield(qa.qualitymeasures,'SurfaceEulerNumber') && ~isempty(qa.qualitymeasures.SurfaceEulerNumber)
    if job.extopts.expertgui
      str{2} = [str{2} struct('name',' Mean surface Euler number:','value',marks2str(qa.qualityratings.SurfaceEulerNumber,...
                sprintf('%g', qa.qualitymeasures.SurfaceEulerNumber)))]; 
              
      if isfield(qa.qualitymeasures,'SurfaceDefectNumber') && ~isempty(qa.qualitymeasures.SurfaceDefectNumber)
        str{2}(end).name  = [str{2}(end).name(1:end-8)  ' / defect number:'];
        str{2}(end).value = [str{2}(end).value ' / ' marks2str(qa.qualityratings.SurfaceDefectNumber,...
                sprintf('%0.2f', qa.qualitymeasures.SurfaceDefectNumber)) ];
      end
            
    else
      str{2} = [str{2} struct('name',' Mean surface Euler number:','value',sprintf('%g', qa.qualitymeasures.SurfaceEulerNumber))]; 
    end
  end
  
  
  if isfield(qa.qualitymeasures,'SurfaceDefectArea') && ~isempty(qa.qualitymeasures.SurfaceDefectArea)
    if job.extopts.expertgui
      str{2} = [str{2} struct('name',' Mean topology defects size:','value',marks2str(qa.qualityratings.SurfaceDefectArea,...
                sprintf('%0.2f%%', qa.qualitymeasures.SurfaceDefectArea)))];

      if isfield(qa.qualitymeasures,'SurfaceSelfIntersections') && ~isempty(qa.qualitymeasures.SurfaceSelfIntersections)
        str{2}(end).name  = [str{2}(end).name(1:end-6)  ' / self-inters. size:'];
        str{2}(end).value = [str{2}(end).value ' / ' marks2str(qa.qualityratings.SurfaceSelfIntersections,...
                sprintf('%0.2f%%', qa.qualitymeasures.SurfaceSelfIntersections)) ];
      end
   
    else
      str{2} = [str{2} struct('name',' Mean size of topology defects:','value',sprintf('%0.2f%%', qa.qualitymeasures.SurfaceDefectArea))];
    end
    
  end
  
  if job.extopts.expertgui && isfield(qa.qualityratings,'SurfaceIntensityRMSE')
      str{2} = [str{2} struct('name',' Surface intensity / position RMSE:','value',[ marks2str( qa.qualityratings.SurfaceIntensityRMSE ,...
        sprintf('%0.2f', qa.qualitymeasures.SurfaceIntensityRMSE)) ' / ' ...
        marks2str( qa.qualityratings.SurfacePositionRMSE ,sprintf('%0.2f', qa.qualitymeasures.SurfacePositionRMSE) ) ] ) ];
  end

  % Subject Measures
  % --------------------------------------------------------------------
  % Volume measures

  % header
  str{3} = struct('name', '\bfVolumes:','value',sprintf('%5s %5s %5s ','CSF','GM','WM')); 
  if job.extopts.WMHC>1, str{3}(end).value = [str{3}(end).value sprintf('%5s ','WMH')]; end
  if job.extopts.SLC>0,  str{3}(end).value = [str{3}(end).value sprintf('%5s ','SL')];  end

  % absolute volumes
  str{3} = [str{3} struct('name', ' Absolute volume:','value',sprintf('%5.0f %5.0f %5.0f ', ...
          qa.subjectmeasures.vol_abs_CGW(1:3)))];
  if job.extopts.WMHC>1,  str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_abs_CGW(4))]; end 
  if job.extopts.SLC>0,   str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_abs_CGW(5))]; end
  str{3}(end).value = [str{3}(end).value 'cm' native2unicode(179, 'latin1')];

  % relative volumes
  str{3} = [str{3} struct('name', ' Relative volume:','value',sprintf('%5.1f %5.1f %5.1f ', ...
          qa.subjectmeasures.vol_rel_CGW(1:3)*100))];
  if job.extopts.WMHC>1,  str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_rel_CGW(4)*100)]; end 
  if job.extopts.SLC>0,   str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_rel_CGW(5)*100)]; end
  str{3}(end).value = [str{3}(end).value '%'];

  % warning if many WMH were found but no correction is active 
  if job.extopts.WMHC<2 && (qa.subjectmeasures.vol_rel_CGW(4)>0.03 || ...
     qa.subjectmeasures.vol_rel_CGW(4)/qa.subjectmeasures.vol_rel_CGW(3)>0.05)
    str{3}(end).value = [str{3}(end).value sprintf('\\bf\\color[rgb]{1 0 1} WMHs!')];  
  end
  
  str{3} = [str{3} struct('name', ' TIV:','value', sprintf(['%0.0f cm' native2unicode(179, 'latin1')],qa.subjectmeasures.vol_TIV))];  
  if isfield(qa.subjectmeasures,'surf_TSA') && job.extopts.expertgui>1
    str{3}(end).name  = [str{3}(end).name  ' / TSA:']; 
    str{3}(end).value = [str{3}(end).value sprintf(' / %0.0f cm%s' ,qa.subjectmeasures.surf_TSA,char(178))];  
  end

  % Surface measures - Thickness, (Curvature, Depth, ...)
  %if cellfun('isempty',strfind({Psurf(:).Pcentral},'ch.')), thstr = 'Cerebral Thickness'; else thstr = 'Thickness'; end
  thstr = 'Thickness';
  if isfield(qa.subjectmeasures,'dist_thickness') && ~isempty(qa.subjectmeasures.dist_thickness)
    str{3} = [str{3} struct('name', ['\bf' thstr ':'],'value',sprintf('%5.2f%s%5.2f mm', ...
           qa.subjectmeasures.dist_thickness{1}(1),native2unicode(177, 'latin1'),qa.subjectmeasures.dist_thickness{1}(2)))];
    if isfield(qa.subjectmeasures,'dist_gyruswidth') && ~isnan(qa.subjectmeasures.dist_gyruswidth{1}(1))
      str{3} = [str{3} struct('name', '\bfGyruswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_gyruswidth{1}(1),native2unicode(177, 'latin1'),qa.subjectmeasures.dist_gyruswidth{1}(2)))];
    end
    if isfield(qa.subjectmeasures,'dist_sulcuswidth') && ~isnan(qa.subjectmeasures.dist_sulcuswidth{1}(1))
      str{3} = [str{3} struct('name', '\bfSulcuswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_sulcuswidth{1}(1),native2unicode(177, 'latin1'),qa.subjectmeasures.dist_sulcuswidth{1}(2)))];
    end
  end

  % Preprocessing Time
  str{2} = [str{2} struct('name','\bfProcessing time:','value',sprintf('%02.0f:%02.0f min', ...
    floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60)))]; 

  % Warnings
  if numel(cat_warnings)>0 && job.extopts.expertgui>0
    str{2} = [str{2} struct('name', '','value','')]; 
    str{2} = [str{2} struct('name', '\bfWarnings:','value','')]; 
    for wi=1:numel(cat_warnings)
      shorter = cat_warnings(wi).identifier;
      % remove leading MATLAB, SPM or CAT elements
      dots    = max([min(strfind(shorter,'MATLAB')+7), ...
                     min(strfind(shorter,'SPM')+4), ...
                     min(strfind(shorter,'CAT')+4)]);
      if ~isempty(dots), shorter = shorter(dots:end); end
      % limit lenght of the string and replace critical character
      shorter = spm_str_manip(shorter,'l40');
      shorter = marks2str(4,shorter);
      shorter = strrep(shorter,'_','\_');
      str{2}    = [str{2} struct('name',shorter,'value','')];  
    end
  end

  % adding one space for correct printing of bold fonts
  for ssi = 1:numel(str)
    for si = 1:numel(str{ssi})
      str{ssi}(si).name   = [str{ssi}(si).name  '  '];   
      str{ssi}(si).value  = [str{ssi}(si).value '  '];
    end
  end

end