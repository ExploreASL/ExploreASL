function cat_main_reportcmd(job,res,qa)
% ______________________________________________________________________
% 
%   Display conclusion of CAT preprocessing in the command window and 
%   cleanup of some figures. 
%
%   cat_main_reportcom(job,res,qa)
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_reportcmd.m 1389 2018-11-11 10:39:41Z dahnke $

  VT0 = res.image0(1);
  [pth,nam] = spm_fileparts(VT0.fname); 

  % command window output
  QMC       = cat_io_colormaps('marks+',17);
  color     = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  mark2rps  = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
  grades    = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  mark2grad = @(mark) grades{min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3)))};
  
  fprintf('\n%s',repmat('-',1,72));
  fprintf(1,'\nCAT12 preprocessing took %0.0f minute(s) and %0.0f second(s).\n', ...
    floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60));
  cat_io_cprintf(color(QMC,qa.qualityratings.IQR), sprintf('Image Quality Rating (IQR):  %5.2f%%%% (%s)\n',...
    mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR)));

  % print subfolders
  %if job.extopts.subfolders
  %  fprintf('Segmentations are saved in %s%s%s\n',pth,filesep,'mri');
  %  fprintf('Reports are saved in %s%s%s\n',pth,filesep,'report');
  %  if job.output.ROI
  %    fprintf('Labels are saved in %s%s%s\n',pth,filesep,'label');
  %  end
  %  if job.output.surface && exist('Psurf','var') && ~isempty(Psurf)
  %    fprintf('Surface measurements are saved in %s%s%s\n',pth,filesep,'surf');
  %  end
  %end

  fprintf('%s\n\n',repmat('-',1,72));

  % finish diary entry of "../report/cmdln_*.txt"
  % read diary and add the command-line output to the *.xml and *.mat file
  % diary off; % ExploreASL fix - disable CAT12 outputs as we write everything to ExploreASL log file.
  try %#ok<TRYNC>
    fid  =fopen(res.catlog);
    txt = fread(fid,200000,'uint8=>char');
    fclose(fid); 
    txt2 = textscan(txt,'%s','Delimiter',''); 
    cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
      'catlog',txt2),'write+'); % here we have to use the write+!
  end    
  
  spm_progress_bar('Clear');

end
