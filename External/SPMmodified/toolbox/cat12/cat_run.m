function varargout = cat_run(job)
% Segment a bunch of images
% ______________________________________________________________________
%
%   FORMAT cat_run(job)
%
%   job.channel(n).vols{m}
%   job.channel(n).biasreg
%   job.channel(n).biasfwhm
%   job.channel(n).write
%   job.tissue(k).tpm
%   job.tissue(k).ngaus
%   job.tissue(k).native
%   job.tissue(k).warped
%
% See the user interface for a description of the fields.
%
% based on John Ashburners version of
% spm_preproc8_run.m 2281 2008-10-01 12:52:50Z john $
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_run.m 1615 2020-05-09 10:19:40Z gaser $

%#ok<*AGROW,*STRIFCND,*STRCL1,*ASGLU,*STREMP>

%rev = '$Rev: 1615 $';

%  -----------------------------------------------------------------
%  Lazy processing (expert feature)
%  -----------------------------------------------------------------
%  If N>10000 files were processed the crash of one of J jobs by 
%  small errors makes it hard to find the unprocess files. 
%  The lazy processing will only process files, if one of the output
%  is missed and if the same preprocessing options were used before.
%  -----------------------------------------------------------------

% disable parallel processing for only one subject
n_subjects = numel(job.data);
if n_subjects == 1, job.nproc = 0; end

%{ send Matlab version to server
if cat_get_defaults('extopts.send_info')
  urlinfo = sprintf('%s%s%s%s%s%s%d',cat_version,'%2F',computer,'%2F','processed',...
     '%2F',n_subjects);
  cat_io_send_to_server(urlinfo);
end
%}

if isfield(job.extopts,'admin') && isfield(job.extopts.admin,'lazy') && job.extopts.admin.lazy && ...
  ~isfield(job,'process_index') && isfield(job,'nproc') && job.nproc>.1 && (~isfield(job,'process_index'))  
  jobo.vout = vout_job(job);    % expected output
  jobl      = update_job(job);
  jobl.vout = vout_job(jobl);   
  job.data  = remove_already_processed(jobl); 
end

% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))  
  % rescue original subjects
  job_data = job.data;
  n_subjects = numel(job.data);
  if job.nproc > n_subjects
    job.nproc = n_subjects;
  end
  job.process_index = cell(job.nproc,1);

  % initial splitting of data
  for i=1:job.nproc
    job.process_index{i} = (1:job.nproc:(n_subjects-job.nproc+1))+(i-1);
  end

  % check if all data are covered
  for i=1:rem(n_subjects,job.nproc)
    job.process_index{i} = [job.process_index{i} n_subjects-i+1];
  end

  tmp_array = cell(job.nproc,1); job.printPID = 1; job.getPID = 2; 
    
  logdate   = datestr(now,'YYYYmmdd_HHMMSS');
  PID       = zeros(1,job.nproc);
  catSID    = zeros(1,job.nproc);
  for i=1:job.nproc
    jobo = job; 
    
    fprintf('\nRunning job %d:\n',i);
    for fi=1:numel(job_data(job.process_index{i}))
      fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
    end
    job.data = job_data(job.process_index{i});
         
    % temporary name for saving job information
    tmp_name = [tempname '.mat'];
    tmp_array{i} = tmp_name; 
    %def = cat_get_defaults; job = cat_io_checkinopt(job,def); % further job update required here to get the latest cat defaults
    spm12def = spm_get_defaults;  %#ok<NASGU>
    cat12def = cat_get_defaults;  %#ok<NASGU>
    save(tmp_name,'job','spm12def','cat12def');
    clear spm12def cat12;
    
    % matlab command, cprintferror=1 for simple printing         
    matlab_cmd = sprintf(...
        ['"global cprintferror; cprintferror=1; addpath %s %s %s %s;load %s; ' ...
         'global defaults; defaults=spm12def; clear defaults; '...
         'global cat; cat=cat12def; clear cat; cat_run(job); "'],...
      spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),...
        fullfile(spm('dir'),'toolbox','OldNorm'),fullfile(spm('dir'),'toolbox','DARTEL'), tmp_name);

    % log-file for output
    log_name{i} = ['catlog_main_' logdate '_log' sprintf('%02d',i) '.txt'];

    % call matlab with command in the background
    if ispc
      % check for spaces in filenames that will not work with windows systems and background jobs
      if strfind(spm('dir'),' ') 
        cat_io_cprintf('warn',...
            ['\nWARNING: No background processes possible because your SPM installation is located in \n' ...
             '         a folder that contains spaces. Please set the number of processes in the GUI \n'...
             '         to ''0''. In order to split your job into different processes,\n' ...
             '         please do not use any spaces in folder names!\n\n']);
         job.nproc = 0; 
         job = update_job(job);
         
         varargout{1} = run_job(job);
         return; 
      end
      % prepare system specific path for matlab
      export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
      [status,result] = system(export_cmd);
      system_cmd = ['start matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name{i}];
    else
      % -nodisplay .. nodisplay is without figure output > problem with CAT report ... was there a server problem with -nodesktop?
      system_cmd = [fullfile(matlabroot,'bin') '/matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name{i} ' 2>&1 & '];
    end
    [status,result] = system(system_cmd); 
    cat_check_system_output(status,result);
    
    
    
    %% look for existing files and extract their PID for later control  
    %  --------------------------------------------------------------------
    test    = 0; lim    = 200; ptime    = 0.5; % exist file?
    testpid = 0; limpid = 400; ptimepid = 2.0; % get PID
    ptimesid = 1 * 30;                        % update every minute? 
    while test<lim
      if ~exist(log_name{i},'file')
        pause(ptime); 
        test = test + ptime; 
        if test>=lim
          cat_io_cprintf('warn',sprintf('"%s" not exist after %d seconds! Proceed! \n',log_name{i},lim));
        end
      else 
        % get PIDs for supervising
        % search for the log entry "CAT parallel processing with MATLAB PID: #####" 
        if job.getPID
          try
            while testpid<limpid
              pause(ptimepid); 
              
              testpid = testpid + ptimepid; 
              FID     = fopen(log_name{i},'r'); 
              txt     = textscan(FID,'%s');
              txt     = txt{1}; 
              PIDi    = find(cellfun('isempty',strfind(txt,'PID:'))==0,1,'first');
              fclose(FID);
              if ~isempty(PIDi)
                PID(i)  = str2double(txt{PIDi+1}); 
                testpid = inf; 
              end
              
              if testpid>=limpid && ~isinf(testpid)
                cat_io_cprintf('warn',sprintf('"%s" no PID information available after %d seconds! Proceed! \n',log_name{i},limpid));
              end
            end
          catch
            cat_io_cprintf('warn',sprintf('No PID information available! Proceed! \n'));
          end
        end
        
        % open file in editor if GUI is available
        test = inf; 
        if ~strcmpi(spm_check_version,'octave') && usejava('jvm') && feature('ShowFigureWindows') && usejava('awt')
          edit(log_name{i});
      end
      end
    end

    % open file in editor if GUI is available
    if ~strcmpi(spm_check_version,'octave') && usejava('jvm') && feature('ShowFigureWindows') && usejava('awt')
      edit(log_name{i});
    end
    if PID(i)>0
      fprintf('\nCheck %s for logging information (PID: ',spm_file(log_name{i},'link','edit(''%s'')')); 
      cat_io_cprintf([1 0 0.5],sprintf('%d',PID(i))); 
    else
      fprintf('\nCheck %s for logging information (',spm_file(log_name{i},'link','edit(''%s'')'));
      cat_io_cprintf([1 0 0.5],'unknown PID'); 
    end
    cat_io_cprintf([0 0 0],sprintf(').\n_______________________________________________________________\n'));

    % starting many large jobs can cause servere MATLAB errors
    pause(1 + rand(1) + job.nproc + numel(job.data)/100);
    jobs(i).data = job.data;
    
    job = jobo; 
  end
  
  job = update_job(job);
  varargout{1} = vout_job(job);
  
  
  
  if job.getPID
    if any(PID==0) 
      cat_io_cprintf('warn',...
        ['\nWARNING: CAT was not able to detect the PIDs of the parallel CAT processes. \n' ...
         '         Please note that no additional modules in the batch can be run \n' ...
         '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
         '         subsequent modules if you split the job into separate processes.\n\n']);
    else
      %% conclusion without filelist
      spm_clf('Interactive'); 
      spm_progress_bar('Init', sum( numel(job_data) ) ,'CAT-Preprocessing','Volumes Complete');      
      
      fprintf('\nStarted %d jobs with the following PIDs:\n',job.nproc);
      for i=1:job.nproc
        fprintf('%3d) %d subjects (PID: ',i,numel(jobs(i).data));
        cat_io_cprintf([1 0 0.5],sprintf('%6d',PID(i))); 
        cat_io_cprintf([0 0 0],sprintf('): ')); 
        cat_io_cprintf([0 0 1],sprintf('%s\n',spm_file(log_name{i},'link','edit(''%s'')')));
      end
      
      
      
      %% supervised pipeline processing 
      %  ------------------------------------------------------------------
      %  This is a "simple" while loop that check if the processes still 
      %  exist and extract information from the log-files, which subject 
      %  was (successfully) processed. 
      %  Finally, a report could be generated and exportet in future that 
      %  e.g. count errors give some suggestions 
      %  ------------------------------------------------------------------
      if job.getPID>1
        cat_io_cprintf('warn',sprintf('\nKilling of this process will not kill the parallel processes!\n'));
        fprintf('_______________________________________________________________\n');
        fprintf('Completed volumes (see catlog files for details!):\n');
        
        % some variables 
        err         = struct('aff',0,'vbm',0,'sbm',0,'else',0,'warn0',0,'warn1',0,'warn2',0); 
        cid         = 0;
        PIDactive   = ones(size(catSID));
        catSIDlast  = zeros(size(catSID));
        [catv,catr] = cat_version;
            
        %% loop as long as data is processed by active tasks
        while cid <= sum( numel(job_data) ) &&  any( PIDactive )
          pause(ptimesid); 
          
          %% get status of each process
          for i=1:job.nproc
            % get FID
            FID = fopen(log_name{i},'r'); 
            txt = textscan(FID,'%s','Delimiter','\n');
            txt = txt{1}; 
            fclose(FID);
            
            % search for the _previous_ start entry "CAT12.# r####: 1/14:   ./MRData/*.nii" 
            catis   = find(cellfun('isempty',strfind(txt,sprintf('%s r%s: ',catv,catr)))==0,2,'last'); 
            catie   = find(cellfun('isempty',strfind(txt,'CAT preprocessing takes'))==0,1,'last');
            if ~isempty(catis) && ( numel(catis)>2 ||  ~isempty(catie) )
              if catis(end)<catie(1)
                cathd = textscan( txt{catis(end)} ,'%s%s%s','Delimiter',':');
              else
                cathd = textscan( txt{catis(1)} ,'%s%s%s','Delimiter',':');
              end
              cathd = textscan( char(cathd{2}) ,'%d','Delimiter','/');
              catSID(i) = cathd{1}(1);
            else 
              catSID(i) = 0; 
            end
            
            % search for the end entry "CAT preprocessing takes ... " to get processing time
            if ~isempty(catie)
              cathd   = textscan( txt{catie} ,'%s%s%s%d%s%s%d%s','Delimiter',' ');
              cattime = [cathd{4}(1) cathd{7}(1)]; 
            else 
              cattime = [0 0];
            end
            
            % search for the end entry "Image Quality Rating (IQR): ... " to get IQR 
            cati = find(cellfun('isempty',strfind(txt,'Image Quality Rating (IQR):'))==0,1,'last');
            if ~isempty(cati) 
              cathd   = textscan( txt{cati} ,'%s%s%s%s%s%s%s','Delimiter',' ');
              catiqr  = [cathd{6} cathd{7}]; 
            else 
              catiqr = {'unknown'};
            end
            
            %% search for preprocessing errors (and differentiate them)
            cati = find(cellfun('isempty',strfind(txt,'CAT Preprocessing error'))==0,1,'last');
            catl = find(cellfun('isempty',strfind(txt,'-----------------------'))==0);
            if ~isempty(cati) && numel(catis)>2 && catie>catis(1) &&  catie<catis(2)
              caterr  = textscan( txt{cati+2} ,'%s','Delimiter','\n');
              caterr  = char(caterr{1});
              caterrcode = ''; 
              for ei = 4:(catl(find(catl>cati,1,'first')+2) - cati)
                catfct{ei-2}  = textscan( txt{cati+ei} ,'%d%s%s','Delimiter',' ');
                if isempty(caterrcode)
                  switch char(catfct{ei-2}{3})
                    % most relevant functions to identify the error 
                    case {'cat_surf_createCS','cat_main','cat_run'}
                      caterrcode = sprintf('%s:%d',char(catfct{ei-2}{3}),double(catfct{ei-2}{1}));
                  end
                end
              end
            else
              caterr     = '';
              caterrcode = ''; 
            end
            %
            % We need some simple error codes that helps the user (check origin)
            % but also us (because they may only send us this line). Hence,
            % the major position of the error (e.g. cat_run/main) is most
            % important.
            %
            % - affreg VBM error > check origin 
            % - R#:cat_run:#:Possible registration error - Check origin! 
            % - R#:cat_main:#:VBM processing error. 
            % - R#:cat_createCS:#:Surface creation error.
            % -----
            % * Handling of warnings?
            %   * yellow light warning vs. orange severe warning  
            %   - low IQR?             < 50%  = yellow warning?
            %   - high topodefarea?    > 1000 = yellow warning? > 5000 = orange warning?
            %   - template corvariance < 0.80 = yellow warning? < 0.60 = orange warning?
            %     this required an update of cat_vol_qa
            % -----

            
            % find out if the current task is still active
            if ispc
              [status,result] = system(sprintf('tasklist /v /fi "PID %d"',PID(i)));  
            else
              [status,result] = system(sprintf('ps %d',PID(i)));
            end
            if isempty( strfind( result , sprintf('%d',PID(i)) ) ) 
              PIDactive(i) = 0; 
            end
            
            
            %% update status
            %  if this tast was not printed before  ( catSIDlast(i) < catSID(i) )  and 
            %  if one subject was successfully or with error processed ( any(cattime>0) || ~isempty(caterr) )
            %fprintf('    %d - %d %d %d %d\n',i,catSIDlast(i), catSID(i), any(cattime>0), ~isempty(caterr) ); 
            if ( catSIDlast(i) < catSID(i) )  &&  ( any(cattime>0) || ~isempty(caterr) )
              cid = cid + 1; 
              catSIDlast(i) = catSID(i);
              
              [pp,ff,ee] = spm_fileparts(jobs(i).data{catSID(i)}); 
              if job.extopts.subfolders
                catlog = fullfile(pp,'report',['catlog_' ff '.txt']); 
              else
                catlog = fullfile(pp,['catlog_' ff '.txt']); 
              end
              if exist(catlog,'file')
                catlogt = ['<a href="matlab:edit(''' catlog ''');">' ...
                  spm_str_manip( catlog , 'k60') ': </a>'];
              else
                catlogt = spm_str_manip( fullfile(pp,[ff ee]), 'k40'); 
              end
              
              
              switch caterr
                case 'Bad SPM-Segmentation. Check image orientation!' % pink error that support user interaction  
                  err.txt   = 'VBM affreg error - Check origin!'; 
                  err.color = [0.9 0 0.9];
                  err.aff   = err.aff + 1;
                case '' % successful processing    
                  % here it would be necessary to differentiate IQR and PQR
                  %if 1 % no warning
                    err.txt   = ''; 
                    err.color = [0 0 0];
                    err.warn0 = err.warn0 + 1; 
                %{
                elseif 0==1 % light yellow warning
                  err.txt   = 'Possible error - Check results!'; 
                  err.color = [0.8 0.6 0];
                  err.warn1 = err.warn1 + 1; 
                elseif 0==2 % severe orange waring 
                  err.txt   = 'Probable error - Check results!'; 
                  err.color = [1 0.3 0];
                  err.warn2 = err.warn2 + 1; 
                end
                %}
                otherwise   
                  err.txt   = errtxt;
                  err.color = [1 0 0];
                  err.vbm   = err.vbm + 1;
              end
              err.txt = sprintf('R%s:%s:%s',catr,caterrcode,err.txt); 
            
              
              % display
              cat_io_cprintf(err.color,sprintf('  %d/%d (job %d: %d/%d): ',...
                cid,sum( numel(job_data) ), i,catSID(i), numel(jobs(i).data) )); 
              cat_io_cprintf([0 0.0 0],catlogt);
              if isempty(caterr)
                cat_io_cprintf(err.color,sprintf(' % 3d.%02d minutes, ',cattime'));
                cat_io_cprintf([0 0.0 0],sprintf('IQR=%s\n',strrep(catiqr{1},'%','%%')));  
              else
                cat_io_cprintf(err.color,sprintf('%s\n',err.txt));  
              end
            end
          end
          spm_progress_bar('Set', cid );
                  
          
        end
      end
    end
    
    %% final report
    fprintf('_______________________________________________________________\n');
    fprintf(['Conclusion: \n' ...
     sprintf('  Processed successfully:% 8d volume(s)\n',err.warn0) ...
     ... sprintf('  Processed with warning:% 8d volume(s)\n',1) ...
     ... sprintf('  Processed with IQR warning:% 8d volume(s)\n',1) ...
     ... sprintf('  Processed with PQR warning:% 8d volume(s)\n',1) ...
     sprintf('  Processed with error:  % 8d volume(s)\n\n',err.aff + err.vbm + err.sbm) ...
             'In case of warnings and errors please check the correct position\n' ...
             'of the AC by using the SPM display function.\n' ...
    ]);   
    fprintf('_______________________________________________________________\n');
    
  else
    cat_io_cprintf('warn',...
      ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
       '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
       '         subsequent modules if you split the job into separate processes.\n\n']);
  end

  spm_progress_bar('Clear');
  return
end

if isfield(job,'printPID') && job.printPID 
  cat_display_matlab_PID
end

job = update_job(job);

varargout{1} = run_job(job);
if isfield(job.extopts,'admin') && isfield(job.extopts.admin,'lazy') && job.extopts.admin.lazy && ...
  ~isfield(job,'process_index') && isfield(job,'nproc') && job.nproc>-1 && (~isfield(job,'process_index'))  
  % set default output even it was not processed this time
  varargout{1} = jobo.vout; 
end

% remove files that do not exist
% EXPLOREASL HACK CREATE DUMMY LOG FILE, SINCE WE DONT USE THE LOGGING OF CAT12
if isfield(varargout{1},'catlog')
    if ~exist(varargout{1}.catlog{1},'file')
        fclose(fopen(varargout{1}.catlog{1}, 'w'));
    end
end

if ~usejava('jvm') % EXPLOREASL HACK TO SKIP THESE WARNINGS WHEN NO JVM
  fprintf('Warning, running CAT12 without JVM, so skipping QC PDF report\n');

  FieldsAre = {'catreport' 'catreportpdf' 'catreportjpg'};
  FilesAre = {'cat_T1.xml' 'catreport_T1.pdf' 'catreportj_T1.jpg'};

  for iFile=1:length(FilesAre)
    if isfield(varargout{1},FieldsAre{iFile})
      if ~exist(varargout{1}.(FieldsAre{iFile}){1},'file')
        fclose(fopen(varargout{1}.(FieldsAre{iFile}){1}, 'w'));
      end
    end
  end
end


varargout{1} = cat_io_checkdepfiles( varargout{1} );
return
%_______________________________________________________________________
function job = update_job(job)

  % set GUI specific parameter if available
  FN = {}; GUIfields = {'registration','segmentation','admin','surface'}; 
  for fnj=1:numel(GUIfields)
    if isfield(job.extopts,GUIfields{fnj})
       FN = [FN;{GUIfields{fnj} fieldnames(job.extopts.(GUIfields{fnj}) )'}];
    end
  end
  for fnj=1:size(FN,1)  
    if isfield(job.extopts,FN{fnj,1})
      for fni=1:numel(FN{fnj,2})
        if isfield(job.extopts.(FN{fnj,1}),FN{fnj,2}{fni})
          job.extopts.(FN{fnj,2}{fni}) = job.extopts.(FN{fnj,1}).(FN{fnj,2}{fni});
        %$else
        %  fprintf('err1: %s\n', FN{fnj,2}{fni});
        end
      end
      job.extopts = rmfield(job.extopts,FN{fnj,1}); % this is just a GUI field! 
    end 
  end
  
  % get defaults
  def = cat_get_defaults;
  
  if isfield(job.extopts,'restypes')
    def.extopts.restype = (char(fieldnames(job.extopts.restypes))); 
    def.extopts.resval  = job.extopts.restypes.(def.extopts.restype);
  end
  
  def.extopts.new_release = 0;
  def.extopts.lazy        = 0;
  def.opts.fwhm           = 1;
  def.nproc               = 0; 
  def.getPID              = 2; % 0 - nothing (old), 1 - get PIDs, 2 - supervise PIDs 
   
  % ROI atlas maps
  if isfield(job.output,'ROImenu') % expert/developer GUI that allows control each atlas map 
    if isfield(job.output.ROImenu,'atlases')
      %% image output
      try, atlases = rmfield(job.output.ROImenu.atlases,'ownatlas'); end
      def.output.atlases = atlases;
      def.output.ROI     = any(cell2mat(struct2cell(atlases))) || ~isempty( job.output.ROImenu.atlases.ownatlas ); 
      
      if ~isempty( job.output.ROImenu.atlases.ownatlas ) && ~isempty( job.output.ROImenu.atlases.ownatlas{1} )
        for i=1:numel( job.output.ROImenu.atlases.ownatlas ) 
          [pp,ff,ee] = spm_fileparts( job.output.ROImenu.atlases.ownatlas{i} ); 
          if any( strcmp( spm_str_manip( def.extopts.atlas( cell2mat(def.extopts.atlas(:,2)) < cat_get_defaults('extopts.expertgui') + 1 ,1) ,'cs') ,ff))
            error('cat_run:ownatlasname', ...
             ['There is a atlas file name conflict. Each atlas name has to be unique. \n' ...
              'Please rename your own atlas map "%s". \n'],fullfile(pp,[ff ee]) ); 
          else
            % add new atlas  
            def.output.atlases.(ff) = 1; 
            def.extopts.atlas = [ def.extopts.atlas; [ {job.output.ROImenu.atlases.ownatlas{i}} {def.extopts.expertgui} {{'gm','wm','csf'}} {0} ] ]; 
          end
        end
      end
    else
      def.output.atlases = struct();
      def.output.ROI     = 0; 
    end
    job = cat_io_checkinopt(job,def);
  end
  % ROI atlas maps
  if isfield(job.output,'sROImenu') % expert/developer GUI that allows control each atlas map 
    if isfield(job.output.sROImenu,'satlases')
      %% image output
      satlases = rmfield(job.output.sROImenu.satlases,'ownatlas'); 
      def.output.satlases = satlases;
      def.output.sROI     = any(cell2mat(struct2cell(satlases))) || ~isempty( job.output.ROImenu.satlases.ownatlas ); 
      
      if ~isempty( job.output.sROImenu.satlases.ownatlas ) && ~isempty( job.output.sROImenu.satlases.ownatlas{1} )
        for i=1:numel( job.output.sROImenu.satlases.ownatlas ) 
          [pp,ff,ee] = spm_fileparts( job.output.sROImenu.satlases.ownatlas{i} ); 
          if any(~cellfun('isempty',strfind( spm_str_manip( def.extopts.satlas(:,1) ,'cs') ,ff)))
            error('cat_run:ownatlasname', ...
             ['There is a surface atlas file name conflict. Each atlas name has to be unique. \n' ...
              'Please rename your own surface atlas map "%s". \n'],fullfile(pp,[ff ee]) ); 
          else
            % add new atlas  
            def.output.satlases.(ff) = 1; 
            def.extopts.satlas = [ def.extopts.satlas; [ {ff} job.output.sROImenu.satlases.ownatlas(i) {def.extopts.expertgui} {0} ] ]; 
          end
        end
      end
    else
      def.output.atlases = struct();
      def.output.sROI    = 0; 
    end
    job = cat_io_checkinopt(job,def);
  else
    def.output.atlases = struct();
    def.output.sROI    = 1; 
    job = cat_io_checkinopt(job,def);
  end
  
  
  if ~isfield(job.output,'atlases') 
    % default GUI that only allow to switch on the settings defined in the default file 
    if ~isfield(job.extopts,'atlas')
      job.extopts.atlas  = def.extopts.atlas;
    end
    
    job.output.atlases   = struct();
    if job.output.ROI 
      % if output, than use the parameter of the default file
      job.output.atlases = cell2struct(job.extopts.atlas(:,4)',spm_str_manip(job.extopts.atlas(:,1),'tr')',2);
      job.output.ROI     = any(cell2mat(struct2cell(job.output.atlases))); 
    end
  end
  if ~isfield(job.output,'satlases') 
    % default GUI that only allow to switch on the settings defined in the default file 
    if ~isfield(job.extopts,'satlas')
      job.extopts.satlas  = def.extopts.satlas;
    end
    
    job.output.satlases   = struct();
    if job.output.sROI 
      % if output, than use the parameter of the default file
      job.output.satlases = cell2struct(job.extopts.satlas(:,4)',spm_str_manip(job.extopts.satlas(:,1),'tr')',2);
      job.output.sROI     = any(cell2mat(struct2cell(job.output.satlases))); 
    end
  end  
  
  if ~isfield(job.output,'atlases')
    if ~isfield(job.extopts,'atlas') && job.output.surface
      job.extopts.satlas = def.extopts.satlas;
    end
  end  
  
  % simplyfied default user GUI input
  if isfield(job.output,'labelnative') 
    job.output.label.native = job.output.labelnative; 
    job.output = rmfield(job.output,'labelnative');
  end

  % simplyfied default user GUI input
  if isfield(job.output,'jacobianwarped') 
    job.output.jacobian.warped = job.output.jacobianwarped; 
    job.output = rmfield(job.output,'jacobianwarped');
  end
  
  % ROI export 
  for ai = 1:size(job.extopts.atlas,1)
    [pp,ff,ee]  = spm_fileparts(job.extopts.atlas{ai,1}); 
    job.extopts.atlas{ai,4} = job.extopts.atlas{ai,2}<=cat_get_defaults('extopts.expertgui') && ...
      exist(job.extopts.atlas{ai,1},'file') && isfield(def.output,'atlases') && isfield(def.output.atlases,ff) && def.output.atlases.(ff);
    % show licence message
    if ~isempty(strfind(ff,'hammers')) && job.extopts.atlas{ai,4}
      disp('--------------------------------------------')
      disp('Free academic end user license agreement for Hammers atlas')
      alert_str = ['For using the Hammers atlas, please fill out license agreement at <a href =',...
      ' "http://brain-development.org/brain-atlases/adult-brain-atlases/adult-brain-maximum-probability-map-hammers-mith-atlas-n30r83-in-mni-space">www.brain-development.org</a>'];
      disp(alert_str);
      disp('--------------------------------------------')
    end
    if ~isempty(strfind(ff,'lpba40')) && job.extopts.atlas{ai,4}
      disp('--------------------------------------------')
      disp('No commercial use of LPBA40 atlas')
      alert_str = ['Permission is granted to use this atlas without charge for non-commercial research purposes only: <a href =',...
      ' "https://www.loni.usc.edu/docs/atlases_methods/Human_Atlas_Methods.pdf">https://www.loni.usc.edu/docs/atlases_methods/Human_Atlas_Methods.pdf</a>'];
      disp(alert_str);
      disp('--------------------------------------------')
    end
  end


  job = cat_io_checkinopt(job,def);
  if ~isfield(job.extopts,'restypes')
    job.extopts.restypes.(def.extopts.restype) = job.extopts.resval;  
  end

  %% handling of SPM biasoptions for specific GUI entry
  if isfield(job.opts,'bias')
    if isfield(job.opts.bias,'spm')
      job.opts.biasstr  = 0; 
      job.opts.biasfwhm = job.opts.bias.spm.biasfwhm; 
      job.opts.biasreg  = job.opts.bias.spm.biasreg; 
    elseif isfield(job.opts.bias,'biasstr')
      job.opts.biasstr  = job.opts.bias.biasstr; 
    end
    job.opts = rmfield(job.opts,'bias'); 
  end
  % the extopts.biasstr controls and overwrites (biasstr>0) the SPM biasreg and biasfwhm parameter
  %   biasstr  = [0.01  0.25  0.50  0.75  1.00] ... result in ?
  %   biasreg  = [0.01  0.0032  0.0010  0.0003  0.0001] ? and ?
  %   biasfwhm = [30 45 60 75 90] for "30 + 60*biasstr? 
  %   biasfwhm = [30.32  42.65  60  84.39 118.71)] for "10^(5/6 + biasstr/3)?  .. allows lower fields 
  if job.opts.biasstr>0 % update biasreg and biasfwhm only if biasreg>0
    % limits only describe the SPM standard range
    job.opts.biasreg	= min(  10 , max(  0 , 10^-(job.opts.biasstr*2 + 2) ));
    job.opts.biasfwhm	= min( inf , max( 30 , 30 + 60*(1-job.opts.biasstr) ));  
  end
  
  % SPM preprocessing accuracy
  if ~isfield(job.opts,'tol')
    job.opts.tol = cat_get_defaults('opts.tol');
  end
  job.opts.tol = min(1e-2,max(1e-6, job.opts.tol));
  
  %% handling of SPM accuracy options for specific GUI entry
  %  Although lower resolution (>3 mm) is not really faster and maybe much 
  %  worse in sense of quality, it is simpler to have a linear decline
  %  rather than describing the other case. 
  %  RD20200130: Takes me a day to figure out that the SPM7771 US failed in 
  %              T1_dargonchow but also single_subjT1 by lower sampl res:
  %                sampval = [3 2.5 2 1.5 1];
  %              Keep in mind that this effects volume resolution (^3), eg
  %              [32 16 8 4 2] .^(1/3) is close to these values. 
  %  RD20200301: However, this setting is really slow and did not solve all
  %              problems, so we go back to previous settings.
  %                sampval =  [5 4 3 2 1]; % that describes volume of 
  %              [125 64 27 8 1] that is also describes the changes in 
  %              processing time roughly 
  sampval               = [5 4 3 2 1]; 
  tolval                = [1e-2 1e-3 1e-4 1e-5 1e-6];
  if isfield(job.opts,'accstr') && ~isfield(job.opts,'acc') 
    job.opts.samp       = sampval( round(job.opts.accstr*4 + 1) );
    job.opts.tol        = tolval(  round(job.opts.accstr*4 + 1) );
  elseif isfield(job.opts,'acc') % developer settings 
    if isfield(job.opts.acc,'accstr')
      job.opts.accstr   = job.opts.acc.accstr; 
      job.opts.samp     = sampval( round(job.opts.acc.accstr*4 + 1));
      job.opts.tol      = tolval(  round(job.opts.acc.accstr*4 + 1));
    elseif isfield(job.opts.acc,'spm')
      job.opts.accstr   = -1; 
      job.opts.samp     = job.opts.acc.spm.samp;
      job.opts.tol      = job.opts.acc.spm.tol;
    end
    job.opts = rmfield(job.opts,'acc'); 
  end
  clear sampval tolval;
  
  %% set Dartel/Shooting templates
  if isfield(job.extopts,'dartel')
    job.extopts.darteltpm   = job.extopts.dartel.darteltpm;
    job.extopts.regstr      = 0; 
  elseif isfield(job.extopts,'shooting')
    job.extopts.shootingtpm = job.extopts.shooting.shootingtpm;
    job.extopts.regstr      = job.extopts.shooting.regstr; 
  end
  
  % find and check the Dartel templates
  [tpp,tff,tee] = spm_fileparts(job.extopts.darteltpm{1});
  job.extopts.darteltpm{1} = fullfile(tpp,[tff,tee]); 
  numpos = min(strfind(tff,'Template_1')) + 8;
  if isempty(numpos)
    error('CAT:cat_main:TemplateNameError', ...
    ['Could not find the string "Template_1" in Dartel template that \n'...
     'indicates the first file of the Dartel template. \n' ...
     'The given filename is "%s.%s" \n'],tff,tee);
  end
  job.extopts.darteltpms = cat_vol_findfiles(tpp,[tff(1:numpos) '*' tff(numpos+2:end) tee],struct('depth',1));
  
  % if we also have found Template_0 we have to remove it from the list
  if numel(job.extopts.darteltpms)==7 
    if ~isempty(strfind(job.extopts.darteltpms{1},'Template_0'))
      for i=1:6, job.extopts.darteltpms{i} = job.extopts.darteltpms{i+1}; end
      job.extopts.darteltpms(7) = [];
    end
  end
  
  job.extopts.darteltpms(cellfun('length',job.extopts.darteltpms)~=length(job.extopts.darteltpm{1}))=[]; % remove to short/long files
  if numel(job.extopts.darteltpms)~=6 && any(job.extopts.regstr==0)
    %%
    files = ''; for di=1:numel(job.extopts.darteltpms), files=sprintf('%s\n  %s',files,job.extopts.darteltpms{di}); end
    error('CAT:cat_main:TemplateFileError', ...
     ['Could not find the expected 6 Dartel template files (Template_1 to Template_6). \n' ...
      'Found %d templates: %s'],numel(job.extopts.darteltpms),files);
  end

  % find and check the Shooting templates
  [tpp,tff,tee] = spm_fileparts(job.extopts.shootingtpm{1});
  job.extopts.shootingtpm{1} = fullfile(tpp,[tff,tee]); 
  numpos = min(strfind(tff,'Template_0')) + 8;
  if isempty(numpos)
    error('CAT:cat_main:TemplateNameError', ...
    ['Could not find the string "Template_0" in Shooting template that \n'...
     'indicates the first file of the Shooting template. \n' ...
     'The given filename is "%s.%s" \n'],tff,tee);
  end
  job.extopts.shootingtpms = cat_vol_findfiles(tpp,[tff(1:numpos) '*' tff(numpos+2:end) tee],struct('depth',1));
  job.extopts.shootingtpms(cellfun('length',job.extopts.shootingtpms)~=length(job.extopts.shootingtpm{1}))=[]; % remove to short/long files
  if numel(job.extopts.shootingtpms)~=5 && any(job.extopts.regstr>0)
    %%
    files = ''; for di=1:numel(job.extopts.shootingtpms), files=sprintf('%s\n  %s',files,job.extopts.shootingtpms{di}); end
    error('CAT:cat_main:TemplateFileError', ...
     ['Could not find the expected 5 Shooting template files (Template_0 to Template_4).\n' ...
      'Found %d templates: %s'],numel(job.extopts.shootingtpms),files);
  end
  
  
  % check range of str variables
  FN = {'WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf'};
  for fni=1:numel(FN)
    if ~isfield(job.extopts,FN{fni})  
      job.extopts.(FN{fni}) = max(0,min(1,job.extopts.(FN{fni})));
    end
  end

  if job.extopts.WMHC<3 && any(cell2mat(struct2cell(job.output.WMH)))
    error('cat_run:bad_WMHC_parameter','Cannot ouput WMH maps if WMHC<3!') 
  end
  if job.extopts.SLC<1 && any(cell2mat(struct2cell(job.output.SL)))
    error('cat_run:bad_SLC_parameter','Cannot ouput stroke lesion maps if SLC is inactive!') 
  end

  
  % deselect ROI output and print warning if ROI output is true and dartel template was changed
  [pth,nam] = spm_fileparts(job.extopts.darteltpm{1});
  if isempty(strfind(nam,'MNI152')) && strcmp(job.extopts.species,'human') && cat_get_defaults('output.ROI')  %~strcmp(nam,'Template_1_IXI555_MNI152')
    warning('DARTEL:template:change',...
      ['Dartel template was changed: Please be aware that ROI analysis \n' ...
       'and other template-specific options cannot be used and ROI \n ' ...
       'output has been deselected.']);
    job.output.ROI = 0;
  end
  
  
  % set boundary box by Template properties if box inf
  if job.extopts.regstr(1)==0
    Vd       = spm_vol([job.extopts.darteltpm{1} ',1']);
  else
    Vd       = spm_vol([job.extopts.shootingtpm{1} ',1']);
  end
  [bb,vox] = spm_get_bbox(Vd, 'old');  
  if bb(1)>bb(2), bbt=bb(1); bb(1)=bb(2); bb(2)=bbt; clear bbt; end
  % Removed BB defintion in GUI and default file in november 2011, because
  % it did not work (need changes in Dartel/Shooting processing) and is not required yet.
  %if job.extopts.bb(1)>job.extopts.bb(2), bbt=job.extopts.bb(1); job.extopts.bb(1)=job.extopts.bb(2); job.extopts.bb(2)=bbt; clear bbt; end
  %job.extopts.bb(isinf(job.extopts.bb))=nan; 
  %job.extopts.bb  = [ min(bb(1,1:3) , job.extopts.bb(1,1:3) ) ; max(bb(2,1:3) , job.extopts.bb(2,1:3) ) ];
  job.extopts.bb = bb; 
  
  job.extopts.vox( isinf(job.extopts.vox) | isnan(job.extopts.vox) ) = []; 
  if isempty( job.extopts.vox ),  job.extopts.vox = cat_get_defaults('extopts.vox'); end 
  job.extopts.vox = abs( job.extopts.vox );
  
  % prepare tissue priors and number of gaussians for all 6 classes
  [pth,nam,ext] = spm_fileparts(job.opts.tpm{1});
  clsn = min(6,numel(spm_vol(fullfile(pth,[nam ext])))); 
  tissue = struct();
  for i=1:clsn
    tissue(i).ngaus = job.opts.ngaus(i);
    tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
  end
  
  tissue(1).warped = [job.output.GM.warped  (job.output.GM.mod==1)        (job.output.GM.mod==2)       ];
  tissue(1).native = [job.output.GM.native  (job.output.GM.dartel==1)     (job.output.GM.dartel==2)    ];
  tissue(2).warped = [job.output.WM.warped  (job.output.WM.mod==1)        (job.output.WM.mod==2)       ];
  tissue(2).native = [job.output.WM.native  (job.output.WM.dartel==1)     (job.output.WM.dartel==2)    ];
  tissue(3).warped = [job.output.CSF.warped (job.output.CSF.mod==1)       (job.output.CSF.mod==2)      ];
  tissue(3).native = [job.output.CSF.native (job.output.CSF.dartel==1)    (job.output.CSF.dartel==2)   ];

  % never write class 4-6
  for i=4:6
    tissue(i).warped = [0 0 0];
    tissue(i).native = [0 0 0];
  end

  job.channel  = struct('vols',{job.data});
  job.tissue   = tissue;

return;

%_______________________________________________________________________
function vout = run_job(job)
  vout   = vout_job(job);

  % load tpm priors 
  fprintf('%s\n','Loading TPM priors');
  tpm = char(cat(1,job.tissue(:).tpm));
  tpm = spm_load_priors8(tpm);

  for subj=1:numel(job.channel(1).vols)
    % __________________________________________________________________
    % Error management with try-catch blocks
    % See also cat_run_newcatch.
    % __________________________________________________________________
    if cat_io_matlabversion>20072 
      cat_run_newcatch(job,tpm,subj); 
    else
      if job.extopts.APP == 1070
        cat_run_job1070(job,tpm,subj); 
      else
        cat_run_job(job,tpm,subj); 
      end
    end
  end

  if usejava('jvm') %% EXPLOREASL HACK TO AVOID ERROR WHEN noJVM
    colormap(gray);
  end
  
  if isfield(job,'nproc') && job.nproc>0 
    fprintf('\n%s',repmat('_',1,72));
    fprintf('\nCAT12 Segmentation job finished.\n');
  end
return
%_______________________________________________________________________

function vout = vout_job(job)
% ----------------------------------------------------------------------
% create output structure for SPM batch mode
% ----------------------------------------------------------------------

n     = numel(job.channel(1).vols);

parts = cell(n,4); % fileparts

biascorr    = {};
wbiascorr   = {};
ibiascorr   = {};
wibiascorr  = {};
ribiascorr  = {};
aibiascorr  = {};
label       = {};
wlabel      = {};
rlabel      = {};
alabel      = {};
catreportjpg= {};
catreportpdf= {};
catreport   = {};
catlog      = {};


%roi         = {};
%fordef      = {};
%invdef      = {};
jacobian    = {};

if job.extopts.subfolders
  roifolder    = 'label';
  surffolder   = 'surf';
  mrifolder    = 'mri';
  reportfolder = 'report';
else
  roifolder    = '';
  surffolder   = '';
  mrifolder    = '';
  reportfolder = '';
end

for j=1:n
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end

% CAT report XML file
% ----------------------------------------------------------------------
roi = cell(n,1);
for j=1:n
    catreport{j}    = fullfile(parts{j,1},reportfolder,['cat_',parts{j,2},'.xml']);
    catreportpdf{j} = fullfile(parts{j,1},reportfolder,['catreport_',parts{j,2},'.pdf']);
    catreportjpg{j} = fullfile(parts{j,1},reportfolder,['catreportj_',parts{j,2},'.jpg']);
end

% 
for j=1:n
    catlog{j} = fullfile(parts{j,1},reportfolder,['catlog_',parts{j,2},'.txt']);
end


% lh/rh/cb central/white/pial/layer4 surface and thickness
% ----------------------------------------------------------------------
surfaceoutput = { % surface texture
  {'central'}                 % no measures - just surfaces
  {}                          % default
  {}                          % expert
  {'pial','white'}            % developer
};
measureoutput = {
  {'thickness'}               % default
  {}                          % no measures
  {}                          % expert
  {'depthWM','depthCSF'}      % developer
};
% no output of intlayer4 or defects in cat_surf_createCS but in cat_surf_createCS2 (but not with fast) 
if isfield(job,'extopts') && isfield(job.extopts,'surface') && ...
   isfield(job.extopts.surface,'collcorr') && job.extopts.surface.collcorr>19 

  surfaceoutput{1} = [surfaceoutput{1},{'pial','white'}];
  surfaceoutput{4} = {}; 
  if any( job.output.surface ~= [ 5 6 ] ) % fast pipeline
    surfaceoutput{3} = {'layer4'}; 
    measureoutput{3} = {'intlayer4','defects'};
  end
end

sides = {'lh','rh'}; 
if any( job.output.surface == [ 2 6 8 ] )
  sides = [sides {'cb'}]; 
end
voutsfields = {};

def.output.surf_measures = 1;
def.extopts.expertgui    = 0;
job = cat_io_checkinopt(job,def); 
% create fields
for si = 1:numel(sides)
  % surfaces
  for soi = 1:numel(surfaceoutput)
    if soi < job.extopts.expertgui + 2
      for soii = 1:numel(surfaceoutput{soi})
        eval( sprintf('%s%s = {};' , sides{si} , surfaceoutput{soi}{soii} ) ); 
        if ~isempty( surfaceoutput{soi} ) && job.output.surface
          eval( sprintf('%s%s = cell(n,1);' , sides{si} , surfaceoutput{soi}{soii} ) ); 
          for j = 1:n
            eval( sprintf('%s%s{j} = fullfile(  parts{j,1} , surffolder , ''%s.%s.%s.gii'' ); ' , ...
              sides{si} , surfaceoutput{soi}{soii} , ...
              sides{si} , surfaceoutput{soi}{soii} , parts{j,2} ) ); 
            voutsfields{end+1} = sprintf('%s%s',  sides{si} , surfaceoutput{soi}{soii} );
          end
        end
      end
    end
  end
  % measures
  for soi = 1:numel(measureoutput)
    if soi < job.extopts.expertgui + 2
      for soii = 1:numel(measureoutput{soi})
        eval( sprintf('%s%s = {};' , sides{si} , measureoutput{soi}{soii} ) ); 
        if ~isempty( measureoutput{soi} ) && job.output.surface
          eval( sprintf('%s%s = cell(n,1);' , sides{si} , measureoutput{soi}{soii} ) ); 
          for j = 1:n
            eval( sprintf('%s%s{j} = fullfile( parts{j,1} , surffolder , ''%s.%s.%s'' ); ' , ...
              sides{si} , measureoutput{soi}{soii} , ...
              sides{si} , measureoutput{soi}{soii} , parts{j,2} ) ); 
            voutsfields{end+1} = sprintf('%s%s',  sides{si} , measureoutput{soi}{soii} );
          end
        end
      end
    end
  end
end


% XML label
% ----------------------------------------------------------------------
if job.output.ROI
    roi = cell(n,1);
    for j=1:n
        roi{j} = fullfile(parts{j,1},roifolder,['catROI_',parts{j,2},'.xml']);
    end
end

% bias
% ----------------------------------------------------------------------
if job.output.bias.native
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},mrifolder,['m',parts{j,2},'.nii']);
    end
end

if job.output.bias.warped
    wbiascorr = cell(n,1);
    for j=1:n
        wbiascorr{j} = fullfile(parts{j,1},mrifolder,['wm',parts{j,2},'.nii']);
    end
end

if job.output.bias.dartel==1
    rbiascorr = cell(n,1);
    for j=1:n
        rbiascorr{j} = fullfile(parts{j,1},mrifolder,['rm',parts{j,2},'_rigid.nii']);
    end
end

if job.output.bias.dartel==2
    abiascorr = cell(n,1);
    for j=1:n
        abiascorr{j} = fullfile(parts{j,1},mrifolder,['rm',parts{j,2},'_affine.nii']);
    end
end

% intensity corrected bias
% ----------------------------------------------------------------------
if job.output.las.native
    ibiascorr = cell(n,1);
    for j=1:n
        ibiascorr{j} = fullfile(parts{j,1},mrifolder,['mi',parts{j,2},'.nii']);
    end
end

if job.output.las.warped
    wibiascorr = cell(n,1);
    for j=1:n
        wibiascorr{j} = fullfile(parts{j,1},mrifolder,['wmi',parts{j,2},'.nii']);
    end
end

if job.output.las.dartel==1
    ribiascorr = cell(n,1);
    for j=1:n
        ribiascorr{j} = fullfile(parts{j,1},mrifolder,['rmi',parts{j,2},'_rigid.nii']);
    end
end

if job.output.las.dartel==2
    aibiascorr = cell(n,1);
    for j=1:n
        aibiascorr{j} = fullfile(parts{j,1},mrifolder,['rmi',parts{j,2},'_affine.nii']);
    end
end


% label
% ----------------------------------------------------------------------
if job.output.label.native
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},mrifolder,['p0',parts{j,2},'.nii']);
    end
end

if job.output.label.warped
    wlabel = cell(n,1);
    for j=1:n
        wlabel{j} = fullfile(parts{j,1},mrifolder,['wp0',parts{j,2},'.nii']);
    end
end

if job.output.label.dartel==1
    rlabel = cell(n,1);
    for j=1:n
        rlabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'_rigid.nii']);
    end
end

if job.output.label.dartel==2
    alabel = cell(n,1);
    for j=1:n
        alabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'_affine.nii']);
    end
end


% tissues
% ----------------------------------------------------------------------
tiss = struct('p',{},'rp',{},'rpa',{},'wp',{},'mwp',{},'m0wp',{});
for i=1:numel(job.tissue)
    if job.tissue(i).native(1)
        tiss(i).p = cell(n,1);
        for j=1:n
            tiss(i).p{j} = fullfile(parts{j,1},mrifolder,['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2)
        tiss(i).rp = cell(n,1);
        for j=1:n
            tiss(i).rp{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'_rigid.nii']);
        end
    end
    if job.tissue(i).native(3)
        tiss(i).rpa = cell(n,1);
        for j=1:n
            tiss(i).rpa{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'_affine.nii']);
        end
    end
    if job.tissue(i).warped(1)
        tiss(i).wp = cell(n,1);
        for j=1:n
            tiss(i).wp{j} = fullfile(parts{j,1},mrifolder,['wp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2)
        tiss(i).mwp = cell(n,1);
        for j=1:n
            tiss(i).mwp{j} = fullfile(parts{j,1},mrifolder,['mwp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(3)
        tiss(i).m0wp = cell(n,1);
        for j=1:n
            tiss(i).m0wp{j} = fullfile(parts{j,1},mrifolder,['m0wp',num2str(i),parts{j,2},'.nii']);
        end
    end
end


% deformation fields
% ----------------------------------------------------------------------
if job.output.warps(1)
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},mrifolder,['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

if job.output.warps(2)
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},mrifolder,['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end


% jacobian
% ----------------------------------------------------------------------
if job.output.jacobian.warped
    jacobian = cell(n,1);
    for j=1:n
        jacobian{j} = fullfile(parts{j,1},mrifolder,['wj_',parts{j,2},'.nii']);
    end
end


% ----------------------------------------------------------------------
vout  = struct('tiss',tiss,'label',{label},'wlabel',{wlabel},'rlabel',{rlabel},'alabel',{alabel},...
               'biascorr',{biascorr},'wbiascorr',{wbiascorr},'roi',{roi},'ibiascorr',{ibiascorr},...
               'wibiascorr',{wibiascorr},'ribiascorr',{ribiascorr},'aibiascorr',{aibiascorr},...
               'invdef',{invdef},'fordef',{fordef},'jacobian',{jacobian},'catreport',{catreport},...
               'catlog',{catlog},'catreportpdf',{catreportpdf},'catreportjpg',{catreportjpg});
             
% add surface fields            
for fi=1:numel(voutsfields)
  eval( sprintf( 'vout.(voutsfields{fi}) = %s;', voutsfields{fi} ) ); 
end

%_______________________________________________________________________
return

%=======================================================================
function [data,err] = remove_already_processed(job,verb)
  if ~exist('verb','var'), verb=0; end
  remove = []; err = zeros(size(job));
  cat_io_cprintf('warn','Lazy processing: \n');
  for subj = 1:numel(job.data)
    [lazy,err(subj)] = checklazy(job,subj,verb); 
    if lazy
      remove = [remove subj];
    end
  end
  cat_io_cprintf('warn','  Skip %d subjects!\n',numel(remove));
  data = job.data(setxor(1:numel(job.data),remove)); 
  cat_io_cprintf([0 0.4 0.6],'\n\nProcess:\n');
  for subj = 1:numel(data)
    cat_io_cprintf([0 0.4 0.6],sprintf(' Code%3d: "%s"\n',err(subj),data{subj}));
  end
  cat_io_cprintf('warn',sprintf('  Process %d subjects!\n',numel(data)));
return
%=======================================================================
function [lazy,FNok] = checklazy(job,subj,verb) %#ok<INUSD>
  if job.extopts.subfolders
    roifolder    = 'label';
    surffolder   = 'surf';
    reportfolder = 'report';
  else
    roifolder    = '';
    surffolder   = '';
    reportfolder = '';
  end

  lazy = 0;
  
  [pp,ff] = spm_fileparts(job.data{subj}); 
  catxml  = fullfile(pp,reportfolder,['cat_' ff '.xml']);
  
  FNok = 0;
  if exist(catxml,'file')

    xml         = cat_io_xml(catxml);
    
    FNopts      = fieldnames(job.opts); 
    FNextopts   = fieldnames(job.extopts);
    FNok        = 1; 
    FNextopts   = setxor(FNextopts,{'LAB','lazy','mrf','atlas','NCstr','resval'});
   
   
    %% check opts
    if isempty(FNopts) || isempty(FNextopts) || ...
       ~isfield(xml.parameter,'opts') || ~isfield(xml.parameter,'extopts')
      return
    end
    for fni=1:numel(FNopts)
      if ~isfield(xml.parameter.opts,FNopts{fni})
        FNok = 2; break
      end
      if ischar(xml.parameter.opts.(FNopts{fni}))
        if ischar(job.opts.(FNopts{fni}))
          if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni}))
            FNok = 3; break
          end
        else
          if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni}){1})
            FNok = 4; break
          end
        end
      else
        if isnumeric(job.opts.(FNopts{fni}))
          if strcmp(FNopts{fni},'ngaus') && numel(xml.parameter.opts.(FNopts{fni}))==4
            % nothing to do (skull-stripped case)
          else
            if xml.parameter.opts.(FNopts{fni}) ~= job.opts.(FNopts{fni})
              FNok = 5; break
            end
          end
        elseif ischar(job.opts.(FNopts{fni}))
          if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni})) 
            FNok = 5; break
          end
        end
      end
    end
    if FNok~=1 % different opts
      return
    end

    %% check extopts
    for fni=1:numel(FNextopts)
      if ~isfield(xml.parameter.extopts,FNextopts{fni})
        FNok = 6; break
      end
      if ischar(xml.parameter.extopts.(FNextopts{fni}))
        if ischar(job.extopts.(FNextopts{fni}))
          if ~strcmp(xml.parameter.extopts.(FNextopts{fni}),job.extopts.(FNextopts{fni}))
            FNok = 7; break
          end
        else
          if ~strcmp(xml.parameter.extopts.(FNextopts{fni}),job.extopts.(FNextopts{fni}){1})
            FNok = 8; break
          end
        end
      elseif iscell(xml.parameter.extopts.(FNextopts{fni}))
        if numel(xml.parameter.extopts.(FNextopts{fni}))~=numel(job.extopts.(FNextopts{fni}))
          FNok = 9; break
        end
        for fnic = 1:numel(xml.parameter.extopts.(FNextopts{fni}))
          if iscell(xml.parameter.extopts.(FNextopts{fni}){fnic})
            for fnicc = 1:numel(xml.parameter.extopts.(FNextopts{fni}){fnic})
              if xml.parameter.extopts.(FNextopts{fni}){fnic}{fnicc} ~= job.extopts.(FNextopts{fni}){fnic}{fnicc}
                FNok = 10; break
              end
            end
            if FNok==10; break; end
          else
            try
              if any(xml.parameter.extopts.(FNextopts{fni}){fnic} ~= job.extopts.(FNextopts{fni}){fnic})
                FNok = 11; break
              end
            catch
                FNok = 11;
            end
            if FNok==11; break; end
          end
          if FNok==11 || FNok==10; break; end
        end
      elseif isstruct(xml.parameter.extopts.(FNextopts{fni}))
        FNX = fieldnames(xml.parameter.extopts.(FNextopts{fni}));
        for fnic = 1:numel(FNX)
          if any(xml.parameter.extopts.(FNextopts{fni}).(FNX{fnic}) ~= job.extopts.(FNextopts{fni}).(FNX{fnic}))
            FNok = 12; break
          end
          if FNok==12; break; end
        end
      else
        % this did not work anymore due to the GUI subfields :/
        %if any(xml.parameter.extopts.(FNextopts{fni}) ~= job.extopts.(FNextopts{fni}))
        %  FNok = 13; break
        %end
      end
    end
    if FNok~=1 % different extopts
      return
    end
    

    % check output
    
    % surface
    if job.output.surface && exist(fullfile(pp,surffolder),'dir')
      Pcentral = cat_vol_findfiles(fullfile(pp,surffolder),['*h.central.' ff '.gii']);
      if  numel(Pcentral)==1 % miss ROI xml
        return
      end
    end
    
    % rois
    if job.output.ROI && ~exist(fullfile(pp,roifolder,['catROI_' ff '.xml']),'file')  % miss ROI xml
      return
    end
      
    %% volumes 
    FNO = fieldnames(job.vout);
    for fnoi = 1:numel(FNO)
      if isempty(job.vout.(FNO{fnoi}))
        continue
      elseif iscell(job.vout.(FNO{fnoi}))
         if ~isempty(job.vout.(FNO{fnoi}){subj}) && ~exist(job.vout.(FNO{fnoi}){subj},'file')
           FNok = 14; break
         end
      elseif isstruct(job.vout.(FNO{fnoi}))
        for si = numel(job.vout.(FNO{fnoi}))
          FNOS = fieldnames(job.vout.(FNO{fnoi})); 
          for fnosi = 1:numel(FNOS)
            if isempty([job.vout.(FNO{fnoi})(si).(FNOS{fnosi})])
              continue
            elseif ~exist(job.vout.(FNO{fnoi})(si).(FNOS{fnosi}){subj},'file')
              FNok = 14; break
            end
          end
        end
      end
      if FNok==14 % miss volumes
        return
      end
    end
    %%
    
    lazy = FNok==1; 
      
  end
 
  if lazy 
    cat_io_cprintf('warn','  "%s" \n',job.data{subj});
  end
return
