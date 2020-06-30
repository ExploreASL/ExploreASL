function varargout = cat_parallelize(job,func,datafield)
% ______________________________________________________________________
% Function to parallelize other functions with job structure, by the 
% following call:
% 
%   SVNid = '$Rev: 1575 $';
% 
%   ... further initialization code
%  
%   % split job and data into separate processes to save computation time
%   if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
%     if nargout==1
%       varargout{1} = cat_parallelize(job,mfilename,'data_surf');
%     else
%       cat_parallelize(job,mfilename,'data_surf');
%     end
%     return
%   elseif isfield(job,'printPID') && job.printPID 
%     cat_display_matlab_PID
%   end 
%  
%   % new banner
%   if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
%   
%   % add system dependent extension to CAT folder
%   if ispc
%     job.CATDir = [job.CATDir '.w32'];
%   elseif ismac
%     job.CATDir = [job.CATDir '.maci64'];
%   elseif isunix
%     job.CATDir = [job.CATDir '.glnx86'];
%   end  
%
%   ... main code
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id: cat_parallelize.m 1575 2020-03-03 11:49:34Z dahnke $

%#ok<*STRIFCND,*STREMP,*STRCL1> % MATLAB contains function 
%#ok<*ASGLU>

  def.verb      = cat_get_defaults('extopts.verb'); 
  def.lazy      = 0; % reprocess exist results
  def.debug     = cat_get_defaults('extopts.verb')>2;
  def.getPID    = 2;
  job.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
 
  job = cat_io_checkinopt(job,def);

  if ~exist('datafield','var'), datafield = 'data'; end

  % rescue original subjects
  job_data = job.(datafield);
  if isstruct(job.(datafield))
    n_subjects = numel(job.(datafield));
  elseif iscell(job.(datafield){1})
    n_subjects = numel(job.(datafield){1});
  else
    n_subjects = numel(job.(datafield));
  end
  if job.nproc > n_subjects
    job.nproc = n_subjects;
  end
  job.process_index = cell(job.nproc,1);
  job.verb = 1; 

  % initial splitting of data
  for i=1:job.nproc
    job.process_index{i} = (1:job.nproc:(n_subjects-job.nproc+1))+(i-1);
  end

  
  % check if all data are covered
  for i=1:rem(n_subjects,job.nproc)
    job.process_index{i} = [job.process_index{i} n_subjects-i+1];
  end

  tmp_array = cell(job.nproc,1);
  logdate   = datestr(now,'YYYYmmdd_HHMMSS');
  PID       = zeros(1,job.nproc);
  %catSID    = zeros(1,job.nproc);
  SID       = cell(1,job.nproc); 
  for i=1:job.nproc
    jobo = job; 
    
    fprintf('Running job %d (with datafield 1):\n',i);
    if isstruct(job.(datafield))
      job.(datafield) = job_data(job.process_index{i});
    elseif iscell(job.(datafield){1})
      for fi=1:numel(job_data{1}(job.process_index{i}))
        fprintf('  %s\n',spm_str_manip(char(job_data{1}(job.process_index{i}(fi))),'a78')); 
      end
      for ix=1:numel(job_data)
        job.(datafield){ix} = job_data{ix}(job.process_index{i});
      end
    else
      for fi=1:numel(job_data(job.process_index{i}))
        fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
      end
      job.(datafield) = job_data(job.process_index{i});
    end

    
    job.verb        = 1; 
    job.printPID    = 1; 
    % temporary name for saving job information
    tmp_name = [tempname '.mat'];
    tmp_array{i} = tmp_name; 
    global defaults cat;  %#ok<TLEV>
    spm12def = defaults;  %#ok<NASGU>
    cat12def = cat;       %#ok<NASGU>
    save(tmp_name,'job','spm12def','cat12def');
    clear spm12def cat12;
 
    % matlab command, cprintferror=1 for simple printing        
    if nargout 
      %%
      matlab_cmd = sprintf(['" ' ...
        'warning off;global cprintferror; cprintferror=1; addpath %s %s %s; load %s; ' ...
        'global defaults; defaults=spm12def; clear defaults; '...
        'global cat; cat=cat12def; clear cat; cat_display_matlab_PID; ' ...
        'output = %s(job); try, save(''%s'',''output''); end; exit"'],... 
        spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),fileparts(which(func)),tmp_name,func,tmp_name);
    else
      matlab_cmd = sprintf(['" ' ...
        'warning off;global cprintferror; cprintferror=1; addpath %s %s %s; load %s; ' ...
        'global defaults; defaults=spm12def; clear defaults; '...
        'global cat; cat=cat12def; clear cat; cat_display_matlab_PID; ' ...
        '%s(job); exit"'],... 
        spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),fileparts(which(func)),tmp_name,func);
    end
    
    % log-file for output
    log_name{i} = ['log_' func '_' logdate '_' sprintf('%02d',i) '.txt']; %#ok<AGROW>
    
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
         %job = update_job(job);
         %varargout{1} = run_job(job);
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
        
        % open file in editor
        test = inf; 
        edit(log_name{i});
      end
    end

    edit(log_name{i});
    if PID(i)>0
      fprintf('\nCheck %s for logging information (PID: ',spm_file(log_name{i},'link','edit(''%s'')')); 
      cat_io_cprintf([1 0 0.5],sprintf('%d',PID(i))); 
    else
      fprintf('\nCheck %s for logging information (',spm_file(log_name{i},'link','edit(''%s'')'));
      cat_io_cprintf([1 0 0.5],'unknown PID'); 
    end
    cat_io_cprintf([0 0 0],').\n_______________________________________________________________\n');

    % starting many large jobs can cause servere MATLAB errors
    pause(1 + rand(1) + job.nproc + numel(job.(datafield))/100);
    jobs(i).(datafield) = job.(datafield); %#ok<AGROW>
   
    job = jobo; 
  end

  
  %job = update_job(job);
  varargout{1} = job; 
  %vout_job(job);

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
        fprintf('%3d) %d subjects (PID: ',i,numel(jobs(i).(datafield)));
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
        fprintf('Process volumes (see catlog files for details!):\n');
        
        % jobIDs .. variable to handle different processing parameters
        %           [ subID , asignedproc , procID , started , finished , error , printed_started , printed_finished ]  
        %           subID, asignedproc, procID .. identification of the subject in job and jobs(i)
        %           started   .. count all appearance of the subject name)
        %           finished  .. set job as finished (e.g. if the next subject name appears)
        %           error     .. not used yet
        %           printed   .. printed on cmd line (do not want to print every loop)
        %
        % A matrix gives a better overview of the status although the fieldnames of a structure are also nice 
        %   jobIDs = struct('subID',0,'jobID',0,'jobsubID',0,'started',0,'finished',0,'error',0,'pstarted',0,'pfinished',0);
        %
        jobIDs      = zeros(sum( numel(job_data) ),9); 
        for pi = 1:job.nproc
          for spi = 1:numel(job.process_index{pi})
            si = job.process_index{pi}(spi);
            jobIDs(si,1) = si; 
            jobIDs(si,2) = pi;
            jobIDs(si,3) = spi;
          end
        end
        
        % some older variables
        %err         = struct('err',0); 
        catSID      = zeros(1,job.nproc);
        catSIDlast  = zeros(size(catSID));
        
        
        % loop as long as data is processed by active tasks
        PIDactive   = ones(size(catSID)); % active jobs
        activePIDs  = 1;                  % loop variable to have an extra loop
        cid         = 0;
        while ( cid <= sum( numel(job_data) ) + 1 ) && activePIDs % plus 1 because only staring jobs are shown!
          pause(ptimesid); 
          if all( PIDactive==0 ), activePIDs = 0; end % if all jobs are finish, we can also stop looping
          %fprintf('--\n');
          
          % get status of each process
          for i=1:job.nproc
            % get FID and read the processing output 
            FID = fopen(log_name{i},'r'); 
            txt = textscan(FID,'%s','Delimiter','\n');
            txt = txt{1}; 
            fclose(FID);
            
            
            % find out if the current task is still active
            if ispc
              [status,result] = system(sprintf('tasklist /v /fi "PID %d"',PID(i))); 
            else
              [status,result] = system(sprintf('ps %d',PID(i)));
            end
            if isempty( strfind( result , sprintf('%d',PID(i)) ) ) 
              PIDactive(i) = 0; 
            end
            
            
            %% search for the start/end entries of a subject, e.g. "CAT12.# r####: 1/14:   ./MRData/*.nii" 
            %  This is the dirty part that is expected to need adaption for 
            %  each new routine that utilize cat_parallelize and has a new
            %  unique job and log structure. 
            %  If this process failed, we have no progress information but
            %  we the batch can still work by asking if the taskIDs exist.
            %
            %   jobs(i) .. the i-th parallel job (0<i<=nproc)
            %   si      .. subject counter in parallel jobs 
            %try
            
              clear findSID SID;
              
              if isstruct( jobs(i).(datafield) )
                % ---------------------------------------------------------
                % if the datafield is a structure, e.g. in cat_long_main 
                % with "jobs(i).subj.mov(1){files}" than we have to run 
                % though a lot of levels "jobs(i).data(si).data2(fni){ci}"
                % ---------------------------------------------------------
                FN = fieldnames( jobs(i).(datafield) ); 
                for fni=1:numel( FN ) % multiple fields
                  for si=1:numel( jobs(i).(datafield) ) % multiple entries 
                    if iscell( jobs(i).(datafield)(si).(FN{fni}) )
                      for ii=1 %:numel( jobs(i).(datafield)(si).(FN{fni}) ) % just count the first job
                        %%
                        [pp,ff,ee] = spm_fileparts( jobs(i).(datafield)(si).(FN{fni}){ii} );
                        
                        % try to find the filename in the log file
                        SID2 = find( cellfun('isempty', strfind( txt , ff ))==0 ); 
                        
                        % if you found the name ...
                        if ~isempty( SID2 )
                          % ... then get the subject number in the original job structure 
                          jobSID = jobIDs(jobIDs(:,2)==i & jobIDs(:,3)==si,1); 
                          
                          jobIDs( jobSID , 4) = numel( SID2 );  % count the number of filename output 
                          if si>1
                            jobSIDp = jobIDs(jobIDs(:,2)==i & jobIDs(:,3)==(si-1),1); 
                            jobIDs( jobSIDp , 5) = 1; % if work on a new subject that mark the previous one as finished 
                          end
                        
                          if strcmp(FN{fni},'mov') % long
                            % this does not work
                            FIN = find( cellfun('isempty', strfind( txt, 'Finished CAT12 longitudinal processing of '))); 
                            if ~isempty( FIN ) && numel(txt)>FIN(end)
                              SIDFIN = find( cellfun('isempty', strfind( txt( FIN(end)+1 ) , ff ))==0 ); 
                              if ~isempty( SIDFIN )
                                jobIDs( jobSID , 5) = 1;
                                jobIDs( jobSID , 6) = 0;
                              else
                                jobIDs( jobSID , 5) = 1;
                                jobIDs( jobSID , 6) = 1;
                              end
                            end
                          end

                          % if you found something and did not print it before than print that you have started this subject 
                          if jobIDs(jobSID,4)>0 && jobIDs(jobSID,7)==0 && jobIDs(jobSID,8)==0  
                            jobIDs(jobSID,7) = 1;
                            cat_io_cprintf([0 0 0.5],sprintf('  started  %d/%d (pjob %d: %d/%d): %s\n',...
                              jobSID,sum( numel(job_data) ),  i, jobIDs(jobSID,3) , numel(jobs(i).(datafield)), ...
                              spm_str_manip( jobs(i).(datafield)(si).(FN{1}){1}, 'k40') ));
                          end
                          
                          % if you started the next subject or the job is finished then print that you finished the job
                          if ( jobIDs(jobSID,5)>0 || PIDactive(i)==0 ) && jobIDs(jobSID,8)==0
                            jobIDs(jobSID,8) = 1;
                            if jobIDs( jobSID , 6) == 0
                              cat_io_cprintf( [0 0.5 0] ,sprintf('  finished %d/%d (pjob %d: %d/%d): %s\n',...
                                jobSID,sum( numel(job_data) ),  i, jobIDs(jobSID,3) , numel(jobs(i).(datafield)), ...
                                spm_str_manip( jobs(i).(datafield)(si).(FN{1}){1}, 'k40') ));
                            else
                              cat_io_cprintf( [1 0 0] ,sprintf('  failed %d/%d (pjob %d: %d/%d): %s\n',...
                                jobSID,sum( numel(job_data) ),  i, jobIDs(jobSID,3) , numel(jobs(i).(datafield)), ...
                                spm_str_manip( jobs(i).(datafield)(si).(FN{1}){1}, 'k40') ));
                            end
                            cid = cid + 1; 
                          end
                          
                          % if you found something and did not print it before than print that you have started this subject 
                          if jobIDs(jobSID,4)>0 && jobIDs(jobSID,7)==0
                            jobIDs(jobSID,7) = 1;
                            cat_io_cprintf([0 0 0.5],sprintf('  started  %d/%d (pjob %d: %d/%d): %s\n',...
                              jobSID,sum( numel(job_data) ),  i, jobIDs(jobSID,3) , numel(jobs(i).(datafield)), ...
                              spm_str_manip( jobs(i).(datafield)(si).(FN{1}){1}, 'k40') ));
                          end
                          
                        end
                      end
                    end
                  end
                end
              elseif strcmp( datafield , 'data_surf' )
                % ---------------------------------------------------------
                % surfaces ... here the filenames of the processed data
                % change strongly due to side coding ...
                % ---------------------------------------------------------
                for si=1:numel( jobs(i).(datafield) )
                  if iscell( jobs(i).(datafield){si} )
                    for sii=1:numel( jobs(i).(datafield){si} )
                      [pp,ff,ee] = spm_fileparts(jobs(i).(datafield){si}{sii}); 

                      SID{si} = ...
                        find(cellfun('isempty', strfind( txt , pp ))==0,1,'first') & ... 
                        find(cellfun('isempty', strfind( txt , ff ))==0,1,'first') & ...
                        find(cellfun('isempty', strfind( txt , ee ))==0,1,'first'); 
                    end
                  else
                    [pp,ff,ee] = spm_fileparts(jobs(i).(datafield){si}); 

                    SID{si} = ...
                      find(cellfun('isempty', strfind( txt , pp ))==0,1,'first') & ... 
                      find(cellfun('isempty', strfind( txt , ff ))==0,1,'first') & ...
                      find(cellfun('isempty', strfind( txt , ee ))==0,1,'first'); 
                  end
                end
              else % volumes
                for si=1:numel( jobs(i).(datafield) )
                  SID{si} = find(cellfun('isempty', strfind( txt , jobs(i).(datafield){si} ))==0,1,'first');
                end
              end
               
            %{
            catch
              if ~exist('noSID','var') || noSID==0
                noSID = 1; 
                cat_io_cprintf('warn','  Progress bar did not work but still monitoring the tasks.\n'); 
              end
            end
            %}
              
            
            
            %% update status (old version!)
            %  if this task was not printed before  ( catSIDlast(i) < catSID(i) )  and 
            %  if one subject was successfully or with error processed ( any(cattime>0) || ~isempty(caterr) )
            if ~isstruct( jobs(i).(datafield) )
              %%
              findSID = find(cellfun('isempty',SID)==0,1,'last'); 
              if ~isempty(findSID), catSID(i) = findSID; end
              
              if numel(catSID)>1 && ( catSIDlast(i) < catSID(i) ) 
                try
                  catSIDlast(i) = catSID(i);
                  cat_io_cprintf([ 0 0 0 ],sprintf('  %d/%d (pjob %d: %d/%d): %s\n',...
                    sum(catSID) ,sum( numel(job_data) ),  i,  catSID(i), numel(jobs(i).(datafield)), ...
                    spm_str_manip( jobs(i).(datafield){catSID(i)} , 'k40') )); 
                  cid = cid + 1; 
                end
              end
            end
          end
          
          % further error handling of different functions
          % 'Error using MATLABbatch system'
          
          
          %
          spm_progress_bar('Set', cid );
                  
          
        end
        
      end
      %fprintf('done\n');
    end
    
    
    %% Merge output of the subprocesses to support SPM dependencies.
    %  --------------------------------------------------------------------
    %  The results of a function are saved as variable named "output" in 
    %  the temporary matlab files that were also used for data input.
    %  However, currently only the first output can be used and a structure
    %  array is expected that have maximal two sublevels, e.g.
    %    output.data = {file1, ...}
    %    output.subject.data = {file1, ...}
    %  --------------------------------------------------------------------
    if nargout>0
      %% load the results of the first job
      load(tmp_array{1});
      
      % add further output results of the other jobs 
      for oi=1:numel(tmp_array)
        clear output
        load(tmp_array{oi}); 

        if exist('output','var')
          if ~exist('varargout','var') 
            % here we can simply set the output
            varargout{1} = output; 
          else
            % here we have to add the data to each field 
            % (cat_struct does therefore not work)
            FN = fieldnames(output);
            for fni=1:numel(FN)
              if ~isfield( varargout{1} , FN{fni} ) 
                % no subfield ( = new subfield ) > just add it
                varargout{1}.(FN{fni}) = output.(FN{fni}); 
              else
                % existing subfield > merge it
                if ~isstruct( output.(FN{fni}) ) 
                  if isfield(varargout{1},FN{fni}) && ~isstruct( varargout{1}.(FN{fni}) )
                    if size(varargout{1}.(FN{fni}),1) > size(varargout{1}.(FN{fni}),2) || ...
                      size(output.(FN{fni}),1) > size(output.(FN{fni}),2) 
                      varargout{1}.(FN{fni}) = [varargout{1}.(FN{fni}); output.(FN{fni})]; 
                    else
                      varargout{1}.(FN{fni}) = [varargout{1}.(FN{fni}), output.(FN{fni})]; 
                    end

                    % cleanup (remove empty entries of failed processings)
                    if iscell(varargout{1}.(FN{fni}))
                      for ffni=numel( varargout{1}.(FN{fni}) ):-1:1
                        if isempty(  varargout{1}.(FN{fni}){ffni} )
                          varargout{1}.(FN{fni})(ffni) = []; 
                        end
                      end
                    end
                  else
                    varargout{1}.(FN{fni}) = output.(FN{fni}); 
                  end
                else
                  % this is similar to the first level ...
                  FN2 = fieldnames(output.(FN{fni}));
                  for fni2 = 1:numel(FN2)
                    for fnj=1:numel( output.(FN{fni}) ) 
                      if ~isstruct( output.(FN{fni})(fnj).(FN2{fni2}) )
                        if isfield(varargout{1}.(FN{fni})(fnj),FN2{fni2}) && ~isstruct( varargout{1}.(FN{fni})(fnj).(FN2{fni2}) )
                          if size(varargout{1}.(FN{fni})(fnj).(FN2{fni2}),1) > size(varargout{1}.(FN{fni})(fnj).(FN2{fni2}),2) || ...
                             size(output.(FN{fni})(fnj).(FN2{fni2}),1) > size(output.(FN{fni})(fnj).(FN2{fni2}),2)
                            varargout{1}.(FN{fni})(fnj).(FN2{fni2}) = [varargout{1}.(FN{fni})(fnj).(FN2{fni2}); output.(FN{fni})(fnj).(FN2{fni2})]; 
                          else
                            varargout{1}.(FN{fni})(fnj).(FN2{fni2}) = [varargout{1}.(FN{fni})(fnj).(FN2{fni2}), output.(FN{fni})(fnj).(FN2{fni2})]; 
                          end

                          % cleanup (remove empty entries of failed processings)
                          if iscell(varargout{1}.(FN{fni})(fnj).(FN2{fni2}))
                            for ffni=numel( varargout{1}.(FN{fni})(fnj).(FN2{fni2}) ):-1:1
                              if isempty(  varargout{1}.(FN{fni})(fnj).(FN2{fni2}){ffni} )
                                varargout{1}.(FN{fni})(fnj).(FN2{fni2})(ffni) = []; 
                              end
                            end
                          end
                        else
                          varargout{1}.(FN{fni})(fnj).(FN2{fni}) = output.(FN{fni})(fnj).(FN2{fni}); 
                        end
                      else
                        error('Only 2 level in output structure supported.')
                      end
                    end
                  end
                end
              end
            end
          end
        else
          % create an error message?
          oistr = cat_io_strrep( sprintf('%dth',oi) , {'1th'; '2th'; '3th'} , {'1st','2nd','3rd'} ); 
          cat_io_cprintf('error',sprintf(...
            'The %s processes does not contain the output variable for depending jobs. Check log-file for errors: \n  ', oistr));
          %fprintf([spm_file(log_name{oi},'link',sprintf('edit(%s)',log_name{oi})) '\n\n']);
          cat_io_cprintf([0 0 1],sprintf('%s\n\n',spm_file(log_name{i},'link','edit(''%s'')')));
        end
      end
    end
    
    % no final report yet ...
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