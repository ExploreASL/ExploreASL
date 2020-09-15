function cat_run_newcatch1585(job,tpm,subj)
% ______________________________________________________________________
% This function contains the new matlab try-catch block.
% The new try-catch block has to be in a separate file to avoid an error.
%
% See also cat_run_newcatch.
% ______________________________________________________________________
% $Revision: 1563 $  $Date: 2020-02-11 17:00:46 +0100 (Tue, 11 Feb 2020) $

  global cat_err_res;

  [pth,nam,ext] = spm_fileparts(job.channel(1).vols{subj}); 

  try
    cat_run_job1585(job,tpm,subj); % the cat_run_job1070 is only called by older functions
  catch caterr 
    %% add further information for special errors
    if isempty(caterr.identifier)
      switch caterr.message
        case 'insufficient image overlap'
          adderr = MException('SPM:AlignmentError','There is not enough overlap in the images to obtain a solution.');
        otherwise
          adderr = MException('SPM:CAT:cat_main',strrep(caterr.message,'\','\\'));
      end
      caterr = addCause(caterr,adderr);
    end
    
    if job.extopts.subfolders
      mrifolder = 'mri';
    else
      mrifolder = '';
    end

    cat_io_cprintf('err',sprintf('\n%s\nCAT Preprocessing error for %s:\n%s\n%s\n%s\n', ...
      repmat('-',1,72),nam,repmat('-',1,72),caterr.message,repmat('-',1,72)));  
    
    
    %% check for filenames that are usually indicated by '"'
    testmessage = 0; % just for tests 
    if testmessage
     caterr.identifier =  'CAT:error0815:BadFileInput'; 
     caterr.message    = ['Bad values variable "xyz" in "C:\private\patient\name" ' ...
       'and \n also in the second line "/private2/patient/name".\n ' ...
       'But also without /private3/pat3/nam3 and C:\blub\bla. ' ...
       'Don''t forget special character for    area (33 mm' char(178) ...
       ') or volume (3 mm' char(179) ') or 33' char(177) '0.33%.\n' ...
       'Ignore/Replace unclear chars :' char(200:210)];  
    end
    
    % anonymize by removing filename within '"' characters were found
    ind_str  = strfind(caterr.message,'"');
    ind_str2 = [strfind(caterr.message,'/') strfind(caterr.message,'\') ...
      strfind(caterr.message,'.nii') strfind(caterr.message,'.gii')];
    
    caterr_id = caterr.identifier;
    caterr_message_str = caterr.message;
    if mod(length(ind_str),2) == 0
      for i = length(ind_str):-2:1
        if any( ind_str2>ind_str(i-1) & ind_str2<ind_str(i) ) % replace only files "a/b" "c:\c" but not variables "var" 
          caterr_message_str = [ caterr_message_str(1:ind_str(i-1)) 'FILE' caterr_message_str(ind_str(i):end) ];
        end
      end
    end
    
    % again check for filenames indicated by slashes/backslashes ... , 
    caterr_id          = cat_io_strrep(caterr_id         ,{'\\','\n','\t','\a','\f','\r','\v'},' ');  
    caterr_message_str = cat_io_strrep(caterr_message_str,{'\\','\n','\t','\a','\f','\r','\v'},' ');  
    caterr_message_str = cat_io_strrep(caterr_message_str,{'   '},' ');  
    caterr_message_str = cat_io_strrep(caterr_message_str,{'  '},' ');  
    
    % replace special character may used in some messages
    rc = { {char(177) char(178) char(179) } {'+-' '2' '3'} };
    caterr_id           = char( cat_io_strrep(caterr_id         , rc{1} , rc{2}) ); 
    caterr_message_str  = char( cat_io_strrep(caterr_message_str, rc{1} , rc{2}) ); 
    
    
    % further anonymize by removing filename if slashes/backslashes were found
    % seperate the string into words, find path-strings and replace them 
    words = textscan(caterr_message_str,'%s'); words = words{1}; 
    files = find( ~cellfun('isempty',strfind(words,'\')) | ~cellfun('isempty',strfind(words,'/')) == true);
    for wi = numel(files):-1:1
      if any(words{files(wi)}(end) == '.,;'), words{files(wi)}(end) = []; end  
      caterr_message_str = char( strrep( caterr_message_str , words(files(wi)) , '"FILE"') ); 
    end
    clear words files
    
    % finally remove/replace bad characters
    whitelist           = [ char(48:57) char(65:89) char(97:122) ' +-\/_;:."<=>|()&!?%?`''?^[]'];
    caterr_idi          = false( size( caterr_id ) );
    caterr_message_stri = false( size( caterr_message_str) ); 
    for wi=1:numel(whitelist)
      caterr_idi( caterr_id == whitelist(wi) ) = true; 
      caterr_message_stri( caterr_message_str == whitelist(wi) ) = true; 
    end
    if 0 % remove 
      caterr_id          = caterr_id(caterr_idi);
      caterr_message_str = caterr_message_str(caterr_message_stri);
    else % replace ... maybe better do avoid empty messages
      caterr_id(~caterr_idi) = 'X';
      caterr_message_str(~caterr_message_stri) = 'X';
    end
    
    % Or a private/general part that is useful for users but not for us. 
    
    
    % We may can use a general limitation for the error message?
    maxlength = 256;
    if numel(caterr_message_str)>maxlength
      caterr_message_str = [spm_str_manip(caterr_message_str,sprintf('f%d',maxlength)) ' ...']; 
    end
    
    if testmessage
      disp(caterr_message_str)
    end
    
    %% remove uninteresting messages
    ignore_message = 0; 
    keywords = {
      ... possible orientation and resolution errors issues
      'insufficient image overlap'
      'Image does not have 3 dimensions.'
      'Voxel resolution has to be better than 5 mm'
      'Out of memory.' 
      'cat_run_job:restype'
      ... file reading/writing messages
      '** failed to open'
      'Access is denied.'
      'Cant open file'
      'Cant create file mapping. ' 
      'Cannot create output file '
      'cp: cannot create regular file'
      'File too small'
      'Invalid file identifier.'
      'Permission denied'
      'The process cannot access the file '
      'There was a problem while generating '
      'Unable to write file'
      };
    for ki = 1:numel(keywords)
      if strfind( caterr_message_str , keywords{ki} ) 
        ignore_message = 1; 
      end
    end

    %% send error information, CAT12 version and computer system
    if cat_get_defaults1585('extopts.send_info') && ~ignore_message && job.extopts.expertgui<2
      [v,rev] = cat_version; expertguistr = ' ed';
      str_err = sprintf('%s%s|',rev,deblank(expertguistr(job.extopts.expertgui + 1))); % revision and guilevel
      for si=1:numel(caterr.stack)
        str_err = [str_err '|' caterr.stack(si).name ':' num2str(caterr.stack(si).line)];
      end      
      str_err = str_err(2:end); % remove first "|"
      urlinfo = sprintf('%s%s%s%s%s%s%s%s%s%s',cat_version,'%2F',computer,'%2F','errors','%2F',caterr_id,'%2F',caterr_message_str,str_err);
      cat_io_send_to_server(urlinfo);
    end

    % write error report
    caterrtxt = cell(numel(caterr.stack)+2,1);
    caterrtxt{1} = sprintf('%s\n',caterr.identifier);
    caterrtxt{2} = sprintf('%s\n',caterr.message); 
    for si=1:numel(caterr.stack)
      cat_io_cprintf('err',sprintf('% 5d - %s\n',caterr.stack(si).line,caterr.stack(si).name));  
      caterrtxt{si+2} = sprintf('% 5d - %s\n',caterr.stack(si).line,caterr.stack(si).name); 
    end
    cat_io_cprintf('err',sprintf('%s\n',repmat('-',1,72)));  

    % save cat xml file
    caterrstruct = struct();
    for si=1:numel(caterr.stack)
      caterrstruct(si).line = caterr.stack(si).line;
      caterrstruct(si).name = caterr.stack(si).name;  
      caterrstruct(si).file = caterr.stack(si).file;  
    end
    
    % better to have the res that the opt field
    if isfield(cat_err_res,'res')
      job.SPM.res = cat_err_res.res;
    elseif isfield(cat_err_res,'obj')
      job.SPM.opt = cat_err_res.obj;
    end
    
    qa = cat_vol_qa('cat12err',struct('write_csv',0,'write_xml',1,'caterrtxt',{caterrtxt},'caterr',caterrstruct,'job',job,'subj',subj));
    cat_io_report(job,qa,subj)
    
    % delete noise corrected image
    if exist(fullfile(pth,mrifolder,['n' nam ext]),'file')
      try %#ok<TRYNC>
        delete(fullfile(pth,mrifolder,['n' nam ext]));
      end
    end
    
    if job.extopts.subfolders
      reportfolder = 'report';
    else
      reportfolder = '';
    end
    % create an error directory with errortype subdirectory for all failed datasets
    % copy the cat*.xml and catreport_*pdf 
    % create a symbolic link of the original file
    if job.extopts.subfolders
      %%
      errfolder    = 'err';
      [ppe,ffe]    = spm_fileparts(caterr.stack(1).file); 
      suberrfolder = sprintf('%s.line%d.%s',ffe,caterr.stack(1).line,caterr.identifier); 
      suberrfolder = char(regexp(strrep(suberrfolder,':','.'),'[A-Za-z0-9_.\- ]','match'))'; % remove bad chars
      suberrfolder = strrep(suberrfolder,' ','_');
      if ~exist(fullfile(pth,errfolder,suberrfolder),'dir'), mkdir(fullfile(pth,errfolder,suberrfolder)); end
      catfile = fullfile(pth,reportfolder,['cat_' nam '.xml']);
      logfile = fullfile(pth,reportfolder,['catlog_' nam '.txt']);
      repfile = fullfile(pth,reportfolder,['catreport_' nam '.pdf']);
      if exist(catfile,'file'), copyfile(catfile,fullfile(pth,errfolder,suberrfolder)); end
      if exist(logfile,'file'), copyfile(catfile,fullfile(pth,errfolder,suberrfolder)); end
      if exist(repfile,'file'), copyfile(repfile,fullfile(pth,errfolder,suberrfolder)); end
      if ismac || isunix
        [ST, RS] = system(sprintf('ln -s -F "%s" "%s"',...
          fullfile(pth,[nam ext]),fullfile(pth,errfolder,suberrfolder,[nam ext])));
          cat_check_system_output(ST,RS,job.extopts.verb>2);
      end  
    end
    
    %% rethrow error 
    if ~job.extopts.ignoreErrors
      rethrow(caterr); 
    end 
  end
end
