function [S,Stype,removed] = cat_io_checkdepfiles(S,usedummy)
% _________________________________________________________________________
% Remove non-existing files from the (SPM dependency) variable S. 
% If there is only one file than create a dummy file to avoid batch errors.
% However, this may cause larger problems and is false by default and print
% an error message (just print a message, no real error). 
% Moreover, it is possible that removing an entry can cause problems in 
% related datasets that count on a specific file order, so a message is 
% printed even in this case  ...
%
%   S = cat_io_checkdepfiles(S);
%
%   [S,Stype] = cat_io_checkdepfiles(S,usedummy)
% _________________________________________________________________________
% Robert Dahnke
% $Id: cat_io_checkdepfiles.m 1614 2020-05-08 14:39:49Z gaser $


  if ~exist('usedummy','var')
    usedummy = 0;
  end

  Stype   = '';
  removed = 0;
  
  if ischar(S)
    for i = 1:size(S,1)
      if ~exist(S(i,:),'file')
        % Here, I remove just one element that is missing (e.g. one failed
        % preprocessing). However, a following job (e.g. smoothing) does
        % not count on it.
        [pp,ff,ee] = spm_fileparts(S(i,:));
        S(i,:)  = '';
        Stype   = ee;
        removed = 1; 
        % prevent interpreting backslash as escape character
        removeFile = strrep(fullfile(pp,[ff ee]), '\', '\\');
        cat_io_cprintf('warn',sprintf('  Remove "%s" from dependency list because it does not exist!\n',removeFile));
      end
    end
  elseif iscell(S)
    for i = 1:numel(S)
      S{i} = cat_io_checkdepfiles( S{i} , usedummy );
    end
  elseif isstruct(S)
    %%
    FN = fieldnames(S);
    for i = 1:numel(FN)
      SFN = S.(FN{i}); 
      for ii = 1:numel( S.(FN{i}) )
        SFNi    = SFN(ii);
        [SFNi,SFNtype,removed] = cat_io_checkdepfiles( SFNi ,usedummy );
        SFN(ii) = SFNi; 
        clear SFNi;
      end    
      if isstruct(SFN)==0 
       if numel(SFN)==1 && size(SFN{1},1)==0 && ~isequal(SFN,S.(FN{i}))
        if usedummy
          SFN = {create_dummy_volume(SFNtype)};
        else
          % Here, the full entry of a class of files becomes empty and this 
          % of cause will trouble an following job (empty imput).
          % I give not further file information because these were given by
          % the upper warning!
          cat_io_cprintf('err',sprintf(['One or multiple files do not exist and were removed from the dependency list \n' ...
            'and following batches will may not work correctly!\n\n'])); 
          SFN = {};
        end
       elseif removed
        cat_io_cprintf('warn',sprintf([
           'One or multiple files do not exist and were removed from the dependency list' \n ...
           'and following batches those file input number and order is relevant may do not work properly!\n\n'])); 
       end
      end
      S.(FN{i}) = SFN; 
      
      clear SFN;      
    end
  end
end
function Pdummy = create_dummy_volume(type)
  switch type
    case '.nii'
      Pvol   = fullfile(spm('dir'),'toolbox','cat12','templates_volumes','cat.nii');
      Pdummy = fullfile(spm('dir'),'toolbox','cat12','cattest','batchdummy.nii');
      if ~exist( fileparts(Pdummy) , 'dir')
        mkdir( fileparts(Pdummy) ); 
      end
      copyfile(Pvol,Pdummy);
    case {'.xml','.csv','.txt',''}
      Pvol   = fullfile(spm('dir'),'toolbox','cat12','templates_volumes','mori.csv');
      Pdummy = fullfile(spm('dir'),'toolbox','cat12','cattest',['batchdummy' type]);
      if ~exist( fileparts(Pdummy) , 'dir')
        mkdir( fileparts(Pdummy) ); 
      end
      copyfile(Pvol,Pdummy);
  end
end