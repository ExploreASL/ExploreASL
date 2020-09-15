function varargout = cat_get_defaults1173(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = cat_get_defaults1173(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT cat_get_defaults1173(defstr, defval)
% Sets the cat value associated with identifier "defstr". The new
% cat value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit cat_defaults.m.
%
% FORMAT cat_get_defaults1173(defstr, 'rmfield')
% Removes last field of defstr eg. defstr = 'opts.sopt.myfield' will remove 
% 'myfield'. 
%
% FORMAT cat_get_defaults1173(defstr, 'rmentry')
% Removes last field of defstr eg. defstr = 'opts.sopt.myfield' will remove 
% 'sopt' with all all subfield. 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% based on Volkmar Glauches version of
% spm_get_defaults
% $Id: cat_get_defaults1173.m 1389 2018-11-11 10:39:41Z dahnke $

global cat1173;
if isempty(cat1173) || ~isfield(cat1173,'version') || cat1173.version~=1173;
  cat_defaults1173;
end

if nargin == 0
    varargout{1} = cat1173;
    return
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    % default output
    try
      varargout{1} = subsref(cat1173, subs);
    catch
      varargout{1} = []; 
    end
    return;
elseif nargin == 2
    switch varargin{1}
        case 'rmfield'
          % remove the last field of the given defstr
            mainfield = tags{1}{1}; 
            for ti=2:numel(tags{1})-1
                mainfield = [mainfield '.' tags{1}{ti}]; %#ok<AGROW>
            end
            subfield  = tags{1}{end};  
            %fprintf('Remove field "%s" in "cat.%s"!\n',subfield,mainfield);
            eval(sprintf('cat.%s = rmfield(cat.%s,subfield);',mainfield,mainfield));
        case 'rmentry'
          % removes the complete entry of the given defstr
            %fprintf('Remove entry "%s" "cat"!\n',tags{1}{1});
            cat1173 = rmfield(cat1173,defstr); 
        otherwise
          % add an new entry
            cat1173 = subsasgn(cat1173, subs, varargin{1});
    end
end
if nargout == 1
  % output in case changes in cat
    varargout{1} = cat1173;
end