function varargout = ps_LST_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = cg_vbm_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT cg_vbm_get_defaults(defstr, defval)
% Sets the vbm value associated with identifier "defstr". The new
% vbm value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit cg_vbm_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% based on Volkmar Glauches version of
% spm_get_defaults
% $Id: cg_vbm_get_defaults.m 477 2013-04-22 09:08:41Z gaser $

global lst_obj;
if isempty(lst_obj)
    ps_LST_defaults;
end

if nargin == 0
    varargout{1} = lst_obj;
    return
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(lst_obj, subs);
else
    lst_obj = subsasgn(lst_obj, subs, varargin{1});
end
