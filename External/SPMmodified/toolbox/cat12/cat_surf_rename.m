function [PO,sinfo] = cat_surf_rename(P,varargin)
% ______________________________________________________________________
% Rename parts of a cat12 surface filename.
%
%   [PO,sinfo] = cat_surf_rename(P,varargin)
% 
%   P   = 'lh.central.test.gii';
%   Pth = cat_surf_rename(P,'dataname','s3tickness');
%
%   Check the help for cat_surf_info for information about fields that 
%   can be renamed.
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_surf_rename.m 1590 2020-03-23 07:38:58Z dahnke $

% Todo: update sinfo! also for other fields!


  PN = struct();
  if mod(nargin-1,2)==1
    error('paired input');
  else
    for i=1:2:numel(varargin)
      PN.(varargin{i}) = varargin{i+1}; 
    end
  end
  
  if ~isstruct(P)
    sinfo = cat_surf_info(P);
  else
    sinfo = P;
  end  
  
  
  if nargin>0
    FN = fieldnames(PN);
    PO = cell(size(sinfo));
    for i=1:numel(sinfo)

      if ~isempty(PN)
        for fni=1:numel(FN)
          sinfo(i).(FN{fni}) = PN.(FN{fni}); 
        end
      end

      if any(~cellfun('isempty',strfind(FN,'templateresampled')))
        sinfo(i).resampled = 1;
        templateresampled  = PN.templateresampled; 
      else
        if sinfo(i).resampled==1
          if sinfo(i).template==1
            templateresampled=''; %.template';
          else
            templateresampled='.resampled';
          end
        else
          templateresampled='';
        end
      end

      if isempty(sinfo(i).name), namedot=''; else namedot='.'; end
      if isempty(sinfo(i).side), sidedot=''; else sidedot='.'; end
      if isempty(templateresampled), tempdot=''; else tempdot='.'; end

      PO{i} = fullfile(sinfo(i).pp,sprintf('%s%s%s%s%s%s%s%s',...
        sinfo(i).preside,...
        sinfo(i).side,...
        sidedot, ...
        sinfo(i).dataname,...
        tempdot, ...
        templateresampled,...
        namedot,...
        sinfo(i).name,...
        sinfo(i).ee));

      if isempty(strfind(sinfo(i).ff,'..')), PO{i} = strrep(PO{i},'..','.'); end
    end
  else
    PO = P;
  end
end
    