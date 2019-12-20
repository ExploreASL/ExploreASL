function varargout = ps_LST_update(varargin)
%ps_LST_tlv   Check if the newest version of LST is installed.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html

installed = importdata(fullfile(spm('dir'), 'toolbox', 'LST', 'lst-version.txt'));
proceed = 1;
try
    online = urlread('http://www.statistical-modelling.de/lst-version.txt');    
catch
    proceed = 0;
end

if proceed
    if strcmp(installed, online)
        if varargin{1} < 2
            fprintf(['You have installed the latest version of LST (', installed{1}, ').\n']);
        end
    else
        fprintf(['\nThere is a newer version of LST available (', online, ').\n']);
        fprintf('Visit www.statistical-modelling.de/lst.html for more information.\n\n');
    end
else
    if varargin{1} < 2
        fprintf('Internet connection required.\n')
    end
end
return;

