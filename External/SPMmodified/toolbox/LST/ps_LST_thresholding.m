function varargout = ps_LST_thresholding(varargin)
%ps_LST_thresholding   Thresholds probability lesion maps.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
% 
%   ps_LST_thresholding thresholds probability lesion maps in order to
%   obtain a binary lesion map.
%
%   ps_LST_thresholding asks the user for the lesion probability maps and
%   writes corresponding binary maps thresholded by 0.5.
%
%   ps_LST_thresholding(V, thr) thresholds the images given in V by the
%   threshold thr. V must be a character like a call from spm_select.
%

% Check input
if ~isempty(varargin) && isfield(varargin{1}, 'data_plm')
    viajob = 1;
else
    viajob = 0;
end
void = 0;
thr = 0.5;
if ~viajob
    if nargin == 0
        % Select files by function
        Vles = spm_select(Inf, 'image', 'Select probability lesion maps.');        
    end    
    if nargin > 0
        if ischar(varargin{1})
            Vles = varargin{1};
        else            
            fprintf('Input for Vles must be a character, like a call from spm_select.\n');            
            return;
        end        
    end
    if nargin == 2
        void = varargin{2};
        thr = varargin{2};
    end    
else    
    job = varargin{1};    
    Vles = job.data_plm;
end

if ~void
    % Welcome text
    fprintf('\n')
    fprintf(repmat('-', 1, 72));
    fprintf('\n')
    %fprintf('\nThis is LST')
    strout = 'This is LST';
    fprintf(strout)
    tt = 'www.statistical-modeling.de/lst.html\n';
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
    fprintf(strout)
    
    strout = 'Algorithm:';
    fprintf(strout)
    tt = 'Create binary lesion maps\n';
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
    fprintf(strout)
        
end

% Load volumes
Vles = spm_vol(Vles); 

% Check input
if numel(Vles) == 0
    fprintf('No images selected.\n');
    return;
end

% Summarize input
if ~void
    strout = 'Number of images to process:';
    fprintf(strout)
    tt = [num2str(numel(Vles)), '\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
    fprintf(strout)
end


if ~void
    strout = 'Progress';
    fprintf(strout)
    tt = '0.0%%';
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 1), tt];
    fprintf(strout)
end

for i = 1:numel(Vles)
    
    if viajob        
        Vles_tmp = Vles{i};    
    else            
        Vles_tmp = Vles(i);
    end       
    
    les = 1 .* (spm_read_vols(Vles_tmp) > thr);
    Vles_tmp.fname = strrep(Vles_tmp.fname, 'ples_', ['bles_', num2str(thr), '_']);
    Vles_tmp.descrip = strrep(Vles_tmp.descrip, 'Probability', 'Binary');
    spm_write_vol(Vles_tmp, les);
    
    if ~void
        fprintf(repmat('\b', 1, numel(tt)-1))        
        tt = [num2str(round(i / numel(Vles) * 1000) / 10), '%%'];        
        strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt];
        fprintf(tt)
    end
end

if ~void
    fprintf('\n')
    c = clock();
    strout = 'Finished successfully on';
    fprintf(strout)
    tt = [datestr(c), '\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
    fprintf(strout)
    fprintf(repmat('-', 1, 72));
    fprintf('\n')
end
varargout{:} = [];

return

