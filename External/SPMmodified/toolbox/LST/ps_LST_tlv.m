function varargout = ps_LST_tlv(varargin)
%ps_LST_tlv   Compute total lesion volume.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
%
%   ps_LST_tlv computes the total lesion volume (TLV) and the number of
%   lesions for a given set of lesion probability maps. The result
%   of this function is a CSV-file with the name LST_tlv_[date].csv
%   with the following four columns: 
%      Path     ... path of image
%      FileName ... name of image
%      LGA      ... dummy that indicate weather the lesion map was obtained
%                   by LGA (1) or LPA (0)
%      TLV      ... total lesion volume in ml
%      N        ... Number of lesions. This is probably not the best way to
%                   describe lesion pattern.
%
%   tlv = ps_LST_tlv(V) computes the TLV for the files given in V. V must 
%   be a character like a call from spm_select. Binary threshold is set to
%   its default value (0.5). The CSV-file is saved in MATLAB's curren
%   directory.
%
%   tlv = ps_LST_tlv(V, void, thr) calculates the TLV for the files given 
%   in V. V must be a character like a call from spm_select. If void is set
%   to 1 the algorithm does not produce any command line output during the
%   calculation. Binary threshold can be set by thr. The CSV-file is saved 
%   in MATLAB's curren directory.
%
%   tlv = ps_LST_tlv(job) Same as above but with job being a harvested job 
%   data structure (see matlabbatch help).
%
%   tlv = ps_LST_tlv asks for file names by spm_select beforehand and uses
%   void = 0 and thr = 0.5.
%
% ExploreASL Hack to allow setting MinimalLesionSize:
% lines 46, 67-71, 170
% ExploreASL Hacks with comments as explanation of the volumetric
% procedure: lines 164-174
%
% EXPLANATION FOR EXPLOREASL:
% FROM WMH SEGMENTATION (AFTER CLEANUP)
% 1) THIS SCRIPT EXCLUDES ALL VOXELS <thr
% THIS VOXELWISE probability THRESHOLD IS 0.5 BY DEFAULT
% BUT CAN BE DECREASED BECAUSE OF OUR CLEANUP (e.g. to 0.05)
% 
% VOLUME FACTOR = 0.001 for 1x1x1 mm ^3 voxels
% SO 
% 2) EXCLUSION OF LESION CLUSTERS < MINIMALLESIONSIZE
% WOULD REMOVE THOSE SMALLER THAN 15 voxels with pWMH==1 for default
% MinimalLesionSize = 0.015

% Check input
if ~isempty(varargin) && isfield(varargin{1}, 'data_lm')
    viajob = 1;
else
    viajob = 0;
end
% set defaults
void = 0;
thr = 0.5;
MinimalLesionSize = 0.015;
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
    end    
    if nargin == 3
        void = varargin{2};
        thr = varargin{3};
    end    
    if nargin == 4 %% ExploreASL Hack
        void = varargin{2};
        thr = varargin{3};
        MinimalLesionSize = varargin{4};
    end
else    
    job = varargin{1};    
    Vles = job.data_lm;
    thr = job.bin_thresh;
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
    tt = 'Compute total lesion volume\n';
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

% Check weather images come from lga or lpa
lga = zeros(numel(Vles), 1);
for i = 1:numel(Vles)
    if viajob
        name_tmp = Vles{i}.fname;
    else
        name_tmp = Vles(i).fname;
    end
    lga_tmp = regexp(name_tmp, '_lga_', 'once');
    if isempty(lga_tmp); lga_tmp = 0; else lga_tmp = 1; end
    lga(i) = lga_tmp;
end
if ~(sum(lga) == 0 || sum(lga) == numel(Vles))
    fprintf('\nWarning:\n   Total lesion volume may not be comparable between different lesion\n');
    fprintf('   segmentation algorithms! I hope you know what you are doing!\n\n');
end

% Create TLV file
nameTLV = ['LST_tlv_', num2str(thr), '_', ps_create_timestamp, '.csv'];
fileID = fopen(nameTLV, 'wt');
fprintf(fileID, 'Path,FileName,LGA,TLV,N\n');

if ~void
    strout = 'Results are written to';
    fprintf(strout)
    tt = [nameTLV, '\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
    fprintf(strout)
    strout = 'Progress';
    fprintf(strout)
    tt = '0.0%%';
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 1), tt];
    fprintf(strout)
end

for i = 1:numel(Vles) % iterate over images
    
    if viajob        
        Vles_tmp = Vles{i};    
    else            
        Vles_tmp = Vles(i);
    end   
    [pth, name, ext] = fileparts(Vles_tmp.fname);
    pth = strrep(pth, '\', '/');
    
    les = spm_read_vols(Vles_tmp);
    volfactor = abs(det(Vles_tmp.mat(1:3,1:3))) /  1000; % volume per voxel
    % thr = voxel-wise threshold
    % MinimalLesionSize = minimal volume of each lesions
    if any(les(:) > thr) % if any voxel survives the voxel-wise threshold
        bw = ps_bwlabeln(1 .* (les > thr)); % create binary mask above threshold
        % and label each region with connected WMH separately -> lesions
        les_size = zeros(max(bw(:)), 1); % create empty LesionSize column
        for j = 1:max(bw(:)) % for all lesions
            les_size(j) = sum(bw(:) == j) * volfactor; % calculate the volume of each lesion
        end
        les_size = les_size(les_size > MinimalLesionSize); % remove lesions below volume of MinimalLesionSize
    else
        les_size = [];
    end    
    fprintf(fileID, [pth, ',', name, ext, ',', num2str(lga(i)), ',', ...
         num2str(round(sum(les_size) * 1000) / 1000), ',', num2str(numel(les_size)), '\n']);
    
        
    if ~void
        fprintf(repmat('\b', 1, numel(tt)-1))        
        tt = [num2str(round(i / numel(Vles) * 1000) / 10), '%%'];        
        strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt];
        fprintf(tt)
    end
end

fclose(fileID);

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
varargout{:} = les_size;

return

