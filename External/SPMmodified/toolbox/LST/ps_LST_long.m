function varargout = ps_LST_long(varargin)
%ps_LST_long   Longitudinal lesion segmentation.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
%
%   ps_LST_long quantifies the change in lesion structure between different
%   lesion probability maps. The pipeline proceeds by comparing all
%   consecutive time points in an iterative manner. It decides if changes
%   in lesion structure are significant or due to natural variations of the
%   FLAIR signal. Non-significant changes are labeled as lesions in both 
%   probability maps, thus, probability lesion maps are corrected within 
%   this procedure and may differ from the ones that served as input. As a 
%   final result, lesion change labels are produced for all consecutive 
%   time points. In these images the three possible cases decrease, no 
%   change and increase are labeled by the numbers 1, 2, and 3, 
%   repsectively.
%   
%   The longitudinal pipeline works through the following steps:
%       1) Filling of 
%           (a) T1 images for all time points if LGA was used.
%           (b) FLAIR images for all time points if LGA was used without 
%               reference images.
%           (c) reference images for all time points if LGA was used with
%               reference images.
%       2) Coregistration of the filled image for time point t + 1 to time
%          point t, for t = 1, ..., m - 1. The estimated coregistrations 
%          are then applied to (FLAIR(t+1)) and (ples(t+1)).
%       3) Analyzation of relative change in FLAIR intensity 
%          (FLAIR(t) - FLAIR(t+1)) / ((FLAIR(t) + FLAIR(t+1)) / 2) in 
%          healthy WM in order to obtain the amount of normal/random deviation.
%       4) Information of 3) is used to threshold relative change in FLAIR 
%          intensity along observed lesions.
%
%   As a result this pipeline produces the following output:
%   lples_[].nii ... probability lesion map for all time points after 
%                    including all information in the segmentation process.
%                    Here, [] is either lga_[kappa]_FLAIR(t) for LGA or l
%                    pa_namFLAIR(t+1) for LPA. All lesion maps are in the 
%                    same space as the image for t = 1. The 'l' indicates 
%                    'longitudinal'.
%   LCL_[].nii   ... lesion change label, i.e. an image in the same space 
%                    as the first image that indicates if lesions 
%                    disapeared (1), remain constant (2), or appeared new 
%                    (3). Here, [] is something like 
%                    [LGA/LPA]_FLAIR(t)_FLAIR(t+1)
%   report_LST_long_[] ... A HTML report of the segmentation if sepcified
%                          by the user.
%   LST_long_[]  ... a folder that contains a report (HTML) along with some
%                    images. 
%   
%   ps_LST_long(m) asks the user for the lesion probability (lesion 
%   probability maps computed by LGA or LPA, but not both) maps of the m 
%   time points.

%   ps_LST_long(m, html) same as above but with the possibility to
%   deactivate the creation of the HTML report.
%
%   ps_LST_long(job) Same as above but with job being a harvested job data
%   structure (see matlabbatch help).
%

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
tt = 'Longitudinal segmentation\n';
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)

% Create log file
pthor = cd;
nameLog = ['LST_log_', ps_create_timestamp, '.txt'];
fileID = fopen(nameLog, 'wt');
strout = 'If anything goes wrong:';
fprintf(strout)
tt = [nameLog, '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)

% Check input
if ~isempty(varargin) && isfield(varargin{1}, 'data_long_tmp')
    viajob = 1;
else
    viajob = 0;
end
if ~viajob
    if nargin == 0
        fprintf('I need at least the number of time points!\n');
        fclose(fileID);
        spm_unlink(nameLog);
        return;
    end
    if nargin > 2        
        fprintf('Wrong number of input parameters!\n');
        fclose(fileID);
        spm_unlink(nameLog);
        return;
    else
        % If function is called without a harvested job data structure input
        % must be a scalar that gives the number of time points.        
        if isnumeric(varargin{1})
            fprintf(fileID, 'Select files by function ... ');
            m = varargin{1};
            Vles = cell(m, 1);
            for j = 1:m
                Vles{j} = spm_select(Inf, 'image', ['Select lesion probability maps for t = ', num2str(j)]);
            end
            fprintf(fileID, 'ok.\n');
        else
            fprintf(['Input for this function must either be a harvested job ', ...
                'data structure, a numeric value, or a cell with cells of ', ...
                'image names.\n']);                
            fclose(fileID);
            spm_unlink(nameLog);
            return;
        end        
        if nargin > 1
            html_report = varargin{2};
        else
            html_report = 1;
        end
    end    
else
    job = varargin{1};
    Vles = job.data_long_tmp;
    html_report = job.html_report;
    m = numel(Vles);    
end

fprintf(fileID, 'Load volume header ... ');
Vles = spm_vol(Vles);
fprintf(fileID, 'ok.\n');
%Vles{1}
%numel(Vles{1})
%arrayfun(@(x) numel(Vles{x}), (1:m))
% Check if all time points have the same number of images
if ~all(cellfun(@numel, Vles) == numel(Vles{1}))
    fprintf('Number of images for all time points must match.\n');
    fclose(fileID);
    spm_unlink(nameLog);
    return;
end
if numel(Vles{1}) == 0
    fprintf('No images selected.\n');
    fclose(fileID);
    spm_unlink(nameLog);
    return;
end

fprintf(fileID, 'Longitudinal segmentation\n');

% Summarize input
fprintf(fileID, 'Input summary:\n');
fprintf(fileID, ['Jobs: ', num2str(numel(Vles{1})), '\n']);
strout = 'Number of jobs to process:';
fprintf(strout)
tt = [num2str(numel(Vles{1})), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)
fprintf(fileID, ['Time points: ', num2str(m), '\n']);
strout = 'Number of time points:';
fprintf(strout)
tt = [num2str(m), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)

% Loop over all subjects
for i = 1:numel(Vles{1})
    
    fprintf(fileID, '\n******************************\n');
    fprintf(fileID, ['Job ', num2str(i), ' of ', num2str(numel(Vles{1})), '\n']);
    fprintf(fileID, '*******************************\n');
    
    % Which subject?    
    strout = '\nWorking on job';
    fprintf(strout)
    tt = [num2str(i), ' out of ', num2str(numel(Vles{1})), ' (', num2str(i/numel(Vles{1})*100), '%%)\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 5), tt];
    fprintf(strout)        
       
    % Extract names of the ith subject
    Vles_tmp = cell(m, 1);
    pthles = Vles_tmp; namles = pthles; extles = pthles;
    for j = 1:m        
        if viajob
            Vles_tmp{j} = Vles{j}{i};
        else
            Vles_tmp{j} = Vles{j}(i);
        end        
        [pthles{j}, namles{j}, extles{j}] = fileparts(Vles_tmp{j}.fname);
    end
    
    % ---------------------------------------------------------------------
    % Preprocessing
    % ---------------------------------------------------------------------
        
    % Check weather lga or lpa was used for lesion segmentation.
    ch_lga = cellfun(@(x) ~isempty(regexp(x, '_lga_') > 0), namles);
    if ~all(ch_lga == ch_lga(1))
        fprintf('Images are not from the same lesion segmentation algorithm.\n');
        fclose(fileID);
        spm_unlink(nameLog);
        return;
    end
    
    % Get some information about the lesion segmentation    
    if ch_lga(1)
                        
        % get FLAIR names
        alg_text = 'LGA';
        us_pos = cellfun(@(x) regexp(x, '_'), namles, 'UniformOutput', false);
        namf2 = arrayfun(@(x) namles{x}((us_pos{x}(3)+1):end), 1:m, 'UniformOutput', false);
        
        % lga.mat files
        lst = cell(m, 1);
        for j = 1:m
            cd(pthles{j})
            load(['LST_lga_', namf2{j}, '.mat'])
            lst{j} = lga;
        end               
        cd(pthles{1})
        
        % name of images that need to be filled/coregistered
        namcoreg = cellfun(@(x) x.Vt1.fname, lst, 'UniformOutput', false);
        %for j = 1:numel(lst)
        %    namcoreg{j} = [ps_fileparts(namcoreg{j}, 1:2), '_', num2str(j), ps_fileparts(namcoreg{j}, 3)];
        %end
        
        % Create folder for report and stuff        
        nameFolder = ['LST_long_lga_', num2str(m), '_', ps_create_timestamp];
        mkdir(nameFolder)
        
        % write p0 for t = 2 , ..., t = m
        for j = 1:m
            p0 = zeros(lst{j}.Vt1.dim); p0(lst{j}.indx_brain) = lst{j}.p0_vec;
            Vp0 = lst{j}.Vt1; Vp0.fname = fullfile(nameFolder, ['p0_', num2str(j), '.nii']);
            spm_write_vol(Vp0, p0);
        end
    else
        
        % get FLAIR names
        alg_text = 'LPA';
        us_pos = cellfun(@(x) regexp(x, '_'), namles, 'UniformOutput', false);
        namf2 = arrayfun(@(x) namles{x}((us_pos{x}(2)+1):end), 1:m, 'UniformOutput', false);
        
        % lpa.mt files
        lst = cell(m, 1);
        for j = 1:m
            cd(pthles{j})
            load(['LST_lpa_', namf2{j}, '.mat'])
            lst{j} = lpa;
        end
        cd(pthles{1})
        
        % name of images that need to be filled/coregistered
        if isfield(lst{1}, 'Vref')
            namcoreg = cellfun(@(x) x.Vref.fname, lst, 'UniformOutput', false);            
        else
            namcoreg = cellfun(@(x) x.Vf2.fname, lst, 'UniformOutput', false);            
        end
        %for j = 1:numel(lst)
        %    namcoreg{j} = [ps_fileparts(namcoreg{j}, 1:2), '_', num2str(j), ps_fileparts(namcoreg{j}, 3)];
        %end
        
        % Create folder for report and stuff
        nameFolder = ['LST_long_lpa_', num2str(m), '_', ps_create_timestamp];
        mkdir(nameFolder)        
        
        % write p0 for t = 2, ..., t = m
        for j = 1:m
            p0 = zeros(lst{j}.Vf2.dim); p0(lst{j}.indx_brain) = lst{j}.I;
            Vp0 = lst{j}.Vf2; Vp0.fname = fullfile(nameFolder, ['p0_', num2str(j), '.nii']);
            spm_write_vol(Vp0, p0);
        end
            
    end
    or = lst{1}.or;
    fl = lst{1}.fl;
    
    % ---------------------------------------------------------------------
    % Filling
    % ---------------------------------------------------------------------
    
    fprintf(fileID, 'Filling ...');
    strout = 'Fill images ';
    fprintf(strout)
    tic
    
    % Copy images to temporary folder
    %cellfun(@(x) copyfile(x, fullfile(nameFolder, '.')), namcoreg);    
    for j = 1:m        
        copyfile(namcoreg{j}, fullfile(nameFolder, [ps_fileparts(namcoreg{j}, 2), '_', num2str(j), ps_fileparts(namcoreg{j}, 3)]));
        namcoreg{j} = [ps_fileparts(namcoreg{j}, 1:2), '_', num2str(j), ps_fileparts(namcoreg{j}, 3)];
    end
    for j = 1:m
        if strcmp(ps_fileparts(namcoreg{j}, 3), '.img')
            copyfile([ps_fileparts(namcoreg{j}, 1:2), '.hdr'], fullfile(nameFolder, '.'))
        end
    end
    namcoreg = cellfun(@(x) fullfile(pthles{1}, nameFolder, ps_fileparts(x, 2:3)), namcoreg, 'UniformOutput', false);
        
    % fill images
    nam_del = cell(m, 1);
    for j = 1:m
        nam_del{j} = namcoreg{j};
        namcoreg{j} = fullfile(pthles{1}, nameFolder, ps_LST_lesfill(namcoreg{j}, Vles_tmp{j}.fname, 1, 0));
        spm_unlink(nam_del{j});        
    end    
                
    tt = toc; tt = [num2str(round(tt)), 's'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
    fprintf(strout)
    fprintf(fileID, ' ok.\n');    
    
    % ---------------------------------------------------------------------
    % Coregistration
    % ---------------------------------------------------------------------
    
    [~, job2] = ps_LST_lpa_preproc_default;
    fprintf(fileID, 'Coregistration ...');
    strout = 'Coregistration ';
    fprintf(strout)
    tic
        
    job2.roptions.prefix = 'rl';    
    for j= 2:m                
        cd(pthles{1})
        copyfile(fullfile(pthles{j}, [namf2{j}, '.nii']), fullfile(nameFolder, [namf2{j}, '_', num2str(j), '.nii']))
        copyfile(Vles_tmp{j}.fname, fullfile(nameFolder, [ps_fileparts(Vles_tmp{j}.fname, 2), '_', num2str(j), '.nii']))
        job2.ref = {namcoreg{1}};
        job2.source = {namcoreg{j}};
        %job2.other = {fullfile(pthles{j}, [namf2{j}, '.nii']), ...
        %               Vles_tmp{j}.fname, ...
        %               fullfile(nameFolder, ['p0_', num2str(j), '.nii'])};
        job2.other = {fullfile(nameFolder, [namf2{j}, '_', num2str(j), '.nii']), ...
                       fullfile(nameFolder, [namles{j}, '_', num2str(j), '.nii']), ...
                       fullfile(nameFolder, ['p0_', num2str(j), '.nii'])};
        ps_LST_spm_run_coreg(job2);
        spm_unlink(namcoreg{j})
        spm_unlink(fullfile(nameFolder, ['p0_', num2str(j), '.nii']))
        copyfile(fullfile(nameFolder, ['rl', namf2{j}, '_', num2str(j), '.nii']), pthles{j})
        copyfile(fullfile(nameFolder, ['rl', namles{j}, '_', num2str(j), '.nii']), pthles{j})
    end
    
    copyfile(fullfile(pthles{1}, [namles{1}, extles{1}]), fullfile(pthles{1}, ['rl', namles{1}, '_1', extles{1}]));
    copyfile(fullfile(pthles{1}, [namf2{1}, '.nii']), fullfile(pthles{1}, ['rl', namf2{1}, '_1', '.nii']));
    copyfile(fullfile(pthles{1}, nameFolder, 'p0_1.nii'), fullfile(pthles{1}, nameFolder, 'rlp0_1.nii'));
    spm_unlink(fullfile(pthles{1}, nameFolder, 'p0_1.nii'));
    
    tt = toc; tt = [num2str(round(tt)), 's'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
    fprintf(strout)
    fprintf(fileID, ' ok.\n');
    
    % ---------------------------------------------------------------------
    % Longitudinal
    % ---------------------------------------------------------------------
    
    fprintf(fileID, 'Load data ...');
    strout = 'Compare time points ';
    fprintf(strout)
    tic
    
    changes = zeros([Vles_tmp{1}.dim m - 1]);    
    counter = 0;
    ch = zeros(1, m - 1);
    while any(ch == 0) && counter < 30
        counter = counter + 1;
        fprintf(fileID, '\n***********************\n');
        fprintf(fileID, ['Iteration ', num2str(counter), '\n']);
        fprintf(fileID, '***********************');
        for j = 1:(m-1)
            fprintf(fileID, ['\n*** Compare time point ', num2str(j), ' with time point ', num2str(j+1), '\n']);
            fprintf(fileID, 'Load images ...');
            if counter == 1
                les_01 = spm_read_vols(spm_vol(fullfile(pthles{j}, ['rl', namles{j}, '_', num2str(j), '.nii'])));
                les_02 = spm_read_vols(spm_vol(fullfile(pthles{j+1}, ['rl', namles{j+1}, '_', num2str(j+1), '.nii'])));
            else
                les_01 = spm_read_vols(spm_vol(fullfile(pthles{j}, ['l', namles{j}, '_', num2str(j), '.nii'])));
                les_02 = spm_read_vols(spm_vol(fullfile(pthles{j+1}, ['l', namles{j+1}, '_', num2str(j+1), '.nii'])));
            end
            les_01 = 1 .* (les_01 > 0.1);
            les_02 = 1 .* (les_02 > 0.1);
            les_0102 = les_02 - les_01;

            joint = 1 .* (les_01 > 0 | les_02 > 0);
            Vf2 = spm_vol(fullfile(pthles{j}, ['rl', namf2{j}, '_', num2str(j), '.nii']));
            f2_01 = spm_read_vols(Vf2);
            f2_02 = spm_read_vols(spm_vol(fullfile(pthles{j+1}, ['rl', namf2{j+1}, '_', num2str(j+1), '.nii'])));
            p0_01 = spm_read_vols(spm_vol(fullfile(nameFolder, ['rlp0_', num2str(j), '.nii'])));
            p0_02 = spm_read_vols(spm_vol(fullfile(nameFolder, ['rlp0_', num2str(j+1), '.nii'])));
            p0_01(isnan(f2_01) | isnan(f2_02)) = 0;
            p0_02(isnan(f2_01) | isnan(f2_02)) = 0;                
            fprintf(fileID, ' ok.\n');

            fprintf(fileID, 'Normalize FLAIR images ...');
            % Normalize FLAIR    
            tmp = f2_01((p0_01 > 1.5 & p0_01 < 2.5) | (p0_02 > 1.5 & p0_02 < 2.5) & joint < 1);
            f2_01 = f2_01 ./ mean(tmp(~isnan(tmp)));
            tmp = f2_02((p0_01 > 1.5 & p0_01 < 2.5) | (p0_02 > 1.5 & p0_02 < 2.5) & joint < 1);
            f2_02 = f2_02 ./ mean(tmp(~isnan(tmp)));
            fprintf(fileID, ' ok.\n');

            fprintf(fileID, 'Relative change ...');
            % relative change
            f2_relchange = (f2_02 - f2_01) ./ ((f2_01 + f2_02) ./ 2);
            f2_relchange(isinf(f2_relchange)) = NaN;
            f2_relchange(p0_01 == 0) = 0; f2_relchange(p0_02 == 0) = 0;
            fprintf(fileID, ' ok.\n');

            fprintf(fileID, 'Healthy WM ...');
            % Healthy WM
            tmp = f2_relchange(p0_01 > 2.9 & p0_01 < 3.1 & (p0_02 > 2.9 & p0_02 < 3.1) & joint == 0);    
            %tmp = f2_relchange(p0_01 == 3 & (p0_02 > 2.9 & p0_02 < 3.1) & joint == 0);    
            tmp(isnan(tmp)) = [];
            m3 = mean(tmp);
            sd3 = std(tmp);
            fprintf(fileID, [' ok, mean = ', num2str(m3), ', sd = ', num2str(sd3), '.\n']);

            % Use the smoothed relchange for thresholding
            fprintf(fileID, 'Smoothing and thresholding ...');
            sf2_relchange_les = joint .* f2_relchange;
            spm_smooth(sf2_relchange_les, sf2_relchange_les, [1, 1, 1] .* 2);            
            %thr1 = norminv(.1, m3, sd3);
            thr1 = ps_qnorm(.1, m3, sd3);
            %thr2 = norminv(.9, m3, sd3);
            thr2 = ps_qnorm(.9, m3, sd3);
            change = les_0102 .* (sf2_relchange_les < thr1 | sf2_relchange_les > thr2);
            fprintf(fileID, [' ok, thr1 = ', num2str(thr1), ', thr2 = ', num2str(thr2), '.\n']);

            % delete all change voxels that lie directly on CSF
            indx_tmp = find(change ~= 0);
            csf = 1 .* (p0_01 < 1.5 | p0_02 < 1.5);
            csf = ps_set_border_zero(csf);
            nh = getNeighborhood2(csf, indx_tmp, 3);
            change(indx_tmp(sum(nh > 0) > 0)) = 0;

            % Process all voxels that were identified as change
            fprintf(fileID, 'Postprocessing ...');
            st = 0;
            while ~st
                indx_les_02 = find(les_02 > 0 & les_01 == 0 & change == 0);
                nh = getNeighborhood2(change, indx_les_02, 1);
                indx_tmp = indx_les_02(sum(nh > 0) > 0 & abs(sf2_relchange_les(indx_les_02))' > thr2*.5 & joint(indx_les_02)' > 0);
                if isempty(indx_tmp)
                    st = 1;
                else
                    change(indx_tmp) = 1;
                end
            end
            fprintf(fileID, ' ok ... ');
            st = 0;
            while ~st
                indx_les_01 = find(les_01 > 0 & les_02 == 0 & change == 0);
                nh = getNeighborhood2(change, indx_les_01, 1);
                indx_tmp = indx_les_01(sum(nh < 0) > 0 & abs(sf2_relchange_les(indx_les_01))' > thr1*.5 & joint(indx_les_01)' > 0);
                if isempty(indx_tmp)
                    st = 1;
                else
                    change(indx_tmp) = -1;
                end
            end
            fprintf(fileID, ' ok.\n');
            
            % Delete all lesions that are smaller than 0.01 ml
            fprintf(fileID, 'Delete voxels that are smaller than 0.01 ml ...');
            volfactor = abs(det(Vf2.mat(1:3,1:3))) /  1000;            
            for k = [-1,1]
                %b = bwconncomp(1 .* (change == k), 6);
                %change(cell2mat(b.PixelIdxList(cellfun(@numel, b.PixelIdxList) .* volfactor < 0.01)')) = 0;
                if sum(change(:) == k) > 0
                    b = ps_bwlabeln(1 .* (change == k));
                    if max(b(:)) > 0
                        c_tmp = ps_count(b(b > 0));
                        for kk = 1:size(c_tmp, 2)
                            if (c_tmp(2,kk) * volfactor)  < 0.01
                                change(b == c_tmp(1,kk)) = 0;
                            end
                        end
                    end
                end
            end
            clear change2;
            fprintf(fileID, ' ok.\n');

            fprintf(fileID, 'Create new lesion maps ...');    

            les_01_new = 0.*les_01;
            les_02_new = 0.*les_02;
            les_01_new(change < 0 | (les_01 > 0 & les_02 > 0)) = 1;
            les_02_new(change > 0 | (les_01 > 0 & les_02 > 0)) = 1;

            change = 0 .* les_01;
            change(les_01_new > 0 & les_02_new == 0) = 1;
            change(les_01_new > 0 & les_02_new > 0) = 2;
            change(les_01_new == 0 & les_02_new > 0) = 3; 
            tmp = changes(:,:,:,j) - change;
            ch(j) = numel(tmp(tmp ~= 0)) < 3;
            %ch(j) = isequal(changes(:,:,:,j), change);
            changes(:,:,:,j) = change;
            
            % Adjust old lesion maps
            if counter == 1
                les_01_or = spm_read_vols(spm_vol(fullfile(pthles{j}, ['rl', namles{j}, '_', num2str(j), '.nii'])));
                Vles_tmp{j+1} = spm_vol(fullfile(pthles{j+1}, ['rl', namles{j+1}, '_', num2str(j+1), '.nii']));
                les_02_or = spm_read_vols(Vles_tmp{j+1});
            else
                les_01_or = spm_read_vols(spm_vol(fullfile(pthles{j}, ['l', namles{j}, '_', num2str(j), '.nii'])));
                les_02_or = spm_read_vols(spm_vol(fullfile(pthles{j+1}, ['l', namles{j+1}, '_', num2str(j+1), '.nii'])));
            end
            les_02_or(les_02_or > 1) = 1;
            les_or = [les_01_or(:), les_02_or(:)];
            les_01_or(change == 2) = max(les_or(change == 2,:), [], 2);
            les_02_or(change == 2) = max(les_or(change == 2,:), [], 2);
            les_01_or(change == 3) = 0;
            les_02_or(change == 1) = 0;
            les_01_or(change == 0) = 0;
            les_02_or(change == 0) = 0;
            Vles_tmp{j}.fname = fullfile(pthles{j}, ['l', namles{j}, '_', num2str(j), extles{j}]);
            Vles_tmp{j}.descrip = 'Probability lesion map obtained by longitudinal pipeline within LST toolbox';
            spm_write_vol(Vles_tmp{j}, les_01_or);
            Vles_tmp{j+1}.fname = fullfile(pthles{j+1}, ['l', namles{j+1}, '_', num2str(j+1), extles{j+1}]);
            Vles_tmp{j+1}.descrip = 'Probability lesion map obtained by longitudinal pipeline within LST toolbox';            
            spm_write_vol(Vles_tmp{j+1}, les_02_or);
        end
    end
    
    % save lesion change label
    for j = 1:(m-1)
        Vles_tmp{j}.fname = fullfile(pthles{j}, ['LCL_', namles{j}, '_', namles{j+1}, '.nii']);
        Vles_tmp{j}.descrip = ['Lesion change label for timepoint ', num2str(j), ' and ', num2str(j+1)];
        spm_write_vol(Vles_tmp{j}, changes(:,:,:,j));
    end
    
    % Change name of images
    for j = 1:m
        movefile(fullfile(pthles{j}, ['l', namles{j}, '_', num2str(j), extles{j}]), ...
            fullfile(pthles{j}, ['l', namles{j}, extles{j}]))
        movefile(fullfile(pthles{j}, ['rl', namf2{j}, '_', num2str(j), '.nii']), ...
            fullfile(pthles{j}, ['rl', namf2{j}, '.nii']))
    end
    
        
    tt = toc; tt = ['finished after ', num2str(counter) ' iterations, ', num2str(round(tt)), 's'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
    fprintf(strout)
    fprintf(fileID, ' ok.\n\n');               
    
    if html_report
        
        % -----------------------------------------------------------------
        % Create HTML report
        % -----------------------------------------------------------------
        
        stroutHTML = 'Create HTML report';
        fprintf(stroutHTML);
        
        fprintf(fileID, 'Create main HTML file ...');
        % Main HTML file
        copyfile(fullfile(spm('dir'), 'toolbox', 'LST', 'LST_main_html.html'), ['report_', nameFolder, '.html'])
        HTMLid = fopen(['report_', nameFolder, '.html'], 'at');
        strout = ['    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'raphael.js'), '\"></script>\n', ...
                  '    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery.min.js'), '\"></script>\n', ...
                  '    <link href=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.css'), '\" rel=\"stylesheet\"></script>\n', ...
                  '    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.js'), '\"></script>\n', ...
                  '  </head>\n  <body>\n'];
        fprintf(HTMLid, strout);
        fprintf(fileID, ' ok.\n');      
        strout = '';
        %strout = '';
            
        fprintf(fileID, 'Create PNGs ...');
        for j = 1:(m-1)
            if j == 1
                Vimg1 = fullfile(pthles{j}, [namf2{j}, '.nii']);
            else
                Vimg1 = fullfile(pthles{j}, ['rl', namf2{j}, '.nii']);
            end
            Vimg2 = fullfile(pthles{j+1}, ['rl', namf2{j+1}, '.nii']);
            if numel(Vimg1) < numel(Vimg2)
                Vimg1 = [Vimg1, repmat(' ', 1, numel(Vimg2) - numel(Vimg1))];
            else
                Vimg2 = [Vimg2, repmat(' ', 1, numel(Vimg1) - numel(Vimg2))];
            end
            Vimg = [Vimg1; Vimg2];

            p0_01 = spm_read_vols(spm_vol(fullfile(nameFolder, ['rlp0_', num2str(j), '.nii'])));
            p0_02 = spm_read_vols(spm_vol(fullfile(nameFolder, ['rlp0_', num2str(j+1), '.nii'])));
            c_tmp = indx2coord(find(p0_01 > .9 & p0_02 > .9), size(p0_01, 1), size(p0_01, 2));
            r = [min(c_tmp(:,3)),max(c_tmp(:,3))];
            pngFailed = '';
            try
                [~, r] = ps_LST_create_gif(Vimg, fullfile(pthles{j}, ['LCL_', namles{j}, '_', namles{j+1}, '.nii']), fullfile(pthles{1}, nameFolder, ['overlay_', num2str(j), '_', num2str(j+1)]), r, [or, fl]);
            catch ME
                fprintf(fileID, ' failed!.\n');             
                pngFailed = ME.message;
            end            
            
            fprintf(fileID, 'Create glass brains ...');  
            change = changes(:,:,:,j);
            % Create images for glass brains
            ps_LST_create_glass_brain(1 .* (change == 1), ...
                1 .* (p0_01 > 0.5 & p0_02 > 0.5), ...
                fullfile(nameFolder, ['gb_decreased_', num2str(j), '_', num2str(j+1)]), [or, fl]);
            ps_LST_create_glass_brain(1 .* (change == 2), ...
                1 .* (p0_01 > 0.5 & p0_02 > 0.5), ...
                fullfile(nameFolder, ['gb_unchanged_', num2str(j), '_', num2str(j+1)]), [or, fl]);
            ps_LST_create_glass_brain(1 .* (change == 3), ...
                1 .* (p0_01 > 0.5 & p0_02 > 0.5), ...
                fullfile(nameFolder, ['gb_increased_', num2str(j), '_', num2str(j+1)]), [or, fl]);
                        
            fprintf(fileID, ' ok.\n');
            
            
            fprintf(fileID, 'Create subject specific HTML file ...');
            les_01_or = spm_read_vols(spm_vol(fullfile(pthles{j}, ['l', namles{j}, '.nii'])));
            les_02_or = spm_read_vols(spm_vol(fullfile(pthles{j+1}, ['l', namles{j+1}, '.nii'])));
            
            if any(les_01_or(:) > .5)
                b = ps_bwlabeln(1*(les_01_or > 0.5));
                if any(b(:) > 0)
                    for k = 1:max(b(:))
                        if (sum(b(:) == k) * volfactor) <= 0.015
                            b(find(b) == k) = 0;
                        end
                    end
                end
            else
                b = 0 .* les_01_or;
            end
            les_01_or = les_01_or .* (b > 0);
            
            if any(les_02_or(:) > .5)
                b = ps_bwlabeln(1*(les_02_or > 0.5));
                if any(b(:) > 0)
                    for k = 1:max(b(:))
                        if (sum(b(:) == k) * volfactor) <= 0.015
                            b(find(b) == k) = 0;
                        end
                    end
                end
            else
                b = 0 .* les_02_or;
            end
            les_02_or = les_02_or .* (b > 0);
            
            tlv1 = sum(les_01_or(:) > .5) * volfactor;
            tlv2 = sum(les_02_or(:) > .5) * volfactor;    
            joint = 0    .* les_01_or;
            joint(:) = max([les_01_or(:), les_02_or(:)], [], 2);
            tlv_joint = sum(joint(:) > .5) * volfactor;
            tlv_unch = sum(joint(change(:) == 2) > .5) * volfactor;
            tlv_decr = sum(joint(change(:) == 1) > .5) * volfactor;
            tlv_incr = sum(joint(change(:) == 3) > .5) * volfactor; 
                
            if any(les_01_or(:) > .5)
                b = ps_bwlabeln(1 .* (les_01_or > .5));
            else
                b = 0 .* les_01_or;
            end
            numles_01 = max(b(:));
            if any(les_02_or(:) > .5)
                b = ps_bwlabeln(1 .* (les_02_or > .5));
            else
                b = 0 .* les_02_or;
            end
            numles_02 = max(b(:));
            
            jsid = [nameFolder, '_', num2str(j), num2str(j+1)];
            jsid(regexp(jsid, '\.')) = [];
                
            strout = [strout, ...
                      '\n<div class=\"container\">\n', ...              
                      '  <h1>Longitudinal lesion segmentation by LST</h1>\n', ...
                      '  <div class=\"column-01\">\n', ...
                      '    <h2>Input summary</h2>\n', ...
                      '      <table style=\"min-width: 500px;\">\n', ...
                      '        <tr>\n', ...
                      '            <td>Date of analysis</td>\n', ...
                      '            <td class=\"ta_right\">', datestr(clock()), '</td>\n', ...
                      '        </tr>\n', ...
                      '        <tr>\n', ...
                      '          <td>Directory for t = ' num2str(j), '</td>\n', ...
                      '          <td class=\"ta_right\">', ps_shorten_string(pthles{1}, 28), '</td>\n', ...
                      '       </tr>\n', ...
                      '       <tr>\n', ...
                      '           <td>Directory for t = ' num2str(j+1), '</td>\n', ...
                      '           <td class=\"ta_right\">', ps_shorten_string(pthles{2}, 28), '</td>\n', ...
                      '       </tr>\n', ...
                      '       <tr>\n', ...
                      '           <td>Lesion map for t = ' num2str(j), '</td>\n', ...
                      '           <td class=\"ta_right\">', namles{j}, extles{j}, '</td>\n', ...
                      '       </tr>\n', ...
                      '       <tr>\n', ...
                      '           <td>Lesion map for t = ' num2str(j+1), '</td>\n', ...
                      '           <td class=\"ta_right\">', namles{j+1}, extles{j+1}, '</td>\n', ...
                      '       </tr>\n', ...
                      '       <tr>\n', ...
                      '           <td>Algorithm used for segmentation</td>\n', ...
                      '           <td class=\"ta_right\">', alg_text, '</td>\n', ...
                      '       </tr>\n', ...
                      '   </table>\n', ...
                      '   </div>\n', ...              
                      '   <div class=\"column-02\">\n', ...
                      '     <h2>Results</h2>\n', ...
                      '         <table style=\"width: 500px\">\n', ...
                      '             <tr>\n', ...
                      '                 <td>Lesion volume for t = ' num2str(j), '</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(tlv1), ' ml</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Number of lesions for t = ' num2str(j), '</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(numles_01), '</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Lesion volume for t = ' num2str(j+1), '</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(tlv2), ' ml</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Number of lesions for t = ' num2str(j+1), '</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(numles_02), '</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Joint lesion volume</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(tlv_joint), ' ml</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Unchanged lesion volume</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(tlv_unch), ' ml</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Decreased lesion volume</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(tlv_decr), ' ml (', num2str(tlv_decr / tlv1 * 100, 3), '%%)</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Increased lesion volume</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(tlv_incr), ' ml (', num2str(tlv_incr / tlv1 * 100, 3), '%%)</td>\n', ...
                      '             </tr>\n', ...
                      '             <tr>\n', ...
                      '                 <td>Lesion volume change</td>\n', ...
                      '                 <td class=\"ta_right\">', num2str(-tlv_decr + tlv_incr), ' ml (', (num2str((-tlv_decr + tlv_incr) / tlv1 * 100, 3)), '%%)</td>\n', ...
                      '             </tr>\n', ...
                      '        </table>\n', ...
                      '    </div>\n', ...
                      '    <div style=\"clear:both\"></div>\n', ...
                      '    <div class=\"column-01\">\n', ...
                      '        <h2>Lesion change plot (LCP)</h2>\n', ... %'        <script src=\"', fullfile(cd, nameFolder, ['lcp_', num2str(id), '.js']), '\" type=\"text/javascript\"></script>\n', ...
                      '        <script src=\"', fullfile(cd, nameFolder, ['lcp_', num2str(j), num2str(j+1), '.js']), '\" type=\"text/javascript\"></script>\n', ...
                      '        <div id=\"canvas_', nameFolder, '_', num2str(j), num2str(j+1), '\" class=\"canvas\"></div>\n', ...
                      '    </div>\n', ...
                      '    <div class=\"column-02\" style=\"vertical-align: top;\">\n', ...
                      '        <h2>Overlay</h2>\n'];
                  if strcmp(pngFailed, '')
                      strout = [strout, ...
                          '        <img width=\"450px\" id=\"overlay', jsid, '\" src=\"', fullfile(cd, nameFolder, ['overlay_', num2str(j), '_', num2str(j+1), '_', num2str(round(mean(r))), '.png']), '\" />\n', ...
                          '        <div id=\"slider_', jsid, '\" style=\"width: 450px; text-align: center;\"></div>\n', ...
                          '        <div style=\"width: 450px; text-align: center;\">\n', ...
                          '           <button id=\"button-left', jsid, '\">\n', ...
                          '              <\n', ...
                          '           </button>\n', ...
                          '           <span id="slice', jsid, '">Slice ', num2str(round(mean(r))), '</span>\n', ...
                          '           <button id=\"button-right', jsid, '\">\n', ...
                          '              >\n', ...
                          '           </button>\n', ...
                          '        </div>\n'];
                  else
                      strout = [strout, ...
                        'Sorry, there was a problem when creating the PNG images. MATLAB said: ',  pngFailed, '\n'];
                  end                          
                  strout = [strout, ...    
                      '    </div>\n', ...
                      '    <div style=\"clear:both\"></div>\n', ...
                      '    <div>\n', ...
                      '  <h2>Change location</h2>\n', ...
                      '  <div class=\"column-center\">\n', ...
                      '    <h4 style=\"text-align: center\">Unchanged ', num2str(tlv_unch), ' ml</h4>\n', ...%'    <img src=\"', fullfile(cd, nameFolder, ['c22_', id, '.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_unchanged_', num2str(j), '_', num2str(j+1), '_2.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_unchanged_', num2str(j), '_', num2str(j+1), '_1.png']), '\" width=\"169.41px\" style=\"vertical-align: top;\"><br>\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_unchanged_', num2str(j), '_', num2str(j+1), '_3.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
                      '  </div>\n', ...
                      '  <div class=\"column-left\">\n', ...
                      '    <h4 style=\"text-align: center\">Decrease ', num2str(tlv_decr), ' ml (', num2str(tlv_decr / tlv1 * 100, 3), '%%)</h4>\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_decreased_', num2str(j), '_', num2str(j+1), '_2.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_decreased_', num2str(j), '_', num2str(j+1), '_1.png']), '\" width=\"169.41px\" style=\"vertical-align: top;\"><br>\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_decreased_', num2str(j), '_', num2str(j+1), '_3.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
                      '  </div>\n', ...
                      '  <div class=\"column-right\">\n', ...
                      '    <h4 style=\"text-align: center\">Increase ', num2str(tlv_incr), ' ml (', num2str(tlv_incr / tlv1 * 100, 3), '%%)</h4>\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_increased_', num2str(j), '_', num2str(j+1), '_2.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_increased_', num2str(j), '_', num2str(j+1), '_1.png']), '\" width=\"169.41px\" style=\"vertical-align: top;\"><br>\n', ...
                      '    <img src=\"', fullfile(cd, nameFolder, ['gb_increased_', num2str(j), '_', num2str(j+1), '_3.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...               
                      '    </div>\n', ...
                      '  </div>\n', ...
                      '</div>\n']; %% !!
                
    
    
            
            % Lesion change plot
            fprintf(fileID, 'Add to lcp.js ...');
            %b = bwconncomp(1*(change > 0), 6);
            if any(change(:) > 0)
                b = ps_bwlabeln(1*(change > 0));
                if any(b(:) > 0)
                    for k = 1:max(b(:))
                        if (sum(b(:) == k) * volfactor) <= 0.015
                            b(find(b) == k) = 0;
                        end
                    end
                end
                b = ps_bwlabeln(1*(b > 0));
            else
                b = 0 .* change;
            end
            
            if any(b(:) > 0)
                change_les = zeros(max(b(:)), 3);
                for k = 1:size(change_les, 1)
                    %tmp = change(b.PixelIdxList{k});
                    tmp = change(b == k);
                    change_les(k,1) = sum(tmp == 1);
                    change_les(k,2) = sum(tmp == 2);
                    change_les(k,3) = sum(tmp == 3);
                end
                change_les = change_les .* volfactor;
                change_les = [change_les(change_les(:,1) == 0 & change_les(:,3) == 0,:);
                    change_les(change_les(:,1) > 0 | change_les(:,3) > 0,:)];
                joint_les = sum(change_les, 2);
                change_les = change_les(joint_les > 0.002,:);
                joint_les = joint_les(joint_les > 0.002,:);
                strout1 = ['joint = [', num2str(joint_les(1))];
                strout2 = ['decr = [', num2str(change_les(1,1))];
                strout3 = ['unch = [', num2str(change_les(1,2))];
                strout4 = ['incr = [', num2str(change_les(1,3))];
                strout5 = ['tlv1 = [', num2str(change_les(1,1) + change_les(1,2))];
                strout6 = ['tlv2 = [', num2str(change_les(1,2) + change_les(1,3))];
                for k = 2:numel(joint_les)
                    strout1 = [strout1, ',', num2str(joint_les(k))];
                    strout2 = [strout2, ',', num2str(change_les(k,1))];
                    strout3 = [strout3, ',', num2str(change_les(k,2))];
                    strout4 = [strout4, ',', num2str(change_les(k,3))];
                    strout5 = [strout5, ',', num2str(change_les(k,1) + change_les(k,2))];
                    strout6 = [strout6, ',', num2str(change_les(k,2) + change_les(k,3))];
                end
                strout1 = [strout1, '],\n']; strout2 = [strout2, '],\n'];
                strout3 = [strout3, '],\n']; strout4 = [strout4, '],\n'];
                strout5 = [strout5, '],\n']; strout6 = [strout6, '],\n'];
                
                grid = linspace(0, max([max(change_les(:,1) + change_les(:,2)), max(change_les(:,2) + change_les(:,3))]), 7);
                grid = round(grid * 100) / 100;
                %grid = 0:.4:max([max(change_les(:,1) + change_les(:,2)), max(change_les(:,2) + change_les(:,3))]);
                strgrid = num2str(grid(1));
                for k = 2:numel(grid)
                    strgrid = [strgrid, ', ', num2str(grid(k))];
                end    
            else
                strout1 = 'joint = 0';
                strout2 = 'decr = 0';
                strout3 = 'unch = 0';
                strout4 = 'incr = 0';
                strout5 = 'tlv1 = 0';
                strout6 = 'tlv2 = 0';
                strgrid = '[0,0.25,0.5,0.75,1]';
            end

            
                
            JSid = fopen(fullfile(nameFolder, ['lcp_', num2str(j), num2str(j+1), '.js']), 'wt');
            js_strout = ['$(function () {\n', ...             
              '  var min_slice', jsid, ' = ', num2str(r(1)), ',\n', ...
              '  max_slice', jsid, ' = ', num2str(r(2)), ',\n', ...
              '  slice', jsid, ' = ', num2str(round(mean(r))), ';\n', ...
              '  $( \"#slider_', jsid, '\" ).slider({\n', ...
              '    min: min_slice', jsid, ',\n', ...
              '    max: max_slice', jsid, ',\n', ...
              '    value: slice', jsid, ',\n', ...
              '    slide: function( event, ui ) {\n', ...
              '      slice', jsid, ' = ui.value;\n', ...
              '      $( \"#overlay', jsid, '\" ).attr(\"src\", \"', fullfile(cd, nameFolder, ['overlay_', num2str(j), '_', num2str(j+1), '_\" + ui.value + \".png\"']), ');\n', ...
              '      $( \"#slice', jsid, '\" ).text(\"Slice \" + ui.value);\n', ...
              '    }\n', ...
              '  });\n', ...
              '  $( \"#button-left', jsid, '\" ).button({\n', ...
              '    icons: {\n', ...
              '      primary: \"ui-icon-carat-1-w\"\n', ...
              '    },\n', ...
              '    text: false,\n', ...
              '    }).click(function(event, ui){\n', ...
              '        if(slice', jsid, ' > min_slice', jsid, '){\n', ...
              '            slice', jsid, ' = slice', jsid, ' - 1;\n', ...
              '            $( \"#overlay', jsid, '\" ).attr(\"src\", \"', fullfile(cd, nameFolder, ['overlay_', num2str(j), '_', num2str(j+1), '_\" + slice', jsid, ' + \".png\"']), ');\n', ...
              '            $( \"#slice', jsid, '\" ).text(\"Slice \" + slice', jsid, ');\n', ...
              '            $(\"#slider_', jsid, '\").slider(\"option\", "value", slice', jsid, ');\n', ...
              '        }\n', ...
              '  });\n', ...
              '  $(\"#button-right', jsid, '\").button({\n', ...
              '    icons: {\n', ...
              '      primary: \"ui-icon-carat-1-e\"\n', ...
              '    },\n', ...
              '    text: false,\n', ...
              '    }).click(function(event, ui){\n', ...
              '        if(slice', jsid, ' < max_slice', jsid, '){\n', ...
              '            slice', jsid, ' = slice', jsid, ' + 1;\n', ...
              '            $( \"#overlay', jsid, '\").attr(\"src\", \"', fullfile(cd, nameFolder, ['overlay_', num2str(j), '_', num2str(j+1), '_\" + slice', jsid, ' + \".png\"']), ');\n', ...
              '            $( \"#slice', jsid, '\").text(\"Slice \" + slice', jsid, ');\n', ...
              '            $(\"#slider_', jsid, '\").slider(\"option\", \"value\", slice', jsid, ');\n', ...
              '        }\n', ...
              '  });\n', ...
              '  var paper_', jsid, ' = Raphael(\"canvas_', jsid, '\", 475, 325),\n', ...
              '       wframe = 475,\n', ...
              '       hframe = 325,\n', ...
              '       wbox1 = 275,\n', ...
              '       hbox1 = 275,\n', ...
              '       wbox2 = 85,\n', ...
              '       m_left = 60,\n', ...
              '       m_top = 5,\n', ...
              '       m_left2 = 20,\n', ...
              '       m_top2 = 20,\n', ...
              '       ', strout5, ...
              '       ', strout6, ...
              '       ', strout1, ...
              '       ', strout4, ...
              '       ', strout3, ...
              '       ', strout2, ...                    
              '       global_change = [', num2str(tlv_joint), ', ', num2str(tlv_decr), ', ',  num2str(tlv_unch), ', ',  num2str(tlv_incr), '],\n', ...
              '       tlv1_max = Math.max.apply(null, tlv1) * 1.05,\n', ...
              '       tlv2_max = Math.max.apply(null, tlv2) * 1.05,\n', ...
              '       tlv_max = Math.max(tlv1_max, tlv2_max),\n', ...
              '       tlv_min = 0,\n', ...
              '       xgrid = [', strgrid, '],\n', ...
              '       ygrid = xgrid,\n', ...
              '       joint_max = Math.max.apply(null, joint),\n', ...
              '       w_max = 40;\n\n', ...
              '  // ******************\n', ...
              '  // Lesion change plot\n', ...
              '  // ******************\n', ...           
              '  paper_', jsid, '.rect(m_left, m_top, wbox1, hbox1).attr({fill: \"#E1E1E1\", stroke: \"none\"});\n\n', ...                
              '  // x-axis\n', ...
              '  paper_', jsid, '.setStart();\n', ...
              '  for(var i = -1; i++ < (xgrid.length-1);){\n', ...
	          '     paper_', jsid, '.path(\"M\" + ((xgrid[i] / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2) + \" \" + (m_top) + \"L\" + ((xgrid[i] / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2) + \" \" + (hbox1 + m_top)).attr({stroke: \"#fff\", \"stroke-dasharray\": \"-\"});\n', ...
              '  }\n', ...
              '  var xgrid_code = paper_', jsid, '.setFinish();\n\n', ...                
              '  // ticks\n', ...
              '  paper_', jsid, '.setStart();\n', ...
              '  for(var i = -1; i++ < (xgrid.length-1);){\n', ...
	          '     paper_', jsid, '.text(((xgrid[i] / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2), (m_top * 2 + hbox1), xgrid[i]).attr({\"font-family\": \"Courier New\", fill: \"#606060\"});\n', ...
              '  }\n', ...
              '  var xticks_code = paper_', jsid, '.setFinish();\n', ...
              '  paper_', jsid, '.text(m_left + (wbox1 / 2), hbox1+25, \"Lesion volume (ml) for t = ', num2str(j), '\").attr({\"font-family\": \"Courier New\", fill: \"#606060\"});\n\n', ...                
              '  // y-axis\n', ...
              '  paper_', jsid, '.setStart();\n', ...
              '  for(var i = -1; i++ < (xgrid.length-1);){\n', ...
	          '     paper_', jsid, '.path(\"M\" + m_left + \" \" + ((hbox1 + m_top) - ((xgrid[i] / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2)) + \"L\" + (m_left + wbox1) + \" \" +  ((hbox1 + m_top) - ((xgrid[i] / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2))).attr({stroke: \"#fff\", \"stroke-dasharray\": \"-\"});\n', ...
              '  }\n', ...
              '  var ygrid_code = paper_', jsid, '.setFinish();\n\n', ...                
              '  // ticks\n', ...
              '  paper_', jsid, '.setStart();\n', ...
              '  for(var i = -1; i++ < (xgrid.length-1);){\n', ...
	          '     paper_', jsid, '.text(m_left * .75, (hbox1 + m_top) - ((xgrid[i] / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2), xgrid[i]).attr({\"font-family\": \"Courier New\", fill: \"#606060\"});\n', ...
              '  }\n', ...
              '  var yticks_code = paper_', jsid, '.setFinish();\n', ...
              '  paper_', jsid, '.text(m_left * .4, hbox1/2, \"Lesion volume (ml) for t = ', num2str(j+1), '\").transform(\"r270\").attr({\"font-family\": \"Courier New\", fill: \"#606060\"});\n', ...                
              '  // line\n', ...
              '  var diagonal = paper_', jsid, '.path(\"M\" + ((-0 / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2 - m_left2) + \" \" + ((hbox1 + m_top) - ((0 / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2) + m_left2) + \"L\" + ((tlv_max*1.05 / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2) + \" \" + ((hbox1 + m_top) - ((tlv_max*1.05 / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2))).attr({stroke: \"#A8A8A8\"});\n\n', ...
              '  // Rectangles\n', ...
              '  for(var i = -1; i++ < (tlv1.length-1);){\n', ...
              '      if((joint[i] / joint_max) * w_max < 5){\n', ...
              '          var w_tmp = 5;\n', ...
              '      } else {\n', ...
              '          var w_tmp = (joint[i] / joint_max) * w_max;\n', ...
              '      }\n', ...
              '      paper_', jsid, '.rect((tlv1[i] / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2 - w_tmp / 2, \n', ...
              '                (hbox1 + m_top) - ((tlv2[i] / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2) + w_tmp/2 - w_tmp * (decr[i] / joint[i]),\n', ...
              '                w_tmp, w_tmp * (decr[i] / joint[i])).attr({fill: \"#00FF66\", stroke: \"none\", opacity: .5});\n', ...
              '      paper_', jsid, '.rect((tlv1[i] / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2 - w_tmp / 2, \n', ...
              '                (hbox1 + m_top) - ((tlv2[i] / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2) + w_tmp/2 - w_tmp * (decr[i] / joint[i]) - w_tmp * (unch[i] / joint[i]),\n', ...
              '                w_tmp, w_tmp * (unch[i] / joint[i])).attr({fill: \"#909090\", stroke: \"none\", opacity: .5});\n', ...
              '      paper_', jsid, '.rect((tlv1[i] / tlv_max) * (wbox1 - 2*m_left2) + m_left + m_left2 - w_tmp / 2, \n', ...
              '                (hbox1 + m_top) - ((tlv2[i] / tlv_max) * (hbox1 - 2*m_top2) + m_top + m_top2) + w_tmp/2 - w_tmp * (decr[i] / joint[i]) - w_tmp * (unch[i] / joint[i]) - w_tmp * (incr[i] / joint[i]),\n', ...
              '                w_tmp, w_tmp * (incr[i] / joint[i])).attr({fill: \"#D00000\", stroke: \"none\", opacity: .5});\n', ...
              '  }\n\n', ...                
              '  // ********\n', ...
              '  // barchart\n', ...
              '  // ********\n\n', ...              
              '  // green\n', ...
              '  paper_', jsid, '.rect(m_left + wbox1 + 10 + m_left2/2, m_top + m_top2/2, wbox2 - m_left2, hbox1 - m_top2).attr({fill: \"#00FF66\", stroke: \"#none\"});\n', ...
              '  // grey\n', ...
              '  paper_', jsid, '.rect(m_left + wbox1 + 10 + m_left2/2, m_top + m_top2/2, wbox2 - m_left2, (hbox1 - m_top2) * (global_change[3]/global_change[0])).attr({fill: \"#D00000\", stroke: \"#none\"});\n', ...
              '  // red\n', ...
              '  paper_', jsid, '.rect(m_left + wbox1 + 10 + m_left2/2, m_top + m_top2/2 + (hbox1 - m_top2) * (global_change[3]/global_change[0]), wbox2 - m_left2, (hbox1 - m_top2) * (global_change[2]/global_change[0])).attr({fill: \"#909090\", stroke: \"#none\"});\n\n', ...                
              '  // text\n', ...
              '  paper_', jsid, '.text(m_left + wbox1 + 10 + .9*m_left2 + wbox2, m_top + m_top2/2 + (hbox1 - m_top2) * (global_change[3]/global_change[0])/2, Math.round(global_change[3]*100)/100 + \" ml\").attr({\"font-family\": \"Courier New\", fill: \"#606060\"});\n', ...
              '  paper_', jsid, '.text(m_left + wbox1 + 10 + .9*m_left2 + wbox2, m_top + m_top2/2 + (hbox1 - m_top2) * (global_change[3]/global_change[0]) + (hbox1 - m_top2) * (global_change[2]/global_change[0])/2, Math.round(global_change[2]*100)/100 + \" ml\").attr({\"font-family\": \"Courier New\", fill: \"#606060\"});\n', ...
              '  paper_', jsid, '.text(m_left + wbox1 + 10 + .9*m_left2 + wbox2, m_top + m_top2/2 + (hbox1 - m_top2) * (global_change[3]/global_change[0]) + (hbox1 - m_top2) * (global_change[2]/global_change[0]) + (hbox1 - m_top2) * (global_change[1]/global_change[0])/2, Math.round(global_change[1]*100)/100 + \" ml\").attr({\"font-family\": \"Courier New\", fill: \"#606060\"});\n', ...                
            '});'];
        
            
            fprintf(JSid, js_strout);
            fclose(JSid);
            
        end
           
        strout = [strout, '\n<br><hr><br>'];
        HTMLid2 = fopen(fullfile(nameFolder, [nameFolder, '.html']), 'wt');
        fprintf(HTMLid2, strout);
        fclose(HTMLid2);
            
        fprintf(HTMLid, strout);
        fclose(HTMLid);            
        fprintf(fileID, ' ok.\n');
        
        tt = toc; tt = [num2str(round(tt)), 's'];
        strout = [repmat(' ', 1, 72 - numel(tt) - numel(stroutHTML)), tt, '\n'];
        fprintf(strout)
        fprintf(fileID, ' ok.\n');
        
    end
    
    % delete some images
    spm_unlink(fullfile(pthles{1}, ['rl', namf2{1}, '.nii']))    
    for j = 1:m
        spm_unlink(fullfile(pthles{j}, ['rl', namles{j}, '_', num2str(j), extles{j}]))
        spm_unlink(fullfile(nameFolder, ['rlp0_', num2str(j), '.nii']))
        spm_unlink(namcoreg{j})
        spm_unlink(fullfile(nameFolder, ['rl', ps_fileparts(namcoreg{j}, 2:3)]))
    end        
    cd(nameFolder)
    ls = dir;
    for k = 1:numel(ls)
        if regexp(ls(k).name, '.nii', 'once')
           spm_unlink(ls(k).name) 
        end
        if regexp(ls(k).name, '.hdr', 'once')
           spm_unlink(ls(k).name) 
        end
        if regexp(ls(k).name, '.img', 'once')
           spm_unlink(ls(k).name) 
        end
    end
    
end

fclose(fileID);
cd(pthor)
spm_unlink(nameLog);
    
c = clock();
strout = 'Finished successfully ';
fprintf(strout)
tt = datestr(c);
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n\n'];
fprintf(strout)
fprintf(repmat('-', 1, 72));
fprintf('\n')

varargout{:} = 'Don''t forget to cite the toolbox.';%['Finished successfully on ', datestr(c)];