function varargout = ps_LST_lga(varargin)
%ps_LST_lga   Lesion segmentation by a lesion growth algorithm (LGA)
%   Part of the LST toolbox, www.statistical-modelling.de/lst.html
%
%   ps_LST_lga Lesion segmentation by the LGA requires a T1 and a FLAIR 
%   image. Furthermore, the user has to specify an initial threshold
%   (kappa). See Schmidt et al. (2012) for details.
%
%   This routine produces lesion probability maps (ples...), coregistered
%   bias corrected versions of the FLAIR inputs, a .mat-file with
%   information about the segmentation that is needed for a re-run of the 
%   algorithm, and a HTML report along with a folder if this option has 
%   been chosen desired.
%   
%   ps_LST_lga asks the user for the input images (T1 and FLAIR) and sets
%   kappa to its default value (0.3). MRF parameter is set to 1 and, number 
%   of maximum iterations are 50 and a HTML report is produced.
%
%   ps_LST_lga(Vt1, Vf2, kappa, phi, mxiter, html) performs lesion 
%   segmentation for the image volumes given in Vt1 and Vf2. Both must be 
%   characters like a call from from spm_select. Initial thresholds are 
%   colectedInitial in kappa, MRF parameter in phi, number of maximum 
%   iterations in maxiter and a dummy for the HTML report in html. If the
%   last four arguments are missing they are replaced by their default
%   values, see above.
%
%   ps_LST_lga(job) Same as above but with job being a harvested job data
%   structure (see matlabbatch help).
%
%   References
%   P. Schmidt, Gaser C., Arsic M., Buck D., F?rschler A., Berthele A., 
%   Hoshi M., Ilg R., Schmid V.J., Zimmer C., Hemmer B., and M?hlau M. An 
%   automated tool for detec- tion of FLAIR-hyperintense white-matter 
%   lesions in Multiple Sclerosis. NeuroImage, 59:3774?3783, 2012.
%

%addpath(fullfile('/Users/paul/Software/spm8','toolbox','LST'));

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
tt = 'Lesion growth algorithm (LGA)\n';
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
if ~isempty(varargin) && isfield(varargin{1}, 'data_T1')
    viajob = 1;
else
    viajob = 0;
end
if ~viajob
    if nargin == 0
        % call to ps_LST_lga
        fprintf(fileID, 'Select files by function ... ');
        Vt1 = spm_select(Inf, 'image', 'Select T1 images.');
        Vf2 = spm_select(Inf, 'image', 'Select FLAIR images.');
        fprintf(fileID, 'ok.\n');    
        kappa = 0.3;
        maxiter = 50;
        phi = 1;
        html_report = 1;
    end
    if nargin == 1
        fprintf(fileID, 'Only one argument.');
        fprintf('Please give me at least two arguments.\n');
        fprintf('See ?ps_LST_lga for help.\n');
        fclose(fileID);
        return;        
    end
    if nargin > 1
        if isempty(varargin{1})
            fprintf(fileID, 'Select files by function ... ');
            Vt1 = spm_select(Inf, 'image', 'Select T1 images.');
            fprintf(fileID, 'ok.\n'); 
        else
            if ischar(varargin{1})
                Vt1 = varargin{1};
            else
                fprintf(fileID, 'Vt1 is not a character!');
                fprintf('Input for Vt1 must be a character, like from spm_select.\n');
                fclose(fileID);
                return;
            end
        end
        if isempty(varargin{2})
            fprintf(fileID, 'Select files by function ... ');
            Vf2 = spm_select(Inf, 'image', 'Select FLAIR images.');
            fprintf(fileID, 'ok.\n'); 
        else
            if ischar(varargin{2})
                Vf2 = varargin{2};
            else
                fprintf(fileID, 'Vf2 is not a character!');
                fprintf('Input for Vf2 must be a character, like from spm_select.\n');
                fclose(fileID);
                return;
            end
        end
        if nargin > 2
            if isempty(varargin{3})
                kappa = 0.3;
            else
                kappa = varargin{3};
            end
        else
            kappa = 0.3;
        end
        if nargin > 3
            if isempty(varargin{4})
                maxiter = 50;
            else                
                maxiter = varargin{4};
            end
        else
            maxiter = 50;
        end
        if nargin > 4
            if isempty(varargin{5})
                phi = 1;
            else                
                phi = varargin{5};
            end
        else
            phi = 1;
        end
        if nargin > 5
            if isempty(varargin{6})
                html_report = 1;
            else                
                html_report = varargin{6};
            end
        else
            html_report = 1;
        end
    end    
else    
    job = varargin{1};    
    Vt1 = job.data_T1;
    Vf2 = job.data_F2;
    kappa = job.opts_lga.initial;
    maxiter = job.opts_lga.maxiter;
    phi = job.opts_lga.mrf;
    html_report = job.html_report;
	xasl_quality = job.xasl_quality;
end

if ~exist('xasl_quality','var')
	xasl_quality = 1;
end

fprintf(fileID, 'Load volumes ... ');
Vt1 = spm_vol(Vt1);
Vf2 = spm_vol(Vf2);
fprintf(fileID, 'ok.\n');

if ~isequal(numel(Vt1), numel(Vf2))
    fprintf(fileID, 'numel(T1) != numel(FLAIR)');
    fclose(fileID);
    error('Number of T1 images must match number of FLAIR images.\n'); 
end
if numel(Vt1) == 0 || numel(Vf2) == 0
    fprintf(fileID, 'No images selected.');
    fprintf('No images selected.\n');
    fclose(fileID);
    return;
end

% Summarize input
fprintf(fileID, 'Input summary:\n');
fprintf(fileID, ['Jobs: ', num2str(numel(Vf2)), '\nkappa: ', num2str(kappa(kappa > 0)), ...
    '\nMaxiter: ', num2str(maxiter), '\nPhi: ', num2str(phi)]);
strout = 'Number of jobs to process:';
fprintf(strout)
tt = [num2str(numel(Vf2)), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)
if numel(kappa(kappa > 0)) == 1
    strout = 'One kappa value: ';
    tt = num2str(kappa(kappa > 0));
else
    strout = [num2str(numel(kappa(kappa > 0))), ' different kappa values: '];
    if numel(kappa(kappa > 0)) > 3
        tt = ['ranging from ', num2str(min(kappa(kappa > 0))), ' to ', num2str(max(kappa(kappa > 0)))];        
    else
        tmp = kappa(kappa > 0);
        tt = '(';
        for j = 1:(numel(tmp) - 1)
           tt = [tt, num2str(tmp(j)), ', ']; 
        end
        tt = [tt, num2str(tmp(j+1)), ')'];
    end
end
fprintf(strout)
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
fprintf(strout)
strout = 'Maxiter:';
fprintf(strout)
tt = [num2str(maxiter), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)
strout = 'Phi:';
fprintf(strout)
tt = [num2str(phi), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)

% Loop over all subjects
for i = 1:numel(Vt1)
    
    fprintf(fileID, '\n--------------\n');
    fprintf(fileID, ['Job ', num2str(i), ' of ', num2str(numel(Vt1)), '\n']);
    fprintf(fileID, '--------------\n');
    
    % Extract file information
    if viajob
        Vt1_tmp = Vt1{i};
        Vf2_tmp = Vf2{i};    
    else    
        Vt1_tmp = Vt1(i);
        Vf2_tmp = Vf2(i);
    end
    [ptht1, namt1, extt1] = fileparts(Vt1_tmp.fname);
    [pthf2, namf2, extf2] = fileparts(Vf2_tmp.fname);
    cd(ptht1)
        
    % Which subject?    
    strout = '\nWorking on job';
    fprintf(strout)
    tt = [num2str(i), ' out of ', num2str(numel(Vf2)), ' (', num2str(i/numel(Vf2)*100), '%%)\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 5), tt];
    fprintf(strout)
        
    % Correct directory?
    cd_tmp = cd;
    cd_tmp  = xASL_adm_ConvertSlash(cd_tmp,1); %% ExploreASL bugfix slashes
    strout = 'Current directory ';
    if numel(cd_tmp) > (72 - numel(strout) - 5)
        fprintf([strout, '...', cd_tmp((end - (72 - numel(strout) - 4)):end), '\n'])
    else
        fprintf([strout, repmat(' ', 1, 72 - numel(strout) - numel(cd_tmp)), cd_tmp, '\n'])
    end
    fprintf(fileID, ['Current directory is ', cd_tmp, '\n']);
    
    % Check if user has writing permissions
    [~, struc] = fileattrib;
    if ~struc.UserWrite 
        fprintf(fileID, 'User has no writing permissions!');
        error('You do not have writing permissions for the current folder.\n'); 
    end
    
    % Check if results from a previous run are avaliable
    if exist(['LST_lga_rm', namf2, '.mat'], 'file')
        fprintf(fileID, 'Found results from a previous run.\n');
        fprintf('Load results from a previous run.\n');
        pr = 1;
        load(['LST_lga_rm', namf2, '.mat'])
        indx_brain = lga.indx_brain;
        f2 = zeros(lga.dim); f2(:) = lga.f2_vec;
        p0 = zeros(lga.dim); p0(indx_brain) = lga.p0_vec;
        atlas_wm = zeros(lga.dim); atlas_wm(indx_brain) = lga.atlas_wm_vec;
        nx = lga.dim(1); ny = lga.dim(2); nz = lga.dim(3);
        noles = 0 .* f2; noles(indx_brain) = lga.noles_vec;
        or = lga.or;
        fl = lga.fl;
        %seg = zeros(lga.dim); seg(indx_brain) = lga.I;        
    else
        pr = 0;
    end
    
    if ~pr
        
        % Create temporary folder
        % ID
        id = [num2str(round(rand(1, 1) * 1e4)), '_', num2str(round(rand(1, 1) * 1e4)), '_', num2str(round(rand(1, 1) * 1e4))];        
        tmpFolder = ['LST_tmp_', id];
        if ~isdir(tmpFolder)  % ExploreASL bugfix
            mkdir(tmpFolder);
        end

        % Preprocessing for T1
        % --------------------

        
        
        % get job defaults, ExploreASL hack to speed up for reproducibility tests
        %if  exist(fullfile(cd_tmp,'LowQ.status'),'file')
		if xasl_quality
			[job1, job2, job3] = ps_LST_lga_preproc_default; % normal, high quality
		else
			[job1, job2, job3] = ps_LST_lga_preproc_default_LowQ; % low quality for fast repro testing
		end
        
        if ~(exist(['m', namt1, '.nii'], 'file') && exist(['c1', namt1, '.nii'], 'file') ...
                && exist(['c2', namt1, '.nii'], 'file') && exist(['c3', namt1, '.nii'], 'file') ...
                && exist(['iy_', namt1, '.nii'], 'file'))
            fprintf(fileID, ['Preprocessing of T1 (', namt1, extt1, ') ...']);
            strout = ['Preprocessing of ', namt1, extt1, ' '];
            fprintf(strout)
            tic        
            % Copy t1 into temporary folder
            copyfile(Vt1_tmp.fname, fullfile(tmpFolder, '.'))
            if strcmp(extt1, '.img')
                copyfile([namt1, '.hdr'], fullfile(tmpFolder, '.'))
            end
            job1.channel.vols = {fullfile(tmpFolder, [namt1, extt1])};
            ps_LST_spm_preproc_run(job1);
            tt = toc; tt = [num2str(round(tt)), 's'];
            strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
            fprintf(strout)
            %t1del = 1;
            fprintf(fileID, ' ok.\n');
        else
            fprintf(fileID, 'Skipped preprocessing for T1.\n');
            fprintf('All images of T1 preprocessing exist\n');
            copyfile(['m', namt1, '.nii'], fullfile(tmpFolder, '.'))
            copyfile(['c1', namt1, '.nii'], fullfile(tmpFolder, '.'))
            copyfile(['c2', namt1, '.nii'], fullfile(tmpFolder, '.'))
            copyfile(['c3', namt1, '.nii'], fullfile(tmpFolder, '.'))
            copyfile(['iy_', namt1, '.nii'], fullfile(tmpFolder, '.'))
            %t1del = 0;
        end

        % Preprocessing for FLAIR
        % -----------------------

        %if ~(exist(['m', namf2, '.nii'], 'file') || exist(['rm', namf2, '.nii'], 'file'))
            fprintf(fileID, ['Preprocessing of FLAIR (', namf2, extf2, ') ...']);
            strout = ['Preprocessing of ', namf2, extf2, ' '];
            fprintf(strout)
            tic
            % Copy FLAIR into temporary folder
            copyfile(Vf2_tmp.fname, fullfile(tmpFolder, '.'))
            if strcmp(extf2, '.img')
                copyfile([namf2, '.hdr'], fullfile(tmpFolder, '.'))
            end
            job2.channel.vols = {fullfile(tmpFolder, [namf2, extf2])};
            ps_LST_spm_preproc_run(job2);
            tt = toc; tt = [num2str(round(tt)), 's'];
            strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
            fprintf(strout)
            f2del = 1;
            fprintf(fileID, ' ok.\n');
        %else
        %    fprintf(fileID, 'Skipped preprocessing for FLAIR.\n');
        %    fprintf('Bias corrected FLAIR image exists\n');
        %    f2del = 0;
        %end


        % Coregister FLAIR to T1
        % ----------------------
	%% ExploreASL Hack: commented lines 360-371, to disable coregistration (already performed by ExploreASL)
        %if ~exist(['rm', namf2, '.nii'], 'file')
%%            fprintf(fileID, 'Coregistration ...');
%%            strout = 'Coregister FLAIR to T1 ';
%%            fprintf(strout)
%%            tic            
%%            %job.ref = {Vref{i}.fname};
%%            job3.ref = {fullfile(tmpFolder, ['m', namt1, '.nii'])};
%%            job3.source = {fullfile(tmpFolder, ['m', namf2, '.nii'])};
%%            ps_LST_spm_run_coreg(job3);
%%            tt = toc; tt = [num2str(round(tt)), 's'];
%%            strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
%%            fprintf(strout)
%%            fprintf(fileID, ' ok.\n');
            copyfile(fullfile(tmpFolder, [ 'm', namf2, '.nii']), fullfile(tmpFolder, [ 'rm', namf2, '.nii']));
            copyfile(fullfile(tmpFolder, ['rm', namf2, '.nii']), '.');
        %else
        %    fprintf(fileID, 'Skipped coregistration.\n');
        %    fprintf('Coregistered bias corrected FLAIR image exists\n');
        %end

        %if f2del
        %    spm_unlink(['m', namf2, '.nii']);
        %end
    
    
        % Initialization
        % --------------

        strout = 'Initialize lesions ';
        fprintf(strout)
        tic

        % 1) PVE label
        fprintf(fileID, 'Load tissue probabilities ...');
        % load tissue probabilities
        p1 = spm_read_vols(spm_vol(fullfile(tmpFolder, ['c1', namt1, '.nii'])));
        p2 = spm_read_vols(spm_vol(fullfile(tmpFolder, ['c2', namt1, '.nii'])));
        p3 = spm_read_vols(spm_vol(fullfile(tmpFolder, ['c3', namt1, '.nii'])));
        t1 = spm_read_vols(spm_vol(fullfile(tmpFolder, ['m', namt1, '.nii'])));
        %if t1del
        %    spm_unlink(['c1', namt1, '.nii'])
        %    spm_unlink(['c2', namt1, '.nii'])
        %    spm_unlink(['c3', namt1, '.nii'])
        %    spm_unlink(['m', namt1, '.nii'])
        %end
        fprintf(fileID, ' ok.\n');

        fprintf(fileID, 'Hard segmentation ...');
        prob = [p3(:), p1(:), p2(:)];        
        prob = bsxfun(@times, prob, 1 ./ sum(prob, 2));
        indx_brain = find(sum(prob, 2) > 0);
        prob = prob(indx_brain,:);                        

        % Hard segmentation
        seg = 0 .* t1;
        [~, I] = max(prob, [], 2);
        seg(indx_brain) = I;
        seg(:,:,1) = 0 .* seg(:,:,1); seg(:,:,end) = 0 .* seg(:,:,end);
        seg(:,1,:) = 0 .* seg(:,1,:); seg(:,end,:) = 0 .* seg(:,end,:);
        seg(1,:,:) = 0 .* seg(1,:,:); seg(end,:,:) = 0 .* seg(end,:,:);
        clear I;
        fprintf(fileID, ' ok.\n');
        
        % No les
        fprintf(fileID, 'No lesion voxel ...');
        seg2 = seg; 
        seg2(1,:,:) = 1; seg2(end,:,:) = 1;
        seg2(:,1,:) = 1; seg2(:,end,:) = 1;
        seg2(:,:,1) = 1; seg2(:,:,end) = 1;
        if exist('bwlabeln', 'builtin')
            b = bwlabeln(1 .* (seg2 == 0));
        else
            b = ps_bwlabeln(1 .* (seg2 == 0));
        end
        clear seg2;
        c_tmp = ps_count(b(b > 0));
        indx_tmp = find(b == c_tmp(1,find(c_tmp(2,:) == max(c_tmp(2,:)))));
        tmp = 1 .* (seg == 1); tmp(indx_tmp) = NaN;
        for j = 1:5
            indx_tmp = find(seg == 1 & ~isnan(tmp));
            nh = getNeighborhood2(tmp, indx_tmp, 3);
            tmp(indx_tmp(sum(isnan(nh) > 0) > 0)) = NaN;
        end
        noles = 1 .* (isnan(tmp) & seg > 0);
        clear tmp;
        fprintf(fileID, ' ok.\n');
        
        % Generate PVE label
        fprintf(fileID, 'PVE label ...');
        m = [mean(t1(seg == 1)), mean(t1(seg == 2)), mean(t1(seg == 3))];
        p0 = 0 .* t1;
        p0(seg == 1) = (-1) .* ps_scale(p3(seg == 1), -1.5, -1);
        p0(seg == 2 & t1 < m(2)) = ps_scale(p1(seg == 2 & t1 < m(2)), 1.5, 2);
        p0(seg == 2 & t1 > m(2)) = (-1) .* ps_scale(p1(seg == 2 & t1 > m(2)), -2.5, -2);
        p0(seg == 3) = ps_scale(p2(seg == 3), 2.5, 3);

        f2 = spm_read_vols(spm_vol(['rm', namf2, '.nii']));    
        f2(isnan(f2)) = 0;
        indx_tmp = find(p0 < 1.5 & p0 > 0 & f2 > 1 .* mean(f2(seg == 2)));
        p0(indx_tmp) = 2.4;

        tmp = 0 .* t1; tmp(indx_tmp) = 1;
        st = 0;
        counter = 0;
        while ~st
            counter = counter + 1;
            indx_tmp_tmp = find(seg == 2 & tmp < 1);
            nh = getNeighborhood2(tmp, indx_tmp_tmp, 1);
            indx_tmp_tmp_tmp = indx_tmp_tmp(sum(nh > 0) > 0 & f2(indx_tmp_tmp)' > ps_quantile(f2(seg == 2), .95));
            if isempty(indx_tmp_tmp_tmp) || counter > 50
                st = 1;
                if counter > 50
                   fprintf(fileID, 'Counter exceded 50 iterations ...'); 
                end
            else
                tmp(indx_tmp_tmp_tmp) = 1;
            end
        end
        p0(tmp > 0) = 2.4;
        p0(:,:,1) = 0 .* p0(:,:,1); p0(:,:,end) = 0 .* p0(:,:,end);
        p0(:,1,:) = 0 .* p0(:,1,:); p0(:,end,:) = 0 .* p0(:,end,:);
        p0(1,:,:) = 0 .* p0(1,:,:); p0(end,:,:) = 0 .* p0(end,:,:);
        indx_brain = find(p0 > 0);
        fprintf(fileID, [' ok with ', num2str(counter), ' iterations.\n']);

        % 2) Warp WM atlas in native space
        % copy WM TPM to image folder    
        fprintf(fileID, 'Copy atlas_wm.nii ...');
        copyfile(fullfile(spm('dir'), 'toolbox', 'LST', 'atlas_wm.nii'), fullfile(tmpFolder, '.'))
        fprintf(fileID, ' ok.\n');

        % Apply inverse deformation field
        fprintf(fileID, 'Warp atlas_wm.nii to native space ...');
        cd(tmpFolder)
        clear job    
        job.comp{1}.def = {['iy_', namt1, '.nii']};
        job.out{1}.pull.fnames = {'atlas_wm.nii'};
        job.out{1}.pull.savedir.savepwd = 1;
        job.out{1}.pull.interp = 0;
        job.out{1}.pull.mask = 1;
        job.out{1}.pull.fwhm = [0 0 0];
        spm_deformations(job);
        
        load(fullfile(spm('dir'), 'toolbox', 'LST', 'LST_lpa_stuff.mat'), 'noles', 'bp_mni')
        Vatlas = spm_vol('atlas_wm.nii');
        Vatlas.fname = 'noles.nii';
        tmp = zeros(121, 145, 121);
        tmp(bp_mni) = noles;
        spm_write_vol(Vatlas, tmp);
        clear job    
        job.comp{1}.def = {['iy_', namt1, '.nii']};
        job.out{1}.pull.fnames = {'noles.nii'};
        job.out{1}.pull.savedir.savepwd = 1;
        job.out{1}.pull.interp = 0;
        job.out{1}.pull.mask = 1;
        job.out{1}.pull.fwhm = [0 0 0];
        spm_deformations(job);
        noles = spm_read_vols(spm_vol('wnoles.nii'));
        
        cd ..
        fprintf(fileID, ' ok.\n');

        % Load images    
        fprintf(fileID, 'Load watlas_wm.nii ...');
        atlas_wm = spm_read_vols(spm_vol(fullfile(tmpFolder, 'watlas_wm.nii'))); 
        atlas_wm(noles > 0) = 0;
        %spm_unlink('atlas_wm.nii');
        %spm_unlink('watlas_wm.nii');
        %if t1del
        %    spm_unlink(['iy_', namt1, '.nii']);
        %end
        atlas_wm(isnan(atlas_wm)) = 0;
        atlas_wm = atlas_wm .* (p0 > 0);
        atlas_wm(atlas_wm < 0) = 0;
        nx = size(f2, 1); ny = size(f2, 2); nz = size(f2, 3);
        fprintf(fileID, ' ok.\n');
        
        % Save all relevant information for lesion filling or rerun
        lga.indx_brain = indx_brain;
        lga.f2_vec = f2(:);
        lga.p0_vec = p0(indx_brain);
        lga.atlas_wm_vec = atlas_wm(indx_brain);
        lga.I = seg(indx_brain);
        lga.dim = [nx, ny, nz];
        lga.Vt1 = Vt1_tmp;
        lga.Vf2 = Vf2_tmp;
        lga.noles_vec = noles(indx_brain);
        % Should we flip the images?
        V = spm_vol(fullfile(spm('dir'), 'toolbox', 'LST', 'atlas_wm.nii'));
        or_ch = zeros(V.dim);
        for j = 1:121
            or_ch(60,70,j) = 1;
            or_ch(j,70,60) = 2;
        end
        for j = 1:72
            or_ch(60,j,60) = 3;
            or_ch(60,146-j,60) = 4;
        end
        cd(tmpFolder)
        V.fname = 'or_ch.nii';
        spm_write_vol(V, or_ch);    
        clear job    
        job.comp{1}.def = {['iy_', namt1, '.nii']};
        job.out{1}.pull.fnames = {'or_ch.nii'};
        job.out{1}.pull.savedir.savepwd = 1;
        job.out{1}.pull.interp = 0;
        job.out{1}.pull.mask = 1;
        job.out{1}.pull.fwhm = [0 0 0];
        spm_deformations(job);    
        or_ch = spm_read_vols(spm_vol('wor_ch.nii'));
        spm_unlink('wor_ch.nii');
        spm_unlink('or_ch.nii');
        or_ch(isnan(or_ch)) = 0;
        cd ..
        % Obtain voxel size
        vs = zeros(1, 3);
        point1 = cor2mni([1 1 1], Vf2_tmp.mat);
        % z-direction
        point2 = cor2mni([1 1 2], Vf2_tmp.mat);
        vs(3) = sqrt(sum((point2 - point1).^2));
        % y-direction    
        point2 = cor2mni([1 2 1], Vf2_tmp.mat);
        vs(2) = sqrt(sum((point2 - point1).^2));
        % x-direction    
        point2 = cor2mni([2 1 1], Vf2_tmp.mat);
        vs(1) = sqrt(sum((point2 - point1).^2));
        c_tmp = indx2coord(find(or_ch == 1), size(or_ch, 1), size(or_ch, 2));
        %z = find(ps_range(c_tmp) == max(ps_range(c_tmp)));
        z = find(ps_range(c_tmp).*vs == max(ps_range(c_tmp).*vs));
        c_tmp = indx2coord(find(or_ch == 2), size(or_ch, 1), size(or_ch, 2));    
        %x = find(ps_range(c_tmp) == max(ps_range(c_tmp)));   
        x = find(ps_range(c_tmp).*vs == max(ps_range(c_tmp).*vs));   
        c_tmp = indx2coord(find(or_ch == 3 | or_ch == 4), size(or_ch, 1), size(or_ch, 2));    
        %y = find(ps_range(c_tmp) == max(ps_range(c_tmp)));    
        y = find(ps_range(c_tmp).*vs == max(ps_range(c_tmp).*vs));
        or = [x y z];
        if any(~ismember(1:3, or))            
            or = [1, 2, 3];            
        end
        % Flip?
        or_ch = permute(or_ch, or);
        c_tmp3 = indx2coord(find(or_ch == 3), size(or_ch, 1), size(or_ch, 2));
        c_tmp4 = indx2coord(find(or_ch == 4), size(or_ch, 1), size(or_ch, 2));
        fl =  max(c_tmp3(:,2)) > max(c_tmp4(:,2));    
        lga.or = or;
        lga.fl = fl;
        save(['LST_lga_rm', namf2, '.mat'], 'lga');
        
        % delete temporary files
        rmdir(tmpFolder, 's')
        
        %spm_unlink([namt1, '_seg8.mat']);
        %spm_unlink([namf2, '_seg8.mat']);
    end % END if(~pr)
    
    if pr
        strout = 'Initialize lesions ';
        fprintf(strout)
        tic
    end
    fprintf(fileID, 'Scale FLAIR image ...');
    if 0
        f2_norm = f2 ./ mean(f2(p0 > 0));
        f2_norm = f2_norm .* (p0 > 0);
    else
        %[f, x] = ksdensity(f2(p0 > 0), 0:max(f2(p0 > 0)));
        x = 0:max(f2(p0 > 0 & f2 > 0));
        f = histc(f2(p0 > 0 & f2 > 0), x);
        xmax = x(f == max(f));
        f2_norm = f2 ./ mean(xmax);
    end
    clear f2;
    fprintf(fileID, ' ok.\n');

    % 3) Lesion belief maps    
    % lesion belief map for GM
    fprintf(fileID, 'Calculate lesion belief maps ...');
    mean_gm = mean(f2_norm((p0 > 1.5) & (p0 < 2.5)));
    B_gm = p0 .* (p0 > 1.5 & p0 < 2.5) .* (f2_norm - mean_gm) .* atlas_wm;
    B_gm = B_gm.*(B_gm > 0);
    B_gm(noles > 0) = 0;

    % lesion belief map for CSF
    mean_csf = mean(f2_norm(p0 <= 1.5 & p0 > 0));
    B_csf = p0 .* (p0 < 1.5 & p0 > 0) .* (f2_norm - mean_csf) .* atlas_wm;
    B_csf(B_csf < 0) = 0;
    B_csf(noles > 0) = 0;

    % lesion belief map for WM
    mean_wm = mean(f2_norm(p0 > 2.5));
    B_wm = p0 .* (p0 > 2.5) .* (f2_norm - mean_wm) .* atlas_wm;
    B_wm(B_wm < 0) = 0;
    B_wm(noles > 0) = 0;
    clear atlas_wm;
    
    % complete lesion belief map
    B = B_gm + B_wm + B_csf;
    B_init = B_gm;
    fprintf(fileID, ' ok.\n');
    
    % smooth lesion belief maps with a simple mean filter
    fprintf(fileID, 'Smooth lesion belief maps ...');
    indx_les = find(B_init > 0);
    neighborhood = getNeighborhood2(B_init, indx_les, 1);
    neighborhood_mean = mean(neighborhood, 1);
    B_init_mean = B_init .* 0;
    B_init_mean(indx_les) = neighborhood_mean;
    B_init_mean = B_init_mean .* B_init;
    clear B_gm; clear B_wm; clear B_csf;
    
    indx_les = find(B > 0);
    neighborhood = getNeighborhood2(B, indx_les, 1);
    neighborhood_mean = mean(neighborhood, 1);
    B_mean = B .* 0;
    B_mean(indx_les) = neighborhood_mean;
    clear [B, neighborhood_mean, neighborhood];
    fprintf(fileID, ' ok.\n');
    
    % create independence structure
    fprintf(fileID, 'Create independence structure ...');
    indi_struct_brain = uint8(createIndependenceStructure(nx, ny, nz, 1) .* (f2_norm > 0));
    fprintf(fileID, ' ok.\n');
    
    tt = toc; tt = [num2str(round(tt)), 's'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
    fprintf(strout)
    
    % Lesion growing
    % --------------    
        
    for kappa_tmp = kappa(kappa > 0)
        
        fprintf(fileID, ['Grow lesions (kappa = ', num2str(kappa_tmp), ') ']);
        strout = ['Grow lesions (kappa = ', num2str(kappa_tmp), ') '];
        fprintf(strout)
        tic                
        % binary initial lesion map
        Lesion_init = 1 .* (B_init_mean > kappa_tmp);

        % delete initial voxels that have no neighboring Lesion        
        neighborhood = getNeighborhood2(Lesion_init, indx_brain, 1);
        neighborhood_img = B_mean .* 0;
        neighborhood_img(indx_brain) = sum(neighborhood > 0, 1);
        Lesion_init(neighborhood_img < 2) = 0;
        clear [neighborhood_img, neighborhood];
        iter = uint8(0);
        
        if sum(Lesion_init(indx_brain) > 0) % check if any lesions are initialized
            
            indx = uint32(find(B_mean > 0));
            Lesion_iter = Lesion_init;
            clear Lesion_init;            
            max_lesion_iter = 1;            

            while iter < maxiter && max_lesion_iter > 0.01

                % display progress
                iter = iter + 1;                

                % fit a gamma distribution to the lesion class
                [a,b] = ps_LST_fitgamma(double(f2_norm(Lesion_iter > 0.5)));
                % and a gaussian mixture model to the other classes
                dens_csfgmwm = ps_LST_calc_mixture(p0, f2_norm, Lesion_iter);

                % get the voxels that are neighbors to lesion voxels
                neighborhood = getNeighborhood2(Lesion_iter, indx, 1);
                neighborhood_img = B_mean .* 0;
                neighborhood_img(indx) = sum(neighborhood, 1);
                neighborhood_img(Lesion_iter > 0) = 0;
                indx_tmp1 = find(neighborhood_img > 0 & indi_struct_brain == 1);
                indx_tmp2 = find(neighborhood_img > 0 & indi_struct_brain == 2);
                indx_tmp = find(neighborhood_img > 0);

                norm_les = 0 .* B_mean;
                norm_les(indx_tmp) = ps_LST_dgamma(f2_norm(indx_tmp), a, b);
                prob = (B_mean .* norm_les) ./ dens_csfgmwm;

                % independence stucture = 1
                neighborhood = getNeighborhood2(Lesion_iter, indx_tmp1, 1);
                neighborhood = sum(neighborhood, 1);

                mf_img = 0 .* B_mean;
                mf_img(indx_tmp1) = exp(- phi .* (6 - neighborhood) + phi .* neighborhood);
                Lesion_iter = Lesion_iter + prob .* mf_img .* (neighborhood_img > 1);
                Lesion_iter(Lesion_iter > 1) = 1;

                % independence stucture = 2
                neighborhood = getNeighborhood2(Lesion_iter, indx_tmp2, 1);
                neighborhood = sum(neighborhood, 1);

                mf_img = 0 .* B_mean;
                mf_img(indx_tmp2) = exp(- phi .* (6 - neighborhood) + phi .* neighborhood);
                Lesion_iter = Lesion_iter + prob .* mf_img .* (neighborhood_img > 1);
                Lesion_iter(Lesion_iter > 1) = 1;

                max_lesion_iter = max([Lesion_iter(indx_tmp1); Lesion_iter(indx_tmp2)]);
                if isempty(max_lesion_iter)
                    max_lesion_iter = 0;
                end

            end
            tt = toc; tt = ['finished after ', num2str(iter), ' iterations, ', num2str(round(tt)), 's'];
            strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
            fprintf(strout)
            fprintf(fileID, [' ok, finished after ', num2str(iter), ' iterations.\n']);
            
            % Which voxels are surrounded by lesions only?            
            neighborhood = getNeighborhood2(Lesion_iter, indx_brain, 3);
            neighborhood_img = zeros(nx, ny, nz);
            neighborhood_img(indx_brain) = sum(neighborhood > 0.01, 1);
            Lesion_iter(Lesion_iter < 1 & neighborhood_img > 18) = 1;

            % Delete all voxels that do not have a neighbor in their first order
            % neighborhood            
            neighborhood = getNeighborhood2(Lesion_iter, indx_brain, 1);
            neighborhood_img = zeros(nx, ny, nz);
            neighborhood_img(indx_brain) = sum(neighborhood > 0.01, 1);
            Lesion_iter(neighborhood_img == 0) = 0;

            clear neighborhood_img;
            clear norm_les;
            clear prob;
            clear mf_img;
            
        else 
            tt = toc; tt = ['finished after ', num2str(iter), ' iterations, ', num2str(round(tt)), 's'];
            strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
            fprintf(strout)
            fprintf(fileID, ' ok, no lesions found.\n');
            Lesion_iter = Lesion_init;
        end % END sum(Lesion_init(indx) > 0)
        Lesion_iter(isnan(Lesion_iter)) = 0;
        
        fprintf(fileID, 'Write results ...');
        Vles = spm_vol(['rm', namf2, '.nii']);
        Vles.fname = ['ples_lga_', num2str(kappa_tmp), '_rm', namf2, '.nii'];
        Vles.descrip = 'Probability lesion map obtained by LGA within LST toolbox';
        spm_write_vol(Vles, Lesion_iter);
        fprintf(fileID, ' ok.\n');
        
        if html_report
            
            % HTML report
            % -------------------------------------------------------------
            
            stroutHTML = 'Create HTML report';
            fprintf(stroutHTML);
            tic
    
            % create HTML report
            nameFolder = ['LST_lga_', num2str(kappa_tmp), '_rm', namf2];
            warning('off');
            mkdir(nameFolder)
            warning('on');
            
            % Create PNGs
            fprintf(fileID, 'Create PNGs ...');
            pngFailed = '';
            try
                [~, r] = ps_LST_create_gif(fullfile(cd, ['rm', namf2, '.nii']), ...
                    Vles.fname, nameFolder, 0, [or, fl]);
                fprintf(fileID, ' ok.\n');
            catch ME
                fprintf(fileID, ' failed!.\n');             
                pngFailed = ME.message;
                r = 0:1;
            end
            
            % Create images for glass brains
            ps_LST_create_glass_brain(1 .* (Lesion_iter > .5), ...
                1 .* (p0 > 0), ...
                fullfile(nameFolder, ['gb']), [or, fl]);
            
            fprintf(fileID, 'Create main HTML file ...');
            % Main HTML file#
            nameHTML = ['report_LST_lga_', num2str(kappa_tmp), '_rm', namf2, '.html'];
            copyfile(fullfile(spm('dir'), 'toolbox', 'LST', 'LST_main_html.html'), fullfile(cd, nameHTML))
            HTMLid = fopen(nameHTML, 'at');
            strout = ['    <script src=\"', ps_fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'raphael.js'), '\"></script>\n', ...
                  '    <script src=\"', ps_fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery.min.js'), '\"></script>\n', ...
                  '    <link href=\"', ps_fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.css'), '\" rel=\"stylesheet\"></script>\n', ...
                  '    <script src=\"', ps_fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.js'), '\"></script>\n', ...
                  '  </head>\n  <body>\n'];
            fprintf(HTMLid, strout);
            fprintf(fileID, ' ok.\n');
            
            % create subject specific html file
            fprintf(fileID, 'Create subject specific HTML file ...');
            volfactor = abs(det(Vles.mat(1:3,1:3))) /  1000;
            
            if any(Lesion_iter(:) > .5)
                bw = ps_bwlabeln(1. * (Lesion_iter > .5));
                les_size = zeros(max(bw(:)), 1);
                for j = 1:max(bw(:))
                    les_size(j) = sum(bw(:) == j) * volfactor;
                end
                les_size = les_size(les_size > 0.015);
            else
                les_size = [];
            end 

            tlv = sum(les_size);        
            numles = numel(les_size);            
            jsid = [nameFolder, '_', ps_create_timestamp];
            jsid(regexp(jsid, '\.')) = [];
                        
            strout = ['\n<div class=\"container\">\n', ...              
              '  <h1>Lesion segmentation by LST</h1>\n', ...
              '  <div class=\"column-01\">\n', ...
              '    <h2>Input summary</h2>\n', ...
              '      <table style=\"min-width: 500px;\">\n', ...
              '        <tr>\n', ...
              '            <td>Date of analysis</td>\n', ...
              '            <td class=\"ta_right\">', datestr(clock()), '</td>\n', ...
              '        </tr>\n', ...
              '       <tr>\n', ...
              '           <td>Algorithm used for segmentation</td>\n', ...
              '           <td class=\"ta_right\">LGA</td>\n', ...
              '       </tr>\n', ...
              '        <tr>\n', ...
              '          <td>T1 image</td>\n', ...
              '          <td class=\"ta_right\">', ps_fileparts(ps_shorten_string(Vt1_tmp.fname, 28), 2:3), '</td>\n', ...
              '       </tr>\n', ...
              '        <tr>\n', ...
              '          <td>FLAIR image</td>\n', ...
              '          <td class=\"ta_right\">', ps_fileparts(ps_shorten_string(Vf2_tmp.fname, 28), 2:3), '</td>\n', ...
              '       </tr>\n', ...
              '        <tr>\n', ...
              '          <td>Initial threshold (kappa)</td>\n', ...
              '          <td class=\"ta_right\">', num2str(kappa_tmp), '</td>\n', ...
              '       </tr>\n', ...
              '        <tr>\n', ...
              '          <td>MRF parameter</td>\n', ...
              '          <td class=\"ta_right\">', num2str(phi), '</td>\n', ...
              '       </tr>\n', ...
              '        <tr>\n', ...
              '          <td>Maximum iterations</td>\n', ...
              '          <td class=\"ta_right\">', num2str(maxiter), '</td>\n', ...
              '       </tr>\n', ...              
              '   </table>\n', ...
              '   </div>\n', ...              
              '   <div class=\"column-02\">\n', ...
              '     <h2>Results</h2>\n', ...
              '         <table style=\"width: 500px\">\n', ...
              '             <tr>\n', ...
              '                 <td>Number of iterations</td>\n', ...
              '                 <td class=\"ta_right\">', num2str(iter), '</td>\n', ...
              '             </tr>\n', ...
              '             <tr>\n', ...
              '                 <td>Lesion map</td>\n', ...
              '                 <td class=\"ta_right\">', Vles.fname, '</td>\n', ...
              '             </tr>\n', ...
              '             <tr>\n', ...
              '                 <td>Lesion volume</td>\n', ...
              '                 <td class=\"ta_right\">', num2str(tlv), ' ml</td>\n', ...
              '             </tr>\n', ...
              '             <tr>\n', ...
              '                 <td>Number of lesions</td>\n', ...
              '                 <td class=\"ta_right\">', num2str(numles), '</td>\n', ...
              '             </tr>\n', ...
              '        </table>\n', ...
              '    </div>\n', ...
              '    <div style=\"clear:both\"></div>\n', ...
              '    <div class=\"column-01\">\n', ...
              '      <h2>Lesion location</h2>\n', ...%'    <img src=\"', fullfile(cd, nameFolder, ['c22_', id, '.png']), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
              '      <img src=\"', ps_fullfile(cd, nameFolder, 'gb_2.png'), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
              '      <img src=\"', ps_fullfile(cd, nameFolder, 'gb_1.png'), '\" width=\"169.41px\" style=\"vertical-align: top;\"><br>\n', ...
              '      <img src=\"', ps_fullfile(cd, nameFolder, 'gb_3.png'), '\" width=\"120px\" style=\"vertical-align: top;\">\n', ...
              '    </div>\n', ...
              '    <div class=\"column-02\" style=\"vertical-align: top;\">\n', ...
              '        <h2>Overlay</h2>\n'];
          if strcmp(pngFailed, '')
               strout = [strout, ...
                  '        <script src=\"', ps_fullfile(cd, nameFolder, 'lga.js'), '\" type=\"text/javascript\"></script>\n', ...              
                  '        <img width=\"450px\" id=\"overlay', jsid, '\" src=\"', ps_fullfile(cd, nameFolder, ['overlay_', num2str(round(mean(r))), '.png']), '\" />\n', ...
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
              '  </div>\n', ...
              '<br><hr>\n']; %% !!
          
          JSid = fopen(fullfile(nameFolder, 'lga.js'), 'wt');
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
              '      $( \"#overlay', jsid, '\" ).attr(\"src\", \"', ps_fullfile(cd, nameFolder), '/overlay_\" + ui.value + \".png\"', ');\n', ...
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
              '            $( \"#overlay', jsid, '\" ).attr(\"src\", \"', ps_fullfile(cd, nameFolder), ['/overlay_\" + slice', jsid, ' + \".png\"'], ');\n', ...
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
              '            $( \"#overlay', jsid, '\").attr(\"src\", \"', ps_fullfile(cd, nameFolder), ['/overlay_\" + slice', jsid, ' + \".png\"'], ');\n', ...
              '            $( \"#slice', jsid, '\").text(\"Slice \" + slice', jsid, ');\n', ...
              '            $(\"#slider_', jsid, '\").slider(\"option\", \"value\", slice', jsid, ');\n', ...
              '        }\n', ...
              '  });', ...
              '});'];
    
            fprintf(JSid, js_strout);
            fclose(JSid);
            
            %HTMLid = fopen(fullfile(nameFolder, [id, '.html']), 'wt');
            fprintf(HTMLid, strout);
                
            HTMLid2 = fopen(fullfile(nameFolder, ['LST_lga_', num2str(kappa_tmp), '_rm', namf2, '.html']), 'wt');
            fprintf(HTMLid2, strout);
            fclose(HTMLid2);

            strout = '  </head>\n  <body>\n';
            fprintf(HTMLid, strout);
            fclose(HTMLid);
            fprintf(fileID, ' ok.\n');

            tt = toc; tt = [num2str(round(tt)), 's'];
            stroutHTML = [repmat(' ', 1, 72 - numel(tt) - numel(stroutHTML)), tt, '\n'];
            fprintf(stroutHTML)
    
            
        end
    end
    
end

% delete log file if segmentation terminated successfully
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