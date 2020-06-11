function varargout = ps_LST_lpa(varargin)
%ps_LST_lpa   Lesion segmentation by a lesion predicialgorithm (LGA)
%   Part of the LST toolbox, www.statistical-modelling.de/lst.html
%
%   ps_LST_lpa Lesion segmentation by the LPA requires a FLAIR image only.
%   However, the user is free to choose an additional image that serves as
%   a reference image during a coregistration step before the main lesion
%   segmentation. This may be helpful if the dimension of the FLAIR image
%   is low or if the goal of the lesion segmentation is to fill lesions in
%   T1 images. Beside that no additional parameters need to be set. No
%   other parameters need to be set by the user.
%
%   This routine produces lesion probability maps (ples...), [coregistered]
%   bias corrected versions of the FLAIR inputs, a .mat-file with
%   information about the segmentation that is needed for a re-run of the
%   algorithm, and a HTML report along with a folder if this option has
%   been chosen desired. See the documentation for details.
%
%   ps_LST_lpa asks the user for the input images (FLAIR and optional
%   reference images). A HTML report is produced by default.
%
%   ps_LST_lpa(Vf2, Vref, html) performs lesion
%   segmentation for the FLAIR volumes given in Vf2 and reference volumes
%   given in Vref. The letter can be left empty. If specified, both must be
%   characters like a call from from spm_select. The last argument is  a
%   dummy for the HTML reports. If the last argument ismissing it is
%   replaced by its default value, see above.
%
%   ps_LST_lpa(job) Same as above but with job being a harvested job data
%   structure (see matlabbatch help).
%


%spm_jobman('initcfg')
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
tt = 'Lesion prediction algorithm (LPA)\n';
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
if ~isempty(varargin) && isfield(varargin{1}, 'data_F2')
    viajob = 1;
else
    viajob = 0;
end
if ~viajob
    if nargin == 0
        % call to ps_LST_lpa
        fprintf(fileID, 'Select files by function ... ');
        Vf2 = spm_select(Inf, 'image', 'Select FLAIR images.');
        Vref = spm_select(Inf, 'image', 'Select reference images (optional).');
        html_report = 1;
        fprintf(fileID, 'ok.\n');
    end
    if nargin > 0
        html_report = 1;
        if ischar(varargin{1})
            Vf2 = varargin{1};
        else
            fprintf(fileID, 'Vf2 is not a character!');
            fprintf('Input for Vf2 must be a character, like from spm_select.\n');
            fclose(fileID);
            return;
        end
        if nargin > 1
            if ischar(varargin{2})
                Vref = varargin{2};
            else
                fprintf(fileID, 'Vref is not a character!');
                fprintf('Input for Vref must be a character, like from spm_select.\n');
                fclose(fileID);
                return;
            end
        end
        if nargin == 3
            html_report = varargin{3};
        end
        if nargin > 3
            fprintf(fileID, 'Too many input arguments!');
            fprintf('To many input arguments. Did you mean ps_LS_lga?\n');
            fclose(fileID);
            return;
        end
    end
else
    job = varargin{1};
    Vf2 = job.data_F2;
    Vref = job.data_coreg;
    if isempty(Vref{1})
        Vref = [];
    end
    html_report = job.html_report;
	xasl_quality = job.xasl_quality;
end

if ~exist('xasl_quality','var')
	xasl_quality = 1;
end

fprintf(fileID, 'Load volumes ... ');
coreg = ~isempty(Vref);
Vf2 = spm_vol(Vf2);
Vref = spm_vol(Vref);
fprintf(fileID, 'ok.\n');

if coreg && ~isequal(numel(Vf2), numel(Vref))
    fprintf(fileID, 'numel(FLAIR) != numel(ref)');
    error('Number of FLAIR images must match number of reference images.\n');
end
if numel(Vf2) == 0
    fprintf(fileID, 'No images selected.');
    fprintf('No images selected.\n');
    return;
end

% Summarize input
fprintf(fileID, 'Input summary:\n');
fprintf(fileID, ['Jobs: ', num2str(numel(Vf2))]);
strout = 'Number of jobs to process:';
fprintf(strout)
tt = [num2str(numel(Vf2)), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)
fprintf(fileID, ['\nCoregister to reference images: ', num2str(coreg)]);
strout = 'Coregister to reference images:';
fprintf(strout)
if coreg
    tt = 'yes\n';
else
    tt = 'no\n';
end
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)

% Loop over all subjects
for i = 1:numel(Vf2)

    fprintf(fileID, '\n--------------\n');
    fprintf(fileID, ['Job ', num2str(i), ' of ', num2str(numel(Vf2)), '\n']);
    fprintf(fileID, '--------------\n');

    % Extract file information
    if viajob
        Vf2_tmp = Vf2{i};
    else
        Vf2_tmp = Vf2(i);
    end
    Vf2_tmp_or = Vf2_tmp;
    if coreg
        if viajob
            Vref_tmp = Vref{i};
        else
            Vref_tmp = Vref(i);
        end
    end
    [pthf2, namf2, extf2] = fileparts(Vf2_tmp.fname);
    cd(pthf2)

    % Which subject?
    strout = '\nWorking on job';
    fprintf(strout)
    tt = [num2str(i), ' out of ', num2str(numel(Vf2)), ' (', num2str(i/numel(Vf2)*100), '%%)\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 5), tt];
    fprintf(strout)

    % Correct directory?
    cd_tmp = cd;
    cd_tmp = xASL_adm_ConvertSlash( cd_tmp, 1); %% ExploreASL bugfix slashes
    strout = 'Current directory ';
    if numel(cd_tmp) > (72 - numel(strout) - 5)
        fprintf([strout, '...', cd_tmp((end - (72 - numel(strout) - 5)):end), '\n'])
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

    % Create temporary folder
    % ID
    id = [num2str(round(rand(1, 1) * 1e4)), '_', num2str(round(rand(1, 1) * 1e4)), '_', num2str(round(rand(1, 1) * 1e4))];
    tmpFolder = ['LST_tmp_', id];
    xASL_adm_CreateDir(tmpFolder); % EXPLOREASL hack: bugfix

    % get job defaults, ExploreASL hack to speed up
	if xasl_quality==1
		[job1, job2] = ps_LST_lpa_preproc_default;
	elseif xasl_quality==2 % ultralow quality for essentialy skipping preprocessing
        [job1, job2] = ps_LST_lpa_preproc_default_UltraLowQ;
    else
		[job1, job2] = ps_LST_lpa_preproc_default_LowQ;
	end

    % Coregister images
    % -----------------
    % Copy FLAIR into temporary folder
    copyfile(Vf2_tmp.fname, fullfile(tmpFolder, '.'))
    if strcmp(extf2, '.img')
        copyfile([namf2, '.hdr'], fullfile(tmpFolder, '.'))
    end

    coreg = 0; % ExploreASL hack, disable coregistration, for better reproducibility. coreg already performed in ExploreASL

    if coreg% && ~exist(['mr', namf2, '.nii'], 'file')% && ~exist(['r', namf2, extf2], 'file')
        fprintf(fileID, 'Coregistration ...');
        strout = 'Coregister FLAIR to reference ';
        fprintf(strout)
        tic
        job2.ref = {Vref_tmp.fname};
        job2.source = {fullfile(tmpFolder, [namf2, extf2])};
        ps_LST_spm_run_coreg(job2);
        tt = toc; tt = [num2str(round(tt)), 's'];
        strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
        fprintf(strout)
        fprintf(fileID, ' ok.\n');
        %f2del = 1;
        Vf2_tmp = spm_vol(fullfile(tmpFolder, ['r', namf2, extf2]));
        [~, namf2, ~] = fileparts(Vf2_tmp.fname);
        Vf2_tmp_or = Vf2_tmp;
    else
        %f2del = 0;
        fprintf(fileID, 'Skipped coregistration.\n');
        if coreg
            namf2 = ['r', namf2];
        end
    end

    % Calculate bias correction and deformation field
    % -----------------------------------------------

    %if ~(exist(['m', namf2, '.nii'], 'file') && exist(['iy_', namf2, '.nii'], 'file'))
        fprintf(fileID, ['Preprocessing of FLAIR (', namf2, extf2, ') ...']);
        strout = ['Preprocessing of ', namf2, extf2, ' '];
        fprintf(strout)
        tic
        job1.channel.vols = {fullfile(tmpFolder, [namf2, extf2])};
        ps_LST_spm_preproc_run(job1);
        tt = toc; tt = [num2str(round(tt)), 's'];
        strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
        fprintf(strout)
        fprintf(fileID, ' ok.\n');
        copyfile(fullfile(tmpFolder, ['m', namf2, '.nii']), '.')
        %copyfile(fullfile(tmpFolder, [namf2, '.nii']), ['m', namf2, '.nii'])
    %else
    %    fprintf(fileID, 'Skipped preprocessing for FLAIR.\n');
    %    fprintf('Bias corrected FLAIR and inverse deformation field exist\n');
    %end
    %if f2del
    %    spm_unlink(Vf2_tmp.fname);
    %end
    Vf2_tmp = spm_vol(['m', namf2, '.nii']);
    [~, namf2, ~] = fileparts(Vf2_tmp.fname);

    % Load some stuff
    load(fullfile(spm('dir'), 'toolbox', 'LST', 'LST_lpa_stuff.mat'))


    % Warp TPMs and spatial effect to native space
    % --------------------------------------------

    strout = 'Rough segmentation of FLAIR:   \n'; % EXPLOREASL HACK, cosmetic
    fprintf(strout)
    tic

    % Warp brain position of MNI space in native space
    fprintf(fileID, 'Warp brain position ... ');
    cd(tmpFolder)
    bp = ps_LST_lpa_mni2ns(bp_mni, 'bp', Vf2_tmp, [], [], []);
    indx_brain = find(bp > 0);
    fprintf(fileID, 'ok.\n');
    cd ..

    % Exclude some background voxels
    fprintf(fileID, 'Exclude background voxels ... ');
    bg = ps_LST_bc_mni2ns(bg_mni, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);
    indx_brain = find(bp > 0 & bg < 0.3);
    fprintf(fileID, 'ok.\n');
    clear bg;

    % Get TPMs in native space
    fprintf(fileID, 'TPMs in native space ... ');
    csf = ps_LST_bc_mni2ns(csf_mni, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);
    gm = ps_LST_bc_mni2ns(gm_mni, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);
    wm = ps_LST_bc_mni2ns(wm_mni, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);
    fprintf(fileID, 'ok.\n');

    % Rough segmentation of FLAIR into CSF, GM and WM
    % -----------------------------------------------

    fprintf(fileID, 'Rough segmentation of FLAIR ... ');

    % For the calculation of the lesion belief maps we need a rough
    % segmentation of the FLAIR image. As the segmentation by SPM is
    % quite poor we use our one one. Here, we set up a simple mixture
    % model with one Gaussian per tissue class and keep the prior
    % probabilities (the TPMs) for these classes fixed.

    % For the segmentation we smooth the FLAIR image with a Gaussian
    % kernel with FWHM at 1.5 mm
    f2 = spm_read_vols(Vf2_tmp);
    sf2 = 0 .* f2;
    spm_smooth(f2, sf2, [1, 1, 1] .* 1.5);
    f2_vec = sf2(indx_brain);
    clear sf2;

    % Prior probabilities are the TPMs
    prior = [csf(indx_brain), gm(indx_brain), wm(indx_brain)];
    clear gm;
    prior = bsxfun(@times, prior, 1 ./ sum(prior, 2));
    post = 0 .* prior;
    [m, I] = max(prior, [], 2);
    st = 0;
    counter = 0;

    % Main loop of rough segmentation
    while ~st
        counter = counter + 1;

        % Mean and SD of Gaussians
        m = [mean(f2_vec(I == 1)), ...
             mean(f2_vec(I == 2)), ...
             mean(f2_vec(I == 3))];
        s = [std(f2_vec(I == 1)), ...
             std(f2_vec(I == 2)), ...
             std(f2_vec(I == 3))];

        % Posterior probabilities
        for j = 1:3
            %post(:,j) = normpdf(f2_vec, m(j), s(j)) .* prior(:,j);
            post(:,j) = ps_dnorm(f2_vec, m(j), s(j)) .* prior(:,j);
        end
        post = bsxfun(@times, post, 1 ./ sum(post, 2));

        % New class labels
        [~, I_new] = max(post, [], 2);
        % Correct missclassified voxels, i.e. voxels that are segmented
        % as CSF but have intensities larger than the mean of GM
        I_new(f2_vec > m(2) & I_new == 1) = 2;

        % Check if the label of any voxels have been changed.
        if isequal(I, I_new)
            st = 1;
        else
            I = I_new;
        end
    end % END while ~st

    % One last iteration. Here we use the last posterior probability as
    % the prior probabilities. This helps in identifying CSF voxels
    % better.
    for j = 1:3
        %post(:,j) = normpdf(f2_vec, m(j), s(j)) .* mean(I_new == j);
        post(:,j) = ps_dnorm(f2_vec, m(j), s(j)) .* mean(I_new == j);
    end
    post = bsxfun(@times, post, 1 ./ sum(post, 2));
    [~, I_new] = max(post, [], 2);
    I_new(f2_vec > m(2) & I_new == 1) = 2;

    seg = 0.*f2;
    seg(indx_brain) = I_new;
    % Correct voxels in the TPM for CSF
    csf(seg > 1) = 0;

    fprintf(fileID, ['ok, with ' num2str(counter), ' iterations.\n']);

    tt = toc; tt = [num2str(round(tt)), 's'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
    fprintf(strout)

    %% Predict probabilities
    %% ---------------------


    strout = 'Predict lesion probabilities ';
    fprintf(strout)
    tic

    fprintf(fileID, 'Lesion belief maps ... ');
    % For the lesion belief map we use the original FLAIR image
    f2_vec = f2(indx_brain);
    tmp = f2_vec(I == 2);
    f2_vec = f2_vec ./ mean(tmp(~isnan(tmp)));

    % Lesion belief map for WM
    Bf2_wm = 0 .* f2;
    tmp = f2_vec(I == 3);
    Bf2_wm(indx_brain) = (f2_vec - mean(tmp(~isnan(tmp))));
    Bf2_wm(Bf2_wm < 0) = 0;

    % Lesion belief map for GM
    Bf2_gm = 0 .* f2;
    tmp = f2_vec(I == 2);
    Bf2_gm(indx_brain) = (f2_vec - mean(tmp(~isnan(tmp))));
    Bf2_gm(Bf2_gm < 0) = 0;

    % The final lesion belief map is weightes by the prior probability
    % of WM.
    gm = ps_LST_bc_mni2ns(gm2_mni, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);
    wmprob = 1 - (csf + gm);
    wmprob(wm < .05) = 0;
    clear csf; clear gm; clear wm;

    % Final lesion belief map
    Bf2 = (Bf2_gm + Bf2_wm) .* wmprob;
    Bf2_2 = (Bf2_gm + Bf2_wm);
    clear wmprob; clear Bf2_gm; clear Bf2_wm;
    fprintf(fileID, 'ok.\n');
    %Vf2.fname = 'Bf2.nii';
    %spm_write_vol(Vf2, Bf2);

    %{
    gm_mask = spm_read_vols(spm_vol('/Users/paul/201503_NewSegmentation/Data/gm-mask4.img'));
    wm_mask = spm_read_vols(spm_vol('/Users/paul/201503_NewSegmentation/Data/wm-mask.img'));
    for z = 10:111
        gm_mask(z,:,:) = gm_mask(121-z,:,:);
    end
    gm2_img = zeros(121, 145, 121); gm2_img(bp_mni) = gm_mni;
    gm2_img(gm_mask > 0) = gm2_img(gm_mask > 0) ./ 2;
    gm2_mni = gm2_img(bp_mni);

    sp = zeros(121,145,121);
    sp(bp_mni) = sp_mni_Bf2;
     for z = 1:60
            tmp = (sp(z,:,:) + sp(121-z,:,:))./2;
            sp(z,:,:) = tmp;
            sp(121-z,:,:) = tmp;
        end

    bp_img = 0 .* sp; bp_img(bp_mni) = 1;
    sp(abs(sp) < 0.1 & bp_img > 0) = 100;
    sp(sp == 0) = NaN;

    st = 0;
    while ~st
        indx_tmp = find(sp > 99);
        if isempty(indx_tmp)
           st = 1;
        else
            nh = getNeighborhood2(sp, indx_tmp, 1);
            f = find(sum(nh < 99) > 0);
            mnh = max(nh .* (nh < 99 & abs(nh) > .05));
            sp(indx_tmp(f)) = mnh(f);
        end
    end
    sp(isnan(sp)) = 0;

    sp_mni2 = sp(bp_mni);

    %}


    % Calculate lesion probability
    fprintf(fileID, 'Lesion probability maps ... ');

    % Warp spatial effect to native space
    sp_Bf2 = ps_LST_bc_mni2ns(sp_mni2_Bf2, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);

    % Linear predictor
    eta = intercept + Bf2_eff .* Bf2(indx_brain) + sp_Bf2(indx_brain);
    prob = 0 .* f2;
    % Probabilities
    prob(indx_brain) = 1 ./ (1+ exp(-eta));

    % In order to obtain a smoother version we smooth all voxels > 0.1
    % with a Gaussian kernel with FWHM at 1 mm
    prob(prob < 0.1) = NaN;
    sprob = 0 .* f2;
    spm_smooth(prob, sprob, [1, 1, 1]);
    prob(prob > 0) = sprob(prob > 0);
    prob(isnan(prob)) = 0;
    nl = ps_LST_bc_mni2ns(noles, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);
    prob(nl > 0) = 0;

    fprintf(fileID, 'ok.\n');
    tt = toc; tt = [num2str(round(tt)), 's'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
    fprintf(strout)

    %% ExploreASL hack
    % If WMH_SEGM.nii already exists, replace prob by this
    % Allows to follow the same Clean Up procedure as below
    % And also add some extra cleaning here

%     [Fpath, ~, ~] = xASL_fileparts(Vf2{1}.fname);
%     WMH_SEGMpath = fullfile(Fpath,'WMH_SEGM.nii');
%     if xASL_exist(WMH_SEGMpath, 'file')
%         prob = xASL_io_Nifti2Im(WMH_SEGMpath);
%         % Convert NaNs to zeros & positive & cutoff extremely low probabilities
%         prob(~isfinite(prob)) = 0;
%         prob(prob<0.001) = 0;
%         prob(prob>1) = 1;
%         fprintf('\n-> -> LST segmentation replaced by pre-existing WMH_SEGM <- <-\n\n');

%         PriorPath = fullfile(Fpath,'atlas_wm.nii');
%         bFunction = exist('xASL_spm_deformations');
%         if isdeployed || (bFunction == 2) || (bFunction == 5)
%         
%         
%             % Multiply with MNI WM prior probability
%             PriorWM = fullfile(spm('dir'),'toolbox','LST','atlas_wm.nii');
%             if ~exist(PriorPath,'file')
%                 ExistPrior = 0;
%                 xASL_spm_deformations([],PriorWM,PriorPath,1,[],[],fullfile(Fpath,tmpFolder,'iy_rFLAIR.nii'));
%             else
%                 ExistPrior  = 1;
%             end
% 
%             PriorIM = xASL_io_Nifti2Im(PriorPath);
%             PriorIM(PriorIM>0.5) = 1;
%             PriorIM(PriorIM<0.5) = PriorIM(PriorIM<0.5).^0.33;
% 
%             prob = PriorIM.*prob;
% 
%             if ~ExistPrior
%                 xASL_delete(PriorPath);
%             end
%         end
%     end

    %% Clean up
    %% --------

    strout = 'Clean up ';
    fprintf(strout)
    tic

    % Delete all voxels that have no direct neighbor
    indx_les = find(prob > 0);
    nh = getNeighborhood2(1 .* (prob > 0), indx_les, 1);
    prob(indx_les(sum(nh) == 0)) = 0;

    % Delete voxels that are too close to outer CSF voxels

    fprintf(fileID, 'Select outer voxels ... ');
    % Outer voxels
    outerCSF = 0 .* f2;
    nh = getNeighborhood2(seg, indx_brain, 3);
    outerCSF(indx_brain(sum(nh == 0) > 0)) = 1;
    % Add all CSF voxels except for the ventricles
    vm = ps_LST_bc_mni2ns(vm_mni, 'tissue', Vf2_tmp, bp_mni, indx_brain, bp);
    outerCSF(seg == 1 & vm == 0) = 1;
    fprintf(fileID, 'ok.\n');

    % In order to calculate the correct distance (in mm) we need to
    % correct for different voxel dimensions. We do this by artifically
    % blowing the brain up until all voxels have a natural number as
    % voxel dimension.

    fprintf(fileID, 'Voxel size ...');
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
    fprintf(fileID, 'ok.\n');

    %%% ExploreASL hack - we use our own implementation of bwdist
    %%% if ~exist('bwdist', 'builtin')
    %%%     fprintf(fileID, 'Calculate distance ... ');
    %%%     indx_les = find(prob > 0);
    %%%     coords_csf = indx2coord(find(outerCSF > 0), size(f2, 1), size(f2, 2));
    %%%     coords_les = indx2coord(indx_les, size(f2, 1), size(f2, 2));

    %%%     %tic
    %%%     d = zeros(numel(indx_les), 1);
    %%%     for k = 1:size(coords_les, 1)
    %%%         d_tmp = ((coords_csf(:,1) - coords_les(k,1)) .* vs(1)).^2 + ...
    %%%                   ((coords_csf(:,2) - coords_les(k,2)) .* vs(2)).^2 + ...
    %%%                   ((coords_csf(:,3) - coords_les(k,3)) .* vs(3)).^2;
    %%%         d(k) = min(d_tmp);
    %%%     end
    %%%     outerCSF_dist = 0 .* f2;
    %%%     outerCSF_dist(indx_les) = d;
    %%%     %toc
    %%%     fprintf(fileID, 'ok.\n');
    %%%
    %%% else

        fprintf(fileID, 'Rounded voxel size ...');
        % Round voxel size to nearest .5-value
        for j = 1:3
            ch = .5:.5:5;
            [~, mi] = min(abs(vs(j) - ch));
            vs(j) = ch(mi);
        end
        fprintf(fileID, [' ok, ', num2str(vs), '\n']);

        fprintf(fileID, 'Factor ...');
        % Find factor that gives a natural number for each voxel size
        fac_rs = 1;
        vs_rs = vs;
        while any(mod(vs_rs, 1) > 0)
            fac_rs = fac_rs + 1;
            vs_rs = vs .* fac_rs;
        end
        fprintf(fileID, [' ok, ', num2str(fac_rs), '\n']);

        fprintf(fileID, 'Translate coordinates ... ');
        % Create an empty image with the new voxel size
        outerCSF_rs = zeros(size(f2, 1) .* vs_rs(1), ...
            size(f2, 2) .* vs_rs(2), ...
            size(f2, 3) .* vs_rs(3));

        % Translate coordinate of the original image into coordinates of
        % the new image
        coords = indx2coord(indx_brain, size(f2, 1), size(f2, 2));
        coords_rs = coords;
        for j = 1:3
            if vs_rs(j) > 0
                coords_rs(:,j) = round(ps_scale(coords(:,j), ...
                    min(coords(:,j)), ...
                    max(coords(:,j)) .* vs_rs(j)));
            end
        end
        % Translate coordinates back to indx
        indx_rs = coord2indx(coords_rs(:,1), coords_rs(:,2), coords_rs(:,3), ...
                    size(f2, 1) .* vs_rs(1), size(f2, 2) .* vs_rs(2));
        fprintf(fileID, 'ok.\n');

        fprintf(fileID, 'Calculate distance ... ');
        % Fill new image
        outerCSF_rs(indx_rs) = outerCSF(indx_brain);
        % Calculate distance of all voxels to non-zero voxels

        %%% ExploreASL hack - use our own function to calculate Euclidean distance
        %%% outerCSF_rs_dist = bwdist(outerCSF_rs);
        [outerCSF_rs_dist,~,~,~] = xASL_im_DistanceTransform( outerCSF_rs );

        clear outerCSF_rs;
        % Put calculated distances back to original image
        outerCSF_dist = 0 .* f2;
        outerCSF_dist(indx_brain) = outerCSF_rs_dist(indx_rs) ./ fac_rs;
        clear outerCSF_rs_dist;
        fprintf(fileID, 'ok.\n');
    %%% end

    % Finally: delete all voxels whose distance to outer CSF is smaller
    % or equal than 4 mm
    fprintf(fileID, 'Delete voxels ... ');
    prob2 = prob;
    prob2(outerCSF_dist <= 4) = 0;
    fprintf(fileID, 'ok.\n');

    % Correct for big lesions that are too close at outer CSF
    fprintf(fileID, 'Correct for big lesions that are too close at outer CSF ... ');
    indx_les = find(prob > 0);
    bw = ps_bwlabeln(prob > 0);
    
    if ~isempty(indx_les) && sum(bw(:))>0 % Skip this part when no lesions found close to outer CSF %% EXPLOREASL HACK

        bw_tmp = bw(indx_les);
        p_tmp = prob(indx_les);
        p2_tmp = prob2(indx_les);
        gs1 = arrayfun(@(x) mean(p_tmp(bw_tmp == x)), (1:max(bw_tmp))');
        gs2 = arrayfun(@(x) mean(p2_tmp(bw_tmp == x)), (1:max(bw_tmp))');
        %gs1 = grpstats(p_tmp, bw_tmp, 'mean');
        %gs2 = grpstats(p2_tmp, bw_tmp, 'mean');
        indx_correct = find(gs1 > .5 & gs2 < .5 & gs2 > 0);
        if ~isempty(indx_correct)
            for j = indx_correct'
                prob2(bw == j) = prob(bw == j);
            end
        end
    end
    prob = prob2; clear prob2;
    fprintf(fileID, 'ok.\n');

    % Delete small lesions with low probability
    fprintf(fileID, 'Delete small lesions with low probability ... ');
    bw = ps_bwlabeln(prob > 0);
    volfactor = abs(det(Vf2_tmp.mat(1:3,1:3))) /  1000;
    
    if sum(bw(:))>0 % Skip this part when no small lesions found with low probability %% EXPLOREASL HACK
        for j = 1:max(bw(bw > 0))
            indx_tmp = find(bw == j);
            if (sum(prob(indx_tmp) > .1)*volfactor < 3*0.003)
                prob(indx_tmp) = 0;
            end
        end
    end
    fprintf(fileID, 'ok.\n');

    % Delete holes
    fprintf(fileID, 'Delete holes ... ');
    ch = 0;
    iter = 0;
    while ~ch
        iter = iter + 1;
        indx_tmp = find(prob == 0 & bp > 0);
        nh = getNeighborhood2(prob, indx_tmp, 1);
        indx_tmp = indx_tmp(sum(nh > 0) > 4);
        if isempty(indx_tmp) || iter > 10
            ch = 1;
        else
            prob(indx_tmp) = mean(nh(:,sum(nh > 0) > 4));
        end
    end
    fprintf(fileID, ['ok, with ', num2str(iter), ' iterations.\n']);

    % A bit of growing
    fprintf(fileID, 'Grow a bit ... ');
    st = 0;
    iter = 0;
    while ~st
        iter = iter + 1;
        indx_tmp = find(Bf2_2 > 0.5 & prob == 0);
        nh = getNeighborhood2(prob, indx_tmp, 3);
        tmp_tmp = find(mean(nh > 0) > 0.25);
        if numel(tmp_tmp) > 10
            tmp = nh(:,tmp_tmp); tmp(tmp == 0) = NaN;
            prob(indx_tmp(tmp_tmp)) = ps_quantile(tmp, .85);
        else
            st = 1;
        end
    end
    fprintf(fileID, ['ok, with ', num2str(iter), ' iterations.\n']);

    tt = toc; tt = [num2str(round(tt)), 's'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n'];
    fprintf(strout)

    % Write results
    % -----------------------------------------------------------------

    % Save results in LST.mat
    fprintf(fileID, 'Save LST.mat ... ');
    lpa.indx_brain = indx_brain;
    lpa.bp_indx = bp(indx_brain);
    lpa.seg_vec = I_new;
    lpa.I = I;
	% lpa.Bf2 = Bf2;  %% ExploreASL hack: comments added, which were here in removed/unused function
    % lpa.Bf2_2 = Bf2_2; %% ps_LST_lpa_2.m. otherwise, these functions were identical
    lpa.f2_vec = f2(:);
    if coreg
        lpa.Vref = Vref_tmp;
    end
    %Vf2_tmp.dt = [16 0];
    lpa.Vf2 = Vf2_tmp;
    % Should we flip the images?
    V = spm_vol(fullfile(spm('dir'), 'toolbox', 'LST', 'atlas_wm.nii'));
    or_ch = zeros(V.dim);
    for j = 1:121
        or_ch(60,70,j) = 1;
        or_ch(61,71,j) = 1;
        or_ch(59,69,j) = 1;
        or_ch(j,70,60) = 2;
        or_ch(j,71,61) = 2;
        or_ch(j,69,59) = 2;
    end
    for j = 1:72
        or_ch(60,j,60) = 3;
        or_ch(61,j,61) = 3;
        or_ch(59,j,59) = 3;
        or_ch(60,146-j,60) = 4;
        or_ch(61,146-j,61) = 4;
        or_ch(59,146-j,59) = 4;
    end
    cd(tmpFolder)
    V.fname = 'or_ch.nii';
    spm_write_vol(V, or_ch);
    clear job
    job.comp{1}.def = {['iy_', namf2(2:end), '.nii']};
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
    if any(or == 0)
        unique(or(or > 0))
    end
    % Flip?

    %%% ExploreASL hack: flipping bugfix
    if  length(unique(or))<3
        % fix this, by assuming that the repeated value is incorrect
        if  sum(or(1:2) - or(2:3))>0
            or = [3 2 1]; % flip
        else
            or = [1 2 3]; % don't flip
        end
    end

    or_ch = permute(or_ch, or);
    c_tmp3 = indx2coord(find(or_ch == 3), size(or_ch, 1), size(or_ch, 2));
    c_tmp4 = indx2coord(find(or_ch == 4), size(or_ch, 1), size(or_ch, 2));
    fl =  max(c_tmp3(:,2)) > max(c_tmp4(:,2));
    lpa.or = or;
    lpa.fl = fl;
    save(['LST_lpa_', namf2, '.mat'], 'lpa')
    fprintf(fileID, 'ok.\n');

    fprintf(fileID, 'Write results ... ');
    Vles = Vf2_tmp_or;
    Vles.dt = Vf2_tmp.dt;
    Vles.fname = ['ples_lpa_', namf2, '.nii'];
    Vles.descrip = 'Probability lesion map obtained by LPA within LST toolbox';
    spm_write_vol(Vles, prob);
    rmdir(tmpFolder, 's')
    fprintf(fileID, 'ok.\n');

    if html_report

        % HTML report
        % -----------------------------------------------------------------

        stroutHTML = 'Create HTML report';
        fprintf(stroutHTML);
        tic

        % create HTML report
        nameFolder = ['LST_lpa_' namf2];
        warning('off');
        mkdir(nameFolder)
        warning('on');

        % Create PNGs
        fprintf(fileID, 'Create PNGs ...');
        pngFailed = '';
        try
            [~, r] = ps_LST_create_gif(fullfile(cd, [namf2, '.nii']), ...
                Vles.fname, nameFolder, 0, [or, fl]);
            fprintf(fileID, ' ok.\n');
        catch ME
            fprintf(fileID, ' failed!.\n');
            pngFailed = ME.message;
            r = 0:1;
        end

        % Create images for glass brains
        ps_LST_create_glass_brain(1 .* (prob > .5), ...
            1 .* (seg > 0), ...
            fullfile(nameFolder, 'gb'), [or, fl]);

        fprintf(fileID, 'Create main HTML file ...');
        % Main HTML file#
        nameHTML = ['report_LST_lpa_', namf2, '.html'];
        copyfile(fullfile(spm('dir'), 'toolbox', 'LST', 'LST_main_html.html'), nameHTML)
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
        if any(prob(:) > .5)
            bw = ps_bwlabeln(1. * (prob > .5));
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
          '           <td class=\"ta_right\">LPA</td>\n', ...
          '       </tr>\n', ...
          '        <tr>\n', ...
          '          <td>FLAIR image</td>\n', ...
          '          <td class=\"ta_right\">', ps_fileparts(ps_shorten_string(Vf2_tmp.fname, 28), 2:3), '</td>\n', ...
          '       </tr>\n', ...
          '        <tr>\n', ...
          '          <td>Reference image</td>\n'];%, ...
      if coreg
           strout = [strout, ...
               '          <td class=\"ta_right\">', ps_fileparts(ps_shorten_string(Vref_tmp.fname, 28), 2:3), '</td>\n'];
      else
           strout = [strout, ...
               '          <td class=\"ta_right\">none</td>\n'];
      end
      strout = [strout, ...
          '       </tr>\n', ...
          '   </table>\n', ...
          '   </div>\n', ...
          '   <div class=\"column-02\">\n', ...
          '     <h2>Results</h2>\n', ...
          '         <table style=\"width: 500px\">\n', ...
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
                  '        <script src=\"', ps_fullfile(cd, nameFolder, 'lpa.js'), '\" type=\"text/javascript\"></script>\n', ...
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

              %'      $( \"#overlay', jsid, '\" ).attr(\"src\", \"', ps_fullfile(cd, nameFolder, 'overlay_\" + ui.value + \".png\"'), ');\n', ...
              %'            $( \"#overlay', jsid, '\" ).attr(\"src\", \"', ps_fullfile(cd, nameFolder, ['overlay_\" + slice', jsid, ' + \".png\"']), ');\n', ...
              %'            $( \"#overlay', jsid, '\").attr(\"src\", \"', ps_fullfile(cd, nameFolder, ['overlay_\" + slice', jsid, ' + \".png\"']), ');\n', ...
          JSid = fopen(fullfile(nameFolder, 'lpa.js'), 'wt');
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

        HTMLid2 = fopen(fullfile(nameFolder, ['LST_lpa_', namf2, '.html']), 'wt');
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
