function varargout = ps_LST_doit(varargin)
%ps_LST_doit   Determination of the initial threshold for the lesion growth
%   algorithm (LGA).
%   Part of the LST toolbox, www.statistical-modelling.de/lst.html
%
%   ps_LST_doit This module offers the opportunity to determine the optimal 
%   initial threshold (kappa) based on the Dice coefficient, see Schmidt 
%   et al. (2012) for details. This requires the existence of reference 
%   segmentations to be compared with the lesion maps. These reference 
%   images are a binary images in the space of the T1 images where a 1 
%   indicates a lesion. 
%
%   This routine saves a CSV file (LST_doit_[date]_[time].csv) in MATLAB's 
%   current directory. The CSV file contains columns for the folder of the 
%   reference images, the name of the FLAIR image, the value of kappa,
%   the Dice coefficient as well as values for sensitivity and specifity.
%
%   ps_LST_doit asks for the reference images by a call to spm_select.
%   The algorithm automatically searches the folder of each reference image
%   for lesion probability maps obtained by LGA, thus the reference images 
%   need to be in the same folder as the lesion probability maps. The
%   threshold for producing binary lesion maps is automatically set to 0.5.
%
%   ps_LST_doit(Vref, thr) compares the reference images given in Vref with
%   lesion probability maps that are obtained by LGA and are in the same
%   folder as the reference images. Lesion probability maps are thresholded
%   by thr. This argument can left empty.
% 


% Welcome text
fprintf(repmat('=', 1, 72));
fprintf('\nRun determination of the initial threshold. Visit\n')
fprintf('www.statistical-modeling.de/lst.html for updates and more information.\n')
fprintf(repmat('=', 1, 72));

% Check input

if ~isempty(varargin) && isfield(varargin{1}, 'data_ref')
    viajob = 1;
else
    viajob = 0;
end
if ~viajob
    if nargin == 0
        Vref = spm_select(Inf, 'image', 'Select reference images.');
        thr = 0.5;
    end
    if nargin == 1
        Vref = varargin{1};
        thr = 0.5;
    end
    if nargin == 2
        Vref = varargin{1};
        thr = varargin{2};
    end    
else
    job = varargin{1};    
    Vref = job.data_ref;    
    thr = job.bin_thresh;
end
Vref = spm_vol(Vref);
res = zeros(0, 5);

pthor = cd;
nameCSV = ['LST_doit_', ps_create_timestamp, '.csv'];
fileID = fopen(nameCSV, 'wt');
fprintf(fileID, 'Folder,FLAIR,kappa,DC,SE,SP\n');

% Summarize input
strout = '\n\nNumber of subjects to process:';
fprintf(strout)
tt = [num2str(numel(Vref)), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 6), tt];
fprintf(strout)
strout = 'Threshold for binary lesion map:';
fprintf(strout)
tt = [num2str(thr), '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
fprintf(strout)

for i = 1:numel(Vref)
    
    strout = '\nWorking on subject';
    fprintf(strout)
    tt = [num2str(i), ' out of ', num2str(numel(Vref)), ' (', num2str(i/numel(Vref)*100), '%%)\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 5), tt];
    fprintf(strout)

    % Load data
    if viajob
        Vref_tmp = Vref{i};
    else    
        Vref_tmp = Vref(i);
    end
    [pthref, namref, extref] = fileparts(Vref_tmp.fname);
    cd(pthref)
    ref = spm_read_vols(Vref_tmp);
    
    % How many lpm's are there?
    files = dir;
    files = arrayfun(@(x) files(x).name, (1:numel(files))', 'UniformOutput', false);
    lpms = files(cellfun(@(x) ~isempty(x), arrayfun(@(x) regexp(x, 'ples_lga_'), files)));    
    tmp = cellfun(@(x) x(2:3), cellfun(@(x) regexp(x, '_'), lpms, 'UniformOutput', false), 'UniformOutput', false);
    kappas = arrayfun(@(i) str2double(lpms{i}((tmp{i}(1)+1):(tmp{i}(2)-1))), (1:numel(lpms))');
    [~,I] = sort(kappas);
    lpms = lpms(I);
    Vlpms = spm_vol(lpms);
    matloaded = 0;
    
    strout = 'Number of lesion maps found';
    fprintf(strout)
    tt = [num2str(numel(lpms)), '\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 2), tt];
    fprintf(strout)
    pthref = strrep(pthref, '\', '/');
        
    for k = 1:numel(lpms)
        
        % Load LST.mat
        us_pos = regexp(Vlpms{k}.fname, '_');        
        namf2 = Vlpms{k}.fname((us_pos(3)+1):end);
        namf2 = namf2(1:(regexp(namf2, '.nii')) - 1);
        kappa_tmp = str2double(Vlpms{k}.fname((us_pos(2)+1):(us_pos(3)-1)));
        if ~matloaded
            load(['LST_lga_', namf2, '.mat'])
            matloaded = 1;
            indx_brain = lga.indx_brain;
        end
        lpm = spm_read_vols(Vlpms{k});
        
        % calculate ...
        % ... FP ...
        FP = sum(ref(indx_brain) < thr & lpm(indx_brain) >= thr);
        % ... TP ...
        TP = sum(ref(indx_brain) >= thr & lpm(indx_brain) >= thr);
        % ... FN ...
        FN = sum(ref(indx_brain) >= thr & lpm(indx_brain) < thr);
        % ... TN ...
        TN = sum(ref(indx_brain) < thr & lpm(indx_brain) < thr);
        % ... sensitivity ...
        SE = TP / (TP + FN);
        % ... specificity ...
        SP = TN / (TN + FP);
        % ... Dice coefficient
        DC = 2 * TP / (2 * TP + FP + FN);
                
        fprintf(fileID, [pthref, ',', namf2, ',', num2str(kappa_tmp), ',', num2str(DC), ...
            ',', num2str(SE), ',', num2str(SP), '\n']);
    end

end

% write results to current folder
cd(pthor)
fclose(fileID);
%csvwrite(fullfile(pwd, nameCSV), res)

strout = '\nResults have been written to the file';
fprintf(strout)
tt = [nameCSV, '\n'];
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 4), tt];
fprintf(strout)

c = clock();
strout = 'Finished successfully ';
fprintf(strout)
tt = datestr(c);
strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout)), tt, '\n\n'];
fprintf(strout)
fprintf(repmat('-', 1, 72));
fprintf('\n')

varargout{:} = 'Don''t forget to cite the toolbox.';%['Finished successfully on ', datestr(c)];

