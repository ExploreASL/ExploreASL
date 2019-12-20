%% xASL_qc_TestExploreASL

% Run ExploreASL to get directories
x = ExploreASL_Master('',0);

% Get username so we can locate the paths
if isunix()
	CurrentUser = getenv('USER');
else
	CurrentUser = getenv('username');
end

% Option 1 = run ExploreASL serially
% Option 2 = run ExploreASl parallel (start new MATLAB instances)
% Option 3 = run ExploreASL compilation serially
% Option 4 = run ExploreASL compilation parallel

if isunix && strcmp(CurrentUser,'hjmutsaerts')
    TestDirOrig = '/scratch/hjmutsaerts/TestDataSet/ExploreASL_TestCases';
    TestDirDest = '/scratch/hjmutsaerts/TestDataSet/TempTestResults';
    RunMethod = 2;
elseif ismac && strcmp(CurrentUser,'henk')
    TestDirOrig = '/Users/henk/ExploreASL/TestDataOrig';
    TestDirDest = '/Users/henk/ExploreASL/TestDataDest';
    RunMethod = 1;
elseif ispc && strcmp(CurrentUser,'henk')
    TestDirOrig = 'C:\Users\kyrav\Desktop\SurfDrive\HolidayPics\ExploreASL_TestCases';
    TestDirDest = 'c:\TempTestDir';
    RunMethod = 1;
elseif isunix && strcmp(CurrentUser,'janpetr')
	TestDirOrig = '/pet/projekte/asl/data/ExploreASL_TestCases';
	TestDirDest = '/pet/projekte/asl/data/ExploreASL_TempRes';
	RunMethod = 1;
end

%% Initialize SPM

% Reset path
path(pathdef);
% Remove ExploreASL paths
warning('off','MATLAB:rmpath:DirNotFound');
rmpath(genpath(x.MyPath));
warning('on','MATLAB:rmpath:DirNotFound');

% Add SPM path
addpath(fullfile(x.MyPath,'External','SPMmodified'));

% Initialize SPM, but only SPM
spm('defaults','FMRI');
spm('asciiwelcome');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

%% Copy all data for testing
if exist(TestDirDest,'dir')
    fprintf('Deleting previous results...\n');
    system(['rm -rf ' TestDirDest]);
end
xASL_Copy(TestDirOrig, TestDirDest);

%% Test standalone SPM on low quality
x.Quality = false;
% Find the first directory and copy out the first T1 and ASL just for SPM testing

Dlist = xASL_adm_GetFileList(TestDirDest,'^.*$','List',[0 Inf], true);

xASL_Copy(fullfile(TestDirDest,Dlist{1},'001DM_1','T1.nii'),fullfile(TestDirDest,'T1.nii'));
xASL_Copy(fullfile(TestDirDest,Dlist{1},'001DM_1','ASL_1','ASL4D.nii'),fullfile(TestDirDest,'ASL4D.nii'));

xASL_io_ReadNifti(fullfile(TestDirDest,'ASL4D.nii'));
xASL_io_ReadNifti(fullfile(TestDirDest,'T1.nii'));

% Test spm_realign on low quality
matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.write.data               = {fullfile(TestDirDest,'ASL4D.nii')};
matlabbatch{1}.spm.spatial.realign.write.roptions.which     = [2 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.interp    = 0; % low quality
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap      = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask      = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix    = 'r';
spm_jobman('run',matlabbatch);


% Test spm_reslice on low quality
xASL_spm_reslice(fullfile(TestDirDest,'ASL4D.nii'), fullfile(TestDirDest,'T1.nii'), [], [], 0);

% Test spm_smooth on low quality (small kernel)
xASL_spm_smooth(fullfile(TestDirDest,'rT1.nii'), [2 2 2], fullfile(TestDirDest,'rsT1.nii'));

% Test coreg on low quality
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(TestDirDest,'ASL4D.nii,1')};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {fullfile(TestDirDest,'rT1.nii,1')};
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'mi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = 9;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);

% Test CAT12

matlabbatch = [];
SPMTemplateNII    = fullfile(x.MyPath,'External','SPMmodified', 'tpm', 'TPM.nii');
DartelTemplateNII = fullfile(x.MyPath,'External','SPMmodified', 'toolbox', 'cat12', 'templates_1.50mm', 'Template_1_IXI555_MNI152.nii');
GSTemplateNII     = fullfile(x.MyPath,'External','SPMmodified', 'toolbox', 'cat12', 'templates_1.50mm', 'Template_0_IXI555_MNI152_GS.nii');
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm                         = {SPMTemplateNII};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm   = {DartelTemplateNII}; % Runs DARTEL to this n=555 subjects template
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = {GSTemplateNII}; % Runs Geodesic Shooting to this n=555 subjects template
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 4;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP           = 0; % low quality
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr        = 0; % low quality
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr       = 0; % low quality
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox           = 3; % low quality
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr          = eps; % low quality
matlabbatch{1}.spm.tools.cat.estwrite.opts.samp             = 9;   % low quality
matlabbatch{1}.spm.tools.cat.estwrite.data                  = {fullfile(TestDirDest,'T1.nii')}; % T1.nii
matlabbatch{1}.spm.tools.cat.estwrite.nproc                 = 0; 
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg           = 'mni'; % regularize affine registration for MNI European brains
matlabbatch{1}.spm.tools.cat.estwrite.output.surface        = 0;   % don't do surface modeling
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI  = struct([]); % don't do ROI estimations
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 1;   % save c1T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 0;   % don't save DARTEL space c1T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 1;   % save c2T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 0;   % don't save DARTEL space c2T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native      = 1;   % save c3T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel      = 0;   % don't save DARTEL space c2T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobian.warped= 0;   % don't save Jacobians
matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [1 0]; % save warp to MNI
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 0;   % don't save bias-corrected T1.nii
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native      = 1;   % save c3T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel      = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed= [1 0.1]; % process everything on 1 mm fixed resolution (default)
spm_jobman('run',matlabbatch);

% Test deformations after CAT12
xASL_spm_deformations(x, fullfile(TestDirDest,'rT1.nii'), fullfile(TestDirDest,'rsT1.nii'), 0);

% Remove temporary derivatives
xASL_adm_DeleteFileList(TestDirDest, '.*', false, [0 Inf]);

%% Test ExploreASL itself

% Get list of data to test
Dlist = xASL_adm_GetFileList(TestDirDest,'^.*$','List',[0 Inf], true);
LogFiles = cellfun(@(y) fullfile(TestDirDest,y,'Population','xASL_module_Population.log'), Dlist, 'UniformOutput',false);

for iList=1:length(Dlist)
    AnalysisDir = fullfile(TestDirDest,Dlist{iList});
    DataParFile{iList} = xASL_adm_GetFileList(AnalysisDir,'(DAT|dat|Dat).*\.(m|json)');
    xASL_delete(LogFiles{iList}); % useful for rerun when debugging
    
    % Run ExploreASL
    cd(x.MyPath);
    switch RunMethod
        case 1 % run ExploreASL serially
            ExploreASL_Master(DataParFile{iList}{1}, true, true); % can we run screen from here? or run matlab in background, linux easy
        case 2 % run ExploreASl parallel (start new MATLAB instances)
            MatlabRunString = 'matlab -nodesktop -nosplash -r ';
            RunExploreASLString = ['"cd(''' x.MyPath ''');ExploreASL_Master(''' DataParFile{iList}{1} ''',true,true);exit"'];
            system([MatlabRunString RunExploreASLString ' &']);
        case 3 % run ExploreASL compilation serially
        case 4 % run ExploreASL compilation parallel
        otherwise
    end
end

% Wait until all *.log files exist (which will surely be created, even with a crash)

%% Pause until all results exist
LogsDontExist = max(cellfun(@(y) ~exist(y,'file'), LogFiles));

CountTime = 0;
TimeStepSeconds = 30;
fprintf(xASL_adm_ConvertSeconds2TimeString(CountTime));

while max(cellfun(@(y) ~exist(y,'file'), LogFiles)) % logs don't exist
    pause(TimeStepSeconds);
    CountTime = CountTime+TimeStepSeconds;
    TimeString = xASL_adm_ConvertSeconds2TimeString(CountTime);
    fprintf(['\b\b\b\b\b\b' TimeString]);
end
fprintf('\n');

%% Read results
ResultsTable = {'Data', 'mean_qCBF_TotalGM' 'median_qCBF_TotalGM' 'median_qCBF_DeepWM' 'CoV_qCBF_TotalGM' 'GMvol' 'WMvol' 'CSFvol' 'PipelineCompleted'};
for iList=1:length(Dlist)
    ResultsTable{1+iList,1} = Dlist{iList};
    AnalysisDir = fullfile(TestDirDest,Dlist{iList});
    StatsDir = fullfile(AnalysisDir, 'Population','Stats');
    clear ResultsFile
    ResultFile{1} = xASL_adm_GetFileList(StatsDir,'^mean_qCBF_TotalGM.*PVC2\.tsv','FPList');
    ResultFile{2} = xASL_adm_GetFileList(StatsDir,'^median_qCBF_TotalGM.*PVC0\.tsv','FPList');
    ResultFile{3} = xASL_adm_GetFileList(StatsDir,'^median_qCBF_DeepWM.*PVC0\.tsv','FPList');    
    ResultFile{4} = xASL_adm_GetFileList(StatsDir,'^CoV_qCBF_TotalGM.*PVC0\.tsv','FPList');
    ResultFile{5} = xASL_adm_GetFileList(StatsDir,'^CoV_qCBF_TotalGM.*PVC0\.tsv','FPList');

    for iFile=1:length(ResultFile) % iterate over example datasets
        if length(ResultFile{iFile})<1
            ResultsTable{1+iList,1+iFile} = 'empty';
            ResultsTable{1+iList,2+iFile} = 'empty';
            ResultsTable{1+iList,3+iFile} = 'empty';
        elseif iFile<5 % check the ASL parameters
            [~, TempTable] = xASL_adm_csv2tsv(ResultFile{iFile}{1});
            ResultsTable{1+iList,1+iFile} = TempTable{3,end-2};
        else % check the volumetric parameters
            IndexGM = find(cellfun(@(y) strcmp(y,'GM_vol'), TempTable(1,:)));
            IndexWM = find(cellfun(@(y) strcmp(y,'WM_vol'), TempTable(1,:)));
            IndexCSF = find(cellfun(@(y) strcmp(y,'CSF_vol'), TempTable(1,:)));
            ResultsTable{1+iList,1+iFile} = TempTable{3,IndexGM};
            ResultsTable{1+iList,2+iFile} = TempTable{3,IndexWM};
            ResultsTable{1+iList,3+iFile} = TempTable{3,IndexCSF};
        end
        % check if there are missing lock files
        if exist(fullfile(AnalysisDir,'Missing_Lock_files.csv'),'file')
            ResultsTable{1+iList,4+iFile} = 0; % pipeline not completed
        else
            ResultsTable{1+iList,4+iFile} = 1; % pipeline completed
        end
    end
end
