function [ResultsTable] = xASL_qc_TestExploreASL(TestDirOrig, TestDirDest, RunMethod, bTestSPM, MatlabPath, EmailAddress, Password, bOverwrite, testDataUsed, RunTimePath, bPull)
%xASL_qc_TestExploreASL Do a thorough test of the validity and reproducibility of ExploreASL
%
% FORMAT: [ResultsTable] = xASL_qc_TestExploreASL(TestDirOrig, RunMethod)
% 
% INPUT:
%   TestDirOrig - path to root folder containing all datasets to test processing on (REQUIRED)
%   TestDirDest - path to folder where the results are stored (this is temporarily, automatically deleted upon restart, REQUIRED)
%   RunMethod   - Choose a value for how to test ExploreASL (REQUIRED)
%                 Option 1 = run ExploreASL serially
%                 Option 2 = run ExploreASl parallel (start new MATLAB instances)
%                 FUTURE Option 3 = run ExploreASL compilation serially
%                 FUTURE Option 4 = run ExploreASL compilation parallel
%   bTestSPM    - boolean for testing if SPM standalone with xASL modifications works (DEFAULT=true)
%   MatlabPath  - path to matlab executable or compilation bash script (OPTIONAL, required in some 
%                 cases)
%   EmailAddress- string with e-mail address for gmail account to use (OPTIONAL, DEFAULT = skip e-mailing results)
%   Password    - string with password for this gmail account (REQUIRED when EmailAddress provided)
%   bOverwrite  - Overwrite existing test results (OPTIONAL, DEFAULT=true);
%   testDataUsed- Option 1: TestDataSet as an input
%                 Option 0: Other (DEFAULT)
%   RunTimePath - When using a compiled version, the location of the
%                 Matlab RunTime libraries (e.g. '/usr/local/MATLAB/MATLAB_Runtime/v96')
%   bPull       - pull new version of the software (OPTIONAL, DEFAULT=true)
%                 
% OUTPUT:
%   ResultsTable - Table containing all results from the test runs
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function will run ExploreASL on several different
%              datasets, do perform a full and thorough test of ExploreASL.
%              Different environment variables include:
%              x.Quality 0 and 1
%              
%              Different data setups include:
%              ASL readouts (3D spiral, 3D GRASE, 2D EPI)
%              ASL vendors (GE, Philips, Siemens)
%              With/without background suppression
%              With/without FLAIR processing (LST LGA or LPA)
%              With/without lesion masking from tumor
%              Performed as first or second run of ExploreASL, fully or partly done before
%              Longitudinal or cross-sectional data
%              FEAST
%              With/without disabling quantification
%              With/without saving 4D PWI/CBF maps
%              With/without M0
%              With single or multiple ASL sessions (i.e. runs)
%
%              This function performs the following steps:
%              1) Pull latest GitHub version
%              2) Initialize SPM
%              3) Copy all data for testing
%              4) Test standalone SPM on low quality
%              5) Test ExploreASL itself
%              6) Pause until all results exist (if running parallel in background)
%              7) Compile results table
%              8) Compare table with reference table
%              9) E-mail results
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE for Jan: [ResultsTable] = xASL_qc_TestExploreASL('/pet/projekte/asl/data/ExploreASL_TestCases', '/pet/projekte/asl/data/ExploreASL_TempRes', 1);
% EXAMPLE for Henk on MacOS: [ResultsTable] = xASL_qc_TestExploreASL('/Users/henk/surfdrive/HolidayPics/ExploreASL_TestCases', '/Users/henk/ExploreASL/ASL/ExploreASL_TestCasesProcessed', 1, 0,[],'henkjanmutsaerts@gmail.com');
% EXAMPLE for VUmc server: [ResultsTable] = xASL_qc_TestExploreASL('/radshare/ExploreASL_Test/ExploreASL_TestCases', '/radshare/ExploreASL_Test/ExploreASL_TestCasesProcessed', 1);
% __________________________________
% Copyright 2015-2020 ExploreASL
% ToDo: split this function in subfunctions, each function/subfunction
% should be short, making it easier to read what this function does

% ============================================================
%% Admin
% Run ExploreASL to get directories
if isempty(which('ExploreASL_Master'))
    cd ..;
else
    cd(fileparts(which('ExploreASL_Master')));
end
if nargin<3
    RunMethod = 1; % Set 'serial' to be default
end
if nargin<4 || isempty(bTestSPM)
    bTestSPM = true;
end
if nargin<5 || isempty(MatlabPath)
    MatlabPath = 'matlab';
end
if nargin<6 || isempty(EmailAddress)
    EmailAddress = [];
    Password = [];
end

if ~isempty(EmailAddress) && (nargin<7 || isempty(Password))
    warning('Please provide password for g-mail account!');
end

if nargin<8 || isempty(bOverwrite)
    bOverwrite = true;
end

if nargin<9 || isempty(testDataUsed)
    testDataUsed = 0;
end

if RunMethod>2
    % We will test the compiled version, but do some checks first
    if isempty(MatlabPath) || ~exist(MatlabPath, 'file') || ~strcmp(MatlabPath(end-2:end),'.sh')
        warning('Please provide the path to the bash script calling the compiled ExploreASL, skipping');
        return;
    elseif nargin<10 || isempty(RunTimePath) || ~exist(RunTimePath, 'dir')
        warning('Please provide the path to the Matlab Runtime installation, skipping');
        return;        
    end
end

if nargin<11 || isempty(bPull)
	bPull = 1;
end

% ============================================================
%% 1) Pull latest GitHub version
% assuming we are in ExploreASL folder
if bPull
	Answer = system('source ~.bashrc;git fetch','-echo');
	Answer = system('source ~.bashrc;git pull','-echo');
end

% Initialize ExploreASL
x = ExploreASL_Initialize;

clc;

% ============================================================
%% 2) Initialize SPM
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

% ============================================================
%% 3) Copy all data for testing

% Ask for directories if they were not defined
if ~exist('TestDirDest','var'), TestDirDest = uigetdir(pwd, 'Select testing directory...'); end
if ~exist('TestDirOrig','var'), TestDirOrig = uigetdir(pwd, 'Select datasets for testing...'); end

% Clone testdataset repository if not detected
if ~exist(TestDirOrig, 'dir')
    TestDirRoot = fileparts(TestDirOrig);
    TestDataSetRepository = 'https://github.com/ExploreASL/TestDataSets.git';
    fprintf('%s\n', ['TestDataSet repository not found in: ' TestDirRoot]);
    fprintf('%s\n', ['Attempting to clone: ' TestDataSetRepository]);
    system(['source ~.bashrc;cd ' TestDirRoot '; git clone ' TestDataSetRepository], '-echo');
end

% Remove previous results
if bOverwrite && exist(TestDirDest,'dir')
    fprintf('Deleting previous results...\n');
    if ispc
        system(['rmdir /s /q ' TestDirDest]);
    else
        system(['rm -rf ' xASL_adm_UnixPath(TestDirDest)]);
    end
end
% Copy data sets into testing directory
xASL_Copy(TestDirOrig, TestDirDest);    



% ============================================================
%% 4) Test standalone SPM on low quality
if bTestSPM

    x.Quality = false;
    % Find the first directory and copy out the first T1 and ASL just for SPM testing

    Dlist = xASL_adm_GetFileList(TestDirDest,'^.*$','List',[0 Inf], true);
    
    if testDataUsed
        % TestDataSet detected
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'T1.nii'),fullfile(TestDirDest,'T1.nii'));
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'ASL_1','ASL4D.nii'),fullfile(TestDirDest,'ASL4D.nii'));
    else
        % Default
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'001DM_1','T1.nii'),fullfile(TestDirDest,'T1.nii'));
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'001DM_1','ASL_1','ASL4D.nii'),fullfile(TestDirDest,'ASL4D.nii'));
    end
    
    % Read Nifti files
    try
        xASL_io_ReadNifti(fullfile(TestDirDest,'ASL4D.nii'));
        xASL_io_ReadNifti(fullfile(TestDirDest,'T1.nii'));
    catch
        error('ASL and T1 data set import failed...')
    end
    
    

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
	[~,catVer] = cat_version();
	if str2double(catVer) > 1500
		catTempDir = 'templates_volumes';
	else
		catTempDir = 'templates_1.50mm';
	end
    matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm              = {SPMTemplateNII};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP           = 0; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr        = 0; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr       = 0; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox           = 6; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr          = eps; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.opts.samp             = 9;   % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.data                  = {fullfile(TestDirDest,'T1.nii')}; % T1.nii
    matlabbatch{1}.spm.tools.cat.estwrite.nproc                 = 0; 
    matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg           = 'mni'; % regularize affine registration for MNI European brains
    matlabbatch{1}.spm.tools.cat.estwrite.output.surface        = 0;   % don't do surface modeling
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI  = struct([]); % don't do ROI estimations
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 0;   % save c1T1 in native space
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 0;   % don't save modulation
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 0;   % don't save DARTEL space c1T1, this happens below in the reslice part
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 0;   % save c2T1 in native space
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 0;   % don't save modulation
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 0;   % don't save DARTEL space c2T1, this happens below in the reslice part
    matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [1 0]; % save warp to MNI
    matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 0;   % don't save bias-corrected T1.nii
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed= [1 0.1]; % process everything on 1 mm fixed resolution (default)
    spm_jobman('run',matlabbatch);

    % Test deformations after CAT12
    xASL_Copy(fullfile(TestDirDest,'mri','y_T1.nii'), fullfile(TestDirDest,'y_T1.nii'));
    matlabbatch = [];
    matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(TestDirDest,'y_T1.nii')};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(TestDirDest,'rT1.nii')};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
    spm_jobman('run',matlabbatch);

    % Remove temporary derivatives
    xASL_adm_DeleteFileList(TestDirDest, '.*', false, [0 Inf]);
    DirsAre = {'err' 'ASL_1' 'mri' 'report'};
    for iDir=1:length(DirsAre)
        DirIs = fullfile(TestDirDest,DirsAre{iDir});
        xASL_adm_DeleteFileList(DirIs, '.*', true, [0 Inf]);
        xASL_delete(DirIs);
    end
end

% ============================================================
%% 5) Test ExploreASL itself
x = ExploreASL_Initialize; % here we return the ExploreASL paths, which we removed above for testing SPM

% Remove lock folders, useful for rerun when debugging
LockFolders = xASL_adm_GetFileList(TestDirDest, '^locked$', 'FPListRec', [0 Inf], true);
if ~isempty(LockFolders)
    cellfun(@(y) xASL_delete(y), LockFolders, 'UniformOutput', false);
end

% Evaluate parallelization in linux
if isunix && (RunMethod==2 || RunMethod==4)
    [Result1, Result2] = system('screen -dmS TryOut exit');
    if Result1~=0
        warning('Please install screen for testing ExploreASL parallel in a mac/linux');
        fprintf('%s\n', Result2);
        error('Skipping');
    end
end

% Get list of data to test
Dlist = xASL_adm_GetFileList(TestDirDest,'^.*$','List',[0 Inf], true);
LogFiles = cellfun(@(y) fullfile(TestDirDest,y,'Population','xASL_module_Population.log'), Dlist, 'UniformOutput',false);

for iList=1:length(Dlist)
    AnalysisDir = fullfile(TestDirDest,Dlist{iList});
    DataParFile{iList} = xASL_adm_GetFileList(AnalysisDir,'(DAT|dat|Dat).*\.(m|json)');
    xASL_delete(LogFiles{iList}); % useful for rerun when debugging
    ScreenName = ['TestxASL_' num2str(iList)];
    
    if ~isempty(DataParFile{iList})
        try
            % Run ExploreASL
            cd(x.MyPath);
            
            if RunMethod>2 % prepare compilation testing
                [Fpath, Ffile, Fext] = fileparts(MatlabPath);
                if isunix
                    CompilationString = ['cd ' Fpath ';bash ' Ffile Fext ' ' RunTimePath ' ' DataParFile{iList}{1}];
                else
                    CompilationString = ['cd ' Fpath '; ' Ffile Fext ' ' RunTimePath ' ' DataParFile{iList}{1}];
                end
            end
            
            switch RunMethod
                case 1 % run ExploreASL serially
                    ExploreASL_Master(DataParFile{iList}{1}, 0, 1, false); % can we run screen from here? or run matlab in background, linux easy
                case 2 % run ExploreASl parallel (start new MATLAB instances)
                    if isunix
                        ScreenString = ['screen -dmS ' ScreenName ' nice -n 10 ' MatlabPath ' -nodesktop -nosplash -r '];
                        RunExploreASLString = ['"cd(''' x.MyPath ''');ExploreASL_Master(''' DataParFile{iList}{1} ''',0,1,0);system([''screen -SX ' ScreenName ' kill'']);"'];
                    else
                        ScreenString = [MatlabPath ' -nodesktop -nosplash -r '];
                        RunExploreASLString = ['"cd(''' x.MyPath ''');ExploreASL_Master(''' DataParFile{iList}{1} ''',0,1,0);system([''exit'']);"'];
                    end
                    system([ScreenString RunExploreASLString ' &']);
                case 3 % run ExploreASL compilation serially
                    system(CompilationString);

                case 4 % run ExploreASL compilation parallel
                    if isunix
                        ScreenString = ['screen -dmS ' ScreenName ' nice -n 10 '];
                    else
                        ScreenString = [];
                    end
                    
                    system([ScreenString CompilationString ' &']);
                otherwise
            end
        catch ME
            warning('Something went wrong:');
            fprintf('%s\', ME.message);
        end
    end
end

% Wait until all *.log files exist (which will surely be created, even with a crash)

% ============================================================
%% 6) Pause until all results exist (if running parallel in background)
if RunMethod==2
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
end

% ============================================================
%% 7) Compile results table
ResultsTable = {'Data', 'mean_qCBF_TotalGM' 'median_qCBF_TotalGM' 'median_qCBF_DeepWM' 'CoV_qCBF_TotalGM' 'GMvol' 'WMvol' 'CSFvol' 'PipelineCompleted' 'TC_ASL_Registration' 'TC_M0_Registration'};
fprintf('Reading & parsing results:   ');
for iList=1:length(Dlist) % iterate over example datasets
    xASL_TrackProgress(iList, length(Dlist));
    ResultsTable{1+iList,1} = Dlist{iList};
    AnalysisDir = fullfile(TestDirDest, Dlist{iList});
    PopulationDir = fullfile(AnalysisDir, 'Population');
    StatsDir = fullfile(PopulationDir, 'Stats');
	VolumeDir = fullfile(PopulationDir, 'TissueVolume');
    
    clear ResultsFile
    ResultFile{1} = xASL_adm_GetFileList(StatsDir,'^mean_qCBF.*TotalGM.*PVC2\.tsv','FPList');
    ResultFile{2} = xASL_adm_GetFileList(StatsDir,'^median_qCBF.*TotalGM.*PVC0\.tsv','FPList');
    ResultFile{3} = xASL_adm_GetFileList(StatsDir,'^median_qCBF.*DeepWM.*PVC0\.tsv','FPList');
    ResultFile{4} = xASL_adm_GetFileList(StatsDir,'^CoV_qCBF.*TotalGM.*PVC0\.tsv','FPList');
    ResultFile{5} = xASL_adm_GetFileList(VolumeDir,'^TissueVolume.*\.tsv','FPList');

	for iFile=1:length(ResultFile) % iterate over ROI results
		if length(ResultFile{iFile})<1
            ResultsTable{1+iList,1+iFile} = 'empty';
            ResultsTable{1+iList,2+iFile} = 'empty';
            ResultsTable{1+iList,3+iFile} = 'empty';
		elseif iFile<5 % check the ASL parameters
            [~, TempTable] = xASL_bids_csv2tsvReadWrite(ResultFile{iFile}{end});
            ResultsTable{1+iList,1+iFile} = TempTable{3,end-2};
		else % check the volumetric parameters
			[~, TempTable] = xASL_bids_csv2tsvReadWrite(ResultFile{iFile}{end});
			% Backward compatibility:
			% Volumetrics used to be saved as '_(L)', but this is converted by spm_jsonread to its HEX counterpart '_0x28L0x29'
			% Now we always save to _L to avoid this. For backward compatibility we still check the old options here
            IndexGM = find(cellfun(@(y) ~isempty(regexpi(y,'(GM_volume_L|GM_volume_(L)|GM_volume_0x28L0x29)')), TempTable(1,:)));
            IndexWM = find(cellfun(@(y) ~isempty(regexpi(y,'(WM_volume_L|WM_volume_(L)|WM_volume_0x28L0x29)')), TempTable(1,:)));
            IndexCSF = find(cellfun(@(y) ~isempty(regexpi(y,'(CSF_volume_L|CSF_volume_(L)|CSF_volume_0x28L0x29)')), TempTable(1,:)));
			if ~isempty(IndexGM)
				ResultsTable{1+iList,1+iFile} = TempTable{2, IndexGM};
			else
				ResultsTable{1+iList,1+iFile} = 'n/a';
			end
			if ~isempty(IndexWM)
				ResultsTable{1+iList,2+iFile} = TempTable{2, IndexWM};
			else
				ResultsTable{1+iList,2+iFile} = 'n/a';
			end
			if ~isempty(IndexCSF)
				ResultsTable{1+iList,3+iFile} = TempTable{2, IndexCSF};
			else
				ResultsTable{1+iList,2+iFile} = 'n/a';
			end
		end
	end
	% check if there are missing lock files
	if exist(fullfile(AnalysisDir,'Missing_Lock_files.csv'),'file')
		ResultsTable{1+iList,4+length(ResultFile)} = 0; % pipeline not completed
	else
		ResultsTable{1+iList,4+length(ResultFile)} = 1; % pipeline completed
	end
	
    % Get registration performance
    PathTemplateASL = fullfile(x.D.TemplateDir, 'Philips_2DEPI_Bsup_CBF.nii');
    PathTemplateM0 = fullfile(x.D.TemplateDir, 'Philips_2DEPI_noBsup_Control.nii');
    PathCBF = xASL_adm_GetFileList(PopulationDir,'^qCBF(?!.*(4D|masked|Visual2DICOM)).*\.nii$', 'FPList');
    PathM0 = xASL_adm_GetFileList(PopulationDir,'^(noSmooth_M0|mean_control).*\.nii$', 'FPList');
    if ~isempty(PathCBF)
        ResultsTable{1+iList,5+length(ResultFile)} = xASL_qc_TanimotoCoeff(PathCBF{1}, PathTemplateASL, x.S.masks.WBmask, 3, 0.975, [4 0]); % Tanimoto Coefficient, Similarity index
    end
    if ~isempty(PathM0)
        ResultsTable{1+iList,6+length(ResultFile)} = xASL_qc_TanimotoCoeff(PathM0{1}, PathTemplateM0, x.S.masks.WBmask, 3, 0.975, [4 0]); % Tanimoto Coefficient, Similarity index
    end
end
fprintf('\n');

% Save results
ResultTableName = datestr(now,'yyyy-mm-dd_HH_MM');
ResultTableName(end-2) = 'h';
ResultTableFile = [ResultTableName,'_ResultsTable.mat'];
SaveFile = fullfile(TestDirOrig, ResultTableFile);
save(SaveFile, 'ResultsTable');

% ============================================================
%% 8) Compare table with reference table

% Comparison with tsv file
[ReferenceTables,ReferenceTable] = xASL_qc_LoadRefTable(fullfile(x.MyPath,'Testing','Reference','ReferenceValues.tsv'));
ResultsComparison = xASL_qc_CompareTables(ReferenceTable,ResultsTable);
save(SaveFile, 'ResultsTable', 'ReferenceTables', 'ReferenceTable', 'ResultsComparison');

% Comparison with mat file
try
    % Find all result tables in directory
    AllResultTables = xASL_adm_GetFsList(fullfile(TestDirOrig),'^.+\_ResultsTable.mat$',false)';
    IndexCurrentTable = find(~isempty(strfind(AllResultTables, ResultTableFile)));
    if ~isempty(IndexCurrentTable)
        AllResultTables(IndexCurrentTable) = [];
    end
    % Previous Table name
    if ~isempty(AllResultTables)
        PreviousTableName = AllResultTables{numel(AllResultTables)};
        PreviousSaveFile = fullfile(TestDirOrig, PreviousTableName);
        PreviousTable = load(PreviousSaveFile, '-mat');

        clear DifferenceTable
        DifferenceTable(1:size(ResultsTable,1),1:size(ResultsTable,2)) = {''};
        DifferenceTable(1,:) = ResultsTable(1,:);
        DifferenceTable(:,1) = ResultsTable(:,1);
        for iX=2:size(ResultsTable,1)
            for iY=2:size(ResultsTable,2)
                A = xASL_str2num(ResultsTable{iX,iY});
                B = xASL_str2num(PreviousTable.ResultsTable{iX,iY});
                AsymmIndex = (A-B)/(A+B);
                DifferenceTable{iX,iY} = [xASL_num2str(AsymmIndex) '%'];
            end
        end
    else
        % Could not find previous tables
        DifferenceTable = [];
    end
    
    % ============================================================
    %% 9) E-mail results
    if ~isempty(EmailAddress)
        % First convert table to string to send by e-mail
        NewTable{1,1} = 'mean_qCBF_TotalGM     median_qCBF_TotalGM     median_qCBF_DeepWM     CoV_qCBF_TotalGM             GMvol                 WMvol                 CSFvol             PipelineCompleted     TC_ASL_Registration    TC_M0_Registration';
        SingleEmptyString1 = repmat(' ',[1 44]);  
        SingleEmptyString2 = repmat(' ',[1 27]);
        for iX=2:size(DifferenceTable,1)
            NewTable{iX,1} = '';
            for iY=2:size(DifferenceTable,2)
                if iY<5
                    NewCell = SingleEmptyString1;
                else
                    NewCell = SingleEmptyString2;
                end

                NewCell(1:length(DifferenceTable{iX,iY})) = DifferenceTable{iX,iY};
                NewTable{iX,1} = [NewTable{iX,1} NewCell];
            end
            NewTable{iX,1} = [NewTable{iX,1} DifferenceTable{iX,1}];
        end


        % See here: https://nl.mathworks.com/help/matlab/import_export/sending-email.html
        fprintf('Sending e-mail with results\n');
        setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
        setpref('Internet', 'E_mail', EmailAddress);
        setpref('Internet', 'SMTP_Username', EmailAddress);
        setpref('Internet', 'SMTP_Password', Password);
        props = java.lang.System.getProperties;
        props.setProperty('mail.smtp.auth','true');
        props.setProperty('mail.smtp.starttls.enable', 'true');
        props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
        props.setProperty('mail.smtp.socketFactory.port','465'); % or port 587
        EmailAddresses = {'Patricia.Clement@ugent.be', 'Pieter.Vandemaele@UZGENT.be', 'j.petr@hzdr.de', 'henkjanmutsaerts@gmail.com'};
        sendmail(EmailAddresses, 'ExploreASL TestRun: %AsymmetryIndexWithTemplateResults (should be <0.01%)', NewTable);
    end
catch ME
    warning('Something went wrong in trying to create difference table & mailing it to receivers');
    fprintf('%s\n', ME.message);
end
    

end

% Load Reference Table
function [ReferenceTables,ReferenceTable] = xASL_qc_LoadRefTable(pathRefTable)

    % Load TSV file
    ReferenceTables = xASL_tsvRead(pathRefTable);
    
    iRow = 1;
    while iRow<=size(ReferenceTables,1)
        if ~isempty(regexp(ReferenceTables{iRow,1},'^xASL_', 'once'))
            versionXASL = ReferenceTables{iRow,1};
            operatingSystem = ReferenceTables{iRow+1,1};
            fprintf('Version: %s\n', versionXASL);
            fprintf('OS:      %s\n', operatingSystem);
            ReferenceTables(iRow,:) = []; % Remove row 1
            ReferenceTables(iRow,:) = []; % Remove row 2
            ReferenceTable.(versionXASL) = ReferenceTables(iRow:iRow+10,:);
        end
        iRow=iRow+1;
    end

end

% Compare Results and Reference Tables
function ResultsComparison = xASL_qc_CompareTables(ReferenceTable,ResultsTable)

    % Compare tables (skip first row)
    ResultsComparison = ReferenceTable;
    
    % Iterate over versions
    versionsXASL = fieldnames(ResultsComparison);
    for iVersion = 1:size(versionsXASL,1)
        % Get current reference table
        currentReferenceTable = ReferenceTable.(versionsXASL{iVersion});
        
        for iRow = 2:size(currentReferenceTable,1)
            ResultsComparison.(versionsXASL{iVersion})(iRow,2:end) = {NaN};
            if strcmp(currentReferenceTable{iRow,1},ResultsTable{iRow,1})
                % Iterate over results
                for iColumn = 2:size(currentReferenceTable,2)
                    refValue = xASL_str2num(currentReferenceTable{iRow,iColumn});
                    resValue = xASL_str2num(ResultsTable{iRow,iColumn});
					if isnan(refValue) && isnan(resValue)
						% NaNs are assigned to strings, if both are strings or NaNs then the comparison is fine
						ResultsComparison.(versionsXASL{iVersion})(iRow,iColumn) = {true};
					else
						% Check if difference is smaller than 0.1% of reference value
						if abs(refValue-resValue)<(abs(refValue)*0.01)
							ResultsComparison.(versionsXASL{iVersion})(iRow,iColumn) = {true};
						else
							ResultsComparison.(versionsXASL{iVersion})(iRow,iColumn) = {false};
						end
					end
                end
            else
                fprintf('Dataset names do not match...\n');
                fprintf('Reference: %s\n', currentReferenceTable{iRow,1});
                fprintf('Results:   %s\n', ResultsTable{iRow,1});
            end
        end
    end
    
    % Check: single value passed/failed
    versionsXASL = fieldnames(ResultsComparison);
    for iVersion = 1:size(versionsXASL,1)
        % Get current reference table
        currentResults = ResultsComparison.(versionsXASL{iVersion});
        resultValues = currentResults(2:end,2:end);
        % Check if there are "failed" values
        passed = true;
        for iRow=1:size(resultValues,1)
            for iColumn=1:size(resultValues,2)
                %fprintf('%d\n', resultValues{iRow,iColumn});
                if resultValues{iRow,iColumn}==0
                    passed = false; % At least one value is "failed"
                end
            end
        end
        % Store value
        ResultsComparison.([versionsXASL{iVersion} '_Passed']) = passed;
        if passed
            fprintf('Version: %s\nPassed:  %s\n', versionsXASL{iVersion}, 'true');
        else
            fprintf('Version: %s\nPassed:  %s\n', versionsXASL{iVersion}, 'false');
        end
    end



end


