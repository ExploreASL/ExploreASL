function xASL_adm_MakeStandalone(outputPath, bCompileSPM, importDCM, markAsLatest)
%xASL_adm_MakeStandalone This function was written to create a compiled "standalone" version of
% ExploreASL using the mcc compiler from Matlab.
%
% FORMAT: xASL_adm_MakeStandalone(outputPath, bCompileSPM, importDCM, markAsLatest);
%
% INPUT:
%   outputPath      - Folder where the compiled version should be saved (REQUIRED)
%   bCompileSPM     - Boolean specifying whether SPM is compiled first
%                     (OPTIONAL, DEFAULT=true)
%   importDCM       - Generate a separate standalone import for DICOM2BIDS.
%                     (OPTIONAL, DEFAULT=true)
%   markAsLatest    - Option to mark the generated compiled versions as
%                     "latest" instead, to simplify the docker integration for example.
%                     (OPTIONAL, DEFAULT=true)
%
% OUTPUT:       Generates a standalone/executable version of ExploreASL.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates an output folder including a
%              standalone version of ExploreASL, which can be used with the Matlab
%              Runtime outside of Matlab itself.
%              
%              A quick fix to solve path dependencies etc. is to first
%              compile SPM (but this can be turned off for speed).
%
% This function performs the following steps:
%
% 1. Manage ExploreASL and compiler code folders
% 2. Capture version/date/time
% 3. File management output folder & starting diary
% 4. Handle SPM Specific Options
% 5. Manage compilation paths
% 6. Run SPM compilation
% 7. Run ExploreASL compilation
% 8. Copy .bat file for Windows compilation
% 9. Save Log-file
%
% EXAMPLE: xASL_adm_MakeStandalone('/Path2/Compilation')
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL


%% 1) Manage ExploreASL and compiler code folders
if nargin<1 || isempty(outputPath),     error('OutputPath input missing');  end
if nargin<2 || isempty(bCompileSPM),    bCompileSPM = true;                 end
if nargin<3,                            importDCM = true;                   end
if nargin<4,                            markAsLatest = true;                end

% Initialize Explore ASL
x = ExploreASL_Initialize([],0);

% Get current directory
ExploreASLPath = x.MyPath;

%% 2) Define Versioning
VersionPath = xASL_adm_GetFileList(ExploreASLPath, '^VERSION_.*', 'List', [0 Inf]);
if ~isempty(VersionPath)
    xASLVersion = VersionPath{1}(9:end);
else
    xASLVersion = '';
end

if markAsLatest
    xASLVersion = 'latest';
end

% Define file name (add version)
Version = xASL_adm_CorrectName(['xASL_' xASLVersion]);

% Version name of standalone DICOM import
if importDCM
    VersionImport = xASL_adm_CorrectName(['xASL_' xASLVersion '_Import']);
    % Output folder
    outputPathImport = fullfile(outputPath,VersionImport);
    if ~exist(outputPathImport, 'dir')
        mkdir(outputPathImport);
    end
end

%% 3) File management output folder & starting diary
outputPath = fullfile(outputPath,Version);
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
else
    warning(['Compilation output folder ' outputPath ' already existed, will now delete all files in this folder recursively']);
    fprintf('%s\n','Press any key to start processing & analyzing');
    fprintf('%s\n','Or press CTRL/command-C to cancel...  ');
    pause;
    xASL_adm_DeleteFileList(outputPath, '.*', true, [0 Inf]);
end
diary(fullfile(outputPath,'compilationLog.txt'));
diary on



%% 4) Handle SPM Specific Options
% Static listing of batch application initialisation files
cfg_util('dumpcfg');

% Duplicate Contents.m in Contents.txt for use in spm('Ver')
sts = copyfile(fullfile(spm('Dir'),'Contents.m'), fullfile(spm('Dir'),'Contents.txt'));
               
if ~sts
    warning('Copy of Contents.m failed');
end


%% 5) Manage compilation paths
% It is important that there aren't any syntax errors in the scripts of the
% folders below (Development, Functions, etc.). To include all files of the
% described folders is necessary, so that no important scripts are missing
% in the standalone version. Functions which are executed by using eval for
% example can be missed by the automatic matlab mcc search trees, which is
% why I chose this method.

% Add all SPM folders
addpath(genpath(fullfile(spm('Dir'))));

% First remove folders that we want to exclude
warning('off','MATLAB:rmpath:DirNotFound');
FoldersNot2Deploy = {'CustomScripts','Development','GitHub','TestDataSet',fullfile('External','DIP'),...
    fullfile(spm('Dir'),'toolbox', 'FieldMap'),fullfile(spm('Dir'),'toolbox', 'SRender')};
for iFolder=1:length(FoldersNot2Deploy)
    rmpath(genpath(fullfile(ExploreASLPath, FoldersNot2Deploy{iFolder})));
end
warning('on','MATLAB:rmpath:DirNotFound');

opts = {'-p',fullfile(matlabroot,'toolbox','signal')};
if ~exist(opts{2},'dir'); opts = {}; end



%% 6) Run SPM compilation
if bCompileSPM
    fprintf('First compiling SPM as test\n');
    spmCompilationName = xASL_adm_CorrectName(['xASL_' xASLVersion '_SPM_Dummy']);
    DummyDir = fullfile(fileparts(outputPath), spmCompilationName);
    xASL_adm_CreateDir(DummyDir);
    spm_make_standalone(DummyDir);
end

if ~isempty(VersionPath)
    AddExploreASLversion = fullfile(ExploreASLPath,VersionPath{1});
else
    AddExploreASLversion = '';
end


%% 7) Run ExploreASL compilation
fprintf('Compiling ExploreASL\n');

% Compilation
mcc('-m', '-C', '-v',... % '-R -nodisplay -R -softwareopengl',... % https://nl.mathworks.com/matlabcentral/answers/315477-how-can-i-compile-a-standalone-matlab-application-with-startup-options-e-g-nojvm
    fullfile(ExploreASLPath,'ExploreASL_Master.m'),...
    '-d', fullfile(outputPath),...
    '-o', Version,...
    '-N', opts{:},...
    '-a', spm('Dir'),... % For SPM support
    '-a', AddExploreASLversion,...
    '-a', fullfile(ExploreASLPath, 'Functions'),...
    '-a', fullfile(ExploreASLPath, 'mex'),...
    '-a', fullfile(ExploreASLPath, 'Maps'),...
    '-a', fullfile(ExploreASLPath, 'Modules'),...
    '-a', fullfile(ExploreASLPath, 'Modules', 'SubModule_Structural'),...
    '-a', fullfile(ExploreASLPath, 'Modules', 'SubModule_ASL'),...
    '-a', fullfile(ExploreASLPath, 'Modules', 'SubModule_Population'),...
    '-a', fullfile(ExploreASLPath, 'External', 'isnear'),...
    '-a', fullfile(ExploreASLPath, 'External', 'Atlases', 'License_CC_BY_4_0', 'Mindboggle-101', 'OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_v2_1_5mm.tsv'),...
    '-a', fullfile(ExploreASLPath, 'External', 'Atlases', 'License_CC_BY_4_0', 'Mindboggle-101', 'OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_v2_1_5mm.nii'),...
    '-a', fullfile(ExploreASLPath, 'External', 'Atlases', 'License_Restricted_Use'));

% Compilation DICOM import (Work in progress -> Meant to be used for docker integration)
if importDCM
    mcc('-m', '-C', '-v',... % '-R -nodisplay -R -softwareopengl',... % https://nl.mathworks.com/matlabcentral/answers/315477-how-can-i-compile-a-standalone-matlab-application-with-startup-options-e-g-nojvm
    fullfile(ExploreASLPath,'External','mediri','xASL_io_mTRIAL'),...
    '-d', fullfile(outputPathImport),...
    '-o', VersionImport,...
    '-N', opts{:},...
    '-a', spm('Dir'),... % For SPM support
    '-a', AddExploreASLversion,...
    '-a', fullfile(ExploreASLPath,'ExploreASL_Master.m'),...
    '-a', fullfile(ExploreASLPath,'ExploreASL_Initialize.m'),...
    '-a', fullfile(ExploreASLPath,'ExploreASL_Import.m'),...
    '-a', fullfile(ExploreASLPath,'ExploreASL_ImportConfig.m'),...
    '-a', fullfile(ExploreASLPath,'Functions'),...
    '-a', fullfile(ExploreASLPath,'mex'),...
    '-a', fullfile(ExploreASLPath,'Modules'),...
    '-a', fullfile(ExploreASLPath,'External', 'MRIcron')); % DICOM import
end

%% 8) Copy .bat file for Windows compilation
createBat = false;
if ispc && createBat
    NewPath = fullfile(outputPath, 'RunExploreASL.bat');
    xASL_Copy(outputPath, NewPath, true);
end

%% Restore cfg file
% if exist(BackupCfgPath, 'file')
%     xASL_Move(BackupCfgPath, CfgPath, true);
% end
% if exist(BackupCfgPath2, 'file')
%     xASL_Move(BackupCfgPath2, CfgPath2, true);
% end

%% 9) Save Log-file
fprintf('Done\n');
diary(fullfile(outputPath,'compilationLog.txt'));
diary off


end

% Back up function
function backUpConfigFile(outputPath)

    % First back up the cfg file
    CfgPath = fullfile(spm('Dir'), 'matlabbatch', 'private', 'cfg_mlbatch_appcfg_master.m');
    CfgPath2 = fullfile(spm('Dir'), 'matlabbatch', 'private', 'cfg_mlbatch_appcfg_1.m');
    BackupCfgPath = fullfile(outputPath, 'cfg_mlbatch_appcfg_master.m');
    BackupCfgPath2 = fullfile(outputPath, 'cfg_mlbatch_appcfg_1.m');
    if exist(CfgPath, 'file')
        xASL_Copy(CfgPath, BackupCfgPath, true);
    end
    if exist(CfgPath2, 'file')
        xASL_Copy(CfgPath2, BackupCfgPath2, true);
    end

    % Static listing of SPM toolboxes
    fid = fopen(fullfile(spm('Dir'),'config','spm_cfg_static_tools.m'),'wt');
    fprintf(fid,'function values = spm_cfg_static_tools\n');
    fprintf(fid,...
        '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
    % Get the list of toolbox directories
    tbxdir = fullfile(spm('Dir'),'toolbox');
    d = [tbxdir; cellstr(spm_select('FPList',tbxdir,'dir'))];
    ft = {};
    % Look for '*_cfg_*.m' files in these directories
    for i=1:numel(d)
        fi = spm_select('List',d{i},'.*_cfg_.*\.m$');
        if ~isempty(fi) && isempty(regexp(fi,'cfg_(fieldmap|render)')) % remove fieldmap & render, we don't use these, to avoid conflicts
            ft = [ft(:); cellstr(fi)];
        end
    end
    % Create code to insert toolbox config
    if isempty(ft)
        ftstr = '';
    else
        ft = spm_file(ft,'basename');
        ftstr = sprintf('%s ', ft{:});
    end
    fprintf(fid,'values = {%s};\n', ftstr);
    fclose(fid);

end

