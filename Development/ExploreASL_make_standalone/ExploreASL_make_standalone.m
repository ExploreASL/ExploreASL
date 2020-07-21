function ExploreASL_make_standalone(outputPath, bCompileSPM, bRelease, importDCM)
%ExploreASL_make_standalone This function was written to create a compiled "standalone" version of
% ExploreASL using the mcc compiler from Matlab.
%
% INPUT:
%   outputPath      - Folder where the compiled version should be saved (REQUIRED)
%   bCompileSPM     - Boolean specifying whether SPM is compiled first
%                     (OPTIONAL, DEFAULT=true)
%   bRelease        - Set to true for release compilations. Changes the filename.
%                     DEFAULT=false
%   importDCM       - Generate a separate standalone import for DICOM2BIDS.
%                     (OPTIONAL, DEFAULT=false)
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
% EXAMPLE: ExploreASL_make_standalone('/Path2/Compilation')
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL




%% 1) Manage ExploreASL and compiler code folders
if nargin<1 || isempty(outputPath)
    error('OutputPath input missing\n');
end
if nargin<2 || isempty(bCompileSPM)
    bCompileSPM = true;
end
if nargin<3
    bRelease = false;
end
if nargin<4
    importDCM = false;
end

try
    ExploreASL_Master('',0);
catch
    cd ../..
    ExploreASL_Master('',0);
end

OldPath = 'RunExploreASL.bat';

if ~exist(OldPath, 'file')
    try
        cd('Development');
        cd('ExploreASL_make_standalone');
    catch ME
        warning('Couldnt access standalone compiler folder');
        fprintf('%s\n', ME.message);
    end
    if ~exist(OldPath, 'file')
        error('Please start this script in its own folder');
    end
end

CurrDir = fileparts(mfilename('fullpath'));
ExploreASLPath = fileparts(fileparts(CurrDir)); % assuming to folder layers


%% 2) Capture version/date/time
Time = clock;
Time = [num2str(round(Time(4))) 'h' num2str(round(Time(5))) 'm'];
VersionPath = xASL_adm_GetFileList(ExploreASLPath, '^VERSION_.*', 'List', [0 Inf]);
if ~isempty(VersionPath)
    xASLVersion = [VersionPath{1}(9:end) '_'];
else
    xASLVersion = '';
end

MVersion = version;
[StartIndex, EndIndex] = regexp(MVersion,'R\d{4}(a|b)');
if ~isempty(StartIndex)
    MVersion = MVersion(StartIndex:EndIndex);
else
    warning('Couldnt find Matlab version');
    MVersion = '';
end

% Different notation for compiled release version
if bRelease
    Version = xASL_adm_CorrectName(['xASL_' xASLVersion '_Release']);
else
    Version = xASL_adm_CorrectName(['xASL_' xASLVersion '_' MVersion '_' date '_' Time]);
end

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
    warning('Copy of Contents.m failed.');
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
    DummyDir = fullfile(fileparts(outputPath), 'DummySPM_CompilationTest');
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
    '-o', strcat('ExploreASL_',Version),...
    '-N', opts{:},...
    '-a', spm('Dir'),... % For SPM support
    '-a', AddExploreASLversion,...
    '-a', fullfile(ExploreASLPath,'Functions'),...
    '-a', fullfile(ExploreASLPath,'mex'),...
    '-a', fullfile(ExploreASLPath,'Maps'),...
    '-a', fullfile(ExploreASLPath,'Modules'),...
    '-a', fullfile(ExploreASLPath,'Modules', 'SubModule_Structural'),...
    '-a', fullfile(ExploreASLPath,'Modules', 'SubModule_ASL'),...
    '-a', fullfile(ExploreASLPath,'Modules', 'SubModule_Population'),...
    '-a', fullfile(ExploreASLPath,'External','isnear'),...
    '-a', fullfile(ExploreASLPath,'External','AtlasesNonCommercial'));

% Compilation DICOM import (Work in progress -> Meant to be used for docker integration)
if importDCM
    mcc('-m', '-C', '-v',... % '-R -nodisplay -R -softwareopengl',... % https://nl.mathworks.com/matlabcentral/answers/315477-how-can-i-compile-a-standalone-matlab-application-with-startup-options-e-g-nojvm
    fullfile(ExploreASLPath,'Functions','xASL_dcm_Import'),...
    '-d', fullfile(outputPathImport),...
    '-o', strcat('ExploreASL_',VersionImport),...
    '-N', opts{:},...
    '-a', spm('Dir'),... % For SPM support
    '-a', AddExploreASLversion,...
    '-a', fullfile(ExploreASLPath,'ExploreASL_Master.m'),...
    '-a', fullfile(ExploreASLPath,'ExploreASL_Initialize.m'),...
    '-a', fullfile(ExploreASLPath,'ExploreASL_Import.m'),...
    '-a', fullfile(ExploreASLPath,'ExploreASL_ImportConfig.m'),...
    '-a', fullfile(ExploreASLPath,'Development','xASL_par_Fix.m'),... % Fix JSON files
    '-a', fullfile(ExploreASLPath,'Functions'),...
    '-a', fullfile(ExploreASLPath,'mex'),...
    '-a', fullfile(ExploreASLPath,'Modules'),...
    '-a', fullfile(ExploreASLPath,'External','MRIcron')); % DICOM import
end

%% 8) Copy .bat file for Windows compilation
if ispc
    NewPath = fullfile(outputPath, 'RunExploreASL.bat');
    xASL_Copy(OldPath, NewPath, true);
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
