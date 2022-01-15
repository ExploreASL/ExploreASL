function xASL_adm_MakeStandalone(outputPath, bCompileSPM, markAsLatest)
%xASL_adm_MakeStandalone This function was written to create a compiled "standalone" version of
% ExploreASL using the mcc compiler from Matlab.
%
% FORMAT: xASL_adm_MakeStandalone(outputPath, bCompileSPM, importDCM, markAsLatest);
%
% INPUT:
%   outputPath      - Folder where the compiled version should be saved (REQUIRED)
%   bCompileSPM     - Boolean specifying whether SPM is compiled first (OPTIONAL, DEFAULT = true)
%   markAsLatest    - Option to mark the generated compiled versions as "latest" instead, to simplify 
%                     the docker integration for example (OPTIONAL, DEFAULT = true)
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
% 8. Print done
%
% EXAMPLE: 
%
% 1. Compilation marked as latest including compiled SPM: xASL_adm_MakeStandalone('Drive/User/Folder/',true,true);
% 2. Compilation marked with version number including compiled SPM: xASL_adm_MakeStandalone('Drive/User/Folder/',true,false);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL


%% 1) Manage ExploreASL and compiler code folders
if nargin<1 || isempty(outputPath);     error('OutputPath missing...');     end
if nargin<2 || isempty(bCompileSPM);    bCompileSPM = true;                 end
if nargin<3 || isempty(markAsLatest);   markAsLatest = true;                end

% Mana file handle issues
close all

% Initialize Explore ASL
x = ExploreASL_Initialize;

%% 2) Define Versioning
VersionPath = xASL_adm_GetFileList(x.opts.MyPath, '^VERSION_.*', 'List', [0 Inf]);
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

%% 3) File management output folder & starting diary
pathLogFile = fullfile(outputPath,'compilationLog.txt');
outputPath = fullfile(outputPath,Version);
if ~exist(outputPath, 'dir')
    xASL_adm_CreateDir(outputPath);
else
    warning(['Compilation output folder ' outputPath ' already existed, will now delete all files in this folder recursively']);
    fprintf('%s\n','Press any key to start processing & analyzing');
    fprintf('%s\n','Or press CTRL/command-C to cancel...  ');
    pause;
    xASL_adm_DeleteFileList(outputPath, '.*', true, [0 Inf]);
end
diary(pathLogFile);
diary on

%% 4) Handle SPM Specific Options
% Static listing of batch application initialisation files
cfg_util('dumpcfg');

% Duplicate Contents.m in Contents.txt for use in spm('Ver')
xASL_Copy(fullfile(spm('Dir'),'Contents.m'), fullfile(spm('Dir'),'Contents.txt'));

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
    rmpath(genpath(fullfile(x.opts.MyPath, FoldersNot2Deploy{iFolder})));
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
    zip([fullfile(DummyDir) '.zip'],fullfile(DummyDir));
    xASL_delete(DummyDir, true);
end

%% 7) Run ExploreASL compilation
fprintf('Compiling ExploreASL\n');

if ~isempty(VersionPath)
    AddExploreASLversion = fullfile(x.opts.MyPath,VersionPath{1});
else
    AddExploreASLversion = '';
end

% Compilation
mcc('-m', '-C', '-v',... % '-R -nodisplay -R -softwareopengl',... % https://nl.mathworks.com/matlabcentral/answers/315477-how-can-i-compile-a-standalone-matlab-application-with-startup-options-e-g-nojvm
    fullfile(x.opts.MyPath,'ExploreASL.m'),...
    '-d', fullfile(outputPath),...
    '-o', Version,...
    '-N', opts{:},...
    '-a', spm('Dir'),... % For SPM support
    '-a', AddExploreASLversion,...
    '-a', fullfile(x.opts.MyPath,'Functions'),...
    '-a', fullfile(x.opts.MyPath,'Functions','Administration'),...
    '-a', fullfile(x.opts.MyPath,'Functions','BIDS'),...
    '-a', fullfile(x.opts.MyPath,'Functions','FSL'),...
    '-a', fullfile(x.opts.MyPath,'Functions','ImageProcessing'),...
    '-a', fullfile(x.opts.MyPath,'Functions','Import'),...
    '-a', fullfile(x.opts.MyPath,'Functions','Initialization'),...
    '-a', fullfile(x.opts.MyPath,'Functions','InputOutput'),...
    '-a', fullfile(x.opts.MyPath,'Functions','mex'),...
    '-a', fullfile(x.opts.MyPath,'Functions','QualityControl'),...
    '-a', fullfile(x.opts.MyPath,'Functions','Quantification'),...
    '-a', fullfile(x.opts.MyPath,'Functions','SPM'),...
    '-a', fullfile(x.opts.MyPath,'Functions','Statistics'),...
    '-a', fullfile(x.opts.MyPath,'Functions','Visualization'),...
    '-a', fullfile(x.opts.MyPath,'Maps'),...
    '-a', fullfile(x.opts.MyPath,'Modules'),... % Modules
    '-a', fullfile(x.opts.MyPath,'Modules', 'SubModule_Import'),...
    '-a', fullfile(x.opts.MyPath,'Modules', 'SubModule_Structural'),...
    '-a', fullfile(x.opts.MyPath,'Modules', 'SubModule_ASL'),...
    '-a', fullfile(x.opts.MyPath,'Modules', 'SubModule_Population'),...
    '-a', fullfile(x.opts.MyPath,'External','MRIcron'),... % DICOM import
    '-a', fullfile(x.opts.MyPath,'External', 'isnear'),...
    '-a', fullfile(x.opts.MyPath,'External', 'Atlases'));

% Copy version file to compilation folder
xASL_Copy(AddExploreASLversion, fullfile(outputPath, VersionPath{1}), 1);

% Zip the compilation
zip([fullfile(outputPath) '.zip'],fullfile(outputPath));

% Close diary
diary(pathLogFile);
diary off

% Delete unzipped folder
xASL_delete(outputPath, true);

%% 8) Print done
fprintf('Done\n');

end


