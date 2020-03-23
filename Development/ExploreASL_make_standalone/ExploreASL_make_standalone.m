function ExploreASL_make_standalone(outputPath, bCompileSPM)
%% ExploreASL_make_standalone
% This function was written to create a compiled "standalone" version of
% ExploreASL using the mcc compiler from Matlab.
%
% INPUT:
%   outputPath      - folder where the compiled version should be saved (REQUIRED)
%   bCompileSPM     - boolean specifying whether SPM is compiled first
%                     (OPTIONAL, DEFAULT=true)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates an output folder including a
%              standalone version of ExploreASL, which can be used with the Matlab
%              Runtime outside of Matlab itself.
%              
%              A quick fix to solve path dependencies etc. is to first
%              compile SPM.
%
% EXAMPLE: ExploreASL_make_standalone('/Path2/Compilation')
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2019 ExploreASL




%% 1) Manage ExploreASL and compiler code folders
if nargin<1 || isempty(outputPath)
    error('OutputPath input missing\n');
end
if nargin<2 || isempty(bCompileSPM)
    bCompileSPM = true;
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

% currentPath = mfilename('fullpath');
% currentPath = strsplit(currentPath,filesep);
% [~,id] = intersect(currentPath,'ExploreASL','stable');
% currentPath = string(currentPath(1,1:id));
% ExploreASLPath = '';
% for n=1:length(currentPath)
%     ExploreASLPath = strcat(ExploreASLPath,currentPath(n),filesep);
% end
% clear n id currentPath



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

Version = xASL_adm_CorrectName(['xASL' xASLVersion '_' MVersion '_' date '_' Time]);



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

% % First back up the cfg file
% CfgPath = fullfile(spm('Dir'), 'matlabbatch', 'private', 'cfg_mlbatch_appcfg_master.m');
% CfgPath2 = fullfile(spm('Dir'), 'matlabbatch', 'private', 'cfg_mlbatch_appcfg_1.m');
% BackupCfgPath = fullfile(outputPath, 'cfg_mlbatch_appcfg_master.m');
% BackupCfgPath2 = fullfile(outputPath, 'cfg_mlbatch_appcfg_1.m');
% if exist(CfgPath, 'file')
%     xASL_Copy(CfgPath, BackupCfgPath, true);
% end
% if exist(CfgPath2, 'file')
%     xASL_Copy(CfgPath2, BackupCfgPath2, true);
% end

% Static listing of SPM toolboxes
% fid = fopen(fullfile(spm('Dir'),'config','spm_cfg_static_tools.m'),'wt');
% fprintf(fid,'function values = spm_cfg_static_tools\n');
% fprintf(fid,...
%     '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
% % Get the list of toolbox directories
% tbxdir = fullfile(spm('Dir'),'toolbox');
% d = [tbxdir; cellstr(spm_select('FPList',tbxdir,'dir'))];
% ft = {};
% % Look for '*_cfg_*.m' files in these directories
% for i=1:numel(d)
%     fi = spm_select('List',d{i},'.*_cfg_.*\.m$');
%     if ~isempty(fi) && isempty(regexp(fi,'cfg_(fieldmap|render)')) % remove fieldmap & render, we don't use these, to avoid conflicts
%         ft = [ft(:); cellstr(fi)];
%     end
% end
% % Create code to insert toolbox config
% if isempty(ft)
%     ftstr = '';
% else
%     ft = spm_file(ft,'basename');
%     ftstr = sprintf('%s ', ft{:});
% end
% fprintf(fid,'values = {%s};\n', ftstr);
% fclose(fid);

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



%% 6) timeFirst run SPM compilation, see if this helps remove any errors (e.g. path dependency-issues)...
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


%% Run ExploreASL compilation
fprintf('Compiling ExploreASL\n');
mcc('-m', '-C', '-v',... % '-R -nodisplay -R -softwareopengl',... % https://nl.mathworks.com/matlabcentral/answers/315477-how-can-i-compile-a-standalone-matlab-application-with-startup-options-e-g-nojvm
    fullfile(ExploreASLPath,'ExploreASL_Master.m'),...
    '-d', fullfile(outputPath),...
    '-o', strcat('ExploreASL_',Version),...
    '-N', opts{:},...                                   % Added for SPM support
    '-a', spm('Dir'),...                                % Added for SPM support
    '-a', AddExploreASLversion,...
    '-a', fullfile(ExploreASLPath,'Functions'),...
    '-a', fullfile(ExploreASLPath,'External','SPMmodified','xASL'),...
    '-a', fullfile(ExploreASLPath,'Maps'),...
    '-a', fullfile(ExploreASLPath,'mex'),...
    '-a', fullfile(ExploreASLPath,'Modules'));

% put opengl('software') in code?

%% Copy .bat file for Windows compilation
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



%% Compilation successful, save Log-file
fprintf('Done\n');
diary(fullfile(outputPath,'compilationLog.txt'));
diary off


end
