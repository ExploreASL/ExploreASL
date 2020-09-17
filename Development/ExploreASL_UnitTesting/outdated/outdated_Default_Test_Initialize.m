function x = xASL_qc_Default_Test_Initialize(x,rootDir,subject,moduleName)
%xASL_qc_Iteration_Initialize Script to run the xASL_Iteration initialization.
%
% FORMAT:       x = xASL_qc_Default_Test_Initialize(x);
% 
% INPUT:        x           structure containing all the initialized fields
%               rootDir     root directory of the test folder / test data set
%               subject     Sub-001 for example
%               moduleName  xASL_module_Structural for example
%
% OUTPUT:       x structure
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  My version of the initialization part of xASL_Iteration for
%               the QC testing of individual modules and submodules.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES:     x = xASL_qc_Default_Test_Initialize(x,TestDataSetPath,x.SUBJECTS{1});
% __________________________________
% Copyright 2015-2020 ExploreASL

% Use the Structural module in the default case
if nargin<4
    moduleName = 'xASL_module_Structural';
end

% General settings
dbSettings.x                     = x;
dbSettings.x.RERUN               = false;
dbSettings.x.MUTEXID             = moduleName;
dbSettings.x.LockDir             = fullfile(rootDir, 'lock', moduleName);

% Extract the ModName
if  strcmp(moduleName(1:12),'xASL_module_')
    ModName     = moduleName(13:end);
else
    ModName     = moduleName;
end

% Specific settings
if  strcmp(ModName,'DARTEL') || strcmp(ModName,'LongReg')
    dbSettings.x.LockDir         = [dbSettings.x.LockDir '_' x.P.STRUCT ];
end

if ~isempty(regexp(ModName,'(Struct|ASL|func|LongReg|dwi)'))
    dbSettings.x.SUBJECTDIR      = fullfile(rootDir,subject);
    dbSettings.x.LockDir         = fullfile(dbSettings.x.LockDir, subject);
end

if ~isempty(regexp(ModName,'(ASL|func|dwi)'))
    dbSettings.x.MUTEXID         = [dbSettings.x.MUTEXID '_<SESSION>'];
    dbSettings.x.SESSIONDIR      = '<ROOT>/<SUBJECT>/<SESSION>';
end

% Write the dbSetting back to the x structure
x = dbSettings.x;


end









