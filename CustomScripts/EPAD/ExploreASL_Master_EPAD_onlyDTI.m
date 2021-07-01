function [x] = ExploreASL_Master_EPAD_onlyDTI(DataParPath, ProcessData, SkipPause, iWorker, nWorkers)

%%%%ADAPTATION OF EXPLOREASL_MASTER_EPAD FOR RUNNING ONLY DTI ANALYSIS %%%%
%%% WE WANT TO DENOISE DTI DATA (MRTRIX) AND THEN RUN xASL, THIS IS A TEST
%%% SCRIPT %%%%%


%ExploreASL_Master_EPAD Run ExploreASL_Master without population module, and add xASL_module_func &
%xASL_module_dwi

    % -----------------------------------------------------------------------------
    %% 1 Initialization
    % -----------------------------------------------------------------------------

    % NB: *.mat files that contain statistics in data-root folder
    % (e.g.\analysis) should contain the following format:

    % 1st column should be subject-id column
    % 2nd column should be parameter values
    % or 2nd column is session-id (e.g. 'ASL_1', then
    % 3rd column contains parameter values

    % using exist(var) here as nargin doesnt work when debugging
    if ~exist('nWorkers','var') || isempty(nWorkers)
        nWorkers = 1;
    end
    if ~exist('iWorker','var') || isempty(iWorker)
        iWorker = 1;
    end
    if ~exist('DataParPath','var') || isempty(DataParPath)
        DataParPath = [];
    end
    if ~exist('ProcessData','var') || isempty(ProcessData)
        ProcessData = [];
    end
    if ~exist('SkipPause','var') || isempty(SkipPause)
        SkipPause = false; % by default we don't skip the pause question below
    end
%     if ~exist('iModules','var') || isempty(iModules)
%         iModules = [1 2 3 4 5]; % structural func dwi ASLS
%     end    

    if regexp(fileparts(pwd), 'CustomScripts')
        cd ..; cd ..;
    end
   
    %% WE DONT RUN NEITHER STRUCTURAL/ASL OR FUNCTIONAL MODULES
   
    % WE ONLY INITIALIZE EXPLOREASL
    x = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers); 
    addpath(genpath(fullfile(x.MyPath,'CustomScripts', 'EPAD'))); % add paths recursively, to be sure it detects the EPAD-specific functions
    
% % %     % 2    xASL_module_func (==fMRI)
% % %     1     Motion correction func
% % %     2     Registration func sessions to T1w
% % %     3     Reslice func data to high resolution standard space
% % %     4     Visualization
% % %     5     QC
% % 
% %     x.nSessions = 1;
% %     x.SESSIONS = {'func'}; % hack ExploreASL to believe only 1 functional session
% %     [~, x] = xASL_Iteration(x,'xASL_module_func');

    %% 3    xASL_module_dwi (==dwi)
    % 1     Motion correction func
    % 2     Registration func sessions to T1w
    % 3     Reslice func data to high resolution standard space
    % 4     Visualization
    % 5     QC

    x.nSessions = 1;
    x.SESSIONS = {'dwi'}; % hack ExploreASL to believe only 1 functional session
    [~, x] = xASL_Iteration(x,'xASL_module_dwi');


end