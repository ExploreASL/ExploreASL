function xASL_wrp_GetROIstatistics(x)
%xASL_wrp_GetROIstatistics Compute statistics for each ROI
%
% FORMAT: xASL_wrp_GetROIstatistics(x)
% 
% INPUT:
%   x                            - struct containing statistical pipeline environment parameters (REQUIRED)
%   x.WBmask                     - WholeBrain mask used to convert image to column &
%                                  vice versa (ExploreASL compression method) (REQUIRED)
%   x.S.InputDataStr             - prefix of data files (i.e. before SubjectSession
%                                  name), to specify which data to compute ROI
%                                  statistics for (e.g. 'qCBF', 'M0' 'mrc1T1') (REQUIRED)
%   x.S.InputAtlasPath           - path to NIfTI file containing atlas to load (REQUIRED)
%   x.S.output_ID                - name of data that is analyzed, to be
%                                  used in the created TSV summary file
%                                  (e.g. CBF, pGM, mean_control) (REQUIRED)
%   x.S.ASL                      - true for ASL data, which will take into account its limited resolution
%                                  & regions of poor SNR: apply PVC, PVC expansion, avoid too
%                                  small ROIs, add vascular & susceptibility artifacts masks (OPTIONAL, DEFAULT=true)
%   x.S.IsVolume                 - true for running volumetric analyses (OPTIONAL, DEFAULT=false)
%   x.D.PopDir                   - path to population folder where standard
%                                  space data should be stored (OPTIONAL, DEFAULT='/analysis/Population')
%   x.S.SubjectWiseVisualization - true to visualize effect of PVC expansion & ROI
%                                  creation. Takes long (OPTIONAL, DEFAULT=off)
%   x.LabEffNorm                 - true to normalize labeling efficiency
%                                  over vascular territories (OPTIONAL, DEFAULT=false)
%   x.KISS                       - true to keep statistics short & simple (OPTIONAL, DEFAULT=true)
%   x.S.GetStats                 - true to determine whether to run extra statistics 
%                                  (permuting over sets/covariates) (OPTIONAL, DEFAULT=false)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This wrapper organizes the computation of statistics for different ROIs
%              in a [1.5 1.5 1.5] mm MNI space:
%              1) Load the atlas: xASL_stat_AtlasForStats
%              2) Organize TSV output name: using x.S.output_ID
%              3) Obtain the ROI statistics: xASL_stat_GetROIstatistics
%              4) Print statistics in TSV files: xASL_stat_PrintStats
%              5) Runs statistics permutation (if requested):
%                using xASL_wrp_PermuteOverSets to permute over
%                sets/covariates & xASL_stat_PrintBasicStats to print
%                several statistics. This last part is historical and not
%                used much anymore.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_wrp_GetROIstatistics(x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% ------------------------------------------------------------------------------------------------------------
%% Admin

if ~isfield(x.S,'IsASL') || isempty(x.S.IsASL)
    x.S.IsASL = true;
end
if ~isfield(x.S,'IsVolume') || isempty(x.S.IsVolume)
    x.S.IsVolume = false; % default
end

if isfield(x.D,'PopDir') || ~isempty(x.D.PopDir)
    x.S.InputDataDir = x.D.PopDir; 
end

x.S.CheckMasksDir = fullfile(x.S.StatsDir,'CheckMasksVisually');

x.S.Input2Check = {'InputDataStr'           'InputAtlasPath'              }; % 'InputDataDir'
x.S.InputFormat = {'string'                 'file'                        }; % 'dir'
x.S.InputNaming = {'name of the input data' 'atlas file (.nii or .nii.gz)'}; % 'input data'

[x] = xASL_adm_uiGetInput(x); % Request missing input fields
[x.S] = xASL_adm_uiGetInput(x.S); % Request missing input fields

if ~isfield(x.S,'SubjectWiseVisualization') || isempty(x.S.SubjectWiseVisualization)
    x.S.SubjectWiseVisualization = false; % this takes lots of computation time, off by default
end    
if iscell(x.S.InputDataStr)
    x.S.InputDataStr = x.S.InputDataStr{1}; 
end

if ~isfield(x,'LabEffNorm') || isempty(x.LabEffNorm)
    x.LabEffNorm = false;
elseif x.LabEffNorm ==true
    x.S.output_ID = [x.S.output_ID '_LabEffNorm'];
end

%% ------------------------------------------------------------------------------------------------------------
%% 1) Load the atlas
x = xASL_stat_AtlasForStats(x); % check all atlases that are requested (see InputAtlasPath above)
% check them, check their ROI names, and make them ready for use below



%% ------------------------------------------------------------------------------------------------------------
%% 2) Organize TSV output name
[~, Ffile] = xASL_fileparts(x.S.InputAtlasPath);
if ~isfield(x.S,'output_ID'); x.S.output_ID  = ''; 
elseif isempty(x.S.output_ID)
    % leave like this
else
    x.S.output_ID = [x.S.output_ID '_']; % add underscore if not empty
end

x.S.output_ID = [x.S.output_ID x.S.InputDataStr '_' Ffile '_n=' num2str(x.nSubjects) '_' date];


%% ------------------------------------------------------------------------------------------------------------
%% 3) Obtain the ROI statistics
x = xASL_stat_GetROIstatistics(x); % calculate median CBF/spatial CoV etc per ROI


%% ------------------------------------------------------------------------------------------------------------
%% 4)  Print statistics in TSV files
if ~isempty(strfind(x.S.InputDataStr,'CBF'))
    x.S.unit = 'mL/100g/min';
elseif ~isempty(strfind(x.S.InputDataStr,'T1'))
    x.S.unit = 'mL';
else
    x.S.unit = 'au';
end

Statistics = {'median' 'mean' 'CoV' 'CoV' 'sum'};
for iStat=1:length(Statistics)
    for iPVC=1:3
        StatField = ['DAT_' Statistics{iStat} '_PVC' num2str(iPVC-1)];
        if isfield(x.S,StatField)
            x.S.DAT = x.S.(StatField);
            FileName = [Statistics{iStat} '_' x.S.output_ID '_PVC' num2str(iPVC-1) '.tsv'];
            x.S.SaveFile = fullfile(x.S.StatsDir, FileName);
            xASL_stat_PrintStats(x); % print individual overview
        end
    end
end



%% -----------------------------------------------------------------------------------------------
%% 5) Runs statistics permutation (if requested)
if ~isfield(x.S,'GetStats') || isempty(x.S.GetStats)
    x.S.GetStats = 0; 
end % This parameter determines whether we run stats

if ~isfield(x,'KISS') || isempty(x.KISS) % keep it short & simply
    x.KISS  = 1;
end
x.S.function2call = @xASL_stat_PrintBasicStats; % if we do statistics, then create a CSV overview with all ...
% types of statistical tests (requires that groups are defined in x.S.SetsID, x.S.SetsNames, x.S.SetsOptions etc)


% if  x.S.GetStats
%     xASL_wrp_PermuteOverSets(x);
% end

    
end
    
%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
%     x.S.HistogramDir      = fullfile( x.HistogramDir,S.output_ID);
%     xASL_adm_CreateDir(x.S.HistogramDir);
%     x.S.oriDIR            = x.S.HistogramDir;
%     x.S.function2call     = @xASL_stat_CreateHistograms;
%     fprintf('%s\n',['Creating ' x.S.output_ID ' histograms']);
%     xASL_wrp_PermuteOverSets(x.S, xASL);
% 
%     if isfield(x.S,'SetsID')
%         x.S.function2call     = @xASL_stat_SpaghettiPlot;
%         fprintf('%s\n',['Creating ' x.S.output_ID ' spaghetti plots']);
%         xASL_wrp_PermuteSetsPer1_2SampleTests( x );
%     end