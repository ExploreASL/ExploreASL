function xASL_wrp_PermuteSets1( x )
%xASL_wrp_PermuteSets1 
% INPUT
% x = from ExploreASL
% masks   = physically loaded masks for each subject
% ASL     = the data to be analyzed. Could be ASL, or e.g. SD or SNR masks
%
%
% By HJMM Mutsaerts, ExploreASL 2016
%
% The idea behind this permutation function is that you want for some functions
% all possible permutations. E.g. for summarizing data in means and comparing sets using t-tests,
% it will not necessarily differ between 1 (sessions) or 2 sample (cohorts/scanners) sets
% e.g. with the use of histogram_creations or printing_ROI_summaries

% This function prints an overview of all data first, the permutation
% happens in the bottom




%% -----------------------------------------------------------------------------------------------
%% Admin

fprintf('%s\n',['Printing tsv-files with ' x.S.output_ID ' statistics']);

if  isfield(x.S,'SetsName')
    if  length(x.S.SetsName)~= size(x.S.SetsID,2)
        error('Number of sets names differs from number of values (SetsID)');
    end

%    if  ~strcmp(x.S.output_ID,'volume') % with volume, if there are multiple sessions the x.S.SetsID is restructured into a single session per subject,
%        % assuming a single T1 per subject
%        x.S.SetsName              = x.SetsName;
%        x.S.SetsID                = x.SetsID;    
%    end
end

if ~isfield(x.S,'GetStats')
    x.S.GetStats  = 0;
end





%% -----------------------------------------------------------------------------------------------
%% Here the main output directory is selected
%  This should initially get the output_ID (e.g. CBF, volume, motion etc.)
%  but not if it already is some levels lower

if  ~isempty(strfind(x.S.StatsDir, 'permute')) || ~isempty(strfind(x.S.StatsDir, 'All scans'))
    x.S.oriDIR                = fullfile( x.S.StatsDir);
else
    x.S.oriDIR                = fullfile( x.S.StatsDir, x.S.output_ID);
end    



xASL_adm_CreateDir(x.S.oriDIR);
x.S.SaveFile                  = fullfile( x.S.oriDIR, [x.S.output_ID '.tsv']);

[x] = xASL_stat_PrintStats(x); % print individual overview





%% -----------------------------------------------------------------------------------------------
%% This runs the permutation of statistical comparisons, if requested

if  x.S.GetStats
    xASL_wrp_PermuteOverSets( x );
end
    
 
 end
