function [result, x] = xASL_module_Population(x)
%xASL_module_Population ExploreASL module for population-based/group-based processing
%
% FORMAT: [result, x] = xASL_module_Population(x)
%
% INPUT:
%   x       - x structure containing all input parameters (REQUIRED)
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL module processes all available images on the
% group level. It assumes that all images were adequately processed in the
% previous modules. It will perform the following group-wise processing and
% checks:
% 
% 010_CreatePopulationTemplates - Create population average images, to compare scanners, cohorts etc without physiological variance
% 020_CreateAnalysisMask        - Generate a group-level mask by combining individuals masks, for ROI-based analysis & VBA
% 030_CreateBiasfield           - When there are multiple scanners, create scanner-specific biasfields (uses Site.mat for this)
% 040_GetDICOMStatistics        - Create TSV file with overview of DICOM parameters
% 050_GetVolumeStatistics       - Create TSV file with overview of volumetric parameters
% 060_GetMotionStatistics       - Create TSV file with overview of motion parameters
% 070_GetROIstatistics          - Create TSV file with overview of regional values (e.g. qCBF, mean control, pGM etc)
% 080_SortBySpatialCoV          - Sort ASL_Check QC images by their spatial CoV in quality bins
% 090_DeleteAndZip              - Delete temporary files and gzip all NIfTIs
%
% EXAMPLE: [~, x] = xASL_module_Population(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2019 ExploreASL



%% ------------------------------------------------------------------------------------------------------------
%% Admin

x = xASL_init_GenericMutexModules(x, 'QC'); % starts mutex locking process to ensure that everything will run only once
x = xASL_init_FileSystem(x);
xASL_adm_CreateDir(x.S.StatsDir);

% Define ASL sequence (quick&dirty, if ASL information exists)
if ~isfield(x,'Sequence') && isfield(x,'readout_dim')
    if strcmp(x.readout_dim,'2D')
       x.Sequence = '2D_EPI'; % assume that 2D is 2D EPI, irrespective of vendor
    elseif strcmp(x.readout_dim,'3D') && ( ~isempty(regexp(x.Vendor,'Philips')) || ~isempty(regexp(x.Vendor,'Siemens')) )
           x.Sequence = '3D_GRASE'; % assume that 3D Philips or Siemens is 3D GRASE
    elseif strcmp(x.readout_dim,'3D') && ~isempty(regexp(x.Vendor,'GE'))
           x.Sequence = '3D_spiral'; % assume that 3D GE is 3D spiral
    end
end

StateName{1} = '010_CreatePopulationTemplates';
StateName{2} = '020_CreateAnalysisMask';
StateName{3} = '030_CreateBiasfield';
StateName{4} = '040_GetDICOMStatistics';
StateName{5} = '050_GetVolumeStatistics';
StateName{6} = '060_GetMotionStatistics';
StateName{7} = '070_GetROIstatistics';
StateName{8} = '080_SortBySpatialCoV';
StateName{9} = '090_DeleteAndZip';



%% ------------------------------------------------------------------------------------------------------------
%% 1    Create template images
if ~x.mutex.HasState(StateName{1})
    
    % Get initial data
    IndexCohort = find(cellfun(@(y) strcmp(y,'Cohort'), x.S.SetsName));
    IndexTP = find(cellfun(@(y) strcmp(y,'LongitudinalTimePoint'), x.S.SetsName));
    CohortIs{1} = x.S.SetsID(:,IndexCohort)==1 & x.S.SetsID(:,IndexTP)==1; % HIV 1
    CohortIs{2} = x.S.SetsID(:,IndexCohort)==1 & x.S.SetsID(:,IndexTP)==2; % HIV 2
    CohortIs{3} = x.S.SetsID(:,IndexCohort)==2 & x.S.SetsID(:,IndexTP)==1; % HC 1
    CohortIs{4} = x.S.SetsID(:,IndexCohort)==2 & x.S.SetsID(:,IndexTP)==2; % HC 2
    NameIs = {'HIV1' 'HIV2' 'HC1' 'HC2'};
    
    % Backup data
    SubjectsBackup = x.SUBJECTS;
    SetsIDBackup = x.S.SetsID;

    PathFLAIR = fullfile(x.D.TemplatesStudyDir,'FLAIR_bs-mean.nii');
    PathWMH = fullfile(x.D.TemplatesStudyDir,'WMH_SEGM_bs-sum.nii');
    
    for iCohort=1:4
        x.SUBJECTS = SubjectsBackup(CohortIs{iCohort});
        x.nSubjects = length(x.SUBJECTS);
        x.S.SetsID = SetsIDBackup(CohortIs{iCohort},:);
        xASL_wrp_CreatePopulationTemplates(x, false, false, {{['r' x.P.FLAIR]} {'FLAIR'} 0}, false, false, {{@xASL_stat_MeanNan} {'mean'}});
        xASL_wrp_CreatePopulationTemplates(x, false, false, {{['r' x.P.WMH_SEGM]} {'WMH_SEGM'} 0}, false, false, {{@xASL_stat_SumNan} {'sum'}});
        Path_TemplateFLAIR = fullfile(x.D.TemplatesStudyDir,['FLAIR_' NameIs{iCohort} '.nii']);
        Path_TemplateWMH = fullfile(x.D.TemplatesStudyDir,['WMH_' NameIs{iCohort} '.nii']);
        xASL_delete(Path_TemplateFLAIR);
        xASL_delete(Path_TemplateWMH);
        xASL_Move(PathFLAIR, Path_TemplateFLAIR);
        xASL_Move(PathWMH, Path_TemplateWMH);
    end
        

%     % 5) Get average FLAIR
%     x.SUBJECTS = SubjectsBackup;
%     x.S.SetsID = SetsIDBackup;
%     xASL_wrp_CreatePopulationTemplates(x, false, false, {{['r' x.P.FLAIR]} {'FLAIR'} 0}, false, false);

    x.mutex.AddState(StateName{1});
    fprintf('%s\n',[StateName{1} ' was performed']);
else
    fprintf('%s\n',[StateName{1} ' has already been performed, skipping...']);
end

    
    
%% -----------------------------------------------------------------------------
%% 9    Reduce data size
if ~x.mutex.HasState(StateName{9})

    % Zip temporary files
    % This way, temporary files that take up a lot of data, will take up have
    % the disc space. These files may be re-used in the event of a
    % re-processing, which is why we zip them instead of deleting them.
    xASL_adm_GzipAllFiles(x.D.ROOT);
    x.mutex.AddState(StateName{9});
    fprintf('%s\n',[StateName{9} ' was performed']);
else
    fprintf('%s\n',[StateName{9} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
result = true;

end
