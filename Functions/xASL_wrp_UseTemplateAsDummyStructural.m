function xASL_wrp_UseTemplateAsDummyStructural(x, Template)
% bUseTemplateAsDummyStructural for copying Structural templates and
% running registration of ASL to templates
% INPUT:
%   x                   - ExploreASL structure containing fields with global information about the pipeline environment
%                         and settings (e.g. x.settings.Quality), useful when you want this script to copy the options of an ExploreASL pipeline run 
%   Template            - Template name string used for copying to subject directory and registration of ASL to this structural template and derivatives (default = MNI_Structural)
%
% OUTPUT: n/a


%% 1. admin
% paths -> may contain all Template paths currently supported by ExploreASL
TemplateList = {'TotalGM', 'TotalWM', 'DeepWM', 'WholeBrain', 'MNI_Structural', 'CortVascTerritoriesTatu', 'TatuICA_PCA', 'LabelingTerritories', 'ATTbasedFlowTerritories'};
if nargin<2 || isempty('Template') % use default MNI template
    Template = 'MNI_Structural';
end

if strcmp(Template, 'MNI_Structural')
    TemplateFolder = x.D.MapsSPMmodifiedDir;
    rc1T1 = 'rc1T1.nii';
    rc2T1 = 'rc2T1.nii';
    rc3T1 = 'rc3T1.nii';
    rT1 = 'rT1.nii';    
    
    TemplatePath = fullfile(TemplateFolder, 'rT1.nii');
elseif strcmp(Template, 'QASPER') % check if template is QASPER
    % implement proper location later
    TemplateFolder = fullfile(x.D.TemplateDir, 'QASPER');
    rc1T1 = 'QASPER_pInlet.nii';
    rc2T1 = 'QASPER_pOutlet.nii';
    rc3T1 = 'QASPER_pPorousMedium.nii';
    rT1 = 'QASPER_Template.nii';
    
    TemplatePath = fullfile(TemplateFolder, rT1);

elseif ismember(TemplateList, Template) % check if supplied template is supported by ExploreASL
    fprintf('Currently, templates other than MNI_Structural (default) or QASPER are not implemented yet')
end


%% 2.1 copy and register structural template files to subject directory and subject space  
fprintf('Missing structural scans, using ASL registration only instead, copying structural template as dummy files\n');

IDmatrixPath = fullfile(x.D.MapsSPMmodifiedDir, 'Identity_Deformation_y_T1.nii');

if xASL_exist(TemplatePath,'file')
    % Copy dummy transformation field
    xASL_Copy(IDmatrixPath, x.P.Path_y_T1, true);

    % Use the GM, WM and T1 from the template to create dummy structural files
    % Create dummy structural derivatives in standard space
    xASL_Copy(fullfile(TemplateFolder, rc1T1), x.P.Pop_Path_rc1T1);
    xASL_Copy(fullfile(TemplateFolder, rc2T1), x.P.Pop_Path_rc2T1);
    xASL_Copy(fullfile(TemplateFolder, rc3T1), x.P.Pop_Path_rc3T1);
    xASL_Copy(fullfile(TemplateFolder, rT1), x.P.Pop_Path_rT1);

    % Create dummy structural derivatives in native space
    xASL_spm_deformations(x, {x.P.Pop_Path_rc1T1, x.P.Pop_Path_rc2T1, x.P.Pop_Path_rc3T1, x.P.Pop_Path_rT1}, {x.P.Path_c1T1, x.P.Path_c2T1, x.P.Path_c3T1, x.P.Path_T1});

else 
    fprintf('Template cannot be found, check if this template is available within ExploreASL');
end
%% 2.2 create dummy and lock files
        % Dummy files
        catVolFile = fullfile(x.D.TissueVolumeDir,['cat_' x.P.STRUCT '_' x.P.SubjectID '.mat']);
        MatFile   = fullfile(x.dir.SUBJECTDIR, [x.P.STRUCT '_seg8.mat']);
        dummyVar = [];
        save(catVolFile,'dummyVar');
        save(MatFile,'dummyVar');

        SaveFile = fullfile(x.D.TissueVolumeDir,['TissueVolume_' x.P.SubjectID '.csv']);
        FileID = fopen(SaveFile,'wt');
        fprintf(FileID,'%s', '0, 0, 0');
        fclose(FileID);
        xASL_bids_csv2tsvReadWrite(SaveFile, true);

        % To lock in the structural part
        % Save the ASL lock and unlock
        jj = strfind(x.dir.LockDir,'xASL_module_ASL');
        jj = jj(1);
        oldRoot = x.mutex.Root;
        oldID = x.mutex.ID;
        newRoot = fullfile(x.dir.LockDir(1:(jj-1)),'xASL_module_Structural',x.dir.LockDir((jj+16):end));
        x.mutex.Unlock();

        % Look the structural part
        x.mutex.Root = newRoot;
        x.mutex.Lock('xASL_module_Structural');

        % Add the correct lock-files
        x.mutex.AddState('010_LinearReg_T1w2MNI');
        x.mutex.AddState('060_Segment_T1w');
        x.mutex.AddState('080_Resample2StandardSpace');
        x.mutex.AddState('090_GetVolumetrics');
        x.mutex.AddState('100_VisualQC_Structural');
        x.mutex.AddState('110_DoWADQCDC');

        % Unlock the structural and lock again the ASL part
        x.mutex.Unlock();
        x.mutex.Root = oldRoot;
        x.mutex.Lock(oldID);


end 