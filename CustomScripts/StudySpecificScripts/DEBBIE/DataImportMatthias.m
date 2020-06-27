%% Import FME Adaptive ASL DEBBIE grant
AnalysisDir = '/Users/henk/ExploreASL/ASL/BBB_ASL/analysis_AdaptiveASL';

Subjects = xASL_adm_GetFileList(AnalysisDir, '^fme_19_\d{4}$', 'List', [0 Inf], 1);

for iSubject=1:length(Subjects)
    xASL_TrackProgress(iSubject,length(Subjects));
    NiftiDir = fullfile(AnalysisDir, Subjects{iSubject}, 'Nifties');
    PathOff = xASL_adm_GetFileList(NiftiDir, '^ASL_adaptive.*off\.nii$', 'FPList');
    PathOn = xASL_adm_GetFileList(NiftiDir, '^ASL_adaptive.*on(|_rep1)\.nii$', 'FPList');
    PathM0 = xASL_adm_GetFileList(NiftiDir, '^M0\.nii$', 'FPList');
    
    % for off
    if isempty(PathOff) || length(PathOff)>1
        fprintf('missing off\n');
        warning(num2str(iSubject));
    else
        % Load ASL
        IM = xASL_io_Nifti2Im(PathOff{1});
        
        for iTimePoint=1:7
            TPDir = fullfile(fileparts(NiftiDir), ['ASL_' num2str(iTimePoint)]);
            xASL_adm_CreateDir(TPDir);
            
            % Create ASL
            PathASL = fullfile(TPDir, 'ASL4D.nii');
            xASL_io_SaveNifti(PathOff{1}, PathASL, IM(:,:,:,iTimePoint), [], 0);
        end
    end
    
    % same for on
    if isempty(PathOn) || length(PathOn)>1
        fprintf('missing on\n');
        warning(num2str(iSubject));
    else
        % Load ASL
        IM = xASL_io_Nifti2Im(PathOn{1});
        
        for iTimePoint=1:7
            TPDir = fullfile(fileparts(NiftiDir), ['ASL_' num2str(iTimePoint+7)]);
            xASL_adm_CreateDir(TPDir);
            
            % Create ASL
            PathASL = fullfile(TPDir, 'ASL4D.nii');
            xASL_io_SaveNifti(PathOn{1}, PathASL, IM(:,:,:,iTimePoint), [], 0);            
        end        
    end
    
    % for M0
    if isempty(PathM0) || length(PathM0)>1
        fprintf('missing M0\n');
        warning(num2str(iSubject));
    else
        for iTimePoint=1:14
            TPDir = fullfile(fileparts(NiftiDir), ['ASL_' num2str(iTimePoint)]);
            if exist(TPDir, 'dir')
                PathM0dest = fullfile(TPDir, 'M0.nii');
                xASL_Copy(PathM0{1}, PathM0dest);
            end
        end
    end
end
    
    