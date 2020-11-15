function extract_perfusion_stats(xASL_path,data_dir,kernel_size,PVC_mask_treshold,quantile_lim)
% This is a script created by Tamas Jozsa to extract statistics from CBF maps processed by ExploreASL
%
% data_dir - string of the folder containing the data
% kernel_size - odd integer, larger than 2: 3, 5, 7, ...
% PVC_mask_treshold - float representing tissue partial volume, so that 0.6 will consider only data points where PV_GM and PV_WM > 0.6
% qunatile_lim - float to extract symmetric quantiles: 0.25 will lead to 25 and 75 percentiles
%
% Example:
% extract_perfusion_stats('./ExploreASL/','./EPAD_stats/',3,0.6,0.25)

clc;

%% 1. define files, folders, and kernel size
% TODO: please specify folders and Kernel size
stat_dir = data_dir;

outfile = ['patient_stats_kernel' num2str(kernel_size) '_trsh' num2str(PVC_mask_treshold,'%2.0e') '_ql' num2str(0.25,'%2.1e') '.csv'];

myKernel = [kernel_size,kernel_size,kernel_size];

%% 2. initialise xASL
CurrentPath = pwd;
cd(xASL_path);
ExploreASL_Master('',0); % initialize ExploreASL, without processing data
cd(CurrentPath); % back to the initial path

%% 3. compute statistical data for patients

load(fullfile(stat_dir,'patient_info.mat'));
ASL_data_folder = 'ASL_1/';

% ID, age, sex, Vol_GM, Vol_WM, mCBF_GM, minCBF_GM, maxCBF_GM, mCBF_WM, minCBF_WM, maxCBF_WM
patient_data = zeros(sum(pass_check),11);

[VGM, VWM, meanCBF_GM_corr, GM05perc, GM95perc, meanCBF_WM_corr, WM05perc, WM95perc] = ...
    deal( 0,0,0,0,0,0,0,0 );

cter = 1;
for ii = 1:n_patients
    if pass_check(ii) && exist(fullfile(stat_dir, patient_folders{ii}, patient_folders{ii}, ASL_data_folder, 'CBF.nii'))>0
        % read age, sex info
        % TODO: please implement reading patients' age and sex -> H: Will do
        age = str2double( Age(ii) );
        if strcmp(Sex{ii},'Male')
            sex = 1;
        elseif strcmp(Sex{ii},'Female')
            sex = 2;
        else
            sex = NaN;
        end

        if isnan( age ) == false && isnan( sex ) == false
            PathCBF = fullfile(stat_dir, patient_folders{ii}, patient_folders{ii}, ASL_data_folder, 'CBF.nii');
            NIfTI = xASL_io_ReadNifti(PathCBF);
            Volume_voxel = prod(NIfTI.hdr.pixdim(2:4));

            % read images
            PVgm = xASL_io_Nifti2Im(fullfile(stat_dir, patient_folders{ii}, patient_folders{ii}, ASL_data_folder, 'PVgm.nii')); % use fullfile (OS independent)
            % note that the xASL_ prefixed scripts all deal transparently with % .nii(.gz)
            PVwm = xASL_io_Nifti2Im(fullfile(stat_dir, patient_folders{ii}, patient_folders{ii}, ASL_data_folder, 'PVwm.nii'));

            CBF = xASL_io_Nifti2Im(PathCBF);
            MaskVascular = xASL_io_Nifti2Im(fullfile(stat_dir, patient_folders{ii}, patient_folders{ii}, ASL_data_folder, 'MaskVascular.nii'));

            % check readings
        %     [nx, ny, nz] = size(MaskVascular);
        %
        %     figure(1); imshow(rot90(squeeze(PVgm(:,round(ny/2),:)),1)) % I can
        %     recommend the DIP (Delft Image Library), with e.g. dip_image()
        %     instead of imshow(), very nice to scroll through slices/volumes etc,
        %     and joking
        %     we do have our own replacements of this image processing toolbox,
        %     e.g. xASL_im_rotate()

            % compute perfusion stats [mL/min/100g]
            [meanCBF_GM_corr, meanCBF_WM_corr] = xASL_stat_ComputeMean(CBF, MaskVascular, [], 2, PVgm, PVwm);

            [CBF_stack_corr,residual] = xASL_im_PVCkernel(CBF, cat(4,PVgm,PVwm), myKernel, 'gauss');
            CBF_GM_corr = CBF_stack_corr(:,:,:,1);
            CBF_WM_corr = CBF_stack_corr(:,:,:,2);

            GM_mask = PVgm>PVC_mask_treshold & MaskVascular==1;
            WM_mask = PVwm>PVC_mask_treshold & MaskVascular==1;

            CBF_GM_array = CBF_GM_corr(GM_mask);
            CBF_WM_array = CBF_WM_corr(WM_mask);

            % compute GM and WM volumes [mL]
            VGM = sum(Volume_voxel*PVgm(PVgm>PVC_mask_treshold)); % if you mask the PVC CBF above with e.g. pGM>0.5, then you want to do the same here for getting the volume right?
            VWM = sum(Volume_voxel*PVwm(PVwm>PVC_mask_treshold));

            if exist('quantile')>0
                GM_stat = num2cell( quantile(CBF_GM_array,[quantile_lim, 1-quantile_lim]) );
                WM_stat = num2cell( quantile(CBF_WM_array,[quantile_lim, 1-quantile_lim]) );
            else
                GM_stat = num2cell( quantile_fct(CBF_GM_array,[quantile_lim, 1-quantile_lim]) );
                WM_stat = num2cell( quantile_fct(CBF_WM_array,[quantile_lim, 1-quantile_lim]) );
            end

            [GM_perc_a,GM_perc_b] = deal( GM_stat{:} );
            [WM_perc_a,WM_perc_b] = deal( WM_stat{:} );

            % create database
            patient_data(cter,:) = [cter, age, sex, VGM, VWM, ...
                                  meanCBF_GM_corr, GM_perc_a, GM_perc_b, ...
                                  meanCBF_WM_corr, WM_perc_a, WM_perc_b];
            cter = cter + 1;
        end
    end
end

%% save database
patient_data = array2table(patient_data,'VariableNames',...
               {'ID','age','sex','Vol_GM','Vol_WM','mCBF_GM','pA_CBF_GM','pB_CBF_GM','mCBF_WM','pA_CBF_WM','pB_CBF_WM'});
writetable(patient_data,outfile);
