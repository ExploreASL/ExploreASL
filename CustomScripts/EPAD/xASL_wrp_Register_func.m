function xASL_wrp_Register_func(x)
% xASL_wrp_Register_func (part of ExploreASL)
% Same as xASL_wrp_Register_fMRI, but all ASL-specific stuff removed
% Kept very simple for now, only rigid-body registration

%% ----------------------------------------------------------------------------------------
%% Administration: manage sequences


%% 1    First create ASL-standard space flow field
%       If no T1 flow field exists, create an identity flowfield.
%       We now assume the structural module hasn't run, and we simply want to run the ASL module to quickly check how the images look like
%       So we only run the automatic Center of Mass ACPC alignment
if ~xASL_exist(x.P.Path_y_T1,'file') && ~x.Quality
    IDmatrixPath = fullfile(x.D.MapsDir,'Identity_Deformation_y_T1.nii');
    xASL_Copy(IDmatrixPath, x.P.Path_y_ASL,1);
    xASL_im_CenterOfMass(x.P.Path_func_bold);

elseif ~xASL_exist(x.P.Path_y_T1,'file')
    error('Structural module did not run correctly yet');
else

    %% We assume 2D EPI for fMRI:
    x.Sequence      = '2D_EPI';
    x.readout_dim   = '2D';

    % when we can use a T1w deformation field, smooth it to the ASL resolution
    % Smooth T1w deformation fields to the ASL resolution
    % First estimate ASL resolution
    EstimatedResolution     = xASL_init_DefaultEffectiveResolution( x.P.Path_func_bold, x );
    sKernel                 = (EstimatedResolution.^2 - [1.5 1.5 1.5].^2).^0.5; % assuming the flow fields are in 1.5x1.5x1.5 mm

    % First fill NaNs, to prevent smoothing artifacts
    tIM                     = xASL_io_Nifti2Im(x.P.Path_y_T1);
    for iM=1:size(tIM,5)
        tIM(:,:,:,:,iM)     = xASL_im_ndnanfilter(tIM(:,:,:,:,iM),'gauss',[16 16 16],2);
    end
    xASL_io_SaveNifti(x.P.Path_y_T1,x.P.Path_y_T1,tIM,[],0);

    % Perform the smoothing
    xASL_spm_smooth(x.P.Path_y_T1, sKernel,x.P.Path_y_ASL);
    % Solve edges
    FieldT1                             = xASL_io_Nifti2Im(x.P.Path_y_T1);
    FieldASL                            = xASL_io_Nifti2Im(x.P.Path_y_ASL);
    FieldASL([1:5 end-5:end],:,:,:,:)   = FieldT1([1:5 end-5:end],:,:,:,:); % x
    FieldASL(:,[1:5 end-5:end],:,:,:)   = FieldT1(:,[1:5 end-5:end],:,:,:); % y
    FieldASL(:,:,[1:5 end-5:end],:,:)   = FieldT1(:,:,[1:5 end-5:end],:,:); % z
    xASL_io_SaveNifti( x.P.Path_y_ASL, x.P.Path_y_ASL, FieldASL,[],0);
end

%% 0    Clip image, optimizes image contrast for registration
xASL_delete(x.P.Path_mean_control);
xASL_io_PairwiseSubtraction(x.P.Path_func_bold, x.P.Path_mean_PWI_Clipped,0,0,x); % create PWI & mean_control
if ~exist(x.P.Path_mean_control, 'file')
    xASL_Move(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_control);
else
    xASL_delete(x.P.Path_mean_PWI_Clipped);
end


%% ----------------------------------------------------------------------------------------
%% 1    Registration fMRI -> anat
% Here a temporary CBF image is created, which will be used for
% registration to T1 GM prob map & can be used for registration to previous ASL sessions.
% PWI-based registration can be preferable over EPI-based registration if
% background suppression and/or other 3D readout techniques are used.

OtherList = '';
OtherList{1,1} = x.P.Path_func_bold;
OtherList{end+1,1} = x.P.Path_func_NormPE;
OtherList{end+1,1} = x.P.Path_func_RevPE;

xASL_spm_coreg(x.P.Path_c1T1, x.P.Path_mean_control, OtherList, x );

%% ----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% if this doesn't work, we can use a 2D EPI template instead of the T1w
%% ----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fprintf('\n%s\n','--------------------------------------------------------------------');



end
