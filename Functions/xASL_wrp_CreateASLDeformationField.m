function xASL_wrp_CreateASLDeformationField(x, bOverwrite, EstimatedResolution)
%xASL_wrp_CreateASLDeformationField When we can use a T1w deformation field, smooth it to the ASL resolution
% Smooth T1w deformation fields to the ASL resolution
% First estimate ASL resolution

if nargin<2 || isempty(bOverwrite)
    bOverwrite = false; % By default, this function will be skipped if the ASL deformation field already exists
end

if nargin<3 || isempty(EstimatedResolution)
	EstimatedResolution = xASL_init_DefaultEffectiveResolution(x.P.Path_ASL4D, x);
end

if ~xASL_exist(x.P.Path_y_ASL,'file') || bOverwrite

    if ~xASL_exist(x.P.Path_y_T1,'file')
        warning([x.P.Path_y_T1 ' didnt exist, skipping...']);
        return;
    end

    TemplateResolution = [1.5 1.5 1.5];
    sKernel = (EstimatedResolution.^2 - TemplateResolution.^2).^0.5; % assuming the flow fields are in 1.5x1.5x1.5 mm
    sKernel(EstimatedResolution<TemplateResolution) = 0;

    xASL_im_FixEdgesFlowfield(x.P.Path_y_T1); % First fill NaNs, to prevent smoothing artifacts

    % Perform the smoothing
    % xASL_spm_smooth(x.P.Path_y_T1, sKernel,x.P.Path_y_ASL);
    xASL_im_PreSmooth(x.P.Path_ASL4D, x.P.Path_y_T1, x.P.Path_y_ASL); % we need to add the effective resolution here still!
    % sKernel, as calculated above, can be used for this. But the major
    % rotations need to be taken into account, between the effective
    % resolution as specified, and the one in the different NIfTIs
    % (e.g. the ASL & T1w are usually acquired transversal & sagittally,
    % respectively

    % Solve edges
    FieldT1 = xASL_io_Nifti2Im(x.P.Path_y_T1);
    FieldASL = xASL_io_Nifti2Im(x.P.Path_y_ASL);
    FieldASL([1:5 end-5:end],:,:,:,:) = FieldT1([1:5 end-5:end],:,:,:,:); % x
    FieldASL(:,[1:5 end-5:end],:,:,:) = FieldT1(:,[1:5 end-5:end],:,:,:); % y
    FieldASL(:,:,[1:5 end-5:end],:,:) = FieldT1(:,:,[1:5 end-5:end],:,:); % z
    xASL_io_SaveNifti(x.P.Path_y_ASL, x.P.Path_y_ASL, FieldASL, [], 0);

    xASL_im_FixEdgesFlowfield(x.P.Path_y_ASL); % Again fill NaNs, if needed
end



end
