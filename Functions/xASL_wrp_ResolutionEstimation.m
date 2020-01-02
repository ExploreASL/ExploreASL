function x = xASL_wrp_ResolutionEstimation(x)
if ~isfield(x,'ResolutionEstimation')
    x.ResolutionEstimation    = 0; % default: use default resolutions based on the acquisition resolution
end

%% If ResolutionEstimation is requested, run this here
if x.ResolutionEstimation

    tIM = xASL_io_ReadNifti(x.P.Path_PWI);
    NativeRes = tIM.hdr.pixdim(2:4);

    % Here we choose a range of GM-WM CBF ratios to search along
    % but they can be surprising, so lets give each sequence the same range

    if      strcmp(x.Sequence,'3D_spiral')
            relPSF = [2.6667 2.6667 2.3750];
            % This assumes the inplane interpolation of GE from [4 4 4] to [2 2 4]
    elseif  strcmp(x.Sequence,'2D_EPI')
            relPSF = [1.2167 1.2167 1.0857];
    elseif  strcmp(x.Sequence,'3D_GRASE')
            relPSF = [1.9417 1.9417 1.7304]; % average of 2D EPI & 3D spiral
    end

    x.S.ExpectedFWHM_res = NativeRes.*relPSF;


    x = xASL_im_ResolutionEstim(x); % estimate effective spatial resolution

    % Store results
    xASL_adm_CreateDir( fullfile(x.D.PopDir, 'ResolutionEstimation') );
    CSVfile = fullfile(x.D.PopDir, 'ResolutionEstimation',[x.P.SubjectID '.mat']);
    save(CSVfile,'S');

else
    %% 3 Determine sequence smoothness, using predefined calculations
    x.S.optimFWHM_Res_mm = xASL_init_DefaultEffectiveResolution(x.P.Path_ASL4D, x);

end

x.S.optimFWHM_mm = (abs((x.S.optimFWHM_Res_mm./1.5).^2-1).^0.5)*1.5; % this removes 1 voxel (1.5 mm in MNI) from the smoothing PSF
end
