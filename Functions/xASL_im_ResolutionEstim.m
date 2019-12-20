function x = xASL_im_ResolutionEstim(x)
%ResolutionEstim Estimation spatial effective resolution
% by smoothing the pGM+pWM images, with a given GM-WM ratio,
% and investing the Sum of Squares with the CBF image


    % Example input for resolution estimation on group level
    % x.S.rpGMname  = 'C:\Backup\ASL\ExploreASL_Example\Example_SingleSubject_Unproc\dartel\Templates\Template_mean_pGM_1.nii';
    % x.S.rpWMname  = 'C:\Backup\ASL\ExploreASL_Example\Example_SingleSubject_Unproc\dartel\Templates\Template_mean_pWM_1.nii';
    % x.S.UpsampASL = 'C:\Backup\ASL\ExploreASL_Example\Example_SingleSubject_Unproc\dartel\Templates\Template_mean_CBF_1.nii';
    % % corrected vascular artifacts is better

    imGM        = xASL_io_Nifti2Im(x.S.rpGMname);
    imWM        = xASL_io_Nifti2Im(x.S.rpWMname);
    imCBF       = xASL_io_Nifti2Im(x.S.UpsampASL);

    if ~isfield(x.S,'UpsampFoV')
        imFoV     = ones(size(imCBF));
    else
        imFoV     = xASL_io_Nifti2Im(x.S.UpsampFoV);
    end

    if ~isfield(x.S,'ExpectedFWHM_res')
        x.S.ExpectedFWHM_res  = ([5 5 9]+[4 4 6]+[3 3 7])./3; % average expected resolution for all sequences
	end

	if ~isfield(x.S,'PSF') % this is the kernel we use to estimate smoothness
		PSFtype       = {'gaussian' 'gaussian' 'gaussian'};
	else
		PSFtype = x.S.PSF;
	end

    %startSigmaVec       = (x.S.ExpectedFWHM_res./2.355)./1.5;

    %% Rescale images
    imCBF               = imCBF./max(imCBF(:));
    imGM                = imGM ./max(imGM(:));
    imWM                = imWM ./max(imWM(:));

    %% Define maximal n iterations
    if ~exist('x', 'var')
        x.Quality     = 1;
    elseif ~isfield(x,'Quality')
        x.Quality     = 1;
	end

	if ~isfield(x.S,'MaxIter')
		if  x.Quality
			MaxIter   = 20;
		else
			MaxIter   = 5;
		end
	else
		MaxIter = x.S.MaxIter;
	end

    %% Clip vascular extremes
    qnt             = sort(imCBF(:),'ascend');
    qnt             = qnt(floor(length(qnt).*0.95));
    imCBF           = imCBF./qnt;
    imCBF(imCBF>1)  = 1;

    imGM(imGM<0) = 0;
    imGM(imGM>1) = 1;
    imWM(imWM<0) = 0;
    imWM(imWM>1) = 1;

    %% Define Mask
    % Constricted to whole brain, and inside the ASL FoV. Otherwise it goes
    %
    structMask                      = (imGM+imWM)>0.1;
    FoVMask                         = imFoV>0;
    imMask                          = structMask.*FoVMask;
    imMask(:,:,[1:44 86:end])      = 0;

    %% Run resolution estimation

    fprintf('%s\n','Estimating ASL resolution...  ');

    [ resFWHM, resSigma,~,~,~] = xASL_im_EstimateResolution( imCBF,imGM,imWM,imMask,PSFtype,MaxIter);

    fprintf('\n');

    x.S.optimSigma_Res_vox          = resSigma;   % the estimated resolution in voxels, in sigma
    x.S.optimFWHM_Res_vox           = resFWHM;    % the estimated resolution in voxels, in FWHM
    x.S.optimSigma_Res_mm           = x.S.optimSigma_Res_vox.*1.5;
    x.S.optimFWHM_Res_mm            = x.S.optimFWHM_Res_vox.*1.5; % assuming 1.5 mm voxel-size in ExploreASL

    % this is what we need to smooth the pGM & pWM with, to get the pseudo-CBF map in correct effective resolution

    fprintf('%s\n',['Spatial effective resolution estimated as FWHM = [' num2str(x.S.optimFWHM_Res_mm(1),3) ' ' num2str(x.S.optimFWHM_Res_mm(2),3) ' ' num2str(x.S.optimFWHM_Res_mm(3),3) '] mm']);


end
