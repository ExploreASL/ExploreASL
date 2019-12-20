function [ Corr_M0 ] = xASL_im_ProcessM0Conventional( ImIn, x )
%xASL_im_M0ConventionalProcessing This function uses the conventional M0 masking,
% and only a little smoothing, following what Philips uses for its 3D
% GRASE. Advantages of the newer M0 processing in ExploreASL are the lack
% of use of M0 threshold-based masking, the removal of high CSF values and
% higher SNR for ASL division.



    %% Mask data below threshold
    LowCutOff                       = 0.2;
    HiCutOff                        = 0.6;

    Data_array                      = sort(ImIn(:));
    % 
    Data_array(Data_array<(LowCutOff.*max(Data_array)))    = [];
    Data_array(Data_array>(HiCutOff.*max(Data_array)))    = [];
    % 
    % [N X]=hist(MaskM0)
    % figure(2);plot(X,N)

    Diff_array                      = Data_array(2:end) - Data_array(1:end-1);
    placemax                        = find(Diff_array == max(Diff_array));
    Threshold                       = (Data_array(placemax)+Data_array(placemax+1))./2;

    MaskM0                          = ImIn;
    MaskM0(MaskM0<Threshold(1))     = NaN;


    %% Mask-constrained smoothing

    if ~isfield(x,'SmoothingFactor')
        x.SmoothingFactor        = 8; % 8 mm = consensus paper default, 6.5 mm = Philips default
    end
    
    N                      			= x.SmoothingFactor./x.VoxelSize;
	if numel(N)==1
		N = repmat(N,[1 3]);
	end

    Corr_M0                         = xASL_im_ndnanfilter(MaskM0,'gauss',double(N),0);
end