function ImOut = xASL_im_M0ErodeSmoothExtrapolate( ImIn, x)
%xASL_im_M0ErodeSmoothExtrapolate Erodes, smooths & extrapolates M0
% for ASL population analyses
% Assumes images are in standard space & that GM & WM probability maps
% are registered

% ExploreASL, HJMM Mutsaerts

% Here, we mask the M0, to remove high CSF signal and low extracranial signal,
% enabling us to smooth the image without inserting wrong signal.





%% ------------------------------------------------------------------------------------------
%% Mask 1) Load segmentations, create structural mask
ExistpGMpWM     = xASL_exist(x.P.Pop_Path_rc1T1) && xASL_exist(x.P.Pop_Path_rc2T1);
if  ExistpGMpWM
    fprintf('%s\n','Masking M0 with structural (pGM+pWM)>0.5 & 70% non-zero sorted intensities');
    GMim                        = xASL_io_Nifti2Im( x.P.Pop_Path_rc1T1);
    WMim                        = xASL_io_Nifti2Im( x.P.Pop_Path_rc2T1);

    GMmask                      = GMim>0.7;
    Mask1                       = (GMim+WMim)>0.5;
else
    fprintf('%s\n','Masking M0 with intensity-based mask only, structural (pGM+pWM) files missing');
end




%% ------------------------------------------------------------------------------------------
%% Mask 2) Create intensity-based mask to remove extracranial signal
SortInt                     = sort(ImIn(:));
SortInt                     = SortInt(~isnan(SortInt));
ThresholdN                  = SortInt(round(0.7*length(SortInt)));

% Remove peak signal as well
ThresholdN2                 = SortInt(round(0.999*length(SortInt)));
% Combine these masks
if  ExistpGMpWM
    Mask2                   = Mask1 & (ImIn>ThresholdN) & (ImIn<ThresholdN2);
else
    Mask2                   = (ImIn>ThresholdN) & (ImIn<ThresholdN2);
    GMmask                  = ImIn>SortInt(round(0.85*length(SortInt))) & ImIn<SortInt(round(0.95*length(SortInt)));
end






%% ------------------------------------------------------------------------------------------
%% Mask 3) Erode the combined masks
fprintf('%s\n','Erode M0 mask with 2-voxel sphere');
Mask3                       = xASL_im_DilateErodeFull(Mask2,'erode',xASL_im_DilateErodeSphere(2));






%% ------------------------------------------------------------------------------------------
%% Mask 4) Determine any odd borders
fprintf('%s\n','Identify & remove non-smooth values at the border of the M0 mask');

% First get median inside the mask
ValuesM0mask                = ImIn(Mask3 & isfinite(ImIn));
MedianN                     = median(ValuesM0mask);
% Fill voxels outside mask with this median
TempIM                      = ImIn;
TempIM(~Mask3)              = MedianN;

% Smooth ImOut, to check the voxels that change most by smoothing
% (especially those with low values)
ImOutSmooth                 = xASL_im_ndnanfilter(TempIM,'gauss',double([5 5 5]),0);
DiffIM                      = abs(ImOutSmooth-TempIM);

DiffIMvalues                = DiffIM(Mask3 & isfinite(DiffIM));

MedianN                     = median(DiffIMvalues);
MadN                        = xASL_stat_MadNan(DiffIMvalues);
HiThresh                    = MedianN+4.5*MadN;
% here we remove too many parts, but this is not bad, as it will be smoothed & extrapolated, anyway.
Mask4                       = Mask3 & ~(DiffIM>HiThresh);



%% ------------------------------------------------------------------------------------------
%% 5)   Smoothing

% Smooth M0 map before division
% Currently, same smoothing kernel for 2D & 3D, to remove e.g. Gibbs
% ringing artifact in 3D data. And smoothing sums quadratically

% This should be considerable large smoothing,
% because of the noise that would be introduced by
% voxel-wise division. See Beaumont's thesis, chapter 4

% Initial large smoothing for smooth brain image

MaxIt                   = 24; % usually around this value (used for counting, not for limitation
                                % the number of iterations)

fprintf('Mask M0 with this mask & smooth...  ');
ImOut                   = ImIn.*Mask4;
ImOut(ImOut==0)         = NaN;

% smoothing with interpolation
% The bigger the kernel size, the smoother the biasfield, but the
% more artifacts will be filtered into the data

VoxelSize 				= [1.5 1.5 1.5];

ImOut                   = xASL_im_ndnanfilter(ImOut,'gauss',double([16 16 16]./VoxelSize),0);
xASL_TrackProgress(1,MaxIt);
Im5                     = ImOut;
ImOut                   = xASL_im_ndnanfilter(ImOut,'gauss',double([16 16 16]./VoxelSize),0);
xASL_TrackProgress(2,MaxIt);




%% ------------------------------------------------------------------------------------------
%% 6) Now only extrapolating
ImOut   = xASL_im_ExtrapolateOverNaNs(ImOut); % this should not be masked,
                                                % the idea is to fill the whole FoV
                                                % to prevent ASL/M0 division artifacts



%% ------------------------------------------------------------------------------------------
%% Now calculate back the GM M0
fprintf('\n%s\n','Rescale the smooth biasfield GM M0 values back to non-smooth GM M0 values');
OldGMM0     = ImIn(GMmask & isfinite(ImIn));
OldGMM0     = median(OldGMM0);

NewGMM0     = ImOut(GMmask & isfinite(ImOut));
NewGMM0     = median(NewGMM0);

RatioN      = OldGMM0/NewGMM0;
ImOut       = ImOut.*RatioN;






%% ------------------------------------------------------------------------------------------
%% Print progress

S2S                     = 53; % slice to show
IM                      = [xASL_im_rotate(ImIn(:,:,S2S),90) xASL_im_rotate(ImIn(:,:,S2S).*Mask2(:,:,S2S),90) xASL_im_rotate(ImIn(:,:,S2S).*Mask3(:,:,S2S),90) ; xASL_im_rotate(ImIn(:,:,S2S).*Mask4(:,:,S2S),90) xASL_im_rotate(Im5(:,:,S2S),90) xASL_im_rotate(ImOut(:,:,S2S),90)];

xASL_adm_CreateDir(x.D.M0regASLdir);
OutputFile              = fullfile(x.D.M0regASLdir,['M0_im_proc_' x.P.SubjectID '.jpg']);
fprintf('%s\n',['Writing ' OutputFile]);
xASL_imwrite(IM, OutputFile);

end
