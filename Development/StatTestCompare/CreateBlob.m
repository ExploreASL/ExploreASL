function Kernel = CreateBlob( FwHm, SliceDim, NormValue, PeakValue, NoiseFactor, CutOffLevel, TranslVox )
%CreateBlob Creates 2D Gaussian kernel blob
% as artificial signal
% SliceDim (e.g. [121 145]) & TranslVox should be vector (e.g. [10 30] == 10 right 30 up)
% + = right & up
% - = left  & down
% SliceDim & TranslVox should be integers

    %% Throw error if SliceDim or TranslVox are no integers
    if      max(SliceDim~=round(SliceDim))
            error('SliceDim should be integers!!!');
    elseif  max(TranslVox~=round(TranslVox))
            error('TranslVox should be integers!!!');
    end

    %% Expand dimensions if there is a shift requested
    if  exist('TranslVox','var')

            SliceDimTemp    = SliceDim + (2.*abs(TranslVox));
    else    SliceDimTemp    = SliceDim;
    end    

    %% Factor from FWHM to std
    FwHm2SD     = (2*(2*reallog(2))^0.5);
    SD          = FwHm/FwHm2SD;

    %% Create kernel
    Kernel      = fspecial('gaussian',SliceDimTemp,SD); 
    Kernel      = Kernel ./ sum(Kernel(:)); % normalization

    if exist('NormValue','var')
       Kernel      = Kernel .* NormValue;
    end
    if exist('PeakValue','var')
        Kernel      = Kernel .* (PeakValue ./ max(Kernel(:)) );
    end

    %% Cut-off to determine number of voxels with significant change
    if exist('CutOffLevel','var')
        Kernel(Kernel <(CutOffLevel.*max(Kernel(:))) ) = 0;
    end

    %% Crop (& implicitly translate) if shift requested
    if  exist('TranslVox','var')
        IndexMin    = abs(TranslVox)+TranslVox+1;
        IndexMax    = SliceDim + abs(TranslVox) + TranslVox;

        Kernel      = Kernel( IndexMin(1):IndexMax(1),IndexMin(2):IndexMax(2) );

    end

    %% Add noise

    randNoise   = rand(size(Kernel,1),size(Kernel,2)) .* NoiseFactor;
    Kernel      = Kernel + randNoise;


end