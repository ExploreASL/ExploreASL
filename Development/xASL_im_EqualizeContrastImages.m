function xASL_im_EqualizeContrastImages( IMfile1, IMfile2 )
%xASL_im_EqualizeContrastImages clips images, removes extreme values, &
%increases/reduces contrast of the two images to make contrast similar,
% based on a sorted intensities mask. Tested on M0/mean control images

%% Load images
IMfile{1}       = IMfile1;
IMfile{2}       = IMfile2;

for ii=1:2
    IM{ii}                      = xASL_io_Nifti2Im(IMfile{ii});
    [Fpath Fname{ii} Fext{ii}]  = fileparts(IMfile{ii});
end

fprintf('\n\n\n');
fprintf('%s\n','***');
fprintf('\n');
fprintf('%s\n',['Clipping, equalizing contrast & intensity range for ' Fname{1} Fext{1} ' & ' Fname{2} Fext{2}]);


ThreshClipping      = 0.7875; % lenient, to avoid clipping intra-brain voxels
ThreshMask          = 0.85; % strict, to not include extra-brain voxels

for ii=1:2
    SortInt{ii}     = sort(IM{ii}(isfinite(IM{ii})));
    Clip0{ii}       = SortInt{ii}(round(ThreshMask*length(SortInt{ii})));
    Clip1{ii}       = SortInt{ii}(round(ThreshClipping*length(SortInt{ii}))); % more lenient
    MaskN{ii}       = IM{ii}>Clip0{ii};

    %% Acquire stats before treatment
    IMvalues{ii}    = IM{ii}(MaskN{ii} & isfinite(IM{ii}));
    MinPre{ii}      = min(IMvalues{ii});
    MaxPre{ii}      = max(IMvalues{ii});
    meanPre{ii}     = xASL_stat_MeanNan(IMvalues{ii});
    medianPre{ii}   = xASL_stat_MedianNan(IMvalues{ii});
    stdPre{ii}      = xASL_stat_StdNan(IMvalues{ii});

    %% Remove extreme values
    IM{ii}          = xASL_im_ClipExtremes(IM{ii},0.9975,1,0);

    %% Clip lenient mask
    IM{ii}(IM{ii}<Clip1{ii})  = Clip1{ii};
    IM{ii}          = IM{ii} -  Clip1{ii};

    %% Set minimal value to 0
    IM{ii}(IM{ii}<min(IM{ii}(:)) )  = min(IM{ii}(:));
    IM{ii}          = IM{ii} -  min(IM{ii}(:));

    %% Re-obtain intensities
    SortInt{ii}     = sort(IM{ii}(isfinite(IM{ii})));
    MaxInt{ii}      = max(SortInt{ii});
end

%% Rescale maximum from both images
MaxN                = max([MaxInt{1} MaxInt{2}]); % max from both images

for ii=1:2
    IM{ii}          = IM{ii} ./ MaxInt{ii} .* MaxN;

    %% Equalize contrast for both images, in equal approximation,
    %% to avoid too much noise increase, and avoid too much contrast decrease

    %% Equalize contrast - initialize
    SortInt{ii}     = sort(IM{ii}(isfinite(IM{ii})));
    Clip0{ii}       = SortInt{ii}(round(ThreshMask*length(SortInt{ii})));
    MaskN{ii}       = IM{ii}>Clip0{ii};
    IMColumn{ii}    = IM{ii}(MaskN{ii} & isfinite(IM{ii}));
    spatCoV{ii}     = std(IMColumn{ii})/mean(IMColumn{ii});
    spatCoVpre{ii}  = spatCoV{ii};

    N{ii}           = 1; % starting power
end

% Parameters
AllowN          = 0.1; % allow 10% change in spatial CoV between images
StepSize        = 0.01;

% determine which image has higher spatial CoV
if  spatCoV{1}>spatCoV{2}
    MaxCovN         = 1;
    MinCovN         = 2;
else
    MaxCovN         = 2;
    MinCovN         = 1;
end

while   spatCoV{MaxCovN}>spatCoV{MinCovN}*(1+AllowN)
        for ii=1:2
            IMColumn{ii}    = IMColumn{ii}.^N{ii};
            spatCoV{ii}     = std(IMColumn{ii})/mean(IMColumn{ii});
        end

        N{MaxCovN}   = N{MaxCovN} - StepSize; % lower power high contrast image
        N{MinCovN}   = N{MinCovN} + StepSize; % higher power low contrast image
end

for ii=1:2
    %% Apply contrast changes
    IM{ii}           = IM{ii}.^N{ii};
    MaxInt{ii}       = max(IM{ii}(:));
end

%% Rescale to max of both images
MaxN                = max([MaxInt{1} MaxInt{2}]); % max from both images

for ii=1:2
    IM{ii}          = IM{ii} ./ MaxInt{ii} .* MaxN;

    %% Obtain stats after treatment
    SortInt{ii}     = sort(IM{ii}(isfinite(IM{ii})));
    Clip0{ii}       = SortInt{ii}(round(ThreshMask*length(SortInt{ii})));
    MaskN{ii}       = IM{ii}>Clip0{ii};

    %% Acquire stats before treatment
    IMvalues{ii}     = IM{ii}(MaskN{ii} & isfinite(IM{ii}));
    MinPost{ii}      = min(IMvalues{ii});
    MaxPost{ii}      = max(IMvalues{ii});
    meanPost{ii}     = xASL_stat_MeanNan(IMvalues{ii});
    medianPost{ii}   = xASL_stat_MedianNan(IMvalues{ii});
    stdPost{ii}      = xASL_stat_StdNan(IMvalues{ii});

    %% Save image
    xASL_io_SaveNifti(IMfile{ii},IMfile{ii},IM{ii});

    %% Report changes

    fprintf('\n');
    fprintf('%s\n',['For ' Fname{ii} Fext{ii} ':      median from ' num2str(medianPre{ii},3) ' to ' num2str(meanPost{ii},3)]);
    fprintf('%s\n',['Mean +/- SD from ' num2str(meanPre{ii},3) ' +/- ' num2str(stdPre{ii},3) ' -> ' num2str(meanPost{ii},3) ' +/- ' num2str(stdPost{ii},3)]);
    fprintf('%s\n',['[min max] from [' num2str(MinPre{ii},3) ':' num2str(MaxPre{ii},3) '] -> [' num2str(MinPost{ii},3) ':' num2str(MaxPost{ii},3) ']']);
    fprintf('%s\n',['Spatial CoV from ' num2str(100*spatCoVpre{ii},3) '% -> ' num2str(100*spatCoV{ii},3)]);
end

end
