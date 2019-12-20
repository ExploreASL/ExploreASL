function IMout = xASL_im_IsoDataThreshold(IMsource)

% Dynamic range
range = (max(IMsource(:)) - min(IMsource(:)));

% Initialize the threshold
threshold = range/2 + min(IMsource(:));

it = 0;
lastThreshold = range;

% Do 10 iterations and then the change has to be more than 1/100 of the
% dynamic range
while (threshold~=lastThreshold) && ((it<5) || ( abs(threshold-lastThreshold)>(range/100)))
    lastThreshold = threshold;
    it = it+1;
    IMout = IMsource > threshold;
    threshold = (xASL_stat_MeanNan(IMsource(IMout)) + xASL_stat_MeanNan(IMsource(~IMout)) )/2;
end;

return;
