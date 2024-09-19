% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

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
