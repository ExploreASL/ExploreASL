function UnitTest = xASL_ut_function_xASL_mex_conv3DsepGauss(TestRepository)
%xASL_ut_function_xASL_mex_conv3DsepGauss Individual unit test for xASL_mex_conv3DsepGauss
%
% INPUT:        TestRepository - Path to test repository.
%
% OUTPUT:       UnitTest  - Test structure
%               name      - Name of tested module or submodule (char array)
%               unit      - Insert one of the following: 'Module', 'Submodule' or 'Function'
%               passed    - Result of all subtests combined (true or false)
%               tests     - Structure with individual subtest results
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Tests the separable three-dimensional Gaussian convolution
%               against a base MATLAB reference implementation.
%
%               The accepted maximum absolute difference is 5% of the
%               maximum absolute signal in the reference output.
%
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_mex_conv3DsepGauss(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________

if nargin < 1
    TestRepository = '';
end


%% Test run 1

UnitTest.tests(1).testname = 'Zero sigma returns unchanged image';

testTime = tic;

imageIn = reshape(1:60, [4, 5, 3]);
sigma = [0, 0, 0];

imageOut = xASL_mex_conv3DsepGauss(imageIn, sigma);
imageReference = imageIn;

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(1).duration = toc(testTime);
UnitTest.tests(1).passed = testCondition;


%% Test run 2

UnitTest.tests(2).testname = 'Gaussian convolution in X dimension';

testTime = tic;

imageIn = reshape(sin(1:120), [5, 6, 4]);
sigma = [0.6, 0, 0];

imageOut = xASL_mex_conv3DsepGauss(imageIn, sigma);
imageReference = xASL_ut_ReferenceConv3DsepGauss(imageIn, sigma);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(2).duration = toc(testTime);
UnitTest.tests(2).passed = testCondition;


%% Test run 3

UnitTest.tests(3).testname = 'Gaussian convolution in Y dimension';

testTime = tic;

imageIn = reshape(cos((1:120) ./ 5), [5, 6, 4]);
sigma = [0, 0.6, 0];

imageOut = xASL_mex_conv3DsepGauss(imageIn, sigma);
imageReference = xASL_ut_ReferenceConv3DsepGauss(imageIn, sigma);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(3).duration = toc(testTime);
UnitTest.tests(3).passed = testCondition;


%% Test run 4

UnitTest.tests(4).testname = 'Gaussian convolution in Z dimension';

testTime = tic;

imageIn = reshape(sin((1:150) ./ 7), [5, 6, 5]);
sigma = [0, 0, 0.6];

imageOut = xASL_mex_conv3DsepGauss(imageIn, sigma);
imageReference = xASL_ut_ReferenceConv3DsepGauss(imageIn, sigma);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(4).duration = toc(testTime);
UnitTest.tests(4).passed = testCondition;


%% Test run 5

UnitTest.tests(5).testname = 'Three-dimensional Gaussian convolution';

testTime = tic;

imageIn = reshape(cos((1:210) ./ 7), [6, 7, 5]);
sigma = [0.6, 0.8, 0.5];

imageOut = xASL_mex_conv3DsepGauss(imageIn, sigma);
imageReference = xASL_ut_ReferenceConv3DsepGauss(imageIn, sigma);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(5).duration = toc(testTime);
UnitTest.tests(5).passed = testCondition;


%% Test run 6

UnitTest.tests(6).testname = 'Constant image and boundary normalization';

testTime = tic;

imageIn = ones(5, 6, 5);
sigma = [0.6, 0.6, 0.6];

imageOut = xASL_mex_conv3DsepGauss(imageIn, sigma);
imageReference = imageIn;

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(6).duration = toc(testTime);
UnitTest.tests(6).passed = testCondition;


%% Test run 7

UnitTest.tests(7).testname = 'Gaussian Z impulse response';

testTime = tic;

imageIn = zeros(3, 3, 9);
imageIn(2, 2, 5) = 1;
sigma = [0, 0, 1];

imageOut = xASL_mex_conv3DsepGauss(imageIn, sigma);
imageReference = xASL_ut_ReferenceConv3DsepGauss(imageIn, sigma);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(7).duration = toc(testTime);
UnitTest.tests(7).passed = testCondition;


%% End of testing

UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


function imageOut = xASL_ut_ReferenceConv3DsepGauss(imageIn, sigma)
%xASL_ut_ReferenceConv3DsepGauss Reference Gaussian convolution

imageOut = imageIn;

for iDimension = 1:3

    if sigma(iDimension) <= 0
        continue;
    end

    kernelRadius = ceil(3 .* sigma(iDimension));
    kernelPosition = -kernelRadius:kernelRadius;

    kernel = exp(-(kernelPosition .^ 2) ./ ...
        (2 .* sigma(iDimension) .^ 2));

    kernel = kernel ./ sum(kernel);

    kernelSize = ones(1, 3);
    kernelSize(iDimension) = numel(kernel);
    kernel = reshape(kernel, kernelSize);

    weightedImage = convn(imageOut, kernel, 'same');
    kernelNormalization = convn(ones(size(imageOut)), kernel, 'same');

    imageOut = weightedImage ./ kernelNormalization;

end

end


function testCondition = xASL_ut_CompareRelative(imageOut, imageReference, relativeTolerance)
%xASL_ut_CompareRelative Compare output using a relative signal tolerance

testCondition = isequal(size(imageOut), size(imageReference));

if ~testCondition
    return;
end

maximumSignal = max(abs(imageReference(:)));
absoluteTolerance = relativeTolerance .* maximumSignal;

if absoluteTolerance == 0
    absoluteTolerance = eps;
end

maximumDifference = max(abs(imageOut(:) - imageReference(:)));

if ~isfinite(maximumDifference) || maximumDifference > absoluteTolerance
    testCondition = false;
end

end