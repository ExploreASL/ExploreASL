function UnitTest = xASL_ut_function_xASL_mex_conv3Dsep(TestRepository)
%xASL_ut_function_xASL_mex_conv3Dsep Individual unit test for xASL_mex_conv3Dsep
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
% DESCRIPTION:  Tests three-dimensional separable convolution with supplied
%               one-dimensional kernels.
%
%               The accepted maximum absolute difference is 5% of the
%               maximum absolute signal in the reference output.
%
%               Symmetric kernels are used so that convolution and
%               correlation conventions produce the same result.
%
% EXAMPLE:      UnitTests(1) = xASL_ut_function_xASL_mex_conv3Dsep(TestRepository);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________

if nargin < 1
    TestRepository = '';
end


%% Test run 1

UnitTest.tests(1).testname = 'Empty kernels return unchanged image';

testTime = tic;

imageIn = reshape(1:60, [4, 5, 3]);

imageOut = xASL_mex_conv3Dsep(imageIn, [], [], []);
imageReference = imageIn;

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(1).duration = toc(testTime);
UnitTest.tests(1).passed = testCondition;


%% Test run 2

UnitTest.tests(2).testname = 'Convolution in X dimension';

testTime = tic;

imageIn = reshape(sin(1:60), [5, 4, 3]);
kernelX = [1, 2, 1] ./ 4;

imageOut = xASL_mex_conv3Dsep(imageIn, kernelX);
imageReference = xASL_ut_ReferenceConv3Dsep(...
    imageIn, kernelX, [], []);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(2).duration = toc(testTime);
UnitTest.tests(2).passed = testCondition;


%% Test run 3

UnitTest.tests(3).testname = 'Convolution in Y dimension';

testTime = tic;

imageIn = reshape(cos((1:72) ./ 5), [4, 6, 3]);
kernelY = [1, 4, 1] ./ 6;

imageOut = xASL_mex_conv3Dsep(imageIn, [], kernelY);
imageReference = xASL_ut_ReferenceConv3Dsep(...
    imageIn, [], kernelY, []);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(3).duration = toc(testTime);
UnitTest.tests(3).passed = testCondition;


%% Test run 4

UnitTest.tests(4).testname = 'Convolution in Z dimension';

testTime = tic;

imageIn = reshape(sin((1:100) ./ 7), [4, 5, 5]);
kernelZ = [1, 2, 1] ./ 4;

imageOut = xASL_mex_conv3Dsep(imageIn, [], [], kernelZ);
imageReference = xASL_ut_ReferenceConv3Dsep(...
    imageIn, [], [], kernelZ);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(4).duration = toc(testTime);
UnitTest.tests(4).passed = testCondition;


%% Test run 5

UnitTest.tests(5).testname = 'Three-dimensional separable convolution';

testTime = tic;

imageIn = reshape(cos((1:210) ./ 9), [6, 7, 5]);

kernelX = [1, 2, 1] ./ 4;
kernelY = [1, 4, 1] ./ 6;
kernelZ = [1, 6, 1] ./ 8;

imageOut = xASL_mex_conv3Dsep(...
    imageIn, kernelX, kernelY, kernelZ);

imageReference = xASL_ut_ReferenceConv3Dsep(...
    imageIn, kernelX, kernelY, kernelZ);

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(5).duration = toc(testTime);
UnitTest.tests(5).passed = testCondition;


%% Test run 6

UnitTest.tests(6).testname = 'Constant image and boundary normalization';

testTime = tic;

imageIn = ones(5, 6, 5);

kernelX = [1, 2, 1] ./ 4;
kernelY = [1, 4, 1] ./ 6;
kernelZ = [1, 2, 1] ./ 4;

imageOut = xASL_mex_conv3Dsep(...
    imageIn, kernelX, kernelY, kernelZ);

imageReference = imageIn;

testCondition = xASL_ut_CompareRelative(imageOut, imageReference, 0.001);

UnitTest.tests(6).duration = toc(testTime);
UnitTest.tests(6).passed = testCondition;


%% Test run 7

UnitTest.tests(7).testname = 'Z impulse response';

testTime = tic;

imageIn = zeros(3, 3, 7);
imageIn(2, 2, 4) = 1;

kernelZ = [1, 2, 3, 2, 1];
kernelZ = kernelZ ./ sum(kernelZ);

imageOut = xASL_mex_conv3Dsep(imageIn, [], [], kernelZ);
imageReference = xASL_ut_ReferenceConv3Dsep(...
    imageIn, [], [], kernelZ);

testCondition = xASL_ut_CompareRelative(...
    imageOut, imageReference, 0.001);

UnitTest.tests(7).duration = toc(testTime);
UnitTest.tests(7).passed = testCondition;


%% End of testing

UnitTest = xASL_ut_CheckSubtests(UnitTest);

end


function imageOut = xASL_ut_ReferenceConv3Dsep(...
    imageIn, kernelX, kernelY, kernelZ)
%xASL_ut_ReferenceConv3Dsep Reference separable convolution

imageOut = imageIn;
kernelCell = {kernelX, kernelY, kernelZ};

for iDimension = 1:3

    kernel = kernelCell{iDimension};

    if isempty(kernel) || numel(kernel) == 1
        continue;
    end

    kernel = double(kernel(:));

    kernelSize = ones(1, 3);
    kernelSize(iDimension) = numel(kernel);
    kernel = reshape(kernel, kernelSize);

    weightedImage = convn(imageOut, kernel, 'same');
    kernelNormalization = convn(ones(size(imageOut)), kernel, 'same');

    imageOut = weightedImage ./ kernelNormalization;

end

end


function testCondition = xASL_ut_CompareRelative(...
    imageOut, imageReference, relativeTolerance)
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