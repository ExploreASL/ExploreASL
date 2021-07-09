function [bSuccess] = xASL_fsl_TopUp(InDir, ScanType, x, OutputPath)
%xASL_fsl_TopUp Submodule of ExploreASL Structural Module, that performs several visualizations for QC
%
% FORMAT: xASL_fsl_TopUp(InDir[, ScanType], x)
%
% INPUT:
%   InDir       - path to folder containing the input & output NIfTIs of TopUp (REQUIRED)
%   ScanType    - type of scan, current options = func, asl, dwi (OPTIONAL, DEFAULT=check automatically)
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   OutputPath  - path to NIfTIfile with TopUp applied (OPTIONAL, DEFAULT = skip TopUpApply)   
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function runs FSL TopUp. It assumes that there are 2
%              TopUp images, i.e. 1 blip up & 1 blip down.
%
%              0. Admin: manage ScanType, NIfTI paths, create TopUp
%                 parameter file for image to apply TopUp to & for the TopUp NIfTIs,
%                 delete files from previous run, define the image with the
%                 same acquisition parameters as TopUp (does the image
%                 we apply TopUp to, have the Blip up or down?)
%              1. Register images to image that we apply TopUp to
%                 (registration between blip up/down images is performed by
%                 TopUp)
%              2. Run TopUp estimate (i.e. estimate the geometric distortion field from B0 NIfTI &
%                 parameters file), this takes quite long. Also has a x.settings.Quality=0 option that is very fast
%                 but inaccurate, to try out this pipeline part. Before
%                 TopUp, NaNs (e.g. from resampling) are removed from the images
%                 TopUp is run with default settings
%              3. Apply TopUp
%
% EXAMPLE: xASL_fsl_TopUp('/analysis/Sub-001/dwi', [], x);
%
% REFERENCE:   Please reference as:
%              "Data was collected with reversed phase-encode blips, resulting in pairs of images with distortions going in opposite directions. From these pairs the susceptibility-induced off-resonance field was estimated using a method similar to that described in [Andersson 2003] as implemented in FSL [Smith 2004] and the two images were combined into a single corrected one."
%              [Andersson 2003] J.L.R. Andersson, S. Skare, J. Ashburner How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 20(2):870-888, 2003.
%              [Smith 2004] S.M. Smith, M. Jenkinson, M.W. Woolrich, C.F. Beckmann, T.E.J. Behrens, H. Johansen-Berg, P.R. Bannister, M. De Luca, I. Drobnjak, D.E. Flitney, R. Niazy, J. Saunders, J. Vickers, Y. Zhang, N. De Stefano, J.M. Brady, and P.M. Matthews. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 23(S1):208-219, 2004.
%              https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup
% __________________________________
% Copyright (C) 2015-2020 ExploreASL



%% Admin: set ScanType

% Default = insuccesfull run:
bSuccess = false;

if nargin<3 || isempty(x)
    x = struct;
end
if nargin<2 || isempty(ScanType)
    [~, ScanDir] = fileparts(InDir);
    if ~isempty(regexp(lower(ScanDir), '.*(func|fmri).*'))
        ScanType = 'func';
    elseif ~isempty(regexp(lower(ScanDir), '.*(asl|m0).*'))
        ScanType = 'asl';
    elseif ~isempty(regexp(lower(ScanDir), '.*(dwi|dti).*'))
        ScanType = 'dwi';
    else
        warning('Unknown ScanType, skipping TopUp');
        return;
    end
end

switch lower(ScanType)
    case 'asl'
        RegExpStr{1} = '^M0\.nii$'; % normal phase encoding direction
        RegExpStr{2} = '^(ASL4D|M0).*RevPE\.nii$'; % reversed phase encoding direction
        RegExpStr{3} = '^ASL4D(|_run-\d)\.nii$'; % file to apply TopUp to
    case 'dwi'
        RegExpStr{1} = '^dwi.*NormPE\.nii$'; % normal phase encoding direction
        RegExpStr{2} = '^dwi.*RevPE\.nii$'; % reversed phase encoding direction
        RegExpStr{3} = '^dwi(|_run-\d)(|_dwi)(?!.*ADC)\.nii$'; % file to apply TopUp to
    case 'func'
        RegExpStr{1} = '^func.*NormPE\.nii$'; % normal phase encoding direction
        RegExpStr{2} = '^func.*RevPE\.nii$'; % reversed phase encoding direction
        RegExpStr{3} = '^func(|_run-\d)(|_bold)\.nii$'; % file to apply TopUp to
    otherwise
        warning('Unknown ScanType for TopUp, skipping');
        return;
end

%% Admin: set paths
if ~isfield(x.external,'bAutomaticallyDetectFSL')
    x.external.bAutomaticallyDetectFSL = 0;
end

[FSLdir, x] = xASL_fsl_SetFSLdir(x, x.external.bAutomaticallyDetectFSL); % Find the FSL directory
% Pathb0cfg = fullfile(FSLdir, 'etc', 'flirtsch', 'b02b0.cnf');
Pathb0cfg = fullfile(x.opts.MyPath, 'CustomScripts', 'EPAD', 'b02b0.cnf'); % use our own one for reproducibility
PathB0 = fullfile(InDir, 'B0.nii');
PathLog = fullfile(InDir, 'B0.topup_log'); % path & file must be same as PathB0
PathResults = fullfile(InDir, 'TopUp');
PathResults1 = fullfile(InDir, 'TopUp_fieldcoef.nii');
PathResults2 = fullfile(InDir, 'TopUp_movpar.txt');
PathField = fullfile(InDir, 'Field.nii');
PathUnwarped = fullfile(InDir, 'Unwarped.nii');



%% Define NIfTI paths
for iRegExp=1:length(RegExpStr)
    FileList = xASL_adm_GetFileList(InDir, RegExpStr{iRegExp}, 'FPList', [0 Inf]);
    if isempty(FileList)
        warning('Could not find TopUp NIfTI, skipping');
        return;
    elseif length(FileList)>1
        warning('Too many TopUp NIfTIs, using the first');
    end
    PathNII{iRegExp} = FileList{1}; % defining the NIfTI

    [Fpath, Ffile] = xASL_fileparts(PathNII{iRegExp});
    MATfile = fullfile(Fpath,[Ffile '.mat']);
    if xASL_exist(MATfile)
        fprintf('%s\n', ['Deleting ' MATfile ' as FSL will not take this into account']);
        xASL_delete(MATfile);
    end
end

%% Create Parms file for NifTI on which to apply TopUp
fclose all;
PathParms2 = fullfile(InDir, ['TopUpAcqParms_' ScanType '.txt']);
xASL_delete(PathParms2);
FID = fopen(PathParms2, 'wt');
ParmsOutput = ObtainTopUpParms(PathNII{end}, x);
fprintf(FID, ParmsOutput); % obtain TopUp acquisition parms
fclose(FID);


%% Delete files from a previous run
xASL_delete(PathB0);
xASL_delete(PathResults1);
xASL_delete(PathResults2);
xASL_delete(PathField);
xASL_delete(PathUnwarped);
xASL_delete(PathLog);

%% Revert previous ORI storage
if strcmp(ScanType, 'asl')
    PathsORI{1} = fullfile(InDir, 'ASL4D.nii');
    PathsORI{2} = fullfile(InDir, 'M0.nii');
else
    PathsORI{1} = PathNII{end};
end

for iORI=1:length(PathsORI)
    [Fpath, Ffile] = xASL_fileparts(PathsORI{iORI});
    PathOrig = fullfile(Fpath, [Ffile '_ORI.nii']);
    if xASL_exist(PathOrig)
        xASL_Move(PathOrig, PathNII{end}, true);
    end
end

%% Create Parms file for TopUp
PathParms = fullfile(InDir, 'TopUpAcqParms.txt');
xASL_delete(PathParms);
FID = fopen(PathParms, 'wt');


for iRegExp=1:length(RegExpStr)-1 % assuming the last one is the output file
    % Here we print the acquisition parameters line by line
    AcqParms = ObtainTopUpParms(PathNII{iRegExp}, x); % obtain TopUp acquisition parms
    CompareParms{iRegExp} = AcqParms;

    fprintf(FID, AcqParms);
    if iRegExp<length(RegExpStr)-1
        fprintf(FID, '\n');
    end

    %% Copy the NIfTIs & make sure they have a single volume
    TopUpNIIPath{iRegExp} = fullfile(InDir,['TopUp' num2str(iRegExp) '.nii']);
    xASL_Copy(PathNII{iRegExp}, TopUpNIIPath{iRegExp}, true);
    tempIM = xASL_io_Nifti2Im(TopUpNIIPath{iRegExp});
    xASL_io_SaveNifti(TopUpNIIPath{iRegExp}, TopUpNIIPath{iRegExp}, tempIM(:,:,:,1,1,1,1,1,1,1), [], false);
end

fclose(FID);

%% Check if there is a TopUp image that has the same acquisition parms as output image
SameParmsInd = 0;
for iC=1:length(CompareParms)
    if ~(isempty(ParmsOutput) || isempty(CompareParms{iC}))
        if strcmp(ParmsOutput(1:6),CompareParms{iC}(1:6)) % same blip (Up/Down)
            N1 = xASL_str2num(ParmsOutput(7:end));
            N2 = xASL_str2num(CompareParms{iC}(7:end));
            
            if isnan(N1)
                warning('No TotalReadoutTime detected for output image, assuming it is equal to the TopUp images');
                SameParmsInd = iC;
            elseif ~(N1>1.1*N2 || N2>1.1*N1) % check that TotalReadoutTime is not more than 10% different
                SameParmsInd = iC;
            end
        end
    end
end
if SameParmsInd==0
    warning('No TopUp found with same acquisition parameters as output image, skipping TopUp');
    return;
end

%% ========================================================================
% =========================================================================
%% Run TopUp Estimate

%% 1) Register blip up/down images to TopUp output image, if we find a TopUp image which similar
%  acquisition parameters as the output image
if SameParmsInd~=0
    RegisterTopUptoOutput(PathNII, InDir, TopUpNIIPath, SameParmsInd, x);
end


%% 2) Run TopUp Estimate
% Concatenate the AP and PA NIfTIs to a B0 NIfTI:
% THIS COMMAND ASSUMES 2 TOPUP FILES ONLY
% xASL_fsl_RunFSL(['/bin/fslmerge -t ' xASL_adm_UnixPath(PathB0, 1)...
%     ' ' xASL_adm_UnixPath(TopUpNIIPath{1}, 1) ' ' xASL_adm_UnixPath(TopUpNIIPath{2}, 1)], x); % direction to concatenate over, t = time, a = auto

% If this NIfTI has multiple volumes, we assume that the first is the B0/M0
tempImage = xASL_io_Nifti2Im(TopUpNIIPath{1});
tIM = tempImage(:,:,:,1);
% If this NIfTI has multiple volumes, we assume that the first is the B0/M0
tempImage = xASL_io_Nifti2Im(TopUpNIIPath{2});
tIM(:,:,:,2) = tempImage(:,:,:,1);

tIM(isnan(tIM)) = 0; % Remove NaNs, TopUp cannot deal with NaNs

% NB we use the orientation of the first image (avoiding image registration
% between NormPE & RevPE

if xASL_stat_SumNan(tIM(:))==0
    warning('Empty TopUp images, skipping TopUp');
    result1 = -1;
    bSuccess = false;
else
    xASL_io_SaveNifti(TopUpNIIPath{1}, PathB0, tIM, 32, false); % we also put the images in single precision    

    fprintf('\n\n=========================================================================\n')
    fprintf('Running FSL TopUp\n');
    fprintf('=========================================================================\n')

    % Run TopUp (i.e. estimate the geometric distortion field from B0 NIfTI &
    % parameters file), this takes quite long:
    TopUpCommand = ['/bin/topup --imain=' xASL_adm_UnixPath(PathB0, 1)...
        ' --datain=' xASL_adm_UnixPath(PathParms, 1) ' --out=' xASL_adm_UnixPath(PathResults, 1)...
        ' --fout=' xASL_adm_UnixPath(PathField, 1) ' --iout=' xASL_adm_UnixPath(PathUnwarped, 1)...
        ' --config=' xASL_adm_UnixPath(Pathb0cfg, 1) ' --verbose=true'];

    if x.settings.Quality
        ActualCommand = TopUpCommand;
        % keep default settings (which are pretty extensive, many iterations)
    else
        %   to speed up while testing:
        %      % --subsamp -> first subsample for estimation, to avoid local
        % minima, and to speed up, but can only be used for matrices that are
        % dividible by the number of subsamples (so 4 or would work for 16
        % slices, 3 or 5 for 15 slices, but with 17 slices there is no
        % possibility for subsampling)
        %
        %     % find number that we can divide matrix by
        %     tIM = xASL_io_ReadNifti(PathB0);
        %     SampleMM = 1;
        %     for iMM=4:-1:1
        %         if SampleMM==1
        %             if min(round(tIM.hdr.dim(2:4)./iMM)==tIM.hdr.dim(2:4)./iMM)
        %                 SampleMM = iMM;
        %             end
        %         end
        %     end
        % 	TopUpCommand = [TopUpCommand ' --subsamp=' SampleMM ' --fwhm=8 --miter=1 --splineorder=2 --interp=linear'];
        %   Above part was disabled, as it didnt work always
        % HERE WE DISABLE SUBSAMPLING, IT 
        % ActualCommand = [TopUpCommand ' --subsamp=4 --fwhm=8 --miter=1 --splineorder=2 --interp=linear'];
        % [~, result1] = xASL_fsl_RunFSL(ActualCommand, x);

        % if result1~=0 % try again without subsampling
            ActualCommand = [TopUpCommand ' --miter=1 --splineorder=3 --interp=linear'];
            % splineorder needs to be cubical (3), if set to square (2) this
            % will create weird results
            % default on normal quality is '--miter=5 --splineorder=3 --interp=spline
        % end
    end

    [~, result1] = xASL_fsl_RunFSL(ActualCommand, x);
end

if result1==0 % successfull run
    bSuccess = true;
end
    
fprintf('\n');
     % --regrid=0 -> don resample into another space for estimation


%% ========================================================================
% =========================================================================
%% 3) Apply TopUp
if ~exist('OutputPath','var') || isempty(OutputPath)
    % Skipping TopUp application, wasn't required
elseif result1~=0
        fprintf('TopUp application skipped as TopUp failed to run\n');
else
    
    if strcmp(ScanType, 'asl')
        PathApplyTopUp{1} = fullfile(InDir, 'ASL4D.nii');
        PathApplyTopUp{2} = fullfile(InDir, 'M0.nii');

        PathOutput = PathApplyTopUp;
    else
        PathApplyTopUp{1} = PathNII{end};

        PathOutput{1} = OutputPath;
    end
    
    for iTopUp=1:length(PathApplyTopUp)
        
        [Fpath, Ffile] = xASL_fileparts(PathApplyTopUp{iTopUp});
        PathOrig = fullfile(Fpath, [Ffile '_ORI.nii']);
        [Fpath2, Ffile2] = xASL_fileparts(PathOutput{iTopUp});
        if strcmp(fullfile(Fpath, Ffile), fullfile(Fpath2, Ffile2))
            % backup input file when output file has same name
            % but don't overwrite previous backup

            xASL_Move(PathApplyTopUp{iTopUp}, PathOrig);
            PathApplyTopUp{iTopUp} = PathOrig; % avoid .nii/.nii.gz issues
            % NB: HERE WE APPLY TOPUP TO THE RENAMED _ORI FILE, but THIS IS THE
            % SAME AS THE file without _ORI
        end

        % Now we convert the images in single precision before running TopUp
        xASL_io_SaveNifti(PathApplyTopUp{iTopUp}, PathApplyTopUp{iTopUp}, single(xASL_io_Nifti2Im(PathApplyTopUp{iTopUp})), 32, false);

        fprintf('\n\n=========================================================================\n')
        fprintf('%s\n',['Applying FSL TopUp to ' PathApplyTopUp{iTopUp}]);
        fprintf('=========================================================================\n')

        % Apply TopUp (method: Use jacobian modulation (jac) or least-squares resampling (lsr, default)
        [~, result1] = xASL_fsl_RunFSL(['/bin/applytopup --imain=' xASL_adm_UnixPath(PathApplyTopUp{iTopUp}, 1) ' --inindex=1'...
                 ' --datain=' xASL_adm_UnixPath(PathParms2, 1) ' --topup=' xASL_adm_UnixPath(PathResults, 1)...
                 ' --out=' xASL_adm_UnixPath(PathOutput{iTopUp}, 1) ' --method=jac'], x); %  --verbose=true

        if result1~=0
            warning('Apply TopUp failed, skipping negative signal detection');
            bSuccess = false;
        else
            % Check if TopUp provided significant negative results
            IM = xASL_io_Nifti2Im(PathOutput{iTopUp});
            PercNonfinite = sum(~isfinite(IM(:))) / numel(IM);
            IM = IM(isfinite(IM(:)));
            NegativeTreshold = -0.001 * (60/xASL_stat_MeanNan(IM(:)));
            % We can allow a little bit negative values from resampling
            PercNegative = sum(IM<NegativeTreshold)/numel(IM);
            
            fprintf('%s\n', ['TopUp resampling resulted in ' xASL_num2str(PercNegative*100, 2) '% negative voxels (thresholded at ' xASL_num2str(NegativeTreshold)]);
            fprintf('Note that we do not remove this negative signal, to keep the noise distribution intact\n');
            if PercNegative>0.1 || PercNonfinite>0.1
                warning('Significant amount of negative voxels in TopUp result, could be a bug');
            end
        end
    end
end

%% Householding
if x.settings.DELETETEMP
    for iP=1:length(TopUpNIIPath)
        xASL_delete(TopUpNIIPath{iP});
    end
end


end





%% ========================================================================
% =========================================================================
function RegisterTopUptoOutput(PathNII, InDir, TopUpNIIPath, SameParmsInd, x)
%RegisterTopUptoOutput Register & resample NormPE & RevPE to target image
%
% FORMAT: RegisterTopUptoOutput(PathNII, InDir, TopUpNIIPath, SameParmsInd, x)
%
% INPUT:
%   PathNII         - cell containing path to input source NifTIs (REQUIRED):
%                     {1} = blip up
%                     {2} = blip down (or reversed)
%                     {3} = image to apply TopUp to
%   InDir           - path to folder with TopUp input images & results (REQUIRED)
%   TopUpNIIPath    - path to copied NIfTIs where TopUp will be computed from (REQUIRED)
%   SameParmsInd    - index of path to NIfTI with same acquisition parameters/geometric distortion as image we
%                     apply TopUp to (REQUIRED)
%   x               - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT: n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function ensures that the NormPE & RevPE images are
%              registered & resampled to the image we apply topup to.
%              Note that the NormPE & RevPE are copied to separate
%              TopUp.nii, both the copies and the NormPE/RevPE are registered.
%              A) Create average image (if the goal NIfTI has multiple
%                 images, they are averaged (e.g. time-series) into a
%                 temporary NIfTI TempRegIm.nii
%              B) Images are registered
%              C) Images are resampled, .mat files are deleted
%
% EXAMPLE: RegisterTopUptoOutput({'../dwi_NormPE.nii' '../dwi_RevPE.nii' '../dwi_dwi.nii'}, '/analysis/Sub-001/dwi', {'../TopUp1.nii' '../TopUp2.nii'}, 1, x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


[~, PrintFile, PrintExt] = xASL_fileparts(PathNII{SameParmsInd});
[~, PrintFile2, PrintExt2] = xASL_fileparts(PathNII{end});

fprintf('\n\n=========================================================================\n')
fprintf([PrintFile PrintExt ' has similar acquisition parms as the output image\n']);
fprintf([PrintFile2 PrintExt2 ', so we register them now & resample them to the output image space\n']);
fprintf('=========================================================================\n')

%% A) Create temporary average image of the output image
tIM = xASL_io_Nifti2Im(PathNII{end});
TempRegPath = fullfile(InDir, 'TempRegIm.nii');
if size(tIM,4)>1
    tIM = xASL_stat_MeanNan(tIM,4);
    xASL_io_SaveNifti(PathNII{end}, TempRegPath, tIM, [], 0);
    refPath = TempRegPath;
else
    refPath = PathNII{end};
end

%% B) Register
srcPath = {TopUpNIIPath{SameParmsInd}};
OtherList = {};
for iC=1:length(TopUpNIIPath)
    if iC~=SameParmsInd
        OtherList{end+1} = TopUpNIIPath{iC};
    end
end
% Also add the original images to the OtherList
for iAdd=1:length(PathNII)-1
    OtherList{end+1} = PathNII{iAdd};
end

xASL_spm_coreg(refPath, srcPath, OtherList, x);

%% C) Resample
%% PM: should we resample here only the temporary TopUp images, or also the images we apply TopUp to?
% Currently, we resample all images
for iO=1:length(OtherList)
    srcPath{end+1} = OtherList{iO};
end
for iS=1:length(srcPath)
    xASL_spm_reslice(refPath, srcPath{iS}, [],[],[], srcPath{iS});
end

xASL_delete(TempRegPath);
[Fpath, Ffile] = xASL_fileparts(TempRegPath);
xASL_delete(fullfile(Fpath, [Ffile '.mat']));


end




%% ========================================================================
% =========================================================================
function [AcqParms] = ObtainTopUpParms(PathIn, x)
%ObtainTopUpParms Extract TopUp parameters from JSON sidecar
%
% FORMAT: [AcqParms] = ObtainTopUpParms(PathIn, x)
%
% INPUT:
%   PathIn    - path to NIfTI from which we want to know the TopUp
%               parameters (assuming that this file has a JSON sidecar with the same
%               filename) (REQUIRED)
%   x         - struct with parameters of the pipeline
%
% OUTPUT:
%    AcqParms - string with acquisition parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function extracts TopUp parameters from JSON sidecar:
%              PhaseEncodingDirection - i j k (for orientation) and - for
%                                       negative (e.g. i = LeftToRight, i- = RightToLeft)
%                                       in TopUp format, this becomes '1 0 0
%              However, it seems that this directionality goes with the
%              DICOM acquisition orientation, so i = x, j = y, etc. And i =
%              LeftRight, J=AP, K=IZ only in case of transversal
%              acquisition (which is true 99% of the case of 2D EPI
%              acquisitions)
%              TotalReadoutTime       -
%              Using the following steps:
%              A) Load the JSON
%              B) Read the JSON
%              C) Print the parameters in TopUp format
%
% EXAMPLE: AcqParms = ObtainTopUpParms('analysis/Sub-001/dwi/dwi_NormPe.nii');
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


    %% A) Load the JSON
    AcqParms = '';
    [Fpath, Ffile] = xASL_fileparts(PathIn);
    JSONin = fullfile(Fpath, [Ffile '.json']);
    if ~exist(JSONin, 'file')
        warning('JSON file missing, skipping TopUp');
        return;
    end
    json = xASL_adm_LoadParms(JSONin, x);

 	if ~isfield(json,'PhaseEncodingDirection')
        warning(['PhaseEncodingDirection JSON field missing: ' PathIn]);
        json.PhaseEncodingDirection = NaN;
    elseif ~isfield(json,'TotalReadoutTime')
        warning(['TotalReadoutTime JSON field missing: ' PathIn]);
        json.TotalReadoutTime = NaN;
    end

    %% B) Read the JSON
    switch json.PhaseEncodingDirection
        case 'i'
            AcqParms = '1 0 0';
        case 'i-'
            AcqParms = '-1 0 0';
        case 'j'
            AcqParms = '0 1 0';
        case 'j-'
            AcqParms = '0 -1 0';
        case 'k'
            AcqParms = '0 0 1';
        case 'k-'
            AcqParms = '0 0 -1';
        otherwise
            warning(['Unknown PhaseEncodingDirection: ' PathIn]);
            AcqParms = 'n/a n/a n/a';
    end

    %% C) Print the parameters in TopUp format
    AcqParms = [AcqParms ' ' xASL_num2str(json.TotalReadoutTime)];


end