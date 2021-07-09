%% Create Dummy TestDataSet (digital phantom) 2

% Here we recreate the average timeseries
x = ExploreASL_Initialize('/Users/henk/ExploreASL/ASL/ExploreASL_Manuscript/analysis_EPAD/DATA_PAR_Philips2DEPI_Bsup_040_Only_All.json');

SumASL4D = double(zeros(121,145,121,60));
Counter = 0;

for iSubject=1:x.nSubjects
    xASL_TrackProgress(iSubject, x.nSubjects);
    Path_ASL4Dnii = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'ASL4D.nii');
    Path_ASL4DniiTemp2Delete = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'ASL4D_Temp2Delete.nii');
    Path_ASL4DniiTemp2Delete_StandardSpace = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'ASL4D_Temp2Delete_StandardSpace.nii');
    Path_ASL4Djson = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'ASL4D.json');
    PathY = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'y_ASL.nii');
    
    if ~exist(Path_ASL4Djson, 'file') || ~xASL_exist(Path_ASL4Dnii, 'file')
        fprintf(['iSubject ' x.SUBJECTS{iSubject} ' is missing JSON or ASL4D file(s)\n']); % here we had 040EPAD00014 00064 00066 00085 missing ASL scans
        fprintf('\n  ');
    continue;
    end
    
    JSON = spm_jsonread(Path_ASL4Djson);

    if ~isfield(JSON, 'PhilipsScaleSlope') || ~isfield(JSON, 'PhilipsRescaleSlope')
        warning(['iSubject ' x.SUBJECTS{iSubject} ' is missing scaleslope(s)']);
        continue;
    end
    
    ASL4D = xASL_io_Nifti2Im(Path_ASL4Dnii);
    ASL4D = ASL4D .* JSON.PhilipsScaleSlope .* JSON.PhilipsRescaleSlope;
    
    % Create rescaled native space copy
    if ~xASL_exist(Path_ASL4DniiTemp2Delete, 'file')
        xASL_io_SaveNifti(Path_ASL4Dnii, Path_ASL4DniiTemp2Delete, ASL4D, 32, 0);
    end
    
    % Create standard space copy
    if  ~xASL_exist(PathY, 'file')
        fprintf(['iSubject ' x.SUBJECTS{iSubject} ' is missing y_ASL.nii file(s)\n']);
        fprintf('\n  ');
    continue;
    end
    
    if ~xASL_exist(Path_ASL4DniiTemp2Delete_StandardSpace, 'file')
        xASL_spm_deformations(x, Path_ASL4DniiTemp2Delete, Path_ASL4DniiTemp2Delete_StandardSpace, 2, [], [], PathY);
    end
    
    % Add this image to our sum
    ASL4D = double(xASL_io_Nifti2Im(Path_ASL4DniiTemp2Delete_StandardSpace));
    ASL4D(isnan(ASL4D)) = 0;
    SumASL4D = SumASL4D + ASL4D;
    Counter = Counter+1;
end
    
MeanASL4D = SumASL4D./Counter;

Path_MeanASL4D = fullfile(x.D.TemplatesStudyDir, 'MeanASL4D.nii');
xASL_io_SaveNifti(Path_ASL4DniiTemp2Delete_StandardSpace, Path_MeanASL4D, MeanASL4D, 32, 0);

% Move to native space
PathASL4D_Native = '/Users/henk/ExploreASL/ExploreASL/External/TestDataSet/Sub-001/ASL_1/ASL4D.nii';
InverseSpaceM0 = '/Users/henk/ExploreASL/ExploreASL/External/TestDataSet/Sub-001/ASL_1/M0.nii';
DeformationPath = '/Users/henk/ExploreASL/ExploreASL/External/TestDataSet/Sub-001/y_T1.nii';

xASL_spm_deformations(x, Path_MeanASL4D, PathASL4D_Native, 2, InverseSpaceM0, [], DeformationPath);

% Correct quantification
ScaleFactor = 2400/8;
ASL4D = xASL_io_Nifti2Im(PathASL4D_Native);
xASL_io_SaveNifti(PathASL4D_Native, PathASL4D_Native, ASL4D.* ScaleFactor, [], 0);

M0 = xASL_io_Nifti2Im(InverseSpaceM0);
xASL_io_SaveNifti(InverseSpaceM0, InverseSpaceM0, M0./ 600, [], 0);


% Delete temporary files
for iSubject=1:x.nSubjects
    Path_ASL4DniiTemp2Delete = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'ASL4D_Temp2Delete.nii');
    Path_ASL4DniiTemp2Delete_StandardSpace = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'ASL_1', 'ASL4D_Temp2Delete_StandardSpace.nii');
    xASL_delete(Path_ASL4DniiTemp2Delete);
    xASL_delete(Path_ASL4DniiTemp2Delete_StandardSpace);
end


    
    
    
    
    


%% 0) Copy analyzed TestDataSet from previous study & open it
RootDir = 'C:\TempxASL\TestDataSet';
NewAnalysisDir = fullfile(fileparts(RootDir),'TestDataSetNew');

PathDataParOld = fullfile(RootDir,'DataParameters_HiQ.json');
PathDataParNew = fullfile(NewAnalysisDir,'DataParameters_HiQ.json');

x = ExploreASL_Initialize();
x.settings.Quality = 1;
xASL_adm_CreateDir(NewAnalysisDir);
SubjectDir = fullfile(x.D.ROOT,'Sub-001');
NewSubjectDir = fullfile(NewAnalysisDir,'Sub-001');
NewASLDir = fullfile(NewSubjectDir,'ASL_1');
xASL_adm_CreateDir(NewSubjectDir);
xASL_adm_CreateDir(NewASLDir);

PathM0 = fullfile(NewASLDir,'M0.nii');
PathM0parms = fullfile(NewASLDir,'M0.json');
PathM0parmsOld = fullfile(SubjectDir,'ASL_1','M0.json');
PathASLparmsOld = fullfile(SubjectDir,'ASL_1','ASL4D.json');
PathPWI = fullfile(NewASLDir,'PWI.nii');
PathCBF = fullfile(NewASLDir,'CBF.nii');
PathCBF2 = fullfile(NewASLDir,'CBF2.nii');
Path_SliceGradient = fullfile(NewASLDir,'SliceGradient.nii');
x.P.Path_despiked_ASL4D = fullfile(NewASLDir,'ASL4D.nii');
Path_ASL4DOriginal = fullfile(SubjectDir, 'ASL_1','ASL4D.nii');
x.P.Path_ASL4D = x.P.Path_despiked_ASL4D;
x.P.Path_ASL4D_parms_mat = fullfile(NewASLDir,'ASL4D_parms.mat');
x.P.Path_ASL4D_parms_json = fullfile(NewASLDir,'ASL4D.json');
x.P.Path_mean_control = fullfile(NewASLDir,'mean_control.nii');
MotionMat = fullfile(NewASLDir,'ASL4D.mat');
MotionMatOld = fullfile(SubjectDir,'ASL_1','ASL4D.mat');

xASL_Copy(PathDataParOld, PathDataParNew); % copy datapar
xASL_Copy(PathM0parmsOld, PathM0parms); % copy M0 parms
xASL_Copy(PathASLparmsOld, x.P.Path_ASL4D_parms_json); % copy ASL parms
xASL_Copy(MotionMatOld, MotionMat); % copy motion mat

%% 1) Create native space files from standard space
TemplateDir = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population\Templates';
TemplateFiles = {'T1_bs-mean_Unmasked.nii' 'FLAIR_bs-mean_unmasked.nii' 'noSmooth_M0_bs-mean.nii' 'mean_control_bs-mean.nii' 'CBF_bs-mean.nii'};
NativeFiles = {'T1.nii' 'FLAIR.nii' fullfile('ASL_1','M0.nii') fullfile('ASL_1','mean_control.nii') fullfile('ASL_1', 'CBF.nii')};
DeformatFiles = {'y_T1.nii' 'y_T1.nii' fullfile('ASL_1','y_ASL.nii') fullfile('ASL_1','y_ASL.nii') fullfile('ASL_1','y_ASL.nii')};

for iFile=1:length(TemplateFiles)
    TemplatePath = fullfile(TemplateDir,TemplateFiles{iFile});
    NativePath = fullfile(NewSubjectDir, NativeFiles{iFile});
    InversePath = fullfile(SubjectDir, NativeFiles{iFile});
    DeformatPath = fullfile(SubjectDir, DeformatFiles{iFile});
    
    xASL_spm_deformations(x, TemplatePath, NativePath, 2, InversePath, [], DeformatPath);
end
    
%% 2) Create slice gradient
tNII = xASL_io_ReadNifti(PathM0);
dim = size(tNII.dat(:,:,:,:,:,:,:,:,:));
dim = dim(1:3);
SGim = zeros(dim);

for iSlice=1:dim(3)
    SGim(:,:,iSlice) = iSlice;
end

% We skip motion correction here, to avoid edge artifacts
xASL_io_SaveNifti(PathM0,Path_SliceGradient,SGim,8,0);

%% 3) Obtain the PWI image by running xASL_wrp_Quantify in reverse
xASL_Copy(PathCBF, PathPWI);
xASL_Copy(PathCBF,x.P.Path_ASL4D);
x.Sequence = '2D_EPI';
x.P.SubjectID = 'Sub-001';
x.P.SessionID = 'ASL_1';
xASL_wrp_Quantify(x, PathPWI, PathCBF2, PathM0, Path_SliceGradient);
PWIim = double(xASL_io_Nifti2Im(PathPWI));
M0im = double(xASL_io_Nifti2Im(PathM0));
CBF2im = double(xASL_io_Nifti2Im(PathCBF2));
CBFim = double(xASL_io_Nifti2Im(PathCBF));
FactorIm = (PWIim./M0im)./CBF2im;
TruePWIim = (CBFim .* FactorIm) .* M0im;
xASL_io_SaveNifti(PathPWI, PathPWI, TruePWIim);

% we set the rescale slope in the json sidecar to 1, only to avoid a
% warning for having a different rescale slope in the nifti and json

xASL_delete(x.P.Path_ASL4D);
xASL_delete(PathCBF2);

%% 4) Scale back the M0 (including DICOM scale slopes)
% Here we simply run xASL_quant_M0 in reverse,
% by running it again on the M0 image, but re-applying its before/after
% ratio
% Note that we dont have voxel-size differences between M0 and ASL here, so
% it is fine to only use this function

M0IMbefore = double(xASL_io_Nifti2Im(PathM0));
x.P.Path_M0_parms_mat = PathM0parms;
x.P.Path_M0 = PathM0;
M0IMafter = double(xASL_quant_M0(M0IMbefore, x));

IMScale = M0IMbefore./M0IMafter;
xASL_io_SaveNifti(PathM0,PathM0, M0IMbefore .* IMScale,32,0);

%% 5) Now we reconstruct the label image from mean control & PWI
% PWI = control-label
ControlIm = xASL_io_Nifti2Im(x.P.Path_mean_control);
LabelIm = ControlIm - xASL_io_Nifti2Im(PathPWI);
% get nVol
OriNii = xASL_io_ReadNifti(Path_ASL4DOriginal);
nVol = size(OriNii.dat,4);
% Create timeseries
ASL4D_im(:,:,:,[1:2:nVol-1]) = repmat(ControlIm,[1 1 1 nVol/2]);
ASL4D_im(:,:,:,[2:2:nVol-0]) = repmat(LabelIm,[1 1 1 nVol/2]);

xASL_io_SaveNifti(PathPWI,x.P.Path_ASL4D,ASL4D_im,32,0);

% Apply motion originally estimated
xASL_spm_reslice(x.P.Path_ASL4D, x.P.Path_ASL4D, [], [], 1, x.P.Path_ASL4D, 2);

%% Housekeeping
xASL_delete(Path_SliceGradient);
xASL_delete(PathPWI);
xASL_delete(x.P.Path_mean_control);
xASL_delete(PathCBF);
xASL_delete(MotionMat);
rmdir(x.D.PopDir,'s');

