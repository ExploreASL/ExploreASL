function [x] = xASL_stat_GetROIstatistics(x)
%xASL_stat_GetROIstatistics Compute statistics for each ROI
%
% FORMAT: [x] = xASL_stat_GetROIstatistics(x)
%
% INPUT:
%   x                            - struct containing statistical pipeline environment parameters (REQUIRED)
%   x.S.InputMasks               - ROI masks to compute statistics for (REQUIRED)
%   x.WBmask                     - WholeBrain mask used to convert image to column &
%                                  vice versa (ExploreASL compression method) (REQUIRED)
%   x.S.InputDataStr             - prefix of data files (i.e. before SubjectSession
%                                  name), to specify which data to compute ROI
%                                  statistics for (e.g. 'qCBF', 'M0' 'mrc1T1') (REQUIRED)
%   x.S.NamesROI                 - names of ROIs (REQUIRED)
%   x.S.SubjectWiseVisualization - true to visualize effect of PVC expansion & ROI
%                                  creation. Takes computational time
%                                  (OPTIONAL, DEFAULT=off)
%   x.S.ASL                      - true for ASL data, which will take into account its limited resolution
%                                  & regions of poor SNR: apply PVC, PVC expansion, avoid too
%                                  small ROIs, add vascular & susceptibility artifacts masks (OPTIONAL, DEFAULT=true)
% OUTPUT:
%   x           - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function computes mean and spatial CoV for each ROI,
%              in a [1.5 1.5 1.5] mm MNI space,
%              with several ASL-specific adaptions:
%              1) Skip ROI masks that are smaller than 1 mL (i.e. 296 voxels)
%                 as this would be too noisy for ASL (skipped when x.S.IsASL==false)
%              2) Expand each ROI mask such that it has sufficient WM
%                 content (skipped when IsASL==false)
%              3) Create for each ROI mask a left, right and bilateral copy
%              4) Iterate over all subjects:
%                 a) Load partial volume maps
%                 b) Correct for WMH SEGM -> IS THIS STILL REQUIRED???
%                 c) Load data
%                 d) Show ROIs projected on ASL image
%                 e) Actual data computations
%                    Partial Volume Correction (PVC) options:
%                    PVC==0 -> perform masking only, no regression
%                    PVC==1 -> single compartment regression, for pGM
%                    PVC==2 -> dual compartment regression for pGM & pWM (i.e. normal
%                               PVC)
%                    Here we mask out susceptibility artifacts (including
%                    outside FoV) for all ASL computations, and also mask
%                    out vascular artifacts for computing the mean.
%              While other artifacts/FoV can be masked out on population
%              level (i.e. >0.95 subjects need to have a valid signal in a
%              certain voxel), vascular artifacts differ too much in their
%              location between subjects, so we mask this on subject-level.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_stat_GetROIstatistics(x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% ------------------------------------------------------------------------------------------------------------
%% Administration
x.LeftMask = x.WBmask;
x.LeftMask(1:60,:,:) = 0;
x.LeftMask = xASL_im_IM2Column(x.LeftMask,x.WBmask);

if iscell(x.S.InputDataStr)
    x.S.InputDataStr = x.S.InputDataStr{1};
end

x.S = RemoveSuffixes(x.S); % Remove previous output_ID suffixes

if ~isfield(x,'LabEffNorm')
    x.LabEffNorm = false; % default
end
x.S.InputMasks = logical(x.S.InputMasks);

% Define number of sessions to use
nSessions = 0;
NoSessions = 0;
SessionList = xASL_adm_GetFileList(x.D.PopDir,['^' x.S.InputDataStr '.*_ASL_' num2str(nSessions+1) '\.nii'], 'FPList', [0 Inf]);

while ~isempty(SessionList)
       nSessions = nSessions+1;
       SessionList = xASL_adm_GetFileList(x.D.PopDir,['^' x.S.InputDataStr '.*_ASL_' num2str(nSessions+1) '\.nii'], 'FPList', [0 Inf]);
end

if nSessions==0 % check if there are "subject-files"
    SubjectList = xASL_adm_GetFileList(x.D.PopDir,['^' x.S.InputDataStr '.*\.nii'], 'FPList', [0 Inf]);
    if isempty(SubjectList)
        fprintf('%s\n','No session or subject files found');
        return;
    end
    nSessions = 1;
    NoSessions = 1;
end

% Determine whether group mask exists
if isfield(x.S,'MaskSusceptibility') && ~min(x.S.MaskSusceptibility == xASL_im_IM2Column(ones(121,145,121),x.WBmask))
    HasGroupSusceptMask = true;
else
    HasGroupSusceptMask = false;
end

if HasGroupSusceptMask
    fprintf('%s\n', 'Using population-based susceptibility mask');
elseif ~strcmp(x.Sequence,'3D_spiral') % for 3D spiral we dont need a susceptibility mask
    fprintf('%s\n', 'Using subject-specific susceptibility mask');
end

fprintf('%s\n',['Preparing ROI-based ' x.S.output_ID ' statistics:']);


if x.S.IsASL
    %% ------------------------------------------------------------------------------------------------------------
    %% 1) For all ROIs, skip ROIs smaller than 1 mL (296 voxels @ 1.5x1.5x1.5 mm)
    VoxelSize = [1.5 1.5 1.5];
    MinVoxels = 1000/prod(VoxelSize);
    SumList = squeeze(sum(x.S.InputMasks,1));
    SumList = SumList>MinVoxels;
    if size(SumList,1)==1
        SumList = SumList';
    end



    %% ------------------------------------------------------------------------------------------------------------
    %% 2) For all ROIs, expand ROIs to contain sufficient pWM for PVEc
    %  And this can also be used for the PVC==0 results, which are masked
    %  by individual pGM anyway. This actually could improve individual accuracy/precision,
    %  since it weighs more to the individual masking, compared to the
    %  original atlas which was more narrow around the GM, weighing more
    %  towards an atlas-based masking. However, these results didn"t seem
    %  to change much

    pGM_MNI = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir, 'rc1T1_ASL_res.nii'));
    pWM_MNI = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir, 'rc2T1_ASL_res.nii'));

    fprintf('Expanding ROIs for PVC:   ');
    for iROI=1:size(x.S.InputMasks,2)
        xASL_TrackProgress(iROI,size(x.S.InputMasks,2));
        for iSubject=1:size(x.S.InputMasks,3)
            if sum(SumList(iROI,iSubject))~=0 % skip empty ROIs
                x.S.InputMasks(:,iROI,iSubject) = xASL_im_CreatePVEcROI(x,x.S.InputMasks(:,iROI,iSubject), pGM_MNI, pWM_MNI);
            end
        end
    end
end


%% 3) Create for each ROI mask a left, right and bilateral copy
%  Time is in file loading, not the computation;
% It assumes atlas symmetry though, which is true for most atlases
% For individual atlases (e.g. Lesion/ROI) this part may be skipped
if ~isfield(x.S,'NamesROI')
    warning('No ROIs found');
    return;
end
NamesROIold = x.S.NamesROI;
x.S.NamesROI = {''};
for iR=1:length(NamesROIold)
    x.S.NamesROI{iR*3-2} = [NamesROIold{iR} '_B'];
    x.S.NamesROI{iR*3-1} = [NamesROIold{iR} '_L'];
    x.S.NamesROI{iR*3-0} = [NamesROIold{iR} '_R'];
end

NewSize = [size(x.S.InputMasks,1) size(x.S.InputMasks,2)*3 size(x.S.InputMasks,3)];
x.S.InputMasksTemp = false(NewSize);
x.S.InputMasksTemp(:,[1:3:end-2],:)           = x.S.InputMasks;
x.S.InputMasksTemp( x.LeftMask,[2:3:end-1],:) = x.S.InputMasks( x.LeftMask,:,:); % left  only
x.S.InputMasksTemp(~x.LeftMask,[3:3:end  ],:) = x.S.InputMasks(~x.LeftMask,:,:); % right only
x.S.InputMasks                                = x.S.InputMasksTemp;


%% 4) Iterate over all subjects

fprintf('%s\n','Computing ROI data:   ');
SubjectSpecificMasks = x.S.InputMasks(:,:,1); % changed below if size(x.S.InputMasks,3)>1
% Here we assume there are no subject-specific masks
for iSubject=1:x.nSubjects
    xASL_TrackProgress(iSubject, x.nSubjects);

    if x.S.IsASL
        %% a) Load partial volume maps
        PathGM = fullfile(x.D.PopDir, ['PV_pGM_' x.SUBJECTS{iSubject} '.nii']);
        if xASL_exist(PathGM,'file')
            pGM = xASL_im_IM2Column(xASL_io_Nifti2Im(PathGM),x.WBmask);
        else
            fprintf('%s\n',[PathGM ' missing...']);
            continue;
        end

        PathWM = fullfile(x.D.PopDir, ['PV_pWM_' x.SUBJECTS{iSubject} '.nii']);
        if xASL_exist(PathWM,'file')
            pWM = xASL_im_IM2Column(xASL_io_Nifti2Im(PathWM),x.WBmask);
        else
            fprintf('%s\n',[PathWM ' missing...']);
            continue;
        end

        %% b) Correct for WMH SEGM -> IS THIS STILL REQUIRED???
        WMHfile = fullfile(x.D.PopDir, ['rWMH_SEGM_' x.SUBJECTS{iSubject} '.nii']);
        if xASL_exist(WMHfile, 'file')
            pWMH = xASL_im_IM2Column(xASL_io_Nifti2Im(WMHfile), x.WBmask);
            % we take sqrt(pWMH) to mimic the effective resolution of ASL (instead of smoothing)
            pWMH(pWMH<0) = 0;
            pGM = max(0, pGM - pWMH.^0.67);
            pWM = max(0, pWM - pWMH.^0.67);
        end
    end

    if size(x.S.InputMasks,3)>1
        SubjectSpecificMasks = x.S.InputMasks(:,:,iSubject);
    end

    for iSess=1:nSessions % NB this is not x.nSession
        % but the one above that can be forced to a single nSessions
        % in case of volume or TT

        % ID (which name, group etc), all for identification
        % Subject_session definition
        DataIm = NaN;
        SubjSess = (iSubject-1)* nSessions +iSess;
        if NoSessions
            x.S.SUBJECTID{SubjSess,1} = x.SUBJECTS{iSubject};
        else
            x.S.SUBJECTID{SubjSess,1} = [x.SUBJECTS{iSubject} '_' x.SESSIONS{iSess}];
        end

        %% c) Load data
        FilePath = fullfile(x.D.PopDir, [x.S.InputDataStr '_' x.S.SUBJECTID{SubjSess,1} '.nii']);
        if xASL_exist(FilePath,'file')
            Data3D = xASL_io_Nifti2Im(FilePath,[121 145 121]);
            DataIm = xASL_im_IM2Column(Data3D,x.WBmask);
        end

        % Load vascular mask (this is done subject-wise)
        FilePath = fullfile(x.D.PopDir, ['MaskVascular_' x.S.SUBJECTID{SubjSess,1} '.nii']);
        if xASL_exist(FilePath,'file')
            VascularMask = xASL_im_IM2Column(logical(xASL_io_Nifti2Im(FilePath)), x.WBmask);
        end

        if numel(DataIm)==1
            DataIm = zeros([sum(x.WBmask(:)) 1],'single');
            VascularMask = logical(DataIm);
            DataIm(:) = NaN;
        end

        %% d) Show ROIs projected on ASL image
        if x.S.SubjectWiseVisualization % this takes extra computation time, hence best switched off
            % Prepare visualization settings
            x.S.TraSlices = x.S.slices;
            x.S.CorSlices = [110 90 x.S.slices(1:2)];
            x.S.SagSlices = x.S.slices;

            LabelIM = xASL_Convert4D_3D_atlas(xASL_im_Column2IM(SubjectSpecificMasks(:,[1:3:end]), x.WBmask));
            LabelIM = xASL_im_TransformData2View(LabelIM);
            DataIM = xASL_im_TransformData2View(Data3D);
            CombiIM = xASL_im_ProjectLabelsOverData(DataIM, LabelIM, x);

            xASL_adm_CreateDir(x.S.CheckMasksDir);
            xASL_imwrite(CombiIM, fullfile(x.S.CheckMasksDir,[x.S.output_ID(1:end-16) '_' x.S.SUBJECTID{SubjSess,1} '.jpg']));
        end

%         % Labeling efficiency normalization
%         if  x.LabEffNorm; temp = xASL_im_NormalizeLabelingTerritories( temp, logical(x.masks.Data.data(:,iSub,1)), x); end

       %% e) Actual data computations

        SusceptibilityMask = xASL_im_IM2Column(x.WBmask, x.WBmask); % default = no susceptibility mask
        if HasGroupSusceptMask % use population-based susceptibility mask
            SusceptibilityMask = x.S.MaskSusceptibility;
        elseif strcmp(x.Sequence,'3D_spiral')
            % if 3D spiral, then we dont need a susceptibility mask
        else % fall back to try subject-wise Susceptibility Masks
            FilePath = fullfile(x.D.PopDir, ['rMaskSusceptibility_' x.S.SUBJECTID{SubjSess,1} '.nii']);
            if xASL_exist(FilePath,'file')
                SusceptibilityMask = xASL_im_IM2Column(xASL_io_Nifti2Im(FilePath), x.WBmask);
            else
                if x.S.IsASL
                    fprintf('%s\n',[FilePath ' missing...']);
                end
                if isfield(x.S,'MaskSusceptibility')
                    SusceptibilityMask = x.S.MaskSusceptibility;
                elseif x.S.IsASL % fall back to legacy/backward compatibility  option
                    warning('No susceptibility mask found, using legacy option, no susceptibility masking');
                end
            end
        end



       for iROI=1:size(SubjectSpecificMasks,2)
            % Swap tissue compartments if needed (this is easier scripting-wise)
            %
            % With ASL normally we are interested in GM only
            % which is why we compute the pGM of a ROI throughout this script
            % For the case of WholeBrain or WM, we swap the tissue compartments
            % such that the pGM represents the WholeBrain or WM

            if ~x.S.IsASL
                pGM_here = ones(size(DataIm)); % no masking
            elseif ~isempty(findstr(lower(x.S.NamesROI{iROI}),'whole brain')) || ~isempty(findstr(lower(x.S.NamesROI{iROI}),'wb')) || ~isempty(findstr(lower(x.S.NamesROI{iROI}),'wholebrain'))
                pGM_here = pGM+pWM;
                pWM_here = ones(size(pGM)) - pGM - pWM; % CSF
            elseif  ~isempty(findstr(lower(x.S.NamesROI{iROI}),'white matter')) || ~isempty(findstr(lower(x.S.NamesROI{iROI}),'wm'))  || ~isempty(findstr(lower(x.S.NamesROI{iROI}),'whitematter'))
                % flip the ROIs:
                % we are interested in WM instead of pGM
                pGM_here = pWM;
                pWM_here = pGM;
            else % keep ROIs as they are
                pGM_here = pGM;
                pWM_here = pWM;
            end

            if ~x.S.IsASL
                CurrentMask = SubjectSpecificMasks(:,iROI) & isfinite(DataIm);
                if x.S.IsVolume
                    VoxelSize = 1.5*1.5*1.5;
                    x.S.DAT_sum_PVC0(SubjSess,iROI) = sum(DataIm(CurrentMask)) .* VoxelSize;
                else
                    x.S.DAT_median_PVC0(SubjSess,iROI) = xASL_stat_ComputeMean(DataIm, CurrentMask, 0, 0);
                    x.S.DAT_CoV_PVC0(SubjSess,iROI) = xASL_stat_ComputeSpatialCoV(DataIm, CurrentMask, 0, 0);
                end
            else
                CurrentMask = logical(bsxfun(@times,single(SubjectSpecificMasks(:,iROI)),pGM_here>0.5 & SusceptibilityMask));
                %% CoV
                x.S.DAT_CoV_PVC0(SubjSess,iROI) = xASL_stat_ComputeSpatialCoV(DataIm, CurrentMask, 296, 0);
                % x.S.DAT_CoV_PVC1(SubjSess,iROI) = xASL_stat_ComputeSpatialCoV(DataIm, CurrentMask, 296, 1, pGM_here); % PVC==1, "single-compartment" PVC (regress pGM only)
                x.S.DAT_CoV_PVC2(SubjSess,iROI) = xASL_stat_ComputeSpatialCoV(DataIm, CurrentMask, 296, 2, pGM_here, pWM_here); % PVC==2, "dual-compartment" (full) PVC (regress pGM & pWM)

                %% CBF (now remove vascular artifacts)
                CurrentMask = CurrentMask & VascularMask;
                x.S.DAT_median_PVC0(SubjSess,iROI) = xASL_stat_ComputeMean(DataIm, CurrentMask, 296, 0);
                % x.S.DAT_mean_PVC1(SubjSess,iROI) = xASL_stat_ComputeMean(DataIm, CurrentMask, 296, 1, pGM_here); % PVC==1, "single-compartment" PVC (regress pGM only)
                x.S.DAT_mean_PVC2(SubjSess,iROI) = xASL_stat_ComputeMean(DataIm, CurrentMask, 296, 2, pGM_here, pWM_here); % PVC==2, "dual-compartment" (full) PVC (regress pGM & pWM)

                %% Diff_CoV, new parameter by Jan Petr
                % x.S.DAT_Diff_CoV_PVC0(SubjSess,iROI) = xASL_stat_ComputeDifferCoV(DataIm, CurrentMask, 0);
                % x.S.DAT_Diff_CoV_PVC1(SubjSess,iROI) = xASL_stat_ComputeDifferCoV(DataIm, CurrentMask, 1, pGM_here, [], 0);
                % x.S.DAT_Diff_CoV_PVC2(SubjSess,iROI) = xASL_stat_ComputeDifferCoV(DataIm, CurrentMask, 2, pGM_here, pWM_here, 0);

            end
        end % for iROI=1:size(SubjectSpecificMasks,2)
    end % for iSess=1:nSessions
end % for iSub=1:x.nSubjects


fprintf('\n');

end







%% ------------------------------------------------------------------------------------------------------------
function [S] = RemoveSuffixes(S)
%RemoveSuffixes Clean output_ID by removing previous suffixes
%
% FORMAT: [S] = RemoveSuffixes(S)
%
% INPUT:
%   S           - struct containing statistical pipeline environment parameters (REQUIRED)
%
% OUTPUT:
%   S           - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function removes strings from previous output_ID usage,
%              allowing to reuse this field
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_im_CreateAnalysisMask(x);
% __________________________________
% Copyright 2015-2019 ExploreASL

% Define here which strings to remove
StrRemove = {'PVEC0' 'PVEC1' 'PVEC2' 'AI' 'paired t-test'};

for iStr=1:length(StrRemove)
    INDEXfind = strfind(S.output_ID, StrRemove{iStr} );
    if ~isempty(INDEXfind)
        S.output_ID = [S.output_ID(1:INDEXfind-2) '_' S.output_ID(INDEXfind+length(StrRemove{iStr})+1:end)];
    end
end

UnderScoreFound = strfind(S.output_ID, '_' );
if ~isempty(UnderScoreFound)
    if UnderScoreFound(length(UnderScoreFound))==length(S.output_ID)
        S.output_ID = S.output_ID(1:end-1);
    end
end


end




%% ------------------------------------------------------------------------------------------------------------
function [ROI] = xASL_im_CreatePVEcROI(x, ROI, pGM, pWM)
%xASL_im_CreatePVEcROI Clean output_ID by removing previous suffixes
%
% FORMAT: [ROI] = xASL_im_CreatePVEcROI(x, ROI, pGM, pWM)
%
% INPUT:
%   x     - struct containing pipeline environment parameters (REQUIRED)
%   ROI   - ROI that we wish to expand (REQUIRED)
%   pGM   - GM probability map (this should be in ASL resolution) (REQUIRED)
%   pWM   - WM probability map (this should be in ASL resolution) (REQUIRED)
% OUTPUT:
%   ROI   - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function runs the ROI expansion for PVC, to ensure that
%              there is sufficient pWM inside the ROI for PVC

%              We also remove pCSF from the mask
%              For a ROI that originated from the GM, voxels that contain CSF are shared with GM only.
%              By excluding pCSF>35% we improve PVEc
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ROI] = xASL_im_CreatePVEcROI(x, ROI, pGM, pWM);
% __________________________________
% Copyright 2017-2019 ExploreASL

ROI = ROI>0;

if sum(ROI(:))<296
    % if the ROI is too small, we don't use it for ASL stats, so we can skip this
    return;
end

ROI = xASL_im_Column2IM(ROI, x.WBmask); % convert to image, decompress

ROI = xASL_im_PVC_ROIexpansion(ROI, pGM, pWM, 17); % Dilate ROI
pCSF = max(0,1-pGM-pWM);
ROI(pCSF>0.35) = 0; % Remove pCSF from the mask
ROI = xASL_im_PVC_ROIexpansion(ROI, pGM, pWM, 17); % Dilate ROI

ROI = xASL_im_IM2Column(ROI, x.WBmask); % convert back to column, compress



end








%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
function [ROI] = xASL_im_PVC_ROIexpansion(ROI, pGM, pWM, WMdistMax)
%xASL_im_PVC_ROIexpansion Dilate a GM ROI inwards to include sufficient WM for PVC
%
% FORMAT: [ROI] = xASL_im_PVC_ROIexpansion(ROI, pGM, pWM, WMdistMax)
%
% INPUT:
%   ROI       - ROI binary mask that we wish to expand (NB: this should not be a map) (REQUIRED)
%   pGM       - GM probability map (this should be in ASL resolution) (REQUIRED)
%   pWM       - WM probability map (this should be in ASL resolution) (REQUIRED)
%   WMdistMax - maximum added distance from ROI mask that we want to include WM (in nVoxels).
%               When WMdistMax==0, we compute the WMdistMax for
%               each ROI mask separately (OPTIONAL, DEFAULT=0)
%
% OUTPUT:
%   ROI       - same as input
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function expands a GM ROI mask to include sufficient WM for
%              robust PVC. The rationale is that we often want to obtain GM CBF
%              for which we define an GM ROI, but we need to be sure that
%              this ROI contains sufficient WM content for PVC to work:
%              a) we create GM+CSF masks inside and outside the ROI mask
%              b) for both we compute the distance transform, i.e. a map
%                 stating for each voxel what its distance to the ROI mask is
%              c) Here we dilate the ROI mask until we reach threshold
%                 WMdistMax: i.e. certain nVoxels distance from the ROI mask
%                 We have two options: either 1) predefine the WMdistMax
%                 OR 2) compute the WMdistMax for each ROI separately.
%                 Both options use the same procedure, they only differ in
%                 the stopping criteria. Option 1 is faster and more consistent,
%                 option 2 is perhaps more fair/accurate. Both are
%                 arbitrary/pragmatic solutions.
%
%                 1) dilate ROI mask until we reach threshold WMdistMax
%                 2) dilate ROI mask until we reach the following:
%                    within the ROI:
%                    WM size > 0.5 GM size OR WM > 0.25 GM size & pWM > pGM
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: ROI = xASL_im_PVC_ROIexpansion(ROI, pGM, pWM, 17);
% __________________________________
% Copyright 2017-2019 ExploreASL


%% Admin
GMmask = pGM>0.7;
WMmask = pWM>0.8;

Iteration = 1;
MaxIt = 100;

ROI = logical(ROI);

nVoxelsGMROI = sum(sum(sum(GMmask(ROI))));
nVoxelsWMROI = sum(sum(sum(WMmask(ROI))));
pGMROIsum = sum(sum(sum(pGM(ROI))));
pWMROIsum = sum(sum(sum(pWM(ROI))));

if  xASL_stat_SumNan(ROI(:))==0
    fprintf('%s','Empty ROI, skipping PVEc expansion');
    return;
end

%% a) Create binary GM masks inside/outside ROI
% Take the GM+CSF that is IN and OUT of the ROI
imMaskIN = ROI.*(1-pWM);
imMaskOUT = (1-ROI).*(1-pWM);

% Calculate median GM value in the ROI per slice
 for iSlice = 1:size(ROI,3)
     GMslice = pGM(:,:,iSlice);
     ROIslice = ROI(:,:,iSlice);
     GMvoxels = GMslice(ROIslice);
     GMvoxels = GMvoxels(GMvoxels>0.1);

     if ~isempty(GMvoxels) % Use median GM inside ROI for thresholding the ROIs
         GMmedian = sort(GMvoxels);
         GMmedian = GMmedian(ceil(length(GMmedian)*0.3));
         imMaskIN(:,:,iSlice)  = imMaskIN(:,:,iSlice) > GMmedian;
         imMaskOUT(:,:,iSlice) = imMaskOUT(:,:,iSlice) > GMmedian;
     else % Else threshold by 0.1
         imMaskIN(:,:,iSlice)  = imMaskIN(:,:,iSlice) > 0.1;
         imMaskOUT(:,:,iSlice) = imMaskOUT(:,:,iSlice) > 0.1;
     end
 end

%% b) Compute distance transform
% Now that we have binary masks of GM inside and outside the ROI
% We calculate the distance transform from both
maskINdist = xASL_im_DistanceTransform(imMaskIN);
maskOUTdist = xASL_im_DistanceTransform(imMaskOUT);

%% c) Add WM layers until we reach desired ratio of GM/WM content
if WMdistMax>0 % Add all reasonable WM voxels within the specified maximum distance

    maskAdd = (maskINdist <= WMdistMax); % Find voxels that are distance within the WMdistMax close to the GMmask
    maskAdd = maskAdd .* (maskINdist < maskOUTdist); % If they are not closer to the non-ROI GM
    maskAdd = maskAdd .* (1-ROI); % and if they are not in the ROI yet
    maskAdd = maskAdd .* (pWM > 0.1); % Check if they have WM voxels inside
    ROI = ROI + maskAdd; % Then add them

else % It adds layers of WM until reaching certain ratio of GM and WM size (nVoxels)/content(p) within ROI:
     % size WM/GM > 0.5 OR size WM/GM > 0.25 & content WM>GM

    while Iteration<MaxIt && nVoxelsGMROI>(2*nVoxelsWMROI) && ~(nVoxelsGMROI<(4*nVoxelsWMROI) && pGMROIsum<pWMROIsum)
        % Only add voxels that are closer to the ROI or the GM part of the ROI than they are to surrounding GM structures outside the ROI

        maskAdd = (maskINdist > (Iteration-1)) .* (maskINdist <= Iteration); % Find voxels that are distance Iteration close to the GMmask
        maskAdd = maskAdd .* (maskINdist < maskOUTdist); % If they are not closer to the non-ROI GM
        maskAdd = maskAdd .* (1-ROI); % and if they are not in the ROI yet
        maskAdd = maskAdd .* (pWM > 0.1); % Check if they have WM voxels inside
        ROI = ROI + maskAdd; % Then add them

        % Recalculate the sums
        nVoxelsGMROI = sum(sum(sum(GMmask(logical(ROI)))));
        nVoxelsWMROI = sum(sum(sum(WMmask(logical(ROI)))));
        pGMROIsum = sum(sum(sum(pGM(logical(ROI)))));
        pWMROIsum = sum(sum(sum(pWM(logical(ROI)))));

        Iteration = Iteration+1;
    end
end


end




%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
function [AtlasOut] = xASL_Convert3D_4D_atlas(AtlasIn)
%xASL_Convert3D_4D_atlas Converts 3D atlas with ordinal integers to 4D binary image
%
% FORMAT: [AtlasOut] = xASL_Convert3D_4D_atlas(AtlasIn)
%
% INPUT:
%   AtlasIn   - 3D image matrix with ordinal integers where each integer defines an ROI mask (REQUIRED)
%
% OUTPUT:
%   AtlasOut  - binary 4D matrix, concatenation of 3D binary masks, where
%               fourth dim is equal to the mask integer of AtlasIn
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts 3D atlas with ordinal integers to 4D binary image
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: AtlasOut = xASL_Convert3D_4D_atlas(AtlasIn);
% __________________________________
% Copyright 2017-2019 ExploreASL

AtlasOut = zeros([size(AtlasIn(:,:,:,1,1,1)) max(AtlasIn(:))],'uint8');
for iL=1:max(AtlasIn(:))
    tempIM = zeros(size(AtlasIn(:,:,:,1,1,1)));
    tempIM(AtlasIn==iL) = 1;
    AtlasOut(:,:,:,iL) = tempIM;
end

end



%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
function [AtlasOut] = xASL_Convert4D_3D_atlas(AtlasIn)
%xASL_Convert4D_3D_atlas Converts 4D binary image to 3D atlas with ordinal integers
%
% FORMAT: [AtlasOut] = xASL_Convert4D_3D_atlas(AtlasIn)
%
% INPUT:
%   AtlasIn   - binary 4D matrix, concatenation of 3D binary masks, where
%               fourth dim is equal to the mask number (REQUIRED)
%
% OUTPUT:
%   AtlasOut  - 3D image matrix with ordinal integers where each integer
%               defines an ROI mask, equal to the number of the fourth dim of the input
%               atlas
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts 4D binary image to 3D atlas with ordinal integers
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: AtlasOut = xASL_Convert4D_3D_atlas(AtlasIn);
% __________________________________
% Copyright 2017-2019 ExploreASL

AtlasOut = zeros(size(AtlasIn(:,:,:,1,1,1)),'uint8');
for iL=1:size(AtlasIn,4)
    AtlasOut(logical(AtlasIn(:,:,:,iL))) = iL;
end

end
