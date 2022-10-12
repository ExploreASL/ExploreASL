function []                 = xASL_wrp_GetClusterCBFstats(x,printTitle,InputDataStr)
%% Create a CSV with CBF info for ROIs from the statistical mask

    % "Labels" here can be understood as "clusters"

    if ~isfield(x.S,'LabelMaskNII')
        x.S.LabelMaskNII      = fullfile(x.S.StatsDir,['Clusters_' printTitle '.nii']);
    end
    if ~isfield(x.S,'unit')
        x.S.unit              = 'mL/100g/min';
    end
    x.S.output_ID             = ['VBA-ROI_' printTitle];

    %% Create cluster "pseudo-ROI" names
    [Fpath Ffile Fext]      = xASL_fileparts(x.S.LabelMaskNII);
    TSVpath                 = fullfile(Fpath,[Ffile '.tsv']);
    FID                     = fopen(TSVpath,'w');
    nROIs                   = length(nonzeros(unique(xASL_io_Nifti2Im(x.S.LabelMaskNII))));
    for iROI=1:nROIs
        fprintf(FID,'%s\t',['Cluster_' num2str(iROI)]);
    end
    fclose(FID);


    %% VBA-derived ROI preparation


    if  exist(x.S.LabelMaskNII,'file')
        % only if regions were found these will be used to report their
        % CBF, otherwise this part is omitted

        %% Here, we make sure the clusters have a minimal size.
		%  First this was done with 512 voxels, but this doesn't take into account the part that will be excluded
		%  by further scripts including only GM. Therefore, now including GM volume
		%
        %% We don't want this by default, only when we have clusters from e.g. a p=0.001 threshold that provides very small clusters
		%  Therefore, this part is included here, rather than in the xASL_wrp_GetROIstatistics.m
        LabelIM     = xASL_io_Nifti2Im(x.S.LabelMaskNII);
        pGM_MNI     = xASL_io_Nifti2Im(fullfile(x.D.TemplateDir,'rc1T1_ASL_res.nii'));
        for iL=1:max(LabelIM(:))
            while  sum(sum(sum((LabelIM==iL).*pGM_MNI)))<768
                   tempL                   = xASL_im_DilateErodeFull(LabelIM==iL,'dilate',xASL_im_DilateErodeSphere(1));
                   LabelIM(logical(tempL)) = iL;
            end
        end
        xASL_io_SaveNifti( x.S.LabelMaskNII, x.S.LabelMaskNII, LabelIM, [], 0);

        if ~exist('InputDataStr','var')
            x.S.InputDataStr            = 'qCBF_untreated'; % 'qCBF' 'SD' 'TT' 'ATT' 'TExch' 'M0' 'R1' 'ASL_HctCohort' 'ASL_HctCorrInd'
        else
            x.S.InputDataStr            = InputDataStr;
        end
        x.S.InputAtlasPath            = x.S.LabelMaskNII;

        fprintf('%s\n',['Printing ' x.S.output_ID ' tsv-files']);
        x.S.output_ID                 = ''; % is already acquired in the next script,
        % by the name of the LabelMaskNII/"pseudo-Atlas" created above
        xASL_wrp_GetROIstatistics( x);

        if  isfield(x.S,'DAT')
            x.S                       = rmfield(x.S,'DAT');
        end
    end
end
