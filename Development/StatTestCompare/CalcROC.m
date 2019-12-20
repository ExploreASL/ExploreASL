function [TPR FPR] = CalcROC( GoldPos, H)
%CalcROC Get true positive rate (TPR) & false positive rate (FPR)
% from GoldTruth & binary statistical result

    GoldNeg     = ~GoldPos;

    TruePos     = H .* GoldPos;
    FalsePos    = H .* GoldNeg;

    TPR         = sum(TruePos(:))  ./ sum(GoldPos(:));
    FPR         = sum(FalsePos(:)) ./ (306440 - sum(GoldPos(:)) ); % -> ROI is not complete image but just GM!


end

%     GMmask    = xASL_io_ReadNifti('E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\DARTEL_T1_template.nii');
%     GMmask    = GMmask.dat(:,:,:)>0.5;
%     sum(GMmask(:)) = 306440
%     dip_image(GMmask)
