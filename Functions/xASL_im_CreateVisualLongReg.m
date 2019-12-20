function xASL_im_CreateVisualLongReg( x, CurrentSub)
% xASL_im_CreateVisualLongReg Creates for each Other TimePoint (TP):
% First row 1) First TP 2) Other TP
% Second row First TP with normalized difference image between TPs without (1) and with (2) Longitudinal Registration

    fprintf('%s\n','Creating comparison image for images with and without longitudinal registration');

	xASL_adm_CreateDir(x.D.LongRegCheckDir);

    x.S.CorSlices = [];
    x.S.SagSlices = [];
    x.S.TraSlices = x.S.slicesLarge;

    % Create brainmask with TimePoint 1 pGM & pWM
    pGMpath  = fullfile(x.D.PopDir,['rc1T1_' CurrentSub{1} '.nii']);
    pWMpath  = fullfile(x.D.PopDir,['rc2T1_' CurrentSub{1} '.nii']);
    Bmask    = (xASL_io_Nifti2Im(pGMpath) + xASL_io_Nifti2Im(pWMpath))>0.1;

    % This is before LongReg for the first TimePoint,
    % which we will compare all other TimePoints against,
    % for both before {1} & after {2} LongReg
    T1wOri = fullfile( x.D.ROOT, CurrentSub{1}, ['r' x.P.STRUCT '.nii']);
    imOri = RobustScaling(xASL_io_Nifti2Im(T1wOri));

    % Load files
    for iC=2:length(CurrentSub)

        % before LongReg = {1}, after LongReg = {2}
        T1wFile{1} = fullfile( x.D.ROOT, CurrentSub{iC}, ['r' x.P.STRUCT '.nii']);
        T1wFile{2} = fullfile( x.D.ROOT, CurrentSub{iC}, ['r2' x.P.STRUCT '.nii']);

        for ii=1:2
            im{ii} = RobustScaling(xASL_io_Nifti2Im(T1wFile{ii})); % scale
            DiffIM{ii} = abs(imOri - im{ii}); % subtract original from other TimePoint
            MeanIM{ii} = (imOri + im{ii})./2; % mean original & other TimePoint
            AI{ii} = DiffIM{ii} ./ MeanIM{ii} .* 100; % normalize difference image into asymmetry image
            AI{ii} = AI{ii} .*single(Bmask); % mask asymmetry image
        end

        % Make sure both AIs will get the same scale (as they are scaled separately in xASL_im_CreateVisualFig below)
        SortValues = max(AI{1},AI{2});
        SortValues = sort(SortValues(isfinite(SortValues) & Bmask));
        ThreshMax = SortValues(round(0.99*length(SortValues)));
        ThreshMin = median(SortValues);
        for ii=1:2
            % Clip both with the same max value
            AI{ii}(AI{ii}>ThreshMax) = ThreshMax;
            % Remove all below the median asymmetries
            AI{ii}(AI{ii}<ThreshMin) = 0;
            % Square asymmetries, this helps to remove the noise and show only important differences
            AI{ii} = AI{ii}.^2;
        end




        % Plot first TimePoint & other TimePoint
        OutIM1 = xASL_im_CreateVisualFig(x, fullfile(x.D.ROOT, CurrentSub{1}, ['r' x.P.STRUCT '.nii']));
        OutIM2 = xASL_im_CreateVisualFig(x, fullfile(x.D.ROOT, CurrentSub{iC}, ['r2' x.P.STRUCT '.nii']));
        % Plot first TimePoint + AI image
        OutIM3 = xASL_im_CreateVisualFig(x, {im{1}, AI{1}}, [], [], [], {x.S.gray x.S.jet256}, false);
        OutIM4 = xASL_im_CreateVisualFig(x, {im{1}, AI{2}}, [], [], [], {x.S.gray x.S.jet256}, false);

        xASL_imwrite([OutIM1,OutIM2;OutIM3,OutIM4], fullfile( x.D.LongRegCheckDir,[CurrentSub{1} '_vs_' CurrentSub{iC} '.jpg']));

    end
end





function [IMout]    = RobustScaling( IMin)

    SortedIM    = sort(IMin(isfinite(IMin)));
    Value95     = SortedIM(round(1*length(SortedIM))); % 0.95*(length....
    IMout       = IMin./Value95.*500;

end
