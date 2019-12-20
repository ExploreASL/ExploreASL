function image_out = xASL_im_NormalizeLabelingTerritories( imageIN, GMmask, x)
%xASL_im_NormalizeLabelingTerritories Normalizes per perfusion territory
% mask should be GM mask

    % Load territorymap
    TerritoryMap    = fullfile( x.D.MapsDir, 'LabelingTerritories.nii');
    TerritoryMap    = xASL_io_ReadNifti(TerritoryMap);
    TerritoryMap    = TerritoryMap.dat(:,:,:);

    GMmedian        = xASL_stat_MedianNan(xASL_stat_MedianNan(xASL_stat_MedianNan(imageIN(GMmask))));
    image_out       = nan(size(imageIN,1),size(imageIN,2),size(imageIN,3));

    for ii=1:3
        % Create 3 vascular ROIs
        mask{ii}                = GMmask.*TerritoryMap==ii;
        % Create imageROIs
        imageSub{ii}            = nan(size(imageIN,1),size(imageIN,2),size(imageIN,3));
        imageSub{ii}(mask{ii})  = imageIN(mask{ii});
%         imageSub{ii}  = imageSub{ii}./xASL_stat_MedianNan(imageSub{ii}(:)).*GMmedian; % this is
%         normalization of labeling efficiency within subject but not between
%         subjects (normalize to median GM of subject, not population)
        imageSub{ii}            = imageSub{ii}./xASL_stat_MedianNan(imageSub{ii}(:)).*50;
        % this is normalization of labeling efficiency within and between
        % subjects (normalization to median GM of population, not subject)


        % Assemble image
        image_out(mask{ii})     = imageSub{ii}(mask{ii});
    end


end
