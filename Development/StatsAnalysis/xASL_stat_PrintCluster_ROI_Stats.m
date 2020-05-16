function []                 = xASL_stat_PrintCluster_ROI_Stats(MaskMap,x,printTitle)
%% Create a CSV with properties of the statistical mask, with HO-atlases

% "Label" & "cluster" are interchangeable in this script
% This function searches for Harvard-Oxford cortical/subcortical and/or
% cerebellar regions that are positioned within the clusters, using
% a probability cutoff of 25%.

    tMap                    = MaskMap;
    RoiMap                  = MaskMap;
    RoiMap(isnan(RoiMap))   = 0;
    RoiMap                  = logical(RoiMap);

    L                       = spm_bwlabel(double(RoiMap));

        % Label (i.e. assign a number) to each subzero clusters

        % 26 = corner criterion for connectivity (being strict in the definition of "clusters")

    x.S.SaveFile              = fullfile(x.S.StatsDir,['Clusters_' printTitle '.nii']);
    x.S.LabelMaskNII          = x.S.SaveFile;
    xASL_io_SaveNifti( x.D.ResliceRef, x.S.SaveFile, L, 16,0 );
    RegionProp              = regionprops(L,'Centroid');


    % First build cell array, then print to csv-file
    x.S.SaveFile              = fullfile(x.S.StatsDir,['ClusterStats_' printTitle '.tsv']);
    x.S.FID                   = fopen(x.S.SaveFile,'w');

    % Print header
    Header1                 = {'Label' 'VoxelsLabel_(n)' 'VolumeLabel_(L)' 'meanTstat' 'peakTstat' 'XvoxelpeakTstat' 'Y' 'Z' 'XmnipeakTstat' 'Y' 'Z' 'XvoxelCentroid' 'Y' 'Z' 'XmniCentroid' 'Y' 'Z' 'ROIsInLabel_(ROI name_Percentage_Overlap)_sorted for_Percentage_Overlap'};
    for iH=1:length(Header1)
        fprintf(x.S.FID,'%s\t',Header1{iH});
    end
    fprintf(x.S.FID,'\n\n\n');


    % Create single row for each single label ROI
    fprintf('%s','Printing ROI stats...  ');
    for iRoi=1:max(L(:))
        xASL_TrackProgress(iRoi,max(L(:)));
        clear CurrentLabel

        CurrentLabel    = L==iRoi;
        nVoxelsLabel    = sum(CurrentLabel(:));
        VolumeLabel     = nVoxelsLabel*1.5^3/100^3; % assuming voxelsize of 1.5x1.5x1.5 mm, 100^3 goes from mm^3 to dm^3 which is Liters

        x.S.TrackOverlapFound   = CurrentLabel; % this is to track which percentage of
                                            % label volume we could fit atlas ROIs in

        fprintf(x.S.FID,'%s\t', num2str( iRoi ) ); % print label number
        fprintf(x.S.FID,'%s\t', num2str( nVoxelsLabel ) ); % print nVoxels
        fprintf(x.S.FID,'%s\t', num2str( VolumeLabel ) ); % print volume
        fprintf(x.S.FID,'%s\t', num2str( mean(tMap(L==iRoi)) ) ); % print mean t-stat
        fprintf(x.S.FID,'%s\t', num2str( max(tMap(L==iRoi)) ) ); % print peak t-stat

        %% print coordinates of peak t-stat
        [X Y Z] = ind2sub( size(tMap) , find( tMap==max(tMap(L==iRoi))) );
        % CAVE: this will not necessarily present the area with highest
        % t-values, rather the single voxel with highest value. perhaps
        % slight smoothing should be performed here to make this more stable

        % QuickFix: if there are multiple voxels with the same highest
        % t-Value, take the mean value for each dimension (not really
        % correct, should be something like Centroid, but is fine for now)

        if  prod(size(X))>1
            X   = round(mean(X));
        end
        if  prod(size(Y))>1
            Y   = round(mean(Y));
        end
        if  prod(size(Z))>1
            Z   = round(mean(Z));
        end

        x             = PrintXYZcoordinates(x,X,Y,Z); % print voxel indices for peak t-value
        [X Y Z]     = xASL_adm_Voxel2RealWorldCoordinates(X,Y,Z);
        x             = PrintXYZcoordinates(x,X,Y,Z); % print MNI-coordinates for peak t-value

        %% print coordinates of centroid
        % CAVE Centroid outputs as [Y X Z] instead of [X Y Z]
        X       = round(RegionProp(iRoi).Centroid(2));
        Y       = round(RegionProp(iRoi).Centroid(1));
        Z       = round(RegionProp(iRoi).Centroid(3));

        x             = PrintXYZcoordinates(x,X,Y,Z); % print voxel indices for centroid
        [X Y Z]       = xASL_adm_Voxel2RealWorldCoordinates(X,Y,Z);
        x             = PrintXYZcoordinates(x,X,Y,Z); % print MNI-coordinates for centroid

        %% Collect info of HO-regions that lie in the label area

        x.S.iNext                     = 1;
        Atlases                     = {'HOcort_CONN' 'HOsub_CONN'};

        for iA=1:length(Atlases)
            AtlasPath                   = fullfile(x.D.AtlasDir,[Atlases{iA} '.tsv']);
            [TSVfile AtlasROIs]         = xASL_bids_csv2tsvReadWrite(AtlasPath);

            x.S.iA                        = iA;
            for iR=1:length(AtlasROIs)
                x.S.AtlasROIs{iR*2-1}    = [AtlasROIs{iR} '_L']; % left
                x.S.AtlasROIs{iR*2-0}    = [AtlasROIs{iR} '_R']; % right
            end
            x.S.AtlasNames{iA}            = x.S.AtlasROIs;
            AtlasPath                     = fullfile(x.D.AtlasDir,[Atlases{iA} '.nii']);
            x                             = CollectRegionsInfo(x, AtlasPath, CurrentLabel, nVoxelsLabel);
            x.S                           = rmfield(x.S,'AtlasROIs');
        end

        %% If regions were found, print them

        if  isfield(x.S,'ROIunits') % if regions were found

            x.S.ROIunits                  = sortrows(x.S.ROIunits,-4);

            for iUN=1:size(x.S.ROIunits,1)
                 % if the CurrentMask lies within the CurrentLabel
                    fprintf(x.S.FID,'%s\t', x.S.AtlasNames{x.S.ROIunits(iUN,5)}{x.S.ROIunits(iUN,1)} );    % print ROI name
        %             fprintf(x.S.FID,'%s\t', num2str( OverlapNVoxels) ); % print nVoxels
        %             fprintf(x.S.FID,'%s\t', num2str( OverlapVolume ) ); % print nVolume
                    fprintf(x.S.FID,'%s\t', num2str( x.S.ROIunits(iUN,4)) );% print overlapPercentage
            end


            nOutsideVoxels              = sum(x.S.TrackOverlapFound(:));
            ExplainedPercentage         = (1-nOutsideVoxels/nVoxelsLabel)*100;

            fprintf(x.S.FID,'%s\t', [num2str( nOutsideVoxels) ' nVoxels outside of any ROI'] ); % print number of unexplained voxels
            fprintf(x.S.FID,'%s\t', [', therefore ' num2str( ExplainedPercentage ) '% of cluster volume was found within ROIs'] ); % print number of unexplained voxels
            x.S                           = rmfield(x.S,'TrackOverlapFound');

            fprintf(x.S.FID,'\n'); % go to next line/row

            %% Collect all names
            if      ~isfield(x.S,'ROIunitsAll')
                    x.S.ROIunitsAll                                                                           = x.S.ROIunits;
            else    x.S.ROIunitsAll( size(x.S.ROIunitsAll,1)+1 : size(x.S.ROIunitsAll,1)+size(x.S.ROIunits,1),: )   = x.S.ROIunits; % if this already exists, add new ROI info to it
            end

            if  isfield(x.S,'ROIunits')
                x.S       = rmfield(x.S,'ROIunits');
            end

        else    fprintf(x.S.FID,'%s\t', 'No HO regions were found within this cluster');
        end
    end

    %% Print overview of all ROIs

    fprintf(x.S.FID,'\n\n\n%s\n', 'ROI name, nVoxels, Percentage Overlap of individual cluster (sorted on ROI size in clusters)');

    if      isfield(x.S,'ROIunitsAll') % if any HO-regions lied within any of the clusters

        x.S.ROIunitsAll   = sortrows(x.S.ROIunitsAll,-2);
        MentionedRois   = 0; % keep track of already mentioned ROIs

        for iUN=1:size(x.S.ROIunitsAll,1) % Print all ROIs
            if  isempty(find(MentionedRois==x.S.ROIunitsAll(iUN,1)))
                fprintf(x.S.FID,'%s\t', x.S.AtlasNames{x.S.ROIunitsAll(iUN,5)}{x.S.ROIunitsAll(iUN,1)} );    % print ROI name
                fprintf(x.S.FID,'%s\t',  num2str( x.S.ROIunitsAll(iUN,2)) );% print nVoxels
                fprintf(x.S.FID,'%s\n', num2str( x.S.ROIunitsAll(iUN,4)) );% print overlapPercentage

                MentionedRois(end+1)    = x.S.ROIunitsAll(iUN,1);
            end
        end

    else    fprintf(x.S.FID,'%s\t', 'No HO regions were found in any of the clusters');
    end

    fclose(x.S.FID); % close csv-file

end





function [x]   = PrintXYZcoordinates(x,X,Y,Z)

    fprintf(x.S.FID,'%s\t', num2str( X ) ); % print X-voxel index
    fprintf(x.S.FID,'%s\t', num2str( Y ) ); % print Y-voxel index
    fprintf(x.S.FID,'%s\t', num2str( Z ) ); % print Z-voxel index
end




function [x]   = CollectRegionsInfo(x,AtlasPath,CurrentLabel,nVoxelsLabel)
% Prints ROI name, nVoxels, nVolume, Percentage Overlap

    AtlasNifti                  = xASL_io_ReadNifti( AtlasPath);

    for iM=1:length(x.S.AtlasROIs)
        if ~isempty(x.S.AtlasROIs{iM}) % regions to skip
            clear CurrentMask OverlapMask OverlapNVoxels
            CurrentMask     = int32(AtlasNifti.dat(:,:,:))==ceil(iM/2); % for left/right

            if      round(iM/2)~=iM/2 %  odd == left
                    CurrentMask( 1: 60,:,:,:)   = 0; % exclude the right part from mask
            elseif  round(iM/2)==iM/2 % even == right
                    CurrentMask(61:end,:,:,:)   = 0; % exclude the left part from mask
            end

            OverlapMask         = CurrentMask & CurrentLabel;
            OverlapNVoxels      = sum(OverlapMask(:));

            if  OverlapNVoxels>1

                x.S.ROIunits(x.S.iNext,1)   = iM; % index for sorting later
                x.S.ROIunits(x.S.iNext,2)   = OverlapNVoxels;
                x.S.ROIunits(x.S.iNext,3)   = OverlapNVoxels*1.5^3/100^3; % OverlapVolume
                x.S.ROIunits(x.S.iNext,4)   = 100*(OverlapNVoxels/nVoxelsLabel); % OverlapPerc
                x.S.ROIunits(x.S.iNext,5)   = x.S.iA; % which AtlasNifti
                x.S.TrackOverlapFound(OverlapMask)  = 0;
                x.S.iNext                 = x.S.iNext+1;
            end
        end
    end

end
