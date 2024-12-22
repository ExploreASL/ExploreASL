IMpath = '/Users/hjmutsaerts/ExploreASL/ExploreASL/External/Atlases4ROIs/LicensePermissive/AAL3v1.nii';
IM = xASL_io_Nifti2Im(IMpath);
TSVpath = '/Users/hjmutsaerts/ExploreASL/ExploreASL/External/Atlases4ROIs/LicensePermissive/AAL3v1.tsv';
TSV = xASL_tsvRead(TSVpath);
TSV = TSV(:,2);

checkNR = unique(IM(:));

for iROI = 1:56
    IM(IM==(iROI*2)-1) = iROI;
    IM(IM==iROI*2) = iROI;
    TSV{iROI} = TSV{iROI*2}(1:end-2);
end
for iROI = 113:120
    IM(IM==iROI) = iROI-56;
    TSV{iROI-56} = TSV{iROI};
end
for iROI = (121:144)-120
    IM(IM==(iROI*2)-1+120) = iROI+64;
    IM(IM==iROI*2+120) = iROI+64;
    TSV{iROI+64} = TSV{iROI*2+120}(1:end-2);
end
for iROI = 169:170
    IM(IM==iROI) = iROI-80;
    TSV{iROI-80} = TSV{iROI};
end

TSV = TSV(1:90);

xASL_io_SaveNifti(IMpath, IMpath, IM, []);
xASL_tsvWrite(TSV', TSVpath, 1);
IM = uint8(IM);
save([IMpath '.mat'],'IM');
gzip(IMpath);
delete(IMpath);

The left-right ROIs are averaged as ExploreASL divides ROIs in left-right-bilateral automatically,
except for ROIs 113-120 (vermis) and 169-170 (raphe medial and dorsal) which are bilateral already).
Hence:
  1:112 -> 1:56
113:120 -> 57:64
121:168 -> 65:88
169:170 -> 89:90

%% The same for WMPM

clear TSVnew

IMpath = '/Users/hjmutsaerts/ExploreASL/ExploreASL/External/Atlases4ROIs/LicenseLimited/WMPM_Type_III.nii';
IM = xASL_io_Nifti2Im(IMpath);
TSVpath = '/Users/hjmutsaerts/ExploreASL/ExploreASL/External/Atlases4ROIs/LicenseLimited/WMPM_Type_III.tsv';
TSV = xASL_tsvRead(TSVpath);

IMnew = zeros(size(IM));

% Remove double _
TSV = strrep(TSV, '__', '_');

% Get left or right
LeftIs = contains(TSV(:,1), 'left');
RightIs = contains(TSV(:,1), 'right');
TSV(LeftIs,2) = {1}; % left
TSV(RightIs,2) = {2}; % right

% Remove trailing left/right
TSV(:,1) = strrep(TSV(:,1), '_left', '');
TSV(:,1) = strrep(TSV(:,1), 'left', '');
TSV(:,1) = strrep(TSV(:,1), '_right', '');
TSV(:,1) = strrep(TSV(:,1), 'right', '');

% Bugfix
TSV{108,1} = strrep(TSV{108,1}, '_(a_part_of_MCP)', '');

iROInew = 1;
reportSingleROI = {''};
for iROI=1:length(TSV)
    % is there a ROI with a similar name
    if ~isempty(TSV{iROI,1})
        indicesAre = find(strcmp(TSV(:,1), TSV{iROI,1}));
        
        if numel(indicesAre)==1
            reportSingleROI{end+1,1} = TSV{iROI,1};
            reportSingleROI{end,2} = iROI;
        end

        % assuming they are always either left or right
        TSVnew{iROInew} = TSV{iROI,1};
        for iIndices=1:length(indicesAre)
            IMnew(IM==indicesAre(iIndices)) = iROInew;
        end
        % Remove current ROI
        TSV(indicesAre,1) = {''};
        iROInew = iROInew+1;
    end
end

xASL_io_SaveNifti(IMpath, IMpath, IMnew, []);
xASL_tsvWrite(TSVnew, TSVpath, 1);
IM = uint8(IMnew);
save([IMpath '.mat'],'IM');
gzip(IMpath);
delete(IMpath);