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

The left-right ROIs are averaged as ExploreASL divides ROIs in left-right-bilateral automatically,
except for ROIs 113-120 (vermis) and 169-170 (raphe medial and dorsal) which are bilateral already).
Hence:
  1:112 -> 1:56
113:120 -> 57:64
121:168 -> 65:88
169:170 -> 89:90

