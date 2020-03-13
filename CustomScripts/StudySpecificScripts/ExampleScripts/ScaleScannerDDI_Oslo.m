ROOT    = 'D:\DementiaASL_data_Atle\DDI_IVS\analysis';

load( fullfile(ROOT, 'scanner.mat') );

GM_mask     = 'D:\DementiaASL_data_Atle\DDI_IVS\analysis\dartel\DARTEL_T1_template.nii';
GM_mask     = xASL_io_ReadNifti(GM_mask);
GM_mask     = GM_mask.dat(:,:,:)>0.5;

% dip_image(GM_mask)

% Load CBF maps
Ingenia     = 1;
Achieva     = 1;

for iS=1:length(scanner)
    clear FileName tnii meanGM
    FileName    = fullfile(ROOT, 'dartel', ['qCBF_' scanner{iS,1} '_ASL_1.nii']);

    if xASL_exist(FileName,'file')

        tnii        = xASL_io_ReadNifti(FileName);
        tnii        = tnii.dat(:,:,:);
        meanGM      = mean(tnii(GM_mask & isfinite(tnii)));

        if      scanner{iS,2}==1
                meanGMall(Achieva,1)    = meanGM;
                Achieva     = Achieva + 1;
        elseif  scanner{iS,2}==2
                meanGMall(Ingenia,2)    = meanGM;
                Ingenia     = Ingenia + 1;
        else    error('Unknown option');
        end
    end
end

meanAchieva     = mean(meanGMall(:,1));
meanIngenia     = mean(meanGMall(1:26,2));
ScaleAchieva    = 50/meanAchieva;
ScaleIngenia    = 50/meanIngenia;


for iS=1:length(scanner)
    clear FileName tnii meanGM
    FileName    = fullfile(ROOT, 'dartel', ['qCBF_' scanner{iS,1} '_ASL_1.nii']);

    if  xASL_exist(FileName,'file')

        tnii        = xASL_io_ReadNifti(FileName);
        tnii        = tnii.dat(:,:,:);

        if      scanner{iS,2}==1
                tnii    = tnii.*ScaleAchieva;
        elseif  scanner{iS,2}==2
                tnii    = tnii.*ScaleIngenia;
        else    error('Unknown option');
        end

        xASL_io_SaveNifti( FileName, FileName, tnii);
    end
end
