function [x] = xASL_wrp_Load4DMemMapping_LesionsROIs(x)
%xASL_wrp_Load4DMemMapping_LesionsROIs Part of ExploreASL analysis module
% Loads data & maps it to memory mapping file on disc, if not done before

%% Admin, specify for each modality how to treat it
ModalitiesList  = {'Lesion' 'ROI'};
SumMask         = sum(x.WBmask(:));

x.StudyAtlasDir   = fullfile(x.D.PopDir,'AtlasesStudySpecific');

for iM=1:2

    % Create single 4D-file for memory mapping, if doesn't exist yet,
    % delete for reprocessing!
    LoadFileDat     = fullfile( x.StudyAtlasDir, ['StudySpecific_' ModalitiesList{iM} 's.dat']);

    %% Find available scans
    % for each subject, T1/FLAIR, get a list of scans, and concatenate
    % them, with NaNs if not existing

    MaskType    = '';

    fprintf(['Find available NIfTI ' ModalitiesList{iM} ' masks:   ']);
    for iS=1:x.nSubjects
        xASL_TrackProgress(iS,x.nSubjects);
        LoadFileList                = xASL_adm_GetFileList(x.D.PopDir, ['r' ModalitiesList{iM} '_(T1|FLAIR)_\d*_' x.SUBJECTS{iS} '\.(nii|nii\.gz)'],'FPList',[0 Inf]);
        for iL=1:length(LoadFileList)
            [Fpath Ffile Fext]      = xASL_fileparts(LoadFileList{iL});
            MaskType{end+1,1}       = Ffile(3+length(ModalitiesList{iM}):end-length(x.SUBJECTS{iS})-1);
        end
    end
    fprintf('   ');

    if  isempty(MaskType)
        fprintf('%s\n',['No subject-specific ' ModalitiesList{iM} ' masks found, skipping...']);
    else
        xASL_adm_CreateDir(x.StudyAtlasDir);

        MaskType                        = unique(MaskType); % collect a list of Lesions/ROIs specified for all subjects
        IMtype                          = {'IntraMask' 'PeriMask' 'Added (Intra+Peri) Mask' 'Contralateral added Mask' 'Ipsilateral hemisphere - added Mask' 'Contralateral hemisphere Mask'};
        for iT=1:length(MaskType)
            for iIM=1:length(IMtype)
                iTiIM   = (iT-1)*length(IMtype)+iIM;
                CSVdata{iTiIM}    = [MaskType{iT} '_' IMtype{iIM}];
            end
        end

        %% Create csv-file
        xASL_adm_CreateCSVfile([LoadFileDat(1:end-4) '.csv'],CSVdata);

        nMasks                          = length(MaskType); % assume 6 masks
        nIMs                            = 6;

       if ~exist( LoadFileDat ,'file' )
            %% Start memory mapping

            tempData                    = zeros(SumMask,nIMs*nMasks,x.nSubjects,'uint8'); % pre-allocation
            fprintf(['Memory mapping ' ModalitiesList{iM} ' masks:   ']);
            for iS=1:x.nSubjects
                xASL_TrackProgress(iS,x.nSubjects);
                for iT=1:length(MaskType)
                    LoadFile                = fullfile(x.D.PopDir, ['r' ModalitiesList{iM} '_' MaskType{iT} '_' x.SUBJECTS{iS} '.nii']);

                    if  xASL_exist(LoadFile,'file')
                        tNII = xASL_im_IM2Column(uint8(xASL_io_Nifti2Im( LoadFile)), x.WBmask, false);

                        if  size(tNII,2)<nIMs % ROI mask creation not happened yet, redo
                            % CAVE: we assume 6 masks here
                            xASL_im_Lesion2Mask( LoadFile, [], [], [], x );
                            tNII = xASL_im_IM2Column(uint8(xASL_io_Nifti2Im( LoadFile)), x.WBmask, false);
                        end

                        % Load masks
                        tempData(:,(iT-1)*nIMs+1:(iT-1)*nIMs+nIMs,iS)   = tNII;
                    end
                end
            end

            fileID = fopen(LoadFileDat,'w');
            fwrite( fileID,tempData,'uint8');
            fclose(fileID);
        end % if ~exist( LoadFile ,'file' )
    end  % if  isempty(MaskType)
end % for iM=1:2


end
