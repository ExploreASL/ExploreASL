function loggingTable = xASL_test_Flavors_DCM2BIDS(testConfig, x, loggingTable)
%xASL_test_Flavors_DCM2BIDS Convert ASL flavors from DICOM to BIDS
%
% FORMAT: loggingTable = xASL_test_Flavors_DCM2BIDS(testConfig, x, loggingTable)
%
% INPUT:
%   testConfig         - Struct which contains the paths to the ExploreASL
%                        installation and the testing/flavor repository (REQUIRED)
%   x                  - ExploreASL x struct (STRUCT, OPTIONAL)
%   loggingTable       - Collect errors of import/processing in this table (REQUIRED)
%
% OUTPUT:              
%   General output     - Outputs the converted data and comparison results are printed on screen
%   loggingTable       - Collect errors of import/processing in this table
%         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% Runs the DICOM to ASL-BIDS import for all data in the baseDirImport directory. Study directories are supposed to be in, containing a 'sourcedata' folder - this folder
% can contain subject directories and also sourceStructure.json and studyPar.json specifying the directory structure and the additional study parameters, respectively.
% The import creates first the 'temp' subfolder with data after dcm2nii and with all tags read and saved to JSON. Then it assembles everything with the
% studyParameters and makes sure all is in BIDS format and saves it correctly in the 'rawdata' subdirectory.
%
% This function runs the following sections:
% 1.  Initialization
% 2. DICOM -> NII+JSON (i.e. dcm2niiX)
% 3. Manual curation for certain flavors
% 3a. Siemens_PCASL_3DGRASE_VD13D_2
% 3b. Philips_PCASL_3DGRASE_5.4.1.0_TopUp_1
% 3c. Siemens_PCASL_3DGRASE_VB17A_TopUp_1
% 3d. Siemens_PCASL_3DGRASE_VB17A_multiPLD_1
% 4. Convert NII+JSON -> BIDS
%
% EXAMPLE: loggingTable = xASL_test_Flavors_DCM2BIDS(testConfig, x, loggingTable);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% 1. Initialization
    if nargin < 2 || isempty(x)
        ExploreASL;
    end

    % Iterate over flavors
    for iFlavor = 1:length(testConfig.flavorList)
        
        % 2. DICOM -> NII+JSON (i.e. dcm2niiX)
        xFlavor = ExploreASL(fullfile(testConfig.pathFlavorDatabase, testConfig.flavorList{iFlavor}), [1 0 0], 0, 0);
        if isfield(xFlavor,'logging')
            loggingTable = xASL_test_AddLoggingEntryToTable(testConfig.flavorList{iFlavor},loggingTable,xFlavor.logging);
        end
        
        % 3. Manual curation for certain flavors
        xASL_test_Flavors_ManualFlavors(testConfig.flavorList, testConfig.pathFlavorDatabase, iFlavor);
        
        % 4. Convert NII+JSON -> BIDS
        xFlavor = ExploreASL(fullfile(testConfig.pathFlavorDatabase, testConfig.flavorList{iFlavor}), [0 1], 0, 0);
        if isfield(xFlavor,'logging')
            loggingTable = xASL_test_AddLoggingEntryToTable(testConfig.flavorList{iFlavor},loggingTable,xFlavor.logging);
        end

    end

end


%% Manual curation for certain flavors
function xASL_test_Flavors_ManualFlavors(flavorList, baseDirImport, iFlavor)

DirASL = fullfile(baseDirImport, flavorList{iFlavor}, 'derivatives', 'ExploreASL', 'temp', 'Sub1', 'ASL_1');

    switch flavorList{iFlavor}

             % 3a. 'Siemens_PCASL_3DGRASE_VD13D_2'
        case 'Siemens_PCASL_3DGRASE_VD13D_2'

            xASL_adm_DeleteFileList(DirASL, '^ASL4D_(32|33|34|35)_00001.*$', 1);

            nii_files = xASL_adm_GetFileList(DirASL, '^.*\.nii$', 'FPList', [], false);
            xASL_bids_MergeNifti(nii_files, 'ASL');

            % 3b. 'Philips_PCASL_3DGRASE_5.4.1.0_TopUp_1'
        case 'Philips_PCASL_3DGRASE_5.4.1.0_TopUp_1'

            % xASL_Move(fullfile(DirASL, 'M0_601_00601.nii'), fullfile(DirASL, 'M0.nii'), 1);
            % xASL_Move(fullfile(DirASL, 'M0_601_00601.json'), fullfile(DirASL, 'M0.json'), 1);
            % xASL_Move(fullfile(DirASL, 'M0_701_00701.nii'), fullfile(DirASL, 'M0PERev.nii'), 1);
            % xASL_Move(fullfile(DirASL, 'M0_701_00701.json'), fullfile(DirASL, 'M0PERev.json'), 1);
            xASL_Move(fullfile(DirASL, 'M0_601_00001.nii'), fullfile(DirASL, 'M0.nii'), 1);
            xASL_Move(fullfile(DirASL, 'M0_601_00001.json'), fullfile(DirASL, 'M0.json'), 1);
            xASL_Move(fullfile(DirASL, 'M0_701_00001.nii'), fullfile(DirASL, 'M0PERev.nii'), 1);
            xASL_Move(fullfile(DirASL, 'M0_701_00001.json'), fullfile(DirASL, 'M0PERev.json'), 1);

            % 3c. 'Siemens_PCASL_3DGRASE_VB17A_TopUp_1'
        case 'Siemens_PCASL_3DGRASE_VB17A_TopUp_1'
            imNS = xASL_io_Nifti2Im(fullfile(DirASL, 'ASL4D_NS.nii'));
            imSS = xASL_io_Nifti2Im(fullfile(DirASL, 'ASL4D_SS.nii'));
            imNS(:,:,:,2) = imSS;
            xASL_io_SaveNifti(fullfile(DirASL, 'ASL4D_NS.nii'), fullfile(DirASL, 'ASL4D.nii'), imNS/10, [], 1);

            xASL_delete(fullfile(DirASL, 'ASL4D_NS.nii'));
            xASL_delete(fullfile(DirASL, 'ASL4D_SS.nii'));

            xASL_delete(fullfile(DirASL, 'ASL4D_NS.json'));
            xASL_Move(fullfile(DirASL, 'ASL4D_SS.json'), fullfile(DirASL, 'ASL4D.json'), 1);

            xASL_Move(fullfile(DirASL, 'M0_2.json'), fullfile(DirASL, 'M0PERev.json'), 1);
            xASL_Move(fullfile(DirASL, 'M0_2.nii'), fullfile(DirASL, 'M0PERev.nii'), 1);

            % 3d. 'Siemens_PCASL_3DGRASE_VB17A_multiPLD_1'
        case 'Siemens_PCASL_3DGRASE_VB17A_multiPLD_1'

            if xASL_exist(fullfile(DirASL, 'ASL4D_NS_300.nii'))
                mTIvector = [300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000];
                for iTI = 1:length(mTIvector)
                    if iTI>1
                        xASL_delete(fullfile(DirASL, ['ASL4D_NS_' xASL_num2str(mTIvector(iTI)) '.json']));
                        xASL_delete(fullfile(DirASL, ['ASL4D_SS_' xASL_num2str(mTIvector(iTI)) '.json']));
                        imNSSS(:,:,:,2*(iTI-1)+1) = xASL_io_Nifti2Im(fullfile(DirASL, ['ASL4D_NS_' xASL_num2str(mTIvector(iTI)) '.nii']));
                        imNSSS(:,:,:,2*(iTI-1)+2) = xASL_io_Nifti2Im(fullfile(DirASL, ['ASL4D_SS_' xASL_num2str(mTIvector(iTI)) '.nii']));
                    else
                        xASL_Move(fullfile(DirASL, ['ASL4D_NS_' xASL_num2str(mTIvector(iTI)) '.json']), fullfile(DirASL, 'ASL4D.json'));
                        xASL_delete(fullfile(DirASL, ['ASL4D_SS_' xASL_num2str(mTIvector(iTI)) '.json']));
                        imNSSS = xASL_io_Nifti2Im(fullfile(DirASL, ['ASL4D_NS_' xASL_num2str(mTIvector(iTI)) '.nii']));
                        imNSSS(:,:,:,2) = xASL_io_Nifti2Im(fullfile(DirASL, ['ASL4D_SS_' xASL_num2str(mTIvector(iTI)) '.nii']));
                    end
                end
                xASL_io_SaveNifti(fullfile(DirASL, ['ASL4D_NS_' xASL_num2str(mTIvector(1)) '.nii']),...
                    fullfile(DirASL, 'ASL4D.nii'), imNSSS/10, [], 1, []);
                for iTI = 1:length(mTIvector)
                    xASL_delete(fullfile(DirASL, ['ASL4D_NS_' xASL_num2str(mTIvector(iTI)) '.nii']));
                    xASL_delete(fullfile(DirASL, ['ASL4D_SS_' xASL_num2str(mTIvector(iTI)) '.nii']));
                end
            end
    end


end




