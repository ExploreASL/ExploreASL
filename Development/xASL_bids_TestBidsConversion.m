function xASL_bids_TestBidsConversion(baseDirImport, baseDirReference, bImport, bComparison)
%xASL_bids_TestBidsConversion Convert ASL flavors from DICOM to BIDS, comparing the results with reference data
%
% FORMAT: xASL_bids_TestBidsConversion(baseDirImport[, baseDirReference, bImport, bComparison])
%
% INPUT:
%   baseDirImport    - Directory for import - sourcedata files are in the 'sourcedata' folder (REQUIRED)
%   baseDirReference - Reference directory with correct ASL-BIDS data (OPTIONAL)
%   bImport          - Specify if import should be performed (OPTIONAL, DEFAULT=TRUE)
%   bComparison      - Specify if comparison should be performed (OPTIONAL, DEFAULT=FALSE)
%
% OUTPUT: n/a        - Outputs the converted data and comparison results are printed on screen
%         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% Runs the DICOM to ASL-BIDS import for all data in the baseDirImport directory. Study directories are supposed to be in, containing a 'sourcedata' folder - this folder
% can contain subject directories and also sourceStructure.json and studyPar.json specifying the directory structure and the additional study parameters, respectively.
% The import creates first the 'analysis' subfolder with data after dcm2nii and with all tags read and saved to JSON. Then it assembles everything with the
% studyParameters and makes sure all is in BIDS format and saves it correctly in the 'rawdata' subdirectory.
%
% This function runs the following sections:
% 1.  Initialization
% 2.  Go through all flavors and import them (i.e. convert them to BIDS)
% 2a. Siemens_PCASL_3DGRASE_vascular
% 2b. Philips_PCASL_3DGRASE_R5.4_TopUp
% 2c. Siemens_PCASL_volunteer
% 2d. Siemens_PCASL_multiTI
% 2e. Flavor not needing manual adjustment
% 3.  Run the comparison with the reference directory
%
% EXAMPLE: xASL_bids_TestBidsConversion('mydir/testImport');
%          xASL_bids_TestBidsConversion('mydir/testImport', 'mydir/testReference',0,1);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% Admin
if nargin<2 || isempty(baseDirReference)
	baseDirReference = ''; % isempty() is from [], so still needs conversion
end

if nargin<3 || isempty(bImport)
	bImport = true;
end

if nargin<4 || isempty(bComparison)
	bComparison = false;
end

if bComparison && isempty(baseDirReference)
	error('Cannot compare to a reference if the reference directory missing');
end

%% 1. Initialization
% Initialize ExploreASL 
ExploreASL_Initialize([], false);

% Load the list of the directories
flavorList = xASL_adm_GetFileList(baseDirImport, [], false, [], true);

%% 2. Go through all flavors and import them (i.e. convert them to BIDS)
if bImport
	for iFlavor = 1:length(flavorList)
        
		switch(flavorList{iFlavor})
			% For certain flavors, we first run the import to NII+JSON (i.e. dcm2niiX),
            % after which the ASL files need to be curated
			% Afterwards we can continue with the NII+JSON->ASL-BIDS
            
            % 2a. 'Siemens_PCASL_3DGRASE_vascular'
			case 'Siemens_PCASL_3DGRASE_vascular'
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [1 0 0], false, true, false, false);
                DirASL = fullfile(baseDirImport, flavorList{iFlavor}, 'analysis', 'Sub1', 'ASL_1');
                
                xASL_adm_DeleteFileList(DirASL, '^ASL4D_(7|8|9|10).*$', 1);

				nii_files = xASL_adm_GetFileList(DirASL, '^.*\.nii$', 'FPList', [], false);
				nii_files = xASL_bids_MergeNifti(nii_files, 'ASL');
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [0 1 0], false, true, false, false);
				
            % 2b. 'Philips_PCASL_3DGRASE_R5.4_TopUp'
			case 'Philips_PCASL_3DGRASE_R5.4_TopUp'
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [1 0 0], false, true, false, false);
                DirASL = fullfile(baseDirImport, flavorList{iFlavor}, 'analysis', 'Sub1', 'ASL_1');
                
				xASL_Move(fullfile(DirASL, 'M0_1.nii'), fullfile(DirASL, 'M0.nii'), 1);
				xASL_Move(fullfile(DirASL, 'M0_1.json'), fullfile(DirASL, 'M0.json'), 1);
				xASL_Move(fullfile(DirASL, 'M0_2.nii'), fullfile(DirASL, 'M0PERev.nii'), 1);
				xASL_Move(fullfile(DirASL, 'M0_2.json'), fullfile(DirASL, 'M0PERev.json'), 1);
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [0 1 0], false, true, false, false);
				
            % 2c. 'Siemens_PCASL_volunteer'
			case 'Siemens_PCASL_volunteer'
				ExploreASL_ImportBIDS(fullfile(baseDirImport,flavorList{iFlavor}), [],[], [1 0 0], false, true, false, false);
                DirASL = fullfile(baseDirImport, flavorList{iFlavor}, 'analysis', 'Sub1', 'ASL_1');
                
				if xASL_exist(fullfile(DirASL, 'ASL4D_NS.nii'))
					xASL_delete(fullfile(DirASL, ASL4D_NS.json'));
					xASL_Move(fullfile(DirASL, ASL4D_SS.json'), fullfile(DirASL, ASL4D.json'), 1);
                              
					imNS = xASL_io_Nifti2Im(fullfile(DirASL, ASL4D_NS.nii'));
					imSS = xASL_io_Nifti2Im(fullfile(DirASL, ASL4D_SS.nii'));
					imNS(:,:,:,2) = imSS;
					xASL_io_SaveNifti(fullfile(DirASL, ASL4D_NS.nii'), fullfile(DirASL, ASL4D.nii'), imNS/10, [], 1);
                        
					xASL_delete(fullfile(DirASL, ASL4D_NS.nii'));
					xASL_delete(fullfile(DirASL, ASL4D_SS.nii'));
                    
					xASL_Move(fullfile(DirASL, M0_2.json'), fullfile(DirASL, M0PERev.json'), 1);
					xASL_Move(fullfile(DirASL, M0_2.nii'), fullfile(DirASL, M0PERev.nii'), 1);
				end
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [0 1 0], false, true, false, false);
				
            % 2d. 'Siemens_PCASL_multiTI'
			case 'Siemens_PCASL_multiTI'
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [1 0 0], false, true, false, false);
                DirASL = fullfile(baseDirImport, flavorList{iFlavor}, 'analysis', 'Sub1', 'ASL_1');
                
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
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [0 1 0], false, true, false, false);
				
            % 2e. 'Flavor not needing manual adjustment'
			otherwise
				ExploreASL_ImportBIDS(fullfile(baseDirImport, flavorList{iFlavor}), [],[], [1 1 0], false, true, false, false);
		end
	end
end

%% 3. Run the comparison with the reference directory
if bComparison
	% List all studies in the import directory
	filenameCompare = xASL_adm_GetFileList(baseDirImport, '^.+$', 'List', [], true);
	for iCompare = 1:length(filenameCompare)
		% Compare the imported data in the 'rawdata' subdirectory with the counterpart
		fprintf('%s\n', ['Dataset: '  filenameCompare{iCompare}]);
		xASL_bids_CompareStructures(fullfile(baseDirImport, filenameCompare{iCompare}, 'rawdata'),...
            fullfile(baseDirReference, filenameCompare{iCompare}, 'rawdata'));
	end
end


end