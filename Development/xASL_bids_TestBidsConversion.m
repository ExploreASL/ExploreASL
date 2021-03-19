function xASL_bids_TestBidsConversion(baseDirImport,baseDirReference,bImport,bComparison,bInitialize)
%xASL_bids_TestBidsConversion Runs the conversion for the DICOM flavors to ASL-BIDS and compares the results with reference data
%
% FORMAT: xASL_bids_TestBidsConversion(baseDirImport[,baseDirReference,bImport,bComparison])
%
% INPUT:
%   baseDirImport    - Directory for import - sourcedata files are in the 'sourcedata' folder (REQUIRED)
%   baseDirReference - Reference directory with correct ASL-BIDS data (OPTIONAL)
%   bImport          - Specify if import should be performed (OPTIONAL, DEFAULT=TRUE)
%   bComparison      - Specify if comparison should be performed (OPTIONAL, DEFAULT=FALSE)
%   bInitialize      - Run ExploreASL_Initialize (OPTIONAL, DEFAULT=TRUE)
%
% OUTPUT: n/a        - Outputs the converted data and comparison results are printed on screen
%         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% Runs the DICOM to ASL-BIDS import for all data in the baseDirImport directory. Study directories are supposed to be in, containing a 'sourcedata' folder - this folder
% can contain subject directories and also sourceStructure.json and studyPar.json specifying the directory structure and the additional study parameters, respectively.
% The import creates first the 'analysis' subfolder with data after dcm2nii and with all tags read and saved to JSON. Then it assembles everything with the
% studyParameters and makes sure all is in BIDS format and saves it correctly in the 'rawdata' subdirectory.
% EXAMPLE: xASL_bids_TestBidsConversion('mydir/testImport');
%          xASL_bids_TestBidsConversion('mydir/testImport','mydir/testReference',0,1);
% __________________________________
% Copyright 2015-2020 ExploreASL

%% Admin
if nargin<2 
	baseDirReference = '';
end

if nargin<3 || isempty(bImport)
	bImport = true;
end

if nargin<4 || isempty(bComparison)
	bComparison = false;
end

if nargin<4 || isempty(bInitialize)
	bInitialize = true;
end

if isempty(baseDirReference) && bComparison
	error('Cannot compare to reference if the reference directory is not specified.');
end

%% Initialization
% Initialize ExploreASL 
if bInitialize
    ExploreASL_Initialize([],false);
end

% Load the list of the directories
fList = xASL_adm_GetFileList(baseDirImport,[],false,[],true);

%% Go through all studies and import them
if bImport
	for ii = 1:length(fList)
		% Import the whole session to JSON and NIFTI
		switch(fList{ii})
			% For certain session, we first run the import to NII+JSON, then the ASL files need to be curated and
			% Afterwards we can continue with the NII+JSON->ASL-BIDS
			case 'Siemens_PCASL_3DGRASE_vascular'
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [1 0 0], false, true, false, false);
				system(['rm ' baseDirImport '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/ASL4D_7*']);
				system(['rm ' baseDirImport '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/ASL4D_8*']);
				system(['rm ' baseDirImport '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/ASL4D_9*']);
				system(['rm ' baseDirImport '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/ASL4D_10*']);
				nii_files = xASL_adm_GetFileList([baseDirImport '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/'],'^*.nii$','FPList',[],false);
				nii_files = xASL_bids_MergeNifti(nii_files, 'ASL');
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [0 1 0], false, true, false, false);
				
			case 'Philips_PCASL_3DGRASE_R5.4_TopUp'
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [1 0 0], false, true, false, false);
				system(['mv ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_1.nii ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0.nii']);
				system(['mv ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_1.json ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0.json']);
				system(['mv ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_2.nii ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0PERev.nii']);
				system(['mv ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_2.json ' baseDirImport '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0PERev.json']);
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [0 1 0], false, true, false, false);
				
			case 'Siemens_PCASL_volunteer'
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [1 0 0], false, true, false, false);
				if xASL_exist([baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii'])
					system(['rm ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.json']);
					system(['mv ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_SS.json ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D.json']);
					imNS = xASL_io_Nifti2Im([baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii']);
					imSS = xASL_io_Nifti2Im([baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_SS.nii']);
					imNS(:,:,:,2) = imSS;
					xASL_io_SaveNifti([baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii'],[baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D.nii'],imNS/10,[],1,[]);
					system(['rm ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii']);
					system(['rm ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_SS.nii']);
					system(['mv ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0_2.json ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0PERev.json']);
					system(['mv ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0_2.nii ' baseDirImport '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0PERev.nii']);
				end
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [0 1 0], false, true, false, false);
				
			case 'Siemens_PCASL_multiTI'
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [1 0 0], false, true, false, false);
				if xASL_exist([baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_300.nii'])
					mTIvec = [300,600,900,1200,1500,1800,2100,2400,2700,3000];
					for jj = 1:length(mTIvec)
						if jj>1
							system(['rm ' baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' xASL_num2str(mTIvec(jj)) '.json']);
							system(['rm ' baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' xASL_num2str(mTIvec(jj)) '.json']);
							imNSSS(:,:,:,2*(jj-1)+1) = xASL_io_Nifti2Im([baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' xASL_num2str(mTIvec(jj)) '.nii']);
							imNSSS(:,:,:,2*(jj-1)+2) = xASL_io_Nifti2Im([baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' xASL_num2str(mTIvec(jj)) '.nii']);
						else
							system(['mv ' baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' xASL_num2str(mTIvec(jj)) '.json ' baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D.json']);
							system(['rm ' baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' xASL_num2str(mTIvec(jj)) '.json']);
							imNSSS = xASL_io_Nifti2Im([baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' xASL_num2str(mTIvec(jj)) '.nii']);
							imNSSS(:,:,:,2) = xASL_io_Nifti2Im([baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' xASL_num2str(mTIvec(jj)) '.nii']);
						end
					end
					xASL_io_SaveNifti([baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' xASL_num2str(mTIvec(1)) '.nii'],[baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D.nii'],imNSSS/10,[],1,[]);
					for jj = 1:length(mTIvec)
						system(['rm ' baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' xASL_num2str(mTIvec(jj)) '.nii']);
						system(['rm ' baseDirImport '/Siemens_PCASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' xASL_num2str(mTIvec(jj)) '.nii']);
					end
				end
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [0 1 0], false, true, false, false);
				
			otherwise
				% But in normal case, no special manipulation is needed
				ExploreASL_ImportBIDS(fullfile(baseDirImport,fList{ii}), [],[], [1 1 0], false, true, false, false);
		end
	end
end

%% Runs the comparison with the reference directory
if bComparison
	% List all studies in the import directory
	filenameCompare = xASL_adm_GetFileList(baseDirImport,'^.+$','List',[],true);
	for iCompare = 1:length(filenameCompare)
		% Compare the imported data in the 'rawdata' subdirectory with the counterpart
		display(['Dataset: '  filenameCompare{iCompare}]);
		xASL_bids_CompareStructures(fullfile(baseDirImport, filenameCompare{iCompare},'rawdata'),fullfile(baseDirReference,filenameCompare{iCompare},'rawdata'));
	end
end

end