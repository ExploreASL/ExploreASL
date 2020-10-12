% Calculate the Tanimoto coefficient between ASL and MNI directly in the MNI space
% That makes it possible to do that for both the ASL-pseudoT1w registration and also for the normal registration via segmented T1w

%studyName = 'analysisASLonlyRigidFinal';
%studyName = 'analysisASLonlyAffineFinal';
%studyName = 'analysisASLonlyDCTFinal';
studyName = 'analysisSegm1y';
DataParPath = ['/pet/projekte/asl/data/Craniosynostosis/' studyName '/Craniosynostosis.json'];

%% 1. Initialize the paths
x = ExploreASL_Initialize(DataParPath);
x.iSession = 1; % Go for session 1 only

%% 2. Go through all subjects
for iSubject = 1:x.nSubjects
	x.iSubject = iSubject;
	x.SESSIONDIR = x.SESSIONS{1};
	x.SUBJECTDIR = x.SUBJECTS{iSubject};
	x = xASL_init_FileSystem(x); % initialize FileSystem, quick & dirty
	
	if xASL_exist(fullfile(x.D.PopDir,['MaskVascular_' x.SUBJECTS{iSubject} '_ASL_1.nii']))
		xASL_TrackProgress(iSubject,x.nSubjects);
		% Get the individual mask and CBF image
		MaskFromTemplate = xASL_io_Nifti2Im(fullfile(x.D.PopDir,['MaskVascular_' x.SUBJECTS{iSubject} '_ASL_1.nii']))>0.5;
		PWIim = xASL_io_Nifti2Im(fullfile(x.D.PopDir,['qCBF_' x.SUBJECTS{iSubject} '_ASL_1.nii']));
		
		% Differs between sequences - get the CBF template
		if      strcmpi(x.readout_dim,'2D') && ~isempty(regexpi(x.Vendor,'Philips'))
			x.Mean_MNI = fullfile(x.D.TemplateDir,'Philips_2DEPI_Bsup_CBF.nii');
		elseif  strcmpi(x.readout_dim,'2D') && ~isempty(regexpi(x.Vendor,'Siemens'))
			x.Mean_MNI = fullfile(x.D.TemplateDir,'Siemens_2DEPI_PCASL_noBsup_CBF.nii');
		elseif  strcmpi(x.readout_dim,'3D')&& ~isempty(regexpi(x.Vendor,'Siemens'))
			x.Mean_MNI = fullfile(x.D.TemplateDir,'Siemens_3DGRASE_PASL_CBF.nii');
		elseif  strcmpi(x.readout_dim,'3D')&& ~isempty(regexpi(x.Vendor,'GE'))
			x.Mean_MNI = fullfile(x.D.TemplateDir,'GE_3Dspiral_Product_CBF.nii');
		else
			error('Unknown sequence/readout');
		end
		
		TemplateIm = xASL_io_Nifti2Im(x.Mean_MNI);
		
		% Compute wholebrain Tanimoto coefficient
		TanimotoCoeff = xASL_qc_TanimotoCoeff(PWIim, TemplateIm, MaskFromTemplate, 3, 0.975);
		results{x.iSubject,2} = TanimotoCoeff;
	else
		results{x.iSubject,2} = 'n/a';
	end
	results{x.iSubject,1} = x.SUBJECTS{iSubject};
end

xASL_tsvWrite(results,['/home/janpetr/tmp/TC-ASL-MNI-' studyName '.tsv'], 1);
 
