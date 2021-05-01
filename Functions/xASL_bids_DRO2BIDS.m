function xASL_bids_DRO2BIDS(droTestPatient,droSubject,deleteGroundTruth)
%xASL_bids_DRO2BIDS Prepare DRO test patient for BIDS2Legacy conversion.
%
% FORMAT: xASL_bids_DRO2BIDS(droTestPatient,[droSubject])
% 
% INPUT:
%   droTestPatient      - Path to the DRO (CHAR ARRAY, REQUIRED)
%   droSubject          - Subject name (CHAR ARRAY, OPTIONAL, DEFAULT = 'sub-Sub1')
%   deleteGroundTruth   - Delete DRO ground truth (BOOLEAN, OPTIONAL, DEFAULT = True)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Prepare DRO test patient for BIDS2RAW conversion.
%               This script uses the output of the asldro python script and
%               converts it into a bids structure that can be read by our
%               xASL_bids_BIDS2Legacy script.
%               An exemplary usage is shown in the unit test called
%               xASL_ut_UnitTest_function_BIDS2Legacy.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      droSubject = 'sub-Sub1';
%               droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient');
%               droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient');
%               xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject));
%               xASL_bids_DRO2BIDS(droTestPatient,droSubject);
%               
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% DRO 2 BIDS

    % Check input arguments
    if nargin<2
        droSubject = 'sub-Sub1'; % DEFAULT
    end
    if nargin<3
        deleteGroundTruth = true; % DEFAULT
    end

    % Directory definitions
    perfDirectory = fullfile(droTestPatient,'rawdata',droSubject,'perf');
    groundTruthDirectory = fullfile(droTestPatient,'rawdata',droSubject,'ground_truth');
    anatDirectory = fullfile(droTestPatient,'rawdata',droSubject,'anat');

    % Rename asl to perf
    xASL_Move(fullfile(droTestPatient,'rawdata',droSubject,'asl'),perfDirectory,1);

    % Rename ASL
    xASL_Move(fullfile(perfDirectory,'001_asl.json'),fullfile(perfDirectory,[droSubject,'_asl.json']),1);
    xASL_Move(fullfile(perfDirectory,'001_asl.nii.gz'),fullfile(perfDirectory,[droSubject,'_asl.nii.gz']),1);
    xASL_Move(fullfile(perfDirectory,'001_aslcontext.tsv'),fullfile(perfDirectory,[droSubject,'_aslcontext.tsv']),1);

    % Rename T1w
    xASL_Move(fullfile(anatDirectory,'003_anat.json'),fullfile(anatDirectory,[droSubject,'_T1w.json']),1);
    xASL_Move(fullfile(anatDirectory,'003_anat.nii.gz'),fullfile(anatDirectory,[droSubject,'_T1w.nii.gz']),1);

    % Remove the ground truth files
    if deleteGroundTruth
        xASL_delete(groundTruthDirectory,true);
    end
    
    % Read ASL JSON file
    jsonASL = spm_jsonread(fullfile(perfDirectory,[droSubject,'_asl.json']));
    
    %% Switch software versions
    softwareVersion = jsonASL.DROSoftwareVersion;
    
    switch softwareVersion
        case '2.2.0'
            %% sub-Sub1_asl.json
            
            % Change field name from LabelingType to ArterialSpinLabelingType
            jsonASL.ArterialSpinLabelingType = jsonASL.LabelingType;
            jsonASL = rmfield(jsonASL,'LabelingType');

            % Change field name from RepetitionTime to RepetitionTimePreparation
            jsonASL.RepetitionTimePreparation = jsonASL.RepetitionTime;
            jsonASL = rmfield(jsonASL,'RepetitionTime');

            % Change field name from MrAcquisitionType to MRAcquisitionType
            jsonASL.MRAcquisitionType = jsonASL.MrAcquisitionType;
            jsonASL = rmfield(jsonASL,'MrAcquisitionType');

            % Add necessary fields
            jsonASL.Manufacturer = 'Philips'; % ExploreASL expects a vendor, so we use GE as our dummy vendor here
            jsonASL.M0Type = 'Included';
            jsonASL.BackgroundSuppression = false;
            jsonASL.TotalAcquiredPairs = 1;

            % Write JSON file
            spm_jsonwrite(fullfile(perfDirectory,[droSubject,'_asl.json']),jsonASL);

            %% dataset_description.json

            % Define required fields
            jsonTemplate.Name = 'DRO_Digital_Reference_Object';

            % Call script to fix missing fields
            [jsonDescription] = xASL_bids_CreateDatasetDescriptionTemplate(jsonTemplate);

            % Write dataset_description file
            spm_jsonwrite(fullfile(droTestPatient,'rawdata','dataset_description.json'),jsonDescription);

            %% sourceStructure.json

            % Define sourceStructure
            sourceStructure.folderHierarchy = {'^(.)+$','^(ASL|T1w|M0|T2|FLAIR)$'};
            sourceStructure.tokenOrdering = [1,0,2];
            sourceStructure.tokenSessionAliases = {'',''};
            sourceStructure.tokenScanAliases = {'^ASL$','ASL4D','^T1w$','T1w','^M0$','M0','^T2$','T2w','^FLAIR$','FLAIR'};
            sourceStructure.bMatchDirectories = true;

            % Write sourceStructure file
            spm_jsonwrite(fullfile(droTestPatient,'sourceStructure.json'),sourceStructure);

            %% studyPar.json

            % Define studyPar
            studyPar.DatasetType = 'raw';
            studyPar.License = 'license';
            studyPar.Authors = {'ASPIRE Project: Gold Standard Phantoms'};
            studyPar.Acknowledgements = 'acknowledgements';
            studyPar.HowToAcknowledge = 'Gold Standard Phantoms';
            studyPar.Funding = {'ASPIRE Project'};
            studyPar.EthicsApprovals = {'ASPIRE Project'};
            studyPar.ReferencesAndLinks = {'https://github.com/gold-standard-phantoms/asldro', 'https://pypi.org/project/asldro/', 'https://asldro.readthedocs.io/'};
            studyPar.DatasetDOI = 'https://pypi.org/project/asldro/';
            studyPar.LabelingType = jsonASL.ArterialSpinLabelingType;
            studyPar.ASLContext = 'm0scan,control,label';

            % Write sourceStructure file
            spm_jsonwrite(fullfile(droTestPatient,'studyPar.json'),studyPar);
            
        otherwise
            warning('Unknown DRO version...');
    end

end



