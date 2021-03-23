function xASL_bids_DRO2BIDS(droTestPatient,droSubject)
%xASL_bids_DRO2BIDS Prepare DRO test patient for BIDS2Legacy conversion.
%
% FORMAT: xASL_bids_DRO2BIDS(droTestPatient,[droSubject])
% 
% INPUT:
%   droTestPatient      - Path to the DRO (CHAR ARRAY, REQUIRED)
%   droSubject          - Subject name (CHAR ARRAY, OPTIONAL, DEFAULT = 'sub-Sub1')
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

    % Check input arguments
    if nargin<2
        droSubject = 'sub-Sub1'; % DEFAULT
    end

    % Directory definitions
    perfDirectory = fullfile(droTestPatient,'rawdata',droSubject,'perf');
    groundTruthDirectory = fullfile(droTestPatient,'rawdata',droSubject,'ground_truth');
    anatDirectory = fullfile(droTestPatient,'rawdata',droSubject,'anat');

    % Rename asl to perf
    xASL_Move(fullfile(droTestPatient,'rawdata',droSubject,'asl'),perfDirectory);

    % Rename ASL
    xASL_Move(fullfile(perfDirectory,'001_asl.json'),fullfile(perfDirectory,[droSubject,'_asl.json']));
    xASL_Move(fullfile(perfDirectory,'001_asl.nii.gz'),fullfile(perfDirectory,[droSubject,'_asl.nii.gz']));
    xASL_Move(fullfile(perfDirectory,'001_aslcontext.tsv'),fullfile(perfDirectory,[droSubject,'_aslcontext.tsv']));

    % Add M0 from ground_truh
    xASL_Move(fullfile(groundTruthDirectory,'002_ground_truth_m0.json'),fullfile(perfDirectory,[droSubject,'_m0scan.json']));
    xASL_Move(fullfile(groundTruthDirectory,'002_ground_truth_m0.nii.gz'),fullfile(perfDirectory,[droSubject,'_m0scan.nii.gz']));

    % Rename T1w
    xASL_Move(fullfile(anatDirectory,'003_anat.json'),fullfile(anatDirectory,[droSubject,'_T1w.json']));
    xASL_Move(fullfile(anatDirectory,'003_anat.nii.gz'),fullfile(anatDirectory,[droSubject,'_T1w.nii.gz']));

    % Remove the ground truth files
    xASL_delete(groundTruthDirectory,true);

    % Read ASL JSON file and add the M0Type field
    jsonASL = spm_jsonread(fullfile(perfDirectory,[droSubject,'_asl.json']));
    jsonASL.M0Type = 'Separate';
    spm_jsonwrite(fullfile(perfDirectory,[droSubject,'_asl.json']),jsonASL);

    % Read M0 JSON file and add the IntendedFor field
    jsonM0 = spm_jsonread(fullfile(perfDirectory,[droSubject,'_m0scan.json']));
    jsonM0.IntendedFor = 'perf/sub-Sub1_asl.nii.gz';
    spm_jsonwrite(fullfile(perfDirectory,[droSubject,'_m0scan.json']),jsonM0);

    nameDRO = 'DRO_Digital_Reference_Object';
    [jsonDescription] = xASL_bids_CreateDatasetDescriptionTemplate(nameDRO);

    % Write file
    spm_jsonwrite(fullfile(droTestPatient,'rawdata','dataset_description.json'),jsonDescription);


end



