function xASL_bids_PrepareDRO(droTestPatient,droSubject)
%xASL_bids_PrepareDRO Prepare DRO test patient for BIDS2Legacy conversion.
%
% FORMAT: xASL_bids_PrepareDRO(droTestPatient,droSubject)
% 
% INPUT:
%   droTestPatient      - Path to the DRO (CHAR ARRAY, REQUIRED)
%   droSubject          - Subject name (CHAR ARRAY, REQUIRED)
%
% OUTPUT:
%   n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Prepare DRO test patient for BIDS2Legacy conversion.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      droSubject = 'sub-Sub1';
%               droTestPatientSource = fullfile(TestRepository,'UnitTesting','dro_files','test_patient');
%               droTestPatient = fullfile(TestRepository,'UnitTesting','working_directory','test_patient');
%               xASL_Copy(droTestPatientSource,fullfile(droTestPatient,'rawdata',droSubject));
%               xASL_bids_PrepareDRO(droTestPatient,droSubject);
%               
% __________________________________
% Copyright 2015-2021 ExploreASL

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
    jsonASL.M0Type = "Separate";
    spm_jsonwrite(fullfile(perfDirectory,[droSubject,'_asl.json']),jsonASL);

    % Read M0 JSON file and add the IntendedFor field
    jsonM0 = spm_jsonread(fullfile(perfDirectory,[droSubject,'_m0scan.json']));
    jsonM0.IntendedFor = "perf/sub-Sub1_asl.nii.gz";
    spm_jsonwrite(fullfile(perfDirectory,[droSubject,'_m0scan.json']),jsonM0);

    nameDRO = "DRO_Digital_Reference_Object";
    [jsonDescription] = createDatasetDescriptionTemplate(nameDRO);

    % Write file
    spm_jsonwrite(fullfile(droTestPatient,'rawdata','dataset_description.json'),jsonDescription);


end

%% Create dataset_description.json template
function [json] = createDatasetDescriptionTemplate(name)

    % Create dummy dataset_description.json
    json = struct;
    if nargin > 0
        json.Name = name;
    else
        json.Name = "RandomText";
    end
    json.BIDSVersion = "1.5.0";
    json.DatasetType = "raw";
    json.License = "RandomText";
    json.Authors = "RandomText";
    json.Acknowledgements = "RandomText";
    json.HowToAcknowledge = "Please cite this paper: https://www.ncbi.nlm.nih.gov/pubmed/001012092119281";
    json.Funding = "RandomText";
    json.EthicsApprovals = "RandomText";
    json.ReferencesAndLinks = "RandomText";
    json.DatasetDOI = "RandomText";

end


