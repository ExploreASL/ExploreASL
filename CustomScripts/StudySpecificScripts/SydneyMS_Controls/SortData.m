%% Sort Sydney Control data

% 1) rename ASL_Controls into "raw"
RootFolder = '/Users/henk/ExploreASL/SydneyMS_Controls';

%% we use these:
% This will search for T1.nii that are a file
% 	case 'SydneyMS_Controls'
% 		imPar.folderHierarchy = {'^(\d{3}.*\d{1})$' '^(T1)\.nii$'};
% 		imPar.tokenOrdering = [1 0 2];
% 		imPar.tokenSessionAliases = {};
% 		imPar.tokenScanAliases = {'^asl$', 'ASL4D'; '^T1$', 'T1'};
% 		imPar.bMatchDirectories = false;     

imPar = ExploreASL_ImportConfig(RootFolder);
ExploreASL_Import(imPar);

%% Then this one
% This will search for ASL folders containing DICOMs
% 	case 'SydneyMS_Controls'
% 		imPar.folderHierarchy = {'^(\d{3}.*\d{1})$' '^(ASL)$'};
% 		imPar.tokenOrdering = [1 0 2];
% 		imPar.tokenSessionAliases = {};
% 		imPar.tokenScanAliases = {'^ASL$', 'ASL4D'; '^T1$', 'T1'};
% 		imPar.bMatchDirectories = true;   


imPar = ExploreASL_ImportConfig(RootFolder);
ExploreASL_Import(imPar);

% Then you use the appropriate DataParFile, that contains the info about
% the Control image being the second in the ASL scan
