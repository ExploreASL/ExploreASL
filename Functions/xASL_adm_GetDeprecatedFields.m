function nameConversionTable = xASL_adm_GetDeprecatedFields()
%xASL_adm_GetDeprecatedFields Get the names of all deprecated fields and
%their corresponding new field names
%
% FORMAT: nameConversionTable = xASL_adm_GetDeprecatedFields()
%
% INPUT:
%   n/a
%
% OUTPUT:
%   nameConversionTable - Conversion table
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script is mainly used to improve backwards
%              compatibility. Check the usage in both xASL_adm_LoadX 
%              and xASL_io_ReadDataPar.
%
% EXAMPLE:     n/a
%
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

nameConversionTable = ...
	{'nSessions' 					'dataset' 	'nSessions'                 '';...
	'nSubjectsSessions' 			'dataset' 	'nSubjectsSessions'         '';...
	'nTotalSubjects' 				'dataset' 	'nTotalSubjects'            '';...
	'TotalSubjects' 				'dataset' 	'TotalSubjects'             '';...
	'nTimePoints' 					'dataset' 	'nTimePoints'               '';...
	'exclusion' 					'dataset' 	'exclusion'                 '';...
	'TimePointSubjects' 			'dataset' 	'TimePointSubjects'         '';...
	'nTimePointSubjects' 			'dataset' 	'nTimePointSubjects'        '';...
	'ExcludedSubjects' 				'dataset' 	'ExcludedSubjects'          '';...
	'nTimePointExcluded' 			'dataset' 	'nTimePointExcluded'        '';...
	'TimePointTotalSubjects' 		'dataset' 	'TimePointTotalSubjects'    '';...
	'TimePointExcluded' 			'dataset' 	'TimePointExcluded'         '';...
	'nTimePointTotalSubjects' 		'dataset' 	'nTimePointTotalSubjects'   '';...
	'nTimePointsTotal' 				'dataset' 	'nTimePointsTotal'          '';...
	'nExcluded' 					'dataset' 	'nExcluded'                 '';...
	'TotalInclusionList' 			'dataset' 	'TotalInclusionList'        '';...
    'name'                          'dataset'   'name'                      '';...
    'Sequence'                      'Q'         'Sequence'                  '';...
    'Vendor'                        'Q'         'Vendor'                    '';...
    'readout_dim'       			'Q'         'readoutDim'                '';...
	'SpaghettiDir' 					'D' 		'SpaghettiDir'              '';...
	'HistogramDir' 					'D' 		'HistogramDir'              '';...
	'StatsMaps'    					'D' 		'StatsMaps'                 '';...
	'SPMDIR'       					'D' 		'SPMDIR'                    '';...
	'SPMpath'      					'D' 		'SPMpath'                   '';...
	'LockDir'      					'dir' 		'LockDir'                   '';...
	'dataParType'  					'opts' 		'dataParType'               '';...
	'SPMVERSION'   					'external' 	'SPMVERSION'                '';...
	'M0'                            'settings' 	'M0'                        '';...
    'subject_regexp'       			'settings' 	'subjectRegexp'             '';...
    'Quality'                       'settings' 	'Q'                         '';...
    'BILAT_FILTER'       			'settings' 	'BILAT_FILTER'              '';...
	'Pediatric_Template' 			'settings' 	'Pediatric_Template'        '';...
	'bReproTesting'      			'settings' 	'bReproTesting'             '';...
	'bOverwrite'         			'settings' 	'bOverwrite'                '';...
	'dryRun'             			'settings' 	'dryRun'                    '';...
	'M0_conventionalProcessing' 	'settings' 	'M0_conventionalProcessing' '';...
	'bAutoACPC'          			'settings' 	'bAutoACPC'                 '';...
	'RERUN'              			'settings' 	'RERUN'                     '';...
	'MUTEXID'            			'settings' 	'MUTEXID'                   '';...
    'DELETETEMP'                    'settings' 	'DELETETEMP'                '';...
	'bLesionFilling'     			'settings' 	'bLesionFilling'            '';...
	'stopAfterErrors'    			'settings' 	'stopAfterErrors'           '';...
	'bSkipIfNoASL'    				'settings' 	'bSkipIfNoASL'              '';...
	'bSkipIfNoM0'    				'settings' 	'bSkipIfNoM0'               '';...
	'bSkipIfNoFlair'    			'settings' 	'bSkipIfNoFlair'            '';...
    'bRunModule_DARTEL'       		'modules'   'bRunDARTEL'                '';...
    'bRunModule_LongReg'       		'modules'   'bRunLongReg'               '';...
    'ForceInclusionList'       		'dataset'   'ForceInclusionList'        '';...
    'bAutomaticallyDetectFS'        'external'  'bAutomaticallyDetectFS'    '';...
    'MakeNIfTI4DICOM'       		'settings'  'MakeNIfTI4DICOM'           '';...
	'SegmentSPM12'       			'modules'   'structural'                'bSegmentSPM12';...
    'T1BiasFieldRegularization' 	'modules'   'structural'                'T1BiasFieldRegularization';...
    'bNativeSpaceAnalysis'       	'modules'   'population'                'bNativeSpaceAnalysis';...
    'bHammersCAT12'       			'modules'   'structural'                'bHammersCAT12';...
    'bFixResolution'       			'modules'   'structural'                'bFixResolution';...
    'bDCTRegistration'              'modules'   'asl'                       'bDCTRegistration';...
    'bAffineRegistration'           'modules'   'asl'                       'bAffineRegistration';...
    'bUseMNIasDummyStructural'      'modules'   'asl'                       'bUseMNIasDummyStructural';...
    'PVCNativeSpaceKernel'          'modules'   'asl'                       'PVCNativeSpaceKernel';...
    'bMakeNIfTI4DICOM'              'modules'   'asl'                       'bMakeNIfTI4DICOM';...
    'M0_GMScaleFactor'              'modules'   'asl'                       'M0_GMScaleFactor';...
    'DummyScanPositionInASL4D'      'modules'   'asl'                       'DummyScanPositionInASL4D';...
    'M0PositionInASL4D'             'modules'   'asl'                       'M0PositionInASL4D'
	};

end



