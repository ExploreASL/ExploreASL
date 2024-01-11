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
% Copyright (C) 2015-2024 ExploreASL

nameConversionTable = ...
    {'nSessions'                    'dataset'   'nSessions'                          '';...
    'nSubjects'                     'dataset'   'nSubjects'                          '';...
    'nTotalSubjects'                'dataset'   'nTotalSubjects'                     '';...
    'TotalSubjects'                 'dataset'   'TotalSubjects'                      '';...
    'nTimePoints'                   'dataset'   'nTimePoints'                        '';...
    'exclusion'                     'dataset'   'exclusion'                          '';...
    'TimePointSubjects'             'dataset'   'TimePointSubjects'                  '';...
    'nTimePointSubjects'            'dataset'   'nTimePointSubjects'                 '';...
    'ExcludedSubjects'              'dataset'   'ExcludedSubjects'                   '';...
    'nTimePointExcluded'            'dataset'   'nTimePointExcluded'                 '';...
    'TimePointTotalSubjects'        'dataset'   'TimePointTotalSubjects'             '';...
    'TimePointExcluded'             'dataset'   'TimePointExcluded'                  '';...
    'nTimePointTotalSubjects'       'dataset'   'nTimePointTotalSubjects'            '';...
    'nTimePointsTotal'              'dataset'   'nTimePointsTotal'                   '';...
    'nExcluded'                     'dataset'   'nExcluded'                          '';...
    'TotalInclusionList'            'dataset'   'TotalInclusionList'                 '';...
    'name'                          'dataset'   'name'                               '';...
    'subject_regexp'                'dataset'   'subjectRegexp'                      '';...
    'subjectRegexp'                 'dataset'   'subjectRegexp'                      '';...
    'Sequence'                      'Q'         'Sequence'                           '';...
    'Vendor'                        'Q'         'Vendor'                             '';...
    'readout_dim'                   'Q'         'readoutDim'                         '';...
    'readoutDim'                    'Q'         'readoutDim'                         '';...
    'BackGrSupprPulses'             'Q'         'BackGrSupprPulses'                  '';...
    'LabelingType'                  'Q'         'LabelingType'                       '';...
    'Initial_PLD'                   'Q'         'Initial_PLD'                        '';...
    'LabelingDuration'              'Q'         'LabelingDuration'                   '';...
    'SliceReadoutTime'              'Q'         'SliceReadoutTime'                   '';...
    'Lambda'                        'Q'         'Lambda'                             '';...
    'T2art'                         'Q'         'T2art'                              '';...
    'BloodT1'                       'Q'         'BloodT1'                            '';...
    'TissueT1'                      'Q'         'TissueT1'                           '';...
    'nCompartments'                 'Q'         'nCompartments'                      '';...
    'NumberOfAverages'              'Q'         'NumberOfAverages'                   '';...
    'LabelingEfficiency'            'Q'         'LabelingEfficiency'                 '';...
    'ATT'                           'Q'         'ATT'                                '';...
    'M0'                            'Q'         'M0'                                 '';...
    'bUseBasilQuantification'       'Q'         'bUseBasilQuantification'            '';...
    'BackgroundSuppressionPulseTime'     'Q'    'BackgroundSuppressionPulseTime'     '';...
    'BackgroundSuppressionNumberPulses'  'Q'    'BackgroundSuppressionNumberPulses'  '';...
    'SpaghettiDir'                  'D'         'SpaghettiDir'                       '';...
    'HistogramDir'                  'D'         'HistogramDir'                       '';...
    'StatsMaps'                     'D'         'StatsMaps'                          '';...
    'SPMDIR'                        'D'         'SPMDIR'                             '';...
    'SPMpath'                       'D'         'SPMpath'                            '';...
    'LockDir'                       'dir'       'LockDir'                            '';...
    'SUBJECTDIR'                    'dir'       'SUBJECTDIR'                         '';...
    'dataParType'                   'opts'      'dataParType'                        '';...
    'SPMVERSION'                    'external'  'SPMVERSION'                         '';...
    'Quality'                       'settings'  'Quality'                            '';...
    'BILAT_FILTER'                  'settings'  'BILAT_FILTER'                       '';...
    'Pediatric_Template'            'settings'  'Pediatric_Template'                 '';...
    'bReproTesting'                 'settings'  'bReproTesting'                      '';...
    'bOverwrite'                    'settings'  'bOverwrite'                         '';...
    'dryRun'                        'settings'  'dryRun'                             '';...
    'bAutoACPC'                     'settings'  'bAutoACPC'                          '';...
    'RERUN'                         'settings'  'RERUN'                              '';...
    'MUTEXID'                       'settings'  'MUTEXID'                            '';...
    'DELETETEMP'                    'settings'  'DELETETEMP'                         '';...
    'bLesionFilling'                'settings'  'bLesionFilling'                     '';...
    'stopAfterErrors'               'settings'  'stopAfterErrors'                    '';...
    'bSkipIfNoASL'                  'settings'  'bSkipIfNoASL'                       '';...
    'bSkipIfNoM0'                   'settings'  'bSkipIfNoM0'                        '';...
    'bSkipIfNoFlair'                'settings'  'bSkipIfNoFlair'                     '';...
    'bRunModule_DARTEL'             'modules'   'bRunDARTEL'                         '';...
    'bRunDARTEL'                    'modules'   'bRunDARTEL'                         '';...
    'bRunModule_LongReg'            'modules'   'bRunLongReg'                        '';...
    'bRunLongReg'                   'modules'   'bRunLongReg'                        '';...
    'ForceInclusionList'            'dataset'   'ForceInclusionList'                 '';...
    'bAutomaticallyDetectFS'        'external'  'bAutomaticallyDetectFS'             '';...
    'MakeNIfTI4DICOM'               'settings'  'MakeNIfTI4DICOM'                    '';...
    'Segment_SPM12'                 'modules'   'structural'                         'bSegmentSPM12';...
    'SegmentSPM12'                  'modules'   'structural'                         'bSegmentSPM12';...
    'bSegmentSPM12'                 'modules'   'structural'                         'bSegmentSPM12';...
    'T1BiasFieldRegularization'     'modules'   'structural'                         'T1BiasFieldRegularization';...
    'bHammersCAT12'                 'modules'   'structural'                         'bHammersCAT12';...
    'bFixResolution'                'modules'   'structural'                         'bFixResolution';...
    'bDCTRegistration'              'modules'   'asl'                                'bDCTRegistration';...
    'bAffineRegistration'           'modules'   'asl'                                'bAffineRegistration';...
    'bUseMNIasDummyStructural'      'modules'   'asl'                                'bUseMNIasDummyStructural';...
    'bPVCNativeSpace'               'modules'   'asl'                                'bPVCNativeSpace';...
	'PVCNativeSpaceKernel'          'modules'   'asl'                                'PVCNativeSpaceKernel';...
	'bPVCGaussianMM'                'modules'   'asl'                                'bPVCGaussianMM';...
    'bMakeNIfTI4DICOM'              'modules'   'asl'                                'bMakeNIfTI4DICOM';...
    'M0_GMScaleFactor'              'modules'   'asl'                                'M0_GMScaleFactor';...
    'DummyScanPositionInASL4D'      'modules'   'asl'                                'DummyScanPositionInASL4D';...
    'M0_conventionalProcessing'     'modules'   'asl'                                'M0_conventionalProcessing';...
    'M0PositionInASL4D'             'modules'   'asl'                                'M0PositionInASL4D';...
    'HadamardType'                  'modules'   'asl'                                'HadamardType';...
	'ApplyQuantification'           'modules'   'asl'                                'ApplyQuantification';...
    'SessionMergingList'            'modules'   'asl'                                'SessionMergingList';...
	'bQuantifyMultiTE'              'modules'   'asl'                                'bQuantifyMultiTE';...
	'bQuantifyMultiPLD'             'modules'   'asl'                                'bQuantifyMultiPLD';...
    'bNativeSpaceAnalysis'          'modules'   'population'                         'bNativeSpaceAnalysis'
    };

end




