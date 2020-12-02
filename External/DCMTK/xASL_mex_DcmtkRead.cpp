/*-----------------------------------------------------------------------------------------------------------------------------------------------------
 * xASL_mex_DcmtkRead.cpp
 * 
 * SHORT DESCRIPTION: It reads the DICOM file and parses its header
 * 
 * FORMAT: parms = xASL_mex_DcmtkRead(path, pixelData)
 *
 * INPUT:
 *     path - path to the dicom file
 *     pixelData (optional, default 0) - if 1 then reads the pixel data, if 0 then skips reading and provides and empty field
 * 
 * OUTPUT:
 *     parms - a structure containing the loaded parameters
 *
 * DESCRIPTION:
 *    Provides a limited interface to read the DICOM files in Matlab using the DCMTK library.
 *    Currently, it reads the following parameters from normal/enhanced DICOM (it uses the DCMTK tag definition for that):
 *    RepetitionTime, EchoTime, RescaleSlope, RescaleIntercept, NumberOfTemporalPositions, NumberOfAverages, AcquisitionTime, 
 *    MediaStorageSOPClassUID, Manufacturer, SeriesDescription, ProtocolName, StudyDate, SeriesTime, StudyInstanceUID, 
 *    SeriesInstanceUID, ImageType, AcquisitionDate, SeriesDate, Rows, Columns, AcquisitionMatrix, InPlanePhaseEncodingDirection,
 *    AcquisitionContrast, ComplexImageComponent, PulseSequenceName, InversionTime, TemporalPositionIdentifier, SoftwareVersions,
 *    RWVIntercept, RWVSlope
 *    
 *    The image data are stored in:
 *    PixelData
 *    NOTE: to get pixel data equivalent to Matlab's dicomread, RESHAPE matrix and TRANSPOSE!
 *    img = reshape( dcm.PixelData, dcm.rows, dcm.columns )';
 *   
 *    Additionally, the following private tags are acquired:
 * 	  AssetRFactor - 0x0043, 0x1083
 *    EffectiveEchoSpacing - 0x0043, 0x192c
 *    MRSeriesWaterFatShift - 0x2001, 0x1022
 *    MRSeriesEPIFactor - 0x2001, 0x1013
 *    BandwidthPerPixelPhaseEncode - 0x0019, 0x1028
 *    MRScaleSlope - 0x2005, 0x100e or 
 *                   0x2005, 0x110e or
 *                   0x2005, 0x120e
 *    RescaleSlopeOriginal - 0x2005, 0x140a or 0x2005, 0x110a
 *    GELabelingType - 0x0019, 0x109C
 *    GELabelingDuration - 0x0043, 0x10A5
 *    PhilipsNumberTemporalScans - 0x2001, 0x1008
 *    PhilipsLabelControl - 0x2005, 0x1429
 *    PhoenixProtocol - 0x0029, 0x1020
 *    SiemensSliceTime - 0x0019, 0x0010
 *
 *    What is read is hard coded - to change that, you need to change the MEX file
 *
 * COMPILATION: 
 *    It requires the DCMTK library installed and compiled with the same compiler
 *    You need to specify the path to the include/lib directories in the attached makefile
 *-----------------------------------------------------------------------------------------------------------------------------------------------------
 * REFERENCES:
 * OFFIS DCMTK https://support.dcmtk.org/redmine/projects/dcmtk/wiki/Overview
 * __________________________________
 * Copyright © 2015-2018 ExploreASL
 * 
 * 2010-07-30 Joost Kuijer - VUmc - original version
 * 2018-12-17 Jan Petr - adapted for xASL use
 * 2020-11-12 Jan Petr - added more private tags
 */

// Matlab MEX stuff
#include "mex.h"
#include "matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <string.h>

#include <iostream>
// DCMTK dicom library
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"


// Static global variable points to the dicom file, declared static in Matlab
// main interface function. This way it can be accessed at a later call 
// from Matlab, and is globally accessible.
// Note: beter way to do this is to implement as a C++ Class....
// another Note: static globals have internal linkage.
static bool bFileRead_g = false;
//static DcmFileFormat *pdcmMyFile_g = NULL;


///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: long int --> double
///////////////////////////////////////////////////////////////////
mxArray *MXAGetLongIntAsDouble(DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
	static char	szModule[] = "MXAGetLongIntAsDouble";
	char        szErrMsgTxt[2048];
	mxArray 	*pmxNumMat;
	double 	*pdmxData;
	long int    diData = 0;
	
	// Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
		
		// get data
		if ( dcmItem->findAndGetLongInt( theTagKey, diData ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxNumMat = NULL;
		}
		else
		{
			// create mx double
			pmxNumMat = mxCreateNumericMatrix( 1, 1, mxDOUBLE_CLASS, mxREAL );
			// assign value
			pdmxData = (double*) mxGetData( pmxNumMat );
			*pdmxData = (double) diData;
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexWarnMsgTxt( szErrMsgTxt );
		pmxNumMat = NULL;//mxCreateNumericMatrix( 0, 0, mxDOUBLE_CLASS, mxREAL );
	}
	
	return pmxNumMat;
}

///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: Float64 --> double
///////////////////////////////////////////////////////////////////
mxArray *MXAGetFloat64AsDouble( DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
    static char	szModule[] = "MXAGetFloat64AsDouble";
    char        szErrMsgTxt[2048];
    mxArray 	*pmxNumMat;
    double 	*pdmxData;
    Float64     f64Data = 0;
    
	// Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
		
		// get data
		if ( dcmItem->findAndGetFloat64( theTagKey, f64Data ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxNumMat = NULL;
		}
		else
		{
			// create mx double
			pmxNumMat = mxCreateNumericMatrix( 1, 1, mxDOUBLE_CLASS, mxREAL );
			// assign value
			pdmxData = (double*) mxGetData( pmxNumMat );
			*pdmxData = (double) f64Data;
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexWarnMsgTxt( szErrMsgTxt );
		pmxNumMat = NULL;//mxCreateNumericMatrix( 0, 0, mxDOUBLE_CLASS, mxREAL );
	}
	return pmxNumMat;
}

///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: Float32 --> double
///////////////////////////////////////////////////////////////////
mxArray *MXAGetFloat32AsDouble( DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
    static char	szModule[] = "MXAGetFloat32AsDouble";
    char        szErrMsgTxt[2048];
    mxArray 	*pmxNumMat;
    double 	*pdmxData;
    Float32    f32Data = 0;
    
	// Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
		// get data
		if ( dcmItem->findAndGetFloat32( theTagKey, f32Data, 0, OFTrue ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxNumMat = NULL;
		}
		else
		{
			// create mx double
			pmxNumMat = mxCreateNumericMatrix( 1, 1, mxDOUBLE_CLASS, mxREAL );
			// assign value
			pdmxData = (double*) mxGetData( pmxNumMat );
			*pdmxData = (double) f32Data;
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexWarnMsgTxt( szErrMsgTxt );
		pmxNumMat = NULL;//mxCreateNumericMatrix( 0, 0, mxDOUBLE_CLASS, mxREAL );
	}
	return pmxNumMat;
}

///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: Float64Array --> double
///////////////////////////////////////////////////////////////////
mxArray *MXAGetFloat32ArrayAsDouble(DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
    static char	szModule[] = "MXAGetFloat32ArrayAsDouble";
    char        szErrMsgTxt[2048];
    mxArray 	*pmxNumMat;
    double 	*pdmxData;
    const Float32 * pf32Data = NULL;
	unsigned long ul_nData;
    
	// Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
		
		// get data
		if ( dcmItem->findAndGetFloat32Array( theTagKey, pf32Data, &ul_nData ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxNumMat = NULL;
		}
		else
		{
			// create mx double
			pmxNumMat = mxCreateNumericMatrix( 1, ul_nData, mxSINGLE_CLASS, mxREAL );
			// assign value
			pdmxData = (double*) mxGetData( pmxNumMat );
			memcpy( pdmxData, (const void*) pf32Data, ul_nData * sizeof(Float32) );
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexWarnMsgTxt( szErrMsgTxt );
		pmxNumMat = NULL;//mxCreateNumericMatrix( 0, 0, mxDOUBLE_CLASS, mxREAL );
	}
	return pmxNumMat;
}

///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: Float64Array --> double
///////////////////////////////////////////////////////////////////
mxArray *MXAGetFloat64ArrayAsDouble(DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
    static char	szModule[] = "MXAGetFloat64ArrayAsDouble";
    char        szErrMsgTxt[2048];
    mxArray 	*pmxNumMat;
    double 	*pdmxData;
    const Float64 * pf64Data = NULL;
	unsigned long ul_nData;
    
	// Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
		
		// get data
		if ( dcmItem->findAndGetFloat64Array( theTagKey, pf64Data, &ul_nData ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxNumMat = NULL;
		}
		else
		{
			// create mx double
			pmxNumMat = mxCreateNumericMatrix( 1, ul_nData, mxDOUBLE_CLASS, mxREAL );
			// assign value
			pdmxData = (double*) mxGetData( pmxNumMat );
			memcpy( pdmxData, (const void*) pf64Data, ul_nData * sizeof(Float64) );
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexWarnMsgTxt( szErrMsgTxt );
		pmxNumMat = NULL;//mxCreateNumericMatrix( 0, 0, mxDOUBLE_CLASS, mxREAL );
	}
	return pmxNumMat;
}

///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: Int16 array
// NOTE: DCMTK gives UNSIGNED int but apparently both GE and Matlab intepret these as SIGNED!!
///////////////////////////////////////////////////////////////////
mxArray *MXAGetInt16Array(DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
    static char	szModule[] = "MXAGetUint16Array";
    char        szErrMsgTxt[2048];
    mxArray 	*pmxNumMat;
    void 	*pmxData;
    const Uint16 *pui16Data = NULL;
    unsigned long ul_nData;
	
	// Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
		// get data
		if ( dcmItem->findAndGetUint16Array( theTagKey, pui16Data, &ul_nData ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxNumMat = NULL;
		}
		else
		{
			// create mx int16
			// NOTE: DCMTK gives UNSIGNED int but apparently both GE and Matlab intepret these as SIGNED!!
			pmxNumMat = mxCreateNumericMatrix( 1, ul_nData, mxINT16_CLASS, mxREAL );
			// assign value
			pmxData = mxGetData( pmxNumMat );
			memcpy( pmxData, (const void*) pui16Data, ul_nData * sizeof(Uint16) );
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexPrintf( szErrMsgTxt );
		pmxNumMat = NULL;//mxCreateNumericMatrix( 0, 0, mxINT16_CLASS, mxREAL );
	}
    
    return pmxNumMat;
}

///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: string array
///////////////////////////////////////////////////////////////////
mxArray *MXAGetStringArray(DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
    static char	szModule[] = "MXAGetString";
    char        szErrMsgTxt[2048];
    mxArray 	*pmxStrMat;
    OFString    strData;
    
    // Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
	
		if ( dcmItem->findAndGetOFStringArray( theTagKey, strData ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxStrMat = NULL;
		}
		else
		{
			// create mx double
			pmxStrMat = mxCreateString( strData.c_str() );
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexPrintf( szErrMsgTxt );
		pmxStrMat = NULL;//mxCreateString("");
	}

    return pmxStrMat;
}

///////////////////////////////////////////////////////////////////
// Get a field from the DICOM header: string
///////////////////////////////////////////////////////////////////
mxArray *MXAGetString(DcmItem * dcmItem, const DcmTagKey &theTagKey )
{
    static char	szModule[] = "MXAGetString";
    char        szErrMsgTxt[2048];
    mxArray 	*pmxStrMat;
    OFString    strData;
    
    // Check if the tag exists and is filled, otherwise issue a warning
	if ( dcmItem->tagExistsWithValue( theTagKey ) == OFTrue )
	{
	
		if ( dcmItem->findAndGetOFString( theTagKey, strData ).bad() )
		{
			/* error getting the element */
			sprintf( szErrMsgTxt, "%s: cannot get element %s\n", szModule, theTagKey.toString().c_str() );
			mexWarnMsgTxt(szErrMsgTxt);
			//mexErrMsgTxt(szErrMsgTxt);
			pmxStrMat = NULL;
		}
		else
		{
			// create mx double
			pmxStrMat = mxCreateString( strData.c_str() );
		}
	}
	else
	{
		/* Warning element missing */
		//sprintf( szErrMsgTxt, "xASL_mex_DcmtkRead:%s: element does not exist %s\n", szModule, theTagKey.toString().c_str() );
		//mexPrintf( szErrMsgTxt );
		pmxStrMat = NULL;//mxCreateString("");
	}

    return pmxStrMat;
}


///////////////////////////////////////////////////////////////////
// Read the DICOM file
///////////////////////////////////////////////////////////////////
void VMatDcmtkRead( DcmFileFormat * DcmMyFile, char *pchFileName, mxArray *pmxOutput, int readPixel )
{
    static char	szModule[] = "VMatDcmtkRead";
    int   	i;
    OFCondition dcmStatus;
    char        szErrMsgTxt[2048];
    mxArray 	*pmxNumMat;
	int parmRead = 0;

	// The meta info and dataset objects
	DcmMetaInfo * metaInfo = NULL;
	DcmDataset * dataset = NULL;
	
	DcmItem * 	rwItem     = NULL;
	DcmItem * 	sharedItem = NULL;
	DcmItem *   perFrmItem = NULL;
	DcmItem * 	timingItem = NULL;
	DcmItem * 	echoItem   = NULL;
	DcmItem *   pixelItem  = NULL;
	DcmItem * 	privatItem = NULL;
        
	// Load the file
	dcmStatus = DcmMyFile->loadFile( pchFileName );
    if ( dcmStatus.bad() )
    {
        sprintf( szErrMsgTxt, "%s: cannot read DICOM file \"%s\" (%s)\n", szModule, pchFileName, dcmStatus.text() );
        mexErrMsgTxt( szErrMsgTxt ); 
    }

	// Get pointers to the meta info and dataset objects
	metaInfo = DcmMyFile->getMetaInfo();
	dataset  = DcmMyFile->getDataset();
	
    // Flag that DICOM file was read
    bFileRead_g = true;

	// Reading the Enhanced DICOM -> timingItem
	// SharedFunctionalGroupsSequence.Item_1.MRTimingAndRelatedParametersSequence.Item_1
	dataset->findAndGetSequenceItem(	DCM_SharedFunctionalGroupsSequence, sharedItem, 0 );
	if ( sharedItem )
	{
		sharedItem->findAndGetSequenceItem(	DCM_MRTimingAndRelatedParametersSequence, timingItem, 0 );
	}
	
	// Reading the Enhanced DICOM -> echoItem, pixelItem, privatItem
	// PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1
	// PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1
	// PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1
	dataset->findAndGetSequenceItem(	DCM_PerFrameFunctionalGroupsSequence, perFrmItem, 0 );
	if ( perFrmItem )
	{
		perFrmItem->findAndGetSequenceItem(	DCM_MREchoSequence , echoItem, 0 );
		perFrmItem->findAndGetSequenceItem(	DCM_PixelValueTransformationSequence, pixelItem, 0 );
		perFrmItem->findAndGetSequenceItem(	DcmTagKey(0x2005, 0x140f), privatItem, 0 );
	}
	
	//RealWorldValueMappingSequence.Item_1
	dataset->findAndGetSequenceItem( DCM_RealWorldValueMappingSequence , rwItem, 0);
		
    //mxSetField( pmxOutput, 0, "PatientID",       MXAGetString         ( DCM_PatientID      ) );  
	
	// Obtain these parameters always from normal DICOM
	mxSetField( pmxOutput, 0, "AcquisitionTime"          , MXAGetStringArray( dataset,        DCM_AcquisitionTime   ) );
	mxSetField( pmxOutput, 0, "Manufacturer"             , MXAGetString( dataset,             DCM_Manufacturer      ) );
	mxSetField( pmxOutput, 0, "SeriesDescription"        , MXAGetString( dataset,             DCM_SeriesDescription ) );
	mxSetField( pmxOutput, 0, "ProtocolName"             , MXAGetString( dataset,             DCM_ProtocolName      ) );
	mxSetField( pmxOutput, 0, "SeriesTime"               , MXAGetStringArray( dataset,        DCM_SeriesTime        ) );
	mxSetField( pmxOutput, 0, "StudyInstanceUID"         , MXAGetString( dataset,             DCM_StudyInstanceUID  ) );	
	mxSetField( pmxOutput, 0, "SeriesInstanceUID"        , MXAGetString( dataset,             DCM_SeriesInstanceUID ) );	
	mxSetField( pmxOutput, 0, "StudyDate"                , MXAGetStringArray( dataset,        DCM_StudyDate         ) );
	mxSetField( pmxOutput, 0, "ImageType"                , MXAGetStringArray( dataset,        DCM_ImageType         ) );	
	mxSetField( pmxOutput, 0, "MediaStorageSOPClassUID"  , MXAGetString( metaInfo,            DCM_MediaStorageSOPClassUID ) );	
	mxSetField( pmxOutput, 0, "AcquisitionDate"          , MXAGetStringArray( dataset,        DCM_AcquisitionDate   ) );
	mxSetField( pmxOutput, 0, "SeriesDate"               , MXAGetStringArray( dataset,        DCM_SeriesDate        ) );
	mxSetField( pmxOutput, 0, "Rows"                     , MXAGetInt16Array(  dataset,        DCM_Rows    ));
	mxSetField( pmxOutput, 0, "Columns"                  , MXAGetInt16Array(  dataset,        DCM_Columns ));
	mxSetField( pmxOutput, 0, "AcquisitionContrast"      , MXAGetString( dataset,             DCM_AcquisitionContrast ) );	
	mxSetField( pmxOutput, 0, "ComplexImageComponent"    , MXAGetString( dataset,             DCM_ComplexImageComponent ) );	
	mxSetField( pmxOutput, 0, "PulseSequenceName"        , MXAGetString( dataset,             DCM_PulseSequenceName ) );	
	mxSetField( pmxOutput, 0, "InversionTime"            , MXAGetFloat64AsDouble( dataset,    DCM_InversionTime       ) );
	mxSetField( pmxOutput, 0, "SoftwareVersions"         , MXAGetStringArray( dataset,        DCM_SoftwareVersions        ) );
	
	
	// Parameters for EPI readout needed for TopUp
	mxSetField( pmxOutput, 0, "AssetRFactor"                , MXAGetStringArray(dataset, DcmTagKey(0x0043, 0x1083)));
	mxSetField( pmxOutput, 0, "EffectiveEchoSpacing"        , MXAGetStringArray(dataset, DcmTagKey(0x0043, 0x192c)));
	mxSetField( pmxOutput, 0, "AcquisitionMatrix"           , MXAGetInt16Array(     dataset, DCM_AcquisitionMatrix    ));
	mxSetField( pmxOutput, 0, "MRSeriesWaterFatShift"       , MXAGetStringArray(dataset, DcmTagKey(0x2001, 0x1022)));
	mxSetField( pmxOutput, 0, "MRSeriesEPIFactor"           , MXAGetStringArray(dataset, DcmTagKey(0x2001, 0x1013)));
	mxSetField( pmxOutput, 0, "BandwidthPerPixelPhaseEncode", MXAGetStringArray(dataset, DcmTagKey(0x0019, 0x1028)));
	mxSetField( pmxOutput, 0, "InPlanePhaseEncodingDirection", MXAGetStringArray(dataset, DCM_InPlanePhaseEncodingDirection));
	 
	// Extra GE and Siemens parameters
	mxSetField( pmxOutput, 0, "GELabelingType"             , MXAGetString( dataset,       DcmTagKey(0x0019, 0x109C) ) );
	mxSetField( pmxOutput, 0, "GELabelingDuration"         , MXAGetLongIntAsDouble( dataset, DcmTagKey(0x0043, 0x10a5)     ) );
	mxSetField( pmxOutput, 0, "PhoenixProtocol"            , MXAGetStringArray(dataset, DcmTagKey(0x0029, 0x1020)));
	mxSetField( pmxOutput, 0, "SiemensSliceTime"           , MXAGetFloat64ArrayAsDouble( dataset,    DcmTagKey(0x0019, 0x1029)  ) );
			
	if ((rwItem) && (rwItem->tagExistsWithValue(DCM_RealWorldValueIntercept) == OFTrue))
		mxSetField( pmxOutput, 0, "RWVIntercept"         , MXAGetFloat64AsDouble( rwItem,    DCM_RealWorldValueIntercept       ) );
	if ((rwItem) && (rwItem->tagExistsWithValue(DCM_RealWorldValueSlope) == OFTrue))
		mxSetField( pmxOutput, 0, "RWVSlope"             , MXAGetFloat64AsDouble( rwItem,    DCM_RealWorldValueSlope       ) );
	
	// Read the Repetition time from either the enhanced or normal DICOM
	if ( ( timingItem ) && ( timingItem->tagExistsWithValue( DCM_RepetitionTime ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "RepetitionTime"       , MXAGetFloat64AsDouble( timingItem, DCM_RepetitionTime ) );
	else
		mxSetField( pmxOutput, 0, "RepetitionTime"       , MXAGetFloat64AsDouble( dataset,    DCM_RepetitionTime ) );
	
	
	// Read the Echo time from either the enhanced or normal DICOM
	if ( ( echoItem ) && ( echoItem->tagExistsWithValue( DCM_EffectiveEchoTime ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "EchoTime"             , MXAGetFloat64AsDouble( echoItem,   DCM_EffectiveEchoTime) );
	else
		mxSetField( pmxOutput, 0, "EchoTime"             , MXAGetFloat64AsDouble( dataset,    DCM_EchoTime       ) );
	
	// Read the rescale slopes and intercept from normal/enhanced DICOM
	if ( ( pixelItem ) && ( pixelItem->tagExistsWithValue( DCM_RescaleSlope ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "RescaleSlope"             , MXAGetFloat64AsDouble( pixelItem,    DCM_RescaleSlope       ) );
	else
		mxSetField( pmxOutput, 0, "RescaleSlope"             , MXAGetFloat64AsDouble( dataset,    DCM_RescaleSlope       ) );
	
	if ( ( pixelItem ) && ( pixelItem->tagExistsWithValue( DCM_RescaleIntercept ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "RescaleIntercept"         , MXAGetFloat64AsDouble( pixelItem,    DCM_RescaleIntercept       ) );
	else
		mxSetField( pmxOutput, 0, "RescaleIntercept"         , MXAGetFloat64AsDouble( dataset,    DCM_RescaleIntercept       ) );

	
	// Read the temporalpositions and private rescale tags from Philips private tags
	if ( ( privatItem ) && ( privatItem->tagExistsWithValue( DCM_NumberOfTemporalPositions ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "NumberOfTemporalPositions", MXAGetLongIntAsDouble( privatItem,    DCM_NumberOfTemporalPositions) );
	else
		mxSetField( pmxOutput, 0, "NumberOfTemporalPositions", MXAGetLongIntAsDouble( dataset,    DCM_NumberOfTemporalPositions) );
	if ( ( privatItem ) && ( privatItem->tagExistsWithValue( DCM_TemporalPositionIdentifier ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "TemporalPositionIdentifier", MXAGetLongIntAsDouble( privatItem,    DCM_TemporalPositionIdentifier) );
	else
		mxSetField( pmxOutput, 0, "TemporalPositionIdentifier", MXAGetLongIntAsDouble( dataset,    DCM_TemporalPositionIdentifier) );
	
	if ( ( privatItem ) && ( privatItem->tagExistsWithValue( DcmTagKey(0x2001, 0x1008) ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "PhilipsNumberTemporalScans", MXAGetStringArray( privatItem,    DcmTagKey(0x2001, 0x1008)) );
	else
		mxSetField( pmxOutput, 0, "PhilipsNumberTemporalScans", MXAGetStringArray( dataset,    DcmTagKey(0x2001, 0x1008)) );
	
		if ( ( privatItem ) && ( privatItem->tagExistsWithValue( DcmTagKey(0x2005, 0x1429) ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "PhilipsLabelControl", MXAGetStringArray( privatItem,    DcmTagKey(0x2005, 0x1429)) );
	else
		mxSetField( pmxOutput, 0, "PhilipsLabelControl", MXAGetStringArray( dataset,    DcmTagKey(0x2005, 0x1429)) );
	
	
	if ((privatItem) && 
		((privatItem->tagExistsWithValue(DcmTagKey(0x2005, 0x100e))==OFTrue)||
		 (privatItem->tagExistsWithValue(DcmTagKey(0x2005, 0x110e))==OFTrue)|| 	
		 (privatItem->tagExistsWithValue(DcmTagKey(0x2005, 0x120e))==OFTrue)))
	{
		if (privatItem->tagExistsWithValue(DcmTagKey(0x2005, 0x120e))==OFTrue)
			mxSetField( pmxOutput, 0, "MRScaleSlope"        , MXAGetStringArray( privatItem,             DcmTagKey(0x2005, 0x120e) ) );
		else
		{
			if (privatItem->tagExistsWithValue(DcmTagKey(0x2005, 0x110e))==OFTrue)
				mxSetField( pmxOutput, 0, "MRScaleSlope"        , MXAGetStringArray( privatItem,             DcmTagKey(0x2005, 0x110e) ) );
			else
				mxSetField( pmxOutput, 0, "MRScaleSlope"        , MXAGetStringArray( privatItem,             DcmTagKey(0x2005, 0x100e) ) );
		}
		
		if (privatItem->tagExistsWithValue(DcmTagKey(0x2005, 0x140a))==OFTrue)
			mxSetField( pmxOutput, 0, "RescaleSlopeOriginal", MXAGetStringArray( privatItem,             DcmTagKey(0x2005, 0x140a) ) );
		else
			mxSetField( pmxOutput, 0, "RescaleSlopeOriginal", MXAGetStringArray( privatItem,             DcmTagKey(0x2005, 0x110a) ) );
	}
	else
	{
		if (dataset->tagExistsWithValue(DcmTagKey(0x2005, 0x120e))==OFTrue)
			mxSetField( pmxOutput, 0, "MRScaleSlope"        , MXAGetStringArray( dataset,             DcmTagKey(0x2005, 0x120e) ) );
		else
		{
			if (dataset->tagExistsWithValue(DcmTagKey(0x2005, 0x110e))==OFTrue)
				mxSetField( pmxOutput, 0, "MRScaleSlope"        , MXAGetStringArray( dataset,             DcmTagKey(0x2005, 0x110e) ) );
			else
				mxSetField( pmxOutput, 0, "MRScaleSlope"        , MXAGetStringArray( dataset,             DcmTagKey(0x2005, 0x100e) ) );
		}
		
		if (dataset->tagExistsWithValue(DcmTagKey(0x2005, 0x140a))==OFTrue)
			mxSetField( pmxOutput, 0, "RescaleSlopeOriginal", MXAGetStringArray( dataset,             DcmTagKey(0x2005, 0x140a) ) );
		else
			mxSetField( pmxOutput, 0, "RescaleSlopeOriginal", MXAGetStringArray( dataset,             DcmTagKey(0x2005, 0x110a) ) );
	}
		
	if ( ( privatItem ) && ( privatItem->tagExistsWithValue( DCM_NumberOfAverages ) == OFTrue ) )
		mxSetField( pmxOutput, 0, "NumberOfAverages"         , MXAGetFloat64AsDouble( privatItem,    DCM_NumberOfAverages ) );
	else
		mxSetField( pmxOutput, 0, "NumberOfAverages"         , MXAGetFloat64AsDouble( dataset,    DCM_NumberOfAverages ) );
	
    // Note: inconsistent naming between Matlab and DCMTK: EchoNumber and DCM_EchoNumbers
    //mxSetField( pmxOutput, 0, "EchoNumber",      MXAGetLongIntAsDouble( DCM_EchoNumbers    ) );
    //mxSetField( pmxOutput, 0, "EchoTrainLength", MXAGetLongIntAsDouble( DCM_EchoTrainLength) );
    //mxSetField( pmxOutput, 0, "ReceiveCoilName", MXAGetString         ( DCM_ReceiveCoilName) );
    //mxSetField( pmxOutput, 0, "SeriesNumber",    MXAGetLongIntAsDouble( DCM_SeriesNumber   ) );
    //mxSetField( pmxOutput, 0, "InstanceNumber",  MXAGetLongIntAsDouble( DCM_InstanceNumber ) );
    // Note: inconsistent naming between Matlab and DCMTK: rows / columns without capital...?
    //mxSetField( pmxOutput, 0, "rows",            MXAGetLongIntAsDouble( DCM_Rows           ) );
    //mxSetField( pmxOutput, 0, "columns",         MXAGetLongIntAsDouble( DCM_Columns        ) );
    //mxSetField( pmxOutput, 0, "BitsAllocated",   MXAGetLongIntAsDouble( DCM_BitsAllocated  ) );
    //mxSetField( pmxOutput, 0, "BitsStored",      MXAGetLongIntAsDouble( DCM_BitsStored     ) );
    //mxSetField( pmxOutput, 0, "WindowCenter",    MXAGetFloat64AsDouble( DCM_WindowCenter   ) );
    //mxSetField( pmxOutput, 0, "WindowWidth",     MXAGetFloat64AsDouble( DCM_WindowWidth    ) );
    //mxSetField( pmxOutput, 0, "ImagesInAcquisition", MXAGetLongIntAsDouble( DCM_ImagesInAcquisition ) );
    //mxSetField( pmxOutput, 0, "SeriesDescription",   MXAGetString         ( DCM_SeriesDescription   ) );  
    //mxSetField( pmxOutput, 0, "SOPInstanceUID",      MXAGetString         ( DCM_SOPInstanceUID      ) ); 
    //mxSetField( pmxOutput, 0, "SeriesInstanceUID",   MXAGetString         ( DCM_SeriesInstanceUID   ) ); 
    
    // Get the pixel data
	if (readPixel)
	{
		mxSetField( pmxOutput, 0, "PixelData",  MXAGetInt16Array( dataset, DCM_PixelData ) );
	}
	else
	{
		mxSetField( pmxOutput, 0, "PixelData", mxCreateNumericMatrix( 0, 0, mxDOUBLE_CLASS, mxREAL ) );
	}
    // NOTE: to get pixel data equivalent to Matlabs dicomread, reshape matrix and TRANSPOSE!
    // img = reshape( a.PixelData, a.rows, a.columns )';
}


///////////////////////////////////////////////////////////////////
// Main function called by Matlab
///////////////////////////////////////////////////////////////////
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    static DcmFileFormat dcmMyFile;
    char *pchFileName;
    int   status;
    const char *pszFieldnames[] = \
    { \
        "RepetitionTime", "EchoTime", "RescaleSlope", "RescaleIntercept", \
        "NumberOfTemporalPositions", "NumberOfAverages", "AcquisitionTime", \
        "PixelData", "MediaStorageSOPClassUID", "Manufacturer",  "MRScaleSlope", \
		"StudyDate", "StudyInstanceUID", "SeriesInstanceUID", "ImageType", \
	    "SeriesDescription", "ProtocolName", "SeriesTime", "AcquisitionDate", "SeriesDate", \
		"AssetRFactor", "EffectiveEchoSpacing", "AcquisitionMatrix", "MRSeriesWaterFatShift", \
	    "MRSeriesEPIFactor", "BandwidthPerPixelPhaseEncode", "InPlanePhaseEncodingDirection", \
	    "Rows", "Columns", "RescaleSlopeOriginal", "RWVIntercept", "RWVSlope", \
	    "AcquisitionContrast", "ComplexImageComponent", "GELabelingType", "PulseSequenceName", \
		"InversionTime", "GELabelingDuration", "PhilipsNumberTemporalScans", \
		"PhilipsLabelControl", "TemporalPositionIdentifier", "PhoenixProtocol", "SoftwareVersions", \
		"SiemensSliceTime"
    };

    const int inFields = 44;
	int readPixel;
	double *tmp;

    // Make the dicom file known globally.
    //pdcmMyFile_g = &dcmMyFile;

    /* check for proper number of arguments */
    if ( nrhs < 1 ) 
        mexErrMsgTxt( "xASL_mex_DcmtkRead::At least one input required." );
	if ( nrhs > 2)
		mexErrMsgTxt( "xASL_mex_DcmtkRead::Two inputs are a maximum.");
	
	if (nrhs == 1)
	{
		readPixel = 0;
	}
	else
	{
		/* Reads the parameter if to read pixel data */
		if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) || mxIsEmpty(prhs[1]))
		{
			mexErrMsgTxt("xASL_mex_DcmtkRead::ima must be numeric, real, full, double, and not empty");
		}
		tmp = mxGetPr(prhs[1]);
		readPixel = tmp[0];
	}
	
    /* first input must be a string: DICOM file name */
    if ( !mxIsChar( prhs[0] ) )
        mexErrMsgTxt( "xASL_mex_DcmtkRead::First input must be a string." );
    
    /* input must be a row vector */
    if ( mxGetM( prhs[0] ) != 1 )
        mexErrMsgTxt( "xASL_mex_DcmtkRead::Input must be a row vector." );
    
    /* copy the string data from prhs[0] into a C string */
    pchFileName = mxArrayToString( prhs[0] );
    if ( pchFileName == NULL ) 
        mexErrMsgTxt( "xASL_mex_DcmtkRead::Could not convert filename to string." );
    
	/* Check output */
	if ( nlhs > 1 )
		mexErrMsgTxt( "xASL_mex_DcmtkRead::Too many output arguments." );
	
	/* create a structure for output */
	mxArray *pmxOutput = mxCreateStructMatrix( 1, 1, inFields, pszFieldnames );
	
	/* read the DICOM file */
	VMatDcmtkRead(&dcmMyFile, pchFileName, pmxOutput, readPixel );
	
	/* pmxOutput to MATLAB mexFunction output*/
	plhs[0] = pmxOutput;

    mxFree( pchFileName );

    return;
}
