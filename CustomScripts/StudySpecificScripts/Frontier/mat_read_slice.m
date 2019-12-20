% MAT_READ_SLICE Read subheader and 2D matrix data from matrix file.
%
%             [ SHEAD DATA ] = MAT_READ_SLICE( FILENAME, MATSPEC [ , SEGMENT ] )
%             reads a 2D matrix from the file FILENAME. 
%
%             The data read is scaled according to the quantification scale
%             factor from the subheader and returned in the 2D double array
%             DATA. The subheader is returned in the subheader structure SHEAD 
%             with the quantification scale factor set to 1. Use MAT_READ_SUBH
%             to read the actual quantification scale factor.
%
%             MATSPEC = [ FRAME PLANE GATE DATA BED ] is a five-element vector
%             specifying the frame, bed, gate, data and bed index of the
%             requested matrix. For single gate studies, set GATE to 1. For
%             single bed studies, set BED to 0. DATA is usually 0. For volume 
%             files, PLANE indicates the slice to be extracted.
%
%             The optional SEGMENT argument determines which segment of
%             a 3D Scan will be read.
%
%             MAT_READ_SLICE will abort in case of errors.
%
%             Bugs and Limitations:
%             Only works for ECAT files (except Norm3D).
%
%             See also:
%             MATFILE, MAT_CREATE, MAT_READ_MAINH, MAT_READ_SUBH,
%             MAT_READ_VOLUME, MAT_WRITE_MAINH, MAT_WRITE_SUBH, MAT_WRITE 
% 
%             Version 1.4 (11/27/03), Harald Fricke
