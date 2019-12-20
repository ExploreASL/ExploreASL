% MAT_READ_MAINH  Read main header from matrix file.
%
%             MHEAD = MAT_READ_MAINH( FILENAME ) returns the main header
%             fields of the file FILENAME in the main header structure MHEAD.
%
%             MAT_READ_MAINH will abort in case of errors.
%
%             Bugs and limitations:
%             Does not read the patient name from Interfiles.
%
%             See also:
%             MATFILE, MAT_CREATE, MAT_READ_SUBH, MAT_READ_SLICE,
%             MAT_READ_VOLUME, MAT_WRITE_MAINH, MAT_WRITE_SUBH, MAT_WRITE 
% 
%             Version 1.3 (11/27/03), Harald Fricke
