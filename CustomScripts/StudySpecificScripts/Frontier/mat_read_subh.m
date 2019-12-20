% MAT_READ_SUBH  Read subheader from matrix file.
%
%             SHEAD = MAT_READ_SUBH( FILENAME, MATSPEC ) reads a subheader
%             from the file FILENAME and returns it in the subheader structure
%             SHEAD.
%
%             MATSPEC = [ FRAME PLANE GATE DATA BED ] is a five-element vector
%             specifying the frame, bed, gate, data and bed index of the
%             requested matrix. For single-gate studies, set GATE to 1. For
%             single bed studies, set BED to 0. DATA usually is 0.
%
%             MAT_READ_SUBH will abort in case of errors.
%
%             Bugs and Limitations:
%             Only works for ECAT files.
%
%             See also:
%             MATFILE, MAT_CREATE, MAT_READ_MAINH, MAT_READ_SLICE,
%             MAT_READ_VOLUME, MAT_WRITE_MAINH, MAT_WRITE_SUBH, MAT_WRITE 
% 
%             Version 1.3 (11/27/03), Harald Fricke
