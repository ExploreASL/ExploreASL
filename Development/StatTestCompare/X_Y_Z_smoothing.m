function output_data_xyz = X_Y_Z_smoothing(input_data,X_FWHM,Y_FWHM,Z_FWHM)
%XY_Z_smoothing This function enables a different kernel in XY compared to
%Z-direction. May be useful to simulate e.g. GE 3D spiral acquisition.

% Smoothing kernel SD
% FWHM = 2*(2*log(2))^0.5 * SD = 2.3548 * SD

% FWHM of Gaussian bell in 3D-GRASE Z-axis is 1.38 times voxel-size (in
% Z-axis); in 3D FSE it is 1.034 times larger (Vidorreta et al.)
% GE_voxel_size=3.75x3.75x4;
% 3.75*1.38*1.034 = 5.35
% 4*1.38*1.034   = 5.7077 

% % default kernels are GE's kernels based on GE 3D FSE voxel sizes in VESPA
% if ~exist('X_FWHM')
%     X_FWHM                                            =2.5685;end
% if ~exist('Y_FWHM')
%     Y_FWHM                                            =2.5685;end
% if ~exist('Z_FWHM')
%     Z_FWHM                                            =5.7077;end

% Divide by 2.5 to obtain SD instead of FWHM (is an estimation)
X_spiral_smoothing_SD                                   =X_FWHM/ 2.3548;
Y_spiral_smoothing_SD                                   =Y_FWHM/ 2.3548;
Z_spiral_smoothing_SD                                   =Z_FWHM/ 2.3548;

myfilter_x                                              =fspecial('gaussian',[size(input_data,1) size(input_data,2)], X_spiral_smoothing_SD);
myfilter_y                                              =fspecial('gaussian',[size(input_data,2) size(input_data,1)], Y_spiral_smoothing_SD);
myfilter_z                                              =fspecial('gaussian',[size(input_data,1) size(input_data,3)], Z_spiral_smoothing_SD);

% Position in 2D volume
myfilter_x_2D                                           =zeros(size(myfilter_x,1),size(myfilter_x,2));
myfilter_x_2D(round(size(myfilter_x,1)/2),:)            =myfilter_x(round(size(myfilter_x,1)/2),:);

myfilter_y_2D                                           =zeros(size(myfilter_y,2),size(myfilter_y,1));
myfilter_y_2D(:,round(size(myfilter_y,1)/2))            =myfilter_y(round(size(myfilter_y,1)/2),:);

myfilter_z_2D                                           =zeros(size(myfilter_z,2),size(myfilter_z,1));
myfilter_z_2D(:,round(size(myfilter_z,1)/2))            =myfilter_z(round(size(myfilter_z,1)/2),:);

% Position in 3D volume
myfilter_x_3D                                           =zeros(size(input_data,1),size(input_data,2),size(input_data,3));
myfilter_y_3D                                           =zeros(size(input_data,1),size(input_data,2),size(input_data,3));
myfilter_z_3D                                           =zeros(size(input_data,1),size(input_data,2),size(input_data,3));

myfilter_x_3D(round(size(myfilter_x,1)/2),:,round(size(myfilter_x_3D,3)/2))         =myfilter_x_2D(round(size(myfilter_x,1)/2),:);
myfilter_y_3D(:,round(size(myfilter_y,1)/2),round(size(myfilter_y_3D,3)/2))         =myfilter_y_2D(:,round(size(myfilter_y,1)/2));
myfilter_z_3D(round(size(myfilter_z_3D,1)/2),round(size(myfilter_z_3D,2)/2),:)      =myfilter_z_2D(:,round(size(myfilter_z,1)/2));

% Convolve xyz PSF
%output_data_x                                           =imfilter(input_data   , myfilter_x_3D, 'replicate');
%output_data_xy                                          =imfilter(output_data_x, myfilter_y_3D ,'replicate');
%output_data_xyz                                         =imfilter(output_data_xy,myfilter_z_3D ,'replicate');
output_data_xyz = xASL_im_conv3Dsep(output_data_xy,myfilter_x_3D(:),myfilter_y_3D(:),myfilter_z_3D(:));


end

