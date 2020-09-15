%cat_vol_interp3f Region-wise statistic.
%  Fast nearest, bi-linear and bi-cubic interpolation for 3D image data on 
%  a regular grid that is used instand of the double-based MATLAB function.
%  This method handles the border by repeating the closest values to the 
%  point accessed.  This is different from matlabs border handling.
% 
%      R = cat_vol_interp3f(F, X, Y, Z, [method])
%      R = cat_vol_interp3f(Fx, Fy, Fz, F, X, Y, Z, [method])
% 
%  Fx, Fy, Fz .. Coordinate system in which F is given. Only the first and
%                last entry in Fx, Fy, Fz are used, and it is assumed that 
%                the inbetween values are linearly interpolated.
%  F          .. WxHxDxC Image with an arbitray number of channels C.
%  X, Y, Z    .. I_1 x ... x I_n matrices with the x and y coordinates to
%                interpolate.
%  R          .. I_1 x ... x I_n x C matrix, which contains the 
%                interpolated image channels.
%  method     .. nearest, linear, or cubic.
% 
% 
%  Examples:
% 
%     %% Interpolation of 3D volumes (e.g. distance transforms)
%     clear
%     sz=5;
% 
%     % Dist 
%     dist1.D = randn(sz,sz,sz);
%     [dist1.x dist1.y dist.z] = ...
%       meshgrid(linspace(-1,1,sz), linspace(-1,1,sz), linspace(-1,1,sz));
%     
%     R = [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
%     RD = R * [Dx(:)'; Dy(:)'] + 250;
%     RDx = reshape(RD(1,:), size(Dx));
%     RDy = reshape(RD(2,:), size(Dy));
%     
%     methods = {'nearest', 'linear', 'cubic'};
%     la=nan(1,3);
%     for i=1:3
%       la(i) = subplot(2,2,i);
%       tic;
%       IMG_R = ba_interp2(IMG, RDx, RDy, methods{i});
%       elapsed=toc;
%       imshow(IMG_R);
%       title(sprintf(['Rotation and zoom using %s interpolation ' ...
%         'took %gs'], methods{i}, elapsed));
%     end
%     linkaxes(la);
% 
%  See also compile.
%  ________________________________________________________________________
%  <a href="www.brian-amberg.de/">Brian Amberg, 2008</a> 
%  $Id: cat_vol_interp3f.m 1523 2019-11-21 23:12:24Z gaser $
