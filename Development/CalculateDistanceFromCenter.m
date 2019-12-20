% Are you trying to calculate the distance from the center? You can use a
% different and fast function rather than dilation or other morphological
% operators...

nx,ny,nz % is the matrix size
cx,cy,cz % is the center voxel coordinates

[x,y,z] = meshgrid(1:nx,1:ny,1:nz)
d = sqrt((x-cx).^2 + (y-cy).^2 + (z-cz).^2);