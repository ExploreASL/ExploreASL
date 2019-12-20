function coord = indx2coord(indx, nx, ny)
    z = ceil(indx / (nx * ny));
    indx2 = indx - (z - 1) * nx * ny;
    y = ceil(indx2 / nx);
    x = indx2 - (y - 1) * nx;
    %x = ceil(indx3 / ny);
    coord = [x, y, z];
end
