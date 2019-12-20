function indx = coord2indx(x, y, z, nx, ny)
    indx = (z - 1) * nx * ny + y * nx - (nx - x);
end

