function neighborhood = getNeighborhood2(img, indx, neighbor_order)

%     getNeighborhood2.m - part of the LST toolbox
%     Copyright (C) 2011  Paul Schmidt
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

    nx = size(img, 2);
    ny = size(img, 1);
    
    if neighbor_order == 0
       indx1 = [indx - ny, indx - 1, indx + 1, indx + ny];       
       indx_complete = [indx1];
    end
    
    if neighbor_order == 1
       indx1 = [indx - ny, indx - 1, indx + 1, indx + ny];
       indx2 = indx - nx * ny;
       indx3 = indx + nx * ny;
       indx_complete = [indx1, indx2, indx3];
    end
    
    if neighbor_order == 3
       indx1 = [indx - ny - 1, indx - ny, indx - ny + 1, indx - 1, indx + 1, indx + ny - 1, indx + ny + 1, indx + ny];
       indx2 = [indx - ny - 1, indx - ny, indx - ny + 1, indx - 1, indx + 1, indx + ny - 1, indx + ny + 1, indx + ny, indx] - nx * ny;
       indx3 = [indx - ny - 1, indx - ny, indx - ny + 1, indx - 1, indx + 1, indx + ny - 1, indx + ny + 1, indx + ny, indx] + nx * ny;
       indx_complete = [indx1, indx2, indx3];
    end
        
    indx_complete(indx_complete <= 0) = NaN;
    indx_complete(indx_complete > numel(img(:))) = NaN;
    neighborhood = img(indx_complete)';

end

