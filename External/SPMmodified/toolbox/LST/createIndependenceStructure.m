function [indi_struct] = createIndependenceStructure(nx, ny, nz, neighbor_order)

%     createIndependenceStructure.m - part of the LST toolbox
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

    if neighbor_order == 1
        indi_struct1 = repmat([1,2;2,1], [ceil(nx/2), ceil(ny/2)]);
        %indi_struct1 = indi_struct1(1:(end - 1),1:(end - 1));
        indi_struct1 = indi_struct1(1:(end),1:(end));
        indi_struct2 = repmat([2,1;1,2], [ceil(nx/2), ceil(ny/2)]);
        %indi_struct2 = indi_struct2(1:(end - 1),1:(end - 1));
        indi_struct2 = indi_struct2(1:(end),1:(end));
        
        indi_mat = zeros(nx, ny, nz);
        
        for i = 1:nz
            if(mod(i, 2) == 0)
                indi_mat(:,:,i) = indi_struct2(1:nx,1:ny);
            else 
                indi_mat(:,:,i) = indi_struct1(1:nx,1:ny);
            end
        end
        
        indi_struct = indi_mat;
        indi_pos1 = find(indi_mat == 1);
        indi_pos2 = find(indi_mat == 2);
    end
    
    if neighbor_order == 3
        indi_struct1 = repmat([1,3;2,4], nx, ny);
        indi_struct2 = repmat([5,7;6,8], nx, ny);
        indi_mat = zeros(nx, ny, nz);
        
        for i = 1:2:nz
             indi_mat(:,:,i) = indi_struct2(1:nx,1:ny);
        for i = 2:2:nz
             indi_mat(:,:,i) = indi_struct1(1:nx,1:ny);
        end
        indi_struct = indi_mat;
    end
          

end
