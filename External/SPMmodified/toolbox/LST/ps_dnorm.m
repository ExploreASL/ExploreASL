function out = dnorm(x, mu, sigma)

%     dnorm.m - Normal probability density function, 
%     part of the LST toolbox
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

lf = - 0.5 .* log(2 * pi .* sigma .* sigma) - 0.5 .* (sigma .* sigma).^(-1) .* (x - mu).*(x - mu);
out = exp(lf);

end

