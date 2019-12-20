function out = dgamma(y, a, b)

%     dgamma.m - Probability density function for the gamma distribution, 
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

log_y = log(y);
out = a * log(b) - gammaln(a) + (a - 1) * log_y - b .* y;
out = exp(out);

end