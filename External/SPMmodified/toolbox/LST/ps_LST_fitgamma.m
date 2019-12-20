function [a, b] = fitgamma(y)

%     fitgamma.m - Maximum likelihood estimation for the parameters of the
%     gamma distribution 
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

mean_log_y = mean(log(y));
log_mean_y = log(mean(y));

% starting point for a
a = 0.5 / (log_mean_y - mean_log_y);
a_old = 0;
ii = 0;

while(abs(a_old - a) > 0.0001)
    ii = ii + 1;
    a_old = a;
    a = 1 / (1/a_old + (mean_log_y - log_mean_y + log(a_old) - psi(a_old)) / (a_old * a_old * (1 / a_old - psi(1, a_old))));
end
b = a / mean(y);

end