function [OutputN] = xASL_round(InputN, PrecisionN)
%xASL_round Wrapper around Matlab's 'round' for older Matlab versions
%
% FORMAT: [OutputN] = xASL_round(InputN[, PrecisionN])
%
% INPUT:
%   InputN 	- input number (or array of numbers) to be rounded (REQUIRED)
%   PrecisionN  - number of decimals (OPTIONAL, DEFAULT = 0)
%
% OUTPUT:
%   OutputN   - output number (or array of numbers) that are rounded
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Recent Matlab versions support a second input that specifies that number of decimals to round at,
% but earlier Matlab versions do not support this. For backward compatibility, use this wrapper instead of round.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: 
% xASL_round(10.5897, 2)
% ans = 10.5900
% 
% xASL_round(10.5897)
% ans = 11
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE

if nargin<2 || isempty(PrecisionN)
    PrecisionN = 0;
end

PrecisionF      = 10^PrecisionN;
OutputN         = round(InputN.*PrecisionF)./PrecisionF;


end

