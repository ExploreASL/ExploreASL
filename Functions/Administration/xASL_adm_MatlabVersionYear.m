function matlab_version_year = xASL_adm_MatlabVersionYear()
% Extracts the release year from the current MATLAB version.
%
% FORMAT: matlab_version_year = xASL_adm_MatlabVersionYear
%
% INPUT:
%   None
%
% OUTPUT:
%   matlab_version_year - MATLAB release year as a numeric scalar
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Extracts the four-digit release year from the MATLAB version string.
%              For example, both {{R2025a}} and {{R2025b}} return the numeric value 2025.
%              The function is independent of MATLAB's numeric versioning scheme and therefore
%              works with releases before and after the version-numbering change in R2023b.
% EXAMPLE: matlab_version_year = xASL_adm_MatlabVersionYear;
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________

% Get the MATLAB version string
matlab_version = version;

% Extract the four-digit release year
matlab_version_token = regexp(matlab_version, '\(R(\d{4})[ab]\)', 'tokens', 'once');

% Check if the release year was detected
if isempty(matlab_version_token)
	error('xASL_adm_MatlabVersionYear:Unable to determine the MATLAB release year from version string [%s]', matlab_version);
end

% Convert the release year to a numeric value
matlab_version_year = str2double(matlab_version_token{1});
end