function BloodT1 = xASL_quant_Hct2BloodT1(Hematocrit, Y, B0, bVerbose)
% xASL_quant_Hct2BloodT1 Predict blood T1 from venous (antecubital) hematocrit
%
% FORMAT: BloodT1 = xASL_quant_Hct2BloodT1(Hematocrit, Y, B0, bVerbose)
%
% INPUT:
%   Hematocrit  - (antecubital) venous hematocrit, as fraction or percentage (REQUIRED)
%   Y           - O2 saturation (fraction or percentage) (OPTIONAL, DEFAULT=0.97)
%   B0          - magnetic field strength used in Tesla (OPTIONAL, DEFAULT=3)
%   bVerbose    - boolean specifying if we want output to the screen (OPTIONAL, DEFAULT=true)
%
% OUTPUT: BloodT1 - longitudinal relaxation time of blood (ms)
% --------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts hematocrit to blood T1, according to
%              calculations defined by Patrick Hales. With courtesy and thanks!
%              Note that we assume a venous O2 saturation of 68% (Yv=0.68)
%
%              This function performs the following steps:
%              1) Check fraction vs percentage hematocrit & Y, should be between 0 and 1
%              2) Specify defaults (Hb, Fe)
%              3) Perform calculation
%              4) Convert s to ms
%              5) Print what we did
% --------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_wrp_Quantify(x);
% REFERENCE: Hales, 2014 JCBFM
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

%% ---------------------------------------------------------
%% Admin
if nargin<44 || isempty(bVerbose)
    bVerbose = true;
end
if nargin<3 || isempty(B0)
    B0 = 3; % 3T field strength
elseif length(B0)~=1
    warning('Incorrect B0 specified, skipping');
    return;
end
if nargin<2 || isempty(Y)
    Y = 0.97;
elseif length(Y)~=1
    warning('Incorrect Y specified, skipping');
    return;
end
if nargin<1 || isempty(Hematocrit) || length(Hematocrit)~=1
    warning('Incorrect Hematocrit specified, skipping');
    return;
end


%% ---------------------------------------------------------
%% 1) Check fraction vs percentage hematocrit & Y, should be between 0 and 1
if Hematocrit>0 && Hematocrit<1
    % this is fine, just continue
elseif Hematocrit>1 && Hematocrit<100
    % we have percentages, correct
    Hematocrit = Hematocrit/100;
else
    warning('Hematocrit was defined in incorrect range, please correct! Skipping...');
    return;
end

if Y>0 && Y<1
    % this is fine, just continue
elseif Y>1 && y<100
    % we have percentages, correct
    Y = Y/100;
else
    warning('Incorrect Y range, skipping');
    return;
end


%% ---------------------------------------------------------
%% 2) Specify defaults (Hb, Fe)
Hb = 5.15;  % mean corpuscular haemoglobin concentration (mmol Hb tetramer / L plasma)
Fe = (0.70*Hematocrit)/((0.70*Hematocrit)+(0.95*(1-Hematocrit))); % (is this iron content?)


%% ---------------------------------------------------------
%% 3) Perform calculation
part1 = 1.099 - (0.057*B0) + ((0.033*Hb)*(1-Y));
part2 = (1-Fe)*(0.496-(0.023*B0));

BloodT1 = (1/((Fe*part1)+part2)) + 0.108;  % (seconds)
% is this offset for venous to arterial?


%% ---------------------------------------------------------
%% 4) Convert s to ms
BloodT1 = BloodT1*1000;

%% ---------------------------------------------------------
%% 5) Print what we did
if bVerbose
    fprintf('%s\n', ['Calculated blood T1 ' num2str(BloodT1) 'ms from Hct ' num2str(Hematocrit)]);
    fprintf('%s\n', ['Assuming field strength ' xASL_num2str(B0) ', arterial O2 saturation ' xASL_num2str(Y)]);
end


end