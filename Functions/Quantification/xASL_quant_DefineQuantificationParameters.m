function [x] = xASL_quant_DefineQuantificationParameters(x)
%xASL_quant_DefineQuantificationParameters Central function for defining ASL quantification parameters
%
% FORMAT: [x] = xASL_quant_DefineQuantificationParameters(x)
%
% INPUT:
%   x       - structure containing fields with all information required to run this function (REQUIRED)
%             with the following information:
% OUTPUT: (note that these are [] empty by default, unless they are created)
%   x       - structure containing fields with all information required to run this function (REQUIRED)
%             with the following information:
%
% DESCRIPTION: This function defines the following ASL quantification parameters
% 1.   Hematocrit
% 2.   Arterial blood T1
% 3.   Arterial blood T2
% 4.   GM and WM T1
% 5.   GM and WM T2(*) and T2 of extravascular compartment for multiTE fitting
%
% REFERENCES: 
%     Gregori, Johannes et al. “T2-based arterial spin labeling measurements of blood to tissue water transfer 
%     in human brain." Journal of magnetic resonance imaging : JMRI vol. 37,2 (2013): 332-42. doi:10.1002/jmri.23822
%  
%     Lee T, Stainsby JA, Hong J, Han E, Brittain J, Wright GA. Blood Relaxation Properties at 3T --Effects 
%     of Blood Oxygen Saturation. Proc Intl Soc Mag Reson Med. 11:131.
%
%     Rooney WD, Johnson G, Li X, Cohen ER, Kim SG, Ugurbil K, Springer Jr CS. 
%     Magnetic field and tissue dependencies of human brain longitudinal 1H2O relaxation in vivo. 
%     Magnetic Resonance in Medicine. 2007 Feb;57(2):308-18.
%     
%     Rooney WD, Lee JH, Li X, Wang GJ, Franceschi D, Springer CS, Volkow ND. 
%     4.0T water proton T1 relaxation times in normal human brain and during acute ethanol intoxication. 
%     Alcohol Clin Exp Res 2000; 24: 830-836.
%   
%     Voelker MN, Kraff O, Goerke S, Laun FB, Hanspach J, Pine KJ, Ehses P, Zaiss M, Liebert A, Straub S, Eckstein K. 
%     The traveling heads 2.0: Multicenter reproducibility of quantitative imaging methods at 7 Tesla. 
%     NeuroImage. 2021 May 15;232:117910.
%
%     Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele PF, Gruetter R. 
%     MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. 
%     Neuroimage. 2010 Jan 15;49(2):1271-81. 
%    
%     Ivanov D, Gardumi A, Haast RAM, Pfeuffer J, Poser BA, Uludag K. 
%     Comparison of 3 T and 7 T ASL techniques for concurrent functional perfusion and BOLD studies
%     Neuroimage. 2017; 156:363-376.
% 
% EXAMPLE: x = xASL_quant_DefineQuantificationParameters(x);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


%% ------------------------------------------------------------------------------------------------
%% 0.   Admin
if ~isfield(x, 'Q')
    x.Q = struct;
end

if ~isfield(x, 'MagneticFieldStrength') || isempty(x.MagneticFieldStrength)
    warning('MagneticFieldStrength was not defined, defaulting to 3T');
    x.MagneticFieldStrength = 3;
end


%% ------------------------------------------------------------------------------------------------
%% 1.   Hematocrit
% Here, we check if the user has provided a hematocrit value for this
% subject_session_run. Only then, we create a x.Q.T1blood.
% Below, at the quantification section, this is only taken into account
% when x.Q.BloodT1 exists, otherwise default Blood T1 values are used based
% on MagneticFieldStrength.

IndexSetsAge = find(strcmpi(x.S.SetsName, 'age'));
IndexSetsSex = find(strcmpi(x.S.SetsName, 'sex'));
indexSetsHct = find(strcmpi(x.S.SetsName, 'hematocrit'));

%% a. We prioritize participants.tsv>Hematocrit (x.S.SetsID > x.Hematocrit)
if ~isempty(indexSetsHct)
    x.Q.Hematocrit = x.S.SetsID(x.iSubjectSession, indexSetsHct);
    fprintf('%s\n', 'Using hematocrit found in participants.tsv for blood T1 correction');
end

if isfield(x, 'Hematocrit') || isfield(x, 'hematocrit')
    warning('x.Hematocrit detected, we ignore this, considering adding hct values to participants.tsv instead');
end


%% a2. Manage hematocrit usage parameter
if isfield(x.modules.asl, 'bHct2BloodT1')
    if x.modules.asl.bHct2BloodT1 == 2 && isfield(x.Q, 'Hematocrit')
        warning('Parameter x.modules.asl.bHct2BloodT1 was set to 2: trying to infer hematocrit from age and sex, but hematocrit data were also found');
        fprintf('%s\n', 'Consider setting x.modules.asl.bHct2BloodT1 to 1 to use the hematocrit data directly');
    end
    if x.modules.asl.bHct2BloodT1 == 2 && (isempty(IndexSetsAge) || isempty(IndexSetsSex))
        warning('Parameter x.modules.asl.bHct2BloodT1 was set to 2: trying to infer hematocrit from age and sex, but age & sex data were incomplete');
        fprintf('%s\n', 'Consider changing x.modules.asl.bHct2BloodT1 or ensure that age & sex data are present in participants.tsv');
        x.modules.asl.bHct2BloodT1 = []; % Setting this to empty here, so it will be dealt with below
    end
end

if ~isfield(x.modules.asl, 'bHct2BloodT1') || isempty(x.modules.asl.bHct2BloodT1)
    if isfield(x.Q, 'Hematocrit')
        warning('Parameter x.modules.asl.bHct2BloodT1 was not set but Hematocrit data were found');
        fprintf('%s\n', 'Setting x.modules.asl.bHct2BloodT1 to 1: converting hematocrit to blood T1 values');
        fprintf('%s\n', 'Consider setting x.modules.asl.bHct2BloodT1 to avoid this warning');
        x.modules.asl.bHct2BloodT1 = 1;
    elseif ~isempty(IndexSetsAge) || ~isempty(IndexSetsSex)
        warning('Parameter x.modules.asl.bHct2BloodT1 was not set but age & sex data were found');
        fprintf('%s\n', 'Setting x.modules.asl.bHct2BloodT1 to option 2: trying to infer hematocrit from age & sex');
        fprintf('%s\n', 'Consider setting x.modules.asl.bHct2BloodT1 to avoid this warning');
        x.modules.asl.bHct2BloodT1 = 2;
    else
        fprintf('\n%s\n', 'No hematocrit data found, disabling x.modules.asl.bHct2BloodT1');
        x.modules.asl.bHct2BloodT1 = 0;
    end
end


%% b. We model the expected hematocrit from age & sex
if x.modules.asl.bHct2BloodT1 == 2
    fprintf('%s\n', 'Trying to infer hematocrit from age & sex');
    
    age = x.S.SetsID(:, IndexSetsAge);
    
    % Convert sex correctly
    sex = x.S.SetsID(:, IndexSetsSex);
    sexOptions = x.S.SetsOptions{:, IndexSetsSex};
    sexN = nan(length(sex), 1); % currently we can only infer hct from male/female
    % So anything that is not detected, will remain NaNs
    for iOption=1:length(sexOptions)
        if ~isempty(regexpi(sexOptions{iOption}, '^(male|m|man|men)$'))
            sexN(sex==iOption) = 1;
        elseif ~isempty(regexpi(sexOptions{iOption}, '^(female|f|woman|women)$'))
            sexN(sex==iOption) = 2;
        end
    end

    if sum(isnan(sexN))>0
        warning('Unknown sex detected in participants.tsv, currently we can only infer hematocrit from male or female');
        fprintf('%s\n', 'So in participants.tsv specify the words "male" and "female" only');
    end

    % CAVE: here sex 1 ==male 2 ==female
    x.Q.Hematocrit = xASL_quant_AgeSex2Hct(age(x.iSubjectSession), sex(x.iSubjectSession));
    % PM: model hematocrit based on age, sex, ethnicity
end


%% ------------------------------------------------------------------------------------------------
%% 2.   Arterial blood T1
% T1 relaxation time of arterial blood
% There are 3 options for x.Q.T1blood (A has the highest priority:
% A) users have provided x.Q.T1blood
% B) users have provided x.Hematocrit (in any of the forms defined in xASL_wrp_Quantify 3.a-c), which is converted to x.Q.T1blood there
% C) it doesn't exist and is defaulted here based on MagneticFieldStrength

% We convert x.Q.Hematocrit -> x.Q.T1blood
if isfield(x.Q, 'Hematocrit') && (~isfield(x.Q, 'T1blood') || isempty(x.Q.T1blood))
    x.Q.T1blood = xASL_quant_Hct2BloodT1(x.Q.Hematocrit, [], x.MagneticFieldStrength);
end

if ~isfield(x.Q, 'T1blood') || isempty(x.Q.T1blood)
	switch(x.MagneticFieldStrength)
		case 0.2
		    x.Q.T1blood = 776; % Rooney 2007 MRM
            fprintf('%s\n', 'Defaulting x.Q.T1blood to 776 ms for 0.2T (Rooney 2007 MRM)');
	    case 1
		    x.Q.T1blood = 1350; % Rooney 2007 MRM
            fprintf('%s\n', 'Defaulting x.Q.T1blood to 1350 ms for 1T (Rooney 2007 MRM)');
	    case 1.5
		    x.Q.T1blood = 1540; % Rooney 2007 MRM
            fprintf('%s\n', 'Defaulting x.Q.T1blood to 1540 ms for 1.5T (Rooney 2007 MRM)');
	    case 3
		    x.Q.T1blood = 1650; % Alsop 2015 MRM
            fprintf('%s\n', 'Defaulting x.Q.T1blood to 1650 ms for 3T (Alsop 2015 MRM)');
	    case 4
		    x.Q.T1blood = 1914; % Rooney 2007
            fprintf('%s\n', 'Defaulting x.Q.T1blood to 1914 ms for 4T (Rooney 2007 MRM)');
	    case 7
		    %x.Q.T1blood = 2578; % Rooney 2007 MRM
		    x.Q.T1blood = 2100; % Ivanov 2017 NeuroImage
            fprintf('%s\n', 'Defaulting x.Q.T1blood to 2100 ms for 7T (Ivanov 2007 NeuroImage)');
	    otherwise
		    x.Q.T1blood = 1650; % Alsop 2015 MRM - assuming default 3 T
		    fprintf('%s\n',['Warning: Unknown T1-blood for ' xASL_num2str(x.MagneticFieldStrength) 'T scanner, using 3T value (Alsop 2015 MRM)']);
            % PM: NOTE that this situation is unlikely, given that we
            % default to x.MagneticFieldStrength = 3 at section 0
            % Administration above
    end
end


%% ------------------------------------------------------------------------------------------------
%% 3.   Arterial blood T2
if ~isfield(x.Q, 'T2art')
	switch(x.MagneticFieldStrength)
		case 3
			x.Q.T2art = 165; % ms Gregori JMRI 2013; Lee ISMRM 2003
			% Jean Chen 2009 MRM, DOI: 10.1002/mrm.21858 175 ms
			% 175 ms for Hct 0.21; 122 ms for Hct 0.44
		case 1.5
			x.Q.T2art = 239; % ms Lee ISMRM 2003
			% Jean Chen 2009 MRM, DOI: 10.1002/mrm.21858 157 ms
		otherwise
			x.Q.T2art = 165;
			fprintf('%s\n',['Warning: Unknown T2 blood for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);	
	end
end
if ~isfield(x.Q,'Lambda')
    x.Q.Lambda = 0.9; % Brain/blood water coefficient (mL 1H/ mL blood)
end


%% ------------------------------------------------------------------------------------------------
%% 4.   Tissue T1
if ~isfield(x.Q,'T1GM')
	switch(x.MagneticFieldStrength)
		% T1 GM tissue
		case 0.2
			x.Q.T1GM =  635; % Rooney 2007
		case 1
			x.Q.T1GM = 1036; % Rooney 2007
		case 1.5
			x.Q.T1GM = 1188; % Rooney 2007
		case 3
			x.Q.T1GM = 1240; % Alsop 2015
		case 4
			x.Q.T1GM = 1723; % Rooney 2007
		case 7
			x.Q.T1GM = 1920; % Marques 2010
		otherwise
			x.Q.T1GM = 1240;
			fprintf('%s\n',['Warning: Unknown T1 GM for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
	end
end

if ~isfield(x.Q,'T1WM')
	switch(x.MagneticFieldStrength)
		% T1 WM tissue
		case 0.2
			x.Q.T1WM =  361; % Rooney 2007
		case 1
			x.Q.T1WM = 555; % Rooney 2007
		case 1.5
			x.Q.T1WM = 656; % Rooney 2007
		case 3
			x.Q.T1WM = 800; % average of frontal & occipital WM from Lu et al., JMRI 2005 & 3 studies they refer to
		case 4
			x.Q.T1WM = 1010; % Rooney 2007
		case 7
			x.Q.T1WM = 1150; % Marques 2010
		otherwise
			x.Q.T1WM = 800;
			fprintf('%s\n',['Warning: Unknown T1 WM for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
	end
end


%% ------------------------------------------------------------------------------------------------
%% 5.   Tissue T2(*) in GM
if ~isfield(x.Q,'T2starGM') || isempty(x.Q.T2starGM)
    switch(x.MagneticFieldStrength)
		case 3
			x.Q.T2starGM = 47.3; % default for 3T; Lu and van Zijl, MRM 2005, DOI: 10.1002/mrm.20379
		case 7
			x.Q.T2starGM = 35.6; % Voelker 2021
		case 1.5
			x.Q.T2starGM = 62.0; % Lu and van Zijl, MRM 2005, DOI: 10.1002/mrm.20379
		otherwise
			x.Q.T2starGM = 47.3;
			fprintf('%s\n',['Warning: Unknown T2starGM for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
    end
end
if ~isfield(x.Q,'T2tissueMultiTE') || isempty(x.Q.T2tissueMultiTE)
	% T2tissueMultiTE is used for 2-compartment fitting with mutli-TE acquisition. By default, we assume GM
    switch(x.MagneticFieldStrength)
		case 3
			x.Q.T2tissueMultiTE = 85; % in ms - default for 3T (ref Johannes Gregori, JMRI 2013) 88 for frontal GM, 79 for occipital GM (Lu et al, 2005 JMRI)
			% Hct specific values are in 10.1002/mrm.21342
		case 1.5
			x.Q.T2tissueMultiTE = 95; % in ms - 99 for frontal GM, 90 for occipital GM (Lu et al, 2005 JMRI).
		otherwise
			x.Q.T2tissueMultiTE = 85;
			fprintf('%s\n',['Warning: Unknown T2tissueMultiTE for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
    end
end

if ~isfield(x.Q,'T2GM') || isempty(x.Q.T2GM)
    switch(x.MagneticFieldStrength)
		case 3
			x.Q.T2GM = 85; % in ms - default for 3T (ref Johannes Gregori, JMRI 2013) 88 for frontal GM, 79 for occipital GM (Lu et al, 2005 JMRI)
			% Hct specific values are in 10.1002/mrm.21342
		case 1.5
			x.Q.T2GM = 95; % in ms - 99 for frontal GM, 90 for occipital GM (Lu et al, 2005 JMRI).
		otherwise
			x.Q.T2GM = 85;
			fprintf('%s\n',['Warning: Unknown T2GM for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
    end
end

if ~isfield(x.Q,'T2WM') || isempty(x.Q.T2WM)
    switch(x.MagneticFieldStrength)
		case 3
			x.Q.T2WM = 75; % in ms - 69 for frontal WM, 81 for occipital WM (Lu et al, 2005 JMRI)
		case 1.5
			x.Q.T2WM = 86; % in ms - 79 for frontal WM, 92 for occipital WM (Lu et al, 2005 JMRI)
		otherwise
			x.Q.T2WM = 75;
			fprintf('%s\n',['Warning: Unknown T2WM for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
    end
end


end