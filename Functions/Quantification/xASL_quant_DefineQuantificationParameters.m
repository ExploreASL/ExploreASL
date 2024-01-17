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
% 1. Define arterial blood T1
% 2. Define arterial blood T2
% 3. Define tissue T1
% 4. Define tissue T2(*)
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
% Copyright (C) 2015-2024 ExploreASL


%% 0. Admin
if ~isfield(x, 'Q')
    x.Q = struct;
end

if ~isfield(x, 'MagneticFieldStrength') || isempty(x.MagneticFieldStrength)
    warning('MagneticFieldStrength was not defined, defaulting to 3T');
    x.MagneticFieldStrength = 3;
end


%% 1. Define arterial blood T1
if ~isfield(x.Q,'BloodT1') || isempty(x.Q.BloodT1)
    % T1 relaxation time of arterial blood
    % There are 3 options for x.Q.BloodT1:
    % A) users have provided x.Q.BloodT1
    % B) users have provided x.Hematocrit which is converted to x.Q.BloodT1 above
    % C) it doesn't exist and is defaulted here based on MagneticFieldStrength
    switch(x.MagneticFieldStrength)
	    case 0.2 
		    x.Q.BloodT1 = 776; % Rooney 2007 MRM
            fprintf('%s\n', 'Defaulting x.Q.BloodT1 to 776 ms for 0.2T (Rooney 2007 MRM)');
	    case 1
		    x.Q.BloodT1 = 1350; % Rooney 2007 MRM
            fprintf('%s\n', 'Defaulting x.Q.BloodT1 to 1350 ms for 1T (Rooney 2007 MRM)');
	    case 1.5
		    x.Q.BloodT1 = 1540; % Rooney 2007 MRM
            fprintf('%s\n', 'Defaulting x.Q.BloodT1 to 1540 ms for 1.5T (Rooney 2007 MRM)');
	    case 3
		    x.Q.BloodT1 = 1650; % Alsop 2015 MRM
            fprintf('%s\n', 'Defaulting x.Q.BloodT1 to 1650 ms for 3T (Alsop 2015 MRM)');
	    case 4
		    x.Q.BloodT1 = 1914; % Rooney 2007
            fprintf('%s\n', 'Defaulting x.Q.BloodT1 to 1914 ms for 4T (Rooney 2007 MRM)');
	    case 7
		    %x.Q.BloodT1 = 2578; % Rooney 2007 MRM
		    x.Q.BloodT1 = 2100; % Ivanov 2017 NeuroImage
            fprintf('%s\n', 'Defaulting x.Q.BloodT1 to 2100 ms for 7T (Ivanov 2007 NeuroImage)');
	    otherwise
		    x.Q.BloodT1 = 1650; % Alsop 2015 MRM - assuming default 3 T
		    fprintf('%s\n',['Warning: Unknown T1-blood for ' xASL_num2str(x.MagneticFieldStrength) 'T scanner, using 3T value (Alsop 2015 MRM)']);
            % PM: NOTE that this situation is unlikely, given that we
            % default to x.MagneticFieldStrength = 3 at section 0
            % Administration above
    end
end


%% 2. Define arterial blood T2
if ~isfield(x.Q, 'T2art')
	if x.MagneticFieldStrength == 3
		x.Q.T2art = 165; % ms Gregori JMRI 2013; Lee ISMRM 2003
		% Jean Chen 2009 MRM, DOI: 10.1002/mrm.21858 175 ms
		% 175 ms for Hct 0.21; 122 ms for Hct 0.44
	else
		x.Q.T2art = 239; % ms Lee ISMRM 2003
		% Jean Chen 2009 MRM, DOI: 10.1002/mrm.21858 157 ms
	end
end
if ~isfield(x.Q,'Lambda')
    x.Q.Lambda = 0.9; % Brain/blood water coefficient (mL 1H/ mL blood)
end


%% 3. Define tissue T1
if ~isfield(x.Q,'TissueT1')
	switch(x.MagneticFieldStrength)
		% T1 GM tissue
		case 0.2
			x.Q.TissueT1 =  635; % Rooney 2007
		case 1
			x.Q.TissueT1 = 1036; % Rooney 2007
		case 1.5
			x.Q.TissueT1 = 1188; % Rooney 2007
		case 3
			x.Q.TissueT1 = 1240; % Alsop 2015
		case 4
			x.Q.TissueT1 = 1530; % Rooney 2000
		case 7
			x.Q.TissueT1 = 1920; % Marques 2010
		otherwise
			x.Q.TissueT1 = 1240;
			fprintf('%s\n',['Warning: Unknown T1 GM for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
	end
end


%% 4. Define tissue T2(*)
if ~isfield(x.Q,'T2star') || isempty(x.Q.T2star)
    if x.MagneticFieldStrength == 3
	    x.Q.T2star = 47.3; % default for 3T; Lu and van Zijl, MRM 2005, DOI: 10.1002/mrm.20379
    elseif x.MagneticFieldStrength == 7
	    x.Q.T2star = 35.6; % Voelker 2021
    elseif x.MagneticFieldStrength == 1.5
	    x.Q.T2star = 62.0; % Lu and van Zijl, MRM 2005, DOI: 10.1002/mrm.20379
    else
	    x.Q.T2star = 47.3;
	    fprintf('%s\n',['Warning: Unknown T2star for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
    end
end
if ~isfield(x.Q,'T2') || isempty(x.Q.T2)
    if x.MagneticFieldStrength == 3
	    x.Q.T2 = 85; % in ms - default for 3T (ref Johannes Gregori, JMRI 2013) 88 for frontal GM, 79 for occipital GM (Lu et al, 2005 JMRI)
	    % Hct specific values are in 10.1002/mrm.21342
    elseif x.MagneticFieldStrength == 1.5
	    x.Q.T2 = 95; % in ms - 99 for frontal GM, 90 for occipital GM (Lu et al, 2005 JMRI).
    else
	    x.Q.T2 = 85;
	    fprintf('%s\n',['Warning: Unknown T2 for ' num2str(x.MagneticFieldStrength) 'T scanners, using 3T value']);
    end
end


end