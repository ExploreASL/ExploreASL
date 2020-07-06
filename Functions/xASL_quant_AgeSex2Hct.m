function [Hematocrit] = xASL_quant_AgeSex2Hct(age, sex)
% xASL_quant_AgeSex2Hct Calculate estimated hematocrit (hct) based on age and sex
%
% FORMAT: [Hematocrit] = xASL_quant_AgeSex2Hct([age, sex])
% 
% INPUT:
%   age        - age in years, in any numerical format (OPTIONAL, defaulted (also by NaN) to average Hematocrit over lifespan for male/female)
%   sex        - 1=male, 2=female (OPTIONAL, defaulted (also by NaN) to average between sexes)
%
% OUTPUT:
%   Hematocrit - hematocrit as fraction (between 0 and 1)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function estimates a participants hematocrit, based on
% literature-based values for age and sex. It performs the following steps:
%
% 1. Warning unknown age/sex
% 2. Imputing unknown age/sex
% 3. Define hematocrit per age for unknown sex
% 4. Define Hematocrit per age for males
% 5. Define Hematocrit per age for females
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_quant_AgeSex2Hct(50, 2); % for a 50 years-old female
% REFERENCES: with courtesy from Patrick Hales:
%  9. Wu et al, In vivo venous blood T1 measurements ... MRM 2010
% 13. Lubin et al, Hematology of Infancy and Childhood 1981 (book chapter)
% 14. Devine et al, Mean blood hematocrit of adults, Natl Cent Health Stat
% Ser Public Health Serv Publ 1967
% __________________________________
% Copyright 2015-2020 ExploreASL

if nargin < 1
	age = [];
end

if nargin < 2 || isempty(sex)
	sex = nan(size(age));
end



    nHct = numel(age);
    Hematocrit = zeros(1,nHct);
    Hematocrit(:) = NaN; % start with missing values as default
    
    for h=1:nHct
        %% 1. Warning unknown age/sex
        if isempty(age(h)) || isnan(age(h))
            warning('Age unknown, imputing this');
            age(h) = NaN;
        elseif isempty(sex(h)) || isnan(sex(h))
            warning('Sex unknown, imputing this');
            sex(h) = NaN;
        end
        
        %% 2. Imputing unknown age/sex
        if isnan(age(h)) && isnan(sex(h))
            Hematocrit(h) = 0.44;      % mean adult value for m/f, Wu et al.
        end
        if isnan(age(h)) && ~isnan(sex(h))
            if sex(h)==1
                Hematocrit(h) = 0.47;      % mean adult male value, Wu et al.
            else
                Hematocrit(h) = 0.41;     % mean adult female value, Wu et al.
            end
        end
        %% 3. Define hematocrit per age for unknown sex
        if ~isnan(age(h)) && isnan(sex(h))
            % mean of male/female over age(h) ranges (from spreadsheet)
            if age(h) >0 && age(h)<=0.01
                Hematocrit(h) = 0.53;
            end
            if age(h) >0.01 && age(h)<=0.02
                Hematocrit(h) = 0.51;
            end
            if age(h) >0.03 && age(h)<=0.04
                Hematocrit(h) = 0.48;
            end
            if age(h) >0.04 && age(h)<=0.08
                Hematocrit(h) = 0.40;
            end
            if age(h) >0.08 && age(h)<=0.50
                Hematocrit(h) = 0.355;
            end
            if age(h) >0.50 && age(h)<=2.0
                Hematocrit(h) = 0.355;
            end
            if age(h) >2.0 && age(h)<=5.0
                Hematocrit(h) = 0.365;
            end
            if age(h) >5.0 && age(h)<=8.0
                Hematocrit(h) = 0.385;
            end
            if age(h) >8.0 && age(h)<=12.0
                Hematocrit(h) = 0.39;
            end
            if age(h) >12.0 && age(h)<=18.0
                Hematocrit(h) = 0.414;
            end
            % interpolated points between 18 and 21...
            if age(h) >18.0 && age(h)<=19.0
                Hematocrit(h) = 0.415;
            end
            if age(h) >19.0 && age(h)<=20.0
                Hematocrit(h) = 0.428;
            end
            if age(h) >20.0 && age(h)<=21.0
                Hematocrit(h) = 0.441;
            end
            %...end of interpolation
            if age(h) >21 && age(h)<=24
                Hematocrit(h) = 0.4405;
            end
            if age(h) >24 && age(h)<=34
                Hematocrit(h) = 0.444;
            end
            if age(h) >34 && age(h)<=44
                Hematocrit(h) = 0.443;
            end
            if age(h) >44 && age(h)<=54
                Hematocrit(h) = 0.445;
            end
            if age(h) >54 && age(h)<=64
                Hematocrit(h) = 0.449;
            end
            if age(h) >64 && age(h)<=74
                Hematocrit(h) = 0.445;
            end
            if (age(h) >74)
                Hematocrit(h) = 0.441;
            end
        end
        %% 4. Define Hematocrit per age for males
        if sex(h)==1  % male
            if age(h) >0 && age(h)<=0.01
                Hematocrit(h) = 0.53;
            end
            if age(h) >0.01 && age(h)<=0.02
                Hematocrit(h) = 0.51;
            end
            if age(h) >0.03 && age(h)<=0.04
                Hematocrit(h) = 0.48;
            end
            if age(h) >0.04 && age(h)<=0.08
                Hematocrit(h) = 0.40;
            end
            if age(h) >0.08 && age(h)<=0.50
                Hematocrit(h) = 0.355;
            end
            if age(h) >0.50 && age(h)<=2.0
                Hematocrit(h) = 0.355;
            end
            if age(h) >2.0 && age(h)<=5.0
                Hematocrit(h) = 0.365;
            end
            if age(h) >5.0 && age(h)<=8.0
                Hematocrit(h) = 0.385;
            end
            if age(h) >8.0 && age(h)<=12.0
                Hematocrit(h) = 0.39;
            end
            if age(h) >12.0 && age(h)<=18.0
                Hematocrit(h) = 0.415;
            end
            % interpolated points between 18 and 21...
            if age(h) >18.0 && age(h)<=19.0
                Hematocrit(h) = 0.415;
            end
            if age(h) >19.0 && age(h)<=20.0
                Hematocrit(h) = 0.432;
            end
            if age(h) >20.0 && age(h)<=21.0
                Hematocrit(h) = 0.450;
            end
            %...end of interpolation
            if age(h) >21 && age(h)<=24
                Hematocrit(h) = 0.467;
            end
            if age(h) >24 && age(h)<=34
                Hematocrit(h) = 0.470;
            end
            if age(h) >34 && age(h)<=44
                Hematocrit(h) = 0.466;
            end
            if age(h) >44 && age(h)<=54
                Hematocrit(h) = 0.465;
            end
            if age(h) >54 && age(h)<=64
                Hematocrit(h) = 0.461;
            end
            if age(h) >64 && age(h)<=74
                Hematocrit(h) = 0.458;
            end
            if age(h) >74
                Hematocrit(h) = 0.451;
            end
        else
            %% 5. Define Hematocrit per age for females
            if age(h) >0 && age(h)<=0.01
                Hematocrit(h) = 0.53;
            end
            if age(h) >0.01 && age(h)<=0.02
                Hematocrit(h) = 0.51;
            end
            if age(h) >0.03 && age(h)<=0.04
                Hematocrit(h) = 0.48;
            end
            if age(h) >0.04 && age(h)<=0.08
                Hematocrit(h) = 0.40;
            end
            if age(h) >0.08 && age(h)<=0.50
                Hematocrit(h) = 0.355;
            end
            if age(h) >0.50 && age(h)<=2.0
                Hematocrit(h) = 0.355;
            end
            if age(h) >2.0 && age(h)<=5.0
                Hematocrit(h) = 0.365;
            end
            if age(h) >5.0 && age(h)<=8.0
                Hematocrit(h) = 0.385;
            end
            if age(h) >8.0 && age(h)<=12.0
                Hematocrit(h) = 0.39;
            end
            if age(h) >12.0 && age(h)<=18.0
                Hematocrit(h) = 0.415;
            end
            if age(h) >18 && age(h)<=24
                Hematocrit(h) = 0.414;
            end
            if age(h) >24 && age(h)<=34
                Hematocrit(h) = 0.418;
            end
            if age(h) >34 && age(h)<=44
                Hematocrit(h) = 0.42;
            end
            if age(h) >44 && age(h)<=54
                Hematocrit(h) = 0.425;
            end
            if age(h) >54 && age(h)<=64
                Hematocrit(h) = 0.437;
            end
            if age(h) >64 && age(h)<=74
                Hematocrit(h) = 0.432;
            end
            if age(h) >74
                Hematocrit(h) = 0.431;
            end
        end % if sex(h)==1
    end % for h=1:nHct


end

% See example code for visualization below
% Age         = [0.005 0.015 0.035 0.06 0.29  1.25 3.5  6.5  10 15   22.5 29   39   49   59   69   79];
% HctMale     = [53    51     48    40   35.5 35.5 36.5 38.5 39 41.5 46.7 47   46.6 46.5 46.1 45.8 45.1];
% HctFemale   = [53    51     48    40   35.5 35.5 36.5 38.5 39 41.5 41.4 41.8 42   42.5 43.7 43.2 43.1];
% 
% figure(1);plot(Age,HctFemale,'b',Age,HctMale,'r')
% 
% xlabel('Age (years)');
% ylabel('Hematocrit (L/L, %)');
% title('Male (red) vs. female (blue) hematocrit');