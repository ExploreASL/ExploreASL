function Hct = xASL_quant_AgeSex2Hct(age, gender)
% xASL_quant_AgeSex2Hct Function to return estimated Hct, based on age and gender
% Enter age in years, and for gender: 0=female, 1=male
% pass NaN is either not known

% References from paper Patrick Hales:
%  9. Wu et al, In vivo venous blood T1 measurements ... MRM 2010
% 13. Lubin et al, Hematology of Infancy and Childhood 1981 (book chapter)
% 14. Devine et al, Mean blood hematocrit of adults, Natl Cent Health Stat
% Ser Public Health Serv Publ 1967


    nHct = numel(age);
    Hct = zeros(1,nHct);
    
    for h=1:nHct
     
        if (isnan(age(h)) && isnan(gender(h)))
            Hct(h) = 0.44;      % mean adult value for m/f, Wu et al.
        end
        if (isnan(age(h)) && ~isnan(gender(h)))
            if (gender(h)==1)
                Hct(h) = 0.47;      % mean adult male value, Wu et al.
            else
                Hct(h) = 0.41;     % mean adult female value, Wu et al.
            end
        end
        if (~isnan(age(h)) && isnan(gender(h)))
            % mean of male/female over age(h) ranges (from spreadsheet)
            if (age(h) >0 && age(h)<=0.01)
                Hct(h) = 0.53;
            end
            if (age(h) >0.01 && age(h)<=0.02)
                Hct(h) = 0.51;
            end
            if (age(h) >0.03 && age(h)<=0.04)
                Hct(h) = 0.48;
            end
            if (age(h) >0.04 && age(h)<=0.08)
                Hct(h) = 0.40;
            end
            if (age(h) >0.08 && age(h)<=0.50)
                Hct(h) = 0.355;
            end
            if (age(h) >0.50 && age(h)<=2.0)
                Hct(h) = 0.355;
            end
            if (age(h) >2.0 && age(h)<=5.0)
                Hct(h) = 0.365;
            end
            if (age(h) >5.0 && age(h)<=8.0)
                Hct(h) = 0.385;
            end
            if (age(h) >8.0 && age(h)<=12.0)
                Hct(h) = 0.39;
            end
            if (age(h) >12.0 && age(h)<=18.0)
                Hct(h) = 0.414;
            end
            % interpolated points between 18 and 21...
            if (age(h) >18.0 && age(h)<=19.0)
                Hct(h) = 0.415;
            end
            if (age(h) >19.0 && age(h)<=20.0)
                Hct(h) = 0.428;
            end
            if (age(h) >20.0 && age(h)<=21.0)
                Hct(h) = 0.441;
            end
            %...end of interpolation
            if (age(h) >21 && age(h)<=24)
                Hct(h) = 0.4405;
            end
            if (age(h) >24 && age(h)<=34)
                Hct(h) = 0.444;
            end
            if (age(h) >34 && age(h)<=44)
                Hct(h) = 0.443;
            end
            if (age(h) >44 && age(h)<=54)
                Hct(h) = 0.445;
            end
            if (age(h) >54 && age(h)<=64)
                Hct(h) = 0.449;
            end
            if (age(h) >64 && age(h)<=74)
                Hct(h) = 0.445;
            end
            if (age(h) >74)
                Hct(h) = 0.441;
            end

        end
        %% male
        if (~isnan(age(h)) && ~isnan(gender(h)))
            if (gender(h)==1)  % male
                if (age(h) >0 && age(h)<=0.01)
                    Hct(h) = 0.53;
                end
                if (age(h) >0.01 && age(h)<=0.02)
                    Hct(h) = 0.51;
                end
                if (age(h) >0.03 && age(h)<=0.04)
                    Hct(h) = 0.48;
                end
                if (age(h) >0.04 && age(h)<=0.08)
                    Hct(h) = 0.40;
                end
                if (age(h) >0.08 && age(h)<=0.50)
                    Hct(h) = 0.355;
                end
                if (age(h) >0.50 && age(h)<=2.0)
                    Hct(h) = 0.355;
                end
                if (age(h) >2.0 && age(h)<=5.0)
                    Hct(h) = 0.365;
                end
                if (age(h) >5.0 && age(h)<=8.0)
                    Hct(h) = 0.385;
                end
                if (age(h) >8.0 && age(h)<=12.0)
                    Hct(h) = 0.39;
                end
                if (age(h) >12.0 && age(h)<=18.0)
                    Hct(h) = 0.415;
                end
                % interpolated points between 18 and 21...
                if (age(h) >18.0 && age(h)<=19.0)
                    Hct(h) = 0.415;
                end
                if (age(h) >19.0 && age(h)<=20.0)
                    Hct(h) = 0.432;
                end
                if (age(h) >20.0 && age(h)<=21.0)
                    Hct(h) = 0.450;
                end
                %...end of interpolation
                if (age(h) >21 && age(h)<=24)
                    Hct(h) = 0.467;
                end
                if (age(h) >24 && age(h)<=34)
                    Hct(h) = 0.470;
                end
                if (age(h) >34 && age(h)<=44)
                    Hct(h) = 0.466;
                end
                if (age(h) >44 && age(h)<=54)
                    Hct(h) = 0.465;
                end
                if (age(h) >54 && age(h)<=64)
                    Hct(h) = 0.461;
                end
                if (age(h) >64 && age(h)<=74)
                    Hct(h) = 0.458;
                end
                if (age(h) >74)
                    Hct(h) = 0.451;
                end
            else
                %% female
                % females
                if (age(h) >0 && age(h)<=0.01)
                    Hct(h) = 0.53;
                end
                if (age(h) >0.01 && age(h)<=0.02)
                    Hct(h) = 0.51;
                end
                if (age(h) >0.03 && age(h)<=0.04)
                    Hct(h) = 0.48;
                end
                if (age(h) >0.04 && age(h)<=0.08)
                    Hct(h) = 0.40;
                end
                if (age(h) >0.08 && age(h)<=0.50)
                    Hct(h) = 0.355;
                end
                if (age(h) >0.50 && age(h)<=2.0)
                    Hct(h) = 0.355;
                end
                if (age(h) >2.0 && age(h)<=5.0)
                    Hct(h) = 0.365;
                end
                if (age(h) >5.0 && age(h)<=8.0)
                    Hct(h) = 0.385;
                end
                if (age(h) >8.0 && age(h)<=12.0)
                    Hct(h) = 0.39;
                end
                if (age(h) >12.0 && age(h)<=18.0)
                    Hct(h) = 0.415;
                end
                if (age(h) >18 && age(h)<=24)
                    Hct(h) = 0.414;
                end
                if (age(h) >24 && age(h)<=34)
                    Hct(h) = 0.418;
                end
                if (age(h) >34 && age(h)<=44)
                    Hct(h) = 0.42;
                end
                if (age(h) >44 && age(h)<=54)
                    Hct(h) = 0.425;
                end
                if (age(h) >54 && age(h)<=64)
                    Hct(h) = 0.437;
                end
                if (age(h) >64 && age(h)<=74)
                    Hct(h) = 0.432;
                end
                if (age(h) >74)
                    Hct(h) = 0.431;
                end
            end
        end
        
    end


end


% Age         = [0.005 0.015 0.035 0.06 0.29  1.25 3.5  6.5  10 15   22.5 29   39   49   59   69   79];
% HctMale     = [53    51     48    40   35.5 35.5 36.5 38.5 39 41.5 46.7 47   46.6 46.5 46.1 45.8 45.1];
% HctFemale   = [53    51     48    40   35.5 35.5 36.5 38.5 39 41.5 41.4 41.8 42   42.5 43.7 43.2 43.1];
% 
% figure(1);plot(Age,HctFemale,'b',Age,HctMale,'r')
% 
% xlabel('Age (years)');
% ylabel('Hematocrit (L/L, %)');
% title('Male (red) vs. female (blue) hematocrit');