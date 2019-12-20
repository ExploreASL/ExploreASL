function t1bl = xASL_quant_Hct2BloodT1(Hct, Y, B0)
% xASL_quant_Hct2BloodT1 Function to return the predicted T1 of blood, based on Hales' JCBFM paper
% assumes we want this for in vivo studies, so adds offset of 108ms
% By Patrick Hales

% If Hct unknown, use Hct=xASL_quant_AgeSex2Hct(age, gender) with 0=female, 1=male
% If Y unknown, assume Y(arterial)=0.97, Y(venous)=0.68


    if ~exist('Y','var')
        Y   = 0.97;
    end
    if ~exist('B0','var')
        B0 = 3; % 3T field strength
    end

    Hb = 5.15;  % mean corpuscular haemoglobin concentration (mmol Hb tetramer / L plasma)
    
    fe = (0.70*Hct)/((0.70*Hct)+(0.95*(1-Hct)));
    
    part1 = 1.099 - (0.057*B0) + ((0.033*Hb)*(1-Y));
    part2 = (1-fe)*(0.496-(0.023*B0));
    
    t1bl =(1/((fe*part1)+part2)) + 0.108;  % (seconds)

    t1bl    = t1bl*1000; % use this in ms in ExploreASL
    fprintf('%s\n',['Calculated blood T1 ' num2str(t1bl) 'ms from Hct ' num2str(Hct)]);
end