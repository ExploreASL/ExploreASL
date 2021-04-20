function [list , QC_Parameters]= xQC_create_parameter_Table(configfile)

%Create Table of parameter where to store the results with just the
%parameters of interest (Visualize == 1)


list = readtable(configfile, 'Sheet', 'ParameterSheet');


QC_Parameters = {};

for iPar = 1:height(list)
    if list.Visualize(iPar)
        QC_Parameters = [QC_Parameters, list.ParameterName(iPar)];
    end
end



end