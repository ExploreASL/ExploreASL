function QC = xQC_missing(QC,sSubject,modality, configfile)


list = readtable(configfile, 'Sheet', 'ParameterSheet');
list(list.Visualize == 1, :);


for iPar = 1:height(list)
    
    if strcmp(list.Scantype{iPar},modality) 
    QC.(sSubject).(list.Scantype{iPar}).(list.Domain{iPar}).(list.Parameter{iPar})= NaN;
    end
end 



end 