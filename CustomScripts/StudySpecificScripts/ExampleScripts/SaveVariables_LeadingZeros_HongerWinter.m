for ID=1:length(IDnr);
    SOE_id{ID,1}={['SOE_'  num2str(IDnr(ID),'%04d')]};
end

exposure    = SOE_id;
sex         = SOE_id;
age         = SOE_id;
cohort2     = SOE_id;

save( fullfile(x.D.ROOT, 'exposure.mat'), 'exposure');
save( fullfile(x.D.ROOT, 'sex.mat'), 'sex');
save( fullfile(x.D.ROOT, 'age.mat'), 'age');
save( fullfile(x.D.ROOT, 'cohort2.mat'), 'cohort2');