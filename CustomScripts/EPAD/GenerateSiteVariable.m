%% Generate Site.mat

for iSubject=1:x.nSubjects
    Site{iSubject,1} = x.SUBJECTS{iSubject};
    Site{iSubject,2} = x.SUBJECTS{iSubject}(1:3);
end

save(fullfile(x.D.ROOT,'Site.mat' ),'Site');