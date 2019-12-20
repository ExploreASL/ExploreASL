%% Create SliceReadoutTime.mat

for iS=1:x.nSubjects
    Site{iS,1}                      = x.SUBJECTS{iS};
    SliceReadoutTime{iS,1}          = x.SUBJECTS{iS};
    if      str2num(x.SUBJECTS{iS}(1))==1
            Site{iS,2}              = 'Intera';
            SliceReadoutTime{iS,2}  = 34.9;
    elseif  str2num(x.SUBJECTS{iS}(1))==2
            Site{iS,2}              = 'Ingenia';
            SliceReadoutTime{iS,2}  = 37.0625;
    else
        error('Wrong site definition');
    end
end

save( fullfile(x.D.ROOT,'Site.mat'),'Site');
save( fullfile(x.D.ROOT,'SliceReadoutTime.mat'),'SliceReadoutTime');

Subj

TEshortest  = {'151005_2' '171042_2' '171110_2' '172017_2' '172062_2' '172067_2' '172069_2' '173032_2' '173049_2' '174002_2' '174024_2' '174055_2' '231069_2' '232001_2' '232004_2' '232042_2' '234017_2'};


for ii=1:length(age)
    SliceReadoutTime{ii*2-1,1}  = age{ii,1};
    SliceReadoutTime{ii*2-0,1}  = age{ii,1};
    
    SliceReadoutTime{ii*2-1,2}  = 'ASL_1';
    SliceReadoutTime{ii*2-0,2}  = 'ASL_2';
    
    if      strcmp(age{ii,1}(end),'1')
            SliceReadoutTime{ii*2-1,3}  = 34.9;
            SliceReadoutTime{ii*2-0,3}  = 34.9;
    elseif  strcmp(age{ii,1}(end),'2')
        
            SliceReadoutTime{ii*2-1,3}  = 37.86667;
            SliceReadoutTime{ii*2-0,3}  = 37.86667;             
        
            for iT=1:length(TEshortest)
                if  strcmp(TEshortest{iT},age{ii,1})
                    SliceReadoutTime{ii*2-1,3}  = 42.6;  
                end
            end
    else    error('Unknown TimePoint');
    end
    
end

save( fullfile(x.D.ROOT,'SliceReadoutTime.mat'),'SliceReadoutTime');

%% list groups here (will be used for statistical "sets")
% Take scanner group from list that Tanja provided
x.group{1}.name                        = 'scanner';
x.group{1}.options 			        = {'Intera' 'Ingenia'};
x.group{1}.sets_1_2_sample  	        = 2; % 2-sample

x.TotalSubjects                        = xASL_adm_GetFsList(x.D.ROOT, x.subject_regexp ,true);
x.nTotalSubjects                       = length(x.TotalSubjects);

for iSubject=1:x.nTotalSubjects
    x.group{1}.code{iSubject,1}        = x.TotalSubjects{iSubject};

	if   	str2num(x.group{1}.code{iSubject,1}(1))==1
			x.group{1}.code{iSubject,2}    = 1;
	elseif	str2num(x.group{1}.code{iSubject,1}(1))==2
			x.group{1}.code{iSubject,2}    = 2;
	else 	error('Invalid x.TotalSubjects{iSubject}(1)!!!');
	end
end
