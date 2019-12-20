clear AcqDate AcqTime
rawROOT     = 'C:\Backup\ASL\BioCog_Repro\raw';
Dlist       = xASL_adm_GetFsList(rawROOT,'^MRU.*$',1);

for iS=1:length(Dlist)
    xASL_TrackProgress(iS,length(Dlist));
    clear ASLname tDCM
    ASLname         = fullfile(rawROOT,Dlist{iS},'M0.dcm');
    tDCM            = dicominfo(ASLname);
    AcqDate(iS,2)   = str2num(tDCM.AcquisitionDateTime(1: 8));
    AcqTime(iS,2)   = str2num(tDCM.AcquisitionDateTime(9:12));
end

Dlist       = xASL_adm_GetFsList(rawROOT,'^BioCog_MRU.*$',1)';
for iS=1:length(Dlist)
    xASL_TrackProgress(iS,length(Dlist));
    clear ASLname tDCM
    Dir2    = xASL_adm_GetFsList(fullfile(rawROOT,Dlist{iS},'mri'),'^WINTERER.*$',1);
    Dir3    = xASL_adm_GetFsList(fullfile(rawROOT,Dlist{iS},'mri',Dir2{1}),'^.*M0_tra_Augen.*$',1);
    File1   = xASL_adm_GetFileList( fullfile(rawROOT,Dlist{iS},'mri',Dir2{1},Dir3{1}), '^.*\.dcm$','FPListRec',[0 Inf]);
    ASLname = File1{1};
    tDCM            = dicominfo(ASLname);
    AcqDate(iS,1)   = str2num(tDCM.AcquisitionDate(1: 8));
    AcqTime(iS,1)   = str2num(tDCM.AcquisitionTime(1:4));
end

AcqTime     = ConvertTimeNr(AcqTime);


[Xb Nb]     = hist(AcqTime(:,1)');
[Xu Nu]     = hist(AcqTime(:,2)');

figure(1);plot(Nb,Xb,'r-',Nu,Xu,'b-')
xlabel('Time (hours)');
ylabel('n Subjects');
title('ASL acquisition time (from DICOM)');
axis([10 14.5 0 4.5]);

RangeColor  = [1:256]./256;
blue        = [RangeColor' zeros(256,1) zeros(256,1)];

figure(1);bar(Nb,Xb); colormap([1 0 0]);
xlabel('Time (hours)'); ylabel('n Subjects');
title('ASL acquisition time (from DICOM)'); axis([10 14.5 0 4]);
figure(2);bar(Nu,Xu); colormap([0 0 1]);
xlabel('Time (hours)'); ylabel('n Subjects');
title('ASL acquisition time (from DICOM)'); axis([10 14.5 0 4]);

% Date
[AcqDate DayInYear]     = ConvertDateNr(AcqDate);
AcqDiff                 = DayInYear(:,1)-DayInYear(:,2);

[Xb Nb]     = hist(AcqDate(:,1)');
[Xu Nu]     = hist(AcqDate(:,2)');



figure(1);bar(Nb,Xb); colormap([1 0 0]);
xlabel('Date (months)');
ylabel('n Subjects');
title('ASL acquisition date (from DICOM)');
axis([10.5 13 0 4]);

figure(2);bar(Nu,Xu); colormap([0 0 1]);
xlabel('Date (months)');
ylabel('n Subjects');
title('ASL acquisition date (from DICOM)');
axis([10.5 13 0 4]);

[XX NN]     = hist(AcqDiff);
figure(3);bar(NN,XX); colormap(gray);
xlabel('Time difference Utrecht-Berlin (days)');
ylabel('n Subjects');
title('ASL acquisition date (from DICOM)');


% xInt    = linspace(10, 14.5, 100);
% SplB    = spline(Nb,Xb);
% SplU    = spline(Nu,Xu);
% figure(2);plot(xInt,ppval(SplB,xInt),'b-',xInt,ppval(SplU,xInt),'r-');



% coffee break & lunch time differs between Utrecht & Berlin, but both in
% the morning, which is good, same hydration & postprandial CO2 status

[Xb Nb]     = hist(AcqDate(:,1)');
[Xu Nu]     = hist(AcqDate(:,2)');

