% xASL_Motion_ZigZagCheck
clear
% x.D.ROOT    = 'C:\Backup\ASL\Harmy\analysis_Harmy';
x.D.ROOT    = 'C:\Backup\ASL\SleepStudy\Analysis3';

MotionMat       = fullfile(x.D.ROOT,'MeanMotion.mat');
load(MotionMat);



% Mean motion
for iM=1:length(MeanMotion); Motion(iM,1)    = MeanMotion{iM,3}; end
Motion(:,2)     = abs(Motion(:,1)-xASL_stat_MeanNan(Motion(:,1)));
Motion(:,3)     = [1:1:length(Motion)];
[A B] = find(Motion(:,2)==min(Motion(:,2))); % find subject closest to mean motion
Subject2Load{1}    = MeanMotion{A,1};

% Min motion
[A B] = find(Motion(:,1)==min(Motion(:,1))); % find subject closest to mean motion
Subject2Load{2}    = MeanMotion{A,1};

% Max motion
Motion              = sortrows(Motion,1);
A                   = Motion(round(0.95*length(Motion)),3);
Subject2Load{3}    = MeanMotion{A,1};

SubLoad             = {'MeanMotionSubject' 'MinMotionSubject' 'MaxMotionSubject'};

for iS=1:3
    clear Txt2Load rp_ASL4D M rp_ASL4D_BeforeSpikeExclusion
    Txt2Load        = fullfile(x.D.ROOT,Subject2Load{iS},'ASL_1','rp_ASL4D.txt');
    if ~exist(Txt2Load,'file')
        Txt2Load        = fullfile(x.D.ROOT,Subject2Load{iS},'ASL_1','rp_ASL4D_BeforeSpikeExclusion.txt');
    end
    load(Txt2Load);
    if  exist('rp_ASL4D','var')
        M               = rp_ASL4D;
    else
        M               = rp_ASL4D_BeforeSpikeExclusion;
    end
    X_Ax            = [1:1:length(M)];

    TitleChoose     = {'Translation_X' 'Translation_Y' 'Translation_Z' 'Rotation_X_pitch' 'Rotation_Y_roll' 'Rotation_Z_Yaw'};

%     %% 6 graphs
%     figure(1);
%     for ii=1:6
%         subplot(2,3,ii);
%         plot(X_Ax,M(:,ii),X_Ax,M(:,ii+6));
%         title(TitleChoose{ii});
%         xlabel('raw unsubtracted EPI volumes');
%         ylabel('position, blue & red = without & with zigzag removal');
%     end
%     x0=10; y0=10; width=1280; height=720;
%     set(gcf,'units','points','position',[x0,y0,width,height]);
    
    %% 2 graphs
    figure(1);
    subplot(2,2,1); plot(X_Ax,M(:,1),X_Ax,M(:,2),X_Ax,M(:,3));
    title('Translation X (blue) Y (red) Z (yellow), no correction');
    xlabel('raw unsubtracted EPI volumes');
    ylabel('position (mm)');
    axis([0 120 -0.3 0.5]);
    
    subplot(2,2,2); plot(X_Ax,M(:,1+6),X_Ax,M(:,2+6),X_Ax,M(:,3+6));
    title('Translation X (blue) Y (red) Z (yellow), zig-zag correction');
    xlabel('raw unsubtracted EPI volumes');
    ylabel('position (mm)');    
    axis([0 120 -0.3 0.5]);
    
    subplot(2,2,3); plot(X_Ax,M(:,4)*57.2958,X_Ax,M(:,5)*57.2958,X_Ax,M(:,6)*57.2958);
    title('Rotation pitch (blue) roll (red) yaw (yellow), no correction');
    xlabel('raw unsubtracted EPI volumes');
    ylabel('position (degrees)');
    axis([0 120 -0.1 0.5]);

    subplot(2,2,4); plot(X_Ax,M(:,4+6)*57.2958,X_Ax,M(:,5+6)*57.2958,X_Ax,M(:,6+6)*57.2958);
    title('Rotation pitch (blue) roll (red) yaw (yellow), zig-zag correction');
    xlabel('raw unsubtracted EPI volumes');
    ylabel('position (degrees)');
    axis([0 120 -0.1 0.5]);
    
    
    x0=10; y0=10; width=1280; height=720;
    set(gcf,'units','points','position',[x0,y0,width,height]);    
    
    SaveDir     = 'C:\Backup\ASL\ExploreASL_paper\MotionZigZag';
    xASL_adm_CreateDir(SaveDir);
    FileName    = fullfile(SaveDir,[Subject2Load{iS} '_' SubLoad{iS} '.eps']); % jpg
    print(FileName,'-depsc'); % djpeg
    close all;
end
