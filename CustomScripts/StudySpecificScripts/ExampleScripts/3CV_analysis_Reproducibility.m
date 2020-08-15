%% 2) Obtain list of mean CBF & check for normality
% Check normality
for vendor=1:12
    for subject=1:size(CBF{vendor},4)
        for k=1:length(segm_mask_total)
            clear temp
            temp                                            =CBF{vendor}(:,:,:,subject);
            temp                                            =temp( isfinite(temp) & segm_mask_total{k}{1,vendor}(:,:,:,subject) );

            if      sum(isfinite(temp(:)))<10
                    swtest_mean{vendor}(subject,k)=0;
            else    [swtest_mean{vendor}(subject,k),p,W]    =swtest(temp);end
        end
    end
end
All were distributed abnormally

% Obtain median numbers
for vendor=1:12
    for subject=1:size(CBF{vendor},4)
        for j=1:8 % pre-defined ROIs
            temp                                          =CBF{vendor}(:,:,:,subject);
            median_CBF_regions{vendor}(subject,j)         =median(median(median( temp(isfinite(temp) & segm_mask_total{j}{1,vendor}(:,:,:,subject) ) )));
            clear temp
        end
        median_CBF_regions{vendor}(subject, 9)            =median_CBF_regions{vendor}(subject,1)        ./ median_CBF_regions{vendor}(subject,2); % GM-WM CBF ratio
        median_CBF_regions{vendor}(subject,10)            =median_CBF_regions{vendor}(subject,3)        ./ median_CBF_regions{vendor}(subject,4); % ACA-MCA CBF ratio
        median_CBF_regions{vendor}(subject,11)            =median_CBF_regions{vendor}(subject,4)        ./ median_CBF_regions{vendor}(subject,5); % MCA-PCA CBF ratio
    end
end

% % Check normality median list
for vendor=1:12
    if  size(median_CBF_regions{vendor},1)>1
        swtest_CBF_vendors(vendor,j)   =swtest(median_CBF_regions{vendor}(:,j));end;end
% Median numbers were normally distributed -> mean can be used

% Obtain mean of median
for vendor=1:12
    for j=1:11
        if  size(median_CBF_regions{vendor},1)>1
            mean_CBF_vendors(vendor,j)      =mean(median_CBF_regions{vendor}(:,j));
        end
    end
end


% Combine sessions
for ii=1:6
    for j=1:11
        median_CBF_2sessions{ii}(1                                     :size(median_CBF_regions{(ii*2)-1},1)  ,j)               =median_CBF_regions{(ii*2)-1}(:,j);
        median_CBF_2sessions{ii}(1+size(median_CBF_regions{(ii*2)-1},1):size(median_CBF_regions{(ii*2)-1},1)*2,j)               =median_CBF_regions{(ii*2)-0}(:,j);
    end
end

% Obtain mean CBF differences between vendors
for subject=1:22
    for ROI=1:11
        vendor_diff{1}(subject,ROI)         =median_CBF_2sessions{2}(subject,ROI) - median_CBF_2sessions{1}(subject,ROI);
    end
end
for subject=1:28
    for ROI=1:11
        vendor_diff{2}(subject,ROI)         =median_CBF_2sessions{5}(subject,ROI) - median_CBF_2sessions{4}(subject,ROI);
    end
end
for ROI=1:11 % for vendor A study 1 vs. vendor A  study 2
    for subject=1:11
        vendor_diff{3}(subject,ROI)         =median_CBF_2sessions{4}(subject,ROI) - median_CBF_2sessions{2}(subject,ROI);
    end
    for subject=12:22
        vendor_diff{3}(subject,ROI)         =median_CBF_2sessions{4}(subject+3,ROI) - median_CBF_2sessions{2}(subject,ROI);
    end
end

% Mean vendordiff
for vendors=1:3
    for ROI=1:size(vendor_diff{vendors},2)
        mean_vendordiff{vendors}(1,ROI)     =mean( vendor_diff{vendors}(:,ROI) );
    end
    mean_vendordiff{vendors}                =mean_vendordiff{vendors}';
end


% % create combined vendor A vendor A
% median_CBF_2sessions{7}             =median_CBF_2sessions{2};
% median_CBF_2sessions{7}(1+size(median_CBF_2sessions{7},1):size(median_CBF_2sessions{7},1)+size(median_CBF_2sessions{4},1),:)             =median_CBF_2sessions{4};

% Obtain mean of median for combined sessions
for ii=1:6
    for j=1:11
        mean_CBF_2sessions(ii,j)        =mean(median_CBF_2sessions{ii}(:,j));
    end
end

% Obtain inter-session paired differences
for  ii=1:6
        paired_diff{ii}             =median_CBF_regions{ii*2}       - median_CBF_regions{(ii*2)-1};
end

% % Check normality inter-session paired differences
% for  ii=1:length(paired_diff)
%     for j=1:9
%      swtest_2_diff(ii,j)            =swtest(paired_diff{ii}(:,j));end;end
% % All normally distributed

% Obtain mean of paired_diff
for ii=1:length(paired_diff)
    for j=1:11
        mean_2_diff(ii,j)            =mean(paired_diff{ii}(:,j));
    end
end

% Obtain SD of paired_diff
for ii=1:length(paired_diff)
    for j=1:11
        SD_2_diff(ii,j)              =mean( ( paired_diff{ii}(:,j)      - mean(paired_diff{ii}(:,j)) ).^2)^0.5;
    end
end

wsCV                            =100.*(SD_2_diff     ./mean_CBF_2sessions);

wsCV_ratio(1,:)                 =wsCV(3,:)./mean(wsCV(1:2,:),1);
wsCV_ratio(2,:)                 =wsCV(6,:)./mean(wsCV(4:5,:),1);

wsCV_diff(:,1)                  =wsCV(2,:)'-wsCV(1,:)';
wsCV_diff(:,2)                  =wsCV(5,:)'-wsCV(4,:)';
wsCV_diff(:,3)                  =wsCV(4,:)'-wsCV(2,:)';

% Obtain SD of mean
for ii=1:6
    for j=1:11
        SD_2_mean(ii,j)         =mean( ( median_CBF_2sessions{ii}(:,j)      - mean(median_CBF_2sessions{ii}(:,j)) ).^2).^0.5;
    end
end
%% 2) CI values
% Paragraph 2: Precision of the estimated limits of agreement of Bland,
% Altman et al. Measuring agreement in method comparison studies (1999)
% We assume normal distribution

% CI of X would be determined by X-1.96*SD and X+1.96*SD in a normal
% distribution. But with a smaller sample size, this normal distribution
% tails will be broader and longer. Replacing 1.96 by the t-statistic
% accounts for this.

% Obtain t-values, two-tailed (p=0.05/2):
% http://www.google.nl/imgres?imgurl=http://faculty.ksu.edu.sa/salghamdi/Statistical%2520Tables
% /T%2520Table.jpg&imgrefurl=http://faculty.ksu.edu.sa/salghamdi/Pages/StatisticalTables.aspx&h=1692&w=2585&sz=621&tbnid=TccUxA5-UcsHtM:&tbnh=90&tbnw=138&zoom=1&usg
% =__nCr5ml-gPVh7SGpqDHyUC1frFx4=&docid=PBSLuedjmB7J0M&sa=X&ei=yjRlUvKfO4eThQfwhIHwCw&sqi=2&ved=0CDkQ9QEwBA

% with n-1 degrees of freedom
n_masks                 =11;        % Composed of 2 1st level GM WM 3 2nd level vascular 3 3rd level WFU Pick 3 ratios
t(1:2,1:n_masks)        =2.228;     % intra-vendor comparison (n=11; n-1=10)
t(  3,1:n_masks)        =2.080;     % intra-vendor comparison (n=22; n-1=21)
t(4:5,1:n_masks)        =2.160;     % intra-vendor comparison (n=14; n-1=13)
t(  6,1:n_masks)        =2.052;     % intra-vendor comparison (n=28; n-1=27)

n(1:2,1:n_masks)        =11;
n(  3,1:n_masks)        =22;
n(4:5,1:n_masks)        =14;
n(  6,1:n_masks)        =28;

% Precision of mean
% CI_lo_mean2session              =t .* (SD_2_mean./((n).^0.5));
% CI_hi_mean2session              =t .* (SD_2_mean./((n).^0.5));

% SD of mean
CI_lo_mean2session              = SD_2_mean;
CI_hi_mean2session              = SD_2_mean;

% Precision of mean diff
CI_lo_mean_diff2session         =t .* (SD_2_diff./((n ).^0.5));
CI_hi_mean_diff2session         =t .* (SD_2_diff./((n ).^0.5));

CI_lo_mean_diff2session         =round(10.*CI_lo_mean_diff2session)./10;
CI_hi_mean_diff2session         =round(10.*CI_hi_mean_diff2session)./10;

% SD of mean diff
CI_lo_mean_diff2session         = SD_2_diff;
CI_hi_mean_diff2session         = SD_2_diff;

CI_lo_mean_diff2session         = round(10.*CI_lo_mean_diff2session)./10;
CI_hi_mean_diff2session         = round(10.*CI_hi_mean_diff2session)./10;

% Precision of SD diff
CI_lo_SD_diff2session           =t .* (SD_2_diff./ ( (2.*(n-1) ).^0.5));
CI_hi_SD_diff2session           =t .* (SD_2_diff./ ( (2.*(n-1) ).^0.5));

CI_lo_SD_diff2session           =round(10.*CI_lo_SD_diff2session)./10;
CI_hi_SD_diff2session           =round(10.*CI_hi_SD_diff2session)./10;

% Precision of wsCV
% Same applies as above. Mean_CBF and RMS_SD_diff are independent, therefore
% variance of wsCV is sum of variances mean_CBF and RMS_SD_diff.

CI_lo_wsCV                      =t .* ( (SD_2_mean./((n).^0.5)).^2 + (SD_2_diff./ ( (2.*(n-1) ).^0.5)).^2 ).^0.5;
CI_hi_wsCV                      =t .* ( (SD_2_mean./((n).^0.5)).^2 + (SD_2_diff./ ( (2.*(n-1) ).^0.5)).^2 ).^0.5;

CI_lo_wsCV                      =round(10.*CI_lo_wsCV)./10;
CI_hi_wsCV                      =round(10.*CI_hi_wsCV)./10;
%% 2) BA values
% % Mean
% for ii=1:6
%         mean_median_CBF_2sessions{ii}                                                   =( median_CBF_regions{(ii*2)-1} + median_CBF_regions{(ii*2)-0} )./2;
%         diff_median_CBF_2sessions{ii}                                                   =median_CBF_regions{(ii*2)-0} - median_CBF_regions{(ii*2)-1};
%         for j=1:9
%             combined_BA{j}(1:length(mean_median_CBF_2sessions{ii}(:,j)),(ii*2)-1)       =mean_median_CBF_2sessions{ii}(:,j);
%             combined_BA{j}(1:length(diff_median_CBF_2sessions{ii}(:,j)),(ii*2)-0)       =diff_median_CBF_2sessions{ii}(:,j);
%             combined_BA2{j}(ii,1)                                                       =mean_2_diff(ii,j);
%             combined_BA2{j}(ii,2)                                                       =LOA_lo(ii,j);
%             combined_BA2{j}(ii,3)                                                       =LOA_hi(ii,j);
%     end
% end
%% 2) Table values
% for ii=1:11
%     for j=1:6
%         table((ii*2)-1,j*3-2)     =mean_CBF_2sessions(j,ii);
%         table((ii*2)-0,j*3-2)     =wsCV(j,ii);
%
%         table((ii*2)-1,j*3-1)     =CI_lo_mean2session(j,ii);
%         table((ii*2)-0,j*3-1)     =CI_lo_wsCV(j,ii);
%
%         table((ii*2)-1,j*3-0)     =CI_hi_mean2session(j,ii);
%         table((ii*2)-0,j*3-0)     =CI_hi_wsCV(j,ii);
%     end
% end
%
% table( 1:16,:)                       =round(table( 1:16,:).*10)./10;
% table(17:22,:)                       =round(table(17:22,:).*100)./100;

diffVendor                      = roundPoint( 100.*abs( mean_CBF_2sessions([1 4],:) - mean_CBF_2sessions([2 5],:) ) ./ mean_CBF_2sessions([3 6],:),0.1);
diffVendor(3,:)                 = roundPoint( 100.*abs( mean_CBF_2sessions(2,:)     - mean_CBF_2sessions(4,:)     ) ./ mean(mean_CBF_2sessions([2 4],:),1),0.1);

ratioVendor(1,:)                = roundPoint( wsCV(3,:) ./ mean(wsCV([1 2],:),1),0.1);
ratioVendor(2,:)                = roundPoint( wsCV(6,:) ./ mean(wsCV([4 5],:),1),0.1);
ratioVendor(3,:)                = roundPoint( wsCV(2,:) ./ wsCV(4,:),0.1);


for ii=1:11
    for j=1:6
        tableMEAN{ii,j*3-2}     = [ num2str(roundPoint(mean_CBF_2sessions(j,ii),0.1)) ' (' num2str(roundPoint(CI_lo_mean2session(j,ii),0.1)) ')'];

        tableWSCV{ii,j*3-2}     = [ num2str(roundPoint(wsCV(j,ii),0.1)) ' (' num2str(roundPoint(CI_lo_wsCV(j,ii),0.1)) ')'];
    end
end

for ii=1:11
    for j=1:3
        tableDIFF{ii,j}         = num2str(diffVendor(j,ii));

%         tableWSCV{ii,j*3-2}     =[ num2str(roundPoint(wsCV(j,ii),0.1)) ' (' num2str(roundPoint(CI_lo_wsCV(j,ii),0.1)) ')'];
    end
end

for ii=1:11
    for j=1:3
        tableRATIO{ii,j}         = num2str(ratioVendor(j,ii));
    end
end

%% 2) Statistics
% Inter-session difference
for ii=1:6
    for j=1:11
        [mean_diff_statistics_GM(ii,(j*2)-1) mean_diff_statistics_GM(ii,j*2)]       =ttestExploreASL (median_CBF_regions{(ii*2)-1}(:,j),median_CBF_regions{(ii*2)-0}(:,j),0.05,'both');
    end
end

for j=1:11 % Vendor A study 1 vs. vendor A study 2 (Philips), cave n=1:11 & 15:25
    temp_study1( 1:11)                                                           =median_CBF_regions{3}(   :,j);
    temp_study1(12:22)                                                           =median_CBF_regions{4}(   :,j);
    temp_study2( 1:11)                                                           =median_CBF_regions{7}(1:11,j);
    temp_study2(12:22)                                                           =median_CBF_regions{8}(1:11,j);

    [mean_diff_statistics_GM(7,(j*2)-1) mean_diff_statistics_GM(7,j*2)]          =ttestExploreASL ( temp_study1, temp_study2, 0.05, 'both');
	clear temp_study1 temp_study2
end



%  diff variation GE vs PH
for j=1:8
    clear X
    X( 1:11,1)                      =paired_diff{1}(1:11,j);
    X(12:22,1)                      =paired_diff{2}(1:11,j);
    X( 1:11,2)                      =1;
    X(12:22,2)                      =2;
    [GE_vs_PH(1,j) GE_vs_PH(2,j)]   =Levenetest_HJM(X,0.05);end

%  diff variation PH vs SI
for j=1:8
    clear X
    X( 1:11,1)                      =paired_diff{4}(1:11,j);
    X(12:22,1)                      =paired_diff{5}(1:11,j);
    X( 1:11,2)                      =1;
    X(12:22,2)                      =2;
    [PH_vs_SI(1,j) PH_vs_SI(2,j)]   =Levenetest_HJM(X,0.05);end

% Compare inter-vendor inter-session differences with mean intra-vendor differences

for j=1:8
    clear X
    X( 1:11,1)       =paired_diff{1}(:,j);
    X(12:22,1)       =paired_diff{2}(:,j);
    X(23:44,1)       =paired_diff{3}(:,j);
    X( 1:22,2)       =1;
    X(23:44,2)       =2;
    if  mean( X( 1:11,1) ) * mean( X( 12:22,1) ) ~= (mean( X( 1:11,1) )^2)^0.5 * (mean( X(12:22,1) )^2)^0.5
        X(12:22,1)          =X(12:22,1) .*-1;end % To correct if mean difference is other side
    [GE_PH_vs_inter(1,j) GE_PH_vs_inter(2,j)]   =Levenetest_HJM(X,0.1);end

for j=1:8
    clear X
    X( 1:14,1)       =paired_diff{4}(:,j);
    X(15:28,1)       =paired_diff{5}(:,j);
    X(29:56,1)       =paired_diff{6}(:,j);
    X( 1:28,2)       =1;
    X(29:56,2)       =2;
    if  mean( X( 1:11,1) ) * mean( X( 12:22,1) ) ~= (mean( X( 1:11,1) )^2)^0.5 * (mean( X(12:22,1) )^2)^0.5
        X(12:22,1)          =X(12:22,1) .*-1;end % To correct if mean difference is other side
    [PH_SI_vs_inter(1,j) PH_SI_vs_inter(2,j)]   =Levenetest_HJM(X,0.1);end


%  diff variation PH study 1 vs PH study 2
for j=1:8
    clear X
    X( 1:11,1)                      =paired_diff{2}(1:11,j);
    X(12:22,1)                      =paired_diff{4}(1:11,j);
    X( 1:11,2)                      =1;
    X(12:22,2)                      =2;
    [PH_vs_PH(1,j) PH_vs_PH(2,j)]   =Levenetest_HJM(X,0.05);end


%% 3) Visualization mean (GE & PH switched)
% for ii=1:length(CBF)
%     for j=1:size(CBF{ii},1)
%         for k=1:size(CBF{ii},2)
%             for l=1:size(CBF{ii},3)
%                 temp                    =CBF{ii}(j,k,l,:);
%                 mean_CBF{ii}(j,k,l)     =mean(temp(isfinite(temp)));end;end;end;end

% save('C:\Backup\ASL\3CV\Analysis_temp\CBF_mean.mat','mean_CBF');
load('E:\Backup\ASL_E\3CV\Analysis_temp\CBF_mean.mat');

for ii=1:(length(mean_CBF)/2)
    mean_CBF_vendor{ii}             =(mean_CBF{(ii*2)-1}+mean_CBF{(ii*2)-0})./2;end

mean_CBF_vendor{7}                  =(mean_CBF_vendor{2}+mean_CBF_vendor{4})./2;

clear slice_mean_CBF slice_view_mean_CBF overview_mean
for ii=1:length(mean_CBF_vendor)
    slice_mean_CBF{1}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(mean_CBF_vendor{ii}(:,:,62).*skull_mask_visual(:,:,62),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_CBF{2}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(mean_CBF_vendor{ii}(:,:,76).*skull_mask_visual(:,:,76),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_CBF{3}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(mean_CBF_vendor{ii}(:,63,:).*skull_mask_visual(:,63,:)),90),sup_crop,size(CBF{ii},2)-inf_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_CBF{4}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(mean_CBF_vendor{ii}(:,95,:).*skull_mask_visual(:,95,:)),90),sup_crop,size(CBF{ii},2)-inf_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_CBF{5}(:,:,ii)         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( mean_CBF_vendor{ii}(68,:,:) .*skull_mask_visual(68,:,:) ),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);
    slice_mean_CBF{6}(:,:,ii)         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( mean_CBF_vendor{ii}(76,:,:) .*skull_mask_visual(76,:,:) ),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);
end

% Switch GE & PH (vendor B & A)
for ii=1:6
    temp                          =slice_mean_CBF{ii}(:,:,1);
    slice_mean_CBF{ii}(:,:,1)     =slice_mean_CBF{ii}(:,:,2);
    slice_mean_CBF{ii}(:,:,2)     =temp;
end

for ii=[1 2 4 5]
    overview_mean(:,:,ii)             =[slice_mean_CBF{1}(:,:,ii) slice_mean_CBF{2}(:,:,ii) slice_mean_CBF{3}(:,:,ii) slice_mean_CBF{4}(:,:,ii) slice_mean_CBF{5}(:,:,ii) slice_mean_CBF{6}(:,:,ii)];
end

% figure(1);imshow([overview_mean(:,:,1);overview_mean(:,:,2)],[0 80],'Colormap',jet256,'InitialMagnification',100);
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_mean1.tiff']);
% figure(2);imshow([overview_mean(:,:,4);overview_mean(:,:,5)],[0 80],'Colormap',jet256,'InitialMagnification',100);
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_mean2.tiff']);
%% 3) Visualization mean_diff (GE & PH switched)
for ii=1:(length(CBF)/2)
    diff{ii}            =CBF{(ii*2)-1}-CBF{(ii*2)-0};end

for ii=1:length(diff)
    mean_diff{ii}       =mean(diff{ii},4);end

for ii=1:length(mean_diff)
    mean_diff{ii}(mean_diff{ii}<-15)    =-15;
    mean_diff{ii}(mean_diff{ii}>15)     =15;
    mean_diff{ii}                       = rescale( mean_diff{ii},0,30,0);end

for ii=1:length(mean_diff)
    slice_mean_diff{1}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(mean_diff{ii}(:,:,62).*skull_mask_visual(:,:,62),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_diff{2}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(mean_diff{ii}(:,:,76).*skull_mask_visual(:,:,76),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_diff{3}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(mean_diff{ii}(:,63,:).*skull_mask_visual(:,63,:)),90),sup_crop,size(CBF{ii},2)-inf_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_diff{4}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(mean_diff{ii}(:,95,:).*skull_mask_visual(:,95,:)),90),sup_crop,size(CBF{ii},2)-inf_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_mean_diff{5}(:,:,ii)         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( mean_diff{ii}(68,:,:) .*skull_mask_visual(68,:,:) ),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);
    slice_mean_diff{6}(:,:,ii)         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( mean_diff{ii}(76,:,:) .*skull_mask_visual(76,:,:) ),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);end

% Switch GE & PH (vendor B & A)
for ii=1:length(mean_diff)
    temp                          =slice_mean_diff{ii}(:,:,1);
    slice_mean_diff{ii}(:,:,1)     =slice_mean_diff{ii}(:,:,2);
    slice_mean_diff{ii}(:,:,2)     =temp;end
%
for ii=1:length(slice_mean_diff)
    slice_view_mean_diff{ii}              =singlesequencesort(slice_mean_diff{ii},1);end

% for ii=1:length(slice_view_mean_diff)
%     figure(1);imshow(slice_view_mean_diff{ii},[0 30],'Colormap',jet256);
%     print(gcf,'-dtiff','-r600',['C:\Backup\ASL\3CV\Results\Mean_CBF\Mean_diff_' num2str(ii) '.tiff']);end
%% 3) Visualization mean_diff parametric (including vendor A vs. A)
% clear CBF_slices vendor subject
% for vendor=[1:4,7:10]
%     for subject=1:size(CBF{vendor},4)
%         CBF_slices{vendor}{1}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(CBF{vendor}(:,:,62,subject).*GM_mask_visual(:,:,62),90),ant_crop,size(CBF{vendor},2)-pos_crop,L_crop,size(CBF{vendor},1)-R_crop);
%         CBF_slices{vendor}{2}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(CBF{vendor}(:,:,76,subject).*GM_mask_visual(:,:,76),90),ant_crop,size(CBF{vendor},2)-pos_crop,L_crop,size(CBF{vendor},1)-R_crop);
%         CBF_slices{vendor}{3}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(CBF{vendor}(:,63,:,subject).*GM_mask_visual(:,63,:)),90),sup_crop,size(CBF{vendor},2)-inf_crop,L_crop,size(CBF{vendor},1)-R_crop);
%         CBF_slices{vendor}{4}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(CBF{vendor}(:,95,:,subject).*GM_mask_visual(:,95,:)),90),sup_crop,size(CBF{vendor},2)-inf_crop,L_crop,size(CBF{vendor},1)-R_crop);
%         CBF_slices{vendor}{5}(:,:,subject)      =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( CBF{vendor}(68,:,:,subject) .*GM_mask_visual(68,:,:) ),sup_crop,size(CBF{vendor},2)-inf_crop,ant_crop,size(CBF{vendor},2)-pos_crop);
%         CBF_slices{vendor}{6}(:,:,subject)      =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( CBF{vendor}(76,:,:,subject) .*GM_mask_visual(76,:,:) ),sup_crop,size(CBF{vendor},2)-inf_crop,ant_crop,size(CBF{vendor},2)-pos_crop);end;end
%
% clear orislice diff_slices
% for orislice=1:6
%     diff_slices{1}{orislice}(:,:,                                1:  size(CBF_slices{1}{orislice},3))        =CBF_slices{1}{orislice} - CBF_slices{ 3}{orislice}; % GE1 vs PH1
%     diff_slices{1}{orislice}(:,:,size(CBF_slices{1}{orislice},3)+1:2*size(CBF_slices{1}{orislice},3))        =CBF_slices{2}{orislice} - CBF_slices{ 4}{orislice}; % GE2 vs PH2
%     diff_slices{2}{orislice}(:,:,                                1:  size(CBF_slices{7}{orislice},3))        =CBF_slices{7}{orislice} - CBF_slices{ 9}{orislice}; % PH1 vs SI1
%     diff_slices{2}{orislice}(:,:,size(CBF_slices{7}{orislice},3)+1:2*size(CBF_slices{7}{orislice},3))        =CBF_slices{8}{orislice} - CBF_slices{10}{orislice}; % PH2 vs SI2
%     diff_slices{3}{orislice}(:,:,                                1:  size(CBF_slices{1}{orislice},3))        =CBF_slices{3}{orislice}(:,:,1:11) - CBF_slices{ 7}{orislice}(:,:,1:11); % PH1-1 vs PH2-1
%     diff_slices{3}{orislice}(:,:,size(CBF_slices{1}{orislice},3)+1:2*size(CBF_slices{1}{orislice},3))        =CBF_slices{4}{orislice}(:,:,1:11) - CBF_slices{ 8}{orislice}(:,:,1:11); % PH1-2 vs PH2-2
% end
%
% % Paired-sample Student's t-test
% clear vendor orislice H P
% for vendor  =1:3
%     for orislice=1:6
%         [H{vendor}{orislice} P{vendor}{orislice}]=ttestExploreASL(diff_slices{vendor}{orislice},zeros(size(diff_slices{vendor}{orislice},1),size(diff_slices{vendor}{orislice},2),size(diff_slices{vendor}{orislice},3)),0.05,'both',3);
% end;end
%
% % Create post&neg color masks
% clear pos neg pos_mask neg_mask pne_mask parametric
% for orislice=1:6
%     for vendor=1:3
%         H{vendor}{orislice}(isnan(H{vendor}{orislice}))=0;
%
%         pos{vendor}{orislice}               =logical(mean(diff_slices{vendor}{orislice},3)>0);
%         neg{vendor}{orislice}               =logical(mean(diff_slices{vendor}{orislice},3)<0);
%
%         for ii=1:3
%             pos_mask{vendor}{orislice}(:,:,ii)          =double(pos{vendor}{orislice}.*H{vendor}{orislice});
%             neg_mask{vendor}{orislice}(:,:,ii)          =double(neg{vendor}{orislice}.*H{vendor}{orislice});end
%
%         if      vendor==1 | vendor==3    % Colors red & blue switched for switching vendor A & B, not vendor A & C
%                 pos_mask{vendor}{orislice}(:,:,1)   =0;
%                 pos_mask{vendor}{orislice}(:,:,2)   =0.5.*pos_mask{vendor}{orislice}(:,:,3);
%                 neg_mask{vendor}{orislice}(:,:,2)   =0;
%                 neg_mask{vendor}{orislice}(:,:,3)   =0;
%         elseif  vendor==2
%                 neg_mask{vendor}{orislice}(:,:,1)   =0;
%                 neg_mask{vendor}{orislice}(:,:,2)   =0.5.*neg_mask{vendor}{orislice}(:,:,3);
%                 pos_mask{vendor}{orislice}(:,:,2)   =0;
%                 pos_mask{vendor}{orislice}(:,:,3)   =0;
%         end
%         pne_mask{vendor}{orislice}          =pos_mask{vendor}{orislice}+ neg_mask{vendor}{orislice};
%         parametric{vendor}{orislice}        =pne_mask{vendor}{orislice}+ (( GM_slices_clr{orislice} - pne_mask{vendor}{orislice} ).*0.5);
% end;end
%
% for vendor=[1 2 3]
%     overview_diff_parametric{vendor}        =[parametric{vendor}{1} parametric{vendor}{2} parametric{vendor}{3} parametric{vendor}{4} parametric{vendor}{5} parametric{vendor}{6}];end
%
% figure(3);imshow(overview_diff_parametric{1})
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_diff_parametric1.tiff']);
% figure(4);imshow(overview_diff_parametric{2})
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_diff_parametric2.tiff']);
% figure(5);imshow(overview_diff_parametric{3})
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_diff_parametric3.tiff']);
%% 3) CBF histograms
clear N_ROI X_ROI bin_nr min_nr max_nr bin_size myfilter

bin_nr      =117;
min_nr      =-40;
max_nr      =140;
bin_size    =(max_nr-min_nr)/bin_nr;
%myfilter    =fspecial('gaussian',[bin_nr,1],0.02*bin_nr);
% This is 117/ 180 = 0.65 bins per mL/100g/min, histograms shown have 120 range which would have been 78 bins

% Mean CBF histograms

for vendor=[5,6,11,12]
        for ii=1:size(CBF{vendor},4)
            clear temp
            temp                                                    =CBF{vendor}(:,:,:,ii);

            for k=1:8
                [N_ROI{vendor}{k}(:,ii)  X_ROI{vendor}{k}(:,ii)]    =hist(temp(logical( segm_mask_total{1,k}{1,vendor}(:,:,:,ii) ) & isfinite(temp)),[min_nr:(max_nr-min_nr)/bin_nr:max_nr]);
                N_ROI{vendor}{k}(:,ii)                              =N_ROI{vendor}{k}(:,ii)./sum(N_ROI{vendor}{k}(:,ii))./bin_size;end
        end

        for k=1:8
%            N_ROI{vendor}{k}        =imfilter(mean(N_ROI{vendor}{k},2), myfilter, 'replicate');
            N_ROI{vendor}{k}        = xASL_im_ndnanfilter(mean(N_ROI{vendor}{k},2),'gauss',[1*0.02*bin_nr*2.335 0 0],0);
		end
end

% % % % CBF_location_1  = xASL_im_rotate(mean(CBF{1},4),90);
% % % % CBF_location_2  = xASL_im_rotate(mean(CBF{2},4),90);
% % % %
% % % % mask_1          = xASL_im_rotate(mean(segm_mask_total{1}{1},4)>0.5,90);
% % % % mask_2          = xASL_im_rotate(mean(segm_mask_total{1}{2},4)>0.5,90);
% % % %
% % % % voxel_1         = CBF_location_1>0 & CBF_location_1<20 & mask_1;
% % % % voxel_2         = CBF_location_2>0 & CBF_location_2<20 & mask_2;
% % % % sum(voxel_1(:))
% % % % sum(voxel_2(:))
% % % %
% % % % dip_image([voxel_1 voxel_2])


% for k=[1 3:8] % GM
%     figure((k*2)-1);plot(X_ROI{ 5}{k}(:,1),100.*N_ROI{ 5}{k},'g',X_ROI{ 6}{k}(:,1),100.*N_ROI{ 6}{k},'r');
%     axis([-10 110 0 2.5]);
%     print(gcf,'-depsc',['C:\Backup\ASL\3CV\Results\Histograms\CBF_' num2str((k*2)-1) '.eps']);
%
%     figure((k*2)-0);plot(X_ROI{11}{k}(:,1),100.*N_ROI{11}{k},'r',X_ROI{12}{k}(:,1),100.*N_ROI{12}{k},'b');
%     axis([-30 110 0 2.5]);
%     print(gcf,'-depsc',['C:\Backup\ASL\3CV\Results\Histograms\CBF_' num2str((k*2)-0) '.eps']);
% end
%
% for k=2 % WM
%     figure((k*2)-1);plot(X_ROI{ 5}{k}(:,1),100.*N_ROI{ 5}{k},'g',X_ROI{ 6}{k}(:,1),100.*N_ROI{ 6}{k},'r');
%     axis([-30 110 0 5.6]);
%     print(gcf,'-depsc',['C:\Backup\ASL\3CV\Results\Histograms\CBF_' num2str((k*2)-1) '.eps']);
%
%     figure((k*2)-0);plot(X_ROI{11}{k}(:,1),100.*N_ROI{11}{k},'r',X_ROI{12}{k}(:,1),100.*N_ROI{12}{k},'b');
%     axis([-30 110 0 5.6]);
%     print(gcf,'-depsc',['C:\Backup\ASL\3CV\Results\Histograms\CBF_' num2str((k*2)-0) '.eps']);
% end


%% 4) Visualization SD_diff              (including vendor A vs. A)
for ii=1:(length(CBF)/2)
    diff{ii}            =CBF{(ii*2)-1}-CBF{(ii*2)-0};end

diff{7}(:,:,:,1:11)     =CBF{3}(:,:,:,1:11) - CBF{7}(:,:,:,1:11);
diff{8}(:,:,:,1:11)     =CBF{4}(:,:,:,1:11) - CBF{8}(:,:,:,1:11);

for j=1:size(mean_diff,2)
    mean_diff_n(j,1)     =mean(mean(mean( mean_diff{j}(isfinite(mean_diff{j})) )));
end

diff{21}(:,:,:, 1:11)   =diff{1};
diff{21}(:,:,:,12:22)   =diff{2};
diff{22}(:,:,:, 1:14)   =diff{4};
diff{22}(:,:,:,15:28)   =diff{5};

% if  ( mean_diff_n(1)<0 &  mean_diff_n(2)>0 ) | ( mean_diff_n(2)<0 &  mean_diff_n(1)>0 )
%     diff{21}(:,:,:, 1:11)   =-diff{21}(:,:,:, 1:11);
% end
% if  ( mean_diff_n(4)<0 &  mean_diff_n(5)>0 ) | ( mean_diff_n(5)<0 &  mean_diff_n(4)>0 )
%     diff{22}(:,:,:, 1:14)   =-diff{22}(:,:,:, 1:14);
% end

for ii=1:length(diff)
    mean_diff{ii}       =mean(diff{ii},4);end

for ii=1:length(diff)
    for j=1:size(diff{ii},4)
        mean_diff_piled{ii}(:,:,:,j)       =mean_diff{ii};end;end

for ii=1:length(diff)
    RMS_SD_diff{ii}                     =(mean((diff{ii} - mean_diff_piled{ii}).^2,4)).^0.5;end

clear mean_diff_piled

for ii=1:length(mean_diff)
    if ~isempty(mean_diff{ii})
        slice_SD_diff{1}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(RMS_SD_diff{ii}(:,:,62).*skull_mask_visual(:,:,62),90),ant_crop,size(CBF{1},2)-pos_crop,L_crop,size(CBF{1},1)-R_crop);
        slice_SD_diff{2}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(RMS_SD_diff{ii}(:,:,76).*skull_mask_visual(:,:,76),90),ant_crop,size(CBF{1},2)-pos_crop,L_crop,size(CBF{1},1)-R_crop);
        slice_SD_diff{3}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(RMS_SD_diff{ii}(:,63,:).*skull_mask_visual(:,63,:)),90),ant_crop,size(CBF{1},2)-pos_crop,L_crop,size(CBF{1},1)-R_crop);
        slice_SD_diff{4}(:,:,ii)         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(RMS_SD_diff{ii}(:,95,:).*skull_mask_visual(:,95,:)),90),ant_crop,size(CBF{1},2)-pos_crop,L_crop,size(CBF{1},1)-R_crop);
        slice_SD_diff{5}(:,:,ii)         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( RMS_SD_diff{ii}(68,:,:) .*skull_mask_visual(68,:,:)),sup_crop,size(CBF{1},2)-inf_crop,ant_crop,size(CBF{1},2)-pos_crop);
        slice_SD_diff{6}(:,:,ii)         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( RMS_SD_diff{ii}(76,:,:) .*skull_mask_visual(76,:,:)),sup_crop,size(CBF{1},2)-inf_crop,ant_crop,size(CBF{1},2)-pos_crop);
    end
end

for ii=1:length(slice_mean_diff)
    slice_view_SD_diff{ii}           =singlesequencesort(slice_SD_diff{ii},1);end


% for ii=1:length(slice_view_mean_diff)
%     figure(1);imshow(slice_view_SD_diff{ii},[0 25],'Colormap',jet256);
%     print(gcf,'-dtiff','-r600',['C:\Backup\ASL\3CV\Results\Mean_CBF\RMS_SD_' num2str(ii) '.tiff']);end
%% 4) Visualization wsCV (GE & PH switched)
for ii=1:6
    wsCV_im{ii}                     =100.*(RMS_SD_diff{ii} ./ mean_CBF_vendor{ii});end % RMS_SD_diff & mean_CBF_vendor weren't switched

for ii=1:6
    slice_wsCV{1}(:,:,ii)           =xASL_vis_CropParmsApply(xASL_im_rotate(wsCV_im{ii}(:,:,62).*GM_mask_visual(:,:,62),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_wsCV{2}(:,:,ii)           =xASL_vis_CropParmsApply(xASL_im_rotate(wsCV_im{ii}(:,:,76).*GM_mask_visual(:,:,76),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_wsCV{3}(:,:,ii)           =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(wsCV_im{ii}(:,63,:).*GM_mask_visual(:,63,:)),90),sup_crop,size(CBF{ii},2)-inf_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_wsCV{4}(:,:,ii)           =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(wsCV_im{ii}(:,95,:).*GM_mask_visual(:,95,:)),90),sup_crop,size(CBF{ii},2)-inf_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_wsCV{5}(:,:,ii)           =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( wsCV_im{ii}(68,:,:) .*GM_mask_visual(68,:,:)),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);
    slice_wsCV{6}(:,:,ii)           =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( wsCV_im{ii}(76,:,:) .*GM_mask_visual(76,:,:)),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);end

% Switch GE & PH (vendor B & A)
for ii=1:6
    temp                            =slice_wsCV{ii}(:,:,1);
    slice_wsCV{ii}(:,:,1)           =slice_wsCV{ii}(:,:,2);
    slice_wsCV{ii}(:,:,2)           =temp;end

for ii=1:length(slice_wsCV)
    slice_view_wsCV{ii}             =singlesequencesort(slice_wsCV{ii}(:,:,1:6),1);end

view_wsCV_total                     =[slice_view_wsCV{1} slice_view_wsCV{2} slice_view_wsCV{3} slice_view_wsCV{4} slice_view_wsCV{5}  slice_view_wsCV{6}];

% figure(1);imshow(view_wsCV_total,[0 50],'Colormap',jet256);

% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\wsCV_overview.tiff']);

% figure(2);imshow(view_wsCV_vendorA_A,[0 50],'Colormap',jet256);
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\wsCV_vendorA.tiff']);
%
% for ii=1:length(slice_view_mean_diff)
%     figure(ii);imshow(slice_view_wsCV{ii},[0 50],'Colormap',jet256);
%     print(gcf,'-dtiff','-r600',['C:\Backup\ASL\3CV\Results\Mean_CBF\wsCV_' num2str(ii) '.tiff']);end
%% 4) Project wsCV over GM masks
clear r_slice g_slice b_slice slice_view_wsCV_clr_prep slice_wsCV_color

slice_view_wsCV_clr_prep                                                =slice_view_wsCV;

% 1) Rescale maps
min_intensity                                                           =0;
max_intensity                                                           =50;

for ii=1:6
    slice_view_wsCV_clr_prep{ii}(slice_view_wsCV_clr_prep{ii}<min_intensity)      =min_intensity;
    slice_view_wsCV_clr_prep{ii}(slice_view_wsCV_clr_prep{ii}>max_intensity)      =max_intensity;

    slice_wsCV_color{ii}                                                =round(double(rescale( slice_view_wsCV_clr_prep{ii},0,255,1)));
    slice_wsCV_color{ii}                                                =max(1,min(slice_wsCV_color{ii},size(jet256,1)));

    r_slice{ii}     =zeros(size(slice_wsCV_color{ii})); r_slice{ii}(:)  =jet256(slice_wsCV_color{ii},1);
    g_slice{ii}     =zeros(size(slice_wsCV_color{ii})); g_slice{ii}(:)  =jet256(slice_wsCV_color{ii},2);
    b_slice{ii}     =zeros(size(slice_wsCV_color{ii})); b_slice{ii}(:)  =jet256(slice_wsCV_color{ii},3);

    b_slice{ii}( (slice_view_GM{ii}==1) & r_slice{ii}==0 & g_slice{ii}==0 & b_slice{ii}==0)    =0.5313;

    rout_slice{ii}        = zeros([size(r_slice{ii}),3]);
    rout_slice{ii}(:,:,1) = r_slice{ii}.*slice_view_skull_clr{ii}(:,:,1);
    rout_slice{ii}(:,:,2) = g_slice{ii}.*slice_view_skull_clr{ii}(:,:,2);
    rout_slice{ii}(:,:,3) = b_slice{ii}.*slice_view_skull_clr{ii}(:,:,3);

%     slice_wsCV_color{ii}                         =(0.5.*rout_slice{ii})+(0.5.*rout_slice{ii}.*slice_view_GM_clr{ii});
    slice_wsCV_color{ii}                         =rout_slice{ii};

    % Put GM map on background
    for i=1:size(slice_wsCV_color{ii},1)
        for j=1:size(slice_wsCV_color{ii},2)
            if  slice_wsCV_color{ii}(i,j,1)==0 & slice_wsCV_color{ii}(i,j,2)==0 & slice_wsCV_color{ii}(i,j,3)==0
                slice_wsCV_color{ii}(i,j,1)=0.8.*slice_view_GM_background_clr{ii}(i,j,1);
                slice_wsCV_color{ii}(i,j,2)=0.8.*slice_view_GM_background_clr{ii}(i,j,2);
                slice_wsCV_color{ii}(i,j,3)=0.8.*slice_view_GM_background_clr{ii}(i,j,3);end;end;end
end

view_wsCV_total     =[slice_wsCV_color{1} slice_wsCV_color{2} slice_wsCV_color{3} slice_wsCV_color{4} slice_wsCV_color{5}  slice_wsCV_color{6}];

% figure(1);imshow([view_wsCV_total]);
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\wsCV_overview.tiff']);
%% 4) wsCV histograms
clear N_ROI X_ROI bin_nr min_nr max_nr bin_size myfilter

bin_nr      =160;
min_nr      =0;
max_nr      =200;
bin_size    =(max_nr-min_nr)/bin_nr;
myfilter    =fspecial('gaussian',[bin_nr,1],1*0.02*bin_nr);
% 117 bins / 200 range = 0.59 bins per % wsCV, for histogram range of 100 % wsCV this whould have been 59 bins

% Mean CBF histograms

for vendor=1:6
    clear temp
    temp                        =wsCV_im{vendor};

    for k=1:8
        [N_ROI{vendor}{k}  X_ROI{vendor}{k}]    =hist(temp( mean_segm_mask{k}>0.5 & isfinite(temp)),[min_nr:(max_nr-min_nr)/bin_nr:max_nr]);
        N_ROI{vendor}{k}                        =N_ROI{vendor}{k}./sum(N_ROI{vendor}{k})./bin_size;end
    for k=1:8
        %N_ROI{vendor}{k}                        =imfilter(N_ROI{vendor}{k}, myfilter', 'replicate');
	    N_ROI{vendor}{k}                        = xASL_im_ndnanfilter(N_ROI{vendor}{k},'gauss',[0 1*0.02*bin_nr*2.335 0],0);
	end
end

% for k=1:8
%     figure((k*2)-1);plot(X_ROI{ 1}{k},100.*N_ROI{ 1}{k},'g',X_ROI{ 2}{k},100.*N_ROI{ 2}{k},'r',X_ROI{ 3}{k},100.*N_ROI{ 3}{k},'k');
%     axis([0 100 0 6]);
%     print(gcf,'-depsc',['C:\Backup\ASL\3CV\Results\Histograms\wsCV_' num2str((k*2)-1) '.eps']);
%
%     figure((k*2)-0);plot(X_ROI{ 4}{k},100.*N_ROI{ 4}{k},'r',X_ROI{ 5}{k},100.*N_ROI{ 5}{k},'b',X_ROI{ 6}{k},100.*N_ROI{ 6}{k},'k');
%     axis([0 100 0 6]);
%     print(gcf,'-depsc',['C:\Backup\ASL\3CV\Results\Histograms\wsCV_' num2str((k*2)-0) '.eps']);
% end
% close all

%% 4) Visualization ratio intra-vendor vs inter-vendor
clear ICC slice_ICC view_ICC_total

ICC{1}          = 100.* ( RMS_SD_diff{3} ./ (RMS_SD_diff{21}) ) ;
ICC{2}          = 100.* ( RMS_SD_diff{6} ./ (RMS_SD_diff{22}) ) ;

ICC{3}          = 100.* ( RMS_SD_diff{7} ./ RMS_SD_diff{8} ) ;


for ii=1:length(ICC)
    slice_ICC{1}{ii}         =xASL_vis_CropParmsApply(xASL_im_rotate(ICC{ii}(:,:,62).*GM_mask_visual(:,:,62),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_ICC{2}{ii}         =xASL_vis_CropParmsApply(xASL_im_rotate(ICC{ii}(:,:,76).*GM_mask_visual(:,:,76),90),ant_crop,size(CBF{ii},2)-pos_crop,L_crop,size(CBF{ii},1)-R_crop);
    slice_ICC{3}{ii}         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(ICC{ii}(:,63,:).*GM_mask_visual(:,63,:)),90),sup_crop,size(mean_GM,2)-inf_crop,L_crop,size(mean_GM,1)-R_crop);
    slice_ICC{4}{ii}         =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(ICC{ii}(:,95,:).*GM_mask_visual(:,95,:)),90),sup_crop,size(mean_GM,2)-inf_crop,L_crop,size(mean_GM,1)-R_crop);
    slice_ICC{5}{ii}         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( ICC{ii}(68,:,:) .*GM_mask_visual(68,:,:)),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);
    slice_ICC{6}{ii}         =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( ICC{ii}(76,:,:) .*GM_mask_visual(76,:,:)),sup_crop,size(CBF{ii},2)-inf_crop,ant_crop,size(CBF{ii},2)-pos_crop);end

for ii=1:2
    view_ICC_total{ii}           =[slice_ICC{1}{ii} slice_ICC{2}{ii} slice_ICC{3}{ii} slice_ICC{4}{ii} slice_ICC{5}{ii}  slice_ICC{6}{ii}];end

% figure(1);imshow(view_ICC_total{1},[-25 50],'Colormap',jet256);
% figure(2);imshow(view_ICC_total{2},[-25 50],'Colormap',jet256);
% % print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\max_ICC_overview.tiff']);
%% 4) Project ICC over GM masks
clear r_slice g_slice b_slice slice_view_ICC_clr_prep slice_ICC_color rout_slice
slice_view_ICC_clr_prep                                                =slice_ICC;

% 1) Rescale maps
min_intensity                                                           =50;
max_intensity                                                           =150;

for ii=1:6
        for j=1:2
        slice_view_ICC_clr_prep{ii}{j}(slice_view_ICC_clr_prep{ii}{j}<min_intensity)      =min_intensity;
        slice_view_ICC_clr_prep{ii}{j}(slice_view_ICC_clr_prep{ii}{j}>max_intensity)      =max_intensity;

        slice_ICC_color{ii}{j}                                                =round(double(rescale( slice_view_ICC_clr_prep{ii}{j},0,255,1)));
        slice_ICC_color{ii}{j}                                                =max(1,min(slice_ICC_color{ii}{j},size(jet256,1)));

        r_slice{ii}{j}     =zeros(size(slice_ICC_color{ii}{j})); r_slice{ii}{j}(:)  =jet256(slice_ICC_color{ii}{j},1);
        g_slice{ii}{j}     =zeros(size(slice_ICC_color{ii}{j})); g_slice{ii}{j}(:)  =jet256(slice_ICC_color{ii}{j},2);
        b_slice{ii}{j}     =zeros(size(slice_ICC_color{ii}{j})); b_slice{ii}{j}(:)  =jet256(slice_ICC_color{ii}{j},3);

        b_slice{ii}{j}( (GM_slices{ii}==1) & r_slice{ii}{j}==0 & g_slice{ii}{j}==0 & b_slice{ii}{j}==0)    =0.5313;

        rout_slice{ii}{j}        = zeros([size(r_slice{ii}{j}),3]);
        rout_slice{ii}{j}(:,:,1) = r_slice{ii}{j}.*slice_view_skull_clr{ii}(1:size(r_slice{ii}{j},1),:,1);
        rout_slice{ii}{j}(:,:,2) = g_slice{ii}{j}.*slice_view_skull_clr{ii}(1:size(r_slice{ii}{j},1),:,1);
        rout_slice{ii}{j}(:,:,3) = b_slice{ii}{j}.*slice_view_skull_clr{ii}(1:size(r_slice{ii}{j},1),:,1);

    %     slice_ICC_color{ii}{j}                         =(0.5.*rout_slice{ii}{j})+(0.5.*rout_slice{ii}{j}.*GM_slices_clr{ii});
        slice_ICC_color{ii}{j}                         =rout_slice{ii}{j};

    % Put GM map on background
    for i=1:size(slice_ICC_color{ii}{j},1)
        for k=1:size(slice_ICC_color{ii}{j},2)
            if  slice_ICC_color{ii}{j}(i,k,1)==0 & slice_ICC_color{ii}{j}(i,k,2)==0 & slice_ICC_color{ii}{j}(i,k,3)==0
                slice_ICC_color{ii}{j}(i,k,1)=0.8.*GM_background_clr{ii}(i,k,1);
                slice_ICC_color{ii}{j}(i,k,2)=0.8.*GM_background_clr{ii}(i,k,2);
                slice_ICC_color{ii}{j}(i,k,3)=0.8.*GM_background_clr{ii}(i,k,3);end;end;end
    end
end

for ii=1:2
    view_ICC_total{ii}     =[slice_ICC_color{1}{ii} slice_ICC_color{2}{ii} slice_ICC_color{3}{ii} slice_ICC_color{4}{ii} slice_ICC_color{5}{ii}  slice_ICC_color{6}{ii}];
    figure(ii);imshow([view_ICC_total{ii}]);
    print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\view_ICC_total_' num2str(ii) '.tiff']);
end
%% 4) Visualization variation_diff parametric, including vendor A study 1 vs. vendor A study 2
clear CBF_slices vendor subject diff_slices
for vendor=1:12
    for subject=1:size(CBF{vendor},4)
        CBF_slices{vendor}{1}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(CBF{vendor}(:,:,62,subject).*GM_mask_visual(:,:,62),90),ant_crop,size(CBF{vendor},2)-pos_crop,L_crop,size(CBF{vendor},1)-R_crop);
        CBF_slices{vendor}{2}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(CBF{vendor}(:,:,76,subject).*GM_mask_visual(:,:,76),90),ant_crop,size(CBF{vendor},2)-pos_crop,L_crop,size(CBF{vendor},1)-R_crop);
        CBF_slices{vendor}{3}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(CBF{vendor}(:,63,:,subject).*GM_mask_visual(:,63,:)),90),sup_crop,size(CBF{vendor},2)-inf_crop,L_crop,size(CBF{vendor},1)-R_crop);
        CBF_slices{vendor}{4}(:,:,subject)      =xASL_vis_CropParmsApply(xASL_im_rotate(FlipOrientation_isotropic(CBF{vendor}(:,95,:,subject).*GM_mask_visual(:,95,:)),90),sup_crop,size(CBF{vendor},2)-inf_crop,L_crop,size(CBF{vendor},1)-R_crop);
        CBF_slices{vendor}{5}(:,:,subject)      =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( CBF{vendor}(68,:,:,subject) .*GM_mask_visual(68,:,:) ),sup_crop,size(CBF{vendor},2)-inf_crop,ant_crop,size(CBF{vendor},2)-pos_crop);
        CBF_slices{vendor}{6}(:,:,subject)      =xASL_vis_CropParmsApply(FlipOrientation2_isotropic( CBF{vendor}(76,:,:,subject) .*GM_mask_visual(76,:,:) ),sup_crop,size(CBF{vendor},2)-inf_crop,ant_crop,size(CBF{vendor},2)-pos_crop);
    end
end

clear orislice diff_slices
for ii=1:6
    for orislice=1:6
        diff_slices{ii}{1,orislice}             =CBF_slices{1,ii*2}{1,orislice} - CBF_slices{ 1,(ii*2)-1}{1,orislice};
    end
end

clear orislice
for orislice=1:6 % vendor A vs. vendor A
    diff_slices{7}{1,orislice}                  =CBF_slices{1,3}{1,orislice}(:,:,11) - CBF_slices{ 1,7}{1,orislice}(:,:,11);
    diff_slices{8}{1,orislice}                  =CBF_slices{1,4}{1,orislice}(:,:,11) - CBF_slices{ 1,8}{1,orislice}(:,:,11);
end

% % Levene's test for comparison inter-vendor against intra-vendor
% clear vendor orislice P
% for orislice=1:6
%     for ii=1:size(diff_slices{1}{orislice},1)
%         for j=1:size(diff_slices{1}{orislice},2)
%             clear X
%             X( 1:11,1)              =squeeze(diff_slices{1}{orislice}(ii,j,:)); X( 1:11,2)=1;
%             X(12:22,1)              =squeeze(diff_slices{2}{orislice}(ii,j,:)); X(12:22,2)=1;
%             X(23:44,1)              =squeeze(diff_slices{3}{orislice}(ii,j,:)); X(23:44,2)=2;
%             [H{1}{orislice}(ii,j) p_value] =Levenetest_HJM(X,0.05);end;end;end
%
%
% % Levene's test
% clear vendor orislice P
% for orislice=1:6
%     for ii=1:size(diff_slices{1}{orislice},1)
%         for j=1:size(diff_slices{1}{orislice},2)
%             clear X
%             X( 1:14,1)              =squeeze(diff_slices{4}{orislice}(ii,j,:)); X( 1:14,2)=1;
%             X(15:28,1)              =squeeze(diff_slices{5}{orislice}(ii,j,:)); X(15:28,2)=1;
%             X(29:56,1)              =squeeze(diff_slices{6}{orislice}(ii,j,:)); X(29:56,2)=2;
%             [H{2}{orislice}(ii,j) p_value] =Levenetest_HJM(X,0.05);end;end;end


% clear vendor orislice P
% for orislice=1:6 % vendor A vs. vendor A
%     for ii=1:size(diff_slices{1}{orislice},1)
%         for j=1:size(diff_slices{1}{orislice},2)
%             clear X
%
%             X( 1:11,1)                      =squeeze(diff_slices{7}{orislice}(ii,j,:));     X( 1:11,2)=1;
%             X(12:22,1)                      =squeeze(diff_slices{8}{orislice}(ii,j,:));     X(12:22,2)=2;
%
%             [H{3}{orislice}(ii,j) p_value]  =Levenetest_HJM(X,0.05);
%         end
%     end
% end


 % save('C:\Backup\ASL\3CV\Analysis_temp\Levene_intra_vs_inter.mat','H');
load('C:\Backup\ASL\3CV\Analysis_temp\Levene_intra_vs_inter.mat');

% Create post&neg color masks
clear pos neg pos_mask neg_mask pne_mask parametric
for orislice=1:6
    for j=1:3
        H{j}{orislice}(isnan(H{j}{orislice}))   =0;

        pos{j}{orislice}                        =logical(slice_ICC{1,orislice}{j}>100);
        neg{j}{orislice}                        =logical(slice_ICC{1,orislice}{j}<100) & slice_ICC{1,orislice}{j}>0;

        for ii=1:3
            pos_mask{j}{orislice}(:,:,ii)       =double(pos{j}{orislice}.*H{j}{orislice});
            neg_mask{j}{orislice}(:,:,ii)       =double(neg{j}{orislice}.*H{j}{orislice});end

        neg_mask{j}{orislice}(:,:,1)            =0;
        neg_mask{j}{orislice}(:,:,2)            =0.5.*neg_mask{j}{orislice}(:,:,3);
        pos_mask{j}{orislice}(:,:,2)            =0;
        pos_mask{j}{orislice}(:,:,3)            =0;

        pne_mask{j}{orislice}                   =pos_mask{j}{orislice}+ neg_mask{j}{orislice};
        parametric{j}{orislice}                 =pne_mask{j}{orislice}+ (( GM_slices_clr{orislice} - pne_mask{j}{orislice} ).*0.5);
    end
end

for j=1:3
    overview_diff_parametric{j}                 =[parametric{j}{1} parametric{j}{2} parametric{j}{3} parametric{j}{4} parametric{j}{5} parametric{j}{6}];end

% figure(3);imshow(overview_diff_parametric{1});
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_diff_parametric_3.tiff']);
% figure(4);imshow(overview_diff_parametric{2});
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_diff_parametric_4.tiff']);
% figure(5);imshow(overview_diff_parametric{3});
% print(gcf,'-dtiff','-r600',['C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\overview_diff_parametric_5.tiff']);
%% 4) ratio histograms
clear bin_nr min_nr max_nr bin_size

bin_nr      =100;
min_nr      =0;
max_nr      =200;
bin_size    =(max_nr-min_nr)/bin_nr;
%myfilter    =fspecial('gaussian',[bin_nr,1],1*0.02*bin_nr);
% 100 bins / 200 ratio % range = 0.5 bins / ratio %, 160 ratio % used which would have been 80 bins

% Ratio histograms

for ii=1:2
    clear temp N_ROI X_ROI
    temp                        =ICC{ii};

    [N_ROI  X_ROI]              =hist(temp( mean_segm_mask{1}>0.5 & isfinite(temp)),[min_nr:(max_nr-min_nr)/bin_nr:max_nr]);
    N_ROI                       =N_ROI./sum(N_ROI)./bin_size;
    %N_ROI                       =imfilter(N_ROI, myfilter', 'replicate');
	N_ROI                       = xASL_im_ndnanfilter(N_ROI,'gauss',[0 1*0.02*bin_nr*2.335 0],0);

    figure(ii);plot(X_ROI,100.*N_ROI,'k');
    axis([20 180 0 2]);
    set( gca, 'Color','none');
%     set(gca,'xtick',[],'ytick',[],'Color','none');
%     eval(['export_fig C:\Backup\Copy\Copy\AMC\Studies\3CV\Article\Figures\hist_ratio_study_' num2str(ii) '.eps -eps -transparent -CMYK -painters']);
%     close all
end