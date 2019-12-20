%% QIBA wsCV analysis

PopulationName = {'NOVICE_childrenHIV' 'Sleep_adultsHealthy' 'EPAD_elderly'};

StudiesDir{1} = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_NOVICE\Population';
StudiesDir{2} = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_Sleep\Population';
StudiesDir{3} = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population';

% 1) Compute volumes HO regions
ROIname = {'HOcort_CONN' 'HOsub_CONN'};
for iRegions=1:2
    ROIpath{iRegions} = fullfile('C:\ExploreASL\Maps\Atlases',[ROIname{iRegions} '.nii']);
end

%% Volumes
for iStudy=1:3
    clear VolROI
    for iRegions=1:2
        PathGM = fullfile(StudiesDir{iStudy}, 'Templates', 'pGM_bs-mean.nii');
        GMim = xASL_io_Nifti2Im(PathGM);
        GMmask = GMim>0.7.*max(GMim(:));
        
        ROIim = xASL_io_Nifti2Im(ROIpath{iRegions});
        for iROI=1:max(ROIim(:))
            tIM = ROIim==iROI;
            tIM = tIM.*GMmask;
            VolROI{iRegions}(iROI,1) = sum(sum(sum(tIM))).*1.5^3 / 1000; % in mL
        end
    end
    % Combine HOcort & HOsubcort
    VolFull(1:48) = VolROI{1};
%     VolFull(49:48+8) = VolROI{2}; % WITHOUT SUBCORTICAL
%     !!!!!!!!!!!!!!!!!!!!!

    %% wsCV
    for iRegions=1:2
        clear ReproValues ReproNum
        for iPart=1:2
            PathPWI{iPart} = xASL_adm_GetFileList(fullfile(StudiesDir{iStudy},'Stats'),['^median_PWI_part' num2str(iPart) '_' ROIname{iRegions} '.*\.tsv$'],'FPList');
            if length(PathPWI{iPart})~=1
                error('Wrong length PathPWI{iPart}');
            end
            [~, tempCell] = xASL_adm_csv2tsv(PathPWI{iPart}{1});
            nVols = size(VolROI{iRegions},1);
            ReproValues{iPart} = tempCell(3:end,end-nVols*3+1:end);
            % take left & right separately here
            WhichColumns = repmat([0 1 1],1,nVols);
            ReproValues{iPart} = ReproValues{iPart}(:,logical(WhichColumns));
            for iX=1:size(ReproValues{iPart},1)
                for iY=1:size(ReproValues{iPart},2)
                    ReproNum{iPart}(iX,iY) = xASL_str2num(ReproValues{iPart}(iX,iY));
                end
            end
        end
        % Compute wsCV for left & right separately here
        meanValues{iRegions} = xASL_stat_MeanNan((ReproNum{1}+ReproNum{2})./2,1);
        SDdiff{iRegions} = xASL_stat_StdNan(ReproNum{1}-ReproNum{2},[],1);
        wsCV{iRegions} = SDdiff{iRegions}./meanValues{iRegions}.*100;
        wsCV{iRegions} = (wsCV{iRegions}(1:2:end-1) + wsCV{iRegions}(2:2:end))./2; % average LeftRight
    end
    % Combine HOcort & HOsubcort
    wsCVfull(1:48) = wsCV{1};
%     wsCVfull(49:48+8) = wsCV{2}; % WITHOUT SUBCORTICAL
%     !!!!!!!!!!!!!!!!!!!!!
    
    wsCVfull = wsCVfull./(2^0.5); % correct for having half scan-time
    
    figure(iStudy);
    plot(VolFull,wsCVfull,'.');
    xlabel('Volume (mL)');
    ylabel('wsCV (%)');
    title(PopulationName{iStudy});
end
        
        