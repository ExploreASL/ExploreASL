function [diff_view_mean H_ttestCONTRAST x] = xASL_spm_GLMcreateStatsFigures(x)
%xASL_spm_GLMcreateStatsFigures Takes result from SPM and puts it in nice figures


%% First convert datasets into correct format

for iM=1:length(x.S.DATASETS)
    if  size(x.S.DATASETS{iM},3)==1
        x.S.DATASETS{iM}     = xASL_im_Column2IM(x.S.DATASETS{iM},x.S.VBAmask);
    end
end

view_ASL = xASL_im_TransformData2View( x.S.DATASETS, x );

for iSet=1:x.S.nSets
    view_ASL{iSet}(view_ASL{iSet}<0)        = 0; % clip values <0
end

% Save contrast map for later
diff_view_mean{1}       = (x.S.ContrastMap.^2).^0.5; % abs T-stat map
diff_view_mean{1}       = xASL_im_TransformData2View( diff_view_mean{1}, x );
H_ttestCONTRAST         = xASL_im_TransformData2View( x.S.MaskMap>0, x );

if  x.S.GlobalNormalization
    for ii=1:length(view_ASL)
        for iSubj=1:size(view_ASL{ii},3)
            view_ASL{ii}(:,:,iSubj)     = 100.* (view_ASL{ii}(:,:,iSubj) ./ xASL_stat_MeanNan( xASL_stat_MeanNan( view_ASL{ii}(:,:,iSubj) ) ) );
        end
    end
end

% For selected significant voxels, get mean difference or CoV
% for diff_view_mean{1}. diff_view_mean{1} is t-stat or f-stat (see
% above)
if      x.S.nSets==1 && x.S.RegressionCOVAR
        FList                       = xASL_adm_GetFileList( x.S.SPMdir, '^beta_\d{4}\.(nii|nii\.gz)$');
        for FL=1:length(FList)
            DIFFMap                 = xASL_io_Nifti2Im( FList{FL} ); % First diff_view_mean is t-stat or f-stat map
            diff_view_mean{FL+1}    = xASL_im_TransformData2View( DIFFMap, x );

            clear DIFFMap
        end
        clear FList FL

elseif  x.S.nSets==1 % 1-sample t-test
        diff_view_mean{2}           = xASL_stat_MeanNan(view_ASL{1},3); % simply show mean CBF map

elseif  x.S.nSets>2
        AllSet(:,:,1:size(view_ASL{1},3))   = view_ASL{1};
        for iSet3=2:x.S.nSets
            AllSet(:,:,size(AllSet,3)+1:size(AllSet,3)+size(view_ASL{iSet3},3))    = view_ASL{iSet3};
        end
        diff_view_mean{2}   = 100.*(xASL_stat_StdNan(AllSet,0,3)./xASL_stat_MeanNan(AllSet,3)); % Coefficient of Variance
else
        diff_view_mean{2}   = xASL_stat_MeanNan(view_ASL{2},3) - xASL_stat_MeanNan(view_ASL{1},3);
end




H_ttestCONTRAST(isnan(H_ttestCONTRAST))             = 0;
for iD=1:length(diff_view_mean)
    diff_view_mean{iD}(isnan(diff_view_mean{iD}))   = 0;
end

end
