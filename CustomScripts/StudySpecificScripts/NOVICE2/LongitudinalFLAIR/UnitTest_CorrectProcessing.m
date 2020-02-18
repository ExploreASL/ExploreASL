% Unittest

SubjSetInd = find(strcmp(x.S.SetsName,'SubjectNList'));
SubjectList = x.S.SetsID(:,SubjSetInd);

TPSetInd = find(strcmp(x.S.SetsName,'LongitudinalTimePoint'));
TPList = x.S.SetsID(:,TPSetInd);

for iSubject=1:length(unique(SubjectList))
    xASL_TrackProgress(iSubject,length(unique(SubjectList)));
    % find TPs
    CurrentSubjIndex = find(SubjectList==iSubject);
    CurrentTPIndex = TPList(CurrentSubjIndex);
    
    % only process if multiple TPs present & contain first timepoint
    if ~(length(CurrentSubjIndex)>1 && length(CurrentTPIndex)>1 && min(CurrentTPIndex)==1)
        warning([num2str(iSubject) ' didnt have 2 TPs']);
        fprintf('\n\n\n');
        continue;
    end
        
    % load images (assuming there are only 2 TPs)
    clear PathIs IM
    for ii=1:2
        PathIs{ii} = fullfile(x.ROOT, x.SUBJECTS{CurrentSubjIndex(ii)},'WMH_SEGM.nii');
        IM{ii} = xASL_io_Nifti2Im(PathIs{ii});
        
        if length(unique(IM{ii}(:)))>2 % should only have dichotomous lesion mask
            warning(['Not binary mask ' PathIs{ii}]);
            fprintf('\n\n\n');
        end
    end
    if sum(sum(sum(IM{1} & ~IM{2})))~=0 % there should not be lesions on TP1 that dont exist on TP2    
        warning(['Lesions on TP1 that dont exist on TP2' PathIs{1}]);
        fprintf('\n\n\n');
    end
end    
    





% sum(sum(ImageIs(33:34,48:50,55)))
% 
% 
% IMcheck = xASL_im_Column2IM(IM,x.WBmask);
% IMcheck = IMcheck(33:34,48:50,55,:);
% IMcheck = squeeze(sum(sum(sum(IMcheck,1),2),3))
% 
% 
% dip_image(squeeze(IMcheck(:,:,55,:))>0)