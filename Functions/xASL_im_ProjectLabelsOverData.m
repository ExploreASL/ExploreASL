function OutputIM = xASL_im_ProjectLabelsOverData(DataIM,LabelIM,x,ScaleFactorData,ScaleFactorLabel)
%xASL_im_ProjectLabelsOverData This script projects labels over an image,
% but works only in 2D. Make sure to make a 2D image from a 3D or 4D image
% using xASL_vis_TransformData2View.m
% can be used in combination with xASL_vis_Imwrite.m

if ~exist('ScaleFactorData','var')
    ScaleFactorData     = 0.3;
end
if ~exist('ScaleFactorLabel','var')
    ScaleFactorLabel     = 1;
end

LabelIM                             = xASL_im_CreateLabelImage(LabelIM,x);

% robust scale DataIM
SortInt                             = sort(DataIM(isfinite(DataIM)));
MaxInt                              = SortInt(round(length(SortInt).*0.95));
DataIM(isnan(DataIM))               = 0;
DataIM(DataIM<0)                    = 0;
DataIM(DataIM>MaxInt)               = MaxInt;
DataIM                              = (DataIM+eps)./max(DataIM(:));
DataIM                              = ind2rgb(round(DataIM.*255),x.S.gray); % the 0.1 & 10 is to avoid clipping
OutputIM                            = DataIM.*ScaleFactorData+LabelIM.*ScaleFactorLabel.*0.25;

end














function [NewIM] = xASL_im_CreateLabelImage( IM, x, ColorShades )
%xASL_im_CreateLabelImage Takes label image from e.g. SPM & turns it into colors
% to be used for Figure creation
% By HJMM Mutsaerts, ExploreASL 2016

%% Third argument ColorShades is useful when you want
% multiple regions with the same colors but different brightness
% This can be useful e.g. when showing bilateral regions that you are
% merging for the post-hoc ROI analysis
% This argument can be omitted

%% Convert image in consequential integer labels
NewMatrix       = zeros(size(IM));
UniqueLabels    = unique(nonzeros(IM(:)));

for iU=1:length(UniqueLabels)
    NewMatrix(IM==UniqueLabels(iU))     = iU;
end

IM              = NewMatrix;

%% If ColorShades are requested, this part will run
% Third argument ColorShades is useful when you want
% multiple regions with the same colors but different brightness
% This can be useful e.g. when showing bilateral regions that you are
% merging for the post-hoc ROI analysis
% This argument can be omitted

ShadeFactor     = 1;

if  exist('ColorShades','var')
    for iC=1:length(ColorShades)

        FirstClusterColor   = x.S.LabelClr(ColorShades{iC}(1),:);
        nShades             = length(ColorShades{iC});
        Shades              = ((1-FirstClusterColor)./nShades)./ShadeFactor;
        
        for iC2=2:length(ColorShades{iC})
            CurrentColor                    = min([1 1 1],FirstClusterColor+(iC2-1).*Shades);
            x.S.LabelClr(ColorShades{iC}(iC2),:)= CurrentColor;
        end
    end
end

%% Administration
if  length(size(IM))~=2
    error('Script currently only works for images with 2 dimensions');
end
    
if  max(IM(:))>size(x.S.LabelClr,1)
    warning([num2str(size(x.S.LabelClr,1)) ' colors available, but ' num2str(max(IM(:))) ' colors required, using colors multiple times']);
    % NB: check whether number of labels is not more than colors specified above
end


%% Create new image

NewIM   = repmat(zeros(size(IM)),[1 1 3]);
OldIM   = repmat(IM,[1 1 3]);

% assign colors to labels
for iLabel=1:max(IM(:))
    % Choose color
    iColor  = iLabel;
    if  iColor>size(x.S.LabelClr,1)
        % start re-using colors, by multiplication
        while   iColor>size(x.S.LabelClr,1)
                % simply each time subtract total number of colors, which
                % will start over, but still keep unique values for how
                % many values there are in x.S.LabelClr
                iColor      = iColor-size(x.S.LabelClr,1);
        end
    end            
                    
    % create full label image
    LabelIM                 = repmat(zeros(size(IM)),[1 1 3]);
    for iL=1:3
        LabelIM(:,:,iL)     = x.S.LabelClr(iColor,iL);
    end
    
    NewIM(OldIM==iLabel)    = LabelIM(OldIM==iLabel);
end

end