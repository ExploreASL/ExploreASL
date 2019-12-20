function [ImOut]   = xASL_im_ExtrapolateOverNaNs(ImIn)
%xASL_im_ExtrapolateOverNaNs Smoothing with a small kernel [2 2 2] goes much faster than big
% kernel [8 8 8], hence for low quality, or outside the mask, use the
% [2 2 2]


fprintf('%s','Extrapolating over NaNs...  ');

% Rough estimation MaxIt for timer:
% average of [8 8 8] & [2 2 2] = [6 6 6] = 6^3
MaxIt   = 5; % sum(sum(sum(isnan(ImIn)))) / 6^3;

iCTotal     = 0;
maxItTotal = MaxIt * size(ImIn,4) * size(ImIn,5) * size(ImIn,6);

VoxelSize 				= [1.5 1.5 1.5];
ImOut = zeros(size(ImIn));

for iX=1:size(ImIn,4)
    for iY=1:size(ImIn,5)
        for iZ=1:size(ImIn,6)
            TempIm = ImIn(:,:,:,iX,iY,iZ);
			iC = iCTotal + 1;
			
            while   sum(sum(sum(isnan(TempIm))))>0 % loop extrapolation until full FoV is filled
                xASL_TrackProgress(iC, maxItTotal);
                TempIm = xASL_im_ndnanfilter(TempIm,'gauss',double([12 12 12]./VoxelSize),2); % 2 at the end means extrapolation only
                iC=iC+1;
			end
			iCTotal = iCTotal + MaxIt;
            ImOut(:,:,:,iX,iY,iZ) = TempIm;
        end
    end
end

fprintf('\n');


end
