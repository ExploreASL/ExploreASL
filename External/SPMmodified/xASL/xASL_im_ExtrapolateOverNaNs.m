function [ImOut]   = xASL_im_ExtrapolateOverNaNs(ImIn)
% It removes the NaNs by extrapolating the non-NaN values over the NaN values.
% In case there are only NaNs in the image - i.e. there's nothing to extrapolate - then replaces the whole image by zeros.

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
			
			if sum(isnan(TempIm(:))) < numel(TempIm)
				while sum(sum(sum(isnan(TempIm))))>0 % loop extrapolation until full FoV is filled
					xASL_TrackProgress(iC, maxItTotal);
					TempIm = xASL_im_ndnanfilter(TempIm,'gauss',double([12 12 12]./VoxelSize),2); % 2 at the end means extrapolation only
					iC=iC+1;
				end
			else
				% Only NaNs - replace all by zeros
				TempIm(isnan(TempIm)) = 0;
			end
			iCTotal = iCTotal + MaxIt;
            ImOut(:,:,:,iX,iY,iZ) = TempIm;
        end
    end
end

fprintf('\n');


end
