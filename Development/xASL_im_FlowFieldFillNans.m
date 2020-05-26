function imOut = xASL_im_FlowFieldFillNans(imIn)
% Fill in the NaNs at the borders of the flow field

imOut = imIn;
% Zeros are also converted to NaNs - because this is sometimes unclear in the flow fields.
imOut(imOut==0) = NaN;

% The first layer next to NaNs is sometimes intrapolated, so rather assigned to Nans as well
[a,~] = xASL_im_DistanceTransform(isnan(imOut));
imOut(a<2) = NaN;

nX = size(imIn,1);
nY = size(imIn,2);
nZ = size(imIn,3);
% Applies the same along the 4th dimension
for tt = 1:size(imIn,4)
	% Applies three times along the remaining dimensions to make sure all the corners are filled
	for ii = 1:3
		if sum(isnan(imOut(:)))
			% Does it along the remaining three dimensions
			for dd = 1:min(3,ndims(imOut))
				% Extracts column wise values
				switch(dd)
					case 1
						imCol = reshape(imOut(:,:,:,tt),[nX,nY*nZ]);
					case 2
						imCol = reshape(shiftdim(imOut(:,:,:,tt),1),[nY,nZ*nX]);
					case 3
						imCol = reshape(imOut(:,:,:,tt),[nX*nY,nZ])';
				end
				
				% Find those columns that incude NaNs, but there are more than 50% non-NaNs (only apply this condition on the first run)
				if ii>1
					indCol = sum(isnan(imCol),1)>0;
				else
					indCol = (sum(isnan(imCol),1)>0 & (sum(~isnan(imCol),1)>(0.5*size(imCol,1))));
				end
					
				nCol = size(imCol,1);
				nColHalf = floor(nCol/2);
				% Find columns starting and ending with nans
				indColLow  = find(indCol & isnan(imCol(1,:)));
				indColHigh = find(indCol & isnan(imCol(end,:)));
				
				% For each column starting with Nans
				for cc = indColLow
					minInd = find(~isnan(imCol(:,cc)),1,'first');
					diff = sum(imCol(minInd:nColHalf,cc) - imCol((minInd+1):(nColHalf+1),cc))/(nColHalf-minInd+1);
					if isnan(diff)
						diff = xASL_stat_MeanNan(imCol(minInd:nColHalf,cc) - imCol((minInd+1):(nColHalf+1),cc));
					end
					imCol(1:(minInd-1),cc) = imCol(minInd,cc) + ((minInd-1):-1:1)*diff;
				end
				
				% For each column ending with NaNs
				for cc = indColHigh
					maxInd = find(~isnan(imCol(:,cc)),1,'last');
					diff = sum(imCol(nColHalf:maxInd,cc)-imCol((nColHalf-1):(maxInd-1),cc))/(maxInd-nColHalf+1);
					if isnan(diff)
						diff = xASL_stat_MeanNan(imCol(nColHalf:maxInd,cc)-imCol((nColHalf-1):(maxInd-1),cc));
					end
					imCol((maxInd+1):end,cc) = imCol(maxInd,cc) + (1:(nCol-maxInd))*diff;
				end
				
				% Put back together to the original image
				switch(dd)
					case 1
						imOut(:,:,:,tt) = reshape(imCol,[nX,nY,nZ]);
					case 2
						if ndims(imOut)<3
							imOut(:,:,:,tt) = shiftdim(reshape(imCol,[nY,nX]),1);
						else
							imOut(:,:,:,tt) = shiftdim(reshape(imCol,[nY,nZ,nX]),2);
						end
					case 3
						imOut(:,:,:,tt) = reshape(imCol',[nX,nY,nZ]);
						
				end
			end
		end
	end
end

% Any remaining NaNs had to be inside, so we convert them back to NaN
imOut(isnan(imOut)) = 0;

end