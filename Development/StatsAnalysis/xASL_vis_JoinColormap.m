function joinedColormap = xASL_io_joinColormap(blackMiddle,colormap1,colormap2)
% Take two colormaps and put them together, creating a colormap with 256 values
% input colormaps 1 and 2 are stacked as colormap1 first, then colormap2
% They can be of anylength [N1,3] and [N2,3]
% joinedColormap has size [256,3]

% Create an empty map
joinedColormap = zeros(256,3);

% Fill the middle values with black, or make a smooth transition
if blackMiddle
	% Specify the coordinates of the output
	Nout1 = 1:127;
	Nout2 = 130:256;
else
	Nout1 = 1:128;
	Nout2 = 129:256;
end

joinedColormap(Nout1,:) = interp2(colormap1,1:3,((Nout1-1)*(255/(Nout1(end)-1)) + 1)');
joinedColormap(Nout2,:) = interp2(colormap2,1:3,((Nout2-Nout2(1))*(255/(Nout2(end)-Nout2(1))) + 1)');

return
