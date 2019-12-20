% Reads all slices and repetitions of an Ecat file

function [im mhead shead] = readEcat(filename)
% Get the information about the size of the volume

infoVol = mat_read_mainh(filename);
infoImg = mat_read_subh(filename,[1,1,1,0,0]);
    
% Create an empty matrix to fill in the image data
im = zeros([infoImg.y_dimension,infoImg.x_dimension,infoImg.z_dimension,infoVol.num_frames,max(infoVol.num_gates,1)]);
    
% Fill the data for all the times and slices
for frame = 1:infoVol.num_frames
	for gate = 1:max(infoVol.num_gates,1)
		[infoTmp, imTmp] = mat_read_volume(filename,[frame,1,gate,0,0]);
		% The first gate is ignored as the mean image was saved there
		im(:,:,:,frame,gate) = imTmp;
	end
end

shead = infoImg;
mhead = infoVol;
return;