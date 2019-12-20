function DiceCoeff = xASL_im_ComputeDice(imA,imB)
%xASL_im_ComputeDice Calculate Dice coefficient of image overlap

DiceCoeff(1)   = (2 .* sum(sum(sum(imA & imB)))) ./ (sum(imA(:)) + sum(imB(:)));


end

