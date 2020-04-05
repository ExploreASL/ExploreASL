function DiceCoeff = xASL_im_ComputeDice(imA, imB)
%xASL_im_ComputeDice Calculate Dice coefficient of image overlap

if ~islogical(imA)
    warning('Input image was not dichotomous, dichotomizing now');
    imA = imA>0;
end
if ~islogical(imB)
    warning('Input image was not dichotomous, dichotomizing now');
    imB = imB>0;
end

DiceCoeff(1)   = (2 .* sum(sum(sum(imA & imB)))) ./ (sum(imA(:)) + sum(imB(:)) );


end

