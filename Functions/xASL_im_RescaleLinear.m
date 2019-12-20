function [ NewMatrix ] = xASL_im_RescaleLinear(OriMatrix,NewMin,NexMax,NonZerosOption)
%rescale Linearly rescales input matrix to output matrix,
% applying a new minimum and new maximum.

% % Test values
% % old_matrix      =round(rand(5,5).*10);
% % new_minimum     =2;
% % new_maximum     =6;
% % 

OriMatrix  = double( OriMatrix );

% If matrix contains zeros only, reinforce nonzerosoption==0
if  max(OriMatrix(:))==0
    NonZerosOption=0;
end

% define minimum old matrix
if      NonZerosOption==1
        OriMin         = min(nonzeros(OriMatrix(:)));
else
        OriMin         = min(OriMatrix(:));
end


OriMax                 = max(OriMatrix(:));
OriRange               = OriMax-OriMin;
NewRange               = NexMax-NewMin;

NewMatrix              = NewMin+NewRange.* (OriMatrix-OriMin)./OriRange;


end

