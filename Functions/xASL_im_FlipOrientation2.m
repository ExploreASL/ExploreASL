function image_out = xASL_im_FlipOrientation2(image_in)
%FlipOrientation This function flips the 3 dimensions from sagittal to
%cor or tra to sag. Leaves other dimensions untouched.


image_out               = xASL_im_rotate(shiftdim(xASL_im_rotate(image_in,90),2),180);
image_out               = image_out(:,:,size(image_out,3):-1:1); % flip


% 
% 
% temp=xASL_im_rotate(single(image_in),90);
% dim=size(temp);
% image_out=zeros(dim(3),dim(1),dim(2),size(temp,4));
% 
% for l=1:size(temp,4)
%     for i=1:dim(3)
%         for j=1:dim(1)
%             for k=1:dim(2)
%                 image_out(i,j,k,l)=temp(j,k,i,l); % new z-direction is magnified/blurred by 2. 2*i-1 = old i & 2*i is average of old i & old (i+1)
%             end
%         end
%     end
% end
% 
% image_out=xASL_im_rotate(image_out,180);
% 
%new_y=old_z dim(3)->dim(2)
%new_x=old_y dim(2)->dim(1)
%new_z=old_x dim(1)->dim(3)



end

