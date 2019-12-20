function image_out = xASL_im_FlipOrientation(image_in)
%FlipOrientation This function flips the 3 dimensions from sagittal to
%transversal or tra to cor. Leaves other dimensions untouched.


image_out   = shiftdim(xASL_im_rotate(image_in,90),1);

% OLD CODE
% temp=single(xASL_im_rotate(image_in,90));
% dim=size(temp);
% image_out=zeros(dim(2),dim(3),dim(1),size(temp,4));
% 
% 
% 
% for l=1:size(temp,4)
%     for i=1:dim(2)
%         for j=1:dim(3)
%             for k=1:dim(1)
%                 image_out(i,j,k,l)=temp(k,i,j,l);
%             end
%         end
%     end
% end


%new_y=old_z dim(3)->dim(2)
%new_x=old_y dim(2)->dim(1)
%new_z=old_x dim(1)->dim(3)



end

