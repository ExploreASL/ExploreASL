function [xmin xmax ymin ymax] = xASL_vis_CropParmsAcquire(temp_image)
%xASL_vis_CropParmsAcquire Goes from outside to inside to acquire crop settings.
% Works with grayscale images (2 dimensions per slice).image position information (2D matrix) should be first 2
% dimensions. Could include colordimension later on.

dim=size(temp_image);
if  length(dim)~=3
    dim(3)=1;
end    


for k=1:dim(3)
    image=temp_image(:,:,k);


    skip=0;

    % xmin(k)
    xmin(k)=1;
    encountered=0;
    for i=1:dim(1)
        if  sum(image(i,:))==0 && encountered==0
            xmin(k)=xmin(k)+1;
        else
            encountered=1;
        end
    end

    % ymin(k)
    ymin(k)=1;
    encountered=0;
    for i=1:dim(2)
        if  sum(image(:,i))==0 && encountered==0
            ymin(k)=ymin(k)+1;
        else
            encountered=1;
        end
    end

    % xmax(k)
    xmax(k)=dim(1);
    encountered=0;
    for i=1:dim(1)
        j=dim(1)-i+1;
        if  sum(image(j,:))==0 && encountered==0
            xmax(k)=xmax(k)-1;
        else
            encountered=1;
        end
    end

    % ymax(k)
    ymax(k)=dim(2);
    encountered=0;
    for i=1:dim(2)
        j=dim(2)-i+1;
        if  sum(image(:,j))==0 && encountered==0
            ymax(k)=ymax(k)-1;
        else
            encountered=1;
        end
    end

    % skip image if empty (only zeros)
    if  xmax(k)==1 && ymax(k)==1 && xmin(k)==dim(1) && ymin(k)==dim(2)
        skip=1;
    end

    % make results into square
    x_width=xmax(k)-xmin(k)+1;
    y_width=ymax(k)-ymin(k)+1;

    if      y_width>x_width
            diff=y_width-x_width;
            xmin(k)=xmin(k)-ceil(diff/2);
            xmax(k)=xmax(k)+(diff-ceil(diff/2));
    elseif  x_width>y_width
            diff=x_width-y_width;
            ymin(k)=ymin(k)-ceil(diff/2);
            ymax(k)=ymax(k)+(diff-ceil(diff/2));
    end


    if  skip==1
        xmin(k)=1;
        ymin(k)=1;
        xmax(k)=dim(1);
        ymax(k)=dim(2);
    end

end

xmin=min(xmin);
ymin=min(ymin);
xmax=max(xmax);
ymax=max(ymax);


end

