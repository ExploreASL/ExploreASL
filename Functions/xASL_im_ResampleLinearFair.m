function [output_res]=xASL_im_ResampleLinearFair(im_input,newsize,showWaitbar)
% xASL_im_ResampleLinearFair
%
% FORMAT:       [output_res]=xASL_im_ResampleLinearFair(im_input,newsize)
% 
% INPUT:        im_input    - Image matrix (REQUIRED, DOUBLE, SINGLE or INT)
%               newsize     - Size of ouput image (REQUIRED, INTEGER ARRAY)
%               showWaitbar - Show waitbar (OPTIONAL, DEFAULT = true)
%
% OUTPUT:       output_res  - Resampled image matrix (SINGLE)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Downsample (or upsample, works similarly) old_res image to
%               low_res image, trilinear.
%
%               {{NB:}} new_res should fit exactly integer fold in old_res
%
%               {{NB:}} all dimensions of new_res should have equal size
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      image = xASL_io_Nifti2Im('.\M0.nii.gz');
%               output = xASL_im_ResampleLinearFair(image,[2,2,2],false);
%
% __________________________________
% Copyright 2015-2021 ExploreASL

% Input check
if nargin < 1 || nargin < 2 || ~exist('im_input','var') || ~exist('newsize','var')
    error('Missing input arguments...');
end
if nargin < 3 || ~exist('showWaitbar','var')
    showWaitbar = true;
end

% Calculation
old_res          =single(im_input);
new_res          =zeros(newsize(1),newsize(2),newsize(3));
 
if newsize(1)~=size(im_input,1) && newsize(2)~=size(im_input,2) && newsize(3)~=size(im_input,3)
if showWaitbar
    h = waitbar(1,'Initializing resample_min...');
end
%% 3-dimensional change
for k=1:size(new_res,3)
    for j=1:size(new_res,2) % prepare new_res blocks
        for i=1:size(new_res,1) % prepare new_res blocks
 
            % Display progress
            if showWaitbar
                percentage      =( (k-1)*size(new_res,2)*size(new_res,1)+(j-1)*size(new_res,1)+i )   /   (size(new_res,3) * size(new_res,2) * size(new_res,1) );
                h = waitbar(100*round(percentage),h,['resample ' num2str(100*round(percentage)) '% completed']);
            end

            start1              =(i-1)* (size(old_res,1) / size(new_res,1));
            end1                =i    * (size(old_res,1) / size(new_res,1));
            start2              =(j-1)* (size(old_res,2) / size(new_res,2));
            end2                =j    * (size(old_res,2) / size(new_res,2));
            start3              =(k-1)* (size(old_res,3) / size(new_res,3));
            end3                =k    * (size(old_res,3) / size(new_res,3));            
            
 
            sum                 =0;
            weigh               =0;
 
        if  floor(start3)<ceil(start3)                                      % if there is a head3 (heading partial voxel)
        h_index3            =ceil(start3);                                  % head index number
        h_weigh3            =ceil(start3)-start3;                           % voxelnr* weighing            
 
 
                if  floor(start2)<ceil(start2)                                      % if there is a head2 (heading partial voxel)
                h_index2            =ceil(start2);                                  % head index number
                h_weigh2            =ceil(start2)-start2;                           % voxelnr* weighing
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),h_index2,h_index3);    % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*h_weigh2*h_weigh3) ;
                        weigh           =weigh      + (h_weigh1*h_weigh2*h_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,h_index2,h_index3);        % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*h_weigh2*h_weigh3) ;
                        weigh           =weigh      + (b_weigh1*h_weigh2*h_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),h_index2,h_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*h_weigh2*h_weigh3);
                        weigh           =weigh      + (t_weigh1*h_weigh2*h_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end
 
 
                if  ceil(start2)<floor(end2)        % if there is a body2
                for b_index2    =ceil(start2)+1:floor(end2)
                b_weigh2        =1;
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),b_index2,h_index3);             % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*b_weigh2*h_weigh3) ;
                        weigh           =weigh      + (h_weigh1*b_weigh2*h_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,b_index2,h_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*b_weigh2*h_weigh3) ;
                        weigh           =weigh      + (b_weigh1*b_weigh2*h_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),b_index2,h_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*b_weigh2*h_weigh3);
                        weigh           =weigh      + (t_weigh1*b_weigh2*h_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end;end
 
                if  floor(end2)<ceil(end2)          % if there is a tail2
                t_index2        =ceil(end2);
                t_weigh2        =(end2-floor(end2));
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),t_index2,h_index3);             % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*t_weigh2*h_weigh3) ;
                        weigh           =weigh      + (h_weigh1*t_weigh2*h_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,t_index2,h_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*t_weigh2*h_weigh3) ;
                        weigh           =weigh      + (b_weigh1*t_weigh2*h_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),t_index2,h_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*t_weigh2*h_weigh3);
                        weigh           =weigh      + (t_weigh1*t_weigh2*h_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end
                clear h_index2 h_weigh2 b_index2 b_weigh2 t_index2 t_weigh2;end
                
        if  ceil(start3)<floor(end3)        % if there is a body3
        for b_index3    =ceil(start3)+1:floor(end3)
        b_weigh3        =1;                
                
 
                if  floor(start2)<ceil(start2)                                      % if there is a head2 (heading partial voxel)
                h_index2            =ceil(start2);                                  % head index number
                h_weigh2            =ceil(start2)-start2;                           % voxelnr* weighing
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),h_index2,b_index3);    % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*h_weigh2*b_weigh3) ;
                        weigh           =weigh      + (h_weigh1*h_weigh2*b_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,h_index2,b_index3);        % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*h_weigh2*b_weigh3) ;
                        weigh           =weigh      + (b_weigh1*h_weigh2*b_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),h_index2,b_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*h_weigh2*b_weigh3);
                        weigh           =weigh      + (t_weigh1*h_weigh2*b_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end
 
 
                if  ceil(start2)<floor(end2)        % if there is a body2
                for b_index2    =ceil(start2)+1:floor(end2)
                b_weigh2        =1;
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),b_index2,b_index3);             % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*b_weigh2*b_weigh3) ;
                        weigh           =weigh      + (h_weigh1*b_weigh2*b_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,b_index2,b_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*b_weigh2*b_weigh3) ;
                        weigh           =weigh      + (b_weigh1*b_weigh2*b_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),b_index2,b_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*b_weigh2*b_weigh3);
                        weigh           =weigh      + (t_weigh1*b_weigh2*b_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end;end
 
                if  floor(end2)<ceil(end2)          % if there is a tail2
                t_index2        =ceil(end2);
                t_weigh2        =(end2-floor(end2));
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),t_index2,b_index3);             % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*t_weigh2*b_weigh3) ;
                        weigh           =weigh      + (h_weigh1*t_weigh2*b_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,t_index2,b_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*t_weigh2*b_weigh3) ;
                        weigh           =weigh      + (b_weigh1*t_weigh2*b_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),t_index2,b_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*t_weigh2*b_weigh3);
                        weigh           =weigh      + (t_weigh1*t_weigh2*b_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end
                clear h_index2 h_weigh2 b_index2 b_weigh2 t_index2 t_weigh2;end;end
            
        if  floor(end3)<ceil(end3)          % if there is a tail3
        t_index3        =ceil(end3);
        t_weigh3        =(end3-floor(end3));            
            
                if  floor(start2)<ceil(start2)                                      % if there is a head2 (heading partial voxel)
                h_index2            =ceil(start2);                                  % head index number
                h_weigh2            =ceil(start2)-start2;                           % voxelnr* weighing
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),h_index2,t_index3);    % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*h_weigh2*t_weigh3) ;
                        weigh           =weigh      + (h_weigh1*h_weigh2*t_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,h_index2,t_index3);        % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*h_weigh2*t_weigh3) ;
                        weigh           =weigh      + (b_weigh1*h_weigh2*t_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),h_index2,t_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*h_weigh2*t_weigh3);
                        weigh           =weigh      + (t_weigh1*h_weigh2*t_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end
 
 
                if  ceil(start2)<floor(end2)        % if there is a body2
                for b_index2    =ceil(start2)+1:floor(end2)
                b_weigh2        =1;
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),b_index2,t_index3);             % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*b_weigh2*t_weigh3) ;
                        weigh           =weigh      + (h_weigh1*b_weigh2*t_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,b_index2,t_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*b_weigh2*t_weigh3) ;
                        weigh           =weigh      + (b_weigh1*b_weigh2*t_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),b_index2,t_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*b_weigh2*t_weigh3);
                        weigh           =weigh      + (t_weigh1*b_weigh2*t_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end;end
 
                if  floor(end2)<ceil(end2)          % if there is a tail2
                t_index2        =ceil(end2);
                t_weigh2        =(end2-floor(end2));
 
                        if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),t_index2,t_index3);             % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*t_weigh2*t_weigh3) ;
                        weigh           =weigh      + (h_weigh1*t_weigh2*t_weigh3);end       % add weight to totalweigh
 
                        if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,t_index2,t_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*t_weigh2*t_weigh3) ;
                        weigh           =weigh      + (b_weigh1*t_weigh2*t_weigh3);end;end
 
                        if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),t_index2,t_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*t_weigh2*t_weigh3);
                        weigh           =weigh      + (t_weigh1*t_weigh2*t_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;end;end
                clear h_index2 h_weigh2 b_index2 b_weigh2 t_index2 t_weigh2;end            
            
            
            output_res(i,j,k)      =sum/weigh;                
            clear h_index3 h_weigh3 b_index3 b_weigh3 t_index3 t_weigh3 sum weigh                   
        end
    end
end
%% or
elseif newsize(1)==size(im_input,1) && newsize(2)==size(im_input,2) && newsize(3)~=size(im_input,3)
%% 1-dimensional change third dimension
for k=1:size(new_res,3)
 
    start3              =(k-1)* (size(old_res,3) / size(new_res,3));
    end3                =k    * (size(old_res,3) / size(new_res,3));            
 
    sum                 =zeros(size(new_res,1),size(new_res,2));
    weight              =0;
 
    if  floor(start3)<ceil(start3)                              % if there is a head3 (heading partial voxel)
        h_value3        =old_res(:,:,ceil(start3));              % head values
        h_weight3       =( ceil(start3)-start3 );               % head weight (partial voxel)
        sum             =sum        + (h_value3.*h_weight3) ;
        weight          =weight      + h_weight3;end            % add weight to totalweight
 
    if  ceil(start3)<floor(end3)                                % if there is a body3
        for b_index3        =ceil(start3)+1:floor(end3)
            b_value3        =old_res(:,:,b_index3);              % body values
            b_weight3       =1;                                 % body weight (or count, are complete voxels)
            sum             =sum        + (b_value3.*b_weight3) ;
            weight          =weight     + (b_weight3);end;end
 
    if  floor(end3)<ceil(end3)                                  % if there is a tail3
        t_value3        =old_res(:,:,ceil(end3));                % tail index value
        t_weight3       =end3-floor(end3);                      % tail weight (partial voxel)
        sum             =sum        + (t_value3.*t_weight3);
        weight          =weight     + (t_weight3);end
 
        output_res(:,:,k)      =sum./weight;
        clear h_value3 h_weight3 b_index3 b_value3 b_weight3 t_value3 t_weight3 sum weight
end


 elseif newsize(1)==size(im_input,1) && newsize(2)~=size(im_input,2) && newsize(3)==size(im_input,3)
%% 1-dimensional change second dimension
for k=1:size(new_res,2)
 
    start2              =(k-1)* (size(old_res,2) / size(new_res,2));
    end2                =k    * (size(old_res,2) / size(new_res,2));            
 
    sum                 =zeros(size(new_res,1),1,size(new_res,3));
    weight              =0;
 
    if  floor(start2)<ceil(start2)                              % if there is a head2 (heading partial voxel)
        h_value2        =old_res(:,ceil(start2),:);             % head values
        h_weight2       =( ceil(start2)-start2 );               % head weight (partial voxel)
        sum             =sum        + (h_value2.*h_weight2) ;
        weight          =weight      + h_weight2;end            % add weight to totalweight
 
    if  ceil(start2)<floor(end2)                                % if there is a body2
        for b_index2        =ceil(start2)+1:floor(end2)
            b_value2        =old_res(:,b_index2,:);             % body values
            b_weight2       =1;                                 % body weight (or count, are complete voxels)
            sum             =sum        + (b_value2.*b_weight2) ;
            weight          =weight     + (b_weight2);end;end
 
    if  floor(end2)<ceil(end2)                                  % if there is a tail2
        t_value2        =old_res(:,ceil(end2),:);               % tail index value
        t_weight2       =end2-floor(end2);                      % tail weight (partial voxel)
        sum             =sum        + (t_value2.*t_weight2);
        weight          =weight     + (t_weight2);end
 
        output_res(:,k,:)      =sum./weight;
        clear h_value2 h_weight2 b_index2 b_value2 b_weight2 t_value2 t_weight2 sum weight
end

 elseif newsize(1)~=size(im_input,1) && newsize(2)==size(im_input,2) && newsize(3)==size(im_input,3)
%% 1-dimensional change first dimension
for k=1:size(new_res,1)
 
    start1              =(k-1)* (size(old_res,1) / size(new_res,1));
    end1                =k    * (size(old_res,1) / size(new_res,1));            
 
    sum                 =zeros(1,size(new_res,2),size(new_res,3));
    weight              =0;
 
    if  floor(start1)<ceil(start1)                              % if there is a head1 (heading partial voxel)
        h_value1        =old_res(ceil(start1),:,:);             % head values
        h_weight1       =( ceil(start1)-start1 );               % head weight (partial voxel)
        sum             =sum        + (h_value1.*h_weight1) ;
        weight          =weight      + h_weight1;end            % add weight to totalweight
 
    if  ceil(start1)<floor(end1)                                % if there is a body1
        for b_index1        =ceil(start1)+1:floor(end1)
            b_value1        =old_res(b_index1,:,:);             % body values
            b_weight1       =1;                                 % body weight (or count, are complete voxels)
            sum             =sum        + (b_value1.*b_weight1) ;
            weight          =weight     + (b_weight1);end;end
 
    if  floor(end1)<ceil(end1)                                  % if there is a tail1
        t_value1        =old_res(ceil(end1),:,:);               % tail index value
        t_weight1       =end1-floor(end1);                      % tail weight (partial voxel)
        sum             =sum        + (t_value1.*t_weight1);
        weight          =weight     + (t_weight1);end
 
        output_res(k,:,:)      =sum./weight;
        clear h_value1 h_weight1 b_index1 b_value1 b_weight1 t_value1 t_weight1 sum weight
end


else output_res=im_input;

end
    % Close waitbar after processing
    if showWaitbar
        try
            close(h);
        catch
            warning('Waitbar manually closed...');
        end
    end
end


