function output_res = xASL_im_ResampleLinearFair(im_input, newsize)
%xASL_im_ResampleLinearFair Downsample or upsample a 1D/2D/3D image
%
% FORMAT:       imOutput = xASL_im_ResampleLinearFair(imInput, newSize)
%
% INPUT:        imInput     - Image matrix (REQUIRED, DOUBLE, SINGLE or INT)
%               newSize     - Size of ouput image (REQUIRED, INTEGER ARRAY)
%
% OUTPUT:       imOutput    - Resampled image matrix (SINGLE)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Downsample or upsample an image from its old to a new resolution.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      image = xASL_io_Nifti2Im('.\M0.nii.gz');
%               output = xASL_im_ResampleLinearFair(image, [2,2,2]);
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Calculation
    old_res = single(im_input);
    new_res = zeros(newsize(1), newsize(2), newsize(3));
    

    %% Get correct case
    if newsize(1)==size(im_input,1) && newsize(2)==size(im_input,2) && newsize(3)==size(im_input,3)
        % Return output as input
        output_res = im_input;
    elseif newsize(1)==size(im_input,1) && newsize(2)==size(im_input,2) && newsize(3)~=size(im_input,3)
        % 1-dimensional change third dimension
        output_res = xASL_im_ResampleLinearFair_1D(old_res,new_res,3);
    elseif newsize(1)==size(im_input,1) && newsize(2)~=size(im_input,2) && newsize(3)==size(im_input,3)
        % 1-dimensional change second dimension
        output_res = xASL_im_ResampleLinearFair_1D(old_res,new_res,2);
    elseif newsize(1)~=size(im_input,1) && newsize(2)==size(im_input,2) && newsize(3)==size(im_input,3)
        % 1-dimensional change first dimension
        output_res = xASL_im_ResampleLinearFair_1D(old_res,new_res,1);
    else
        % 3-dimensional change
        output_res = xASL_im_ResampleLinearFair_3D(old_res,new_res);
    end


end


%% 3-dimensional change
function output_res = xASL_im_ResampleLinearFair_3D(old_res,new_res)

    % Dummy initialization
    output_res = 0;

    % Feedback
    fprintf('Resample:    ');
    
    % Iterate over image
    for k=1:size(new_res,3)
        % Prepare new_res blocks
        for j=1:size(new_res,2)
            % Prepare new_res blocks
            for it=1:size(new_res,1)
                output_res = xASL_im_ResampleLinearFair_Loop3D(it,j,k,old_res,new_res,output_res);
            end
        end
    end
    fprintf('\n');

end


%% Resample 3D
function output_res = xASL_im_ResampleLinearFair_Loop3D(it,j,k,old_res,new_res,output_res)

    % Display progress
    percentage =((k-1)*size(new_res,2)*size(new_res,1)+(j-1)*size(new_res,1)+it) / ...
        (size(new_res,3) * size(new_res,2) * size(new_res,1));
    xASL_TrackProgress(percentage*100);

    % Determine start and end for each dim
    start1  = (it-1)* (size(old_res,1) / size(new_res,1));
    end1    = it    * (size(old_res,1) / size(new_res,1));
    start2  = (j-1) * (size(old_res,2) / size(new_res,2));
    end2    = j     * (size(old_res,2) / size(new_res,2));
    start3  = (k-1) * (size(old_res,3) / size(new_res,3));
    end3    = k     * (size(old_res,3) / size(new_res,3));

    % Initialization
    sum    = 0;
    weigh  = 0;

    % Head-3
    [sum,weigh] = xASL_im_ResampleLinearFair_Head3(sum,weigh,start1,end1,start2,end2,start3,end3,old_res,new_res,output_res);
    
    % Body-3
    [sum,weigh] = xASL_im_ResampleLinearFair_Body3(sum,weigh,start1,end1,start2,end2,start3,end3,old_res,new_res,output_res);
    
    % Tail-3
    [sum,weigh] = xASL_im_ResampleLinearFair_Tail3(sum,weigh,start1,end1,start2,end2,start3,end3,old_res,new_res,output_res);
    
    % Store result voxel
    output_res(it,j,k) = sum/weigh;
    
    % Free up space
    clear sum weigh

end


%% Head-3
function [sum,weigh] = xASL_im_ResampleLinearFair_Head3(sum,weigh,start1,end1,start2,end2,start3,end3,old_res,new_res,output_res)

    if  floor(start3)<ceil(start3)                                          % if there is a head3 (heading partial voxel)
        h_index3            =ceil(start3);                                  % head index number
        h_weigh3            =ceil(start3)-start3;                           % voxelnr* weighing


        if  floor(start2)<ceil(start2)                                      % if there is a head2 (heading partial voxel)
            h_index2            =ceil(start2);                              % head index number
            h_weigh2            =ceil(start2)-start2;                       % voxelnr* weighing

            if  floor(start1)<ceil(start1)                                  % if there is a head1
                h_value1        =old_res(ceil(start1),h_index2,h_index3);   % head index value
                h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                sum             =sum        + (h_value1*h_weigh1*h_weigh2*h_weigh3); % add weight to totalweigh
                weigh           =weigh      + (h_weigh1*h_weigh2*h_weigh3);
            end

            if  ceil(start1)<floor(end1)                                    % if there is a body1
                for b_index1=ceil(start1)+1:floor(end1)
                    b_value1        =old_res(b_index1,h_index2,h_index3);   % body index value
                    b_weigh1        =1;                                     % body weigh (or count, are complete voxels)
                    sum             =sum        + (b_value1*b_weigh1*h_weigh2*h_weigh3) ;
                    weigh           =weigh      + (b_weigh1*h_weigh2*h_weigh3);
                end
            end

            if  floor(end1)<ceil(end1)                                      % if there is a tail1
                t_value1        =old_res(ceil(end1),h_index2,h_index3);     % tail index value
                t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                sum             =sum        + (t_value1*t_weigh1*h_weigh2*h_weigh3);
                weigh           =weigh      + (t_weigh1*h_weigh2*h_weigh3);
                clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
            end
        end


        if  ceil(start2)<floor(end2)                                        % if there is a body2
            for b_index2    =ceil(start2)+1:floor(end2)
                b_weigh2        =1;

                if  floor(start1)<ceil(start1)                              % if there is a head1
                    h_value1        =old_res(ceil(start1),b_index2,h_index3); % head index value
                    h_weigh1        =( ceil(start1)-start1 );               % head weigh (partial voxel)
                    sum             =sum        + (h_value1*h_weigh1*b_weigh2*h_weigh3) ;
                    weigh           =weigh      + (h_weigh1*b_weigh2*h_weigh3); % add weight to totalweigh
                end

                if  ceil(start1)<floor(end1)                                % if there is a body1
                    for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,b_index2,h_index3); % body index value
                        b_weigh1        =1;                                 % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*b_weigh2*h_weigh3) ;
                        weigh           =weigh      + (b_weigh1*b_weigh2*h_weigh3);
                    end
                end

                if  floor(end1)<ceil(end1)                                  % if there is a tail1
                    t_value1        =old_res(ceil(end1),b_index2,h_index3); % tail index value
                    t_weigh1        =end1-floor(end1);                      % tail weigh (partial voxel)
                    sum             =sum        + (t_value1*t_weigh1*b_weigh2*h_weigh3);
                    weigh           =weigh      + (t_weigh1*b_weigh2*h_weigh3);
                    clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
                end
            end
        end

        if  floor(end2)<ceil(end2)          % if there is a tail2
            t_index2        =ceil(end2);
            t_weigh2        =(end2-floor(end2));

            if  floor(start1)<ceil(start1)                              % if there is a head1
                h_value1        =old_res(ceil(start1),t_index2,h_index3);             % head index value
                h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                sum             =sum        + (h_value1*h_weigh1*t_weigh2*h_weigh3) ;
                weigh           =weigh      + (h_weigh1*t_weigh2*h_weigh3);       % add weight to totalweigh
            end

            if  ceil(start1)<floor(end1)                                % if there is a body1
                for b_index1=ceil(start1)+1:floor(end1)
                    b_value1        =old_res(b_index1,t_index2,h_index3);                 % body index value
                    b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                    sum             =sum        + (b_value1*b_weigh1*t_weigh2*h_weigh3) ;
                    weigh           =weigh      + (b_weigh1*t_weigh2*h_weigh3);
                end
            end

            if  floor(end1)<ceil(end1)                                  % if there is a tail1
                t_value1        =old_res(ceil(end1),t_index2,h_index3);               % tail index value
                t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                sum             =sum        + (t_value1*t_weigh1*t_weigh2*h_weigh3);
                weigh           =weigh      + (t_weigh1*t_weigh2*h_weigh3);
                clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
            end
        end
        clear h_index2 h_weigh2 b_index2 b_weigh2 t_index2 t_weigh2;
    end

end



%% Body-3
function [sum,weigh] = xASL_im_ResampleLinearFair_Body3(sum,weigh,start1,end1,start2,end2,start3,end3,old_res,new_res,output_res)

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
                    weigh           =weigh      + (h_weigh1*h_weigh2*b_weigh3);       % add weight to totalweigh
                end

                if  ceil(start1)<floor(end1)                                % if there is a body1
                    for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,h_index2,b_index3);        % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*h_weigh2*b_weigh3) ;
                        weigh           =weigh      + (b_weigh1*h_weigh2*b_weigh3);
                    end
                end

                if  floor(end1)<ceil(end1)                                  % if there is a tail1
                    t_value1        =old_res(ceil(end1),h_index2,b_index3);               % tail index value
                    t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                    sum             =sum        + (t_value1*t_weigh1*h_weigh2*b_weigh3);
                    weigh           =weigh      + (t_weigh1*h_weigh2*b_weigh3);
                    clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
                end
            end


            if  ceil(start2)<floor(end2)        % if there is a body2
                for b_index2    =ceil(start2)+1:floor(end2)
                    b_weigh2        =1;

                    if  floor(start1)<ceil(start1)                              % if there is a head1
                        h_value1        =old_res(ceil(start1),b_index2,b_index3);             % head index value
                        h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                        sum             =sum        + (h_value1*h_weigh1*b_weigh2*b_weigh3) ;
                        weigh           =weigh      + (h_weigh1*b_weigh2*b_weigh3);       % add weight to totalweigh
                    end

                    if  ceil(start1)<floor(end1)                                % if there is a body1
                        for b_index1=ceil(start1)+1:floor(end1)
                            b_value1        =old_res(b_index1,b_index2,b_index3);                 % body index value
                            b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                            sum             =sum        + (b_value1*b_weigh1*b_weigh2*b_weigh3) ;
                            weigh           =weigh      + (b_weigh1*b_weigh2*b_weigh3);
                        end
                    end

                    if  floor(end1)<ceil(end1)                                  % if there is a tail1
                        t_value1        =old_res(ceil(end1),b_index2,b_index3);               % tail index value
                        t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                        sum             =sum        + (t_value1*t_weigh1*b_weigh2*b_weigh3);
                        weigh           =weigh      + (t_weigh1*b_weigh2*b_weigh3);
                        clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
                    end
                end
            end

            if  floor(end2)<ceil(end2)          % if there is a tail2
                t_index2        =ceil(end2);
                t_weigh2        =(end2-floor(end2));

                if  floor(start1)<ceil(start1)                              % if there is a head1
                    h_value1        =old_res(ceil(start1),t_index2,b_index3);             % head index value
                    h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                    sum             =sum        + (h_value1*h_weigh1*t_weigh2*b_weigh3) ;
                    weigh           =weigh      + (h_weigh1*t_weigh2*b_weigh3);       % add weight to totalweigh
                end

                if  ceil(start1)<floor(end1)                                % if there is a body1
                    for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,t_index2,b_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*t_weigh2*b_weigh3) ;
                        weigh           =weigh      + (b_weigh1*t_weigh2*b_weigh3);
                    end
                end

                if  floor(end1)<ceil(end1)                                  % if there is a tail1
                    t_value1        =old_res(ceil(end1),t_index2,b_index3);               % tail index value
                    t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                    sum             =sum        + (t_value1*t_weigh1*t_weigh2*b_weigh3);
                    weigh           =weigh      + (t_weigh1*t_weigh2*b_weigh3);
                    clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
                end
            end
            clear h_index2 h_weigh2 b_index2 b_weigh2 t_index2 t_weigh2;
        end
    end

end


%% Tail-3
function [sum,weigh] = xASL_im_ResampleLinearFair_Tail3(sum,weigh,start1,end1,start2,end2,start3,end3,old_res,new_res,output_res)

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
                weigh           =weigh      + (h_weigh1*h_weigh2*t_weigh3);       % add weight to totalweigh
            end

            if  ceil(start1)<floor(end1)                                % if there is a body1
                for b_index1=ceil(start1)+1:floor(end1)
                    b_value1        =old_res(b_index1,h_index2,t_index3);        % body index value
                    b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                    sum             =sum        + (b_value1*b_weigh1*h_weigh2*t_weigh3) ;
                    weigh           =weigh      + (b_weigh1*h_weigh2*t_weigh3);
                end
            end

            if  floor(end1)<ceil(end1)                                  % if there is a tail1
                t_value1        =old_res(ceil(end1),h_index2,t_index3);               % tail index value
                t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                sum             =sum        + (t_value1*t_weigh1*h_weigh2*t_weigh3);
                weigh           =weigh      + (t_weigh1*h_weigh2*t_weigh3);
                clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
            end
        end


        if  ceil(start2)<floor(end2)        % if there is a body2
            for b_index2    =ceil(start2)+1:floor(end2)
                b_weigh2        =1;

                if  floor(start1)<ceil(start1)                              % if there is a head1
                    h_value1        =old_res(ceil(start1),b_index2,t_index3);             % head index value
                    h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                    sum             =sum        + (h_value1*h_weigh1*b_weigh2*t_weigh3) ;
                    weigh           =weigh      + (h_weigh1*b_weigh2*t_weigh3);       % add weight to totalweigh
                end

                if  ceil(start1)<floor(end1)                                % if there is a body1
                    for b_index1=ceil(start1)+1:floor(end1)
                        b_value1        =old_res(b_index1,b_index2,t_index3);                 % body index value
                        b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                        sum             =sum        + (b_value1*b_weigh1*b_weigh2*t_weigh3) ;
                        weigh           =weigh      + (b_weigh1*b_weigh2*t_weigh3);
                    end
                end

                if  floor(end1)<ceil(end1)                                  % if there is a tail1
                    t_value1        =old_res(ceil(end1),b_index2,t_index3);               % tail index value
                    t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                    sum             =sum        + (t_value1*t_weigh1*b_weigh2*t_weigh3);
                    weigh           =weigh      + (t_weigh1*b_weigh2*t_weigh3);
                    clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
                end
            end
        end

        if  floor(end2)<ceil(end2)          % if there is a tail2
            t_index2        =ceil(end2);
            t_weigh2        =(end2-floor(end2));

            if  floor(start1)<ceil(start1)                              % if there is a head1
                h_value1        =old_res(ceil(start1),t_index2,t_index3);             % head index value
                h_weigh1        =( ceil(start1)-start1 );                   % head weigh (partial voxel)
                sum             =sum        + (h_value1*h_weigh1*t_weigh2*t_weigh3) ;
                weigh           =weigh      + (h_weigh1*t_weigh2*t_weigh3);       % add weight to totalweigh
            end

            if  ceil(start1)<floor(end1)                                % if there is a body1
                for b_index1=ceil(start1)+1:floor(end1)
                    b_value1        =old_res(b_index1,t_index2,t_index3);                 % body index value
                    b_weigh1        =1;                                         % body weigh (or count, are complete voxels)
                    sum             =sum        + (b_value1*b_weigh1*t_weigh2*t_weigh3) ;
                    weigh           =weigh      + (b_weigh1*t_weigh2*t_weigh3);
                end
            end

            if  floor(end1)<ceil(end1)                                  % if there is a tail1
                t_value1        =old_res(ceil(end1),t_index2,t_index3);               % tail index value
                t_weigh1        =end1-floor(end1);                          % tail weigh (partial voxel)
                sum             =sum        + (t_value1*t_weigh1*t_weigh2*t_weigh3);
                weigh           =weigh      + (t_weigh1*t_weigh2*t_weigh3);
                clear h_value1 h_weigh1 b_index1 b_value1 b_weigh1 t_value1 t_weigh1;
            end
        end
        clear h_index2 h_weigh2 b_index2 b_weigh2 t_index2 t_weigh2
    end

end


%% 1-dimensional change n-th dimension
function output_res = xASL_im_ResampleLinearFair_1D(old_res,new_res,dim)

    for k=1:size(new_res,dim)
        % Initialize weights
        weight = 0;
        
        % Define start and end of current dim
        startD = (k-1)* (size(old_res,dim) / size(new_res,dim));
        endD = k * (size(old_res,dim) / size(new_res,dim));

        % Get sum of current dim
        sum = xASL_im_Resample1D_sum(new_res,dim);

        % If there is a head (heading partial voxel)
        [sum,weight] = xASL_im_Resample1D_head(sum,weight,startD,old_res,dim);

        % If there is a body
        [sum,weight] = xASL_im_Resample1D_body(sum,weight,startD,endD,old_res,dim);

        % If there is a tail
        [sum,weight] = xASL_im_Resample1D_tail(sum,weight,endD,old_res,dim);

        % Save resampled voxel
        output_res(:,:,k) = sum./weight;
        
        % Free up space
        clear sum weight
    end

end


%% Determine sum of current dim
function sum = xASL_im_Resample1D_sum(new_res,dim)

    switch dim
        case 3
            sum = zeros(size(new_res,1),size(new_res,2));
        case 2
            sum = zeros(size(new_res,1),1,size(new_res,3));
        case 1
            sum = zeros(1,size(new_res,2),size(new_res,3));
    end

end


%% Determine h value
function h_value = xASL_im_Resample1D_h_value(old_res,startD,dim)

    switch dim
        case 3
            h_value = old_res(:,:,ceil(startD));
        case 2
            h_value = old_res(:,ceil(startD),:);
        case 1
            h_value = old_res(ceil(startD),:,:);
    end

end


%% Determine b value
function b_value = xASL_im_Resample1D_b_value(old_res,b_index,dim)

    switch dim
        case 3
            b_value = old_res(:,:,b_index);
        case 2
            b_value = old_res(:,b_index,:);
        case 1
            b_value = old_res(b_index,:,:);
    end

end


%% Determine t value
function t_value = xASL_im_Resample1D_t_value(old_res,endD,dim)

    switch dim
        case 3
            t_value  = old_res(:,:,ceil(endD));
        case 2
            t_value  = old_res(:,ceil(endD),:);
        case 1
            t_value  = old_res(ceil(endD),:,:);
    end

end


%% Head
function [sum,weight] = xASL_im_Resample1D_head(sum,weight,startD,old_res,dim)

    if floor(startD)<ceil(startD)
        % head values
        h_value = xASL_im_Resample1D_h_value(old_res,startD,dim);

        % head weight (partial voxel)
        h_weight = (ceil(startD)-startD);
        sum = sum + (h_value.*h_weight) ;

        % add weight to totalweight
        weight = weight + h_weight;
    end

end


%% Body
function [sum,weight] = xASL_im_Resample1D_body(sum,weight,startD,endD,old_res,dim)

    if ceil(startD)<floor(endD)
        for b_index = ceil(startD)+1:floor(endD)
            % body values
            b_value = xASL_im_Resample1D_b_value(old_res,b_index,dim);

            % body weight (or count, are complete voxels)
            b_weight = 1;
            sum = sum + (b_value.*b_weight) ;
            weight = weight + (b_weight);
        end
    end

end


%% Tail
function [sum,weight] = xASL_im_Resample1D_tail(sum,weight,endD,old_res,dim)

    if  floor(endD)<ceil(endD)
        % tail index value
        t_value = xASL_im_Resample1D_t_value(old_res,endD,dim);
        % tail weight (partial voxel)
        t_weight = endD-floor(endD);
        sum = sum + (t_value.*t_weight);
        weight = weight + (t_weight);
    end

end

