function [nameGIF, r] = ps_LST_create_gif(varargin)
%ps_LST_create_gif   Creates PNG images for given images.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
%
%   Sorry, there is no further documentation at this moment.
%

% Check input
nonewname = 0;
r = 0;
or = [1 2 3];
if nargin == 0
    V = spm_select(1, 'image', 'Select image.');
    Vover = spm_select(Inf, 'image', 'Select overlay images (optional).');
end
if nargin > 0
    if ischar(varargin{1})
        V = varargin{1};
    else            
        fprintf('Input for V must be a character, like from spm_select.\n');            
        return;
    end
    if nargin > 1
        if ischar(varargin{2})
            Vover = varargin{2};
        else                
            fprintf('Input for Vover must be a character, like from spm_select.\n');                
            return;
        end
    end
    if nargin > 2
        nameGIF = varargin{3};
        nonewname = 1;                
        %V = varargin{1};        
        %Vover = varargin{2};
        %img_input = 1;
    end
    if nargin > 3
        r = varargin{4};        
    end
    if nargin > 4
        or = varargin{5};
        fl = or(end);
        or = or(1:3);
    end
end

% create file name
orpth = cd;
if numel(V) == 1
    [pth, nam, ~] = fileparts(V);
else
    [pth, nam, ~] = fileparts(V(1,:));
end
if ~nonewname    
    nameGIF = ['LST_', nam, '_', ps_create_timestamp];    
end
namePNG = nameGIF;
nameGIF = [nameGIF, '.gif'];
folder = namePNG;
cd(pth)
if ~nonewname
    mkdir(folder);
end

overlay = ~isempty(Vover);
V = spm_vol(V);
Vover = spm_vol(Vover);


% Remove ugly border around the image
%iptsetpref('ImshowBorder','tight');
fid = figure('visible', 'off');

if overlay
    
    % How many overlays?
    n = numel(Vover);
    m = numel(V);
    if m == 1
        % Do all images have the same dimension?
        dims = zeros(n+1, 3);    
        dims(end,:) = V.dim;
        lcl = zeros(n, 1);
        for i = 1:n        
            dims(i,:) = Vover(i).dim;
            lcl(i) = ~isempty(regexp(Vover(i).fname, 'LCL', 'once'));
        end

        if any(sum(dims) ~= (dims(1,:) .* (n+1)))
            fprintf('Images must all have the same dimensions!\n')
            close(fid)
            return;
        end

        % Load images            
        over = zeros(dims(1,1), dims(1,2), dims(1,3), n);
        for i = 1:n
            tmp = spm_read_vols(Vover(i));
            % Flip image
            if any(or ~= [1 2 3])
                tmp = permute(tmp, or);
                if i == 1
                    over = zeros(size(tmp,1), size(tmp,2), size(tmp,3), n);
                end
            end
            if fl
                if exist('flip', 'builtin')
                    tmp = flip(tmp, 2);
                else
                    tmp = flipdim(tmp, 2);
                end
            end
            over(:,:,:,i) = tmp;
        end
        
        img = spm_read_vols(V);
        %if any(or ~= [1 3 2])
        if any(or ~= [1 2 3])
            img = permute(img, or);            
        end
        if fl
            if exist('flip', 'builtin')
                img = flip(img, 2);
            else
                img = flipdim(img, 2);
            end
        end
        img(isnan(img)) = 0;
        img(:) = (img(:) ./ ps_quantile(img(:), .9999));
        img(img >= 1) = 0.9999;

        % Where to start?
        if numel(r) == 1
            c_tmp = indx2coord(find(img > 0), size(img, 1), size(img, 2));
            r = [min(c_tmp(:,3)), max(c_tmp(:,3))];
        end

        counter = 0;
        for z = r(1):r(2)
            z_tmp = r(2) - counter;
            counter = counter + 1;            
            if exist('flip', 'builtin')
                img_tmp = flip(img(:,:,z_tmp)');
                over_tmp = 0 .* flip(over(:,:,z_tmp,1)');
            else
                img_tmp = flipdim(img(:,:,z_tmp)', 1);
                over_tmp = 0 .* flipdim(over(:,:,z_tmp,1)', 1);
            end
            for i = 1:n
                if exist('flip', 'builtin')
                    over_tmp = [over_tmp, flip(over(:,:,z_tmp,i)')];
                    img_tmp = [img_tmp, flip(img(:,:,z_tmp)')];
                else
                    over_tmp = [over_tmp, flipdim(over(:,:,z_tmp,i)', 1)];
                    img_tmp = [img_tmp, flipdim(img(:,:,z_tmp)', 1)];
                end                
            end                

            over_tmp_rgb = zeros(size(over_tmp, 1), size(over_tmp, 2), 3);
            if any(lcl > 0)

                over_tmp2 = over_tmp .* 0;
                over_tmp2(over_tmp == 1) = 0;
                over_tmp2(over_tmp == 2) = 144;
                over_tmp2(over_tmp == 3) = 208;
                over_tmp_rgb(:,:,1) = over_tmp2;

                over_tmp2 = over_tmp .* 0;
                over_tmp2(over_tmp == 1) = 255;
                over_tmp2(over_tmp == 2) = 144;
                over_tmp2(over_tmp == 3) = 0;
                over_tmp_rgb(:,:,2) = over_tmp2;

                over_tmp2 = over_tmp .* 0;
                over_tmp2(over_tmp == 1) = 102;
                over_tmp2(over_tmp == 2) = 144;
                over_tmp2(over_tmp == 3) = 0;
                over_tmp_rgb(:,:,3) = over_tmp2;    
                
                col = zeros(3,3);
                col(:,1) = [42, 253, 110];
                col(:,2) = [252, 189, 65];
                col(:,3) = [192, 7, 19];
                over_tmp_rgb = ind2rgb(over_tmp, col);
                                
            else
                % define jet-colors                
                over_tmp(over_tmp > .9) = .9;
                %[X, ~] = gray2ind(over_tmp, 256);
                X = round(over_tmp .* 255);
                over_tmp_rgb = ind2rgb(X, hot(256));                
            end

            % bring figure to front (pretend that figure is displayed in
            % other figures that may have been clicked by the user)
            figure(fid);
            % make figure invisible
            set(fid, 'visible', 'off')            
            %imshow(over_tmp_rgb);
            imagesc(over_tmp_rgb);

            axis off
            over_tmp(over_tmp > 0) = 1;
            over_tmp = over_tmp == 0;            
            img_tmp_gr = img_tmp;

            hold on                   
            %h = imshow(img_tmp_gr);
            h = imagesc(img_tmp_gr); colormap(gray);
            hold off
            set(h, 'AlphaData', over_tmp);
            
            res = 50; % pixels per inch
            set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 size(img_tmp_gr, 2) size(img_tmp_gr, 1)]/res);
            %saveas(fid, fullfile(namePNG, ['overlay_', num2str(z_tmp), '.png']), 'png') % CHANGE!
            saveas(gcf, fullfile(namePNG, ['overlay_', num2str(z_tmp), '.png']), 'png')

        end
        
    else
        % Do all images have the same dimension?
        dims = zeros(3, 3);    
        dims(end,:) = Vover.dim;
        lcl = zeros(n, 1);
        for i = 1:2
            dims(i,:) = V(i).dim;            
        end
        
        % Load images    
        over = spm_read_vols(Vover);
        if any(or ~= [1 2 3])
            over = permute(over, or);            
            dims = [size(over); size(over)];
        end
        if fl
            if exist('flip', 'builtin')
                over = flip(over, 2);
            else
                over = flipdim(over, 2);
            end
        end
        img = zeros(dims(1,1), dims(1,2), dims(1,3), m);
        for i = 1:m
            tmp = spm_read_vols(V(i));
            if any(or ~= [1 2 3])
                tmp = permute(tmp, or);                
            end
            if fl
                if exist('flip', 'builtin')
                    tmp = flip(tmp, 2);
                else
                    tmp = flipdim(tmp, 2);
                end
            end
            tmp(tmp < 0) = 0;
            tmp(isnan(tmp)) = 0;
            tmp(:) = (tmp(:) ./ ps_quantile(tmp(:), .9999));
            img(:,:,:,i) = tmp;            
        end                        
        img(img > 1) = 1;

        % Where to start?
        if numel(r) == 1
            c_tmp = indx2coord(find(img(:,:,:,1) > 0), size(img, 1), size(img, 2));
            r = [min(c_tmp(:,3)), max(c_tmp(:,3))];
        end

        counter = 0;
        for z = r(1):r(2)
            z_tmp = r(2) - counter;
            counter = counter + 1;
            if exist('flip', 'builtin')
                img_tmp = [flip(img(:,:,z_tmp,1)'), flip(img(:,:,z_tmp,2)'), flip(img(:,:,z_tmp,1)')];
                over_tmp = [0 .* flip(over(:,:,z_tmp)'), 0 .* flip(over(:,:,z_tmp)'), 1 .* flip(over(:,:,z_tmp)')]; 
            else
                img_tmp = [flipdim(img(:,:,z_tmp,1)', 1), flipdim(img(:,:,z_tmp,2)', 1), flipdim(img(:,:,z_tmp,1)', 1)];
                over_tmp = [0 .* flipdim(over(:,:,z_tmp)', 1), 0 .* flipdim(over(:,:,z_tmp)', 1), 1 .* flipdim(over(:,:,z_tmp)', 1)]; 
            end
            
            col = zeros(4,3);
            col(1,:) = [0, 0, 0];
            col(2,:) = [42, 253, 110];
            col(3,:) = [252, 189, 65];
            col(4,:) = [192, 7, 19];
            %over_tmp_rgb = ind2rgb(gray2ind(over_tmp ./ 3, 4), col./255);
            over_tmp_rgb = ind2rgb(uint8(over_tmp), col./255);

            % bring figure to front (pretend that figure is displayed in
            % other figures that may have been clicked by the user)
            figure(fid);
            % make figure invisible
            set(fid, 'visible', 'off')
            % plot image
            %imshow(over_tmp_rgb);
            imagesc(over_tmp_rgb);

            axis off
            over_tmp(over_tmp > 0) = 1;
            over_tmp = over_tmp == 0;
                        
            % Add overlay
            img_tmp_gr = img_tmp;
            hold on       
            %h = imshow(img_tmp_gr);
            h = imagesc(img_tmp_gr); colormap(gray);
            hold off            
            % set transparency
            set(h, 'AlphaData', over_tmp);                        
                        
            res = 50; % pixels per inch
            set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 size(img_tmp_gr, 2) size(img_tmp_gr, 1)]/res);
            if nonewname
                %saveas(h, [namePNG, '_', num2str(z_tmp), '.png'], 'png') % CHANGE!               
                saveas(gcf, [namePNG, '_', num2str(z_tmp), '.png'], 'png')
            else
                %saveas(h, ['overlay_', num2str(z_tmp), '.png'], 'png') % CHANGE!
                saveas(gcf, ['overlay_', num2str(z_tmp), '.png'], 'png')
            end

        end
    end
    
else    
    % Load images    
    img = spm_read_vols(V);
    img(isnan(img)) = 0;
    img(:) = (img(:) ./ ps_quantile(img(:), .9999));
    img(img > 1) = 1;
        
    if any(or ~= [1 2 3])
        img = permute(img, [or, 4]);
    end
    if fl
        if exist('flip', 'builtin')
            img = flip(img, 2);
        else
            img = flipdim(img, 2);
        end
    end
    
    % More than one image?
    n = numel(V);
    
    % Where to start?
    if numel(r) == 1
        if n > 1
            c_tmp = indx2coord(find(img(:,:,:,1) > 0), size(img(:,:,:,1), 1), size(img(:,:,:,1), 2));
            r = [min(c_tmp(:,3)), max(c_tmp(:,3))];
        else
            c_tmp = indx2coord(find(img > 0), size(img, 1), size(img, 2));
            r = [min(c_tmp(:,3)), max(c_tmp(:,3))];
        end
    end
    
    counter = 0;
    for z = r(1):r(2)
        z_tmp = r(2) - counter;
        counter = counter + 1;
        if n > 1
            if exist('flip', 'builtin')
                img_tmp = flip(img(:,:,z_tmp,1)');
                for j = 2:n
                    img_tmp = [img_tmp, flip(img(:,:,z_tmp,j)')];
                end
            else
                img_tmp = flipdim(img(:,:,z_tmp,1)', 1);
                for j = 2:n
                    img_tmp = [img_tmp, flipdim(img(:,:,z_tmp,j)', 1)];
                end
            end
        else
            if exist('flip', 'builtin')
                img_tmp = flip(img(:,:,z_tmp)');
            else
                img_tmp = flipdim(img(:,:,z_tmp)', 1);
            end
        end
        
        img_tmp_gr = img_tmp; 
        % bring figure to front (pretend that figure is displayed in
        % other figures that may have been clicked by the user)
        figure(fid);
        % make figure invisible
        set(fid, 'visible', 'off')
        %h = imshow(img_tmp_gr);        
        %h = imagesc(img_tmp_gr); colormap(gray); % CHANGE!
        imagesc(img_tmp_gr); colormap(gray);
        axis off
        
        res = 50; % pixels per inch
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 size(img_tmp_gr, 2) size(img_tmp_gr, 1)]/res);
        %saveas(h, fullfile(namePNG, ['overlay_', num2str(z_tmp), '.png']), 'png') % CHANGE!
        saveas(gcf, fullfile(namePNG, ['overlay_', num2str(z_tmp), '.png']), 'png')
      
    end
    
end

close(fid)
cd(orpth)
nameGIF = nameGIF((regexp(nameGIF, '/')+1):end);

end