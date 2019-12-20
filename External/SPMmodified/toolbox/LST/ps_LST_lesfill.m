function varargout = ps_LST_lesfill(varargin)
%ps_LST_lesfill   Fill lesions
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
%
%   ps_LST_lesfill Lesion filling can be applied to any image that is in 
%   alignment with the lesion probability map. However, it is required that
%   the .mat-files that are saved during the lesion segmentation process 
%   are available in the same folders as the lesions probability maps.
%   
%   ps_LST_lesfill calls spm_select in order to choose the images that need
%   to be filled (images in native space) and the lesion probability maps.
%   
%   ps_LST_lesfill(Vns, Vles, void, html) fills the native space images in
%   Vns along the lesions given in Vles. Both arguments need to be
%   characters, like a call from spm_selct. The last two arguments indicate
%   weather the algorithm should produce output on the command line and if
%   a HTML report should be produced. These arguments can also be left empty.
%
%   ps_LST_lesfill(job) Same as above but with job being a harvested job data
%   structure (see matlabbatch help).
%


% Check input
if ~isempty(varargin) && isfield(varargin{1}, 'data')
    viajob = 1;
else
    viajob = 0;
end
void = 0;
if ~viajob
    if nargin == 0
        % Select files by function
        Vns = spm_select(Inf, 'image', 'Select images in native space.');
        Vles = spm_select(Inf, 'image', 'Select probability lesion maps.');
        html_report = 1;
    end
    if nargin == 1        
        fprintf('Please give me at least two arguments.\n');
        fprintf('See ?ps_LST_lesfill for help.\n');        
        return;        
    end
    if nargin > 1
        if ischar(varargin{1})
            Vns = varargin{1};
        else            
            if isempty(varargin{1})
                Vns = spm_select(Inf, 'image', 'Select images in native space.');
            else
                fprintf('Input for Vns must be a character, like from spm_select.\n');
                return;
            end
        end
        if ischar(varargin{2})
            Vles = varargin{2};
        else            
            if isempty(varargin{2})
                Vles = spm_select(Inf, 'image', 'Select probability lesion maps.');
            else
                fprintf('Input for Vles must be a character, like from spm_select.\n');            
                return;
            end
        end
        html_report = 1;
        void = 0;
    end
    if nargin == 3
        void = varargin{3};
        html_report = 1;
    end
    if nargin == 4
        void = varargin{3};
        html_report = varargin{4};
    end
    if nargin > 4        
        fprintf('To many input arguments. Did you mean ps_LST_lga?\n');        
        return;
    end
else    
    job = varargin{1};    
    Vns = job.data;
    Vles = job.data_plm; 
    html_report = job.html_report;
end

if ~void
    % Welcome text
    fprintf(repmat('=', 1, 72));
    fprintf('\nRun lesion filling. Visit www.statistical-modeling.de/lst.html \n')
    fprintf('for updates and more information.\n')
    fprintf(repmat('=', 1, 72));
end

% Load volumes
Vns = spm_vol(Vns);
Vles = spm_vol(Vles); 

% Check input
if ~isequal(numel(Vns), numel(Vles))
    fprintf('\nNumber of images must match number of lesion maps.\n'); 
    return;
end
if numel(Vns) == 0 || numel(Vles) == 0
    fprintf('\nNo images selected.\n');
    return;
end

% Summarize input
if ~void
    strout = '\n\nNumber of images to process:';
    fprintf(strout)
    tt = [num2str(numel(Vns)), '\n'];
    strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 6), tt];
    fprintf(strout)
end

for i = 1:numel(Vns)
    
    if viajob
        Vns_tmp = Vns{i};
        Vles_tmp = Vles{i};    
    else    
        Vns_tmp = Vns(i);
        Vles_tmp = Vles(i);
    end
    [pthns, namns, extns] = fileparts(Vns_tmp.fname);
    [pthles, namles, ~] = fileparts(Vles_tmp.fname);
    cd(pthles)
    
    % Which subject?   
    if ~void
        strout = '\nWorking on image';
        fprintf(strout)
        tt = [num2str(i), ' out of ', num2str(numel(Vns)), ' (', num2str(i/numel(Vns)*100), '%%)\n'];
        strout = [repmat(' ', 1, 72 - numel(tt) - numel(strout) + 5), tt];
        fprintf(strout)
    end
    
    % lga or lpa?
    lga = regexp(namles, '_lga_')> 0;
    if isempty(lga); lga = 0; end;
    us_pos = regexp(namles, '_');
    %nam = namles((us_pos(1)+1):(us_pos(end)-1));
    nam = namles((us_pos(1)+1):end);
    if lga
        namf2 = namles((us_pos(3)+1):end);
    else
        namf2 = namles((us_pos(2)+1):end);
    end    
    if lga
        load(['LST_lga_', namf2, '.mat'])
        lst = lga;
        lga = 1;
    else
        load(['LST_lpa_', namf2, '.mat'])
        lst = lpa;
    end
    or = [lst.or, lst.fl];
    
    % load data
    img = spm_read_vols(Vns_tmp);
    les = spm_read_vols(Vles_tmp);
    
    % hard segmentation   
    seg = 0 .* img;
    if lga
        seg(lst.indx_brain) = lst.p0_vec;
        seg(seg > 0 & seg < 1.5) = 1;
        seg(seg >= 1.5 & seg < 2.5) = 2;
        seg(seg >= 2.5) = 3;
    else
        seg(lst.indx_brain) = lst.I;
    end
    c_tmp = indx2coord(lst.indx_brain, size(seg, 1), size(seg, 2));
    %r = [min(c_tmp(:,3)), max(c_tmp(:,3))];
    r = [min(c_tmp(:,find(or == 3))), max(c_tmp(:,find(or == 3)))];
    
    % filling
    filled = ps_LST_lesfill_core(img, les, seg, lst.indx_brain);    
    Vns_tmp.fname = fullfile(pthns, [namns, '_filled_', nam, extns]);
    spm_write_vol(Vns_tmp, filled);
    
    if html_report
        
        % HTML report
        % -----------------------------------------------------------------    
                    
        stroutHTML = 'Create HTML report';
        fprintf(stroutHTML);
        tic

        % create HTML report
        nameFolder = ['LST_filled_' namns, '_', nam];
        warning('off');
        mkdir(nameFolder)
        warning('on');
        
        % Create PNGs        
        Vimg1 = fullfile(cd, [namns, extns]);
        Vimg2 = Vns_tmp.fname;
        if numel(Vimg1) < numel(Vimg2)
            Vimg1 = [Vimg1, repmat(' ', 1, numel(Vimg2) - numel(Vimg1))];
        else
            Vimg2 = [Vimg2, repmat(' ', 1, numel(Vimg1) - numel(Vimg2))];
        end
        Vimg = [Vimg1; Vimg2];
        pngFailed = '';
        try
            [~, r] = ps_LST_create_gif(Vimg, '', nameFolder, r, or);
        catch ME
            pngFailed = ME.message;
            r = 0:1;
        end
        
        % Main HTML file
        nameHTML = ['report_LST_filled_', namns, '_', nam, '.html'];
        copyfile(fullfile(spm('dir'), 'toolbox', 'LST', 'LST_main_html.html'), nameHTML)
        HTMLid = fopen(nameHTML, 'at');
        strout = ['    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'raphael.js'), '\"></script>\n', ...
          '    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery.min.js'), '\"></script>\n', ...
          '    <link href=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.css'), '\" rel=\"stylesheet\"></script>\n', ...
          '    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.js'), '\"></script>\n', ...
          '  </head>\n  <body>\n'];
        fprintf(HTMLid, strout);       

        % create subject specific html file
        if lga
            algtext = 'LGA';
        else
            algtext = 'LPA'; 
        end
        jsid = [nameFolder, '_', ps_create_timestamp];
        jsid(regexp(jsid, '\.')) = []; 

        strout = ['\n<div class=\"container\">\n', ...              
          '  <h1>Lesion filling by LST</h1>\n', ...
          '  <div class=\"column-01\">\n', ...
          '    <h2>Input summary</h2>\n', ...
          '      <table style=\"min-width: 500px;\">\n', ...
          '        <tr>\n', ...
          '            <td>Date of analysis</td>\n', ...
          '            <td class=\"ta_right\">', datestr(clock()), '</td>\n', ...
          '        </tr>\n', ...
          '       <tr>\n', ...
          '           <td>Algorithm used for segmentation</td>\n', ...
          '           <td class=\"ta_right\">', algtext, '</td>\n', ...
          '       </tr>\n', ...          
          '        <tr>\n', ...
          '          <td>Image</td>\n', ...
          '          <td class=\"ta_right\">', ps_shorten_string(fullfile(pthns, namns), 28), '</td>\n', ...
          '       </tr>\n', ...
          '        <tr>\n', ...
          '          <td>Lesion probabilty map</td>\n', ...
          '          <td class=\"ta_right\">', ps_shorten_string(fullfile(pthles, namles), 28), '</td>\n', ...
          '       </tr>\n', ...          
          '   </table>\n', ...
          '   </div>\n', ...              
          '   <div class=\"column-02\">\n', ...
          '     <h2>Results</h2>\n', ...
          '         <table style=\"width: 500px\">\n', ...          
          '             <tr>\n', ...
          '                 <td>Filled image</td>\n', ...
          '                 <td class=\"ta_right\">', [namns, '_filled_', nam, extns], '</td>\n', ...
          '             </tr>\n', ...          
          '        </table>\n', ...
          '    </div>\n', ...
          '    <div style=\"clear:both\"></div>\n', ...
          '        <h2>Overlay</h2>\n'];
      if strcmp(pngFailed, '')
          strout = [strout, ...
              '        <script src=\"', fullfile(cd, nameFolder, 'filled.js'), '\" type=\"text/javascript\"></script>\n', ...          
              '        <img width=\"450px\" id=\"overlay', jsid, '\" src=\"', fullfile(cd, nameFolder, ['overlay_', num2str(round(mean(r))), '.png']), '\" />\n', ...
              '        <div id=\"slider_', jsid, '\" style=\"width: 450px; text-align: center;\"></div>\n', ...
              '        <div style=\"width: 450px; text-align: center;\">\n', ...
              '           <button id=\"button-left', jsid, '\">\n', ...
              '              <\n', ...
              '           </button>\n', ...
              '           <span id="slice', jsid, '">Slice ', num2str(round(mean(r))), '</span>\n', ...
              '           <button id=\"button-right', jsid, '\">\n', ...
              '              >\n', ...
              '           </button>\n', ...
              '        </div>\n'];
      else
          strout = [strout, ...
              'Sorry, there was a problem when creating the PNG images. MATLAB said: ',  pngFailed, '\n'];
      end
      strout = [strout, ...
          '  </div>\n', ...
          '<br><hr>\n']; %% !!
      
        JSid = fopen(fullfile(nameFolder, 'filled.js'), 'wt');
        js_strout = ['$(function () {\n', ...             
          '  var min_slice', jsid, ' = ', num2str(r(1)), ',\n', ...
          '  max_slice', jsid, ' = ', num2str(r(2)), ',\n', ...
          '  slice', jsid, ' = ', num2str(round(mean(r))), ';\n', ...
          '  $( \"#slider_', jsid, '\" ).slider({\n', ...
          '    min: min_slice', jsid, ',\n', ...
          '    max: max_slice', jsid, ',\n', ...
          '    value: slice', jsid, ',\n', ...
          '    slide: function( event, ui ) {\n', ...
          '      slice', jsid, ' = ui.value;\n', ...
          '      $( \"#overlay', jsid, '\" ).attr(\"src\", \"', fullfile(cd, nameFolder, 'overlay_\" + ui.value + \".png\"'), ');\n', ...
          '      $( \"#slice', jsid, '\" ).text(\"Slice \" + ui.value);\n', ...
          '    }\n', ...
          '  });\n', ...
          '  $( \"#button-left', jsid, '\" ).button({\n', ...
          '    icons: {\n', ...
          '      primary: \"ui-icon-carat-1-w\"\n', ...
          '    },\n', ...
          '    text: false,\n', ...
          '    }).click(function(event, ui){\n', ...
          '        if(slice', jsid, ' > min_slice', jsid, '){\n', ...
          '            slice', jsid, ' = slice', jsid, ' - 1;\n', ...
          '            $( \"#overlay', jsid, '\" ).attr(\"src\", \"', fullfile(cd, nameFolder, ['overlay_\" + slice', jsid, ' + \".png\"']), ');\n', ...
          '            $( \"#slice', jsid, '\" ).text(\"Slice \" + slice', jsid, ');\n', ...
          '            $(\"#slider_', jsid, '\").slider(\"option\", "value", slice', jsid, ');\n', ...
          '        }\n', ...
          '  });\n', ...
          '  $(\"#button-right', jsid, '\").button({\n', ...
          '    icons: {\n', ...
          '      primary: \"ui-icon-carat-1-e\"\n', ...
          '    },\n', ...
          '    text: false,\n', ...
          '    }).click(function(event, ui){\n', ...
          '        if(slice', jsid, ' < max_slice', jsid, '){\n', ...
          '            slice', jsid, ' = slice', jsid, ' + 1;\n', ...
          '            $( \"#overlay', jsid, '\").attr(\"src\", \"', fullfile(cd, nameFolder, ['overlay_\" + slice', jsid, ' + \".png\"']), ');\n', ...
          '            $( \"#slice', jsid, '\").text(\"Slice \" + slice', jsid, ');\n', ...
          '            $(\"#slider_', jsid, '\").slider(\"option\", \"value\", slice', jsid, ');\n', ...
          '        }\n', ...
          '  });', ...
          '});'];

        fprintf(JSid, js_strout);
        fclose(JSid);

        %HTMLid = fopen(fullfile(nameFolder, [id, '.html']), 'wt');
        fprintf(HTMLid, strout);

        HTMLid2 = fopen(fullfile(nameFolder, ['LST_filled_', namns, '_', nam, '.html']), 'wt');
        fprintf(HTMLid2, strout);
        fclose(HTMLid2);

        strout = '  </head>\n  <body>\n';
        fprintf(HTMLid, strout);
        fclose(HTMLid);        

        tt = toc; tt = [num2str(round(tt)), 's'];
        stroutHTML = [repmat(' ', 1, 72 - numel(tt) - numel(stroutHTML)), tt, '\n'];
        fprintf(stroutHTML)


    end
end

c = clock();
if ~void
    disp(['Finished successfully on ', datestr(c)]);
    fprintf(repmat('-', 1, 72));
    fprintf('\n')
end
varargout{:} = [namns, '_filled_', nam, extns];


