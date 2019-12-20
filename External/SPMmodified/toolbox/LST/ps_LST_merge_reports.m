function varargout = ps_LST_merge_reports(varargin)
%ps_LST_merge_reports   Merge HTML reports of different segmentations or
%   subjects.
%   Part of the LST toolbox, www.statistical-modeling.de/lst.html
%
%   ps_LST_merge_reports This function allows to merge multiple HTML 
%   reports that have been obtained by different functions of this toolbox.
%   The resulting report can be moved anywhere on your desk as long as it 
%   is provided that the folders for the original reports stay where they 
%   are. A report that is readable on different platforms can be obtained 
%   by exporting/printing the HTML report as a PDF document, a function 
%   that is included in all modern browsers.
%   
%   ps_LST_merge_reports(reports) merges the reports given in the argument 
%   reports. This must be a character, like a call from spm_select.
%   
%   ps_LST_merge_reports(job) Same as above but with job being a harvested 
%   job data structure (see matlabbatch help).
%

% Check input
if ~isempty(varargin) && isfield(varargin{1}, 'data_rep')
    viajob = 1;
else
    viajob = 0;
end
if ~viajob
    if nargin == 0                
        Vrep = spm_select(Inf, 'any', 'Select reports to be merged.');        
    end    
    if nargin == 1
        if ischar(varargin{1})
            Vrep = varargin{1};       
        else            
            fprintf('Input must be character, like from spm_select.\n');
            fprintf('See ?ps_LST_merge_reports for help.\n');            
        return;
        end
    end
    if nargin > 2        
        fprintf('Too many input arguments!\n');
        fprintf('See ?ps_LST_long for help.\n');        
        return;        
    end
else
    job = varargin{1};
    Vrep = job.data_rep;    
end

nameHTML = ['report_LST_', num2str(size(Vrep, 1)), '_', ps_create_timestamp, '.html'];

% Main HTML file
copyfile(fullfile(spm('dir'), 'toolbox', 'LST', 'LST_main_html.html'), nameHTML)
HTMLid = fopen(nameHTML, 'at');
strout = ['    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'raphael.js'), '\"></script>\n', ...
      '    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery.min.js'), '\"></script>\n', ...
      '    <link href=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.css'), '\" rel=\"stylesheet\"></script>\n', ...
      '    <script src=\"', fullfile(spm('dir'), 'toolbox', 'LST', 'js', 'jquery-ui.js'), '\"></script>\n', ...
      '  </head>\n  <body>\n'];
fprintf(HTMLid, strout);

for i = 1:size(Vrep, 1)
    if viajob
        [pthtmp, namtmp, ~] = fileparts(Vrep{i,:});
    else        
        [pthtmp, namtmp, ~] = fileparts(Vrep(i,:));
    end
        
    ul = regexp(namtmp, '_');
    id = namtmp((ul(1)+1):end);    
    strout = fileread(fullfile(pthtmp, id, [id, '.html']));
    strout = strrep(strout, '%', '%%');
    fprintf(HTMLid, strout);
end

strout = '  </body>\n</html>';
fprintf(HTMLid, strout);
fclose(HTMLid);

end