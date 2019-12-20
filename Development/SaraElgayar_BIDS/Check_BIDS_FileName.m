function check = Check_BIDS_FileName(FileName)
% all File names should follow the BIDS standards and be in this formate:
%
%                 [_keyword-<value>]_<suffix>.<ext>
%           for example: sub-001_label-WMH_atlas-LPA_probseg.nii
%              keywords: sub, label, atlas, 
%                values: 001, WMH, LPA
%                suffix:  probseg, dseg
%                   ext:  nii
%
% This script checks the validity of BIDS file name
% it takes a string (file name) and returns true or false 
% also it should output a report file with errors
%  

%keywords = '\w{1,}';
%values = '(\d*|\w+)';  %\d* = any number of digits  \w{1,} = contains one or more characters
%suffix = '\w{0,}';  % example =probseg
%ext = '\w{3,}'; 
%pattern = '('keywords'-'values'_')'+(suffix)\.ext' ; % * 0 or more + %1 or more


%pattern = '^sub-(\d*|\w*)(_\w{1,}-(\d*|\w*))*(_\w)*\.(nii|nii\.gz|csv|mat|json)$';
%pattern = '(keyword-(value)_)+(suffix)\.(ext)';
%pattern = '^sub-(\d+|w*)(_label[-](\d+|\w*))*(_atlas-(\d+|\w*))*(_probseg)*\.(nii|nii\.gz|csv|mat|json)$';

%rec = reconstracution 
%space = ACPC
%ce = Contrast Enhancing Agent 

pattern = '^sub-[0-9a-zA-Z]+(_(ses|label|acq|atlas|space|rec|desc)-[0-9a-zA-Z]+)*(_[0-9a-zA-Z]+)*\.(nii|nii\.gz|csv|mat|json|tsv)$';
matchStr = regexp(FileName, pattern, 'match');
if string(matchStr) == FileName
    check = true;
    disp('BIDS valid');
else
    check = false;
    disp('BIDS NOT valid');
end
