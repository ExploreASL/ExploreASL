function FileName = BIDS_FileName( SubID, ses, space, label, atlas, desc, suffix, ext)

BIDSName =  ['sub-' SubID '_ses-' ses '_space-' space '_label-' label '_atlas-' atlas '_desc-' desc '_' suffix '.' ext ];

if isempty(label) 
    BIDSName =  erase(BIDSName,'_label-');
end
if isempty(ses) 
    BIDSName =  erase(BIDSName,'_ses-');
end
if isempty(atlas)
     BIDSName =  erase(BIDSName,'_atlas-');
end
if isempty(desc)
     BIDSName =  erase(BIDSName,'_desc-');
end
if isempty(space)
     BIDSName =  erase(BIDSName,'_space-');
end
if isempty(suffix)
    % erasing the last suffix '-' 
     indices = strfind(BIDSName, '_');
     s1 = BIDSName(1:indices(end)-1);
     s2 = BIDSName(indices(end)+1:end);
     BIDSName = strcat(s1,s2);
end

if (Check_BIDS_FileName(BIDSName))
    FileName = BIDSName;
else 
    FileName = 'Error' ;
end
