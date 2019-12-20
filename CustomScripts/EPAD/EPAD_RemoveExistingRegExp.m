function [ThroughputString] = EPAD_RemoveExistingRegExp(ThroughputString, InputRegExp)
%EPAD_RemoveExistingRegExp % Remove requested regular expression match from string
% This is useful when we want to add a suffix or prefix to a folder or ListName, but make
% sure that we are not getting multiple suffixes or prefixes

if nargin<2 || isempty(InputRegExp)
    return;
end

if nargin<1 || isempty(ThroughputString)
    return;
end

if isnumeric(ThroughputString)
    return;
end
    
[Ind1, Ind2] = regexp(lower(ThroughputString),lower(InputRegExp), 'start', 'end');

while ~isempty(Ind1) && ~isempty(Ind2)
    ThroughputString = [ThroughputString(1:Ind1(1)-1) ThroughputString(Ind2(1)+1:end)];
    [Ind1, Ind2] = regexp(lower(ThroughputString),lower(InputRegExp), 'start', 'end');
end
    

end