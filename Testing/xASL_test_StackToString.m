%% Stack to string
function stackText = xASL_test_StackToString(stack)

    % Fallback
    stackText = '';
    
    % Iterate over elements
    for iElement = 1:size(stack,1)
        thisFile = stack(iElement).file;
        thisName = stack(iElement).name;
        thisLine = stack(iElement).line;
        stackText = [stackText ', ' thisName ': line ' num2str(thisLine)];
    end
    
    % Remove initial ' ,'
    if ~isempty(stackText) && length(stackText)>3
        stackText = stackText(3:end);
    end

end