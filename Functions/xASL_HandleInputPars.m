function [StructIn DidntContain] = xASL_HandleInputPars(StructIn,FieldName,DefaultV)
%xASL_HandleInputPars Summary of this function goes here
%   Detailed explanation goes here

DidntContain        = false;
if     ~isfield(StructIn,FieldName)
        DidntContain    = true;
elseif  isempty(StructIn.(FieldName))
        DidntContain    = true;
end

if  DidntContain
    StructIn        = setfield(StructIn,FieldName,DefaultV); % create field with default value
end


end

