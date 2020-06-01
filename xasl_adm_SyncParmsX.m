function [x] = xasl_adm_SyncParmsX(Parms, x, bVerbose)
%xasl_adm_SyncParmsX Sync Parms.* with x.(Q.)* (overwrite x/x.Q)
% Input all fields from this single subject/session/run Parms into the x structure (inheritance principle)

%% Admin
if nargin<3 || isempty(bVerbose)
    bVerbose = false;
end

%% Define fields to fill or dive in
ParmsNames = fieldnames(Parms);

%% Do it
for iPar=1:length(ParmsNames)
    if isstruct(Parms.(ParmsNames{iPar}))
        % we go into this struct & copy it
        if ~isfield(x, ParmsNames{iPar})
            % first we create empty X field
            x.(ParmsNames{iPar}) = struct;
        end
        
        x.(ParmsNames{iPar}) = xasl_adm_SyncParmsX(Parms.(ParmsNames{iPar}), x.(ParmsNames{iPar}));
    else % overwrite X by Parms
        x.(ParmsNames{iPar}) = Parms.(ParmsNames{iPar});
        if bVerbose
            if iscell(Parms.(ParmsNames{iPar})) % && length(Parms.(ParmsNames{iPar}))>1
                fprintf('%s\n', ['Overwritten: ' ParmsNames{iPar} ' in x by cell contents']);
            elseif isnumeric(Parms.(ParmsNames{iPar})) && length(Parms.(ParmsNames{iPar}))>1
                fprintf('%s\n', ['Overwritten: ' ParmsNames{iPar} ' in x by numerical table']);
            else
                fprintf('%s\n', ['Overwritten: ' ParmsNames{iPar} ' in x by ' xASL_num2str(Parms.(ParmsNames{iPar}))]);
            end
        end
    end
end


end