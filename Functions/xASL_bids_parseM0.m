function xASL_bids_parseM0(PathNifti)
%xASL_bids_parseM0 This function check the M0 possibilities and will
%convert them to the ExploreASL legacy format

%PathNifti should be ASL Nifti

[Fpath, Ffile] = xASL_fileparts(PathNifti);
PathJSON = fullfile(Fpath, [Ffile '.json']);

% if ~isempty(regexpi(Ffile, 'm0'))
%     NiftiIs = 1;
% elseif ~isempty(regexpi(Ffile, 'asl'))
%     NiftiIs = 2;
% else
%     % no m0scan or asl found, skipping
%     return;
% end
%     

%% Parse & process M0 options

JSON = xASL_import_json(PathJSON);
if isfield(JSON, 'M0Type')
    switch JSON.M0Type
        %% Option 1. Included
        case 'Included'
            PathContext = fullfile(Fpath, [Ffile '_aslcontext.tsv']);
            if exist(PathContext, 'file')
               TSV = xASL_tsvRead(PathContext);
               M0Index = find(strcmp(TSV, 'm0scan'))-1;
               if ~isempty(M0Index)
                   xASL_io_SplitASL_M0(PathNifti, M0Index);
                   JSON.M0 = 'separate_scan';
                   spm_jsonwrite(PathJSON, JSON);
               else
                   warning(['M0Index missing in ' PathContext]);
               end
            else
                warning([PathContext ' missing']);
            end
        %% Option 2. Separate
        case 'Separate'
            PathM0 = fullfile(Fpath, 'M0.nii');
            if ~xASL_exist(fullfile(PathM0))
                warning(['Missing: ' PathM0]);
            end
        %% Option 3. Estimate
        case 'Estimate'
            if isfield(JSON, 'M0Estimate')
                JSON.M0 = JSON.M0Estimate;
                spm_jsonwrite(PathJSON, JSON);
            else
                warning(['Field M0_value missing in ' PathJSON]);
            end
        %% Option 4. Absent
        % (in this case we use the control image as pseudo-M0, 
        % which happens in xASL_wrp_Resample
        case {'use_control_as_m0', 'Absent'}
            if isfield(JSON, 'BackgroundSuppression')
                if JSON.BackgroundSuppression
                    warning('Using mean control as M0 but background suppression was present');
                end
                JSON.M0 = 'UseControlAsM0';
                spm_jsonwrite(PathJSON, JSON);
            else
                warning(['Missing field backgroundSuppression in ' PathJSON]);
            end
        otherwise
            error(['Invalid M0Type in ' PathJSON]);
    end
else
    warning(['Field M0Type missing in ' PathJSON]);
end
        


end