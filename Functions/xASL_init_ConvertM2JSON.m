function [PathJSON] = xASL_init_ConvertM2JSON(DataParFile)
%xASL_init_ConvertM2JSON Convert the old DataPar m-file to JSON format
% DataPar is the settings/parameter file, specific to a dataset to be
% processed by ExploreASL

[Fpath, Ffile] = fileparts(DataParFile);
CurD = pwd;
if ~isempty(Fpath)
    cd(Fpath);
end
PathJSON = fullfile(Fpath, [Ffile '.json']);
FunctionHandle = str2func(Ffile);
try % first try normal running of the m-file
    x = FunctionHandle();
% catch ME1
%     try % then try copying the m-file into a separate temporary file
%         % Bypass eval error stuff with long names, spaces etc
%         TempDataParFile = 'TempDataPar';
%         TempDataParPath = fullfile(pathstr,[TempDataParFile '.m']);
%         while exist(TempDataParPath,'file') % find unexisting file
%             TempDataParFile = [TempDataParFile '_2'];
%             TempDataParPath = fullfile(pathstr,[TempDataParFile '.m']);
%         end
%         copyfile(DataParPath, TempDataParPath, 'f' ); % overwrite ,if exist    
    catch ME2
%         fprintf('%s\n',ME1.message);
        fprintf('%s\n',ME2.message);
        error('Couldnt load the data par file');
%     end
end

cd(CurD);

% Escape '\' for subject_regexp field
FieldsX = fields(x);
for iField=1:length(FieldsX)
    tString = x.(FieldsX{iField});
    
    if ~isstruct(tString)
        % fix Boolean conversion issue JSON
        if islogical(tString) && tString
            tString = 1;
        elseif islogical(tString) && ~tString
            tString = 0;
		end
        
        % Escape illegal characters
		if ~iscell(tString)
			Strs = find(tString=='\');
			for iStr=1:length(Strs)
				tString = [tString(1:Strs(iStr)-1) '\\' tString(Strs(iStr)+1:end)];
				Strs = Strs+1;
			end
		else
			for ii=1:length(tString)
				Strs = find(tString{ii}=='\');
				for iStr=1:length(Strs)
					tString{ii} = [tString{ii}(1:Strs(iStr)-1) '\\' tString{ii}(Strs(iStr)+1:end)];
					Strs = Strs+1;
				end
			end
		end
        x.(FieldsX{iField}) = tString;
    end
end

xASL_delete(PathJSON);
spm_jsonwrite(PathJSON, x);

end

