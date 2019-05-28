function [analysis_dir, anats, functs] = get_neuroimaging_files(analysis_dir)

% Note: this function is only compatible with 4D functionals

if nargin==0 || isempty(analysis_dir) || ~isdir(analysis_dir)
    % Select main analysis folder:
    dialog_title = 'Select main analysis folder...';
    analysis_dir = uigetdir(pwd,dialog_title); % get directory
    if analysis_dir==0; error('User cancelled preprocessing.'); end
end

% Get Subject Folders:
folders = dir(fullfile(analysis_dir,'sub*'));
count = 0; subjDirs = {};
for i = 1:length(folders)
    if folders(i).isdir && isdir(fullfile(analysis_dir,folders(i).name))
        count = count+1;
        subjDirs{count} = fullfile(analysis_dir,folders(i).name); %#ok
    end
end

nSubjects = length(subjDirs);

% Get Anat:
anats = cell(1,nSubjects);
for i = 1:nSubjects
    fpath = fullfile(subjDirs{i},'anat');
    listing = dir(fullfile(fpath,'sub*'));
    for j = 1:length(listing)
        if exist(fullfile(fpath,listing(j).name),'file')==2
            anats{i} = fullfile(fpath,listing(j).name);
            break;
        end
    end
end

% Get Functs:
functs = cell(1,nSubjects);
for i = 1:nSubjects
    fpath = fullfile(subjDirs{i},'funct');
    listing = dir(fullfile(fpath,'sub*'));
    for j = 1:length(listing)
        if exist(fullfile(fpath,listing(j).name),'file')==2
            idx = strfind(listing(j).name,'run-');
            runNum = str2double(listing(j).name(idx+4:idx+5));
            functs{i}{runNum} = fullfile(fpath,listing(j).name);
        end
    end
end

end