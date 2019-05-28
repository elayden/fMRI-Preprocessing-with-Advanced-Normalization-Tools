% append_folder_markers.m
% 
% INFO: this script allows a user to select any number of parent folders.
%   Parent folders contain subfolders which should have matching numbers 
%   across all parent folders. Any non-matching subfolders will be marked 
%   with a specified prefix appended to their titles (e.g., this is 
%   useful to check whether different datasets contain the same subject 
%   numbers).
% 
% Author:       Elliot Layden, University of Chicago, September 13, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nParentFolders = 3;
subfolder_prefix = ''; % this should be a char/string denoting any text 
%   that comes before the subfolder numbers in the subfolder names; if 
%   none, leave empty ('')
marker = 'missing_'; % this is what will be appended to folder names

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentDir = pwd;
% Select Parent Folders:
parent_dirs = cell(1,nParentFolders);
for i = 1:nParentFolders
    dialog_title = sprintf('Select parent folder %02g',i);
    parent_dirs{i} = uigetdir(currentDir, dialog_title); % get directory
    if parent_dirs{i}==0; error('User cancelled preprocessing.'); end
    [currentDir,~,~] = fileparts(parent_dirs{i});
end

% Get Subfolders:
subdir_listings = cell(1,nParentFolders);
for i = 1:nParentFolders
    subdir_listings{i} = dir(fullfile(parent_dirs{i},[subfolder_prefix,'*']));
    removeIdx = false(length(subdir_listings),1);
    for j = 1:length(subdir_listings{i})
        if ~subdir_listings{i}(j).isdir || strcmp(subdir_listings{i}(j).name,'.') || strcmp(subdir_listings{i}(j).name,'..')
            removeIdx(j) = true;
        end
    end; subdir_listings{i}(removeIdx) = [];
end

% Get Numbers:
subdir_nums = cell(1,nParentFolders);
subdir_mark = cell(1,nParentFolders);
for i = 1:nParentFolders
    nSubdir = length(subdir_listings{i});
    subdir_nums{i} = zeros(1,nSubdir);
    subdir_mark{i} = false(1,nSubdir);
    for j = 1:nSubdir
        % Find subfolder prefix start
        if ~isempty(subfolder_prefix)
            ix = strfind(subdir_listings{i}(j).name,subfolder_prefix); 
            ix_begin = ix + length(subfolder_prefix);
        else ix_begin = 1;
        end
        store_num = '';
        for k = ix_begin:length(subdir_listings{i}(j).name)
            if ~isnan(str2double(subdir_listings{i}(j).name(k)))
                store_num = [store_num,subdir_listings{i}(j).name(k)]; %#ok
            else break;
            end
        end
        subdir_nums{i}(j) = str2double(store_num);
    end
end
   
% Check to see that subject numbers were retrieved correctly:
% subdir_nums{1}
% subdir_nums{2}
% subdir_nums{3}

% Check for unique numbers:
for i = 1:nParentFolders
    for j = 1:nParentFolders
        if i~=j
            [~,ix] = setdiff(subdir_nums{i},subdir_nums{j});
            if ~isempty(ix) && ~any(isnan(ix))
               subdir_mark{i}(ix) = true; 
            end
        end
    end
end

% Check to make sure the marking indices (boolean,1/0), worked correctly
% (if all are marked, there may be trouble):
%  subdir_mark{1}
%  subdir_mark{2}
%  subdir_mark{3}

% Mark Unique Subdirs:
for i = 1:nParentFolders
    for j = 1:length(subdir_listings{i})
        if subdir_mark{i}(j)
            origName = fullfile(parent_dirs{i},subdir_listings{i}(j).name);
            newName = fullfile(parent_dirs{i},[marker,subdir_listings{i}(j).name]);
            status = copyfile(origName,newName);
            if status
                status2 = rmdir(origName,'s');
                if ~status2
                    error(['ERROR: Removing previous folder ',origName,' failed.']);
                end
            else
                error(['ERROR: Copying folder ',origName,' failed.']);
            end
                
        end
    end
end