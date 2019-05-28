function [BIDS_output_folder, anats, functs] = organize_neuroimaging_files_point_click
% % [BIDS_output_folder, ANATS, FUNCTS] = ORGANIZE_NEUROIMAGING_FILES;
% 
% INFO: Provides an easy way to organize neuroimaging data into rough
% correspondence with the Brain Imaging Data Structure (BIDS) format
% <http://bids.neuroimaging.io/>.
% Users are directed through a series of file explorer menus and dialogue
% boxes which allow automated point-and-click transfer and renaming of 
% files into the BIDS format directories. Additionally, a 
% 'file_correspondence_info.txt' file is automatically written to the main 
% BIDS directory, allowing users to verify that the correct files were 
% selected, renamed, and copied to the appropriate destination folders.
% 
% Note: this function does not delete any of the original files.
% 
% OUTPUTS,
%   BIDS_output_folder,   string: the full path to the main directory
%   anats,      a cell array wherein each cell contains a string denoting
%               the full filepath to each anatomical image in the BIDS 
%               directories 
%   functs,     a cell array wherein each cell contains a string denoting
%               the full filepath to each functional image in the BIDS 
%               directories; the structure of this cell array will vary
%               depending on the study specifics (i.e., number of runs,
%               etc.)
% 
% Example of BIDS directories:
% % Main Folder
%     sub-01
%         funct    - original functionals
%         anat     - original structurals
%         physio   - scan physio data
%         qcreport  - where the results of preprocessing will be written
%             funcProc
%                sub-%02g_run-%02g_bold.nii
%                wsub-%02g_run-%02g_bold.nii
%                sub-%02g_run-%02g_bold.mat      % affine matrix 
%                rp_sub-%02g_run-%02g_bold.txt   % time (rows) x 6 rigid-body parameter (columns) matrix of realignment motion parameters
%             anatProc
%                sub-%02g_anat.nii
%                y_sub-%02g_anat.nii  % deformation field?
%                msub-%02g_anat.nii
%                c1sub-%02g_anat.nii ... c6sub-%02g_anat.nii
%    ...
%    sub-<last>
% 
% Author:       Elliot Layden, University of Chicago, August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anat_suffix = 'T1w'; % 'T2w'; % 'anat'

orig_dir = pwd;

% Select main output folder:
dialog_title = 'Select main analysis folder...';
BIDS_output_folder = uigetdir(cd,dialog_title); % get directory
if BIDS_output_folder==0; error('User cancelled preprocessing.'); end
            
% Input number of subjects:
answer = inputdlg('Number of Subjects: ');
if isempty(answer); error('User cancelled preprocessing.'); end
nSubjects = str2double(answer);
if ~isnumeric(nSubjects) && nSubjects>0
    error('Specify a positive integer for number of subjects.')
end

% Create subject sub-directories:
subdir = {'funct','anat','physio','qcreport'};
subjDirs = cell(1,nSubjects);
for i = 1:nSubjects
    if nSubjects<100
        subjDirs{i} = fullfile(BIDS_output_folder,sprintf('sub-%02g',i));
        if ~isdir(subjDirs{i}); mkdir(subjDirs{i}); end
        nzero = 2;
    elseif nSubjects<1000
        subjDirs{i} = fullfile(BIDS_output_folder,sprintf('sub-%03g',i));
        if ~isdir(subjDirs{i}); mkdir(subjDirs{i}); end
        nzero = 3;
    else
        subjDirs{i} = fullfile(BIDS_output_folder,sprintf('sub-%04g',i));
        if ~isdir(subjDirs{i}); mkdir(subjDirs{i}); end
        nzero = 4;
    end
    % Create subdirectories:
    for j = 1:4
        if ~isdir(fullfile(subjDirs{i},subdir{j})); mkdir(fullfile(subjDirs{i},subdir{j})); end
    end
end

%% Select Structurals:
qstring = 'Are the anatomical images in...';
str1 = 'Different Directories'; str2 = 'Same Directory';
struct_option = questdlg(qstring,'',str1,str2,str1);

structFiles = cell(1,nSubjects);
switch struct_option
    case str1 % multiple directories
        for i = 1:nSubjects
            [structName, structPath] = uigetfile({'*.nii';'*.nii.gz'},['Select anatomical image ',sprintf(['%0',num2str(nzero),'g:'],i)],'MultiSelect','off'); % open files
            if structName==0; error('User cancelled preprocessing.'); end
            structFiles{i} = fullfile(structPath,structName);
        end
    case str2 % single directory
        [structName, structPath] = uigetfile({'*.nii';'*.nii.gz';'*.img'},'Select anatomical images in order:','MultiSelect','on'); % open files
        if structName==0; error('User cancelled preprocessing.'); end
        for i = 1:nSubjects
            structFiles{i} = fullfile(structPath,structName{i});
        end
end

% Write out list of subject-file correspondences:
write_data = cell(nSubjects,1);
for i = 1:nSubjects
    write_data{i} = ['Anatomical ',sprintf(['%0',num2str(nzero),'g:  '],i),structFiles{i}];
end
write_data = char(write_data);
fileID = fopen(fullfile(BIDS_output_folder,'file_correspondence_info.txt'),'w');
fprintf(fileID,'%s\r\n','-------------------------------- Anatomicals ---------------------------------');
for i = 1:nSubjects
    fprintf(fileID,'%s\r\n',write_data(i,:));
end
fprintf(fileID,'%s\r\n',' ');
fclose(fileID);

% Copy over structs and rename:
[~,~,ext] = fileparts(structFiles{1});
anats = cell(1,nSubjects);
for i = 1:nSubjects
    anats{i} = fullfile(subjDirs{i},'anat',sprintf(['sub-%0',num2str(nzero),'g_',anat_suffix,ext],i));
    copyfile(structFiles{i},anats{i})
end

%% Select Functionals:

% Input number of runs:
answer = inputdlg('Number of Runs Per Subjects: ');
if isempty(answer); error('User cancelled preprocessing.'); end
nRuns = str2double(answer);
if ~isnumeric(nRuns) && nRuns>0
    error('Specify a positive integer for number of runs.')
end

qstring = 'Is a functional run...';
f2str1 = 'A single 4D image'; f2str2 = 'A series of 3D images'; 
funct_dim = questdlg(qstring,'',f2str1,f2str2,f2str1);

functFiles = cell(1,nSubjects);
nRuns_count = 0;
switch funct_dim
    case f2str1
        is4D = true;
        for i = 1:nSubjects
            [functName, functPath] = uigetfile({'*.nii';'*.nii.gz'},...
                ['Select all 4D functional runs for subject ',sprintf(['%0',num2str(nzero),'g'],i),...
                ' in order:'],'MultiSelect','on'); % open files
            if functPath==0; error('User cancelled preprocessing.'); end
            if length(functName)~=nRuns
                warning([num2str(nRuns),' runs specified, but ',num2str(length(functName)),' selected.']);
            end
            for j = 1:length(functName)
                functFiles{i}{j} = fullfile(functPath,functName{j});
                nRuns_count = nRuns_count + 1;
            end
        end
    case f2str2
        is4D = false;
        for i = 1:nSubjects
            for j = 1:nRuns
                if j == 1
                    [functName, functPath] = uigetfile({'*.nii';'*.nii.gz'},...
                        ['Select 3D functional files for subject ',sprintf(['%0',num2str(nzero),'g'],i),...
                        ', run ',sprintf('%02g:',j)],'MultiSelect','on'); % open files
                else
                    cd(functPath)
                    [functName, functPath] = uigetfile({'*.nii';'*.nii.gz'},...
                        ['Select 3D functional files for subject ',sprintf(['%0',num2str(nzero),'g'],i),...
                        ', run ',sprintf('%02g:',j)],'MultiSelect','on'); % open files
                end
                if functName==0; error('User cancelled preprocessing.'); end
                for k = 1:length(functName)
                    functFiles{i}{j}{k} = fullfile(functPath,functName{k});
                end
            end
        end
end

% Write out list of subject-file correspondences:
nWrite = nRuns_count;
write_data = cell(nWrite,1);
count = 0;
for i = 1:nSubjects
    for j = 1:length(functFiles{i})
        count = count+1;
        if iscell(functFiles{i}{j})
            write_data{count} = ['Functional, subject ',sprintf(['%0',num2str(nzero),'g, run %02g:   '],[i,j]),functFiles{i}{j}{1}];
        else
            write_data{count} = ['Functional, subject ',sprintf(['%0',num2str(nzero),'g, run %02g:   '],[i,j]),functFiles{i}{j}];
        end
    end
end
write_data = char(write_data);
fileID = fopen(fullfile(BIDS_output_folder,'file_correspondence_info.txt'),'a');
fprintf(fileID,'%s\r\n','-------------------------------- Functionals ---------------------------------');
for i = 1:size(write_data,1)
    fprintf(fileID,'%s\r\n',write_data(i,:));
end

fclose(fileID);

% Copy over functs and rename:
if is4D; 
    fname = functFiles{1}{1};
else
    fname = functFiles{1}{1}{1}; 
end
[~,~,ext] = fileparts(fname);
functs = cell(1,nSubjects);
hWait = waitbar(0,'Copying files...');
for i = 1:nSubjects
    waitbar(i/nSubjects,hWait);
    for j = 1:length(functFiles{i})
        if is4D
            functs{i}{j} = fullfile(subjDirs{i},'funct',sprintf(['sub-%0',num2str(nzero),'g_run-%02g_bold',ext],[i,j]));
            copyfile(functFiles{i}{j},functs{i}{j})
        else
            for k = 1:length(functFiles{i}{j})
                functs{i}{j}{k} = fullfile(subjDirs{i},'funct',sprintf(['sub-%0',num2str(nzero),'g_run-%02g_bold_%04g',ext],[i,j,k]));
                copyfile(functFiles{i}{j}{k},functs{i}{j}{k})
            end
        end
    end
end; close(hWait);

cd(orig_dir)

end