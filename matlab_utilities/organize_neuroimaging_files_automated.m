% ORGANIZE_NEUROIMAGING_FILES_SCRIPT.m
% 
% INFO: Provides an easy way to organize neuroimaging data into the 
% Brain Imaging Data Structure (BIDS) format <http://bids.neuroimaging.io/>
% A 'file_correspondence_info.txt' file is automatically written to the main 
% BIDS directory, allowing users to verify that the correct files were 
% selected, renamed, and copied to the appropriate destination folders.
% 
% Note: this script does not delete any of the original files.
% 
% OUTPUTS to Matlab workspace:
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
%                sub-%02g_T1w.nii
%                y_sub-%02g_T1w.nii  % deformation field?
%                msub-%02g_T1w.nii
%                c1sub-%02g_T1w.nii ... c6sub-%02g_T1w.nii
%    ...
%    sub-<last>
% 
% Author:       Elliot Layden, University of Chicago, September 13, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Important: Inputs are case sensitive, and for file prefix inputs, the
% accurate placement of wild card patterns (e.g., *) is critical

% Specify these using full paths:
main_input_folder = 'C:\Users\Valued Customer\Desktop\Elliot\Hiro_Data\NTD - Copy'; % this is the main directory which contains the original data
BIDS_output_folder = 'C:\Users\Valued Customer\Desktop\Elliot\Hiro_Data\BIDS_test'; % this is the newly created directory to which BIDS outputs will be written

% Subdirectories (subject folders) of 'main_input_folder':
subject_folder_prefix = ''; % this should be a char/string denoting any text that comes before the subject number in the subject folder names; if none, leave empty ('')
exclude_prefix = 'missing_'; % this should specify any folder prefixes that should be avoided 
%  If you want to exclude certain subject folders, then append a
    %   different prefix to the folder names to "mark" them (other than the 
    %   one listed as 'subject_folder_prefix'); include this marking prefix
    %   here as 'exclude_prefix; leave empty if none ('')
%  Note: do not include a wildcard (*) asterisk for either of these inputs
%  Hint: you can use the 'append_folder_markers.m' script to identify any subject folders that are not consistent across a set of sessions

anat_folder = 'Nifti'; % this is the folder name of the folder which contains anatomicals; this must be consistent across subject folders
funct_folder = 'Nifti'; % this is the folder name of the folder which contains functionals; this must be consistent across subject folders
% Note: if both anatomicals and functionals are contained within the same
%   subfolder (e.g., 'Nifti'), then simply use the same folder name for both
%   'anat_folder' and 'funct_folder'

run_folder = 'Run-'; % this is the folder naming convention for any run-based subfolders within 'funct_folder'
% If no run sub-folders, leave this empty ('')
% Do not need wild card at end here (*)

% These prefixes enable the script to automatically locate the appropriate
%   functional and anatomical files within the specified 
%   'funct_folder'/'run_folder' and 'anat_folder' listed above. 
funct_prefix = '*_WIP_NTD*'; 
anat_prefix = '*_WIP_MPRAGE_*'; 
% Note: if all anatomicals and/or functionals are named identically across
%   subjects, 'funct_prefix' and/or 'anat_prefix' can be listed
%   as their full file names (e.g., 'MPRAGE_anat.nii'); however, if they
%   have characteristics specific to each subject, such as a subject
%   number, list only the prefix with a wildcard '*' (asterisk) at the end,
%   denoting that the rest of the file name can vary beyond this point. If
%   subject-specific naming conventions come before the part shared across
%   subjects, a wildcard * can be used at the start of the name as well, or
%   at both the start and at the end (e.g., '*_WIP_MPRAGE_*').

% List the suffixes to be appended to the appropriate BIDS renamed files:
anat_suffix = 'T1w'; % other options: 'T2w', 'anat'
funct_suffix = 'bold'; 

% List number of subjects and runs:
nSubjects = 3;
nRuns = 3;

suppress_warnings = false;
verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: you can simply move the cursor to this section and click 'Run
% Section' in the Editor Tab above (or press Ctrl + Enter)

orig_dir = pwd;

% First obtain sorted list of subject folders:
subj_listing = dir(fullfile(main_input_folder,[subject_folder_prefix,'*']));
removeIdx = false(length(subj_listing),1);
for i = 1:length(subj_listing)
    if ~subj_listing(i).isdir || strcmp(subj_listing(i).name,'.') || strcmp(subj_listing(i).name,'..')
        removeIdx(i) = true;
    end
    if ~isempty(exclude_prefix)
       check = strfind(subj_listing(i).name,exclude_prefix);
       if ~isempty(check)
           removeIdx(i) = true;
       end
    end
end; subj_listing(removeIdx) = [];

% Get # sub
nSubjects_calc = length(subj_listing);
if verbose
    disp(['Found ',num2str(nSubjects_calc),' subject folders.']);
end
if ~suppress_warnings && (nSubjects ~= nSubjects_calc)
    warning(['Found ',num2str(nSubjects_calc),' subject folders but ',num2str(nSubjects),' were specified.']);
end

% Create BIDS output folder if doesn't exist:
if ~isdir(BIDS_output_folder); mkdir(BIDS_output_folder); end

% Create BIDS output subject folders:
subdir = {'funct','anat','physio','qcreport'};
subjDirs = cell(1,nSubjects_calc);
for i = 1:nSubjects_calc
    if nSubjects_calc<100
        subjDirs{i} = fullfile(BIDS_output_folder,sprintf('sub-%02g',i));
        if ~isdir(subjDirs{i}); mkdir(subjDirs{i}); end
        nzero = 2;
    elseif nSubjects_calc<1000
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

% Find Anatomicals: 
structFiles = cell(1,nSubjects_calc); % full paths to original anats
for i = 1:nSubjects_calc
    % Find Anat:
    inpath = fullfile(main_input_folder,subj_listing(i).name,anat_folder);
    listing = dir(fullfile(inpath,anat_prefix));
    if ~suppress_warnings
        if isempty(listing) 
            warning(['Anatomical for subject ',num2str(i),' not found.']);
        elseif length(listing)>1
            warning([num2str(length(listing)),' anatomicals found for subject ',num2str(i),'; using first.']);
        end
    end
    structFiles{i} = fullfile(inpath,listing(1).name);
end

% Write out list of subject-file correspondences:
write_data = cell(nSubjects_calc,1);
for i = 1:nSubjects_calc
    write_data{i} = ['Anatomical ',sprintf(['%0',num2str(nzero),'g:  '],i),structFiles{i}];
end
write_data = char(write_data);
fileID = fopen(fullfile(BIDS_output_folder,'file_correspondence_info.txt'),'w');
fprintf(fileID,'%s\r\n','-------------------------------- Anatomicals ---------------------------------');
for i = 1:nSubjects_calc
    fprintf(fileID,'%s\r\n',write_data(i,:));
end
fprintf(fileID,'%s\r\n',' ');
fclose(fileID);

% Copy over structs and rename:
if verbose; disp('Copying over anatomicals...'); end
[~,~,ext] = fileparts(structFiles{1});
anats = cell(1,nSubjects_calc);
hWait = waitbar(0,'Copying anatomicals...');
for i = 1:nSubjects_calc
    waitbar(i/nSubjects_calc,hWait);
    anats{i} = fullfile(subjDirs{i},'anat',sprintf(['sub-%0',num2str(nzero),'g_',anat_suffix,ext],i));
    copyfile(structFiles{i},anats{i})
end; close(hWait);
if verbose; disp('Done.'); end

% Find Functionals:
functFiles = cell(1,nSubjects_calc); % full paths to original anats
nRuns_calc = zeros(1,nSubjects_calc); % counts # runs per subject
for i = 1:nSubjects_calc
    inpath = fullfile(main_input_folder,subj_listing(i).name,funct_folder);
    if isempty(run_folder)
        listing = dir(fullfile(inpath,funct_prefix));
    else
        listing = dir(fullfile(inpath,[run_folder,'*']));
        removeIdx = false(length(listing),1);
        for j = 1:length(listing)
            if ~listing(j).isdir || strcmp(listing(j).name,'.') || strcmp(listing(j).name,'..')
                removeIdx(j) = true;
            end
        end
        listing(removeIdx) = [];
        funct_run_folders = cell(1,length(listing));
        fileNames = cell(1,length(listing));
        for j = 1:length(listing)
            funct_run_folders{j} = listing(j).name;
            file_listing = dir(fullfile(inpath,funct_run_folders{j},funct_prefix));
            if length(file_listing) > 1
                error(['ERROR: More than one functional found in ',fullfile(inpath,funct_run_folders{j})])
            elseif isempty(file_listing)
                warning(['No functional file found in ',fullfile(inpath,funct_run_folders{j})])
                continue;
                fileNames{j} = '';
            end
            fileNames{j} = file_listing(1).name;
        end
    end
    nRuns_calc(i) = length(listing);
    if ~suppress_warnings
        if isempty(listing) 
            warning(['No functional runs found for subject ',num2str(i),'.']);
        elseif nRuns_calc(i)~=nRuns
            warning([num2str(nRuns_calc(i)),' functional runs found for subject ',num2str(i),', but ',num2str(nRuns),' runs were specified.']);
        end
    end
    functFiles{i} = cell(1,nRuns_calc(i));
    for j = 1:nRuns_calc(i)
        if isempty(run_folder)
            functFiles{i}{j} = fullfile(inpath,listing(j).name);
        else
            if ~isempty(fileNames{j})
                functFiles{i}{j} = fullfile(inpath,funct_run_folders{j},fileNames{j});
            else
                functFiles{i}{j} = [];
            end
        end
    end
end

% Write out list of subject-file correspondences:
write_data = {};
count = 0;
for i = 1:nSubjects_calc
    for j = 1:length(functFiles{i})
        count = count+1;
        if ~isempty(functFiles{i}{j})
            write_data{count} = ['Functional, subject ',sprintf(['%0',num2str(nzero),'g, run %02g:   '],[i,j]),functFiles{i}{j}]; %#ok
        else
            write_data{count} = ['Functional, subject ',sprintf(['%0',num2str(nzero),'g, run %02g:   '],[i,j]),'------------------MISSING------------------']; %#ok
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
if verbose; disp('Copying over functionals...'); end
for i = 1:nSubjects_calc
    for j = 1:nRuns_calc(1)
        try
            fname = functFiles{1}{1};
            [~,~,ext] = fileparts(fname);
            break;
        catch
            ext = '.nii';
            continue;
        end
    end
end
functs = cell(1,nSubjects_calc);
hWait = waitbar(0,'Copying functionals...');
for i = 1:nSubjects_calc
    waitbar(i/nSubjects_calc,hWait);
    for j = 1:length(functFiles{i})
        if ~isempty(functFiles{i}{j})
            functs{i}{j} = fullfile(subjDirs{i},'funct',sprintf(['sub-%0',num2str(nzero),'g_run-%02g_',funct_suffix,ext],[i,j]));
            copyfile(functFiles{i}{j},functs{i}{j})
        else
            functs{i}{j} = 'MISSING';
        end
    end
end; close(hWait);
if verbose; disp('Done.'); end

cd(orig_dir)
if verbose; disp('All file transfers completed successfully!'); end