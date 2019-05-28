function [analysis_dir, anats, functs] = organize4ants_point_click
% % [ANALYSIS_DIR, ANATS, FUNCTS] = ORGANIZE_NEUROIMAGING_FILES;
% 
% INFO: uses the naming conventions of the Brain Imaging Data Structure 
% (BIDS) format <http://bids.neuroimaging.io/>. However, all files are 
% placed into a single director to enable preprocessing using ANTs. Users 
% are directed through a series of file explorer menus and dialogue
% boxes which allow automated point-and-click transfer and renaming of 
% files into the director. Additionally, a 'file_correspondence_info.txt' 
% file is automatically written to this directory, allowing users to verify 
% that the correct files were selected, renamed, and copied.
% 
% Note: this function does not delete any of the original files.
% 
% OUTPUTS,
%   analysis_dir,   string: the full path to the main directory
%   anats,      a cell array wherein each cell contains a string denoting
%               the full filepath to each anatomical image in the BIDS 
%               directories 
%   functs,     a cell array wherein each cell contains a string denoting
%               the full filepath to each functional image in the BIDS 
%               directories; the structure of this cell array will vary
%               depending on the study specifics (i.e., number of runs,
%               etc.)
% 
% Author:       Elliot Layden, University of Chicago, August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    orig_dir = pwd;
    
    % Identify Function Path and Add Helper Scripts:
    script_fullpath = mfilename('fullpath');
    [script_path,~,~] = fileparts(script_fullpath);
    addpath(genpath(script_path))

    % Select main analysis folder:
    dialog_title = 'Select main analysis folder...';
    analysis_dir = uigetdir(cd,dialog_title); % get directory
    if analysis_dir==0; error('User cancelled preprocessing.'); end

    % Input number of subjects:
    answer = inputdlg('Number of Subjects: ');
    if isempty(answer); error('User cancelled preprocessing.'); end
    nSubjects = str2double(answer);
    if ~isnumeric(nSubjects) && nSubjects>0
        error('Specify a positive integer for number of subjects.')
    end

    % Determine significant digits for numbering files:
    if nSubjects<100
        nzero = 2;
    elseif nSubjects<1000
        nzero = 3;
    else
        nzero = 4;
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
    fileID = fopen(fullfile(analysis_dir,'file_correspondence_info.txt'),'w');
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
        anats{i} = fullfile(analysis_dir,sprintf(['sub-%0',num2str(nzero),'g_anat',ext],i));
        copyfile(structFiles{i},anats{i})
    end
    
    % Null Structs:
    nullStructs(anats);

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
                if ischar(functName)
                    functName={functName};
                end
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
    fileID = fopen(fullfile(analysis_dir,'file_correspondence_info.txt'),'a');
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
                functs{i}{j} = fullfile(analysis_dir,sprintf(['sub-%0',num2str(nzero),'g_run-%02g_bold',ext],[i,j]));
                copyfile(functFiles{i}{j},functs{i}{j})
            else
                for k = 1:length(functFiles{i}{j})
                    functs{i}{j}{k} = fullfile(analysis_dir,sprintf(['sub-%0',num2str(nzero),'g_run-%02g_bold_%04g',ext],[i,j,k]));
                    copyfile(functFiles{i}{j}{k},functs{i}{j}{k})
                end
            end
        end
    end; close(hWait);
    
    display('Adjusting functional header information (turn off qform), please wait (could take several minutes)...')
    nullFuncts(functs);
    
    %% Add bash scripts to analysis_dir:
    copyfile(fullfile(script_path,'custom_bash_utilities','ants_preprocessing.sh'),fullfile(analysis_dir,'ants_preprocessing.sh'))
    copyfile(fullfile(script_path,'custom_bash_utilities','ants_preprocessing_midway1.sbatch'),fullfile(analysis_dir,'ants_preprocessing_midway1.sbatch'))
    
    %% Save a study-specific functional dimension mni_temp.nii:
    mni_temp = fullfile(script_path,'MNI_1mm_Template','mni_temp.nii');
    copyfile(mni_temp,fullfile(analysis_dir,'mni_temp.nii'))
    for i = 1:6
        copyfile(fullfile(script_path,'MNI_1mm_Template',sprintf('mni_tissue_%02g.nii',i)),fullfile(analysis_dir,sprintf('mni_tissue_%02g.nii',i)));
    end
    copyfile(fullfile(script_path,'MNI_1mm_Template','mni_mask.nii'),fullfile(analysis_dir,'mni_mask.nii'));
    
    try 
        img = load_untouch_nii(functs{1}{1}); % 3.2500 3.2500 3.5000
    catch
        error('ERROR <organize4ants.m line 213>: Could not locate funct 01 run 01.')
    end
    
    % Interpolate mni_temp.nii to funct dimensions:
    resize_img(fullfile(analysis_dir,'mni_temp.nii'), img.hdr.dime.pixdim(2:4), [nan nan nan; nan nan nan])
    
    cd(orig_dir)
    
    disp('FINISHED: organize4ants.m finished with 0 errors.')

    %% Helper Functions
    function nullStructs(anats)
        for ix = 1:length(anats)
            struct = load_untouch_nii(anats{ix});
            struct.hdr.hist.qform_code = 0;
            struct.hdr.hist.sform_code = 1;
            struct.hdr.hist.quatern_b = 0;
            struct.hdr.hist.quatern_c = 0;
            struct.hdr.hist.quatern_d = 0;
            struct.hdr.hist.qoffset_x = 0;
            struct.hdr.hist.qoffset_y = 0;
            struct.hdr.hist.qoffset_z = 0;
            save_untouch_nii(struct,anats{ix}) 
        end
    end

    function nullFuncts(functs)
        for ix = 1:length(functs)
            for iy = 1:length(functs{ix})
                funct = load_untouch_nii(functs{ix}{iy});
                funct.hdr.dime.pixdim(6:8) = 1;
                funct.hdr.hist.qform_code = 0;
                funct.hdr.hist.sform_code = 1;
                funct.hdr.hist.quatern_b = 0;
                funct.hdr.hist.quatern_c = 0;
                funct.hdr.hist.quatern_d = 0;
                funct.hdr.hist.qoffset_x = 0;
                funct.hdr.hist.qoffset_y = 0;
                funct.hdr.hist.qoffset_z = 0;
                save_untouch_nii(funct,functs{ix}{iy})
            end
        end
    end
end