function smooth_funcs(subject_folder)

    % Define main output folder:
    output_folder = '/project2/bermanm/NUBE_data/procdata';

    % Define subfolders of subject_folder where anats and funcs are located:
    func_folder = fullfile(subject_folder, 'func');

    % Define wildcard pattern to locate anats and funcs:
    func_pattern = 'sub-*_bold.nii'; % should be original name, not warped
    
    % Options:
    smoothing_FWHM = [6, 6, 6];

    % File prefixes for output:
    warped = 'w_';
    smoothed = 's_';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add paths (probably better to skip this, just add paths initially)
    addpath(genpath('/project2/bermanm/ants_preprocessing_automated'))

    % Define ANTs path:
    antspath = '/software/ANTs-2.1-el6-x86_64/bin/';
    setenv('ANTSPATH', antspath) % set as environment variable (required for ANTs bash scripts)

    % Check for valid input:
    if ~isdir(subject_folder); error('Aborting:  could not locate input ''subject_folder.'''); end
    if ~isdir(output_folder); mkdir(output_folder); end

    % Get subject name/number:
    subject_num = regexp(subject_folder, filesep, 'split');
    subject_num = subject_num(~cellfun('isempty', subject_num));
    output_folder = fullfile(output_folder, subject_num{end});

    % Find anat & func files:
    func_list = dir(fullfile(func_folder, func_pattern));
    
    %% Smoothing (if requested):

    smoothing_sigma = smoothing_FWHM ./ sqrt( 8 * log(2)); % convert from FWHM

    for i = 1:length(func_list)
        in_fname = fullfile(output_folder, 'func', [warped, func_list(i).name]);
        out_fname = fullfile(output_folder, 'func', [smoothed, warped, func_list(i).name]);

        command = [fullfile(antspath,'SmoothImage'),' 4 ', in_fname, ...
            ' ', num2str(smoothing_sigma(1)),'x',num2str(smoothing_sigma(2)),...
            'x',num2str(smoothing_sigma(3)),'x0', ...
            ' ', out_fname, ' 1'];

        % Run command:
        [status, cmdout] = system(command,'-echo'); %#ok
    end
        
end