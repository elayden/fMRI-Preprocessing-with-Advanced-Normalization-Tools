function preprocess_subject(subject_folder)
% Note: subject_folder must have leading '/' and use forward slashes:
%  e.g., subject_folder = '/project2/bermanm/NUBE_data/rawdata/sub-1001/'

    % Note: currently does not do skullstripping (could add ANTs program
    % for this)
    
    % Define main output folder:
    output_folder = '/project2/bermanm/NUBE_data/procdata';
    dataset_name = 'NUBE';

    % Define subfolders of subject_folder where anats and funcs are located:
    anat_folder = fullfile(subject_folder, 'anat');
    func_folder = fullfile(subject_folder, 'func');
    
    % Locate MNI template:
    template_path = '/project2/bermanm/ants_preprocessing_automated/MNI_1mm_Template';
    anat_template = fullfile(template_path, 'mni_temp.nii');

    % Define wildcard pattern to locate anats and funcs:
    anat_pattern = 'sub-*T1w.nii';
    func_pattern = 'sub-*_bold.nii';
    
    % Options:
    perform_tissue_segmentation = true;
    
    % Smoothing:
    perform_smoothing = 1; % (0) none, (1) ANTs, (2) AFNI's 3dBlurInMask [only isotropic allowed, will use FWHMx]
    smoothing_FWHM = [6, 6, 6]; % [FWHMx, FWHMy, FWHMz]
    % If using AFNI, must specify the path to a mask file to smooth within (in normalized / template space):
    maskFile = fullfile(template_path,'mask.nii');
    
    % Slice Timing Correction:
    perform_slice_timing = 2; % (0) none, (1) SPM, (2) ANTs
    % used by both:
    TR = 2;
    % used by ANTs:
    nSlices = 15; 
    sliceTime = TR/nSlices;
    % used by SPM:
    TA = TR - (TR / nSlices);
    sliceOrder = [1 3 5 7 9 11 13 15 2 4 6 8 10 12 14];
    refSlice = 15;
    nTime = 180;
    
    cleanup_files = true; % delete large displacement fields for warping images after use
    save_intermediates = true; % coreg, moco warps, etc.

    % File prefixes for output:
    bias_corrected = 'b_';
    motion_corrected = 'moco_';
    averaged = 'avg_';
    coregistered = 'coreg_';
    warped = 'w_';
    smoothed = 's_';
    tissue_segmented = 'seg_';
    realigned = 'rp_'; % prefix for realignment parameters (motion correction)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add paths (probably better to skip this, just add paths initially)
    addpath(genpath('/project2/bermanm/ants_preprocessing_automated'))
    %addpath(genpath('/project2/bermanm/NUBE_data'))

    % Define ANTs path:
    antspath = '/software/ANTs-2.1-el6-x86_64/bin/';
    setenv('ANTSPATH', antspath) % set as environment variable (required for ANTs bash scripts)

    % Check for valid input:
    if ~isdir(subject_folder); error('Aborting:  could not locate input ''subject_folder.'''); end
    if ~isdir(output_folder); mkdir(output_folder); end

    % Get subject name/number:
    template_out = fullfile(output_folder, 'templates');
    if ~isdir(template_out); mkdir(template_out); end
    subject_num = regexp(subject_folder, filesep, 'split');
    subject_num = subject_num(~cellfun('isempty', subject_num));
    output_folder = fullfile(output_folder, subject_num{end});
    if ~isdir(output_folder); mkdir(output_folder); end
    if ~isdir(fullfile(output_folder,'func')); mkdir(fullfile(output_folder,'func')); end
    if ~isdir(fullfile(output_folder,'anat')); mkdir(fullfile(output_folder,'anat')); end

    % Find anat & func files:
    anat_list = dir(fullfile(anat_folder, anat_pattern));
    func_list = dir(fullfile(func_folder, func_pattern));

    [~, anat_fname] = fileparts(anat_list(1).name);
    anat_fname = [anat_fname, '_'];
    
    %% Downsample MNI template:
    
    downsamp_fname = fullfile(template_out, [dataset_name, '_rmni_temp.nii']);
    
    % Check whether dataset-specific downsampled template is available:
    if ~exist(downsamp_fname,'file')
        % Interpolate mni_temp.nii to funct dimensions:
        func = load_untouch_nii(fullfile(func_folder, func_list(1).name));
        downsample_template(fullfile(template_path,'mni_temp.nii'), downsamp_fname, ...
            func.hdr.dime.pixdim(2:4), [nan nan nan; nan nan nan])
    end
    
    % Check whether affine transformation of mni template -> downsampled template is available:
    if ~exist(fullfile(template_out, [dataset_name,'_0GenericAffine.mat']), 'file')
        % Create command
        command = [fullfile(antspath,'antsRegistrationSyN.sh'),' -d 3 ', ...
            '-f ', downsamp_fname, ...
            ' -m ', fullfile(template_path, 'mni_temp.nii'), ...
            ' -t r -o ', fullfile(template_out, [dataset_name,'_'])];

        % Run command
        [status, cmdout] = system(command,'-echo'); %#ok
    end
    
    %% N4 bias correction:
    
    % Transfer funcs to output folder (avoid altering original data):
    for i = 1:length(func_list)
        copyfile(fullfile(func_folder, func_list(i).name), fullfile(output_folder, 'func', func_list(i).name))
        nullFunct(fullfile(output_folder, 'func', func_list(i).name))
    end
    
    % Transfer anats to output folder (avoid altering original data):
    for i = 1:length(anat_list)
        copyfile(fullfile(anat_folder, anat_list(i).name), fullfile(output_folder, 'anat', anat_list(i).name))
        nullStruct(fullfile(output_folder, 'anat', anat_list(i).name))
    end
    
    % Functional bias correction:
    for i = 1:length(func_list)
        
        % Create command
        command = [fullfile(antspath,'N4BiasFieldCorrection'),' -d 4 -i ', ...
            fullfile(output_folder, 'func', func_list(i).name), ...
            ' -o ', fullfile(output_folder, 'func', [bias_corrected, func_list(i).name])];

        % Run command
        [status, cmdout] = system(command,'-echo'); %#ok
        
    end  
   
    % Structural bias correction:
    for i = 1:length(anat_list)
        % Create command
        command = [fullfile(antspath,'N4BiasFieldCorrection'),' -d 3 -i ', ...
            fullfile(output_folder, 'anat', anat_list(i).name), ...
            ' -o ', fullfile(output_folder, 'anat', [bias_corrected, anat_list(i).name])];

        % Run command
        [status, cmdout] = system(command,'-echo'); %#ok
    end
 
    %% Functional Slice-Timing Correction:

    if perform_slice_timing==1 % SPM
        spm('defaults','fmri');
        spm_jobman('initcfg');
        flist = cell(length(func_list),1);
        for i = 1:length(func_list)
            for j = 1:nTime
                flist{i}{j} = [fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]),',',num2str(j)];
            end
            flist{i} = flist{i}';
        end
        matlabbatch = {};
        matlabbatch{1}.spm.temporal.st.scans = flist; % cell of cells (scan sessions / subjects)
        matlabbatch{1}.spm.temporal.st.nslices = nSlices;
        matlabbatch{1}.spm.temporal.st.tr = TR;
        matlabbatch{1}.spm.temporal.st.ta = TA;
        matlabbatch{1}.spm.temporal.st.so = sliceOrder;
        matlabbatch{1}.spm.temporal.st.refslice = refSlice;
        matlabbatch{1}.spm.temporal.st.prefix = 'a_';
        spm_jobman('run', matlabbatch);
        % remove prefix:
        for i = 1:length(func_list)
            movefile(fullfile(output_folder, 'func', ['a_',bias_corrected, func_list(i).name]),...
                fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]))
            img = load_untouch_nii(fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]));
            img.hdr.dime.pixdim(5) = TR;
            save_untouch_nii(img, fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]));
        end
    elseif perform_slice_timing==2 % ANTs
        for i = 1:length(func_list)
            command = [fullfile(antspath,'ImageMath'),' 4 ', ...
                fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]), ...
                ' SliceTimingCorrection ', fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]),...
                ' ',num2str(sliceTime),' bspline']; % sliceTiming = TR / nSlices 

            % Run command
            [status, cmdout] = system(command,'-echo'); %#ok
        end
    end
    
    %% Functional motion correction:
    
    % Temporally average functionals:
    for i = 1:length(func_list)
        % Create command
        command = [fullfile(antspath,'antsMotionCorr'),' -d 3 -a ', ...
            fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]), ...
            ' -o ', fullfile(output_folder, 'func', [averaged, bias_corrected, func_list(i).name])];

        % Run command
        [status, cmdout] = system(command,'-echo'); %#ok
    end

	% Perform affine motion correction:
    for i = 1:length(func_list)
        [fpath, fname] = fileparts(fullfile(output_folder, 'func', ...
            [motion_corrected, bias_corrected, func_list(i).name]));
        fname = [fname, '_']; %#ok
        
        % Create command
        command = [fullfile(antspath,'antsMotionCorr'),' -d 3 -o ', ...
            '[',fullfile(fpath, fname), ', ', ...
            fullfile(output_folder, 'func', [motion_corrected, bias_corrected, func_list(i).name]), ', ', ...
            fullfile(output_folder, 'func', [averaged, motion_corrected, bias_corrected, func_list(i).name]),']', ' ', ...
            ' -m ', 'MI[', fullfile(output_folder, 'func', [averaged, bias_corrected, func_list(i).name]), ', ', ...
            fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]), ', 1, 32, Regular, 0.2]', ...
            ' -t Rigid[0.05] -i 20x10 -u 1 -e 1 -s 1x0 -f 2x1 -n 10 -w 1'];

        % Run command
        [status, cmdout] = system(command,'-echo'); %#ok
    end
    
    %% Co-register anat and avg moco func:
    for i = 1:length(func_list)
        [~, fname] = fileparts(fullfile(output_folder, 'func', ...
            [coregistered, func_list(i).name]));
        fname = [fname, '_']; %#ok
        
        % Create command
        command = [fullfile(antspath,'antsRegistrationSyN.sh'),' -d 3 -n 6 ', ...
            '-f ', fullfile(output_folder, 'anat', [bias_corrected, anat_list(1).name]), ...
            ' -m ', fullfile(output_folder, 'func', [averaged, motion_corrected, bias_corrected, func_list(i).name]), ...
            ' -t a -o ', fullfile(output_folder, 'func', fname)];
        
        % Run command
        [status, cmdout] = system(command,'-echo'); %#ok
    end
      
    %% Warp anat to template (ICBM-152):
    for i = 1:length(anat_list)
        
        % Create command
        command = [fullfile(antspath,'antsRegistrationSyN.sh'),' -d 3 -f ', ...
            anat_template, ' -m ', fullfile(output_folder, 'anat', [bias_corrected, anat_list(i).name]), ...
            ' -n 1 -t s -o ', fullfile(output_folder, 'anat', [warped, anat_fname])];

        % Run command
        [status, cmdout] = system(command,'-echo'); %#ok
    end

    %% Combine/collapse warps, apply to funcs:
    for i = 1:length(func_list)
        [~, fname] = fileparts(fullfile(output_folder, 'func', ...
            [func_list(i).name]));
        fname = [fname, '_']; %#ok
        
        % Preliminary: get TR & # time points from a functional:
        img = load_untouch_nii(fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]));
        timepoints = img.hdr.dime.dim(5);
        tr = img.hdr.dime.pixdim(5);
        
        % 1. Collapse Funct2Template Transforms
        command = [fullfile(antspath,'antsApplyTransforms'),' -d 3 -o [', ...
            fullfile(output_folder, 'func', ['CollapsedWarp_',fname,'.nii.gz']),',1] ', ...
            ' -t ', fullfile(template_out, [dataset_name, '_0GenericAffine.mat']), ...
            ' -t ', fullfile(output_folder, 'anat', [warped, anat_fname, '1Warp.nii.gz']), ...
            ' -t ', fullfile(output_folder, 'anat', [warped, anat_fname, '0GenericAffine.mat']), ...
            ' -t ', fullfile(output_folder, 'func', [coregistered, fname, '0GenericAffine.mat']), ...
            ' -r ', downsamp_fname];
        
        % Run command:
        [status, cmdout] = system(command,'-echo'); %#ok

        % 2. Replicate these warps to 4D along with template
        command = [fullfile(antspath,'ImageMath'),' 3 ',...
            fullfile(output_folder, 'func', ['rep_CollapsedWarp_',fname,'.nii.gz']), ...
            ' ReplicateDisplacement ', ...
            fullfile(output_folder, 'func', ['CollapsedWarp_',fname,'.nii.gz']), ...
            ' ', num2str(timepoints),' ', num2str(round(tr)),' 0'];
        
        % Run command:
        [status, cmdout] = system(command,'-echo'); %#ok
        
        % 3. Replicate downsampled MNI template:
        command = [fullfile(antspath,'ImageMath'),' 3 ',...
            fullfile(output_folder, 'func', ['rep_', dataset_name, '_rmni_temp.nii']), ...
            ' ReplicateImage ', downsamp_fname, ...
            ' ', num2str(timepoints),' ', num2str(round(tr)),' 0'];
        
        % Run command:
        [status, cmdout] = system(command,'-echo'); %#ok
        
        % 4. Apply warps to original funcs:
        command = [fullfile(antspath,'antsApplyTransforms'),' -d 4 -o ', ...
            fullfile(output_folder, 'func', [warped, func_list(i).name]),...
            ' -t ', fullfile(output_folder, 'func', ['rep_CollapsedWarp_',fname,'.nii.gz']), ...
            ' -t ', fullfile(output_folder, 'func', [motion_corrected, bias_corrected, fname, 'Warp.nii.gz']), ...
            ' -r ', fullfile(output_folder, 'func', ['rep_', dataset_name, '_rmni_temp.nii']), ...
            ' -i ', fullfile(output_folder, 'func', [bias_corrected, func_list(i).name])];
        
        % Run command:
        [status, cmdout] = system(command,'-echo'); %#ok
        
    end
    
    %% Smoothing (if requested):
    if perform_smoothing==1 % ANTs
        
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
        
    elseif perform_smoothing==2 % AFNI 3dBlurInMask
        afniPath = '/software/afni-16.2-x86_64/bin/';
        system('module load afni'); 
%         [status, cmdout] = system([afniPath,'3dBlurInMask -help'],'-echo'); 

        for i = 1:length(func_list)
            in_fname = fullfile(output_folder, 'func', [warped, func_list(i).name]);
            [path, fname, ext] = fileparts(in_fname);
            out_fname = fullfile(path, [smoothed, fname, ext]);
            if exist(out_fname,'file'); delete(out_fname); end
            command = [afniPath,'3dBlurInMask -input ',in_fname,' -FWHM ',num2str(smoothing_FWHM(1)),' -mask ',maskFile,' -prefix ',out_fname];

            % Run command:
            [status, cmdout] = system(command,'-echo'); %#ok
        end
        
    end
    
    %% Cleanup step (delete large warps):
    if cleanup_files
        for i = 1:length(func_list)
            [~, fname] = fileparts(fullfile(output_folder, 'func', ...
                [func_list(i).name]));
            fname = [fname, '_']; %#ok
            delete(fullfile(output_folder, 'func', ['CollapsedWarp_',fname,'.nii.gz']))
            delete(fullfile(output_folder, 'func', ['rep_CollapsedWarp_',fname,'.nii.gz']))
            delete(fullfile(output_folder, 'func', func_list(i).name))
        end
        delete(fullfile(output_folder, 'func', ['rep_', dataset_name, '_rmni_temp.nii']))
        delete(fullfile(output_folder, 'anat', anat_list(1).name))
    end
    
    if ~save_intermediates
        for i = 1:length(func_list)
            [~, fname] = fileparts(fullfile(output_folder, 'func', ...
                [func_list(i).name]));
            fname = [fname, '_']; %#ok
            delete(fullfile(output_folder, 'func', [bias_corrected, func_list(i).name]))
            delete(fullfile(output_folder, 'func', [motion_corrected, bias_corrected, func_list(i).name]))
            delete(fullfile(output_folder, 'func', [motion_corrected, bias_corrected, fname, 'Warp.nii.gz']))
            delete(fullfile(output_folder, 'func', [coregistered, fname, '0GenericAffine.mat']))
            delete(fullfile(output_folder, 'func', [coregistered, fname, 'Warped.nii.gz']))
            delete(fullfile(output_folder, 'func', [averaged, bias_corrected, func_list(i).name]))
            delete(fullfile(output_folder, 'func', [averaged, motion_corrected, bias_corrected, func_list(i).name]))
        end
        
        delete(fullfile(output_folder, 'anat', [bias_corrected, anat_list(1).name]))
        delete(fullfile(output_folder, 'anat', [warped, anat_fname, '1InverseWarp.nii.gz']))
    end
    
    %% Convert ANTs motion parameters to 6 rigid-body motion params (as in SPM):
    for i = 1:length(func_list)
        [~, fname] = fileparts(fullfile(output_folder, 'func', ...
            [func_list(i).name]));
        rigidBodyParams(fullfile(output_folder,'func',[motion_corrected, bias_corrected, fname, '_MOCOparams.csv']),...
            fullfile(output_folder,'func',[realigned, fname, '.txt']));
    end
    
    %% Tissue segmentation
    if perform_tissue_segmentation
        command = [fullfile(antspath,'Atropos'),' -d 3', ...
            ' -a ', fullfile(output_folder, 'anat', [warped, anat_fname, 'Warped.nii.gz']), ...
            ' -i PriorProbabilityImages[6, ', ...
            fullfile(template_path, 'mni_tissue_%02d.nii'),', .6]', ...
            ' -x ', fullfile(template_path, 'mni_mask.nii'), ...
            ' -o [', fullfile(output_folder, 'anat', [tissue_segmented, anat_list(1).name]), ', ',...
            fullfile(output_folder, 'anat', ['tissueProb_%02d_',anat_list(1).name]),']'];
        
        % Run command:
        [status, cmdout] = system(command,'-echo'); %#ok
        
    end

    %% Helper function for pausing execution (not needed):
    function waitFor(status, timeIncr, maxIncr) %#ok
        for ix = 1:maxIncr
            if ~status
                break;
            else
                pause(timeIncr)
            end
        end
    end
    
end