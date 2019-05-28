% ants_preprocess_NUBE.m

dataPath = '/project2/bermanm/NUBE_data/rawdata/';

% Add the path for ANTs scripts:
addpath(genpath('/project2/bermanm/ants_preprocessing_automated/'))
% cd('/project2/bermanm/ants_preprocessing_automated/')

% Add the dataPath and get the subject folders (check this):
addpath(genpath(dataPath))
subFolders = dir(dataPath);
subFolders = subFolders(3:end); % avoid '.', '..'

% Some important notes on parallel processing: 
% 1. ANTs procedures in 'preprocess_subject.m', as called below, should
%    only be run on compute nodes, not on the login node (i.e., you must
%    request resources using either sbatch or sinteractive commands in the
%    Midway terminal). For this script, I recommend starting an
%    sinteractive session using, e.g., the following series of commands:
% 
%    sinteractive --partition=broadwl --nodes=1 --exclusive --time=6:00:00
%    module load matlab
%    matlab
% 
%    At this point, you would be able to navigate to this script and run
%    the parfor loops below in parallel on a single broadwl node. To open
%    multiple sinteractive sessions at once, simply open more Midway
%    terminals and repeat the above steps. The sinteractive session above
%    requests exclusive use of one broadwl node with 28 (i.e., 14) cores
%    for 6 hours.
% 2. can only start pool w/ half the num cores supposedly available on a node (half are logical cores)
%    -e.g., broadwl has 28 cores, but it appears that only 14 can be
%    requested for a parallel pool
% 3. according to the Midway documentation, parallel Matlab jobs are
%    currently limited to running on just one node
% 4. for an unknown reason, ANTs encounters file writing conflicts if
%    parfor loops are run on more than ~7 cores at a time from the parallel
%    pool (possibly memory-related issues, but remains unclear); this is
%    why the processing below is split up into several batchs of parfor
%    loops with only 7 iterations each. This is obviously not ideal, but
%    additional parallel processing can be utilized by starting several
%    sinteractive sessions on different nodes (1 node each) and running
%    separate sets of batches of 7.

% Initialize parallel pool:
pool = parpool('local', str2double(getenv('SLURM_CPUS_ON_NODE'))/2); 

% Run preprocess_subject:
parfor i = 1:7 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 8:14 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 15:21 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 22:28 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 29:35 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 36:42 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 43:49 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 50:56 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 57:63 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 64:65 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

parfor i = 66:69 %length(subFolders)
    preprocess_subject(fullfile(dataPath, subFolders(i).name));
end

%delete(pool)
