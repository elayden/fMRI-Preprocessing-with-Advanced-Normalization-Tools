addpath(genpath('/project2/bermanm/NTD_preprocessed'))

quality_report_funct_linux(...
    'parent_folder','/project2/bermanm/NTD_preprocessed',...
    'output_dir','/project2/bermanm/NTD_preprocessed/QC',...
    'funct_prefix','w_',...
    'funct_ext','.nii',...
    'qc_orientation',[0,1,1],...
    'slice_by_time',[0,1,0],...
    'printRes',100)


parsed.parent_folder = '/project2/bermanm/NTD_preprocessed';
parsed.output_dir = 'qc';
parsed.funct_prefix = 'w_';
parsed.funct_ext = '.nii';
parsed.qc_orientation = [0,1,1];
parsed.slice_by_time = [0,1,0];
parsed.printRes = 100;


% Transfer to an outer directory and just keep sagittal for efficient
% look-through:

parentDir = '/project2/bermanm/NTD_preprocessed';
outDir = '/project2/bermanm/NTD_preprocessed/QC';
if ~isdir(outDir); mkdir(outDir); end

for i = 1:46
    fpath = fullfile(parentDir,sprintf('sub-%02g',i),'qcreport','funcProc','qc');
    listing = dir(fullfile(fpath,'*.png'));
    for j = 1:length(listing)
        copyfile(fullfile(fpath,listing(j).name),fullfile(outDir,listing(j).name))
    end
end

listing = dir(fullfile(outDir,'*Coronal*'));
for i = 1:length(listing)
    delete(fullfile(outDir,listing(i).name))
end

% Check weird image:
img1 = load_nii(fullfile(parentDir,'sub-02','funct','sub-02_run-04_bold.nii'));
img2 = load_nii('/project2/bermanm/NTD_raw/NTD_1003_07202017/Nifti/NTD_1003_WIP_NTD4_SENSE_6_1.nii')

img2 = load_nii('/project2/bermanm/NTD_raw/NTD_1003_07202017/Nifti/NTD_1003_WIP_NTD3_SENSE_5_1.nii')
img2 = load_nii('/project2/bermanm/NTD_raw/NTD_1003_07202017/Nifti/NTD_1003_WIP_NTD5_SENSE_7_1.nii')

corr(img1.img(:),img2.img(:))


img1 = load_nii(fullfile(parentDir,'sub-02','funct','sub-02_run-08_bold.nii'));
img2 = load_nii('/project2/bermanm/NTD_raw/NTD_1003_07202017/Nifti/NTD_1003_WIP_NTD7_SENSE_10_1.nii')

corr(img1.img(:),img2.img(:))


% Rerun Sub2

fpath = '/project2/bermanm/rerun_02/'; mkdir(fpath);
inpath = '/project2/bermanm/NTD_raw/NTD_1003_07202017/Nifti/';
listing = dir(fullfile(inpath,'NTD_1003_WIP_NTD*'));
listing(6) = [];

functs = cell(1,1);
for i = 1:length(listing)
    functs{1}{i} = fullfile(fpath,sprintf('sub-02_run-%02g_bold.nii',i));
%     copyfile(fullfile(inpath,listing(i).name),functs{i})
end


addpath(genpath('/project2/bermanm/ants_preprocessing_automated'))
nullFuncts(functs); % This actually still must be done, even w/ human data


% Rerun QC for sub-02
quality_report_funct_linux(...
    'parent_folder','/project2/bermanm/rerun_02',...
    'output_dir','/project2/bermanm/NTD_preprocessed/QC',...
    'funct_prefix','w_',...
    'funct_ext','.nii',...
    'qc_orientation',[0,1,1],...
    'slice_by_time',[0,1,0],...
    'printRes',100)
