function ants2bids_single(analysis_dir)
% if no analysis_dir, have pop-up menu for use input
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
%                wsub-%02g_anat.nii
%                y_sub-%02g_anat.nii  % deformation field?
%                msub-%02g_anat.nii
%                c1sub-%02g_anat.nii ... c6sub-%02g_anat.nii
%    ...
%    sub-<last>
% 
% Author:       Elliot Layden, University of Chicago, August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get Script Path:
scriptPath = mfilename('fullpath');
[scriptPath,~,~] = fileparts(scriptPath);
addpath(genpath(scriptPath));

% Select main analysis folder:
if nargin<1 || isempty(analysis_dir)
    dialog_title = 'Select main ANTs data folder...';
    analysis_dir = uigetdir(cd,dialog_title); % get directory
    if analysis_dir==0; error('User cancelled preprocessing.'); end
end
             
% Make subject folders:

% First, auto-detect # subjects:
listing = dir(fullfile(analysis_dir,'sub-*_anat.nii'));
nSubjects = length(listing);
sFolders = cell(nSubjects,1);
for i = 1:nSubjects
    sFolders{i} = fullfile(analysis_dir,sprintf('sub-%02g',i));
    if ~isdir(sFolders{i}); mkdir(sFolders{i}); end
    if ~isdir(fullfile(sFolders{i},'anat')); mkdir(fullfile(sFolders{i},'anat')); end
    if ~isdir(fullfile(sFolders{i},'funct')); mkdir(fullfile(sFolders{i},'funct')); end
    if ~isdir(fullfile(sFolders{i},'qcreport')); mkdir(fullfile(sFolders{i},'qcreport')); end
    if ~isdir(fullfile(sFolders{i},'qcreport','anatProc')); mkdir(fullfile(sFolders{i},'qcreport','anatProc')); end
    if ~isdir(fullfile(sFolders{i},'qcreport','funcProc')); mkdir(fullfile(sFolders{i},'qcreport','funcProc')); end
    if ~isdir(fullfile(sFolders{i},'qcreport','funcProc','MoCo')); mkdir(fullfile(sFolders{i},'qcreport','funcProc','MoCo')); end
    try
        transferFile(fullfile(analysis_dir,listing(i).name),fullfile(fullfile(sFolders{i},'anat',listing(i).name)))
    catch
    end
    try
        transferFile(fullfile(analysis_dir,['b_',listing(i).name]),fullfile(fullfile(sFolders{i},'qcreport','anatProc',['b_',listing(i).name])))
    catch
    end
end

% Transfer other structural files:
for i = 1:nSubjects
    
    outPath = fullfile(sFolders{i},'qcreport','anatProc');
    
    % Transfer Normalized Structs
    fnameStruct = fullfile(analysis_dir,sprintf('MNI_sub-%02g_Warped.nii',i));
    if exist(fnameStruct,'file')
        transferFile(fnameStruct,fullfile(outPath,sprintf('w_sub-%02g_T1w.nii',i)))
    end
    if exist([fnameStruct,'.gz'],'file')
        transferFile([fnameStruct,'.gz'],fullfile(outPath,sprintf('w_sub-%02g_T1w.nii.gz',i)))
    end
    
    % Transfer warps
    listing3 = dir(fullfile(analysis_dir,sprintf('MNI_sub-%02g_*',i)));
    for j = 1:length(listing3)
        transferFile(fullfile(analysis_dir,listing3(j).name),fullfile(outPath,listing3(j).name))
    end
    
    % Transfer Tissue Seg Files:
    listing4 = dir(fullfile(analysis_dir,sprintf('labelledTissues_%02g*',i)));
    for j = 1:length(listing4)
        transferFile(fullfile(analysis_dir,listing4(j).name),fullfile(outPath,listing4(j).name))
    end
    
    listingT = dir(fullfile(analysis_dir,sprintf('pTissue_*_s_%02g.nii',i)));
    if numel(listingT)~=0
        for j = 1:length(listingT)
            fname = listingT(j).name;
            transferFile(fullfile(analysis_dir,fname),fullfile(outPath,fname))
        end
    else
        warning('Found no pTissue_* files.')
    end
end

% Transfer functs:
for i = 1:nSubjects
    listing1 = dir(fullfile(analysis_dir,sprintf('sub-%02g_run-*_bold.nii',i)));
    nRuns = length(listing1);
    for j = 1:nRuns
        transferFile(fullfile(analysis_dir,listing1(j).name),fullfile(sFolders{i},'funct',listing1(j).name))
    end
    listing2 = dir(fullfile(analysis_dir,sprintf('w_sub-%02g_run-*_bold.nii',i)));
    for j = 1:length(listing2)
        transferFile(fullfile(analysis_dir,listing2(j).name),fullfile(sFolders{i},'qcreport','funcProc',listing2(j).name))
    end
    listing3 = dir(fullfile(analysis_dir,sprintf('coreg_sub-%02g*',i)));
    for j = 1:length(listing3)
        transferFile(fullfile(analysis_dir,listing3(j).name),fullfile(sFolders{i},'qcreport','funcProc',listing3(j).name))
    end
    listing4 = dir(fullfile(analysis_dir,sprintf('sub-%02g_run-*_bold_avg*',i)));
    for j = 1:length(listing4)
        transferFile(fullfile(analysis_dir,listing4(j).name),fullfile(sFolders{i},'qcreport','funcProc',listing4(j).name))
    end
    listing5 = dir(fullfile(analysis_dir,sprintf('moco_sub-%02g*',i)));
    for j = 1:length(listing5)
        transferFile(fullfile(analysis_dir,listing5(j).name),fullfile(sFolders{i},'qcreport','funcProc','MoCo',listing5(j).name))
    end
    listing6 = dir(fullfile(analysis_dir,sprintf('voxelWiseDisplace_sub-%02g*',i)));
    for j = 1:length(listing6)
        transferFile(fullfile(analysis_dir,listing6(j).name),fullfile(sFolders{i},'qcreport','funcProc','MoCo',listing6(j).name))
    end
    listing7 = dir(fullfile(analysis_dir,sprintf('displaceMap_sub-%02g*',i)));
    for j = 1:length(listing7)
        transferFile(fullfile(analysis_dir,listing7(j).name),fullfile(sFolders{i},'qcreport','funcProc','MoCo',listing7(j).name))
    end
end

% Run conversion of ANTs MOCOparams.csv files to 6 rigid-body motion params:
for i = 1:nSubjects
    currPath = fullfile(sFolders{i},'qcreport','funcProc','MoCo');
    listing = dir(fullfile(currPath,'moco_*MOCOparams.csv'));
    for j = 1:length(listing)
        rigidBodyParams(fullfile(currPath,listing(j).name),fullfile(currPath,sprintf('rp_sub-%02g_run-%02g.txt',[i,j])));
    end
end

% Utilities
    function transferFile(indir,outdir)
        try
            movefile(indir, outdir,'f')
        catch
            warning(['Could not find ',indir,'.'])
        end
    end
    
    function [img,use_hdr] = loadImage(fpath) %#ok
        use_hdr = true;
        try
            img = load_nii(fpath);
        catch
            warning('Possible non-orthogonal shearing detected in affine matrix, trying load_untouch_nii.m')
            try
                img = load_untouch_nii(fpath);
                use_hdr = false;
            catch
                error(['ERROR: Unable to load image',fpath,'.'])
            end
        end
    end

end