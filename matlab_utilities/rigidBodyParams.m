function motion = rigidBodyParams(moco,outName)
% rigidBodyParams accepts either a GenericAffine.mat file output by ANTs 
%   (these store the affine matrices from realignment in vector form) or
%   a *MOCOparams.csv file as output by antsMotionCorr.exe
%   -this function converts inputs to the 6 rigid-body motion parameters 
%   commonly used for nuisance regression.
% 
% INPUT:
%       'moco',     file path to a MOCOparams.csv or a GenericAffine.mat
%                   file as used by ITK / ANTs
%       'outName',  the requested output name and filetype (full file path,
%                   ideally); available extensions are .mat or .txt; a .mat
%                   will include a single variable 'motion' which is a
%                   nTime x 6 matrix of rigid-body motion parameters; a
%                   .txt will simply write these out to a comma-delimited
%                   .txt file
% 
% INFO:
% For reference on how ITK affines are stored, see ANTs manual (page 6), or: 
% https://afni.nimh.nih.gov/afni/community/board/read.php?1,94478,94478#msg-94478
% AffineTransform_double_3_3 = size(1,12) matrix, where affine:
% [1, 2, 3, 10;...
%  4, 5, 6, 11;...
%  7, 8, 9, 12]
% fixed = origin
% 
% For reference on how to extract 6 rigid-body motion parameters from
% affine, see:
% <http://www.fil.ion.ucl.ac.uk/spm/doc/books/hbf2/pdfs/Ch2.pdf> , p. 9
% &
% Smolic, A., Makai, B., Lin, G., & Sikora, T. (1997, September). Estimation 
% of motion parameters of a rigid body from a monocular image sequence for 
% MPEG-4 applications. In International Workshop on Interactive Distributed 
% Multimedia Systems and Telecommunication Services (pp. 11-19). 
% Springer Berlin Heidelberg. p. 12
% <http://download.springer.com/static/pdf/175/chp%253A10.1007%252FBFb0000335.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Fchapter%2F10.1007%2FBFb0000335&token2=exp=1478282521~acl=%2Fstatic%2Fpdf%2F175%2Fchp%25253A10.1007%25252FBFb0000335.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Fchapter%252F10.1007%252FBFb0000335*~hmac=408a3699018b1ced2309c11624686cca701e91607da44382239e678dfcc0f940>
% 
% Author: Elliot Layden, University of Chicago, August 2017
%%
if nargin < 1 || isempty(moco) || ~exist(moco,'file')
    error('ERROR: input ''moco'' not found.')
end

[fpath,fname,ext] = fileparts(moco);

% Output:
if nargin<2 || isempty(outName)
    outName = fullfile(fpath,['rp_',fname,'.mat']);
    outMat = true;
else
    [~,~,outext] = fileparts(outName);
    switch outext
        case '.mat'
            outMat = true;
        case '.txt',
            outMat = false;
        otherwise
            error('ERROR: input ''outName'' should have extension .mat or .txt.')
    end
end
        
switch ext
    case '.mat'
        motion = zeros(1,6);
        load(moco)
        motion(1:3) = AffineTransform_double_3_3(10:12); % translation
        motion(5) = asin(AffineTransform_double_3_3(3)); % y-rotation parameter (phi_y)
        cos_phi_y = cos(motion(5));
        motion(4) = atan2(AffineTransform_double_3_3(6)/cos_phi_y,AffineTransform_double_3_3(9)/cos_phi_y); % x-rotation parameter (phi_x)
        motion(6) = atan2(AffineTransform_double_3_3(2)/cos_phi_y,AffineTransform_double_3_3(1)/cos_phi_y); % z-rotation parameter (phi_z)
        nTime = 1;
    case '.csv'
        moco_data = importdata(moco);
        affine = moco_data.data(:,3:end);
        nTime = size(affine,1);
        motion = zeros(nTime,6);
%         for i = 1:nTime
%             motion(i,1:3) = affine(i,10:12); % translation
%             motion(i,5) = asin(affine(i,3)); % y-rotation parameter (phi_y)
%             cos_phi_y = cos(motion(i,5));
%             motion(i,4) = atan2(affine(i,6)/cos_phi_y,affine(i,9)/cos_phi_y); % x-rotation parameter (phi_x)
%             motion(i,6) = atan2(affine(i,2)/cos_phi_y,affine(i,1)/cos_phi_y); % z-rotation parameter (phi_z)
%         end
        motion = affine;
    otherwise
        error('ERROR: input ''moco'' should have extension .mat or .txt.')
end

if strcmp(ext,'.mat')
    motion(:,3) = -motion(:,3); % makes consistent with SPM
    motion(:,6) = -motion(:,6); % makes consistent with SPM
end

if outMat
    save(outName,'motion')
else
    write_data_cell = cell(nTime,1);
    for i = 1:nTime
        write_data_cell{i} = [num2str(motion(i,1)),',',num2str(motion(i,2)),',',num2str(motion(i,3)),',',num2str(motion(i,4)),',',num2str(motion(i,5)),',',num2str(motion(i,6))];
    end
    write_data = char(write_data_cell);
    fileID = fopen(outName,'w');
    for i = 1:size(write_data,1)
        fprintf(fileID,'%s\r\n',write_data(i,:));
    end
    fclose(fileID);
end

end
