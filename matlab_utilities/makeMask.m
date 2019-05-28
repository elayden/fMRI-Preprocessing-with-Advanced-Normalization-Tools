function makeMask(GM,prob,outName,new_res)
% makeMask creates a thresholded mask from grey matter probability image(s)
% 
% Info: accepts 'GM' as either a file path to an image containing grey
% matter probabilities, or a cell array containing file paths to multiple
% images. In the first case, 'mask' is simply this image thresholded at the
% specified probability 'prob' (default: .2). In the second case, mask is
% thresholded at 'prob' for the average across multiple GM probability
% images.
% 
% Optional Inputs:
%           new_res,    [1 x 3 vector] specifying voxel dimensions of a
%                       desired resampled mask (e.g., resample to
%                       functional dimensions); if option is selected, both
%                       an original space mask and a resampled mask will be
%                       output. To the latter, 'r_' will be prepended 
% 
% Author: Elliot Layden
% U Chicago, August 2017

if nargin < 2 || isempty(prob)
    prob = .2;
end

use_hdr = true;
if ischar(GM)
    try
        img = load_nii(GM);
    catch
        warning('Possible non-orthogonal shearing detected in affine matrix, trying load_untouch_nii.m')
        try
            img = load_untouch_nii(GM);
            use_hdr = false;
        catch
            error('ERROR: Unable to load image ''GM''.')
        end
    end
    mImage = img.img;
    [fpath,~,ext] = fileparts(GM);
elseif iscell(GM)
    N = length(GM); 
    try
        img = load_nii(GM{1});
    catch
        warning('Possible non-orthogonal shearing detected in affine matrix, trying load_untouch_nii.m')
        try
            img = load_untouch_nii(GM{1});
            use_hdr = false;
        catch
            error('ERROR: Unable to load image ''GM''.')
        end
    end
    mImage = img.img;
    for i = 2:N
        if use_hdr
            img = load_nii(GM{i});
        else
            img = load_untouch_nii(GM{i});
        end
        mImage = mImage + img.img;
    end 
    mImage = mImage/N;
    [fpath,~,ext] = fileparts(GM{1});
else
    error('ERROR: Input ''GM'' should either be type char/string or cell array.')
end

maskImg = mImage>=prob;
img.img = maskImg;

if nargin < 3 || isempty(outName)
    outName = fullfile(fpath,['mask',ext]);
end

if use_hdr
    save_nii(img,outName);
else
    save_untouch_nii(img,outName);
end

if nargin>3 && ~isempty(new_res)
    resize_img(outName, new_res, [nan nan nan; nan nan nan])
end

end