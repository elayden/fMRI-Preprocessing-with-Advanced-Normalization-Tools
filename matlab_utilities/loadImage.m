function [img,use_hdr] = loadImage(fpath)
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