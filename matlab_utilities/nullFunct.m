function nullFunct(fname)
    funct = load_untouch_nii(fname);
    funct.hdr.dime.pixdim(6:8) = 1;
    funct.hdr.hist.qform_code = 0;
    funct.hdr.hist.sform_code = 1;
    funct.hdr.hist.quatern_b = 0;
    funct.hdr.hist.quatern_c = 0;
    funct.hdr.hist.quatern_d = 0;
    funct.hdr.hist.qoffset_x = 0;
    funct.hdr.hist.qoffset_y = 0;
    funct.hdr.hist.qoffset_z = 0;
    save_untouch_nii(funct,fname)
end