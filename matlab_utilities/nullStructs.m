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