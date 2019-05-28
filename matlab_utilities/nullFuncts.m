function nullFuncts(functs)
    for ix = 1:length(functs)
        for iy = 1:length(functs{ix})
            funct = load_untouch_nii(functs{ix}{iy});
            funct.hdr.dime.pixdim(6:8) = 1;
            funct.hdr.hist.qform_code = 0;
            funct.hdr.hist.sform_code = 1;
            funct.hdr.hist.quatern_b = 0;
            funct.hdr.hist.quatern_c = 0;
            funct.hdr.hist.quatern_d = 0;
            funct.hdr.hist.qoffset_x = 0;
            funct.hdr.hist.qoffset_y = 0;
            funct.hdr.hist.qoffset_z = 0;
            save_untouch_nii(funct,functs{ix}{iy})
        end
    end
end