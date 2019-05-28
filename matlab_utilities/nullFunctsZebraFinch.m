function nullFunctsZebraFinch(functs)
    for ix = 1:length(functs)
        for iy = 1:length(functs{ix})
            funct = load_untouch_nii(functs{ix}{iy});
            funct.hdr.dime.pixdim(6:8) = 1;
            funct.hdr.dime.pixdim(2:4) = funct.hdr.dime.pixdim([3,2,4]);
            funct.hdr.dime.pixdim(5) = funct.hdr.dime.pixdim(5)*3; % RARE factor 3
            funct.hdr.hist.srow_x = [funct.hdr.dime.pixdim(2),0,0,funct.hdr.dime.pixdim(2)];
            funct.hdr.hist.srow_y = [0,funct.hdr.dime.pixdim(3),0,funct.hdr.dime.pixdim(3)];
            funct.hdr.hist.srow_z = [0,0,funct.hdr.dime.pixdim(4),funct.hdr.dime.pixdim(4)];
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