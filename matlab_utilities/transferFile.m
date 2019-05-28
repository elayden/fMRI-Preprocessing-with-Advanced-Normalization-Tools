    function transferFile(indir,outdir)
        try
            movefile(indir,outdir,'f')
        catch
            warning(['Could not find ',indir,'.'])
        end
    end