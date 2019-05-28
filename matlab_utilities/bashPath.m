function bPath = bashPath(fpath)

    if isempty(strfind(fpath,':'))
        error('Must provide full path.')
    end
    
    bPath = ['/cygdrive/',lower(fpath(1))];
    fpath(1:3) = [];
    
    nFinished = true;
    while (nFinished && numel(fpath)>0)
        ind = strfind(fpath,'\');
        if isempty(ind)
            nFinished = false;
            bPath = [bPath, '/', fpath(1:end)]; %#ok
        else
            bPath = [bPath, '/', fpath(1:ind-1)]; %#ok
            fpath(1:ind) = [];
        end
    end
    
end