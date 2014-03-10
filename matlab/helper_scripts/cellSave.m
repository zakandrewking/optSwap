function cellSave(filename,mycell)

    [nrows,ncols]= size(mycell);

    fid = fopen(filename, 'w');

    for row=1:nrows
        string = '%s';
        for j=1:ncols
            string = [string ' %d'];
        end
        string = [string '\n'];
        
        fprintf(fid, string, mycell{row,:});
    end

    fclose(fid);

end
