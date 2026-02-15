function data = loadDataFromFile(Fname)
    fid=fopen(Fname,'r');
    data = [];
    if(fid==-1)
        waitfor(errordlg({'error in opening file',Fname,'for reading (r)'},'File error')); 
        return;
    else
        data = fscanf(fid,'%f ');
        if(fclose(fid))
            waitfor(errordlg({'error in closing',Fname,'which was opened for reading (r)'},'File error')); 
            return;            
        end
    end
end
