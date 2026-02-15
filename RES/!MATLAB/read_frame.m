function frame = read_frame(filename)
    data = dlmread(filename);
    frame.Nrays = data(1,1);
    data(1,:) = [];
    if(frame.Nrays == 0)
        return;
    end

    frame.type = data(:,1);
    frame.t = data(:,2);
    frame.A = data(:,3);
    frame.c = data(:,4);
    frame.X = data(:,5:7);
    frame.V = data(:,8:10);
end

