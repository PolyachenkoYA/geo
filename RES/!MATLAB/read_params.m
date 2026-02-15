function prm = read_params(filepath)
    data = importdata(filepath);
    %           Nrays           Tmax             dt           Amin
    %               f            tau            eps        use_det
    %        draw_mov      prnt_mode         
    %          Xmin.x         Xmin.y         Xmin.z
    %          Xmax.x         Xmax.y         Xmax.z
    %            dX.x           dX.y           dX.z
    
    k = 0;
    k = k+1; prm.Nrays = data.data(k);
    k = k+1; prm.Tmax = data.data(k);
    k = k+1; prm.dt = data.data(k);
    k = k+1; prm.Amin = data.data(k);
    k = k+1; prm.f = data.data(k); prm.f = prm.f * 10^6;
    k = k+1; prm.tau = data.data(k); prm.tau = prm.tau / 10^6;
    k = k+1; prm.eps = data.data(k);
    k = k+1; prm.use_det = data.data(k);
    k = k+1; prm.draw_mov = data.data(k);
    k = k+1; prm.prnt_mode = data.data(k);
    k = k+1; prm.Xmin(1) = data.data(k); 
    k = k+1; prm.Xmin(2) = data.data(k);   
    k = k+1; prm.Xmin(3) = data.data(k);   
    k = k+1; prm.Xmax(1) = data.data(k);   
    k = k+1; prm.Xmax(2) = data.data(k);   
    k = k+1; prm.Xmax(3) = data.data(k);       
    k = k+1; prm.dX(1) = data.data(k);   
    k = k+1; prm.dX(2) = data.data(k);   
    k = k+1; prm.dX(3) = data.data(k);           
    
    prm.Nfrm = floor(prm.Tmax * (1 + prm.eps) / prm.dt) + 1;
end

