function draw_movie(model_name, N0, N1)
    model_path = fullfile('.', 'DATA', model_name);
    frames_path = fullfile(model_path, 'frames');
    
    close(gcf);
    fig = figure('Name', 'raw rays movie', 'units', 'normalized', 'outerposition', [0 0 0.73 1]);
    axes(fig);
    hold off;
    aviobj = VideoWriter(fullfile(model_path, 'movie.avi'));
    open(aviobj);
    
    prm = read_params(fullfile('.', 'DATA', model_name, [model_name, '.prm']));
    if(prm.Nfrm == 0)
        disp('error: Nfrm == 0');
        return;
    end
    colors.p = [1 0 0];
    colors.s = [0 0 1];
    
    Amax = 0;
    if(~exist('N0', 'var'))
        N0 = 1;
    end
    if(~exist('N1', 'var'))
        N1 = prm.Nfrm - 1;
    end    
    textprogressbar('creating movie: ');
    for i = N0:N1
        frame = read_frame(fullfile(frames_path, [num2str(i) '.frm']));
        if(i == N0)
            Amax = max(frame.A);
        end
        
%         Amax = draw_frame(frame, get(fig,'CurrentAxes'), Amax, 250, ['t = ' num2str(prm.dt*i)],...
%                    [prm.Xmin(1) prm.Xmax(1) prm.Xmin(2) prm.Xmax(2) prm.Xmin(3) prm.Xmax(3)], colors);
        if(frame.Nrays > 1)
            draw_frame(frame, get(fig,'CurrentAxes'), Amax, qop(prm.prnt_mode == 1, 250, 1), ['t = ' num2str(prm.dt*i)],...
                       [prm.Xmin(1) prm.Xmax(1) prm.Xmin(2) prm.Xmax(2) prm.Xmin(3) prm.Xmax(3)], colors);                   
            writeVideo(aviobj,getframe(fig));                   
        end
        textprogressbar((i+1) / prm.Nfrm * 100);        
    end
    textprogressbar(' done');
        
    close(aviobj);
    close(fig);        
end

