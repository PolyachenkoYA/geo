function draw_movie(model_name)
    model_path = fullfile('.', 'DATA', model_name);
    frames_path = fullfile(model_path, 'frames');
    
    close(gcf);
    fig = figure('Name', 'raw rays movie', 'units', 'normalized', 'outerposition', [0 0 0.73 1]);
    axes(fig);
    hold off;
    aviobj = VideoWriter(fullfile(model_path, 'movie.avi'));
    open(aviobj);
    
    prm = read_params(fullfile('.', 'DATA', model_name, [model_name, '.prm']));
    colors.p = [1 0 0];
    colors.sv = [0 0 1];
    
    Amax = 0;
    textprogressbar('creating movie: ');
    for i = 0:(prm.Nfrm - 1)
        frame = read_frame(fullfile(frames_path, [num2str(i) '.frm']), prm.prnt_mode);
        
%        Amax = draw_frame(frame, get(fig,'CurrentAxes'), Amax, ['t = ' num2str(prm.dt*i)],...
%                   [prm.Xmin(1) prm.Xmax(1) prm.Xmin(2) prm.Xmax(2) prm.Xmin(3) prm.Xmax(3)], colors);          
        Amax = draw_frame(frame, get(fig,'CurrentAxes'), Amax, 250, ['t = ' num2str(prm.dt*i)],...
                   [prm.Xmin(1) prm.Xmax(1) prm.Xmin(2) prm.Xmax(2) prm.Xmin(3) prm.Xmax(3)], colors);          
               
        textprogressbar((i+1) / prm.Nfrm * 100);
    
        writeVideo(aviobj,getframe(fig));
    end
    textprogressbar(' done');
        
    close(aviobj);
    close(fig);        
end

