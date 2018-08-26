function Amax = draw_frame(frame, ax, Amax, kA, tit, ax_bounds, colors)
    if(frame.Nrays == 0)        
        cla(ax);        
    else                
        do_norm_v = qop(isempty(Amax), 1, 0);
        if(~do_norm_v)
            Amax = max(max(frame.A), Amax);
        end
        
        draw_ind = ((frame.type == 1) .* (frame.A > 0)) > 0;
        if(do_norm_v)
            scale_k = 1;
        else
            scale_k = frame.A(draw_ind) ./ frame.c(draw_ind) / Amax * kA;
        end        
        quiver3(ax, frame.X(draw_ind,1), frame.X(draw_ind,2), frame.X(draw_ind,3),...
                    frame.V(draw_ind,1) .* scale_k, frame.V(draw_ind,2) .* scale_k, frame.V(draw_ind,3) .* scale_k,...
                    do_norm_v ,'Color', colors.p);
                
                
        hold on;                
        draw_ind = ((frame.type == 2) .* (frame.A > 0)) > 0;
        if(do_norm_v)
            scale_k = 1;
        else            
            scale_k = frame.A(draw_ind) ./ frame.c(draw_ind) / Amax * kA;
        end                
        quiver3(ax, frame.X(draw_ind,1), frame.X(draw_ind,2), frame.X(draw_ind,3),...
                    frame.V(draw_ind,1) .* scale_k, frame.V(draw_ind,2) .* scale_k, frame.V(draw_ind,3) .* scale_k,...
                    do_norm_v,'Color', colors.sv);
        hold off;
        
        axis(ax_bounds);
        legend('P wave','SV wave');
        view(0, 90);        
    end
    
    title(tit);
end

