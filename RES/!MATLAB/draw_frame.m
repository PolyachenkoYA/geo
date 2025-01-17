function Amax = draw_frame(frame, ax, Amax, kA, tit, ax_bounds, colors)
    if(frame.Nrays == 0)        
        cla(ax);        
    else                
        auto_norm = qop(isempty(Amax), 1, 0); % this must be a scalar double because of matlab.
%         if(~auto_norm)
%             Amax = max(max(frame.A), Amax);
%         end
        
        draw_ind = (frame.type == 1);
        if(auto_norm)
            scale_k = 1;
        else
            scale_k = frame.A(draw_ind) ./ frame.c(draw_ind) / Amax * kA;
        end        
        quiver3(ax, frame.X(draw_ind,1), frame.X(draw_ind,2), frame.X(draw_ind,3),...
                    frame.V(draw_ind,1) .* scale_k, frame.V(draw_ind,2) .* scale_k, frame.V(draw_ind,3) .* scale_k,...
                    auto_norm ,'Color', colors.p);
                
                
        hold on;                
        draw_ind = frame.type == 2;
        if(auto_norm)
            scale_k = 1;
        else            
            scale_k = frame.A(draw_ind) ./ frame.c(draw_ind) / Amax * kA;
        end                
        quiver3(ax, frame.X(draw_ind,1), frame.X(draw_ind,2), frame.X(draw_ind,3),...
                    frame.V(draw_ind,1) .* scale_k, frame.V(draw_ind,2) .* scale_k, frame.V(draw_ind,3) .* scale_k,...
                    auto_norm,'Color', colors.s);
        hold off;
        
        axis(ax_bounds);
        legend('P wave','S wave');
        %view(0, 90);
        view(60, 25); % 3D view
    end
    
    title(tit);
end

