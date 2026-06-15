function GPU_layout = RectangularFin(Istart, Iend, Jstart, Jend, NPI, NPJ, l_base_frac, h_base_frac)

global J_fluid_bottom J_fluid_top

GPU_layout = zeros(Iend, Jend);

% Geometry setup
Start_L_base = ceil(l_base_frac * (NPI + 1));
End_limit = ceil((1 - l_base_frac) * (NPI + 1));

H_domain = (NPJ + 1);
Start_H_bottom = J_fluid_bottom;
Start_H_top = J_fluid_top;

H_rectangle = round(0.06 * H_domain);

% Loop over domain
for I = Istart:Iend
    for J = Jstart:Jend
        % Loop over rectangle positions
        if (I >= Start_L_base) && (I <= End_limit)

            Channel_Height = round((Start_H_top - Start_H_bottom - 4*H_rectangle)/3);  
            Rectangle_relative_Height = H_rectangle + Channel_Height;
            
            bottom_wall_rectangle = Start_H_bottom + H_rectangle; 
            top_wall_rectangle = Start_H_top - H_rectangle+1;

            lower_rectangle1 = Start_H_bottom + Rectangle_relative_Height;
            upper_rectangle1 = Start_H_bottom + Rectangle_relative_Height + H_rectangle;

            lower_rectangle2 = Start_H_top - Rectangle_relative_Height - H_rectangle;
            upper_rectangle2 = Start_H_top - Rectangle_relative_Height;

            % Fill regions
            if (J < bottom_wall_rectangle)
                GPU_layout(I, J) = 1;
            end

            if (J > top_wall_rectangle)
                GPU_layout(I, J) = 1;
            end

            if (J > lower_rectangle1 && J < upper_rectangle1)
                GPU_layout(I, J) = 1;
            end

            if (J > lower_rectangle2 && J < upper_rectangle2)
                GPU_layout(I, J) = 1;
            end
        end
    end
end
end
