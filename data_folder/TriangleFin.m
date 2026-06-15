function GPU_layout = TriangleFin(Istart, Iend, Jstart, Jend, NPI, NPJ, l_base_frac, h_base_frac)

GPU_layout = zeros(Iend, Jend);

% Geometry setup
L_triangle = ceil(0.05 * (NPI + 1));
Start_L_base = ceil(l_base_frac * (NPI + 1));
End_limit = ceil((1 - l_base_frac) * (NPI + 1));

H_domain = (NPJ + 1);
Start_H_bottom = ceil(h_base_frac * H_domain);
Start_H_top = H_domain - Start_H_bottom;

H_triangle = ceil(0.08 * H_domain);
slope = H_triangle / L_triangle;

% Loop over domain
for I = Istart:Iend
    for J = Jstart:Jend
k = 0;
        % Loop over triangle positions
        for offset = 0:L_triangle:(End_limit - Start_L_base - L_triangle)

            Start_L_triangle = Start_L_base + offset;
            End_L_triangle   = Start_L_triangle + L_triangle;
             
             if mod(k, 2) == 0
                  s = 1;
                  H_triangle_bottom = ceil(0.08 * H_domain);
                  H_triangle_top = 0;
             else
                  s = -1;
                  H_triangle_bottom = 0;
                  H_triangle_top = ceil(0.08 * H_domain);
             end

            k = k + 1; % update triangle number
            if (I >= Start_L_triangle) && (I <= End_L_triangle)

                i_shift = I - Start_L_triangle;

                lower_line = ceil(s*-i_shift*slope + H_triangle_bottom + Start_H_bottom);
                upper_line = ceil(s*-i_shift*slope + Start_H_top - H_triangle_top);

                Channel_Height = floor((Start_H_top - Start_H_bottom - 3*H_triangle)/3);  
                ZigZag_relative_Height = H_triangle + Channel_Height;
                
                lower_zigzag1 = lower_line + ZigZag_relative_Height - H_triangle;
                upper_zigzag1 = lower_line + ZigZag_relative_Height ;

                lower_zigzag2 = upper_line - ZigZag_relative_Height;
                upper_zigzag2 = upper_line - ZigZag_relative_Height + H_triangle;

                % Fill regions
                if (J < lower_line)
                    GPU_layout(I, J) = 1;
                end

                if (J > upper_line)
                    GPU_layout(I, J) = 1;
                end

                if (J > lower_zigzag1 && J < upper_zigzag1)
                    GPU_layout(I, J) = 1;
                end

                if (J > lower_zigzag2 && J < upper_zigzag2)
                    GPU_layout(I, J) = 1;
                end

            end

        end

    end
end

end