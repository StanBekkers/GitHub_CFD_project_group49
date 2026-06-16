function GPU_layout = ElipseFin(Istart, Iend, Jstart, Jend, NPI, NPJ, l_base_frac, h_base_frac)

GPU_layout = zeros(Iend, Jend);

% Geometry setup
L_Elipse = ceil(0.05 * (NPI + 1));
Start_L_base = ceil(l_base_frac * (NPI + 1));
End_limit = ceil((1 - l_base_frac) * (NPI + 1));

H_domain = (NPJ + 1);
Start_H_bottom = ceil(h_base_frac * H_domain);
Start_H_top = H_domain - Start_H_bottom;

H_Elipse = ceil(0.08 * H_domain);

% Loop over domain
for I = Istart:Iend
    for J = Jstart:Jend

        % Loop over elipse positions
        for offset = 0:L_Elipse:(End_limit - Start_L_base - L_Elipse)

            Start_L_Elipse = Start_L_base + offset;
            End_L_Elipse   = Start_L_Elipse + L_Elipse;
            
            if (I >= Start_L_Elipse) && (I <= End_L_Elipse)
                i_shift = I - Start_L_Elipse - 0.5*L_Elipse;
                
                Channel_Height = floor((Start_H_top - Start_H_bottom - 4*H_Elipse)/3);  
                Elipse_relative_Height = H_Elipse + Channel_Height;
                
                % Elipse:   
                lower_elipse_line = floor(sqrt(1-(i_shift^2)/(0.5*L_Elipse)^2)*0.5*H_Elipse  + Start_H_bottom);
                upper_elipse_line = ceil(-sqrt(1-(i_shift^2)/(0.5*L_Elipse)^2)*0.5*H_Elipse + Start_H_top);

                lower_Elipse1 = upper_elipse_line - 2*Elipse_relative_Height - H_Elipse;
                upper_Elipse1 = lower_elipse_line + Elipse_relative_Height + 0.5*H_Elipse;

                lower_Elipse2 = upper_elipse_line - Elipse_relative_Height - 0.5*H_Elipse;
                upper_Elipse2 = lower_elipse_line + 2*Elipse_relative_Height + H_Elipse;

                % Fill regions
                if (J < lower_elipse_line)
                    GPU_layout(I, J) = 1;
                end

                if (J > upper_elipse_line)
                    GPU_layout(I, J) = 1;
                end

                if (J > lower_Elipse1 && J < upper_Elipse1)
                    GPU_layout(I, J) = 1;
                end

                if (J > lower_Elipse2 && J < upper_Elipse2)
                    GPU_layout(I, J) = 1;
                end

            end

        end

    end
end

end