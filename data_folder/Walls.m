function GPU_layout = Walls(Istart, Iend, Jstart, Jend, NPI, NPJ, h_base_frac)

% Initialize grid
GPU_layout = zeros(Iend, Jend);

for I = Istart:Iend
    for J = Jstart:Jend

        low_wall_face = ceil(h_base_frac * NPJ +1);
        high_wall_face = ceil((1 - h_base_frac) * NPJ + 1);
        
        % lower wall
        if J < low_wall_face 
            GPU_layout(I, J) = 1;
        end

        % upper wall
        if J > high_wall_face
            GPU_layout(I, J) = 1;
        end

    end
end

end



