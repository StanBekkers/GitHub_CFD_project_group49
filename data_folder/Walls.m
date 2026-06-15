function GPU_layout = Walls(Istart, Iend, Jstart, Jend, NPI, NPJ, h_base_frac)

global J_fluid_bottom J_fluid_top

% Initialize grid
GPU_layout = zeros(Iend, Jend);

for I = Istart:Iend
    for J = Jstart:Jend
        % lower wall
        if J < J_fluid_bottom 
            GPU_layout(I, J) = 1;
        end

        % upper wall
        if J > J_fluid_top
            GPU_layout(I, J) = 1;
        end
    end
end
end
