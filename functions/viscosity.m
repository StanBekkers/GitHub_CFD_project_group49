function [] = viscosity()
% Purpose: Calculate turbulent and effective viscosity using k-epsilon model.
% mut  = rho * Cmu * k^2 / (eps + SMALL)
% mueff = mu + mut   (used in momentum, k, eps diffusion)
% Gamma is also updated with turbulent thermal diffusivity:
%   Fluid:  Gamma = 0.6/Cp + mut/Pr_t
%   Solid:  Gamma = 401.0/Cp
 
% constants
global NPI NPJ Cmu SMALL
% variables
global rho k eps mu mut mueff Gamma Cp
global h_base_frac l_base_frac
 
h_base_frac = 2/10;
l_base_frac = 3/10;

Istart = 1;
Iend = NPI+2;  
Jstart = 1;
Jend = NPJ+2; 

Pr_t = 0.9;   % turbulent Prandtl number for water at ~20 deg C

layout_wall = Walls(Istart, Iend, Jstart, Jend, NPI, NPJ, h_base_frac);
layout_fins = TriangleFin(Istart, Iend, Jstart, Jend, NPI, NPJ, l_base_frac, h_base_frac);
cooler_layout = layout_wall | layout_fins;

for I = Istart:Iend 
    for J = Jstart:Jend 
        % 1. Calculate turbulent and effective viscosity
        mut(I,J)   = rho(I,J) * Cmu * k(I,J)^2 / (eps(I,J) + SMALL);
        
        % ---- Add Viscosity Ratio Limiter here ----
        mut_max = 1000.0 * mu(I,J); % Enforce a maximum turbulent viscosity ratio of 1000
        if mut(I,J) > mut_max
            mut(I,J) = mut_max;
        end
        
        mueff(I,J) = mu(I,J) + mut(I,J);
        
        % 2. Identify if cell is in the solid baseplate/fins
        is_solid = false;
        if (cooler_layout(I,J) == 1)
            is_solid = true;
        end
        
        % 3. Apply Gamma according to physical zone
        if is_solid
            Gamma(I,J) = 401.0 / Cp(I,J); % Solid copper
        else
            Gamma(I,J) = 0.6 / Cp(I,J) + mut(I,J) / Pr_t; % Fluid water + turbulent diffusivity
        end
    end
end
end