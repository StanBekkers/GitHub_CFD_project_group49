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
 
Pr_t = 0.9;   % turbulent Prandtl number for water at ~20 deg C
 
% Geometry limits to check solid zones
h_base_frac = 2/10;
l_base_frac = 3/10;
L_triangle = ceil(0.05*(NPI+1));
Start_L_base = ceil(l_base_frac*(NPI+1));
End_limit = ceil((1 - l_base_frac)*(NPI+1));   
H_domain = (NPJ+1);
Start_H_bottom = ceil(h_base_frac*H_domain);
Start_H_top = H_domain - Start_H_bottom;
H_triangle = ceil((1/4)*h_base_frac * H_domain);   
slope = H_triangle / L_triangle;

for I = 1:NPI+2
    for J = 1:NPJ+2
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
        if (J < ceil(h_base_frac*(NPJ+1))) || (J > ceil((1-h_base_frac)*(NPJ+1)))
            is_solid = true;
        end
        
        for offset = 0:L_triangle:(End_limit - Start_L_base - L_triangle)
            Start_L_triangle = Start_L_base + offset;
            End_L_triangle = Start_L_triangle + L_triangle;
            
            if (I >= Start_L_triangle && I <= End_L_triangle)
                i_shift = I - Start_L_triangle;
                lower_line = ceil(-i_shift*slope + H_triangle + Start_H_bottom);
                upper_line = ceil(-i_shift*slope + Start_H_top);
                
                band_half = floor((upper_line - lower_line) / 6);
                lower_zigzag1 = lower_line + band_half;
                upper_zigzag1 = lower_line + 2*band_half;
                lower_zigzag2 = upper_line - 2*band_half;
                upper_zigzag2 = upper_line - band_half;

                if (J < lower_line) || (J > upper_line) || ...
                   (J > lower_zigzag1 && J < upper_zigzag1) || ...
                   (J > lower_zigzag2 && J < upper_zigzag2)
                    is_solid = true;
                end
            end
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