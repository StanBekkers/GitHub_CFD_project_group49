function [] = epscoeff()
% Purpose: Calculate coefficients for the eps (dissipation rate) equation.

% constants
global NPI NPJ Cmu LARGE SMALL sigmaeps kappa C1eps C2eps
% variables
global x x_u y y_v SP Su F_u F_v mut rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP k eps eps_old E2 mu

Istart = 2;
Iend   = NPI+1;
Jstart = 2;
Jend   = NPJ+1;

convect();
viscosity();

% ---- Geometry parameters ----
h_base_frac    = 2/10;
l_base_frac    = 3/10;
L_triangle     = ceil(0.05*(NPI+1));
Start_L_base   = ceil(l_base_frac*(NPI+1));
End_limit      = ceil((1 - l_base_frac)*(NPI+1));
H_domain       = (NPJ+1);
Start_H_bottom = ceil(h_base_frac*H_domain);
Start_H_top    = H_domain - Start_H_bottom;
H_triangle     = ceil((1/4)*h_base_frac * H_domain);
slope          = H_triangle / L_triangle;

J_fluid_bottom = ceil(h_base_frac*(NPJ+1));
J_fluid_top    = ceil((1-h_base_frac)*(NPJ+1));

% Effective turbulent diffusivity for eps (includes molecular viscosity)
Gamma_eps = mu + mut / sigmaeps;

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;

        AREAw = y_v(j+1) - y_v(j);
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i);
        AREAn = AREAs;

        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;

        % ---- Solid check to deactivate eps-transport in solid domains ----
        is_solid = false;
        if (J < J_fluid_bottom) || (J > J_fluid_top)
            is_solid = true;
        end
        
        for offset = 0:L_triangle:(End_limit - Start_L_base - L_triangle)
            Start_L_triangle = Start_L_base + offset;
            End_L_triangle   = Start_L_triangle + L_triangle;

            if (I >= Start_L_triangle) && (I <= End_L_triangle)
                i_shift = I - Start_L_triangle;
                lower_line = ceil(-i_shift*slope + H_triangle + Start_H_bottom);
                upper_line = ceil(-i_shift*slope + Start_H_top);

                band_half     = floor((upper_line - lower_line) / 6);
                lower_zigzag1 = lower_line + band_half;
                upper_zigzag1 = lower_line + 2*band_half;
                lower_zigzag2 = upper_line - 2*band_half;
                upper_zigzag2 = upper_line - band_half;

                isSolidGeometry = (J < lower_line) || (J > upper_line) || ...
                                  (J > lower_zigzag1 && J < upper_zigzag1) || ...
                                  (J > lower_zigzag2 && J < upper_zigzag2);

                if isSolidGeometry
                    is_solid = true;
                end
            end
        end

        if is_solid
            % Enforce eps floor in solids (Bypass relaxation)
            SP(I,J) = -LARGE;
            Su(I,J) = LARGE * 1e-10;
            
            aW(I,J) = 0.0;
            aE(I,J) = 0.0;
            aS(I,J) = 0.0;
            aN(I,J) = 0.0;
            aP(I,J) = LARGE;
            b(I,J)  = Su(I,J);
            continue;
        end

        % Transport by diffusion (harmonic mean of Gamma_eps)
        Dw = (Gamma_eps(I-1,J)*Gamma_eps(I,J)) / ...
             ((Gamma_eps(I-1,J))*(x(I)-x_u(i)) + (Gamma_eps(I,J))*(x_u(i)-x(I-1))) * AREAw;
        De = (Gamma_eps(I,J)*Gamma_eps(I+1,J)) / ...
             ((Gamma_eps(I,J))*(x(I+1)-x_u(i+1)) + (Gamma_eps(I+1,J))*(x_u(i+1)-x(I))) * AREAe;
        Ds = (Gamma_eps(I,J-1)*Gamma_eps(I,J)) / ...
             ((Gamma_eps(I,J-1))*(y(J)-y_v(j)) + (Gamma_eps(I,J))*(y_v(j)-y(J-1))) * AREAs;
        Dn = (Gamma_eps(I,J)*Gamma_eps(I,J+1)) / ...
             ((Gamma_eps(I,J))*(y(J+1)-y_v(j+1)) + (Gamma_eps(I,J+1))*(y_v(j+1)-y(J))) * AREAn;

        % Default: interior source terms
        SP(I,J) = -C2eps * rho(I,J) * eps(I,J) / (k(I,J) + SMALL);
        Su(I,J) =  C1eps * (eps(I,J) / (k(I,J) + SMALL)) * mut(I,J) * E2(I,J); 
        isWall  = false;

        % ---- Channel bottom wall functions ----
        if J == J_fluid_bottom
            y_P      = 0.5 * (y(J_fluid_bottom) - y(J_fluid_bottom - 1)); 
            eps_wall = Cmu^0.75 * k(I,J)^1.5 / (kappa * y_P + SMALL);
            SP(I,J)  = -LARGE;
            Su(I,J)  =  LARGE * eps_wall;
            isWall   = true;
        end

        % ---- Channel top wall functions ----
        if J == J_fluid_top
            y_P      = 0.5 * (y(J_fluid_top + 1) - y(J_fluid_top)); 
            eps_wall = Cmu^0.75 * k(I,J)^1.5 / (kappa * y_P + SMALL);
            SP(I,J)  = -LARGE;
            Su(I,J)  =  LARGE * eps_wall;
            isWall   = true;
        end

        Su(I,J) = Su(I,J) * AREAw * AREAs;
        SP(I,J) = SP(I,J) * AREAw * AREAs;

        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);

        % Standard coefficient formulation
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
        b(I,J)  = Su(I,J);
        
        % ---- Apply under-relaxation for fluid cells ----
        global relax_eps % ensure access to global relax_eps
        aP(I,J) = aP(I,J) / relax_eps;
        b(I,J)  = b(I,J) + (1 - relax_eps) * aP(I,J) * eps(I,J);
    end
end
end