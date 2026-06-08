function [] = kcoeff()
% Purpose: Calculate coefficients for the k (turbulent kinetic energy) equation.
%
% Geometry is read directly from the same parameters used in ucoeff.m and vcoeff.m.
% The solid/fluid mask is computed identically — no separate geometry input needed.
%
% Wall treatment (channel walls and fin surfaces):
%   k_wall = u_tau^2 / sqrt(Cmu)
%   Applied via large-number Dirichlet: SP = -LARGE, Su = LARGE * k_wall
%   This forces k_P → k_wall without floating-point blow-up.
%
% Interior cells: standard k-eps production/destruction source terms.

% constants
global NPI NPJ Cmu sigmak LARGE SMALL
% variables
global x x_u y y_v SP Su F_u F_v mut rho u mu tw uplus Istart Iend ...
    Jstart Jend b aE aW aN aS aP k k_old eps E2

Istart = 2;
Iend   = NPI+1;
Jstart = 2;
Jend   = NPJ+1;

convect();
viscosity();
calculateuplus();

% ---- Geometry parameters (identical to ucoeff.m / vcoeff.m) ----
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

        % Diffusion with harmonic mean of mut (guards against zero with SMALL)
        Dw = (mut(I-1,J)+SMALL)*(mut(I,J)+SMALL)/sigmak / ...
             ((mut(I-1,J)+SMALL)*(x(I)-x_u(i)) + (mut(I,J)+SMALL)*(x_u(i)-x(I-1))) * AREAw;
        De = (mut(I,J)+SMALL)*(mut(I+1,J)+SMALL)/sigmak / ...
             ((mut(I,J)+SMALL)*(x(I+1)-x_u(i+1)) + (mut(I+1,J)+SMALL)*(x_u(i+1)-x(I))) * AREAe;
        Ds = (mut(I,J-1)+SMALL)*(mut(I,J)+SMALL)/sigmak / ...
             ((mut(I,J-1)+SMALL)*(y(J)-y_v(j)) + (mut(I,J)+SMALL)*(y_v(j)-y(J-1))) * AREAs;
        Dn = (mut(I,J)+SMALL)*(mut(I,J+1)+SMALL)/sigmak / ...
             ((mut(I,J)+SMALL)*(y(J+1)-y_v(j+1)) + (mut(I,J+1)+SMALL)*(y_v(j+1)-y(J))) * AREAn;

        % Default: interior production - destruction
        P_k     = 2.0 * mut(I,J) * E2(I,J);
        SP(I,J) = -rho(I,J) * eps(I,J) / (k(I,J) + SMALL);
        Su(I,J) = P_k;
        isWall  = false;

        % ---- Channel bottom wall ----
        if J == J_fluid_bottom
            u_tau   = sqrt(abs(tw(I, J_fluid_bottom)) / (rho(I, J_fluid_bottom) + SMALL));
            k_wall  = u_tau^2 / (sqrt(Cmu) + SMALL);
            SP(I,J) = -LARGE;
            Su(I,J) =  LARGE * k_wall;
            isWall  = true;
        end

        % ---- Channel top wall ----
        if J == J_fluid_top
            u_tau   = sqrt(abs(tw(I, J_fluid_top)) / (rho(I, J_fluid_top) + SMALL));
            k_wall  = u_tau^2 / (sqrt(Cmu) + SMALL);
            SP(I,J) = -LARGE;
            Su(I,J) =  LARGE * k_wall;
            isWall  = true;
        end

        % ---- Fin / solid cells: same geometry as ucoeff.m ----
        for offset = 0:L_triangle:(End_limit - Start_L_base - L_triangle)
            Start_L_triangle = Start_L_base + offset;
            End_L_triangle   = Start_L_triangle + L_triangle;

            if (i >= Start_L_triangle) && (i <= End_L_triangle)
                i_shift = i - Start_L_triangle;
                lower_line = ceil(-i_shift*slope + H_triangle + Start_H_bottom);
                upper_line = ceil(-i_shift*slope + Start_H_top);

                band_half     = floor((upper_line - lower_line) / 6);
                lower_zigzag1 = lower_line + band_half;
                upper_zigzag1 = lower_line + 2*band_half;
                lower_zigzag2 = upper_line - 2*band_half;
                upper_zigzag2 = upper_line - band_half;

                isSolid = (J < lower_line) || (J > upper_line) || ...
                          (J > lower_zigzag1 && J < upper_zigzag1) || ...
                          (J > lower_zigzag2 && J < upper_zigzag2);

                if isSolid && ~isWall
                    % Fin surface: estimate u_tau from local laminar shear
                    % (fin surface cells have near-zero velocity)
                    u_loc   = abs(0.5*(u(i,J) + u(i+1,J)));
                    y_P     = 0.5 * AREAw;
                    u_tau   = sqrt(mu(I,J) * u_loc / (y_P * rho(I,J) + SMALL));
                    k_wall  = u_tau^2 / (sqrt(Cmu) + SMALL);
                    SP(I,J) = -LARGE;
                    Su(I,J) =  LARGE * k_wall;
                    isWall  = true;
                end
            end
        end

        Su(I,J) = Su(I,J) * AREAw * AREAs;
        SP(I,J) = SP(I,J) * AREAw * AREAs;

        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);

        if isWall || J == J_fluid_bottom
            aS(I,J) = 0.;
        else
            aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        end

        if isWall || J == J_fluid_top
            aN(I,J) = 0.;
        else
            aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
        end

        % Steady-state: no time derivative term (matches GPU_Cooler_V1 style)
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
        b(I,J)  = Su(I,J);
    end
end
end