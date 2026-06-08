function [] = epscoeff()
% Purpose: Calculate coefficients for the eps (dissipation rate) equation.
%
% Geometry is identical to ucoeff.m / vcoeff.m / kcoeff.m.
%
% Wall treatment:
%   eps_wall = Cmu^(3/4) * k_P^(3/2) / (kappa * y_P)
%   Applied via large-number Dirichlet: SP = -LARGE, Su = LARGE * eps_wall
%
% Interior: standard C1eps/C2eps production/destruction.

% constants
global NPI NPJ Cmu LARGE SMALL sigmaeps kappa C1eps C2eps
% variables
global x x_u y y_v SP Su F_u F_v mut rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP k eps eps_old E2

Istart = 2;
Iend   = NPI+1;
Jstart = 2;
Jend   = NPJ+1;

convect();
viscosity();

% ---- Geometry parameters (identical to ucoeff.m / vcoeff.m / kcoeff.m) ----
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

        Dw = (mut(I-1,J)+SMALL)*(mut(I,J)+SMALL)/sigmaeps / ...
             ((mut(I-1,J)+SMALL)*(x(I)-x_u(i)) + (mut(I,J)+SMALL)*(x_u(i)-x(I-1))) * AREAw;
        De = (mut(I,J)+SMALL)*(mut(I+1,J)+SMALL)/sigmaeps / ...
             ((mut(I,J)+SMALL)*(x(I+1)-x_u(i+1)) + (mut(I+1,J)+SMALL)*(x_u(i+1)-x(I))) * AREAe;
        Ds = (mut(I,J-1)+SMALL)*(mut(I,J)+SMALL)/sigmaeps / ...
             ((mut(I,J-1)+SMALL)*(y(J)-y_v(j)) + (mut(I,J)+SMALL)*(y_v(j)-y(J-1))) * AREAs;
        Dn = (mut(I,J)+SMALL)*(mut(I,J+1)+SMALL)/sigmaeps / ...
             ((mut(I,J)+SMALL)*(y(J+1)-y_v(j+1)) + (mut(I,J+1)+SMALL)*(y_v(j+1)-y(J))) * AREAn;

        % Default: interior source terms
        SP(I,J) = -C2eps * rho(I,J) * eps(I,J) / (k(I,J) + SMALL);
        Su(I,J) =  C1eps * (eps(I,J) / (k(I,J) + SMALL)) * 2.0 * mut(I,J) * E2(I,J);
        isWall  = false;

        % ---- Channel bottom wall ----
        if J == J_fluid_bottom
            y_P      = y(J_fluid_bottom) - y(J_fluid_bottom - 1);
            eps_wall = Cmu^0.75 * k(I,J)^1.5 / (kappa * y_P + SMALL);
            SP(I,J)  = -LARGE;
            Su(I,J)  =  LARGE * eps_wall;
            isWall   = true;
        end

        % ---- Channel top wall ----
        if J == J_fluid_top
            y_P      = y(J_fluid_top + 1) - y(J_fluid_top);
            eps_wall = Cmu^0.75 * k(I,J)^1.5 / (kappa * y_P + SMALL);
            SP(I,J)  = -LARGE;
            Su(I,J)  =  LARGE * eps_wall;
            isWall   = true;
        end

        % ---- Fin / solid cells ----
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
                    y_P      = 0.5 * AREAw;
                    eps_wall = Cmu^0.75 * k(I,J)^1.5 / (kappa * y_P + SMALL);
                    SP(I,J)  = -LARGE;
                    Su(I,J)  =  LARGE * eps_wall;
                    isWall   = true;
                end
            end
        end

        Su(I,J) = Su(I,J) * AREAw * AREAs;
        SP(I,J) = SP(I,J) * AREAw * AREAs;

        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);

        % Steady-state: no time derivative term
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
        b(I,J)  = Su(I,J);
    end
end
end