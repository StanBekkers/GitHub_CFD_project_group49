function [] = Tcoeff()
% Purpose: To calculate the coefficients for the T equation in a 2D Planar (Top-Down) View.
% constants
global NPI NPJ YMAX XMAX LARGE SMALL t_fluid J_fluid_bottom J_fluid_top
% variables
global x x_u y y_v T Gamma SP Su F_u F_v relax_T Istart Iend Jstart Jend ...
b aE aW aN aS aP Cp l_base_frac P_core

Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;
convect();

% =========================================================================
% --- 2D PLANAR HEAT SOURCE PARAMETERS (Z-DIRECTION) ---
% =========================================================================
% Define chip active footprint centered over the baffle geometry

x_start = l_base_frac * XMAX;
x_end   = (1 - l_base_frac) * XMAX;

W_chip  = 0.015; % m (Width of the chip core, centered vertically)
y_start = 0.5 * YMAX - 0.5 * W_chip;
y_end   = 0.5 * YMAX + 0.5 * W_chip;

A_chip  = (x_end - x_start) * W_chip; % Total GPU area (m^2)
q_flux_core = P_core / A_chip;        % Heat flux from beneath (W/m^2)

Dx = XMAX / NPI;
Dy = YMAX / NPJ;
Area_cell = Dx * Dy; % Area of a single grid cell (m^2)
% =========================================================================

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;
        
        % --- PERFECTLY INSULATED CASING WALL DECOUPLING ---
        % Force the insulated outer casing cells to remain at the inlet temperature (293.15 K)
        % This prevents numerical heat creep in steady-state while keeping them adiabatic.
        if (J < J_fluid_bottom || J > J_fluid_top)
            SP(I,J) = -LARGE;
            Su(I,J) = LARGE * 293.15;
            
            aW(I,J) = 0.0;
            aE(I,J) = 0.0;
            aS(I,J) = 0.0;
            aN(I,J) = 0.0;
            aP(I,J) = LARGE;
            b(I,J)  = Su(I,J);
            continue;
        end

        % Geometrical parameters: Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j);
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i);
        AREAn = AREAs;
        
        % The convective mass flux
        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;
        
        % The transport by diffusion (harmonic mean handles the fluid-solid interface)
        Dw = ((Gamma(I-1,J)*Gamma(I,J))/(Gamma(I-1,J)*(x(I) - x_u(i)) ...
            + Gamma(I,J)*(x_u(i) - x(I-1))))*AREAw;
        De = ((Gamma(I,J)*Gamma(I+1,J))/(Gamma(I,J)*(x(I+1) - x_u(i+1)) ...
            + Gamma(I+1,J)*(x_u(i+1) - x(I))))*AREAe;
        Ds = ((Gamma(I,J-1)*Gamma(I,J))/(Gamma(I,J-1)*(y(J) - y_v(j)) ...
            + Gamma(I,J)*(y_v(j) - y(J-1))))*AREAs;
        Dn = ((Gamma(I,J)*Gamma(I,J+1))/(Gamma(I,J)*(y(J+1) - y_v(j+1)) ...
            + Gamma(I,J+1)*(y_v(j+1) - y(J))))*AREAn;
        
        % --- PLANAR HEAT SOURCE INTEGRATION ---
        SP(I,J) = 0.;
        Su(I,J) = 0.;
        
        % If the cell center lies inside the GPU chip's 2D plan footprint
        if (x(I) >= x_start && x(I) <= x_end && y(J) >= y_start && y(J) <= y_end)
            % Heat rate entering this specific cell from below (Watts)
            % Corrected: Divided by the physical depth of the channel (t_fluid) 
            % to scale the 2D heat source relative to the 2D mass flow rate.
            Q_cell = (q_flux_core / t_fluid) * Area_cell;
            
            % Add to source term, normalized by Cp to match the equations
            Su(I,J) = Q_cell / Cp(I,J);
        end
        
        % The coefficients (hybrid differencing scheme)
        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
        
        % eq. 8.31 without time dependent terms:
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
        
        % Setting the source term equal to b
        b(I,J) = Su(I,J);
        
        % Introducing relaxation
        aP(I,J) = aP(I,J) / relax_T;
        b(I,J) = b(I,J) + (1 - relax_T)*aP(I,J)*T(I,J);
    end
end
end