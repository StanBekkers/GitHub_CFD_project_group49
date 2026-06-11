function [] = Tcoeff()
% Purpose: To calculate the coefficients for the T equation.
% constants
global NPI NPJ YMAX
% variables
global x x_u y y_v T Gamma SP Su F_u F_v relax_T Istart Iend Jstart Jend ...
b aE aW aN aS aP heat_zone Cp
Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;
convect();
h_base_frac = 2/10;
% External air properties
h_air = 15.0; % Convective heat transfer coefficient to outside air [W/m^2K]
T_air = 293.15; % Ambient air temperature [K] (20 degrees C)
k_copper = 401.0; % Solid copper thermal conductivity [W/mK]
Dy = YMAX / NPJ; % Cell height
% Thermal resistance from the boundary cell center to the outside air
R_air = (0.5 * Dy) / k_copper + 1 / h_air;
for I = Istart:Iend
i = I;
for J = Jstart:Jend
j = J;
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
% The source terms
SP(I,J) = 0.;
Su(I,J) = 0.;
% --- Boundary Source Terms (Ambient Air Cooling & Chip Heat) ---
dx = x_u(i+1) - x_u(i);
% 1. Bottom wall boundary cells (J = 2)
if J == 2
is_in_chip_zone = false;
q_chip = 0;
for idx = 1:length(heat_zone)
if (x(I) >= heat_zone(idx).x_start && x(I) <= heat_zone(idx).x_end)
is_in_chip_zone = true;
q_chip = heat_zone(idx).q_wall;
end
end
if is_in_chip_zone
% Bottom is covered by chip: Inject heat flux into bottom of copper baseplate
Su(I,J) = Su(I,J) + (q_chip * dx) / Cp(I,J);
else
% Bottom is exposed to air: Convective cooling
SP(I,J) = SP(I,J) - (dx / (R_air * Cp(I,J)));
Su(I,J) = Su(I,J) + (T_air * dx) / (R_air * Cp(I,J));
end
end
% 2. Top wall boundary cells (J = NPJ + 1)
if J == NPJ+1
% Top is fully exposed to air: Convective cooling
SP(I,J) = SP(I,J) - (dx / (R_air * Cp(I,J)));
Su(I,J) = Su(I,J) + (T_air * dx) / (R_air * Cp(I,J));
end
% The coefficients (hybrid differencing scheme)
aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
aE(I,J) = max([-Fe, De - Fe/2, 0.]);
aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
% (Note: No manual zeroing of coefficients is needed here. Copper wall cells
% now naturally conduct heat in 2D to their neighbors via conjugate heat transfer)
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