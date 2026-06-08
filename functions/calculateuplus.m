function [] = calculateuplus()
% Purpose: Calculate uplus and tw at top and bottom channel walls.
% Uses the log-law wall function. These values feed into kcoeff and epscoeff.
%
% Water properties used (consistent with GPU_Cooler_V1.m):
%   rho = 1000 kg/m^3, mu = 1e-3 Pa.s

% constants
global NPI NPJ Cmu kappa ERough SMALL
% variables
global y rho k tw yplus yplus1 yplus2 uplus mu u

% Geometry — must match ucoeff.m / vcoeff.m
h_base_frac   = 2/10;
J_fluid_bottom = ceil(h_base_frac*(NPJ+1));   % first fluid J (bottom)
J_fluid_top    = ceil((1-h_base_frac)*(NPJ+1)); % last  fluid J (top)

for I = 1:NPI+1
    i = I;

    % --- Bottom channel wall ---
    % y_P = distance from wall face to first fluid cell centre
    y_P = y(J_fluid_bottom) - y(J_fluid_bottom - 1);
    u_P = 0.5*(u(i, J_fluid_bottom) + u(i+1, J_fluid_bottom));

    if yplus1(I, J_fluid_bottom) < 11.63
        % Viscous sublayer: laminar stress
        tw(I, J_fluid_bottom)     = mu(I, J_fluid_bottom) * abs(u_P) / (y_P + SMALL);
        u_tau                     = sqrt(abs(tw(I, J_fluid_bottom)) / (rho(I, J_fluid_bottom) + SMALL));
        yplus1(I, J_fluid_bottom) = rho(I, J_fluid_bottom) * u_tau * y_P / (mu(I, J_fluid_bottom) + SMALL);
        yplus(I, J_fluid_bottom)  = yplus1(I, J_fluid_bottom);
        uplus(I, J_fluid_bottom)  = max(yplus(I, J_fluid_bottom), SMALL);
    else
        % Log-law region
        tw(I, J_fluid_bottom)     = rho(I, J_fluid_bottom) * Cmu^0.25 * sqrt(k(I, J_fluid_bottom) + SMALL) ...
                                    * abs(u_P) / (uplus(I, J_fluid_bottom) + SMALL);
        u_tau                     = sqrt(abs(tw(I, J_fluid_bottom)) / (rho(I, J_fluid_bottom) + SMALL));
        yplus1(I, J_fluid_bottom) = rho(I, J_fluid_bottom) * u_tau * y_P / (mu(I, J_fluid_bottom) + SMALL);
        yplus(I, J_fluid_bottom)  = yplus1(I, J_fluid_bottom);
        uplus(I, J_fluid_bottom)  = log(ERough * max(yplus(I, J_fluid_bottom), 1.0)) / kappa;
    end

    % --- Top channel wall ---
    y_P = y(J_fluid_top + 1) - y(J_fluid_top);
    u_P = 0.5*(u(i, J_fluid_top) + u(i+1, J_fluid_top));

    if yplus2(I, J_fluid_top) < 11.63
        tw(I, J_fluid_top)     = mu(I, J_fluid_top) * abs(u_P) / (y_P + SMALL);
        u_tau                  = sqrt(abs(tw(I, J_fluid_top)) / (rho(I, J_fluid_top) + SMALL));
        yplus2(I, J_fluid_top) = rho(I, J_fluid_top) * u_tau * y_P / (mu(I, J_fluid_top) + SMALL);
        yplus(I, J_fluid_top)  = yplus2(I, J_fluid_top);
        uplus(I, J_fluid_top)  = max(yplus(I, J_fluid_top), SMALL);
    else
        tw(I, J_fluid_top)     = rho(I, J_fluid_top) * Cmu^0.25 * sqrt(k(I, J_fluid_top) + SMALL) ...
                                 * abs(u_P) / (uplus(I, J_fluid_top) + SMALL);
        u_tau                  = sqrt(abs(tw(I, J_fluid_top)) / (rho(I, J_fluid_top) + SMALL));
        yplus2(I, J_fluid_top) = rho(I, J_fluid_top) * u_tau * y_P / (mu(I, J_fluid_top) + SMALL);
        yplus(I, J_fluid_top)  = yplus2(I, J_fluid_top);
        uplus(I, J_fluid_top)  = log(ERough * max(yplus(I, J_fluid_top), 1.0)) / kappa;
    end
end
end