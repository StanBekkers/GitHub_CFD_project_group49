function [] = calculateuplus()
% Purpose: Calculate uplus and tw at top and bottom channel walls.
% Uses the robust k-based y* formulation (Launder-Spalding) to maintain numerical
% stability near separation, stagnation, and recirculation zones.
%
% Water properties used (consistent with GPU_Cooler_V1.m):
%   rho = 1000 kg/m^3, mu = 1e-3 Pa.s

% constants
global NPI NPJ Cmu kappa ERough SMALL
% variables
global y rho k tw yplus yplus1 yplus2 uplus mu u

% --- CONNECT TO GEOMETRY GLOBALS ---
global h_base_frac J_fluid_bottom J_fluid_top

for I = 1:NPI+1
    i = I;

    % ==========================================
    % --- Bottom channel wall ---
    % ==========================================
    y_P = 0.5*(y(J_fluid_bottom) - y(J_fluid_bottom - 1));
    u_P = 0.5*(u(i, J_fluid_bottom) + u(i+1, J_fluid_bottom));

    % Use k-derived wall scale y* to prevent singularity when u_P or tw goes to zero
    k_local = max(k(I, J_fluid_bottom), 0.0);
    yplus_val = rho(I, J_fluid_bottom) * Cmu^0.25 * sqrt(k_local + SMALL) * y_P / (mu(I, J_fluid_bottom) + SMALL);
    
    yplus(I, J_fluid_bottom)  = yplus_val;
    yplus1(I, J_fluid_bottom) = yplus_val; % Keep both synchronized

    if yplus_val < 11.63
        % Viscous sublayer: laminar stress relation
        uplus(I, J_fluid_bottom) = max(yplus_val, SMALL);
        tw(I, J_fluid_bottom)    = mu(I, J_fluid_bottom) * abs(u_P) / (y_P + SMALL);
    else
        % Log-law region: turbulent shear stress
        uplus(I, J_fluid_bottom) = log(ERough * max(yplus_val, 1.0)) / kappa;
        tw(I, J_fluid_bottom)    = rho(I, J_fluid_bottom) * Cmu^0.25 * sqrt(k_local + SMALL) ...
                                   * abs(u_P) / (uplus(I, J_fluid_bottom) + SMALL);
    end

    % ==========================================
    % --- Top channel wall ---
    % ==========================================
    y_P = 0.5*(y(J_fluid_top + 1) - y(J_fluid_top));
    u_P = 0.5*(u(i, J_fluid_top) + u(i+1, J_fluid_top));

    k_local = max(k(I, J_fluid_top), 0.0);
    yplus_val = rho(I, J_fluid_top) * Cmu^0.25 * sqrt(k_local + SMALL) * y_P / (mu(I, J_fluid_top) + SMALL);
    
    yplus(I, J_fluid_top)  = yplus_val;
    yplus2(I, J_fluid_top) = yplus_val; % Keep both synchronized

    if yplus_val < 11.63
        % Viscous sublayer
        uplus(I, J_fluid_top) = max(yplus_val, SMALL);
        tw(I, J_fluid_top)    = mu(I, J_fluid_top) * abs(u_P) / (y_P + SMALL);
    else
        % Log-law region
        uplus(I, J_fluid_top) = log(ERough * max(yplus_val, 1.0)) / kappa;
        tw(I, J_fluid_top)    = rho(I, J_fluid_top) * Cmu^0.25 * sqrt(k_local + SMALL) ...
                                 * abs(u_P) / (uplus(I, J_fluid_top) + SMALL);
    end
end
end
