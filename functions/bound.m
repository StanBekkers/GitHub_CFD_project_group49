function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN YMAX heat_zone x
% variables
global  u v T m_in m_out k eps y_v F_u Ti Cmu

h_base_frac = 2/10;

% Fixed temperature in Kelvin of the incoming fluid (293.15 K = 20°C)
T(1, ceil(h_base_frac*(NPJ+1)):ceil((1-h_base_frac)*(NPJ+1))) = 293.15; 

% Setting the velocity at inlet
u(2, ceil(h_base_frac*(NPJ+1)):ceil((1-h_base_frac)*(NPJ+1))) = U_IN;

% --- Physical Boundary Node Temperature Updates ---
k_copper = 401.0;
h_air = 15.0;
T_air = 293.15;
Dy = YMAX / NPJ;
dy_half = 0.5 * Dy;

for I = 1:NPI+2
    % 1. Bottom boundary surface node (J = 1 at y = 0)
    is_in_chip_zone = false;
    q_chip = 0;
    for idx = 1:length(heat_zone)
        if (x(I) >= heat_zone(idx).x_start && x(I) <= heat_zone(idx).x_end)
            is_in_chip_zone = true;
            q_chip = heat_zone(idx).q_wall;
        end
    end
    
    if is_in_chip_zone
        % Heat enters from chip: boundary surface is warmer than adjacent solid center
        T(I, 1) = T(I, 2) + (q_chip * dy_half) / k_copper;
    else
        % Exposed to air: boundary surface temperature balanced by convection and conduction
        T(I, 1) = (T(I, 2) * (k_copper / dy_half) + h_air * T_air) / ((k_copper / dy_half) + h_air);
    end
    
    % 2. Top boundary surface node (J = NPJ + 2 at y = YMAX)
    % Exposed to air: balanced by convection and conduction
    T(I, NPJ+2) = (T(I, NPJ+1) * (k_copper / dy_half) + h_air * T_air) / ((k_copper / dy_half) + h_air);
end

% L_t based on hydraulic diameter of the channel inlet region
% Define inlet and outlet k and eps
L_t     = 0.07 * (0.2 * 4/6);           % 0.07 * D_h  (D_h ~ channel height)
k_in    = 1.5 * (Ti * U_IN)^2;
eps_in  = Cmu^(3/4) * k_in^(3/2) / L_t;
 
k  (1, ceil(h_base_frac*(NPJ+1)):ceil((1-h_base_frac)*(NPJ+1))) = k_in;
k  (2, ceil(h_base_frac*(NPJ+1)):ceil((1-h_base_frac)*(NPJ+1))) = k_in;
eps(1, ceil(h_base_frac*(NPJ+1)):ceil((1-h_base_frac)*(NPJ+1))) = eps_in;
eps(2, ceil(h_base_frac*(NPJ+1)):ceil((1-h_base_frac)*(NPJ+1))) = eps_in;
 
% Outlet: zero-gradient for k and eps (Neumann)
k  (NPI+2,2:NPJ+1) = k  (NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
 
% Walls: k = 0 (Dirichlet) at solid walls — wall functions handle eps there
k(1:NPI+2,1:ceil(h_base_frac*(NPJ+1))-1)       = 0.;
k(1:NPI+2,ceil((1-h_base_frac)*(NPJ+1))+1:NPJ+2) = 0.;

% begin: globcont()
convect();

m_in = 0.;
m_out = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); 
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out = m_out + F_u(NPI+1,J)*AREAw;
end
% end: globcont()

% correction variables
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1)*m_in/m_out;
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
T(NPI+2,2:NPJ+1) = T(NPI+1,2:NPJ+1);
end