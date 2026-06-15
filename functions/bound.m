function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN YMAX heat_zone x
% variables
global  u v T m_in m_out k eps y_v F_u Ti Cmu

% --- CONNECT TO GEOMETRY GLOBALS ---
global h_base_frac J_fluid_bottom J_fluid_top

% Inlet: fluid cells get 293.15 K
T(1, J_fluid_bottom:J_fluid_top) = 293.15;
% Inlet: wall cells also initialised to 293.15 K (copper will equilibrate)
T(1, 1:J_fluid_bottom-1) = 293.15;
T(1, J_fluid_top+1:NPJ+2) = 293.15;

% Setting the velocity at inlet
u(2, J_fluid_bottom:J_fluid_top) = U_IN;

% --- Physical Boundary Node Temperature Updates ---
k_copper = 401.0;
h_air = 15.0;
T_air = 293.15;
Dy = YMAX / NPJ;
dy_half = 0.5 * Dy;

for I = 1:NPI+2
    % --- TOP-DOWN VIEW MODIFICATION ---
    % Both outer casing walls (J = 1 at y = 0, and J = NPJ + 2 at y = YMAX)
    % are now simply exposed to ambient air.
    
    % 1. Bottom casing wall surface node
    T(I, 1) = (T(I, 2) * (k_copper / dy_half) + h_air * T_air) / ((k_copper / dy_half) + h_air);
    
    % 2. Top casing wall surface node
    T(I, NPJ+2) = (T(I, NPJ+1) * (k_copper / dy_half) + h_air * T_air) / ((k_copper / dy_half) + h_air);
end

% --- CORRECTED TURBULENT LENGTH SCALE TYPO ---
% L_t based on hydraulic diameter of the channel inlet region (uses YMAX instead of U_IN)
L_t     = 0.07 * (YMAX * 4/6);           % 0.07 * D_h  (D_h ~ channel height)
k_in    = 1.5 * (Ti * U_IN)^2;
eps_in  = Cmu^(3/4) * k_in^(3/2) / L_t;
 
k  (1, J_fluid_bottom:J_fluid_top) = k_in;
k  (2, J_fluid_bottom:J_fluid_top) = k_in;
eps(1, J_fluid_bottom:J_fluid_top) = eps_in;
eps(2, J_fluid_bottom:J_fluid_top) = eps_in;
 
% Outlet: zero-gradient for k and eps (Neumann)
k  (NPI+2,2:NPJ+1) = k  (NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
 
% Walls: k = 0 (Dirichlet) at solid walls — wall functions handle eps there
k(1:NPI+2,1:J_fluid_bottom-1)       = 0.;
k(1:NPI+2,J_fluid_top+1:NPJ+2)      = 0.;

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
% Reset outlet solid temperatures to air temperature (no heat source there)
T(NPI+2, 1:J_fluid_bottom-1)   = T(NPI+1, 1:J_fluid_bottom-1);
T(NPI+2, J_fluid_top+1:NPJ+2) = T(NPI+1, J_fluid_top+1:NPJ+2);
end
