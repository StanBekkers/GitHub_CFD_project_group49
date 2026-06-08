function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN
% variables
global  u v T m_in m_out k eps y_v F_u Ti Cmu

global  u v T m_in m_out y_v F_u

% Fixed temperature in Kelvin of the incoming fluid (293.15 K = 20°C)
T(1,ceil((NPJ+1)/6):ceil(5*(NPJ+1)/6)) = 293.15; 

% Setting the velocity at inlet
u(2,ceil((NPJ+1)/6):ceil(5*(NPJ+1)/6)) = U_IN;

% Adiabatic walls: Set solid wall temp equal to adjacent active fluid cell temp
% 1. Lower wall: Set all lower wall rows equal to the first active fluid row
J_fluid_bottom = ceil((NPJ+1)/6);
for J_wall = 1 : J_fluid_bottom - 1
    T(1:NPI+2, J_wall) = T(1:NPI+2, J_fluid_bottom);
end

% 2. Upper wall: Set all upper wall rows equal to the last active fluid row
J_fluid_top = ceil(5*(NPJ+1)/6);
for J_wall = J_fluid_top + 1 : NPJ+2
    T(1:NPI+2, J_wall) = T(1:NPI+2, J_fluid_top);
end

% L_t based on hydraulic diameter of the channel inlet region
% Define inlet and outlet k and eps
L_t     = 0.07 * (0.2 * 4/6);           % 0.07 * D_h  (D_h ~ channel height)
k_in    = 1.5 * (Ti * U_IN)^2;
eps_in  = Cmu^(3/4) * k_in^(3/2) / L_t;
 
k  (1, ceil((NPJ+1)/6):ceil(5*(NPJ+1)/6)) = k_in;
k  (2, ceil((NPJ+1)/6):ceil(5*(NPJ+1)/6)) = k_in;
eps(1, ceil((NPJ+1)/6):ceil(5*(NPJ+1)/6)) = eps_in;
eps(2, ceil((NPJ+1)/6):ceil(5*(NPJ+1)/6)) = eps_in;
 
% Outlet: zero-gradient for k and eps (Neumann)
k  (NPI+2,2:NPJ+1) = k  (NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
 
% Walls: k = 0 (Dirichlet) at solid walls — wall functions handle eps there
k(1:NPI+2,1:ceil((NPJ+1)/6)-1)       = 0.;
k(1:NPI+2,ceil(5*(NPJ+1)/6+1):NPJ+2) = 0.;

% begin: globcont(): Velocity and temperature gradient at outlet = zero:
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
convect();

m_in = 0.;
m_out = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out = m_out + F_u(NPI+1,J)*AREAw;
end
% end: globcont():

% corection varibles: Correction factor m_in/m_out is used to satisfy global continuity
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1)*m_in/m_out;
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
T(NPI+2,2:NPJ+1) = T(NPI+1,2:NPJ+1);
end

