function [] = init()
% Purpose: To initilise all parameters.

% constants
global NPI NPJ LARGE U_IN XMAX YMAX
% variables
global x x_u y y_v u v pc p T rho mu Gamma Cp b SP Su d_u d_v omega SMAX SAVG ...
    m_in m_out relax_u relax_v relax_pc relax_T relax_rho aP aE aW aN aS F_u F_v ...
    k eps uplus yplus yplus1 yplus2 tw ...
    u_old v_old pc_old Ti T_old k_old eps_old Cmu mut mueff
    

% begin: memalloc()=======================================================
% allocate memory for variables
x   = zeros(1,NPI+2);
x_u = zeros(1,NPI+2);
y   = zeros(1,NPJ+2);
y_v = zeros(1,NPJ+2);

u   = zeros(NPI+2,NPJ+2);
v   = zeros(NPI+2,NPJ+2);
pc  = zeros(NPI+2,NPJ+2);
p   = zeros(NPI+2,NPJ+2);
T   = zeros(NPI+2,NPJ+2);
rho = zeros(NPI+2,NPJ+2);
mu  = zeros(NPI+2,NPJ+2);
Gamma = zeros(NPI+2,NPJ+2);
Cp  = zeros(NPI+2,NPJ+2);

aP  = zeros(NPI+2,NPJ+2);
aE  = zeros(NPI+2,NPJ+2);
aW  = zeros(NPI+2,NPJ+2);
aN  = zeros(NPI+2,NPJ+2);
aS  = zeros(NPI+2,NPJ+2);
b   = zeros(NPI+2,NPJ+2);

SP  = zeros(NPI+2,NPJ+2);
Su  = zeros(NPI+2,NPJ+2);

F_u = zeros(NPI+2,NPJ+2);
F_v = zeros(NPI+2,NPJ+2);

d_u = zeros(NPI+2,NPJ+2);
d_v = zeros(NPI+2,NPJ+2);

% Turbulence variables
k      = zeros(NPI+2,NPJ+2);
eps    = zeros(NPI+2,NPJ+2);
uplus  = zeros(NPI+2,NPJ+2);
yplus  = zeros(NPI+2,NPJ+2);
yplus1 = zeros(NPI+2,NPJ+2);
yplus2 = zeros(NPI+2,NPJ+2);
tw     = zeros(NPI+2,NPJ+2);

% Time-level storage arrays
u_old  = zeros(NPI+2,NPJ+2);
v_old  = zeros(NPI+2,NPJ+2);
pc_old = zeros(NPI+2,NPJ+2);
T_old  = zeros(NPI+2,NPJ+2);
k_old  = zeros(NPI+2,NPJ+2);
eps_old= zeros(NPI+2,NPJ+2);

% end of memory allocation=================================================

% begin: grid()===========================================================
% Purpose: Defining the geometrical variables. See fig. 6.2-6.4 in ref. 1
% Length of volume element
Dx = XMAX/NPI;
Dy = YMAX/NPJ;

% Length variable for the scalar points in the x direction
x(1) = 0.;
x(2) = 0.5*Dx;
for I = 3:NPI+1
    x(I) = x(I-1) + Dx;
end
x(NPI+2) = x(NPI+1) + 0.5*Dx;

% Length variable for the scalar points T(I,J) in the y direction
y(1) = 0.;
y(2) = 0.5*Dy;
for J = 3:NPJ+1
    y(J) = y(J-1) + Dy;
end
y(NPJ+2) = y(NPJ+1) + 0.5*Dy;

% Length variable for the velocity components u(i,J) in the x direction
x_u(1) = 0.;
x_u(2) = 0.;
for i = 3:NPI+2
    x_u(i) = x_u(i-1) + Dx;
end

% Length variable for the velocity components v(I,j) in the y direction 
y_v(1) = 0.;
y_v(2) = 0.;
for j = 3:NPJ+2
    y_v(j) = y_v(j-1) + Dy;
end
% end of grid setting======================================================

% begin: init()===========================================================
% Initialising all other variables
omega = 1.0; % Over-relaxation factor for SOR solver

% Initialize convergence parameters at large values
SMAX = LARGE;
SAVG = LARGE;

m_in  = 1.;
m_out = 1.;

u(:,:)   = 0.;    % Velocity in x-direction
v(:,:)   = 0.;    % Velocity in y-direction
p(:,:)   = 0.;    % Relative pressure
pc(:,:)  = 0.;    % Pressure correction (equivalet to p' in ref. 1).
T(:,:)   = 273.;  % Temperature
rho(:,:) = 1.0;   % Density
mu(:,:)  = 2.E-5; % Viscosity
Cp(:,:)  = 1013.; % J/(K*kg) Heat capacity - aSAVGed constant for this problem
Gamma    = 0.0315./Cp; % Thermal conductivity
d_u(:,:) = 0.;    % Variable d(i,j) to calculate pc defined in 6.23
d_v(:,:) = 0.;    % Variable d(i,j) to calculate pc defined in 6.23
b(:,:)   = 0.;	  % The general constant
SP(:,:)  = 0.;    % Source term
Su(:,:)  = 0.;	  % Source term

% Hydraulic diameter ~ channel height (middle region)
L_t   = 0.07 * (YMAX * 4/6);   % turbulent length scale (0.07 * D_h)
k_init  = 1.5 * (Ti * U_IN)^2;
eps_init = Cmu^(3/4) * k_init^(3/2) / L_t;
 
k(:,:)    = k_init;
eps(:,:)  = eps_init;
k_old(:,:)   = k_init;
eps_old(:,:) = eps_init;
 
% Initialise uplus (log-law starting guess, avoid zero-divide in wall functions)
uplus(:,:) = 11.63;
yplus(:,:) = 11.63;
yplus1(:,:) = 11.63;
yplus2(:,:) = 11.63;
 
% Initialise effective viscosity with turbulent contribution
mut(:,:)   = rho(1,1)*Cmu*k_init^2 / eps_init;
mueff(:,:) = mu(1,1) + mut(1,1);
 
u(NPI+1,2:NPJ+1) = 0.5*U_IN;
% Important to avoid crash!! Othervise m_out calculated in subroutine globcont
% would be zero at first iteration=>m_in/m_out =INF

% Setting the relaxation parameters
relax_u   = 0.8;            % See eq. 6.36
relax_v   = relax_u;        % See eq. 6.37
relax_pc  = 1.1 - relax_u;  % See eq. 6.33
relax_T   = 1.0;            % Relaxation factor for temperature
relax_rho = 0.1;            % Relaxation factor for density
% end of initilization=====================================================
end

