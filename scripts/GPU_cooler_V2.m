%% Solves: Steady, compressible convection-diffusion problems.
% Description:
% This program solves steady convection-diffusion problems using the simple algorithm 
% described in ch. 6.4 in "Computational Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. 
% Symbols and variables follow exactly the notations in this reference, and all
% equations cited are from this reference unless mentioned otherwise.

% Converted from C to Matlab by YTANG
% References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. Malalasekera, Longman Group Ltd, 1995

clear
close all
clc
global NPI NPJ XMAX YMAX LARGE SMALL U_IN

% variables
global x x_u y y_v u v pc p T rho mu Gamma Cp aP aE aW aN aS b d_u d_v  SMAX SAVG relax_rho 
global Q_chip x_chip_start x_chip_end y_chip_start y_chip_end heat_zone Ti Cmu
global k eps mut mueff uplus yplus yplus1 yplus2 tw k_old eps_old
global dudx dudy dvdx dvdy E E2
global sigmak sigmaeps C1eps C2eps kappa ERough relax_k relax_eps

heat_zone = struct('x_start', {}, 'x_end', {}, 'q_wall', {}, 'R_copper', {});
    
% constants
NPI        = 4*48;        % number of grid cells in x-direction [-]
NPJ        = 4*24;        % number of grid cells in y-direction [-]
XMAX       = 0.15;      % width of the domain [m]
YMAX       = 0.05;      % height of the domain [m]
MAX_ITER   = 500;      % maximum number of outer iterations [-]
U_ITER     = 1;         % number of Newton iterations for u equation [-]
V_ITER     = 1;         % number of Newton iterations for u equation [-]
PC_ITER    = 200;       % number of Newton iterations for pc equation [-]
T_ITER     = 1;         % number of Newton iterations for T equation [-]
K_ITER     = 1;         % number of Newton iterations for k equation [-]
EPS_ITER   = 1;         % number of Newton iterations for eps equation [-]
SMAXneeded = 1E-6;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-7;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
P_ATM      = 101000.;   % atmospheric pressure [Pa]
U_IN       = 0.2;      % in flow velocity [m/s]
NPRINT     = 1;         % number of iterations between printing output to screen

% k-epsilon constants (standard)
Cmu        = 0.09;
Ti         = 0.04;      % turbulence intensity at inlet [-]
SMALL      = 1E-30;     % arbitrary very small value to prevent division by zero
sigmak     = 1.0;       % turbulent Prandtl number for k
sigmaeps   = 1.3;       % turbulent Prandtl number for eps
C1eps      = 1.44;      % k-eps model constant
C2eps      = 1.92;      % k-eps model constant
kappa      = 0.4187;    % von Karman constant
ERough     = 9.793;     % roughness constant (smooth wall)

% Copper plate properties
k_copper = 401;          % W/m·K
t_copper = 0.003;        % m - 3mm thick copper baseplate

% Contact area per zone (2D: width * unit depth of 1m)
A_left   = (0.3 * XMAX) * 1;     % m²
A_core   = (0.4 * XMAX) * 1;     % m²
A_right  = (0.3 * XMAX) * 1;     % m²

% Total power per zone
P_left   = 75;           % W - VRM/memory left
P_core   = 350;          % W - GPU core
P_right  = 75;           % W - VRM/memory right

% Heat flux at copper surface [W/m²]
% This is what actually enters the fluid after passing through copper
q_flux_left  = P_left  / A_left;     % ~174  W/m²
q_flux_core  = P_core  / A_core;     % ~781  W/m²
q_flux_right = P_right / A_right;    % ~174  W/m²

% Copper thermal resistance per zone [K/W]
% R = t / (k * A) - tells you temperature drop across copper
R_left   = t_copper / (k_copper * A_left);
R_core   = t_copper / (k_copper * A_core);
R_right  = t_copper / (k_copper * A_right);

% Expected temperature drop across copper plate [K]
% Just for reference/sanity check - printed but not used in solver
dT_copper_left  = P_left  * R_left;
dT_copper_core  = P_core  * R_core;
dT_copper_right = P_right * R_right;

fprintf('Copper dT - Left: %.3f K, Core: %.2f K, Right: %.3f K\n', ...
        dT_copper_left, dT_copper_core, dT_copper_right);

% Store zones - note q_wall not q_vol, entering as boundary flux
heat_zone(1) = struct('x_start', 0.0,         'x_end', 0.3*XMAX, ...
                      'q_wall',  q_flux_left,  'R_copper', R_left);

heat_zone(2) = struct('x_start', 0.3*XMAX,    'x_end', 0.7*XMAX, ...
                      'q_wall',  q_flux_core,  'R_copper', R_core);

heat_zone(3) = struct('x_start', 0.7*XMAX,    'x_end', 1.0*XMAX, ...
                      'q_wall',  q_flux_right, 'R_copper', R_right);
%% main calculations
init();  %call initialization function

iter = 1;
% outer iteration loop
while (iter <= MAX_ITER && (SMAX > SMAXneeded || SAVG > SAVGneeded))
    
    bound(); %call boundary function
    % Turbulence: solve k then eps, clip to physical bounds
    derivatives();
    kcoeff();
    for iter_k = 1:K_ITER
        k = solve(k, b, aE, aW, aN, aS, aP);
    end
    k = max(k, 1e-10);   % prevent negative k

    epscoeff();
    for iter_eps = 1:EPS_ITER
        eps = solve(eps, b, aE, aW, aN, aS, aP);
    end
    eps = max(eps, 1e-10); % prevent negative eps

    rho(:,:) = 1000.0;
    mu(1:NPI+2, 2:NPJ+1) = 1.0E-3;
    viscosity();   % now builds mueff and Gamma on fresh mu
    
    ucoeff(); %call ucoeffe.m function to calculate the coefficients for u function
    for iter_u = 1:U_ITER
        u = solve(u, b, aE, aW, aN, aS, aP); %solve u function
    end
    
    vcoeff(); %call vcoeffe.m function to calculate the coefficients for v function
    for iter_v = 1:V_ITER
        v = solve(v, b, aE, aW, aN, aS, aP); %solve v function
    end
    
    bound(); %apply boundary conditions again
    
    pccoeff(); %call pccoeffe.m function to calculate the coefficients for p function
    for iter_pc = 1:PC_ITER
        pc = solve(pc, b, aE, aW, aN, aS, aP); %solve p function
    end
    
    velcorr(); % Correct pressure and velocity

    Tcoeff(); %call Tcoeffe.m function to calculate the coefficients for T function
    for iter_T = 1:T_ITER
        T = solve(T, b, aE, aW, aN, aS, aP); %solve T function
    end
   
    % begin: printConv(iter)========================================================================
    % print temporary results
    if iter == 1
        fprintf ('Iter.\t d_u/u\t\t d_v/v\t\t SMAX\t\t SAVG\n');
    end
    if mod(iter,NPRINT) == 0
        I = round((NPI+1)/2);
        J = round((NPJ+1)/2);
        du = d_u(I,J)*(pc(I-1,J) - pc(I,J));
        dv = d_v(I,J)*(pc(I,J-1) - pc(I,J));
        fprintf ('%3d\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n', iter,du/u(I,J), dv/v(I,J), SMAX, SAVG);
    end
    % end of print temporaty results=================================================================
    
    % increase interation number
    iter = iter + 1;   
end
%% begin: output()
% print out results in files
fp   = fopen('output.dat','w');
str  = fopen('str.dat','w');
velu = fopen('velu.dat','w');
velv = fopen('velv.dat','w');

for I = 1:NPI+1
    i = I;
    for J = 2:NPJ+1
        j = J;
        ugrid = 0.5*(u(i,J)+u(i+1,J));
        vgrid = 0.5*(v(I,j)+v(I,j+1));
        fprintf(fp, '%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n',...
            x(I), y(J), ugrid, vgrid, p(I,J), T(I,J), rho(I,J), mu(I,J), Gamma(I,J));
    end
    fprintf(fp, '\n');
end
fclose(fp);

for I = 1:NPI+1
    i = I;
    for J = 2:NPJ+1
        j = J;
        stream = -(u(i,J+1)-u(i,J))/(y(J+1)-y(J))+(v(I+1,j)-v(I,j))/(x(I+1)-x(I));
        fprintf(str, '%10.2e\t%10.2e\t%10.5e\n',x_u(i), y_v(j), stream);
        fprintf(velu,'%10.2e\t%10.2e\t%10.5e\n',x_u(i), y(J)  , u(i,J));
        fprintf(velv,'%10.2e\t%10.2e\t%10.5e\n',x(I)  , y_v(j), v(I,j));
    end
    fprintf(str, '\n');
    fprintf(velu,'\n');
    fprintf(velv,'\n');
end

fclose(str);
fclose(velu);
fclose(velv);
%% visulize the velocity profile
[X,Y]=meshgrid(x,y);
quiver(X,Y, u', v',1.5);

%% visulize the temperature profile
[X,Y] = meshgrid(x,y);

figure
imagesc(x, y, T')      % transpose because of MATLAB column-major order
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Temperature')

%% visulize the pressure profile
[X,Y] = meshgrid(x,y);

figure
imagesc(x, y, p')      % transpose because of MATLAB column-major order
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Pressure')

%% Visualise turbulence quantities

% --- Turbulent kinetic energy k ---
figure
imagesc(x, y, k')
set(gca, 'YDir', 'normal')
colorbar
xlabel('x [m]'); ylabel('y [m]')
title('Turbulent Kinetic Energy k [m^2/s^2]')
colormap(jet)

% --- Turbulent dissipation rate epsilon ---
figure
imagesc(x, y, eps')
set(gca, 'YDir', 'normal')
colorbar
xlabel('x [m]'); ylabel('y [m]')
title('Turbulent Dissipation Rate \epsilon [m^2/s^3]')
colormap(jet)

% --- Turbulent viscosity ratio mut/mu ---
figure
imagesc(x, y, (mut ./ (mu + 1e-30))')
set(gca, 'YDir', 'normal')
colorbar
xlabel('x [m]'); ylabel('y [m]')
title('Turbulent Viscosity Ratio \mu_t / \mu [-]')
colormap(hot)

% --- y+ distribution (bottom wall) ---
h_base_frac = 2/10;
J_bot = ceil(h_base_frac*(NPJ+1));
figure
plot(x(2:NPI+1), yplus(2:NPI+1, J_bot), 'b-', 'LineWidth', 1.5)
hold on
yline(11.63, 'r--', 'Sublayer limit y^+=11.63')
yline(300,   'k--', 'Log-law upper limit y^+=300')
xlabel('x [m]'); ylabel('y^+')
title('Wall y^+ at Bottom Channel Wall')
legend('y^+', 'Sublayer limit', 'Log-law upper limit')
grid on

%% Visualise velocity magnitude contour
u_grid = zeros(NPI+1, NPJ);
v_grid = zeros(NPI+1, NPJ);

for I = 1:NPI+1
    for J = 2:NPJ+1
        u_grid(I, J-1) = 0.5*(u(I,J) + u(I+1,J));
        v_grid(I, J-1) = 0.5*(v(I,J) + v(I,J+1));
    end
end

V_mag = sqrt(u_grid.^2 + v_grid.^2);

figure
contourf(x(1:NPI+1), y(2:NPJ+1), V_mag', 30, 'LineColor', 'none')
colorbar
colormap(jet)
xlabel('x [m]')
ylabel('y [m]')
title('Velocity Magnitude [m/s]')
set(gca, 'YDir', 'normal')




