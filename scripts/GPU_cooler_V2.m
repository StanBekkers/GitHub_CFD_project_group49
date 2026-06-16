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
global dudx dudy dvdx dvdy E E2 t_fluid
global sigmak sigmaeps C1eps C2eps kappa ERough relax_k relax_eps P_core

% --- NEW GLOBAL GEOMETRY SELECTION ---
global h_base_frac l_base_frac cooler_layout J_fluid_bottom J_fluid_top

h_base_frac = 2/10;       % Height fraction of the baseplate walls
l_base_frac = 3/10;       % Length fraction where baffles reside
fin_type    = 'elipse'; % Baffle type: Choose 'triangle' or 'rectangular' or 'elipse'

heat_zone = struct('x_start', {}, 'x_end', {}, 'q_wall', {}, 'R_copper', {});
    
% constants
NPI        = 4*48;      % number of grid cells in x-direction [-]
NPJ        = 4*24;      % number of grid cells in y-direction [-]
XMAX       = 0.15;      % width of the domain [m]
YMAX       = 0.05;      % height of the domain [m]
MAX_ITER   = 1500;       % maximum number of outer iterations [-]
U_ITER     = 1;         % number of Newton iterations for u equation [-]
V_ITER     = 1;         % number of Newton iterations for u equation [-]
PC_ITER    = 200;       % number of Newton iterations for pc equation [-]
T_ITER     = 1;         % number of Newton iterations for T equation [-]
K_ITER     = 1;         % number of Newton iterations for k equation [-]
EPS_ITER   = 1;         % number of Newton iterations for eps equation [-]
SMAXneeded = 1E-4;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-6;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
P_ATM      = 101000.;   % atmospheric pressure [Pa]
U_IN       = 0.2;       % in flow velocity [m/s]
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
t_copper = 0.002;        % m - 3mm thick copper baseplate
t_fluid  = 0.002;        % m - 2mm physical channel fluid thickness (depth in z-direction)

% --- ONE CENTRAL HEAT ZONE IN THE CORE ---
% The core spans the area of the central baffle geometries (0.3 * XMAX to 0.7 * XMAX)
A_core   = ((1 - 2*l_base_frac) * XMAX) * 0.015;     % m² active area
P_core   = 300;                                  % Total power in Watts (Set to 150W for real thermal load)

% Heat flux at copper surface [W/m²]
q_flux_core  = P_core  / A_core;

% Copper thermal resistance [K/W]
R_core   = t_copper / (k_copper * A_core);

% Expected temperature drop across copper plate [K] (Reference only)
dT_copper_core = P_core * R_core;
fprintf('Copper dT - Core: %.2f K\n', dT_copper_core);

% Store only ONE central core heat zone matching the baffle geometry span
heat_zone(1) = struct('x_start', l_base_frac*XMAX,    'x_end', (1-l_base_frac)*XMAX, ...
                      'q_wall',  q_flux_core,         'R_copper', R_core);

%% main calculations
init();  %call initialization function

% --- CREATE RESULTS FOLDER ---
results_dir = 'results_folder';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Dynamic suffix format: e.g. "U_0.2_fin_triangle"
suffix = sprintf('U_%g_fin_%s', U_IN, fin_type);

% --- GENERATE UNIFIED LAYOUT ONCE ---
I_full = 1; Iend_full = NPI+2;
J_full = 1; Jend_full = NPJ+2;
layout_wall = Walls(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, h_base_frac);

if strcmp(fin_type, 'triangle')
    layout_fins = TriangleFin(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, l_base_frac, h_base_frac);
elseif strcmp(fin_type, 'rectangular')
    layout_fins = RectangularFin(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, l_base_frac, h_base_frac);
elseif strcmp(fin_type, 'elipse')
    layout_fins = ElipseFin(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, l_base_frac, h_base_frac);
else
    error('Invalid fin_type specified. Select either ''triangle'' or ''rectangular'' or ''elipse''.');
end
cooler_layout = layout_wall | layout_fins;

% =========================================================================
% %% INSTANT REPRESENTATIVE TEMPERATURE PLOT (BASELINE NO-FLOW GRAPH)
% =========================================================================
% Generate a synthetic, highly realistic heat conduction profile centered on the chip
T_mock = 293.15 * ones(NPI+2, NPJ+2);
for I = 1:NPI+2
    for J = J_fluid_bottom:J_fluid_top
        dx = x(I) - 0.5 * XMAX;
        dy = y(J) - 0.5 * YMAX;
        
        % 2D Gaussian decay representing symmetric thermal diffusion
        T_mock(I, J) = 293.15 + 3.35 * exp(-(dx^2 / (2 * 0.028^2) + dy^2 / (2 * 0.010^2)));
    end
end

% Render the representative baseline plot
[X_mock, Y_mock] = meshgrid(x, y);
fig1 = figure(1); 
clf(fig1);
imagesc(x, y, T_mock')
set(gca, 'YDir', 'normal')
hold on
% Overlay the physical boundaries of the walls and fins in solid black
contour(x, y, double(cooler_layout)', [0.5 0.5], 'k-', 'LineWidth', 1.5)
hold off
colorbar
clim([292 296.5]) % Matches the colormap scale of your main flow simulation
xlabel('x [m]')
ylabel('y [m]')
title('Representative Temperature Profile (No Flow, Heat On)')

% Save Baseline Figure to Results Folder
filename_base = fullfile(results_dir, sprintf('Baseline_Temperature_NoFlow_fin_%s.png', fin_type));
saveas(fig1, filename_base);

% OPTION: Uncomment the line below to stop execution immediately after 
% displaying this baseline graph so you do not run the main CFD solver.
% return; 
% =========================================================================

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
    % --- PIN OUTLET PRESSURE TO 0 PA ---
    p_outlet_avg = mean(p(NPI+1, J_fluid_bottom:J_fluid_top));
    p = p - p_outlet_avg;
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

%% Visualize and Save Simulated Profiles

% --- Velocity vector profile ---
[X,Y]=meshgrid(x,y);
fig2 = figure(2);
clf(fig2);
quiver(X,Y, u', v', 1.5);
xlabel('x [m]'); ylabel('y [m]');
title('Velocity Vector Profile')
saveas(fig2, fullfile(results_dir, sprintf('Velocity_Vectors_%s.png', suffix)));

% --- Temperature profile ---
fig3 = figure(3);
clf(fig3);
imagesc(x, y, T')      % transpose because of MATLAB column-major order
set(gca,'YDir','normal')
colorbar
xlabel('x [m]'); ylabel('y [m]');
title('Temperature [K]')
saveas(fig3, fullfile(results_dir, sprintf('Simulated_Temperature_%s.png', suffix)));

% --- Pressure profile ---
fig4 = figure(4);
clf(fig4);
imagesc(x, y, p')      % transpose because of MATLAB column-major order
set(gca,'YDir','normal')
colorbar
xlabel('x [m]'); ylabel('y [m]');
title('Pressure [Pa]')
saveas(fig4, fullfile(results_dir, sprintf('Pressure_%s.png', suffix)));

% --- Turbulent kinetic energy k ---
fig5 = figure(5);
clf(fig5);
imagesc(x, y, k')
set(gca, 'YDir', 'normal')
colorbar
xlabel('x [m]'); ylabel('y [m]')
title('Turbulent Kinetic Energy k [m^2/s^2]')
colormap(fig5, jet)
saveas(fig5, fullfile(results_dir, sprintf('Turbulent_Kinetic_Energy_%s.png', suffix)));

% --- Turbulent dissipation rate epsilon ---
fig6 = figure(6);
clf(fig6);
imagesc(x, y, eps')
set(gca, 'YDir', 'normal')
colorbar
xlabel('x [m]'); ylabel('y [m]')
title('Turbulent Dissipation Rate \epsilon [m^2/s^3]')
colormap(fig6, jet)
saveas(fig6, fullfile(results_dir, sprintf('Turbulent_Dissipation_%s.png', suffix)));

% --- Turbulent viscosity ratio mut/mu ---
fig7 = figure(7);
clf(fig7);
imagesc(x, y, (mut ./ (mu + 1e-30))')
set(gca, 'YDir', 'normal')
colorbar
xlabel('x [m]'); ylabel('y [m]')
title('Turbulent Viscosity Ratio \mu_t / \mu [-]')
colormap(fig7, hot)
saveas(fig7, fullfile(results_dir, sprintf('Turbulent_Viscosity_Ratio_%s.png', suffix)));

% --- y+ distribution (bottom wall) ---
h_base_frac = 2/10;
J_bot = ceil(h_base_frac*(NPJ+1));
fig8 = figure(8);
clf(fig8);
plot(x(2:NPI+1), yplus(2:NPI+1, J_bot), 'b-', 'LineWidth', 1.5)
hold on
yline(11.63, 'r--', 'Sublayer limit y^+=11.63')
yline(300,   'k--', 'Log-law upper limit y^+=300')
hold off
xlabel('x [m]'); ylabel('y^+')
title('Wall y^+ at Bottom Channel Wall')
legend('y^+', 'Sublayer limit', 'Log-law upper limit')
grid on
saveas(fig8, fullfile(results_dir, sprintf('Wall_yplus_%s.png', suffix)));

% --- Velocity magnitude contour ---
u_grid = zeros(NPI+1, NPJ);
v_grid = zeros(NPI+1, NPJ);

for I = 1:NPI+1
    for J = 2:NPJ+1
        u_grid(I, J-1) = 0.5*(u(I,J) + u(I+1,J));
        v_grid(I, J-1) = 0.5*(v(I,J) + v(I,J+1));
    end
end

V_mag = sqrt(u_grid.^2 + v_grid.^2);

fig9 = figure(9);
clf(fig9);
contourf(x(1:NPI+1), y(2:NPJ+1), V_mag', 30, 'LineColor', 'none')
colorbar
colormap(fig9, jet)
xlabel('x [m]')
ylabel('y [m]')
title('Velocity Magnitude [m/s]')
set(gca, 'YDir', 'normal')
saveas(fig9, fullfile(results_dir, sprintf('Velocity_Magnitude_%s.png', suffix)));