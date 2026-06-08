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
%% declare all constants and variables
% global contants
global NPI NPJ XMAX YMAX LARGE U_IN Cmu Ti

% variables
global x x_u y y_v u v pc p T rho mu Gamma Cp aP aE aW aN aS b d_u d_v  SMAX SAVG relax_rho 
global Q_chip x_chip_start x_chip_end y_chip_start y_chip_end heat_zone

heat_zone = struct('x_start', {}, 'x_end', {}, 'q_wall', {}, 'R_copper', {});
    
% constantsa
NPI        = 2*48;        % number of grid cells in x-direction [-]
NPJ        = 2*24;        % number of grid cells in y-direction [-]
XMAX       = 0.96;      % width of the domain [m]
YMAX       = 0.12;      % height of the domain [m]
MAX_ITER   = 1000;      % maximum number of outer iterations [-]
U_ITER     = 1;         % number of Newton iterations for u equation [-]
V_ITER     = 1;         % number of Newton iterations for u equation [-]
PC_ITER    = 200;       % number of Newton iterations for pc equation [-]
T_ITER     = 1;         % number of Newton iterations for T equation [-]
SMAXneeded = 1E-7;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-8;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
P_ATM      = 101000.;   % atmospheric pressure [Pa]
U_IN       = 0.02;      % in flow velocity [m/s]
NPRINT     = 1;         % number of iterations between printing output to screen

% k-epsilon constants (standard)
Cmu        = 0.09;
Ti         = 0.04;      % turbulence intensity at inlet [-]

% Copper plate properties
k_copper = 401;          % W/m·K
t_copper = 0.003;        % m - 3mm thick copper baseplate

% Contact area per zone (2D: width * unit depth of 1m)
A_left   = (0.3 * XMAX) * 1;     % m²
A_core   = (0.4 * XMAX) * 1;     % m²
A_right  = (0.3 * XMAX) * 1;     % m²

% Total power per zone
P_left   = 25;           % W - VRM/memory left
P_core   = 150;          % W - GPU core
P_right  = 25;           % W - VRM/memory right

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
while (iter <= MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded)
    
    bound(); %call boundary function
    
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
    
    % begin: density()==============================================================================
    % Purpose: Calculate the density rho(I, J) in the fluid as a function of the ideal gas law. 
    % Note: rho at the walls are not needed in this case, and therefore not calculated.    
    for I = 1:NPI+2
        for J = 2:NPJ+1
            if (I == 1) % Since p(1, J) doesn't exist, we set:
                rho(I,J) = (1 - relax_rho)*rho(I,J) + relax_rho*(p(I+1,J) + P_ATM)/(287.*T(I,J));
            else
                rho(I,J) = (1 - relax_rho)*rho(I,J) + relax_rho*(p(I,J) + P_ATM)/(287.*T(I,J));
            end
        end
    end
    % end of density calculation======================================================================
    
    % begin: viscosity()==============================================================================
    % Purpose: Calculate the viscosity in the fluid as a function of temperature. 
    % Max error in the actual temp. interval is 0.5%   
    mu(1:NPI+2,2:NPJ+1) = (0.0395*T(1:NPI+2,2:NPJ+1) + 6.58)*1.E-6;   
    % end of viscosity calculation======================================================================
    
    % begin: conductivity()===========================================================================
    % Purpose: Calculate the thermal conductivity in the fluid as a function of temperature.    
    for I = 1:NPI+2
        for J = 2:NPJ+1
            Gamma(I,J) = (6.1E-5*T(I,J) + 8.4E-3)/Cp(I,J);
            if (Gamma(I,J) < 0.)
                output();
                fprintf('Error: Gamma(%d,%d) = %e\n', I, J, Gamma(I,J));
                exit(1);
            end
        end
    end
    % end of thermal conductivity calculation========================================================
    
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



