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
global NPI NPJ XMAX YMAX LARGE U_IN

% variables
global x x_u y y_v u v pc p T rho mu Gamma Cp aP aE aW aN aS b d_u d_v  SMAX SAVG relax_rho 
global Q_chip x_chip_start x_chip_end y_chip_start y_chip_end heat_zone Ti Cmu

heat_zone = struct('x_start', {}, 'x_end', {}, 'q_wall', {}, 'R_copper', {});
    
% constants
NPI        = 4*48;        % number of grid cells in x-direction [-]
NPJ        = 4*24;        % number of grid cells in y-direction [-]
XMAX       = 0.15;      % width of the domain [m]
YMAX       = 0.05;      % height of the domain [m]
MAX_ITER   = 100;      % maximum number of outer iterations [-]
U_ITER     = 1;         % number of Newton iterations for u equation [-]
V_ITER     = 1;         % number of Newton iterations for u equation [-]
PC_ITER    = 200;       % number of Newton iterations for pc equation [-]
T_ITER     = 1;         % number of Newton iterations for T equation [-]
SMAXneeded = 1E-7;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-8;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
P_ATM      = 101000.;   % atmospheric pressure [Pa]
U_IN       = 0.05;      % in flow velocity [m/s]
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
P_left   = 50;           % W - VRM/memory left
P_core   = 300;          % W - GPU core
P_right  = 50;           % W - VRM/memory right

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
    % For liquid water, density is constant at 1000 kg/m^3
    rho(:,:) = 1000.0;
    % end of density calculation======================================================================
    
    % begin: viscosity()==============================================================================
    % Constant viscosity for water
    mu(1:NPI+2,2:NPJ+1) = 1.0E-3;   
    % end of viscosity calculation======================================================================
    
% begin: conductivity()===========================================================================
    % Purpose: Calculate thermal conductivity (Water in channel, Solid Copper in walls and fins)
    h_base_frac = 2/10;
    l_base_frac = 3/10;
    L_triangle = ceil(0.05*(NPI+1));
    Start_L_base = ceil(l_base_frac*(NPI+1));
    End_limit = ceil((1 - l_base_frac)*(NPI+1));   
    H_domain = (NPJ+1);
    Start_H_bottom = ceil(h_base_frac*H_domain);
    Start_H_top = H_domain - Start_H_bottom;
    H_triangle = ceil((1/4)*h_base_frac * H_domain);   % = 1/10 * H_domain
    slope = H_triangle / L_triangle;

    for I = 1:NPI+2
        for J = 2:NPJ+1
            % Default: Fluid (water) conductivity
            Gamma(I,J) = 0.6 / Cp(I,J); 
            
            is_solid = false;
            
            % Check if cell falls inside the solid upper/lower wall boundaries
            if (J < ceil(h_base_frac*(NPJ+1))) || (J > ceil((1-h_base_frac)*(NPJ+1)))
                is_solid = true;
            end
            
            % Check if cell falls inside the solid triangle mesh/fins
            for offset = 0:L_triangle:(End_limit - Start_L_base - L_triangle)
                Start_L_triangle = Start_L_base + offset;
                End_L_triangle = Start_L_triangle + L_triangle;
                
                if (I >= Start_L_triangle && I <= End_L_triangle)
                    i_shift = I - Start_L_triangle;
                    lower_line = ceil(-i_shift*slope + H_triangle + Start_H_bottom);
                    upper_line = ceil(-i_shift*slope + Start_H_top);
                    
                    mid       = floor((lower_line + upper_line) / 2);
                    band_half = floor((upper_line - lower_line) / 6);

                    lower_zigzag1 = lower_line + band_half;
                    upper_zigzag1 = lower_line + 2*band_half;

                    lower_zigzag2 = upper_line - 2*band_half;
                    upper_zigzag2 = upper_line - band_half;

                    if (J < lower_line) || (J > upper_line) || ...
                       (J > lower_zigzag1 && J < upper_zigzag1) || ...
                       (J > lower_zigzag2 && J < upper_zigzag2)
                        is_solid = true;
                    end
                end
            end
            
            % If cell is solid (wall or mesh fin), overwrite with Copper thermal properties
            if is_solid
                Gamma(I,J) = 401.0 / Cp(I,J); % Copper thermal conductivity: 401.0 W/m·K
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



%%
function [] = Tcoeff()
% Purpose: To calculate the coefficients for the T equation.

% constants
global NPI NPJ 
% variables
global x x_u y y_v T Gamma SP Su F_u F_v relax_T Istart Iend Jstart Jend ...
    b aE aW aN aS aP heat_zone Cp

Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;
        % Geometrical parameters: Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % = A(i,J) See fig. 6.2 or fig. 6.5
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); % = A(I,j)
        AREAn = AREAs;
        
        % The convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.    
        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;
        
        % The transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition        
        % The conductivity, Gamma, at the interface is calculated with the use of a harmonic mean.        
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
        
        % --- Volumetric heat source: CPU/GPU chip ---
       % Apply ONLY to the bottom-most fluid cell (J = ceil((NPJ+1)/6))
        if J == ceil((NPJ+1)/6)
            for idx = 1:length(heat_zone)
                if (x(I) >= heat_zone(idx).x_start && x(I) <= heat_zone(idx).x_end)
                    % Use cell width (dx) for a boundary surface flux, not cell volume
                    dx = x_u(i+1) - x_u(i);
                    Su(I,J) = Su(I,J) + (heat_zone(idx).q_wall * dx) / Cp(I,J);
                end
            end
        end
        % The coefficients (hybrid differencing scheme)
        aW(I,j) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,j) = max([-Fe, De - Fe/2, 0.]);
        aS(I,j) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,j) = max([-Fn, Dn - Fn/2, 0.]);
        
        
        % transport of T through the baffles can be switched off by setting the coefficients to zero

        %lower walls: 
         if (J < ceil((NPJ+1)/6)) 
            aE(I,j) = 0;
            
            aW(I,j) = 0;
            aN(I,j) = 0;
         end
        %upper walls:
        if (J > ceil(5*(NPJ+1)/6)) 
            aE(I,j) = 0;
            aS(I,j) = 0;
            aW(I,j) = 0;
            
         end
        
        % =================================================================
        % --- Conduction through internal Baffles ---

        k_baffle = 401.0; % Copper: 401.0 W/m·K (or 205.0 for Aluminum)
        Gamma_baffle = k_baffle / Cp(I,J);
        
        % Baffle #1 (Lower baffle)
        if (I == ceil((NPI+1)/5)-1 && J < ceil((NPJ+1)/3))     % left of baffle #1
            aE(I,J) = (Gamma_baffle * AREAe) / (x(I+1) - x(I));
        end       
        if (I == ceil((NPI+1)/5)   && J < ceil((NPJ+1)/3))     % right of baffle #1
            aW(I,J) = (Gamma_baffle * AREAw) / (x(I) - x(I-1));
        end
        
        % Baffle #2 (Upper baffle)
        if (I == ceil(2*(NPI+1)/5)-1 && J > ceil(2*(NPJ+1)/3)) % left of baffle #2
            aE(I,J) = (Gamma_baffle * AREAe) / (x(I+1) - x(I));
        end       
        if (I == ceil(2*(NPI+1)/5)   && J > ceil(2*(NPJ+1)/3)) % right of baffle #2
            aW(I,J) = (Gamma_baffle * AREAw) / (x(I) - x(I-1));
        end
        % =================================================================
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
        
        % Setting the source term equal to b       
        b(I,J) = Su(I,J);
        
        % Introducing relaxation by eq. 6.36 . and putting also the last
        % term on the right side into the source term b(i,J)        
        aP(I,J) = aP(I,J) / relax_T;
        b (I,J) = b (I,J) + (1 - relax_T)*aP(I,J)*T(I,J);
        
        % now the TDMA algorithm can be called to solve the equation.
        % This is done in the next step of the main program.
    end
end
end

