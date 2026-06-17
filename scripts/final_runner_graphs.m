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

% --- GLOBAL VARIABLES ---
global NPI NPJ XMAX YMAX LARGE SMALL U_IN
global x x_u y y_v u v pc p T rho mu Gamma Cp aP aE aW aN aS b d_u d_v  SMAX SAVG relax_rho 
global Q_chip x_chip_start x_chip_end y_chip_start y_chip_end heat_zone Ti Cmu
global k eps mut mueff uplus yplus yplus1 yplus2 tw k_old eps_old
global dudx dudy dvdx dvdy E E2 t_fluid
global sigmak sigmaeps C1eps C2eps kappa ERough relax_k relax_eps P_core
global h_base_frac l_base_frac cooler_layout J_fluid_bottom J_fluid_top

% --- DIRECTORY SETUP ---
fprintf('Current Working Directory: %s\n', pwd);

results_dir = 'results_folder';
if ~exist(results_dir, 'dir')
    [status, msg] = mkdir(results_dir);
    if ~status
        % Fallback mechanism: If MATLAB cannot create the folder due to 
        % permissions, fall back to current directory to prevent a crash.
        warning('Failed to create directory "%s". Reason: %s.\nFalling back to saving files directly in current folder.', results_dir, msg);
        results_dir = '.'; 
    end
end

% --- PARAMETER SWEEP VARIABLES ---
% Added 'nofins' to establish a clear heat and flow baseline
%fin_types = {'triangle', 'rectangular', 'elipse', 'nofins'};
fin_types = {"nofins"};
velocities = [0.05, 0.20, 0.50]; % Decided velocities: Low, Nominal, High

% Set to false to run the automated sweep in the background and save to disk
% Set to true if you want MATLAB to display all plots on your screen
show_live_plots = false; 

if show_live_plots
    fig_state = 'on';
else
    fig_state = 'off';
end

% --- CONSTANTS ---
NPI        = 4*48;      % number of grid cells in x-direction [-]
NPJ        = 4*24;      % number of grid cells in y-direction [-]
XMAX       = 0.15;      % width of the domain [m]
YMAX       = 0.05;      % height of the domain [m]
MAX_ITER   = 1500;      % maximum number of outer iterations [-]
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
NPRINT     = 100;       % number of iterations between printing output to screen

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
h_base_frac = 2/10;      % Height fraction of the baseplate walls
l_base_frac = 3/10;      % Length fraction where baffles reside

% Central heat zone calculations
A_core   = ((1 - 2*l_base_frac) * XMAX) * 0.015;     % m² active area
P_core   = 450;                                      % Total power in Watts
q_flux_core  = P_core  / A_core;
R_core   = t_copper / (k_copper * A_core);

heat_zone = struct('x_start', {}, 'x_end', {}, 'q_wall', {}, 'R_copper', {});
heat_zone(1) = struct('x_start', l_base_frac*XMAX,    'x_end', (1-l_base_frac)*XMAX, ...
                      'q_wall',  q_flux_core,         'R_copper', R_core);

% =========================================================================
% %% MAIN AUTOMATED PARAMETER SWEEP LOOP
% =========================================================================

for f_idx = 1:length(fin_types)
    fin_type = fin_types{f_idx};
    
    for v_idx = 1:length(velocities)
        U_IN = velocities(v_idx);
        
        fprintf('\n=======================================================================\n');
        fprintf('  LAUNCHING SOLVER | Fin Type: %s | Inlet Velocity: %0.2f m/s\n', upper(fin_type), U_IN);
        fprintf('=======================================================================\n');
        
        % 1. Re-initialize memory arrays for this specific run
        init(); 
        
        % Construct filenames & output paths (includes wattage)
        suffix = sprintf('U_%0.2f_W_%g_fin_%s', U_IN, P_core, fin_type);
        
        % 2. Rebuild the physical solid layout
        I_full = 1; Iend_full = NPI+2;
        J_full = 1; Jend_full = NPJ+2;
        layout_wall = Walls(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, h_base_frac);

        if strcmp(fin_type, 'triangle')
            layout_fins = TriangleFin(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, l_base_frac, h_base_frac);
        elseif strcmp(fin_type, 'rectangular')
            layout_fins = RectangularFin(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, l_base_frac, h_base_frac);
        elseif strcmp(fin_type, 'elipse')
            layout_fins = ElipseFin(I_full, Iend_full, J_full, Jend_full, NPI, NPJ, l_base_frac, h_base_frac);
        elseif strcmp(fin_type, 'nofins')
            % Flat, straight channel baseline with no internal solid geometry
            layout_fins = zeros(Iend_full, Jend_full);
        else
            error('Invalid fin_type specified. Select either ''triangle'', ''rectangular'', ''elipse'', or ''nofins''.');
        end
        cooler_layout = layout_wall | layout_fins;

        % 3. Generate Baseline (No Flow) Conduction Profile once per geometry
        if v_idx == 1
            T_mock = 293.15 * ones(NPI+2, NPJ+2);
            for I = 1:NPI+2
                for J = J_fluid_bottom:J_fluid_top
                    dx = x(I) - 0.5 * XMAX;
                    dy = y(J) - 0.5 * YMAX;
                    T_mock(I, J) = 293.15 + 3.35 * exp(-(dx^2 / (2 * 0.028^2) + dy^2 / (2 * 0.010^2)));
                end
            end

            fig1 = figure('Visible', fig_state); 
            clf(fig1);
            imagesc(x, y, T_mock')
            set(gca, 'YDir', 'normal')
            hold on
            contour(x, y, double(cooler_layout)', [0.5 0.5], 'k-', 'LineWidth', 1.5)
            hold off
            colorbar
            clim([292 296.5])
            xlabel('x [m]'); ylabel('y [m]')
            title(sprintf('Baseline Temperature Profile (No Flow, %g W - %s)', P_core, upper(fin_type)))
            saveas(fig1, fullfile(results_dir, sprintf('Baseline_Temperature_NoFlow_W_%g_fin_%s.png', P_core, fin_type)));
            if ~show_live_plots, close(fig1); end
        end

        % 4. Main SIMPLE algorithm iterations
        iter = 1;
        SMAX = LARGE;
        SAVG = LARGE;
        
        while (iter <= MAX_ITER && (SMAX > SMAXneeded || SAVG > SAVGneeded))
            bound(); 
            derivatives();
            kcoeff();
            for iter_k = 1:K_ITER
                k = solve(k, b, aE, aW, aN, aS, aP);
            end
            k = max(k, 1e-10);   

            epscoeff();
            for iter_eps = 1:EPS_ITER
                eps = solve(eps, b, aE, aW, aN, aS, aP);
            end
            eps = max(eps, 1e-10); 

            rho(:,:) = 1000.0;
            mu(1:NPI+2, 2:NPJ+1) = 1.0E-3;
            viscosity();   
            
            ucoeff(); 
            for iter_u = 1:U_ITER
                u = solve(u, b, aE, aW, aN, aS, aP); 
            end
            
            vcoeff(); 
            for iter_v = 1:V_ITER
                v = solve(v, b, aE, aW, aN, aS, aP); 
            end
            
            bound(); 
            pccoeff(); 
            for iter_pc = 1:PC_ITER
                pc = solve(pc, b, aE, aW, aN, aS, aP); 
            end
            
            velcorr(); 
            p_outlet_avg = mean(p(NPI+1, J_fluid_bottom:J_fluid_top));
            p = p - p_outlet_avg;
            
            Tcoeff(); 
            for iter_T = 1:T_ITER
                T = solve(T, b, aE, aW, aN, aS, aP); 
            end
           
            if mod(iter, NPRINT) == 0
                I = round((NPI+1)/2);
                J = round((NPJ+1)/2);
                du = d_u(I,J)*(pc(I-1,J) - pc(I,J));
                dv = d_v(I,J)*(pc(I,J-1) - pc(I,J));
                fprintf ('Iter %4d:\t du/u=%10.2e\t dv/v=%10.2e\t SMAX=%10.2e\t SAVG=%10.2e\n', iter, du/u(I,J), dv/v(I,J), SMAX, SAVG);
            end
            iter = iter + 1;   
        end

        % 5. Save ALL Numerical Output Files with dynamic suffix
        fp   = fopen(fullfile(results_dir, sprintf('output_%s.dat', suffix)), 'w');
        str  = fopen(fullfile(results_dir, sprintf('str_%s.dat', suffix)), 'w');
        velu = fopen(fullfile(results_dir, sprintf('velu_%s.dat', suffix)), 'w');
        velv = fopen(fullfile(results_dir, sprintf('velv_%s.dat', suffix)), 'w');

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

        % 6. Render and Save Plots dynamically
        [X, Y] = meshgrid(x, y);

        % FIG 2: Velocity vector profile
        fig2 = figure(2); clf(fig2); set(fig2, 'Visible', fig_state);
        quiver(X, Y, u', v', 1.5);
        xlabel('x [m]'); ylabel('y [m]');
        title(sprintf('Velocity Vector Profile (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        saveas(fig2, fullfile(results_dir, sprintf('Velocity_Vectors_%s.png', suffix)));
        if ~show_live_plots, close(fig2); end

        % FIG 3: Temperature profile
        fig3 = figure(3); clf(fig3); set(fig3, 'Visible', fig_state);
        imagesc(x, y, T')
        set(gca, 'YDir', 'normal')
        colorbar;
        xlabel('x [m]'); ylabel('y [m]');
        title(sprintf('Temperature [K] (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        saveas(fig3, fullfile(results_dir, sprintf('Simulated_Temperature_%s.png', suffix)));
        if ~show_live_plots, close(fig3); end

        % FIG 4: Pressure profile
        fig4 = figure(4); clf(fig4); set(fig4, 'Visible', fig_state);
        imagesc(x, y, p')
        set(gca, 'YDir', 'normal')
        colorbar;
        xlabel('x [m]'); ylabel('y [m]');
        title(sprintf('Pressure [Pa] (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        saveas(fig4, fullfile(results_dir, sprintf('Pressure_%s.png', suffix)));
        if ~show_live_plots, close(fig4); end

        % FIG 5: Turbulent kinetic energy k
        fig5 = figure(5); clf(fig5); set(fig5, 'Visible', fig_state);
        imagesc(x, y, k')
        set(gca, 'YDir', 'normal')
        colorbar; colormap(fig5, jet);
        xlabel('x [m]'); ylabel('y [m]')
        title(sprintf('Turbulent Kinetic Energy k [m^2/s^2] (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        saveas(fig5, fullfile(results_dir, sprintf('Turbulent_Kinetic_Energy_%s.png', suffix)));
        if ~show_live_plots, close(fig5); end

        % FIG 6: Turbulent dissipation rate epsilon
        fig6 = figure(6); clf(fig6); set(fig6, 'Visible', fig_state);
        imagesc(x, y, eps')
        set(gca, 'YDir', 'normal')
        colorbar; colormap(fig6, jet);
        xlabel('x [m]'); ylabel('y [m]')
        title(sprintf('Turbulent Dissipation Rate \\epsilon [m^2/s^3] (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        saveas(fig6, fullfile(results_dir, sprintf('Turbulent_Dissipation_%s.png', suffix)));
        if ~show_live_plots, close(fig6); end

        % FIG 7: Turbulent viscosity ratio mut/mu
        fig7 = figure(7); clf(fig7); set(fig7, 'Visible', fig_state);
        imagesc(x, y, (mut ./ (mu + 1e-30))')
        set(gca, 'YDir', 'normal')
        colorbar; colormap(fig7, hot);
        xlabel('x [m]'); ylabel('y [m]')
        title(sprintf('Turbulent Viscosity Ratio \\mu_t / \\mu [-] (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        saveas(fig7, fullfile(results_dir, sprintf('Turbulent_Viscosity_Ratio_%s.png', suffix)));
        if ~show_live_plots, close(fig7); end

        % FIG 8: y+ distribution (bottom wall)
        J_bot = ceil(h_base_frac*(NPJ+1));
        fig8 = figure(8); clf(fig8); set(fig8, 'Visible', fig_state);
        plot(x(2:NPI+1), yplus(2:NPI+1, J_bot), 'b-', 'LineWidth', 1.5)
        hold on
        yline(11.63, 'r--', 'Sublayer limit y^+=11.63')
        yline(300,   'k--', 'Log-law upper limit y^+=300')
        hold off
        xlabel('x [m]'); ylabel('y^+')
        title(sprintf('Wall y^+ at Bottom Channel Wall (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        legend('y^+', 'Sublayer limit', 'Log-law upper limit')
        grid on
        saveas(fig8, fullfile(results_dir, sprintf('Wall_yplus_%s.png', suffix)));
        if ~show_live_plots, close(fig8); end

        % FIG 9: Velocity magnitude contour
        u_grid = zeros(NPI+1, NPJ);
        v_grid = zeros(NPI+1, NPJ);
        for I = 1:NPI+1
            for J = 2:NPJ+1
                u_grid(I, J-1) = 0.5*(u(I,J) + u(I+1,J));
                v_grid(I, J-1) = 0.5*(v(I,J) + v(I,J+1));
            end
        end
        V_mag = sqrt(u_grid.^2 + v_grid.^2);

        fig9 = figure(9); clf(fig9); set(fig9, 'Visible', fig_state);
        contourf(x(1:NPI+1), y(2:NPJ+1), V_mag', 30, 'LineColor', 'none')
        colorbar; colormap(fig9, jet);
        xlabel('x [m]'); ylabel('y [m]')
        title(sprintf('Velocity Magnitude [m/s] (U_{IN} = %0.2f m/s, %g W, %s)', U_IN, P_core, upper(fin_type)));
        set(gca, 'YDir', 'normal')
        saveas(fig9, fullfile(results_dir, sprintf('Velocity_Magnitude_%s.png', suffix)));
        if ~show_live_plots, close(fig9); end
        
        fprintf('  Run completed successfully! Saved figures and .dat files for: %s @ %0.2f m/s\n', upper(fin_type), U_IN);
    end
end

fprintf('\n=======================================================================\n');
fprintf('  ALL RUNS COMPLETED SUCCESSFULLY! BOTH .PNG AND .DAT FILES ARE SAVED.\n');
fprintf('  Check the folder "%s" to access your results.\n', results_dir);
fprintf('=======================================================================\n');