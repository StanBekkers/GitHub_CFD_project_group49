function [] = epscoeff()
% Purpose: Calculate coefficients for the eps (dissipation rate) equation.

% constants
global NPI NPJ Cmu LARGE SMALL sigmaeps kappa C1eps C2eps
% variables
global x x_u y y_v SP Su F_u F_v mut rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP k eps eps_old E2 mu
global h_base_frac l_base_frac

Istart = 2;
Iend   = NPI+1;
Jstart = 2;
Jend   = NPJ+1;

convect();
viscosity();

layout_wall = Walls(Istart, Iend, Jstart, Jend, NPI, NPJ, h_base_frac);
layout_fins = TriangleFin(Istart, Iend, Jstart, Jend, NPI, NPJ, l_base_frac, h_base_frac);
cooler_layout = layout_wall | layout_fins;

%Find interfaces 
di = diff(cooler_layout,1,1);
dj = diff(cooler_layout,1,2);
solid_above = (di == -1);
solid_below = (di == 1);
solid_left  = (dj == -1);
solid_right = (dj == 1);

interfaces = false(size(cooler_layout));
interfaces(1:end-1,:) = (di ~= 0);
interfaces(:,1:end-1) = interfaces(:,1:end-1) | (dj ~= 0);


% Effective turbulent diffusivity for eps (includes molecular viscosity)
Gamma_eps = mu + mut / sigmaeps;

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;

        AREAw = y_v(j+1) - y_v(j);
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i);
        AREAn = AREAs;

        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;

          if (cooler_layout(i,j) == 1)
            is_solid = true;
          end

        if is_solid
            % Enforce eps floor in solids (Bypass relaxation)
            SP(I,J) = -LARGE;
            Su(I,J) = LARGE * 1e-10;
            
            aW(I,J) = 0.0;
            aE(I,J) = 0.0;
            aS(I,J) = 0.0;
            aN(I,J) = 0.0;
            aP(I,J) = LARGE;
            b(I,J)  = Su(I,J);
            continue;
        end

        % Transport by diffusion (harmonic mean of Gamma_eps)
        Dw = (Gamma_eps(I-1,J)*Gamma_eps(I,J)) / ...
             ((Gamma_eps(I-1,J))*(x(I)-x_u(i)) + (Gamma_eps(I,J))*(x_u(i)-x(I-1))) * AREAw;
        De = (Gamma_eps(I,J)*Gamma_eps(I+1,J)) / ...
             ((Gamma_eps(I,J))*(x(I+1)-x_u(i+1)) + (Gamma_eps(I+1,J))*(x_u(i+1)-x(I))) * AREAe;
        Ds = (Gamma_eps(I,J-1)*Gamma_eps(I,J)) / ...
             ((Gamma_eps(I,J-1))*(y(J)-y_v(j)) + (Gamma_eps(I,J))*(y_v(j)-y(J-1))) * AREAs;
        Dn = (Gamma_eps(I,J)*Gamma_eps(I,J+1)) / ...
             ((Gamma_eps(I,J))*(y(J+1)-y_v(j+1)) + (Gamma_eps(I,J+1))*(y_v(j+1)-y(J))) * AREAn;

        % Default: interior source terms
        SP(I,J) = -C2eps * rho(I,J) * eps(I,J) / (k(I,J) + SMALL);
        Su(I,J) =  C1eps * (eps(I,J) / (k(I,J) + SMALL)) * mut(I,J) * E2(I,J); 
        isWall  = false;

        % ---- Channel bottom wall functions ----
        if (solid_above(i,j) == true)
            y_P      = 0.5 * (y(J_fluid_bottom) - y(J_fluid_bottom - 1)); 
            eps_wall = Cmu^0.75 * k(I,J)^1.5 / (kappa * y_P + SMALL);
            SP(I,J)  = -LARGE;
            Su(I,J)  =  LARGE * eps_wall;
            isWall   = true;
        end

        % ---- Channel top wall functions ----
        if (solid_below(i,j) == true)
            y_P      = 0.5 * (y(J_fluid_top + 1) - y(J_fluid_top)); 
            eps_wall = Cmu^0.75 * k(I,J)^1.5 / (kappa * y_P + SMALL);
            SP(I,J)  = -LARGE;
            Su(I,J)  =  LARGE * eps_wall;
            isWall   = true;
        end

        Su(I,J) = Su(I,J) * AREAw * AREAs;
        SP(I,J) = SP(I,J) * AREAw * AREAs;

        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);

        % Standard coefficient formulation
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
        b(I,J)  = Su(I,J);
        
        % ---- Apply under-relaxation for fluid cells ----
        global relax_eps % ensure access to global relax_eps
        aP(I,J) = aP(I,J) / relax_eps;
        b(I,J)  = b(I,J) + (1 - relax_eps) * aP(I,J) * eps(I,J);
    end
end
end