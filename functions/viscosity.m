function [] = viscosity()
% Purpose: Calculate turbulent and effective viscosity using k-epsilon model.
% mut  = rho * Cmu * k^2 / (eps + SMALL)
% mueff = mu + mut   (used in momentum, k, eps diffusion)
% Gamma is also updated with turbulent thermal diffusivity: (mu + mut/Pr_t)/Cp

% constants
global NPI NPJ Cmu SMALL
% variables
global rho k eps mu mut mueff Gamma Cp

Pr_t = 0.9;   % turbulent Prandtl number for air

for I = 1:NPI+2
    for J = 1:NPJ+2
        mut(I,J)   = rho(I,J) * Cmu * k(I,J)^2 / (eps(I,J) + SMALL);
        mueff(I,J) = mu(I,J) + mut(I,J);
        % Effective thermal conductivity = (lambda + mu_t*Cp/Pr_t) / Cp
        % Since Gamma = lambda/Cp in the base code, update consistently:
        Gamma(I,J) = (mu(I,J) + mut(I,J)/Pr_t) / Cp(I,J);
    end
end
end