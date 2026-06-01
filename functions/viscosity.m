function [] = viscosity()
% Purpose: To calculate the viscosity in the fluid as a function of temperature.

% constants
global NPI NPJ Cmu SMALL
% variables
global rho k mu mut mueff eps

for I = 1:NPI+1
    for J = 2: NPJ+2
        mut(I,J)   = rho(I,J)*Cmu*k(I,J)^2 /(eps(I,J)+SMALL);
        mueff(I,J)  = mu(I,J) + mut(I,J);
    end
end
end