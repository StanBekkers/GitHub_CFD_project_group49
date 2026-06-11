function [] = Tcoeff()
% Purpose: To calculate the coefficients for the T equation.

% constants
global NPI NPJ YMAX
% variables
global x x_u y y_v T Gamma SP Su F_u F_v relax_T Istart Iend Jstart Jend ...
    b aE aW aN aS aP Cp

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
        AREAw = y_v(j+1) - y_v(j); 
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); 
        AREAn = AREAs;
        
        % The convective mass flux
        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;
        
        % The transport by diffusion (harmonic mean handles the fluid-solid interface)
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
        
        % Note: Manual SP and Su boundary overrides at J=2 and J=NPJ+1 have 
        % been removed. The boundary temperatures are implicitly balanced 
        % and updated in bound.m, allowing conduction to carry physical heat 
        % transfer natively through neighboring cell coupling.
        
        % The coefficients (hybrid differencing scheme)
        aW(I,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,J) = max([-Fe, De - Fe/2, 0.]);
        aS(I,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,J) = max([-Fn, Dn - Fn/2, 0.]);
        
        % eq. 8.31 without time dependent terms:
        mass_imbalance = Fe - Fw + Fn - Fs;
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + max(0, mass_imbalance) - SP(I,J);
        b(I,J)  = Su(I,J) + max(0, -mass_imbalance) * T(I,J);
        
        % Introducing relaxation
        aP(I,J) = aP(I,J) / relax_T;
        b(I,J)  = b(I,J) + (1 - relax_T)*aP(I,J)*T(I,J);
    end
end
end