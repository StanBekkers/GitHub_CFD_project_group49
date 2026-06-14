function [] = vcoeff()
% Purpose: To calculate the coefficients for the v equation.

% constants
global NPI NPJ LARGE
% variables
global x x_u y y_v v p mueff SP Su F_u F_v d_v relax_v Istart Iend Jstart Jend ...
    b aE aW aN aS aP
global h_base_frac l_base_frac

Istart = 2;
Iend = NPI+1;
Jstart = 3;
Jend = NPJ+1;

h_base_frac = 2/10;
l_base_frac = 3/10;

layout_wall = Walls(Istart, Iend, Jstart, Jend, NPI, NPJ, h_base_frac);
layout_fins = TriangleFin(Istart, Iend, Jstart, Jend, NPI, NPJ, l_base_frac, h_base_frac);
cooler_layout = layout_wall | layout_fins;

convect();
for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;        
        % Geometrical parameters: Areas of the cell faces
        AREAw = y(J) - y(J-1);
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i);
        AREAn = AREAs;
        
        % eq. 6.11a-6.11d
        Fw = ((F_u(i,J)   + F_u(i,J-1))/2)*AREAw;
        Fe = ((F_u(i+1,J) + F_u(i+1,J-1))/2)*AREAe;
        Fs = ((F_v(I,j)   + F_v(I,j-1))/2)*AREAs;
        Fn = ((F_v(I,j)   + F_v(I,j+1))/2)*AREAn;
        
        % eq. 6.11e-6.11h
        Dw = ((mueff(I-1,J-1) + mueff(I,J-1) + mueff(I-1,J) + mueff(I,J))/(4*(x(I) - x(I-1))))*AREAw;
        De = ((mueff(I,J-1) + mueff(I+1,J-1) + mueff(I,J) + mueff(I+1,J))/(4*(x(I+1) - x(I))))*AREAe;
        Ds =  (mueff(I,J-1)/(y_v(j) - y_v(j-1)))*AREAs;
        Dn =  (mueff(I,J)/(y_v(j+1) - y_v(j)))*AREAn;
        
        % The source terms
        SP(I,j) = 0.;
        Su(I,j) = 0.;
        
        % The coefficients (hybrid differencing scheme)
        aW(I,j) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,j) = max([-Fe, De - Fe/2, 0.]);
        aS(I,j) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,j) = max([-Fn, Dn - Fn/2, 0.]);

        if (cooler_layout(i,j) == 1)
           aW(I,j) = 0; aE(I,j) = 0;
           aS(I,j) = 0; aN(I,j) = 0;
           SP(I,j) = -LARGE;
        end

        % eq. 8.31 without time dependent terms
        aP(I,j) = aW(I,j) + aE(I,j) + aS(I,j) + aN(I,j) + Fe - Fw + Fn - Fs - SP(I,J);
        
        d_v(I,j) = AREAs*relax_v/aP(I,j);
        
        b(I,j) = (p(I,J-1) - p(I,J))*AREAs + Su(I,j);
        
        aP(I,j) = aP(I,j)/relax_v;
        b(I,j)  = b(I,j) + (1 - relax_v)*aP(I,j)*v(I,j);
    end
end
end