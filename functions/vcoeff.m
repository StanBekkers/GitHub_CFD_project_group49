function [] = vcoeff()
% Purpose: To calculate the coefficients for the v equation.

% constants
global NPI NPJ LARGE
% variables
global x x_u y y_v v p mueff SP Su F_u F_v d_v relax_v Istart Iend Jstart Jend ...
    b aE aW aN aS aP

Istart = 2;
Iend = NPI+1;
Jstart = 3;
Jend = NPJ+1;

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
        
        % lower walls: 
        h_base_frac = 2/10;
        l_base_frac = 3/10;
        if (J < ceil(h_base_frac*(NPJ+1))) 
            aW(I,j) = 0; aE(I,j) = 0;
            aS(I,j) = 0; aN(I,j) = 0;
            SP(I,j) = -LARGE;
        end
        % upper walls:
        if (J > ceil((1-h_base_frac)*(NPJ+1))) 
            aW(I,j) = 0; aE(I,j) = 0;
            aS(I,j) = 0; aN(I,j) = 0;
            SP(I,j) = -LARGE;
        end

       L_triangle = ceil(0.05*(NPI+1));
Start_L_base = ceil(l_base_frac*(NPI+1));
End_limit = ceil((1 - l_base_frac)*(NPI+1));   

H_domain = (NPJ+1);
Start_H_bottom = ceil(h_base_frac*H_domain);
Start_H_top = H_domain - Start_H_bottom;
H_triangle = ceil((1/4)*h_base_frac * H_domain);   % = 1/10 * H_domain
slope = H_triangle / L_triangle;

        for offset = 0:L_triangle:(End_limit - Start_L_base - L_triangle)
            Start_L_triangle = Start_L_base + offset;
            End_L_triangle   = Start_L_triangle + L_triangle;

            if (i >= Start_L_triangle) && (i <= End_L_triangle)

    i_shift = i - Start_L_triangle;

    lower_line = ceil(-i_shift*slope + H_triangle + Start_H_bottom);
    upper_line = ceil(-i_shift*slope + Start_H_top);

mid       = floor((lower_line + upper_line) / 2);
band_half = floor((upper_line - lower_line) / 6);

lower_zigzag1 = lower_line + band_half;
upper_zigzag1 = lower_line + 2*band_half;

lower_zigzag2 = upper_line - 2*band_half;
upper_zigzag2 = upper_line - band_half;

    if (J < lower_line)
        aW(I,j) = 0; aE(I,j) = 0;
        aS(I,j) = 0; aN(I,j) = 0;
        SP(I,j) = -LARGE;
    end
    if (J > upper_line)
        aW(I,j) = 0; aE(I,j) = 0;
        aS(I,j) = 0; aN(I,j) = 0;
        SP(I,j) = -LARGE;
    end
    if (J > lower_zigzag1 && J < upper_zigzag1)
        aW(I,j) = 0; aE(I,j) = 0;
        aS(I,j) = 0; aN(I,j) = 0;
        SP(I,j) = -LARGE;
    end
    if (J > lower_zigzag2 && J < upper_zigzag2)
        aW(I,j) = 0; aE(I,j) = 0;
        aS(I,j) = 0; aN(I,j) = 0;
        SP(I,j) = -LARGE;
    end

end
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