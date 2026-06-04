function [] = ucoeff()
% Purpose: To calculate the coefficients for the u equation.

% constants
global NPI NPJ LARGE 
% variables
global x x_u y y_v u p mu SP Su F_u F_v d_u relax_u Istart Iend Jstart Jend ...
      b aE aW aN aS aP
 
Istart = 3;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;        
        % Geometrical parameters: Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
        AREAe = AREAw;
        AREAs = x(I) - x(I-1);
        AREAn = AREAs;
        
        % eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
        Fw = ((F_u(i,J)   + F_u(i-1,J))/2)*AREAw;
        Fe = ((F_u(i+1,J) + F_u(i,J))/2)*AREAe;
        Fs = ((F_v(I,j)   + F_v(I-1,j))/2)*AREAs;
        Fn = ((F_v(I,j+1) + F_v(I-1,j+1))/2)*AREAn;
        
        % eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        Dw = (mu(I-1,J)/(x_u(i) - x_u(i-1)))*AREAw;
        De = (mu(I,J)/(x_u(i+1) - x_u(i)))*AREAe;
        Ds = ((mu(I-1,J) + mu(I,J) + mu(I-1,J-1) + mu(I,J-1))/(4*(y(J) - y(J-1))))*AREAs;
        Dn = ((mu(I-1,J+1) + mu(I,J+1) + mu(I-1,J) + mu(I,J))/(4*(y(J+1) - y(J))))*AREAn;
        
        % The source terms
        SP(i,J) = 0.;
        Su(i,J) = 0.;
        
         % u can be fixed to zero by setting SP to a very large value

         % lower walls: 
         base_frac = 2/10;
          if (J < ceil(base_frac*(NPJ+1)))
            SP(i,J) = -LARGE;
          end
         % upper walls: 
          if (J > ceil((1-base_frac)*(NPJ+1)))
            SP(i,J) = -LARGE;
          end

L_triangle = ceil(0.05*(NPI+1));
Start_L_base = ceil(base_frac*(NPI+1));
End_limit = ceil((1 - base_frac)*(NPI+1));   

H_domain = (NPJ+1);
Start_H_bottom = ceil(base_frac*H_domain);
Start_H_top = H_domain - Start_H_bottom;
H_triangle = ceil((1/3)*base_frac * H_domain);   % = 1/10 * H_domain
%zigzag_bottom = ceil(0.5*H_triangle);

% --- Loop over repeated triangles ---
for offset = 0:L_triangle:(End_limit - Start_L_base)
    Start_L_triangle = Start_L_base + offset;
    End_L_triangle = Start_L_triangle + L_triangle;
    % local coordinate
    i_shift = i - Start_L_triangle;
    % --- triangle mask ---
    if ((Start_L_triangle <= i) && (i <= End_L_triangle))
      lower_line = ceil((-i_shift*H_triangle/L_triangle) + H_triangle + Start_H_bottom);
      upper_line = ceil((-i_shift*H_triangle/L_triangle) + Start_H_top);
      center_line = ceil(0.5 * (lower_line + upper_line));
      lower_zigzag = center_line - H_triangle;
      upper_zigzag = center_line + H_triangle;

        if (J < lower_line)
            SP(i,J) = -LARGE;
        end
        if (J <  upper_zigzag && J > lower_zigzag)
            SP(i,J) = -LARGE;
        end
        if (J > upper_line)
            SP(i,J) = -LARGE;
        end
        
    end
end

        % The coefficients (hybrid differencing scheme)
        aW(i,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(i,J) = max([-Fe, De - Fe/2, 0.]);
        aS(i,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(i,J) = max([-Fn, Dn - Fn/2, 0.]);
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(i,J) = aW(i,J) + aE(i,J) + aS(i,J) + aN(i,J) + Fe - Fw + Fn - Fs - SP(I,J);
        
        % Calculation of d(i)(J) = d_u(i)(J) defined in eq. 6.23 for use in the
        % equation for pression correction (eq. 6.32). See subroutine pccoeff.
        d_u(i,J) = AREAw*relax_u/aP(i,J);
        
        % Putting the integrated pressure gradient into the source term b(i)(J)
        % The reason is to get an equation on the generalised form
        % (eq. 7.7 ) to be solved by the TDMA algorithm.
        % note: In reality b = a0p*fiP + Su = 0.     
        b(i,J) = (p(I-1,J) - p(I,J))*AREAw + Su(I,J);
        
        % Introducing relaxation by eq. 6.36 . and putting also the last
        % term on the right side into the source term b(i)(J)
        aP(i,J) = aP(i,J) / relax_u;
        b (i,J) = b(i,J) + (1 - relax_u)*aP(i,J)*u(i,J);
        
        % now we have implemented eq. 6.36 in the form of eq. 7.7
        % and the TDMA algorithm can be called to solve it. This is done
        % in the next step of the main program. 
    end
end
end

